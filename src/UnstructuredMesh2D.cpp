#include "pic/UnstructuredMesh2D.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace pic {
namespace {

std::array<double, 4> quadrilateral_shape_weights(double xi, double eta) {
    return {
        0.25 * (1.0 - xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 + eta),
        0.25 * (1.0 - xi) * (1.0 + eta),
    };
}

std::array<Vec2, 4> quadrilateral_shape_derivatives(double xi, double eta) {
    return {
        Vec2{-0.25 * (1.0 - eta), -0.25 * (1.0 - xi)},
        Vec2{0.25 * (1.0 - eta), -0.25 * (1.0 + xi)},
        Vec2{0.25 * (1.0 + eta), 0.25 * (1.0 + xi)},
        Vec2{-0.25 * (1.0 + eta), 0.25 * (1.0 - xi)},
    };
}

std::size_t spatial_coordinate(double value, double minimum, double maximum, std::size_t bins) {
    if (value <= minimum) return 0;
    if (value >= maximum) return bins - 1;
    const double fraction = (value - minimum) / (maximum - minimum);
    return std::min(static_cast<std::size_t>(fraction * static_cast<double>(bins)), bins - 1);
}

} // namespace

UnstructuredMesh2D::UnstructuredMesh2D(ImportedMesh2D topology)
    : topology_(std::move(topology)),
      node_control_areas_(topology_.nodes().size(), 0.0),
      rho_(topology_.nodes().size(), 0.0),
      phi_(topology_.nodes().size(), 0.0),
      electric_(topology_.nodes().size(), Vec2{}) {
    topology_.validate();
    minimum_corner_ = topology_.min_corner();
    maximum_corner_ = topology_.max_corner();
    for (std::size_t i = 0; i < topology_.nodes().size(); ++i) {
        const bool inserted = node_indices_.emplace(topology_.nodes()[i].id, i).second;
        if (!inserted) throw std::runtime_error("unstructured mesh contains duplicate node ids");
    }

    for (const auto& cell : topology_.cells()) {
        if (cell.shape == ImportedCellShape2D::Triangle) {
            const double contribution = topology_.cell_area(cell.id) / 3.0;
            for (const auto node_id : cell.node_ids) {
                node_control_areas_[node_index(node_id)] += contribution;
            }
            continue;
        }

        std::array<Vec2, 4> nodes{};
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            nodes[i] = topology_.node_by_id(cell.node_ids[i]).position;
        }
        const double point = 1.0 / std::sqrt(3.0);
        for (const double xi : {-point, point}) {
            for (const double eta : {-point, point}) {
                const auto weights = quadrilateral_shape_weights(xi, eta);
                const auto derivatives = quadrilateral_shape_derivatives(xi, eta);
                double dx_dxi = 0.0;
                double dx_deta = 0.0;
                double dy_dxi = 0.0;
                double dy_deta = 0.0;
                for (std::size_t i = 0; i < nodes.size(); ++i) {
                    dx_dxi += derivatives[i].x * nodes[i].x;
                    dx_deta += derivatives[i].y * nodes[i].x;
                    dy_dxi += derivatives[i].x * nodes[i].y;
                    dy_deta += derivatives[i].y * nodes[i].y;
                }
                const double jacobian =
                    std::abs(dx_dxi * dy_deta - dx_deta * dy_dxi);
                if (!(jacobian > 0.0) || !std::isfinite(jacobian)) {
                    throw std::runtime_error("invalid quadrilateral Jacobian in imported cell " +
                                             std::to_string(cell.id));
                }
                for (std::size_t i = 0; i < nodes.size(); ++i) {
                    node_control_areas_[node_index(cell.node_ids[i])] += weights[i] * jacobian;
                }
            }
        }
    }

    for (const double area : node_control_areas_) {
        if (!(area > 0.0) || !std::isfinite(area)) {
            throw std::runtime_error("unstructured mesh node control areas must be positive and finite");
        }
    }

    const std::size_t cell_count = topology_.cells().size();
    spatial_bins_x_ = std::max<std::size_t>(
        1, static_cast<std::size_t>(std::ceil(std::sqrt(static_cast<double>(cell_count)))));
    spatial_bins_y_ = 1 + (cell_count - 1) / spatial_bins_x_;
    if (spatial_bins_x_ > std::numeric_limits<std::size_t>::max() / spatial_bins_y_) {
        throw std::overflow_error("unstructured spatial index size overflow");
    }
    spatial_cell_ids_.resize(spatial_bins_x_ * spatial_bins_y_);
    for (const auto& cell : topology_.cells()) {
        Vec2 cell_minimum{std::numeric_limits<double>::infinity(),
                          std::numeric_limits<double>::infinity()};
        Vec2 cell_maximum{-std::numeric_limits<double>::infinity(),
                          -std::numeric_limits<double>::infinity()};
        for (const auto node_id : cell.node_ids) {
            const Vec2 position = topology_.node_by_id(node_id).position;
            cell_minimum.x = std::min(cell_minimum.x, position.x);
            cell_minimum.y = std::min(cell_minimum.y, position.y);
            cell_maximum.x = std::max(cell_maximum.x, position.x);
            cell_maximum.y = std::max(cell_maximum.y, position.y);
        }
        const std::size_t first_i = spatial_coordinate(
            cell_minimum.x, minimum_corner_.x, maximum_corner_.x, spatial_bins_x_);
        const std::size_t last_i = spatial_coordinate(
            cell_maximum.x, minimum_corner_.x, maximum_corner_.x, spatial_bins_x_);
        const std::size_t first_j = spatial_coordinate(
            cell_minimum.y, minimum_corner_.y, maximum_corner_.y, spatial_bins_y_);
        const std::size_t last_j = spatial_coordinate(
            cell_maximum.y, minimum_corner_.y, maximum_corner_.y, spatial_bins_y_);
        const std::size_t covered_i = last_i - first_i + 1;
        const std::size_t covered_j = last_j - first_j + 1;
        constexpr std::size_t MAX_BINS_PER_CELL = 64;
        if (covered_i > MAX_BINS_PER_CELL / covered_j) {
            spatial_global_cell_ids_.push_back(cell.id);
            continue;
        }
        for (std::size_t j = first_j; j <= last_j; ++j) {
            for (std::size_t i = first_i; i <= last_i; ++i) {
                spatial_cell_ids_[j * spatial_bins_x_ + i].push_back(cell.id);
            }
        }
    }
}

std::size_t UnstructuredMesh2D::node_index(std::size_t node_id) const {
    const auto it = node_indices_.find(node_id);
    if (it == node_indices_.end()) {
        throw std::out_of_range("unstructured mesh node id not found: " + std::to_string(node_id));
    }
    return it->second;
}

double UnstructuredMesh2D::node_control_area(std::size_t node_id) const {
    return node_control_areas_[node_index(node_id)];
}

std::optional<ImportedPointLocation2D> UnstructuredMesh2D::locate_point(
    Vec2 point, double relative_tolerance) const {
    if (!std::isfinite(point.x) || !std::isfinite(point.y)) {
        throw std::invalid_argument("unstructured mesh point coordinates must be finite");
    }
    if (!std::isfinite(relative_tolerance) || relative_tolerance < 0.0) {
        throw std::invalid_argument("unstructured point-location tolerance must be finite and non-negative");
    }
    const double scale = std::max(maximum_corner_.x - minimum_corner_.x,
                                  maximum_corner_.y - minimum_corner_.y);
    const double absolute_tolerance =
        std::max(relative_tolerance, 128.0 * std::numeric_limits<double>::epsilon()) * scale;
    if (point.x < minimum_corner_.x - absolute_tolerance ||
        point.x > maximum_corner_.x + absolute_tolerance ||
        point.y < minimum_corner_.y - absolute_tolerance ||
        point.y > maximum_corner_.y + absolute_tolerance) {
        return std::nullopt;
    }
    const std::size_t i = spatial_coordinate(
        std::clamp(point.x, minimum_corner_.x, maximum_corner_.x),
        minimum_corner_.x, maximum_corner_.x, spatial_bins_x_);
    const std::size_t j = spatial_coordinate(
        std::clamp(point.y, minimum_corner_.y, maximum_corner_.y),
        minimum_corner_.y, maximum_corner_.y, spatial_bins_y_);
    for (const auto cell_id : spatial_cell_ids_[j * spatial_bins_x_ + i]) {
        if (auto location = topology_.cell_coordinates(cell_id, point, relative_tolerance)) {
            return location;
        }
    }
    for (const auto cell_id : spatial_global_cell_ids_) {
        if (auto location = topology_.cell_coordinates(cell_id, point, relative_tolerance)) {
            return location;
        }
    }
    return std::nullopt;
}

void UnstructuredMesh2D::clear_charge() {
    std::fill(rho_.begin(), rho_.end(), 0.0);
}

UnstructuredDepositSummary2D deposit_charge_shape(UnstructuredMesh2D& mesh,
                                                   const std::vector<Particle2D>& particles,
                                                   double charge, double weight) {
    return deposit_charge_shape(
        mesh, particles, charge, weight, RuntimePolicy{});
}

UnstructuredDepositSummary2D deposit_charge_shape(
    UnstructuredMesh2D& mesh,
    const std::vector<Particle2D>& particles,
    double charge, double weight,
    const RuntimePolicy& runtime) {
    std::vector<UnstructuredParticleLocation2D> locations(particles.size());
    return deposit_charge_shape(
        mesh, particles, charge, weight, runtime, locations);
}

UnstructuredDepositSummary2D deposit_charge_shape(
    UnstructuredMesh2D& mesh,
    const std::vector<Particle2D>& particles,
    double charge, double weight,
    const RuntimePolicy& runtime,
    std::vector<UnstructuredParticleLocation2D>& locations) {
    if (!std::isfinite(charge)) throw std::invalid_argument("unstructured deposit charge must be finite");
    if (!std::isfinite(weight) || weight <= 0.0) {
        throw std::invalid_argument("unstructured deposit weight must be positive and finite");
    }
    const double particle_charge = charge * weight;
    if (!std::isfinite(particle_charge)) {
        throw std::invalid_argument("unstructured particle charge product must be finite");
    }

    validate_runtime_policy(runtime);
    if (locations.size() != particles.size()) {
        throw std::invalid_argument(
            "unstructured particle-location cache size mismatch");
    }
    if (particles.empty()) return {};
    for (const auto& particle : particles) {
        if (particle.alive &&
            (!std::isfinite(particle.position.x) ||
             !std::isfinite(particle.position.y))) {
            throw std::invalid_argument(
                "unstructured deposited particle position must be finite");
        }
    }
    const std::size_t worker_count =
        std::min(runtime_info(runtime).active_threads, particles.size());
    std::vector<std::vector<double>> local_density(
        worker_count, std::vector<double>(mesh.size(), 0.0));
    std::vector<UnstructuredDepositSummary2D> local_summaries(worker_count);
    std::vector<unsigned char> local_overflow(worker_count, 0);
    runtime_parallel_for(std::size_t{0}, worker_count, runtime,
                         [&](std::size_t worker) {
        const std::size_t begin = particles.size() * worker / worker_count;
        const std::size_t end = particles.size() * (worker + 1) / worker_count;
        auto& density = local_density[worker];
        auto& summary = local_summaries[worker];
        for (std::size_t particle_id = begin; particle_id < end; ++particle_id) {
            const auto& particle = particles[particle_id];
            if (!particle.alive) continue;
            auto& cached = locations[particle_id];
            std::optional<ImportedPointLocation2D> resolved;
            if (cached.valid) {
                resolved = mesh.topology().cell_coordinates(
                    cached.location.cell_id, particle.position);
                if (resolved) ++summary.location_cache_hits;
            }
            if (!resolved) {
                ++summary.location_searches;
                resolved = mesh.locate_point(particle.position);
            }
            if (!resolved) {
                cached.valid = false;
                ++summary.outside_particles;
                continue;
            }
            cached.location = std::move(*resolved);
            cached.valid = true;
            for (std::size_t i = 0; i < cached.location.node_ids.size(); ++i) {
                const std::size_t node_index =
                    mesh.node_index(cached.location.node_ids[i]);
                const double increment =
                    particle_charge * cached.location.shape_weights[i] /
                    mesh.node_control_areas()[node_index];
                const double updated = density[node_index] + increment;
                if (!std::isfinite(updated)) {
                    local_overflow[worker] = 1;
                    continue;
                }
                density[node_index] = updated;
            }
            ++summary.deposited_particles;
            summary.deposited_charge += particle_charge;
            if (!std::isfinite(summary.deposited_charge)) {
                local_overflow[worker] = 1;
            }
        }
    });

    UnstructuredDepositSummary2D summary;
    for (std::size_t worker = 0; worker < worker_count; ++worker) {
        if (local_overflow[worker] != 0) {
            throw std::overflow_error("unstructured deposited charge overflow");
        }
        summary.deposited_particles += local_summaries[worker].deposited_particles;
        summary.outside_particles += local_summaries[worker].outside_particles;
        summary.location_cache_hits += local_summaries[worker].location_cache_hits;
        summary.location_searches += local_summaries[worker].location_searches;
        summary.deposited_charge += local_summaries[worker].deposited_charge;
        if (!std::isfinite(summary.deposited_charge)) {
            throw std::overflow_error("unstructured deposited charge total overflow");
        }
    }
    for (std::size_t node = 0; node < mesh.size(); ++node) {
        double updated = mesh.rho()[node];
        for (std::size_t worker = 0; worker < worker_count; ++worker) {
            updated += local_density[worker][node];
        }
        if (!std::isfinite(updated)) {
            throw std::overflow_error("unstructured deposited charge density overflow");
        }
        mesh.rho()[node] = updated;
    }
    return summary;
}

std::optional<Vec2> interpolate_electric(const UnstructuredMesh2D& mesh, Vec2 position) {
    UnstructuredParticleLocation2D location;
    return interpolate_electric(mesh, position, location);
}

std::optional<Vec2> interpolate_electric(
    const UnstructuredMesh2D& mesh, Vec2 position,
    UnstructuredParticleLocation2D& cached, bool* cache_hit) {
    std::optional<ImportedPointLocation2D> location;
    if (cached.valid) {
        location = mesh.topology().cell_coordinates(
            cached.location.cell_id, position);
    }
    if (cache_hit) *cache_hit = location.has_value();
    if (!location) location = mesh.locate_point(position);
    if (!location) {
        cached.valid = false;
        return std::nullopt;
    }
    cached.location = std::move(*location);
    cached.valid = true;
    Vec2 result{};
    for (std::size_t i = 0; i < cached.location.node_ids.size(); ++i) {
        const Vec2 value =
            mesh.electric()[mesh.node_index(cached.location.node_ids[i])];
        if (!std::isfinite(value.x) || !std::isfinite(value.y)) {
            throw std::runtime_error("unstructured electric field values must be finite");
        }
        result.x += cached.location.shape_weights[i] * value.x;
        result.y += cached.location.shape_weights[i] * value.y;
    }
    return result;
}

} // namespace pic
