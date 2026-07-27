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
    if (!std::isfinite(charge)) throw std::invalid_argument("unstructured deposit charge must be finite");
    if (!std::isfinite(weight) || weight <= 0.0) {
        throw std::invalid_argument("unstructured deposit weight must be positive and finite");
    }
    const double particle_charge = charge * weight;
    if (!std::isfinite(particle_charge)) {
        throw std::invalid_argument("unstructured particle charge product must be finite");
    }

    UnstructuredDepositSummary2D summary;
    for (const auto& particle : particles) {
        if (!particle.alive) continue;
        const auto location = mesh.locate_point(particle.position);
        if (!location) {
            ++summary.outside_particles;
            continue;
        }
        for (std::size_t i = 0; i < location->node_ids.size(); ++i) {
            const std::size_t node_index = mesh.node_index(location->node_ids[i]);
            const double increment =
                particle_charge * location->shape_weights[i] / mesh.node_control_areas()[node_index];
            const double updated = mesh.rho()[node_index] + increment;
            if (!std::isfinite(updated)) {
                throw std::overflow_error("unstructured deposited charge density overflow");
            }
            mesh.rho()[node_index] = updated;
        }
        ++summary.deposited_particles;
        summary.deposited_charge += particle_charge;
        if (!std::isfinite(summary.deposited_charge)) {
            throw std::overflow_error("unstructured deposited charge total overflow");
        }
    }
    return summary;
}

std::optional<Vec2> interpolate_electric(const UnstructuredMesh2D& mesh, Vec2 position) {
    const auto location = mesh.locate_point(position);
    if (!location) return std::nullopt;

    Vec2 result{};
    for (std::size_t i = 0; i < location->node_ids.size(); ++i) {
        const Vec2 value = mesh.electric()[mesh.node_index(location->node_ids[i])];
        if (!std::isfinite(value.x) || !std::isfinite(value.y)) {
            throw std::runtime_error("unstructured electric field values must be finite");
        }
        result.x += location->shape_weights[i] * value.x;
        result.y += location->shape_weights[i] * value.y;
    }
    return result;
}

} // namespace pic
