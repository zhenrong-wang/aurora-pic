#include "pic/UnstructuredSimulation2D.hpp"

#include "pic/Convergence.hpp"
#include "pic/Pusher.hpp"
#include "pic/VTKWriter.hpp"

#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <set>
#include <stdexcept>
#include <utility>

namespace pic {
namespace {

using SteadyClock = std::chrono::steady_clock;
constexpr const char* CHECKPOINT_MAGIC = "AuroraPIC-unstructured-2D-checkpoint-v1";

double cross(Vec2 first, Vec2 second) {
    return first.x * second.y - first.y * second.x;
}

Vec2 subtract(Vec2 first, Vec2 second) {
    return {first.x - second.x, first.y - second.y};
}

Vec2 add(Vec2 first, Vec2 second) {
    return {first.x + second.x, first.y + second.y};
}

Vec2 scale(Vec2 value, double factor) {
    return {factor * value.x, factor * value.y};
}

double dot(Vec2 first, Vec2 second) {
    return first.x * second.x + first.y * second.y;
}

double triangle_area(Vec2 first, Vec2 second, Vec2 third) {
    return 0.5 * std::abs(cross(subtract(second, first), subtract(third, first)));
}

bool cell_has_edge(const ImportedCell2D& cell, std::size_t first, std::size_t second) {
    for (std::size_t i = 0; i < cell.node_ids.size(); ++i) {
        const auto a = cell.node_ids[i];
        const auto b = cell.node_ids[(i + 1) % cell.node_ids.size()];
        if ((a == first && b == second) || (a == second && b == first)) return true;
    }
    return false;
}

Species2DConfig particle_storage_config(const UnstructuredSpecies2DConfig& config) {
    Species2DConfig result;
    result.name = config.name;
    result.charge = config.charge;
    result.mass = config.mass;
    result.weight = config.weight;
    result.particles = config.particles;
    result.drift_velocity_x = config.drift_velocity_x;
    result.drift_velocity_y = config.drift_velocity_y;
    result.thermal_velocity = config.thermal_velocity;
    return result;
}

void initialize_particle_pusher(Particle2D& particle, Vec2 electric, double charge_to_mass,
                                double magnetic_field_z, double dt) {
    if (magnetic_field_z == 0.0) {
        initialize_leapfrog_half_step(particle, electric, charge_to_mass, dt);
    } else {
        initialize_boris_half_step(particle, electric, magnetic_field_z, charge_to_mass, dt);
    }
}

void kick_particle(Particle2D& particle, Vec2 electric, double charge_to_mass,
                   double magnetic_field_z, double dt) {
    if (magnetic_field_z == 0.0) {
        kick_leapfrog(particle, electric, charge_to_mass, dt);
    } else {
        kick_boris(particle, electric, magnetic_field_z, charge_to_mass, dt);
    }
}

void synchronize_particle(Particle2D& particle, Vec2 electric, double charge_to_mass,
                          double magnetic_field_z, double dt) {
    if (magnetic_field_z == 0.0) {
        synchronize_leapfrog(particle, electric, charge_to_mass, dt);
    } else {
        synchronize_boris(particle, electric, magnetic_field_z, charge_to_mass, dt);
    }
}

std::string csv_quote(const std::string& value) {
    std::string result{"\""};
    for (const char character : value) {
        if (character == '"') result += '"';
        result += character;
    }
    result += '"';
    return result;
}

ImportedMesh2D load_configured_mesh(const std::filesystem::path& path) {
    if (path.empty()) {
        throw std::invalid_argument("unstructured simulation mesh path must not be empty");
    }
    return load_gmsh2_ascii_mesh2d(path);
}

std::size_t particle_worker_count(
    std::size_t particle_count, const RuntimePolicy& runtime) {
    return particle_count == 0
               ? 0
               : std::min(runtime_info(runtime).active_threads, particle_count);
}

template <typename Body>
void parallel_particle_chunks(std::vector<Particle2D>& particles,
                              const RuntimePolicy& runtime,
                              Body&& body) {
    const std::size_t workers = particle_worker_count(particles.size(), runtime);
    if (workers == 0) return;
    std::vector<std::exception_ptr> failures(workers);
    runtime_parallel_for(std::size_t{0}, workers, runtime,
                         [&](std::size_t worker) {
        try {
            const std::size_t begin = particles.size() * worker / workers;
            const std::size_t end = particles.size() * (worker + 1) / workers;
            for (std::size_t particle = begin; particle < end; ++particle) {
                body(worker, particles[particle]);
            }
        } catch (...) {
            failures[worker] = std::current_exception();
        }
    });
    for (const auto& failure : failures) {
        if (failure) std::rethrow_exception(failure);
    }
}

} // namespace

UnstructuredSimulation2D::UnstructuredSimulation2D(UnstructuredSimulation2DConfig config)
    : config_(std::move(config)),
      mesh_(load_configured_mesh(config_.mesh_path)),
      rng_(config_.seed) {
    if (!std::isfinite(config_.dt) || config_.dt <= 0.0) {
        throw std::invalid_argument("unstructured simulation dt must be positive and finite");
    }
    if (config_.output_interval == 0) {
        throw std::invalid_argument("unstructured simulation output_interval must be positive");
    }
    if (config_.particle_output_stride == 0) {
        throw std::invalid_argument("unstructured particle output stride must be positive");
    }
    if (config_.checkpoint_output && config_.checkpoint_interval == 0) {
        config_.checkpoint_interval = config_.output_interval;
    }
    if (!std::isfinite(config_.magnetic_field_z)) {
        throw std::invalid_argument("unstructured magnetic_field_z must be finite");
    }
    if (config_.mode == RunMode::SteadyState) {
        if (config_.max_steps == 0 || config_.steady_window == 0 ||
            !std::isfinite(config_.steady_tolerance) || config_.steady_tolerance <= 0.0) {
            throw std::invalid_argument("invalid unstructured steady-state convergence configuration");
        }
    }
    validate_runtime_policy(config_.runtime);

    const auto labels = mesh_.topology().boundary_labels();
    const std::set<std::string> label_set(labels.begin(), labels.end());
    if (config_.particle_boundaries.size() != label_set.size()) {
        throw std::invalid_argument("particle boundary policies must cover every imported boundary label");
    }
    for (const auto& [label, policy] : config_.particle_boundaries) {
        if (!label_set.contains(label)) {
            throw std::invalid_argument("particle boundary label not found in imported mesh: " + label);
        }
        if (policy != ParticleBoundary::Absorbing && policy != ParticleBoundary::Reflecting) {
            throw std::invalid_argument(
                "unstructured particle boundaries currently support absorbing or reflecting policies");
        }
        absorbed_by_label_[label] = 0;
    }

    for (const auto& face : mesh_.topology().boundary_faces()) {
        const Vec2 first = mesh_.topology().node_by_id(face.node_ids[0]).position;
        const Vec2 second = mesh_.topology().node_by_id(face.node_ids[1]).position;
        const ImportedCell2D* adjacent = nullptr;
        for (const auto& cell : mesh_.topology().cells()) {
            if (cell_has_edge(cell, face.node_ids[0], face.node_ids[1])) {
                adjacent = &cell;
                break;
            }
        }
        if (!adjacent) throw std::runtime_error("boundary face has no adjacent imported cell");
        const Vec2 midpoint = scale(add(first, second), 0.5);
        const Vec2 toward_cell = subtract(mesh_.topology().cell_centroid(adjacent->id), midpoint);
        const Vec2 tangent = subtract(second, first);
        Vec2 normal{-tangent.y, tangent.x};
        const double length = std::hypot(normal.x, normal.y);
        if (!(length > 0.0)) throw std::runtime_error("degenerate imported boundary segment");
        normal = scale(normal, 1.0 / length);
        if (dot(normal, toward_cell) < 0.0) normal = scale(normal, -1.0);
        boundary_segments_.push_back({first, second, normal, face.label});
    }

    double cumulative_area = 0.0;
    for (const auto& cell : mesh_.topology().cells()) {
        const auto position = [&](std::size_t local) {
            return mesh_.topology().node_by_id(cell.node_ids[local]).position;
        };
        const auto append_triangle = [&](Vec2 first, Vec2 second, Vec2 third) {
            const double area = triangle_area(first, second, third);
            if (!(area > 0.0) || !std::isfinite(area)) {
                throw std::runtime_error("invalid sampling triangle in imported cell");
            }
            cumulative_area += area;
            if (!std::isfinite(cumulative_area)) {
                throw std::overflow_error("unstructured sampling area overflow");
            }
            sampling_triangles_.push_back({{first, second, third}, cumulative_area});
        };
        append_triangle(position(0), position(1), position(2));
        if (cell.shape == ImportedCellShape2D::Quadrilateral) {
            append_triangle(position(0), position(2), position(3));
        }
    }

    if (config_.species.empty()) {
        config_.species.push_back(UnstructuredSpecies2DConfig{});
    }
    std::set<std::string> species_names;
    for (const auto& species_config : config_.species) {
        if (!species_names.insert(species_config.name).second) {
            throw std::invalid_argument("unstructured species names must be unique");
        }
        if (species_config.initialization_minimum.has_value() !=
            species_config.initialization_maximum.has_value()) {
            throw std::invalid_argument(
                "unstructured species initialization bounds require both minimum and maximum");
        }
        if (species_config.initialization_minimum) {
            const Vec2 minimum = *species_config.initialization_minimum;
            const Vec2 maximum = *species_config.initialization_maximum;
            if (!std::isfinite(minimum.x) || !std::isfinite(minimum.y) ||
                !std::isfinite(maximum.x) || !std::isfinite(maximum.y) ||
                !(minimum.x < maximum.x) || !(minimum.y < maximum.y)) {
                throw std::invalid_argument("unstructured species initialization bounds are invalid");
            }
        }
        species_configs_.push_back(species_config);
        species_.emplace_back(particle_storage_config(species_config));
    }
    poisson_solver_ = std::make_unique<UnstructuredPoissonSolver2D>(
        mesh_, config_.dirichlet_potentials,
        config_.neumann_normal_derivatives, config_.poisson);
}

Vec2 UnstructuredSimulation2D::sample_position(const UnstructuredSpecies2DConfig& config) {
    if (sampling_triangles_.empty()) throw std::runtime_error("unstructured domain has no sampling area");
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    const double total_area = sampling_triangles_.back().cumulative_area;
    for (std::size_t attempt = 0; attempt < 100000; ++attempt) {
        const double target = unit(rng_) * total_area;
        const auto triangle = std::lower_bound(
            sampling_triangles_.begin(), sampling_triangles_.end(), target,
            [](const SamplingTriangle& candidate, double value) {
                return candidate.cumulative_area < value;
            });
        const auto selected =
            triangle == sampling_triangles_.end() ? std::prev(sampling_triangles_.end()) : triangle;
        const double root = std::sqrt(unit(rng_));
        const double second = unit(rng_);
        const std::array<double, 3> weights{
            1.0 - root,
            root * (1.0 - second),
            root * second,
        };
        Vec2 point{};
        for (std::size_t i = 0; i < weights.size(); ++i) {
            point.x += weights[i] * selected->vertices[i].x;
            point.y += weights[i] * selected->vertices[i].y;
        }
        if (config.initialization_minimum) {
            const Vec2 minimum = *config.initialization_minimum;
            const Vec2 maximum = *config.initialization_maximum;
            if (point.x < minimum.x || point.x > maximum.x ||
                point.y < minimum.y || point.y > maximum.y) {
                continue;
            }
        }
        return point;
    }
    throw std::runtime_error(
        "could not sample the requested unstructured species initialization region");
}

void UnstructuredSimulation2D::initialize() {
    time_ = 0.0;
    step_ = 0;
    timing_ = {};
    for (auto& [label, count] : absorbed_by_label_) count = 0;
    for (std::size_t species_id = 0; species_id < species_.size(); ++species_id) {
        auto& particles = species_[species_id].particles();
        const auto& species_config = species_configs_[species_id];
        particles.assign(species_config.particles, Particle2D{});
        std::normal_distribution<double> velocity_x(
            species_config.drift_velocity_x, species_config.thermal_velocity);
        std::normal_distribution<double> velocity_y(
            species_config.drift_velocity_y, species_config.thermal_velocity);
        for (auto& particle : particles) {
            particle.position = sample_position(species_config);
            particle.velocity = {velocity_x(rng_), velocity_y(rng_)};
            particle.velocity_half = particle.velocity;
            particle.alive = true;
        }
    }
    deposit_and_solve();
    const auto particle_start = SteadyClock::now();
    for (auto& species : species_) {
        const double charge_to_mass = species.charge() / species.mass();
        parallel_particle_chunks(
            species.particles(), config_.runtime,
            [&](std::size_t, Particle2D& particle) {
                const auto electric =
                    interpolate_electric(mesh_, particle.position);
                if (!electric) {
                    throw std::runtime_error(
                        "initialized particle is outside imported mesh");
                }
                initialize_particle_pusher(
                    particle, *electric, charge_to_mass,
                    config_.magnetic_field_z, config_.dt);
            });
    }
    timing_.particle_seconds +=
        std::chrono::duration<double>(SteadyClock::now() - particle_start).count();
    initialized_ = true;
}

void UnstructuredSimulation2D::deposit_and_solve() {
    const auto deposition_start = SteadyClock::now();
    mesh_.clear_charge();
    for (const auto& species : species_) {
        const auto deposit =
            deposit_charge_shape(
                mesh_, species.particles(), species.charge(),
                species.weight(), config_.runtime);
        if (deposit.outside_particles != 0) {
            throw std::runtime_error("live particles remain outside the imported mesh during deposition");
        }
    }
    timing_.deposition_seconds +=
        std::chrono::duration<double>(SteadyClock::now() - deposition_start).count();
    const auto solve_start = SteadyClock::now();
    last_poisson_ = poisson_solver_->solve(mesh_);
    timing_.field_solve_seconds +=
        std::chrono::duration<double>(SteadyClock::now() - solve_start).count();
    if (!last_poisson_.converged) {
        throw std::runtime_error("unstructured Poisson solver did not converge");
    }
}

void UnstructuredSimulation2D::advance_with_boundaries(
    Particle2D& particle, Vec2 previous_position,
    std::map<std::string, std::size_t>& absorbed_by_label) {
    Vec2 start = previous_position;
    Vec2 end = particle.position;
    const Vec2 minimum = mesh_.topology().min_corner();
    const Vec2 maximum = mesh_.topology().max_corner();
    const double domain_scale = std::max(maximum.x - minimum.x, maximum.y - minimum.y);
    const double tolerance =
        512.0 * std::numeric_limits<double>::epsilon() * domain_scale;

    for (int collision = 0; collision < 32; ++collision) {
        const Vec2 displacement = subtract(end, start);
        double earliest = std::numeric_limits<double>::infinity();
        const BoundarySegment* hit_segment = nullptr;
        for (const auto& segment : boundary_segments_) {
            if (dot(displacement, segment.inward_normal) >= 0.0) continue;
            const Vec2 edge = subtract(segment.second, segment.first);
            const double denominator = cross(displacement, edge);
            if (std::abs(denominator) <=
                128.0 * std::numeric_limits<double>::epsilon() *
                    std::max(1.0, std::hypot(displacement.x, displacement.y) *
                                      std::hypot(edge.x, edge.y))) {
                continue;
            }
            const Vec2 offset = subtract(segment.first, start);
            double along_path = cross(offset, edge) / denominator;
            const double along_edge = cross(offset, displacement) / denominator;
            if (along_path < -1e-12 || along_path > 1.0 + 1e-12 ||
                along_edge < -1e-12 || along_edge > 1.0 + 1e-12) {
                continue;
            }
            along_path = std::clamp(along_path, 0.0, 1.0);
            if (along_path < earliest) {
                earliest = along_path;
                hit_segment = &segment;
            }
        }
        if (!hit_segment) {
            if (mesh_.locate_point(end)) {
                particle.position = end;
                return;
            }
            throw std::runtime_error(
                "particle left imported domain without a resolvable boundary intersection");
        }

        const ParticleBoundary policy = config_.particle_boundaries.at(hit_segment->label);
        if (policy == ParticleBoundary::Absorbing) {
            particle.position = add(start, scale(displacement, earliest));
            particle.alive = false;
            ++absorbed_by_label.at(hit_segment->label);
            return;
        }

        const Vec2 hit = add(start, scale(displacement, earliest));
        const Vec2 remaining = scale(displacement, 1.0 - earliest);
        const Vec2 reflected_remaining =
            subtract(remaining, scale(hit_segment->inward_normal,
                                      2.0 * dot(remaining, hit_segment->inward_normal)));
        particle.velocity_half =
            subtract(particle.velocity_half,
                     scale(hit_segment->inward_normal,
                           2.0 * dot(particle.velocity_half, hit_segment->inward_normal)));
        start = add(hit, scale(hit_segment->inward_normal, tolerance));
        end = add(start, reflected_remaining);
    }
    throw std::runtime_error("particle exceeded the imported-boundary reflection limit");
}

void UnstructuredSimulation2D::step() {
    if (!initialized_) initialize();
    const auto first_particle_start = SteadyClock::now();
    for (auto& species : species_) {
        const double charge_to_mass = species.charge() / species.mass();
        const std::size_t workers =
            particle_worker_count(species.particles().size(), config_.runtime);
        std::vector<std::map<std::string, std::size_t>> local_absorbed(workers);
        for (auto& counts : local_absorbed) {
            for (const auto& [label, total] : absorbed_by_label_) {
                (void)total;
                counts.emplace(label, 0);
            }
        }
        parallel_particle_chunks(
            species.particles(), config_.runtime,
            [&](std::size_t worker, Particle2D& particle) {
                if (!particle.alive) return;
                const auto electric =
                    interpolate_electric(mesh_, particle.position);
                if (!electric) {
                    throw std::runtime_error(
                        "live particle is outside imported mesh before push");
                }
                kick_particle(
                    particle, *electric, charge_to_mass,
                    config_.magnetic_field_z, config_.dt);
                const Vec2 previous = particle.position;
                drift_leapfrog(particle, config_.dt);
                advance_with_boundaries(
                    particle, previous, local_absorbed[worker]);
            });
        for (const auto& counts : local_absorbed) {
            for (const auto& [label, count] : counts) {
                absorbed_by_label_.at(label) += count;
            }
        }
    }
    timing_.particle_seconds +=
        std::chrono::duration<double>(SteadyClock::now() - first_particle_start).count();

    deposit_and_solve();
    const auto second_particle_start = SteadyClock::now();
    for (auto& species : species_) {
        const double charge_to_mass = species.charge() / species.mass();
        parallel_particle_chunks(
            species.particles(), config_.runtime,
            [&](std::size_t, Particle2D& particle) {
                if (!particle.alive) return;
                const auto electric =
                    interpolate_electric(mesh_, particle.position);
                if (!electric) {
                    throw std::runtime_error(
                        "live particle is outside imported mesh after push");
                }
                synchronize_particle(
                    particle, *electric, charge_to_mass,
                    config_.magnetic_field_z, config_.dt);
            });
    }
    timing_.particle_seconds +=
        std::chrono::duration<double>(SteadyClock::now() - second_particle_start).count();
    ++step_;
    time_ += config_.dt;
}

UnstructuredDiagnosticSample2D UnstructuredSimulation2D::sample() const {
    UnstructuredDiagnosticSample2D result;
    result.step = step_;
    result.time = time_;
    result.absorbed_by_label = absorbed_by_label_;
    result.poisson = last_poisson_;
    for (const auto& species : species_) {
        result.kinetic_energy += species.kinetic_energy();
        result.live_particles += species.live_count();
    }
    for (std::size_t i = 0; i < mesh_.size(); ++i) {
        const double area = mesh_.node_control_areas()[i];
        const Vec2 electric = mesh_.electric()[i];
        result.field_energy +=
            0.5 * EPS0 * (electric.x * electric.x + electric.y * electric.y) * area;
        result.charge_l1 += std::abs(mesh_.rho()[i]) * area;
    }
    result.total_energy = result.kinetic_energy + result.field_energy;
    return result;
}

void UnstructuredSimulation2D::write_diagnostics_header(std::ofstream& output) const {
    output << "step,time,kinetic_energy,field_energy,total_energy,charge_l1,live_particles"
              ",poisson_iterations,poisson_initial_residual,poisson_final_residual"
              ",particle_seconds,deposition_seconds,field_solve_seconds";
    for (const auto& [label, count] : absorbed_by_label_) {
        (void)count;
        output << ',' << csv_quote("absorbed_" + label);
    }
    output << '\n';
}

void UnstructuredSimulation2D::write_diagnostics_sample(
    std::ofstream& output, const UnstructuredDiagnosticSample2D& value) const {
    output << value.step << ',' << std::setprecision(17) << value.time << ','
           << value.kinetic_energy << ',' << value.field_energy << ','
           << value.total_energy << ',' << value.charge_l1 << ','
           << value.live_particles << ',' << value.poisson.iterations << ','
           << value.poisson.initial_residual << ',' << value.poisson.final_residual << ','
           << timing_.particle_seconds << ',' << timing_.deposition_seconds << ','
           << timing_.field_solve_seconds;
    for (const auto& [label, count] : value.absorbed_by_label) {
        (void)label;
        output << ',' << count;
    }
    output << '\n';
    output.flush();
}

void UnstructuredSimulation2D::write_particle_sample(std::size_t step) const {
    std::ofstream output(config_.output_dir / ("particles_" + std::to_string(step) + ".csv"));
    if (!output) throw std::runtime_error("cannot open unstructured particle output");
    output << "species_id,species,x,y,vx,vy,alive\n";
    std::size_t written = 0;
    for (std::size_t species_id = 0; species_id < species_.size(); ++species_id) {
        const auto& species = species_[species_id];
        for (std::size_t particle_id = 0;
             particle_id < species.particles().size(); ++particle_id) {
            if (particle_id % config_.particle_output_stride != 0) continue;
            const auto& particle = species.particles()[particle_id];
            output << species_id << ',' << csv_quote(species.name()) << ','
                   << std::setprecision(17)
                   << particle.position.x << ',' << particle.position.y << ','
                   << particle.velocity.x << ',' << particle.velocity.y << ','
                   << (particle.alive ? 1 : 0) << '\n';
            ++written;
            if (config_.particle_sample_count != 0 &&
                written >= config_.particle_sample_count) {
                return;
            }
        }
    }
}

std::filesystem::path UnstructuredSimulation2D::checkpoint_path_for_step(
    std::size_t step) const {
    if (!config_.checkpoint_path.empty()) return config_.checkpoint_path;
    return config_.output_dir / ("checkpoint_" + std::to_string(step) + ".apc");
}

std::uint64_t UnstructuredSimulation2D::mesh_signature() const {
    constexpr std::uint64_t offset_basis = 14695981039346656037ULL;
    constexpr std::uint64_t prime = 1099511628211ULL;
    std::uint64_t hash = offset_basis;
    const auto append = [&](std::uint64_t value) {
        for (int byte = 0; byte < 8; ++byte) {
            hash ^= (value >> (8 * byte)) & 0xffU;
            hash *= prime;
        }
    };
    const auto append_string = [&](const std::string& value) {
        for (const unsigned char character : value) {
            hash ^= character;
            hash *= prime;
        }
        hash ^= 0xffU;
        hash *= prime;
    };
    for (const auto& node : mesh_.topology().nodes()) {
        append(node.id);
        append(std::bit_cast<std::uint64_t>(node.position.x));
        append(std::bit_cast<std::uint64_t>(node.position.y));
    }
    for (const auto& cell : mesh_.topology().cells()) {
        append(cell.id);
        append(static_cast<std::uint64_t>(cell.shape));
        append(static_cast<std::uint64_t>(cell.physical_tag));
        append_string(cell.label);
        for (const auto node_id : cell.node_ids) append(node_id);
    }
    for (const auto& face : mesh_.topology().boundary_faces()) {
        append(face.id);
        append(static_cast<std::uint64_t>(face.physical_tag));
        append_string(face.label);
        append(face.node_ids[0]);
        append(face.node_ids[1]);
    }
    return hash;
}

void UnstructuredSimulation2D::save_checkpoint(const std::filesystem::path& path) const {
    const auto parent = path.parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("cannot open unstructured checkpoint for writing: " +
                                 path.string());
    }
    output << std::setprecision(17);
    output << CHECKPOINT_MAGIC << '\n';
    output << "mesh_signature " << mesh_signature() << '\n';
    output << "step " << step_ << '\n';
    output << "time " << time_ << '\n';
    output << "absorbed_count " << absorbed_by_label_.size() << '\n';
    for (const auto& [label, count] : absorbed_by_label_) {
        output << "absorbed " << std::quoted(label) << ' ' << count << '\n';
    }
    output << "species_count " << species_.size() << '\n';
    output << "rng " << rng_ << '\n';
    for (std::size_t species_id = 0; species_id < species_.size(); ++species_id) {
        const auto& species = species_[species_id];
        output << "species " << species_id << ' ' << std::quoted(species.name()) << ' '
               << species.particles().size() << '\n';
        for (const auto& particle : species.particles()) {
            output << particle.position.x << ' ' << particle.position.y << ' '
                   << particle.velocity.x << ' ' << particle.velocity.y << ' '
                   << particle.velocity_half.x << ' ' << particle.velocity_half.y << ' '
                   << (particle.alive ? 1 : 0) << '\n';
        }
    }
    if (!output) {
        throw std::runtime_error("failed while writing unstructured checkpoint: " +
                                 path.string());
    }
}

void UnstructuredSimulation2D::load_checkpoint(const std::filesystem::path& path) {
    timing_ = {};
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open unstructured checkpoint for reading: " +
                                 path.string());
    }
    std::string magic;
    std::getline(input, magic);
    if (magic != CHECKPOINT_MAGIC) {
        throw std::runtime_error("invalid unstructured checkpoint magic");
    }
    std::string key;
    std::uint64_t signature = 0;
    input >> key >> signature;
    if (key != "mesh_signature" || signature != mesh_signature()) {
        throw std::runtime_error("unstructured checkpoint mesh does not match configured geometry");
    }
    input >> key >> step_;
    if (key != "step") throw std::runtime_error("unstructured checkpoint missing step");
    input >> key >> time_;
    if (key != "time" || !std::isfinite(time_) || time_ < 0.0) {
        throw std::runtime_error("unstructured checkpoint has invalid time");
    }
    std::size_t absorbed_count = 0;
    input >> key >> absorbed_count;
    if (key != "absorbed_count" || absorbed_count != absorbed_by_label_.size()) {
        throw std::runtime_error("unstructured checkpoint boundary-label count mismatch");
    }
    std::set<std::string> loaded_labels;
    for (std::size_t i = 0; i < absorbed_count; ++i) {
        std::string label;
        std::size_t count = 0;
        input >> key >> std::quoted(label) >> count;
        if (key != "absorbed" || !absorbed_by_label_.contains(label) ||
            !loaded_labels.insert(label).second) {
            throw std::runtime_error("unstructured checkpoint boundary labels do not match");
        }
        absorbed_by_label_[label] = count;
    }
    std::size_t species_count = 0;
    input >> key >> species_count;
    if (key != "species_count" || species_count != species_.size()) {
        throw std::runtime_error("unstructured checkpoint species count mismatch");
    }
    input >> key;
    if (key != "rng") throw std::runtime_error("unstructured checkpoint missing RNG state");
    input >> rng_;
    for (std::size_t species_id = 0; species_id < species_.size(); ++species_id) {
        std::size_t stored_id = 0;
        std::string stored_name;
        std::size_t particle_count = 0;
        input >> key >> stored_id >> std::quoted(stored_name) >> particle_count;
        if (key != "species" || stored_id != species_id ||
            stored_name != species_[species_id].name() ||
            particle_count != species_configs_[species_id].particles) {
            throw std::runtime_error("unstructured checkpoint species metadata mismatch");
        }
        auto& particles = species_[species_id].particles();
        particles.resize(particle_count);
        for (auto& particle : particles) {
            int alive = 0;
            input >> particle.position.x >> particle.position.y
                  >> particle.velocity.x >> particle.velocity.y
                  >> particle.velocity_half.x >> particle.velocity_half.y >> alive;
            particle.alive = alive != 0;
            if (!std::isfinite(particle.position.x) ||
                !std::isfinite(particle.position.y) ||
                !std::isfinite(particle.velocity.x) ||
                !std::isfinite(particle.velocity.y) ||
                !std::isfinite(particle.velocity_half.x) ||
                !std::isfinite(particle.velocity_half.y) ||
                (particle.alive && !mesh_.locate_point(particle.position))) {
                throw std::runtime_error("unstructured checkpoint contains invalid particle state");
            }
        }
    }
    if (!input) throw std::runtime_error("truncated unstructured checkpoint");
    deposit_and_solve();
    initialized_ = true;
}

UnstructuredRunSummary2D UnstructuredSimulation2D::run() {
    if (config_.restart_path.empty()) {
        initialize();
    } else {
        load_checkpoint(config_.restart_path);
    }
    std::filesystem::create_directories(config_.output_dir);
    std::ofstream diagnostics(config_.output_dir / "scalars.csv");
    if (!diagnostics) throw std::runtime_error("cannot open unstructured diagnostics output");
    write_diagnostics_header(diagnostics);

    std::vector<UnstructuredDiagnosticSample2D> history;
    UnstructuredRunSummary2D summary;
    auto initial = sample();
    history.push_back(initial);
    write_diagnostics_sample(diagnostics, initial);
    summary.final_sample = initial;
    if (config_.vtk_output) {
        write_vtk_xml(
            mesh_, config_.output_dir / ("fields_" + std::to_string(step_) + ".vtu"));
    }
    if (config_.particle_output) write_particle_sample(step_);
    if (config_.checkpoint_output) save_checkpoint(checkpoint_path_for_step(step_));

    const std::size_t limit =
        config_.mode == RunMode::SteadyState ? config_.max_steps : config_.steps;
    const std::size_t particle_interval =
        config_.particle_output_interval == 0
            ? config_.output_interval
            : config_.particle_output_interval;
    while (step_ < limit) {
        step();
        if (step_ % config_.output_interval != 0 && step_ != limit) continue;
        auto current = sample();
        history.push_back(current);
        write_diagnostics_sample(diagnostics, current);
        summary.final_sample = current;
        if (config_.vtk_output) {
            write_vtk_xml(
                mesh_, config_.output_dir / ("fields_" + std::to_string(step_) + ".vtu"));
        }
        const bool reached_steady =
            config_.mode == RunMode::SteadyState &&
            adjacent_energy_windows_converged(
                history, config_.steady_window, config_.steady_tolerance);
        if (config_.particle_output &&
            (step_ % particle_interval == 0 || step_ == limit || reached_steady)) {
            write_particle_sample(step_);
        }
        if (config_.checkpoint_output &&
            (step_ % config_.checkpoint_interval == 0 ||
             step_ == limit || reached_steady)) {
            save_checkpoint(checkpoint_path_for_step(step_));
        }
        if (reached_steady) {
            summary.steady_state_reached = true;
            break;
        }
    }
    summary.steps_completed = step_;
    summary.final_time = time_;
    if (summary.final_sample.step != step_) summary.final_sample = sample();
    return summary;
}

} // namespace pic
