#include "pic/UnstructuredSimulation2D.hpp"

#include "pic/Convergence.hpp"
#include "pic/ParticleState.hpp"
#include "pic/Pusher.hpp"
#include "pic/Units.hpp"
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
#include <tuple>
#include <utility>

namespace pic {
namespace {

using SteadyClock = std::chrono::steady_clock;
constexpr const char* CHECKPOINT_MAGIC_V1 =
    "AuroraPIC-unstructured-2D-checkpoint-v1";
constexpr const char* CHECKPOINT_MAGIC_V2 =
    "AuroraPIC-unstructured-2D-checkpoint-v2";
constexpr const char* CHECKPOINT_MAGIC_V3 =
    "AuroraPIC-unstructured-2D-checkpoint-v3";
constexpr const char* CHECKPOINT_MAGIC_V4 =
    "AuroraPIC-unstructured-2D-checkpoint-v4";
constexpr const char* CHECKPOINT_MAGIC_V5 =
    "AuroraPIC-unstructured-2D-checkpoint-v5";
constexpr const char* CHECKPOINT_MAGIC_V6 =
    "AuroraPIC-unstructured-2D-checkpoint-v6";

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
    // Species2D is reused as particle storage and requires a nonzero configured
    // capacity, while imported source-driven species may intentionally start empty.
    result.particles = std::max<std::size_t>(1, config.particles);
    result.drift_velocity_x = config.drift_velocity_x;
    result.drift_velocity_y = config.drift_velocity_y;
    result.drift_velocity_z = config.drift_velocity_z;
    result.thermal_velocity = config.thermal_velocity;
    result.initialization = config.initialization;
    return result;
}

bool has_magnetic_field(Vec3 magnetic_field) {
    return magnetic_field.x != 0.0 ||
           magnetic_field.y != 0.0 ||
           magnetic_field.z != 0.0;
}

Vec3 magnetic_field(
    const UnstructuredSimulation2DConfig& config,
    Vec2 position) {
    if (config.magnetic_field_profile) {
        return config.magnetic_field_profile->evaluate(
            {position.x, position.y, 0.0});
    }
    return {
        config.magnetic_field_x,
        config.magnetic_field_y,
        config.magnetic_field_z};
}

std::string csv_quote(const std::string& value);

void add_collision_statistics(
    CollisionDiagnostics& destination,
    const CollisionStepStatistics& source) {
    destination.candidates += source.candidates;
    destination.null_collisions += source.null_collisions;
    if (destination.channel_collisions.size() !=
        source.channel_collisions.size()) {
        throw std::runtime_error(
            "imported MCC channel statistics size mismatch");
    }
    for (std::size_t channel = 0;
         channel < source.channel_collisions.size(); ++channel) {
        destination.channel_collisions[channel] +=
            source.channel_collisions[channel];
    }
}

void clear_collision_counts(CollisionDiagnostics& diagnostics) {
    diagnostics.candidates = 0;
    diagnostics.null_collisions = 0;
    std::fill(
        diagnostics.channel_collisions.begin(),
        diagnostics.channel_collisions.end(), 0);
}

void write_collision_header(
    std::ofstream& output,
    const CollisionDiagnostics& diagnostics) {
    output << "step,time,candidates,null_collisions";
    for (const auto& name : diagnostics.channel_names) {
        output << ',' << csv_quote(name);
    }
    output << ",cumulative_candidates,cumulative_null_collisions";
    for (const auto& name : diagnostics.channel_names) {
        output << ',' << csv_quote("cumulative_" + name);
    }
    output << '\n';
}

void write_collision_sample(
    std::ofstream& output,
    std::size_t step,
    double time,
    const CollisionDiagnostics& interval,
    const CollisionDiagnostics& totals) {
    output << step << ',' << std::setprecision(17) << time << ','
           << interval.candidates << ',' << interval.null_collisions;
    for (const auto count : interval.channel_collisions) {
        output << ',' << count;
    }
    output << ',' << totals.candidates << ','
           << totals.null_collisions;
    for (const auto count : totals.channel_collisions) {
        output << ',' << count;
    }
    output << '\n';
    output.flush();
}

void initialize_particle_pusher(Particle2D& particle, Vec2 electric, double charge_to_mass,
                                Vec3 magnetic_field, double dt) {
    if (!has_magnetic_field(magnetic_field)) {
        initialize_leapfrog_half_step(particle, electric, charge_to_mass, dt);
    } else {
        initialize_boris_half_step(
            particle, electric, magnetic_field, charge_to_mass, dt);
    }
}

void kick_particle(Particle2D& particle, Vec2 electric, double charge_to_mass,
                   Vec3 magnetic_field, double dt) {
    if (!has_magnetic_field(magnetic_field)) {
        kick_leapfrog(particle, electric, charge_to_mass, dt);
    } else {
        kick_boris(
            particle, electric, magnetic_field, charge_to_mass, dt);
    }
}

void synchronize_particle(Particle2D& particle, Vec2 electric, double charge_to_mass,
                          Vec3 magnetic_field, double dt) {
    if (!has_magnetic_field(magnetic_field)) {
        synchronize_leapfrog(particle, electric, charge_to_mass, dt);
    } else {
        synchronize_boris(
            particle, electric, magnetic_field, charge_to_mass, dt);
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

void write_collision_metadata(
    const std::filesystem::path& output_dir,
    const CollisionConfig& config,
    std::uint64_t signature,
    double neutral_velocity_stddev,
    double neutral_speed_limit_sigma) {
    std::ofstream output(output_dir / "collision_data.txt");
    if (!output) {
        throw std::runtime_error(
            "cannot open imported collision metadata output");
    }
    output << std::setprecision(17);
    output << "format 7\n";
    output << "gas " << std::quoted(config.gas_name) << '\n';
    output << "neutral_mass " << config.neutral_mass << '\n';
    output << "neutral_density " << config.neutral_density << '\n';
    output << "neutral_temperature "
           << config.neutral_temperature << '\n';
    output << "neutral_velocity_stddev "
           << neutral_velocity_stddev << '\n';
    output << "neutral_speed_limit_sigma "
           << neutral_speed_limit_sigma << '\n';
    output << "gas_data_file "
           << std::quoted(config.gas_data_file.string()) << '\n';
    output << "gas_data_version "
           << config.gas_data_version << '\n';
    output << "gas_data_units "
           << std::quoted(to_string(config.gas_data_units)) << '\n';
    output << "dataset_id "
           << std::quoted(config.dataset_id) << '\n';
    output << "dataset_version "
           << std::quoted(config.dataset_version) << '\n';
    output << "data_provenance "
           << std::quoted(config.data_provenance) << '\n';
    output << "citation " << std::quoted(config.citation) << '\n';
    output << "retrieved " << std::quoted(config.retrieved) << '\n';
    output << "license " << std::quoted(config.license) << '\n';
    output << "model_signature " << signature << '\n';
    output << "channel_count " << config.channels.size() << '\n';
    for (const auto& channel : config.channels) {
        output << "channel "
               << std::quoted(channel.name) << ' '
               << std::quoted(to_string(channel.process)) << ' '
               << channel.threshold_energy << ' '
               << channel.energy_scale << ' '
               << channel.cross_section_scale << ' '
               << std::quoted(
                      to_string(channel.energy_frame)) << ' '
               << std::quoted(channel.cross_section_file.string()) << ' '
               << std::quoted(channel.secondary_species) << ' '
               << std::quoted(channel.ion_species) << ' '
               << std::quoted(channel.attachment_species) << ' '
               << std::quoted(
                      to_string(channel.angular_scattering)) << ' '
               << channel.mean_cosine_energy_scale << ' '
               << std::quoted(
                      channel.mean_cosine_file.string()) << ' '
               << std::quoted(
                      to_string(channel.ionization_kinematics)) << ' '
               << channel.ionization_ejected_energy_scale << ' '
               << std::quoted(
                      to_string(
                          channel.cross_section_interpolation)) << ' '
               << std::quoted(
                      to_string(channel.inelastic_transform)) << '\n';
    }
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
                body(worker, particle, particles[particle]);
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
    if (!config_.restart_path.empty() &&
        !config_.initial_state_path.empty()) {
        throw std::invalid_argument(
            "unstructured restart_path and initial_state_path are mutually exclusive");
    }
    if (config_.initial_state_signature &&
        config_.initial_state_path.empty()) {
        throw std::invalid_argument(
            "unstructured initial_state_signature requires initial_state_path");
    }
    validate_initialization_acceptance(
        config_.initialization_acceptance,
        "unstructured initialization acceptance config");
    if (!std::isfinite(config_.dt) || config_.dt <= 0.0) {
        throw std::invalid_argument("unstructured simulation dt must be positive and finite");
    }
    if (config_.output_interval == 0) {
        throw std::invalid_argument("unstructured simulation output_interval must be positive");
    }
    if (config_.particle_output_stride == 0) {
        throw std::invalid_argument("unstructured particle output stride must be positive");
    }
    if (config_.max_particles_per_species == 0) {
        throw std::invalid_argument(
            "unstructured max_particles_per_species must be positive");
    }
    if (config_.checkpoint_output && config_.checkpoint_interval == 0) {
        config_.checkpoint_interval = config_.output_interval;
    }
    if (!std::isfinite(config_.magnetic_field_x) ||
        !std::isfinite(config_.magnetic_field_y) ||
        !std::isfinite(config_.magnetic_field_z)) {
        throw std::invalid_argument(
            "unstructured magnetic_field components must be finite");
    }
    if (config_.magnetic_field_profile) {
        if (config_.magnetic_field_x != 0.0 ||
            config_.magnetic_field_y != 0.0 ||
            config_.magnetic_field_z != 0.0) {
            throw std::invalid_argument(
                "unstructured uniform magnetic_field components and magnetic field profile are mutually exclusive");
        }
        const Vec2 minimum = mesh_.topology().min_corner();
        const Vec2 maximum = mesh_.topology().max_corner();
        config_.magnetic_field_profile->validate_domain(
            {minimum.x, minimum.y, 0.0},
            {maximum.x, maximum.y, 0.0},
            "unstructured 2D simulation");
    }
    if (config_.mode == RunMode::SteadyState) {
        if (config_.max_steps == 0 || config_.steady_window == 0 ||
            !std::isfinite(config_.steady_tolerance) || config_.steady_tolerance <= 0.0) {
            throw std::invalid_argument("invalid unstructured steady-state convergence configuration");
        }
    }
    validate_runtime_policy(config_.runtime);
    config_.poisson.permittivity = config_.units.permittivity();

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
        boundary_lengths_[face.label] += length;
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
            auto& region_triangles =
                region_sampling_triangles_[cell.label];
            const double region_area =
                (region_triangles.empty()
                     ? 0.0
                     : region_triangles.back().cumulative_area) +
                area;
            if (!std::isfinite(region_area)) {
                throw std::overflow_error(
                    "unstructured region sampling area overflow");
            }
            region_triangles.push_back(
                {sampling_triangles_.size() - 1, region_area});
            auto bounds = region_sampling_bounds_.try_emplace(
                cell.label,
                std::pair<Vec2, Vec2>{first, first}).first;
            for (const Vec2 vertex : {first, second, third}) {
                bounds->second.first.x =
                    std::min(bounds->second.first.x, vertex.x);
                bounds->second.first.y =
                    std::min(bounds->second.first.y, vertex.y);
                bounds->second.second.x =
                    std::max(bounds->second.second.x, vertex.x);
                bounds->second.second.y =
                    std::max(bounds->second.second.y, vertex.y);
            }
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
        if (species_config.particles > config_.max_particles_per_species) {
            throw std::invalid_argument(
                "unstructured initial particle count exceeds max_particles_per_species");
        }
        if (species_config.initialization_minimum.has_value() !=
            species_config.initialization_maximum.has_value()) {
            throw std::invalid_argument(
                "unstructured species initialization bounds require both minimum and maximum");
        }
        if (!species_config.initialization_region.empty()) {
            if (species_config.initialization_minimum) {
                throw std::invalid_argument(
                    "unstructured species initialization_region cannot be combined with rectangular initialization bounds");
            }
            if (!region_sampling_triangles_.contains(
                    species_config.initialization_region)) {
                throw std::invalid_argument(
                    "unstructured species initialization region not found in imported mesh: " +
                    species_config.initialization_region);
            }
        }
        if (species_config.initialization_minimum) {
            const Vec2 minimum = *species_config.initialization_minimum;
            const Vec2 maximum = *species_config.initialization_maximum;
            if (!std::isfinite(minimum.x) || !std::isfinite(minimum.y) ||
                !std::isfinite(maximum.x) || !std::isfinite(maximum.y) ||
                !(minimum.x < maximum.x) || !(minimum.y < maximum.y)) {
                throw std::invalid_argument("unstructured species initialization bounds are invalid");
            }
            if (species_config.initialization.loading ==
                ParticleLoading::QuietStart) {
                throw std::invalid_argument(
                    "unstructured quiet_start loading does not yet support rectangular initialization bounds");
            }
        }
        validate_particle_initialization(
            species_config.initialization, 3,
            species_config.thermal_velocity,
            "unstructured species '" + species_config.name + "'");
        validate_density_profile(
            species_config.initialization, 2,
            species_config.particles,
            "unstructured species '" + species_config.name + "'");
        species_configs_.push_back(species_config);
        species_.emplace_back(particle_storage_config(species_config));
        particle_locations_.emplace_back();
    }
    for (const auto& species_config : species_configs_) {
        for (const auto& [label, policy] : config_.particle_boundaries) {
            if (policy == ParticleBoundary::Absorbing) {
                impact_flux_[species_config.name].emplace(
                    label, UnstructuredBoundaryFlux2D{});
            }
        }
    }
    if (config_.collisions.enabled) {
        if (config_.collisions.model !=
            CollisionModelKind::NullCollision) {
            throw std::invalid_argument(
                "imported 2D supports only null-collision MCC");
        }
        if (config_.collisions.gas_name.empty() ||
            config_.collisions.data_provenance.empty() ||
            !std::isfinite(config_.collisions.neutral_mass) ||
            !(config_.collisions.neutral_mass > 0.0) ||
            !std::isfinite(config_.collisions.neutral_temperature) ||
            config_.collisions.neutral_temperature < 0.0) {
            throw std::invalid_argument(
                "imported MCC requires gas name, provenance, positive "
                "neutral mass, and non-negative neutral temperature");
        }
        const auto target = std::find_if(
            species_configs_.begin(), species_configs_.end(),
            [&](const auto& species) {
                return species.name == config_.collisions.species;
            });
        if (target == species_configs_.end()) {
            throw std::invalid_argument(
                "imported MCC target species was not found: " +
                config_.collisions.species);
        }
        mcc_species_id_ = static_cast<std::size_t>(
            target - species_configs_.begin());
        if (target->charge == 0.0 &&
            std::any_of(
                config_.collisions.channels.begin(),
                config_.collisions.channels.end(),
                [](const auto& channel) {
                    return channel.process ==
                           CollisionProcessKind::ChargeExchange;
                })) {
            throw std::invalid_argument(
                "charge exchange requires a charged target species");
        }
        mcc_model_ = std::make_unique<NullCollisionModel>(
            config_.collisions, target->mass);
        collision_totals_.channel_names =
            mcc_model_->channel_names();
        collision_totals_.channel_collisions.assign(
            collision_totals_.channel_names.size(), 0);
        collision_interval_.channel_names =
            collision_totals_.channel_names;
        collision_interval_.channel_collisions.assign(
            collision_totals_.channel_names.size(), 0);
        ionization_channels_.resize(
            config_.collisions.channels.size());
        attachment_channels_.resize(
            config_.collisions.channels.size());
        const auto species_id = [&](const std::string& name) {
            return static_cast<std::size_t>(
                std::find_if(
                    species_configs_.begin(), species_configs_.end(),
                    [&](const auto& species) {
                        return species.name == name;
                    }) - species_configs_.begin());
        };
        for (std::size_t channel = 0;
             channel < config_.collisions.channels.size(); ++channel) {
            const auto& channel_config =
                config_.collisions.channels[channel];
            const auto& target = species_configs_[mcc_species_id_];
            if (channel_config.process ==
                CollisionProcessKind::Ionization) {
                const std::size_t secondary =
                    species_id(channel_config.secondary_species);
                const std::size_t ion =
                    species_id(channel_config.ion_species);
                if (secondary >= species_.size() ||
                    ion >= species_.size()) {
                    throw std::invalid_argument(
                        "ionization product species was not found");
                }
                const auto& secondary_config =
                    species_configs_[secondary];
                const auto& ion_config = species_configs_[ion];
                if (target.charge == 0.0 ||
                    target.weight != secondary_config.weight ||
                    target.weight != ion_config.weight ||
                    secondary_config.mass != target.mass ||
                    secondary_config.charge != target.charge ||
                    ion_config.charge != -target.charge) {
                    throw std::invalid_argument(
                        "ionization currently requires a charged target, "
                        "equal macro weights, secondary mass and charge "
                        "equal to the target, and opposite ion charge");
                }
                ionization_channels_[channel] =
                    IonizationChannelRuntime{secondary, ion};
            } else if (channel_config.process ==
                       CollisionProcessKind::Attachment) {
                const std::size_t product =
                    species_id(channel_config.attachment_species);
                if (product >= species_.size()) {
                    throw std::invalid_argument(
                        "attachment product species was not found");
                }
                const auto& product_config =
                    species_configs_[product];
                if (product == mcc_species_id_ ||
                    !(target.charge < 0.0) ||
                    target.weight != product_config.weight ||
                    target.charge != product_config.charge ||
                    !(product_config.mass > target.mass)) {
                    throw std::invalid_argument(
                        "attachment currently requires a negatively charged "
                        "target and a distinct heavier product with target "
                        "charge and equal macro weight");
                }
                attachment_channels_[channel] =
                    AttachmentChannelRuntime{product};
            }
        }
    }
    std::set<std::string> source_names;
    for (const auto& source_config : config_.sources) {
        if (source_config.name.empty() ||
            !source_names.insert(source_config.name).second) {
            throw std::invalid_argument(
                "unstructured boundary source names must be non-empty and unique");
        }
        if (!species_names.contains(source_config.species)) {
            throw std::invalid_argument(
                "unstructured boundary source species not found: " +
                source_config.species);
        }
        if (!label_set.contains(source_config.boundary)) {
            throw std::invalid_argument(
                "unstructured boundary source label not found: " +
                source_config.boundary);
        }
        if (source_config.particles_per_step == 0) {
            throw std::invalid_argument(
                "unstructured boundary source particles_per_step must be positive");
        }
        if (source_config.end_step != 0 &&
            source_config.end_step <= source_config.start_step) {
            throw std::invalid_argument(
                "unstructured boundary source end_step must exceed start_step");
        }
        if (!std::isfinite(source_config.normal_velocity) ||
            source_config.normal_velocity < 0.0 ||
            !std::isfinite(source_config.tangential_velocity) ||
            !std::isfinite(source_config.thermal_velocity) ||
            !std::isfinite(source_config.out_of_plane_velocity) ||
            source_config.thermal_velocity < 0.0) {
            throw std::invalid_argument(
                "unstructured boundary source velocities are invalid");
        }

        BoundarySourceRuntime source;
        source.config = source_config;
        source.species_id = static_cast<std::size_t>(
            std::find_if(
                species_configs_.begin(), species_configs_.end(),
                [&](const auto& species) {
                    return species.name == source_config.species;
                }) - species_configs_.begin());
        double cumulative_length = 0.0;
        for (std::size_t segment_id = 0;
             segment_id < boundary_segments_.size(); ++segment_id) {
            const auto& segment = boundary_segments_[segment_id];
            if (segment.label != source_config.boundary) continue;
            cumulative_length += std::hypot(
                segment.second.x - segment.first.x,
                segment.second.y - segment.first.y);
            source.segment_indices.push_back(segment_id);
            source.cumulative_lengths.push_back(cumulative_length);
        }
        if (!(cumulative_length > 0.0) || !std::isfinite(cumulative_length)) {
            throw std::runtime_error(
                "unstructured boundary source has no finite boundary length");
        }
        sources_.push_back(std::move(source));
    }
    std::set<std::string> emission_names;
    for (const auto& emission_config : config_.emissions) {
        if (emission_config.name.empty() ||
            !emission_names.insert(emission_config.name).second) {
            throw std::invalid_argument(
                "unstructured emission names must be non-empty and unique");
        }
        if (!label_set.contains(emission_config.boundary) ||
            config_.particle_boundaries.at(emission_config.boundary) !=
                ParticleBoundary::Absorbing) {
            throw std::invalid_argument(
                "unstructured emission boundary must exist and be absorbing");
        }
        if (!species_names.contains(emission_config.incident_species) ||
            !species_names.contains(emission_config.emitted_species)) {
            throw std::invalid_argument(
                "unstructured emission species reference was not found");
        }
        if (!std::isfinite(emission_config.yield) ||
            !(emission_config.yield > 0.0) ||
            emission_config.max_particles_per_impact == 0 ||
            !std::isfinite(emission_config.normal_velocity) ||
            emission_config.normal_velocity < 0.0 ||
            !std::isfinite(emission_config.tangential_velocity) ||
            !std::isfinite(emission_config.thermal_velocity) ||
            !std::isfinite(emission_config.out_of_plane_velocity) ||
            emission_config.thermal_velocity < 0.0) {
            throw std::invalid_argument(
                "unstructured emission parameters are invalid");
        }
        SecondaryEmissionRuntime emission;
        emission.config = emission_config;
        const auto species_id = [&](const std::string& name) {
            return static_cast<std::size_t>(
                std::find_if(
                    species_configs_.begin(), species_configs_.end(),
                    [&](const auto& species) { return species.name == name; }) -
                species_configs_.begin());
        };
        emission.incident_species_id =
            species_id(emission_config.incident_species);
        emission.emitted_species_id =
            species_id(emission_config.emitted_species);
        const double expected_macroparticles =
            emission_config.yield *
            species_configs_[emission.incident_species_id].weight /
            species_configs_[emission.emitted_species_id].weight;
        if (!std::isfinite(expected_macroparticles) ||
            std::ceil(expected_macroparticles) >
                static_cast<double>(
                    emission_config.max_particles_per_impact)) {
            throw std::invalid_argument(
                "unstructured emission macro-particle yield exceeds "
                "max_particles_per_impact");
        }
        emissions_.push_back(std::move(emission));
    }
    poisson_solver_ = std::make_unique<UnstructuredPoissonSolver2D>(
        mesh_, config_.dirichlet_potentials,
        config_.neumann_normal_derivatives, config_.poisson);
}

Vec2 UnstructuredSimulation2D::sample_position(
    const UnstructuredSpecies2DConfig& config,
    std::size_t particle_index,
    std::size_t particle_count,
    std::size_t& profile_attempts) {
    const bool uses_region =
        !config.initialization_region.empty();
    const std::vector<RegionSamplingTriangle>* region_triangles =
        uses_region
            ? &region_sampling_triangles_.at(
                  config.initialization_region)
            : nullptr;
    if (sampling_triangles_.empty() ||
        (region_triangles && region_triangles->empty())) {
        throw std::runtime_error(
            "unstructured initialization region has no sampling area");
    }
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    const double total_area =
        region_triangles
            ? region_triangles->back().cumulative_area
            : sampling_triangles_.back().cumulative_area;
    const bool profiled =
        config.initialization.density_profile !=
        DensityProfileKind::Uniform;
    const auto region_bounds =
        uses_region
            ? std::optional{region_sampling_bounds_.at(
                  config.initialization_region)}
            : std::nullopt;
    const Vec2 profile_minimum =
        config.initialization_minimum.value_or(
            region_bounds
                ? region_bounds->first
                : mesh_.topology().min_corner());
    const Vec2 profile_maximum =
        config.initialization_maximum.value_or(
            region_bounds
                ? region_bounds->second
                : mesh_.topology().max_corner());
    const std::size_t per_particle_limit = profiled
        ? config.initialization.max_profile_sampling_attempts
        : std::size_t{100000};
    for (std::size_t attempt = 0;
         attempt < per_particle_limit; ++attempt) {
        const bool quiet =
            config.initialization.loading == ParticleLoading::QuietStart;
        std::size_t sequence = particle_index;
        if (profiled) {
            if (profile_attempts >=
                config.initialization
                    .max_profile_sampling_attempts) {
                throw std::runtime_error(
                    "unstructured species density-profile sampling exceeded max_profile_sampling_attempts");
            }
            sequence = profile_attempts++;
        }
        const double target =
            (quiet
                 ? (profiled
                        ? quiet_sequence_coordinate(sequence, 0)
                        : quiet_unit_coordinate(
                              particle_index,
                              particle_count, 0))
                 : unit(rng_)) *
            total_area;
        const SamplingTriangle* selected = nullptr;
        if (region_triangles) {
            const auto triangle = std::lower_bound(
                region_triangles->begin(),
                region_triangles->end(), target,
                [](const RegionSamplingTriangle& candidate,
                   double value) {
                    return candidate.cumulative_area < value;
                });
            const auto entry =
                triangle == region_triangles->end()
                    ? std::prev(region_triangles->end())
                    : triangle;
            selected = &sampling_triangles_[
                entry->sampling_triangle_index];
        } else {
            const auto triangle = std::lower_bound(
                sampling_triangles_.begin(),
                sampling_triangles_.end(), target,
                [](const SamplingTriangle& candidate,
                   double value) {
                    return candidate.cumulative_area < value;
                });
            selected =
                triangle == sampling_triangles_.end()
                    ? &sampling_triangles_.back()
                    : &*triangle;
        }
        const double root = std::sqrt(
            quiet
                ? (profiled
                       ? quiet_sequence_coordinate(sequence, 1)
                       : quiet_unit_coordinate(
                             particle_index,
                             particle_count, 1))
                : unit(rng_));
        const double second =
            quiet
                ? (profiled
                       ? quiet_sequence_coordinate(sequence, 2)
                       : quiet_unit_coordinate(
                             particle_index,
                             particle_count, 2))
                : unit(rng_);
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
        if (profiled) {
            const double threshold = quiet
                ? quiet_sequence_coordinate(sequence, 3)
                : unit(rng_);
            if (threshold > density_profile_acceptance(
                    config.initialization,
                    {point.x, point.y, 0.0},
                    {profile_minimum.x, profile_minimum.y, 0.0},
                    {profile_maximum.x, profile_maximum.y, 1.0})) {
                continue;
            }
        }
        return point;
    }
    throw std::runtime_error(
        profiled
            ? "could not sample the requested unstructured density profile within max_profile_sampling_attempts"
            : "could not sample the requested unstructured species initialization region");
}

void UnstructuredSimulation2D::inject_boundary_sources() {
    if (sources_.empty()) return;
    std::vector<std::vector<std::size_t>> reusable(species_.size());
    for (std::size_t species_id = 0; species_id < species_.size(); ++species_id) {
        const auto& particles = species_[species_id].particles();
        for (std::size_t particle_id = 0;
             particle_id < particles.size(); ++particle_id) {
            if (!particles[particle_id].alive) {
                reusable[species_id].push_back(particle_id);
            }
        }
    }

    const Vec2 minimum = mesh_.topology().min_corner();
    const Vec2 maximum = mesh_.topology().max_corner();
    const double domain_scale =
        std::max(maximum.x - minimum.x, maximum.y - minimum.y);
    const double inset =
        std::max(1e-12 * domain_scale,
                 1024.0 * std::numeric_limits<double>::epsilon() * domain_scale);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::normal_distribution<double> thermal(0.0, 1.0);

    for (auto& source : sources_) {
        if (step_ < source.config.start_step ||
            (source.config.end_step != 0 &&
             step_ >= source.config.end_step)) {
            continue;
        }
        auto& species = species_[source.species_id];
        auto& particles = species.particles();
        auto& locations = particle_locations_[source.species_id];
        for (std::size_t injected = 0;
             injected < source.config.particles_per_step; ++injected) {
            std::size_t particle_id = 0;
            if (!reusable[source.species_id].empty()) {
                particle_id = reusable[source.species_id].back();
                reusable[source.species_id].pop_back();
                particles[particle_id] = Particle2D{};
                locations[particle_id] = UnstructuredParticleLocation2D{};
            } else {
                if (particles.size() >= config_.max_particles_per_species) {
                    throw std::runtime_error(
                        "unstructured boundary source exceeded "
                        "max_particles_per_species");
                }
                particle_id = particles.size();
                particles.emplace_back();
                locations.emplace_back();
            }

            const double target =
                unit(rng_) * source.cumulative_lengths.back();
            const auto selected = std::lower_bound(
                source.cumulative_lengths.begin(),
                source.cumulative_lengths.end(), target);
            const std::size_t local_segment =
                selected == source.cumulative_lengths.end()
                    ? source.cumulative_lengths.size() - 1
                    : static_cast<std::size_t>(
                          selected - source.cumulative_lengths.begin());
            const auto& segment =
                boundary_segments_[source.segment_indices[local_segment]];
            const Vec2 edge = subtract(segment.second, segment.first);
            const Vec2 tangent{
                segment.inward_normal.y, -segment.inward_normal.x};
            constexpr double endpoint_margin = 1e-10;
            const double along =
                endpoint_margin +
                (1.0 - 2.0 * endpoint_margin) * unit(rng_);

            auto& particle = particles[particle_id];
            particle.position = add(
                add(segment.first, scale(edge, along)),
                scale(segment.inward_normal, inset));
            const double normal_speed =
                source.config.normal_velocity +
                source.config.thermal_velocity * std::abs(thermal(rng_));
            const double tangent_speed =
                source.config.tangential_velocity +
                source.config.thermal_velocity * thermal(rng_);
            particle.velocity = add(
                scale(segment.inward_normal, normal_speed),
                scale(tangent, tangent_speed));
            particle.velocity_z =
                source.config.out_of_plane_velocity +
                source.config.thermal_velocity * thermal(rng_);
            particle.velocity_half = particle.velocity;
            particle.velocity_half_z = particle.velocity_z;
            particle.alive = true;

            bool cache_hit = false;
            const auto electric = interpolate_electric(
                mesh_, particle.position, locations[particle_id], &cache_hit);
            cache_hit ? ++timing_.location_cache_hits
                      : ++timing_.location_searches;
            if (!electric) {
                throw std::runtime_error(
                    "unstructured boundary source generated an exterior particle");
            }
            initialize_particle_pusher(
                particle, *electric, species.charge() / species.mass(),
                magnetic_field(config_, particle.position),
                config_.dt);
            if (source.injected_particles ==
                std::numeric_limits<std::size_t>::max()) {
                throw std::overflow_error(
                    "unstructured source diagnostic counter overflow");
            }
            ++source.injected_particles;
        }
    }
}

void UnstructuredSimulation2D::initialize() {
    time_ = 0.0;
    step_ = 0;
    timing_ = {};
    for (auto& source : sources_) source.injected_particles = 0;
    for (auto& emission : emissions_) emission.emitted_particles = 0;
    for (auto& [label, count] : absorbed_by_label_) count = 0;
    for (auto& [species, boundaries] : impact_flux_) {
        (void)species;
        for (auto& [label, flux] : boundaries) {
            (void)label;
            flux = {};
        }
    }
    const bool external_state =
        !config_.initial_state_path.empty();
    for (std::size_t species_id = 0;
         species_id < species_.size(); ++species_id) {
        const auto particle_count =
            species_configs_[species_id].particles;
        species_[species_id].particles().assign(
            particle_count, Particle2D{});
        particle_locations_[species_id].assign(
            particle_count,
            UnstructuredParticleLocation2D{});
    }
    if (external_state) {
        std::vector<ExternalSpeciesExpectation> expected;
        expected.reserve(species_.size());
        for (const auto& species : species_) {
            expected.push_back({
                species.name(),
                species.config().particles});
        }
        initial_state_metadata_ =
            load_validated_external_particle_state_bounded(
                config_.initial_state_path, 2, 3,
                config_.units.system, expected,
                "unstructured 2D simulation",
                [&](std::size_t species_id,
                    std::size_t particle_index,
                    const ExternalParticleRecord& record) {
                    auto& particles =
                        species_[species_id].particles();
                    auto& locations =
                        particle_locations_[species_id];
                    const auto& species_config =
                        species_configs_[species_id];
                    const Vec2 position{
                        record.position.x,
                        record.position.y};
                    const auto location =
                        mesh_.locate_point(position);
                    if (!location) {
                        throw std::runtime_error(
                            "external particle for species '" +
                            species_config.name +
                            "' lies outside the imported mesh");
                    }
                    if (species_config.initialization_minimum) {
                        const Vec2 minimum =
                            *species_config.initialization_minimum;
                        const Vec2 maximum =
                            *species_config.initialization_maximum;
                        if (position.x < minimum.x ||
                            position.x > maximum.x ||
                            position.y < minimum.y ||
                            position.y > maximum.y) {
                            throw std::runtime_error(
                                "external particle for species '" +
                                species_config.name +
                                "' lies outside its configured initialization bounds");
                        }
                    }
                    if (!species_config.initialization_region.empty() &&
                        mesh_.topology()
                                .cell_by_id(location->cell_id)
                                .label !=
                            species_config.initialization_region) {
                        throw std::runtime_error(
                            "external particle for species '" +
                            species_config.name +
                            "' lies outside its configured initialization region");
                    }
                    auto& particle =
                        particles.at(particle_index);
                    particle.position = position;
                    particle.velocity = {
                        record.velocity.x,
                        record.velocity.y};
                    particle.velocity_z =
                        record.velocity.z;
                    particle.velocity_half =
                        particle.velocity;
                    particle.velocity_half_z =
                        particle.velocity_z;
                    particle.alive = true;
                    locations.at(particle_index) = {
                        *location, true};
                },
                config_.initial_state_signature);
    }
    for (std::size_t species_id = 0; species_id < species_.size(); ++species_id) {
        auto& particles = species_[species_id].particles();
        const auto& species_config = species_configs_[species_id];
        if (external_state) {
            continue;
        }
        std::size_t profile_attempts = 0;
        const double thermal_x = resolved_thermal_velocity(
            species_config.initialization, 0,
            species_config.thermal_velocity);
        const double thermal_y = resolved_thermal_velocity(
            species_config.initialization, 1,
            species_config.thermal_velocity);
        const double thermal_z = resolved_thermal_velocity(
            species_config.initialization, 2,
            species_config.thermal_velocity);
        if (species_config.initialization.loading ==
            ParticleLoading::Random) {
            std::normal_distribution<double> velocity_x(
                species_config.drift_velocity_x, thermal_x);
            std::normal_distribution<double> velocity_y(
                species_config.drift_velocity_y, thermal_y);
            std::normal_distribution<double> velocity_z(
                species_config.drift_velocity_z, thermal_z);
            for (std::size_t particle_index = 0;
                 particle_index < particles.size(); ++particle_index) {
                auto& particle = particles[particle_index];
                particle.position = sample_position(
                    species_config, particle_index,
                    particles.size(), profile_attempts);
                particle.velocity = {
                    velocity_x(rng_), velocity_y(rng_)};
                particle.velocity_z = velocity_z(rng_);
                particle.velocity_half = particle.velocity;
                particle.velocity_half_z = particle.velocity_z;
                particle.alive = true;
            }
            continue;
        }

        const auto velocity_x = initialize_velocity_component(
            particles.size(), species_config.drift_velocity_x,
            thermal_x, species_config.initialization.loading, rng_);
        const auto velocity_y = initialize_velocity_component(
            particles.size(), species_config.drift_velocity_y,
            thermal_y, species_config.initialization.loading, rng_);
        const auto velocity_z = initialize_velocity_component(
            particles.size(), species_config.drift_velocity_z,
            thermal_z, species_config.initialization.loading, rng_);
        for (std::size_t particle_index = 0;
             particle_index < particles.size(); ++particle_index) {
            auto& particle = particles[particle_index];
            particle.position = sample_position(
                species_config, particle_index,
                particles.size(), profile_attempts);
            particle.velocity = {
                velocity_x[particle_index],
                velocity_y[particle_index]};
            particle.velocity_z = velocity_z[particle_index];
            particle.velocity_half = particle.velocity;
            particle.velocity_half_z = particle.velocity_z;
            particle.alive = true;
        }
    }
    deposit_and_solve();
    const auto particle_start = SteadyClock::now();
    for (std::size_t species_id = 0; species_id < species_.size(); ++species_id) {
        auto& species = species_[species_id];
        const double charge_to_mass = species.charge() / species.mass();
        const std::size_t workers =
            particle_worker_count(species.particles().size(), config_.runtime);
        std::vector<std::size_t> local_hits(workers, 0);
        std::vector<std::size_t> local_searches(workers, 0);
        parallel_particle_chunks(
            species.particles(), config_.runtime,
            [&](std::size_t worker, std::size_t particle_id,
                Particle2D& particle) {
                bool cache_hit = false;
                const auto electric =
                    interpolate_electric(
                        mesh_, particle.position,
                        particle_locations_[species_id][particle_id],
                        &cache_hit);
                cache_hit ? ++local_hits[worker] : ++local_searches[worker];
                if (!electric) {
                    throw std::runtime_error(
                        "initialized particle is outside imported mesh");
                }
                initialize_particle_pusher(
                    particle, *electric, charge_to_mass,
                    magnetic_field(config_, particle.position),
                    config_.dt);
            });
        for (std::size_t worker = 0; worker < workers; ++worker) {
            timing_.location_cache_hits += local_hits[worker];
            timing_.location_searches += local_searches[worker];
        }
    }
    timing_.particle_seconds +=
        std::chrono::duration<double>(SteadyClock::now() - particle_start).count();
    initialized_ = true;
}

void UnstructuredSimulation2D::deposit_and_solve() {
    const auto deposition_start = SteadyClock::now();
    mesh_.clear_charge();
    for (std::size_t species_id = 0; species_id < species_.size(); ++species_id) {
        const auto& species = species_[species_id];
        const auto deposit =
            deposit_charge_shape(
                mesh_, species.particles(), species.charge(),
                species.weight(), config_.runtime,
                particle_locations_[species_id]);
        timing_.location_cache_hits += deposit.location_cache_hits;
        timing_.location_searches += deposit.location_searches;
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

std::optional<UnstructuredSimulation2D::BoundaryImpact>
UnstructuredSimulation2D::advance_with_boundaries(
    Particle2D& particle, Vec2 previous_position) {
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
                return std::nullopt;
            }
            throw std::runtime_error(
                "particle left imported domain without a resolvable boundary intersection");
        }

        const ParticleBoundary policy = config_.particle_boundaries.at(hit_segment->label);
        if (policy == ParticleBoundary::Absorbing) {
            particle.position = add(start, scale(displacement, earliest));
            particle.alive = false;
            return BoundaryImpact{
                0,
                0,
                static_cast<std::size_t>(
                    hit_segment - boundary_segments_.data()),
                particle.position,
                particle.velocity_half,
                particle.velocity_half_z,
            };
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

void UnstructuredSimulation2D::process_boundary_impacts(
    std::vector<BoundaryImpact> impacts) {
    for (auto& [species, boundaries] : impact_flux_) {
        (void)species;
        for (auto& [label, flux] : boundaries) {
            (void)label;
            flux.last_step_macroparticles = 0;
            flux.last_step_physical_particles = 0.0;
            flux.physical_particle_rate = 0.0;
            flux.physical_particle_flux = 0.0;
        }
    }
    if (impacts.empty()) return;
    std::sort(
        impacts.begin(), impacts.end(),
        [](const BoundaryImpact& first, const BoundaryImpact& second) {
            return std::tie(first.species_id, first.particle_id) <
                   std::tie(second.species_id, second.particle_id);
        });

    std::vector<std::vector<std::size_t>> reusable(species_.size());
    for (std::size_t species_id = 0; species_id < species_.size(); ++species_id) {
        const auto& particles = species_[species_id].particles();
        for (std::size_t particle_id = 0;
             particle_id < particles.size(); ++particle_id) {
            if (!particles[particle_id].alive) {
                reusable[species_id].push_back(particle_id);
            }
        }
    }
    const Vec2 minimum = mesh_.topology().min_corner();
    const Vec2 maximum = mesh_.topology().max_corner();
    const double domain_scale =
        std::max(maximum.x - minimum.x, maximum.y - minimum.y);
    const double inset =
        std::max(1e-12 * domain_scale,
                 1024.0 * std::numeric_limits<double>::epsilon() *
                     domain_scale);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::normal_distribution<double> thermal(0.0, 1.0);

    for (const auto& impact : impacts) {
        const auto& segment = boundary_segments_[impact.segment_id];
        const auto& incident = species_[impact.species_id];
        if (absorbed_by_label_.at(segment.label) ==
                std::numeric_limits<std::size_t>::max() ||
            impact_flux_
                    .at(species_configs_[impact.species_id].name)
                    .at(segment.label)
                    .macroparticles ==
                std::numeric_limits<std::size_t>::max()) {
            throw std::overflow_error(
                "unstructured boundary diagnostic counter overflow");
        }
        ++absorbed_by_label_.at(segment.label);
        auto& flux = impact_flux_
                         .at(species_configs_[impact.species_id].name)
                         .at(segment.label);
        const double physical_particles = incident.weight();
        const double speed_squared =
            dot(impact.incident_velocity, impact.incident_velocity) +
            impact.incident_velocity_z * impact.incident_velocity_z;
        ++flux.macroparticles;
        flux.physical_particles += physical_particles;
        flux.charge += incident.charge() * physical_particles;
        flux.kinetic_energy +=
            0.5 * incident.mass() * physical_particles * speed_squared;
        if (!std::isfinite(flux.physical_particles) ||
            !std::isfinite(flux.charge) ||
            !std::isfinite(flux.kinetic_energy)) {
            throw std::overflow_error(
                "unstructured boundary physical diagnostic overflow");
        }
        ++flux.last_step_macroparticles;
        flux.last_step_physical_particles += physical_particles;
        flux.physical_particle_rate =
            flux.last_step_physical_particles / config_.dt;
        flux.physical_particle_flux =
            flux.physical_particle_rate /
            boundary_lengths_.at(segment.label);
        if (!std::isfinite(flux.physical_particle_rate) ||
            !std::isfinite(flux.physical_particle_flux)) {
            throw std::overflow_error(
                "unstructured boundary rate diagnostic overflow");
        }

        for (auto& emission : emissions_) {
            if (emission.incident_species_id != impact.species_id ||
                emission.config.boundary != segment.label) {
                continue;
            }
            const auto& emitted_config =
                species_configs_[emission.emitted_species_id];
            const double expected =
                emission.config.yield * incident.weight() /
                emitted_config.weight;
            const double integral = std::floor(expected);
            std::size_t count = static_cast<std::size_t>(integral);
            if (unit(rng_) < expected - integral) ++count;
            if (count > emission.config.max_particles_per_impact) {
                throw std::runtime_error(
                    "unstructured emission exceeded max_particles_per_impact");
            }

            auto& emitted_species = species_[emission.emitted_species_id];
            auto& particles = emitted_species.particles();
            auto& locations =
                particle_locations_[emission.emitted_species_id];
            const Vec2 tangent{
                segment.inward_normal.y, -segment.inward_normal.x};
            for (std::size_t emitted = 0; emitted < count; ++emitted) {
                std::size_t particle_id = 0;
                if (!reusable[emission.emitted_species_id].empty()) {
                    particle_id =
                        reusable[emission.emitted_species_id].back();
                    reusable[emission.emitted_species_id].pop_back();
                    particles[particle_id] = Particle2D{};
                    locations[particle_id] =
                        UnstructuredParticleLocation2D{};
                } else {
                    if (particles.size() >=
                        config_.max_particles_per_species) {
                        throw std::runtime_error(
                            "unstructured emission exceeded "
                            "max_particles_per_species");
                    }
                    particle_id = particles.size();
                    particles.emplace_back();
                    locations.emplace_back();
                }
                auto& particle = particles[particle_id];
                const Vec2 edge =
                    subtract(segment.second, segment.first);
                const double edge_length_squared = dot(edge, edge);
                const double edge_length = std::sqrt(edge_length_squared);
                const double endpoint_margin =
                    std::min(0.25, 4.0 * inset / edge_length);
                const double along = std::clamp(
                    dot(subtract(impact.position, segment.first), edge) /
                        edge_length_squared,
                    endpoint_margin, 1.0 - endpoint_margin);
                particle.position = add(
                    add(segment.first, scale(edge, along)),
                    scale(segment.inward_normal, inset));
                const double normal_speed =
                    emission.config.normal_velocity +
                    emission.config.thermal_velocity *
                        std::abs(thermal(rng_));
                const double tangential_speed =
                    emission.config.tangential_velocity +
                    emission.config.thermal_velocity * thermal(rng_);
                particle.velocity = add(
                    scale(segment.inward_normal, normal_speed),
                    scale(tangent, tangential_speed));
                particle.velocity_z =
                    emission.config.out_of_plane_velocity +
                    emission.config.thermal_velocity * thermal(rng_);
                particle.velocity_half = particle.velocity;
                particle.velocity_half_z = particle.velocity_z;
                particle.alive = true;
                bool cache_hit = false;
                const auto electric = interpolate_electric(
                    mesh_, particle.position, locations[particle_id],
                    &cache_hit);
                cache_hit ? ++timing_.location_cache_hits
                          : ++timing_.location_searches;
                if (!electric) {
                    throw std::runtime_error(
                        "unstructured emission generated an exterior particle");
                }
                initialize_particle_pusher(
                    particle, *electric,
                    emitted_species.charge() / emitted_species.mass(),
                    magnetic_field(config_, particle.position),
                    config_.dt);
                if (emission.emitted_particles ==
                    std::numeric_limits<std::size_t>::max()) {
                    throw std::overflow_error(
                        "unstructured emission diagnostic counter overflow");
                }
                ++emission.emitted_particles;
            }
        }
    }
}

void UnstructuredSimulation2D::step() {
    if (!initialized_) initialize();
    const auto first_particle_start = SteadyClock::now();
    inject_boundary_sources();
    std::vector<BoundaryImpact> impacts;
    for (std::size_t species_id = 0; species_id < species_.size(); ++species_id) {
        auto& species = species_[species_id];
        const double charge_to_mass = species.charge() / species.mass();
        const std::size_t workers =
            particle_worker_count(species.particles().size(), config_.runtime);
        std::vector<std::vector<BoundaryImpact>> local_impacts(workers);
        std::vector<std::size_t> local_hits(workers, 0);
        std::vector<std::size_t> local_searches(workers, 0);
        parallel_particle_chunks(
            species.particles(), config_.runtime,
            [&](std::size_t worker, std::size_t particle_id,
                Particle2D& particle) {
                if (!particle.alive) return;
                bool cache_hit = false;
                const auto electric =
                    interpolate_electric(
                        mesh_, particle.position,
                        particle_locations_[species_id][particle_id],
                        &cache_hit);
                cache_hit ? ++local_hits[worker] : ++local_searches[worker];
                if (!electric) {
                    throw std::runtime_error(
                        "live particle is outside imported mesh before push");
                }
                kick_particle(
                    particle, *electric, charge_to_mass,
                    magnetic_field(config_, particle.position),
                    config_.dt);
                const Vec2 previous = particle.position;
                drift_leapfrog(particle, config_.dt);
                auto impact = advance_with_boundaries(particle, previous);
                if (impact) {
                    impact->species_id = species_id;
                    impact->particle_id = particle_id;
                    local_impacts[worker].push_back(*impact);
                }
            });
        for (auto& local : local_impacts) {
            impacts.insert(
                impacts.end(),
                std::make_move_iterator(local.begin()),
                std::make_move_iterator(local.end()));
        }
        for (std::size_t worker = 0; worker < workers; ++worker) {
            timing_.location_cache_hits += local_hits[worker];
            timing_.location_searches += local_searches[worker];
        }
    }
    process_boundary_impacts(std::move(impacts));
    timing_.particle_seconds +=
        std::chrono::duration<double>(SteadyClock::now() - first_particle_start).count();

    deposit_and_solve();
    const auto second_particle_start = SteadyClock::now();
    for (std::size_t species_id = 0; species_id < species_.size(); ++species_id) {
        auto& species = species_[species_id];
        const double charge_to_mass = species.charge() / species.mass();
        const std::size_t workers =
            particle_worker_count(species.particles().size(), config_.runtime);
        std::vector<std::size_t> local_hits(workers, 0);
        std::vector<std::size_t> local_searches(workers, 0);
        parallel_particle_chunks(
            species.particles(), config_.runtime,
            [&](std::size_t worker, std::size_t particle_id,
                Particle2D& particle) {
                if (!particle.alive) return;
                bool cache_hit = false;
                const auto electric =
                    interpolate_electric(
                        mesh_, particle.position,
                        particle_locations_[species_id][particle_id],
                        &cache_hit);
                cache_hit ? ++local_hits[worker] : ++local_searches[worker];
                if (!electric) {
                    throw std::runtime_error(
                        "live particle is outside imported mesh after push");
                }
                synchronize_particle(
                    particle, *electric, charge_to_mass,
                    magnetic_field(config_, particle.position),
                    config_.dt);
            });
        for (std::size_t worker = 0; worker < workers; ++worker) {
            timing_.location_cache_hits += local_hits[worker];
            timing_.location_searches += local_searches[worker];
        }
    }
    timing_.particle_seconds +=
        std::chrono::duration<double>(SteadyClock::now() - second_particle_start).count();
    const auto collision_start = SteadyClock::now();
    apply_collisions();
    timing_.particle_seconds +=
        std::chrono::duration<double>(
            SteadyClock::now() - collision_start).count();
    ++step_;
    time_ += config_.dt;
}

void UnstructuredSimulation2D::apply_collisions() {
    if (!mcc_model_) return;
    auto& species = species_[mcc_species_id_];
    const double charge_to_mass =
        species.charge() / species.mass();
    struct IonizationProduct {
        Vec2 position{};
        Vec3 secondary_velocity{};
        Vec3 ion_velocity{};
        IonizationChannelRuntime channel{};
    };
    struct AttachmentProduct {
        Vec2 position{};
        std::size_t consumed_particle_id{0};
        Vec3 product_velocity{};
        AttachmentChannelRuntime channel{};
    };
    std::vector<IonizationProduct> ionization_products;
    std::vector<AttachmentProduct> attachment_products;
    const std::size_t initial_particle_count =
        species.particles().size();
    for (std::size_t particle_id = 0;
         particle_id < initial_particle_count; ++particle_id) {
        auto& particle = species.particles()[particle_id];
        if (!particle.alive) continue;
        Vec3 velocity{
            particle.velocity.x,
            particle.velocity.y,
            particle.velocity_z};
        const auto statistics =
            mcc_model_->collide(velocity, config_.dt, rng_);
        add_collision_statistics(collision_totals_, statistics);
        add_collision_statistics(collision_interval_, statistics);
        for (const auto& secondary : statistics.secondaries) {
            if (secondary.channel >= ionization_channels_.size() ||
                !ionization_channels_[secondary.channel]) {
                throw std::logic_error(
                    "MCC produced an unmapped ionization channel");
            }
            ionization_products.push_back({
                particle.position,
                secondary.velocity,
                secondary.ion_velocity,
                *ionization_channels_[secondary.channel]});
        }
        if (statistics.primary_removal_channel) {
            const std::size_t channel =
                *statistics.primary_removal_channel;
            if (channel >= attachment_channels_.size() ||
                !attachment_channels_[channel]) {
                throw std::logic_error(
                    "MCC produced an unmapped attachment channel");
            }
            attachment_products.push_back({
                particle.position,
                particle_id,
                statistics.primary_removal_product_velocity.value_or(
                    Vec3{}),
                *attachment_channels_[channel]});
            continue;
        }
        particle.velocity = {velocity.x, velocity.y};
        particle.velocity_z = velocity.z;
        bool cache_hit = false;
        const auto electric = interpolate_electric(
            mesh_, particle.position,
            particle_locations_[mcc_species_id_][particle_id],
            &cache_hit);
        cache_hit ? ++timing_.location_cache_hits
                  : ++timing_.location_searches;
        if (!electric) {
            throw std::runtime_error(
                "MCC particle is outside imported mesh");
        }
        initialize_particle_pusher(
            particle, *electric, charge_to_mass,
            magnetic_field(config_, particle.position),
            config_.dt);
    }
    std::vector<std::size_t> required_products(species_.size(), 0);
    for (const auto& product : ionization_products) {
        ++required_products[product.channel.secondary_species_id];
        ++required_products[product.channel.ion_species_id];
    }
    for (const auto& product : attachment_products) {
        ++required_products[product.channel.product_species_id];
    }
    for (std::size_t species_id = 0;
         species_id < species_.size(); ++species_id) {
        const auto& particles = species_[species_id].particles();
        const std::size_t reusable = static_cast<std::size_t>(
            std::count_if(
                particles.begin(), particles.end(),
                [](const Particle2D& particle) {
                    return !particle.alive;
                }));
        const std::size_t growth =
            required_products[species_id] > reusable
                ? required_products[species_id] - reusable
                : 0;
        if (growth >
            config_.max_particles_per_species - particles.size()) {
            throw std::runtime_error(
                "reactive collisions exceeded max_particles_per_species "
                "for species '" +
                species_[species_id].name() + "'");
        }
    }
    for (const auto& product : attachment_products) {
        species.particles()[product.consumed_particle_id].alive = false;
    }
    const auto append_product = [&](std::size_t species_id,
                                    Vec2 position,
                                    Vec3 velocity) {
        auto& product_species = species_[species_id];
        auto& particles = product_species.particles();
        auto& locations = particle_locations_[species_id];
        auto dead = std::find_if(
            particles.begin(), particles.end(),
            [](const Particle2D& particle) {
                return !particle.alive;
            });
        std::size_t particle_id = 0;
        if (dead == particles.end()) {
            if (particles.size() >=
                config_.max_particles_per_species) {
                throw std::runtime_error(
                    "reactive collisions exceeded "
                    "max_particles_per_species");
            }
            particle_id = particles.size();
            particles.emplace_back();
            locations.emplace_back();
        } else {
            particle_id = static_cast<std::size_t>(
                dead - particles.begin());
        }
        auto& particle = particles[particle_id];
        particle = {};
        particle.position = position;
        particle.velocity = {velocity.x, velocity.y};
        particle.velocity_z = velocity.z;
        particle.velocity_half = particle.velocity;
        particle.velocity_half_z = particle.velocity_z;
        particle.alive = true;
        bool cache_hit = false;
        const auto electric = interpolate_electric(
            mesh_, position, locations[particle_id], &cache_hit);
        cache_hit ? ++timing_.location_cache_hits
                  : ++timing_.location_searches;
        if (!electric) {
            throw std::runtime_error(
                "collision product is outside imported mesh");
        }
        initialize_particle_pusher(
            particle, *electric,
            product_species.charge() / product_species.mass(),
            magnetic_field(config_, particle.position),
            config_.dt);
    };
    for (const auto& product : ionization_products) {
        append_product(
            product.channel.secondary_species_id,
            product.position, product.secondary_velocity);
        append_product(
            product.channel.ion_species_id,
            product.position, product.ion_velocity);
    }
    for (const auto& product : attachment_products) {
        append_product(
            product.channel.product_species_id,
            product.position, product.product_velocity);
    }
}

UnstructuredDiagnosticSample2D UnstructuredSimulation2D::sample() const {
    UnstructuredDiagnosticSample2D result;
    result.step = step_;
    result.time = time_;
    result.absorbed_by_label = absorbed_by_label_;
    for (const auto& source : sources_) {
        result.injected_by_source.emplace(
            source.config.name, source.injected_particles);
    }
    for (const auto& emission : emissions_) {
        result.emitted_by_rule.emplace(
            emission.config.name, emission.emitted_particles);
    }
    result.impact_flux = impact_flux_;
    result.poisson = last_poisson_;
    for (const auto& species : species_) {
        result.kinetic_energy += species.kinetic_energy();
        result.live_particles += species.live_count();
    }
    for (std::size_t i = 0; i < mesh_.size(); ++i) {
        const double area = mesh_.node_control_areas()[i];
        const Vec2 electric = mesh_.electric()[i];
        result.field_energy +=
            0.5 * config_.units.permittivity() *
            (electric.x * electric.x + electric.y * electric.y) * area;
        result.charge_l1 += std::abs(mesh_.rho()[i]) * area;
    }
    result.total_energy = result.kinetic_energy + result.field_energy;
    return result;
}

void UnstructuredSimulation2D::write_diagnostics_header(std::ofstream& output) const {
    output << "step,time,kinetic_energy,field_energy,total_energy,charge_l1,live_particles"
              ",poisson_iterations,poisson_initial_residual,poisson_final_residual"
              ",particle_seconds,deposition_seconds,field_solve_seconds"
              ",location_cache_hits,location_searches";
    for (const auto& [label, count] : absorbed_by_label_) {
        (void)count;
        output << ',' << csv_quote("absorbed_" + label);
    }
    std::map<std::string, std::size_t> injected;
    for (const auto& source : sources_) {
        injected.emplace(source.config.name, source.injected_particles);
    }
    for (const auto& [name, count] : injected) {
        (void)count;
        output << ',' << csv_quote("injected_" + name);
    }
    std::map<std::string, std::size_t> emitted;
    for (const auto& emission : emissions_) {
        emitted.emplace(
            emission.config.name, emission.emitted_particles);
    }
    for (const auto& [name, count] : emitted) {
        (void)count;
        output << ',' << csv_quote("emitted_" + name);
    }
    for (const auto& [species, boundaries] : impact_flux_) {
        for (const auto& [label, flux] : boundaries) {
            (void)flux;
            const std::string key = species + "@" + label;
            output << ',' << csv_quote("impact_macroparticles_" + key)
                   << ',' << csv_quote("impact_physical_particles_" + key)
                   << ',' << csv_quote("impact_charge_" + key)
                   << ',' << csv_quote("impact_kinetic_energy_" + key)
                   << ',' << csv_quote("impact_rate_" + key)
                   << ',' << csv_quote("impact_flux_" + key);
        }
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
           << timing_.field_solve_seconds << ','
           << timing_.location_cache_hits << ',' << timing_.location_searches;
    for (const auto& [label, count] : value.absorbed_by_label) {
        (void)label;
        output << ',' << count;
    }
    for (const auto& [name, count] : value.injected_by_source) {
        (void)name;
        output << ',' << count;
    }
    for (const auto& [name, count] : value.emitted_by_rule) {
        (void)name;
        output << ',' << count;
    }
    for (const auto& [species, boundaries] : value.impact_flux) {
        (void)species;
        for (const auto& [label, flux] : boundaries) {
            (void)label;
            output << ',' << flux.macroparticles
                   << ',' << flux.physical_particles
                   << ',' << flux.charge
                   << ',' << flux.kinetic_energy
                   << ',' << flux.physical_particle_rate
                   << ',' << flux.physical_particle_flux;
        }
    }
    output << '\n';
    output.flush();
}

void UnstructuredSimulation2D::write_particle_sample(std::size_t step) const {
    std::ofstream output(config_.output_dir / ("particles_" + std::to_string(step) + ".csv"));
    if (!output) throw std::runtime_error("cannot open unstructured particle output");
    output << "species_id,species,x,y,vx,vy,vz,alive\n";
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
                   << particle.velocity_z << ','
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
    output << CHECKPOINT_MAGIC_V6 << '\n';
    output << "mesh_signature " << mesh_signature() << '\n';
    output << "units " << to_string(config_.units.system) << ' '
           << config_.units.relative_permittivity << ' '
           << config_.units.permittivity() << '\n';
    output << "step " << step_ << '\n';
    output << "time " << time_ << '\n';
    output << "collision_model "
           << (mcc_model_ ? "null_collision" : "off") << ' '
           << (mcc_model_ ? mcc_model_->signature() : 0) << '\n';
    output << "collision_totals " << collision_totals_.candidates
           << ' ' << collision_totals_.null_collisions
           << ' ' << collision_totals_.channel_collisions.size();
    for (const auto count : collision_totals_.channel_collisions) {
        output << ' ' << count;
    }
    output << '\n';
    output << "absorbed_count " << absorbed_by_label_.size() << '\n';
    for (const auto& [label, count] : absorbed_by_label_) {
        output << "absorbed " << std::quoted(label) << ' ' << count << '\n';
    }
    output << "source_count " << sources_.size() << '\n';
    for (const auto& source : sources_) {
        output << "source " << std::quoted(source.config.name) << ' '
               << std::quoted(source.config.species) << ' '
               << std::quoted(source.config.boundary) << ' '
               << source.config.particles_per_step << ' '
               << source.config.start_step << ' '
               << source.config.end_step << ' '
               << source.config.normal_velocity << ' '
               << source.config.tangential_velocity << ' '
               << source.config.thermal_velocity << ' '
               << source.config.out_of_plane_velocity << ' '
               << source.injected_particles << '\n';
    }
    output << "emission_count " << emissions_.size() << '\n';
    for (const auto& emission : emissions_) {
        output << "emission " << std::quoted(emission.config.name) << ' '
               << std::quoted(emission.config.boundary) << ' '
               << std::quoted(emission.config.incident_species) << ' '
               << std::quoted(emission.config.emitted_species) << ' '
               << emission.config.yield << ' '
               << emission.config.max_particles_per_impact << ' '
               << emission.config.normal_velocity << ' '
               << emission.config.tangential_velocity << ' '
               << emission.config.thermal_velocity << ' '
               << emission.config.out_of_plane_velocity << ' '
               << emission.emitted_particles << '\n';
    }
    std::size_t impact_flux_count = 0;
    for (const auto& [species, boundaries] : impact_flux_) {
        (void)species;
        impact_flux_count += boundaries.size();
    }
    output << "impact_flux_count " << impact_flux_count << '\n';
    for (const auto& [species, boundaries] : impact_flux_) {
        for (const auto& [boundary, flux] : boundaries) {
            output << "impact_flux " << std::quoted(species) << ' '
                   << std::quoted(boundary) << ' '
                   << flux.macroparticles << ' '
                   << flux.physical_particles << ' '
                   << flux.charge << ' '
                   << flux.kinetic_energy << ' '
                   << flux.last_step_macroparticles << ' '
                   << flux.last_step_physical_particles << ' '
                   << flux.physical_particle_rate << ' '
                   << flux.physical_particle_flux << '\n';
        }
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
                   << particle.velocity_z << ' '
                   << particle.velocity_half.x << ' ' << particle.velocity_half.y << ' '
                   << particle.velocity_half_z << ' '
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
    const bool checkpoint_v1 = magic == CHECKPOINT_MAGIC_V1;
    const bool checkpoint_v2 = magic == CHECKPOINT_MAGIC_V2;
    const bool checkpoint_v3 = magic == CHECKPOINT_MAGIC_V3;
    const bool checkpoint_v4 = magic == CHECKPOINT_MAGIC_V4;
    const bool checkpoint_v5 = magic == CHECKPOINT_MAGIC_V5;
    const bool checkpoint_v6 = magic == CHECKPOINT_MAGIC_V6;
    if (!checkpoint_v1 && !checkpoint_v2 && !checkpoint_v3 &&
        !checkpoint_v4 && !checkpoint_v5 && !checkpoint_v6) {
        throw std::runtime_error("invalid unstructured checkpoint magic");
    }
    std::string key;
    std::uint64_t signature = 0;
    input >> key >> signature;
    if (key != "mesh_signature" || signature != mesh_signature()) {
        throw std::runtime_error("unstructured checkpoint mesh does not match configured geometry");
    }
    input >> key;
    if (key == "units") {
        std::string unit_system;
        double relative_permittivity = 0.0;
        double permittivity = 0.0;
        input >> unit_system >> relative_permittivity >> permittivity;
        if ((!checkpoint_v4 && !checkpoint_v5 && !checkpoint_v6) ||
            unit_system != to_string(config_.units.system) ||
            relative_permittivity !=
                config_.units.relative_permittivity ||
            permittivity != config_.units.permittivity()) {
            throw std::runtime_error(
                "unstructured checkpoint unit system mismatch");
        }
        input >> key;
    } else if (checkpoint_v4 || checkpoint_v5 || checkpoint_v6 ||
               config_.units.system != UnitSystem::Normalized ||
               config_.units.relative_permittivity != 1.0) {
        throw std::runtime_error(
            "legacy unstructured checkpoint without unit metadata "
            "requires normalized units");
    }
    input >> step_;
    if (key != "step") throw std::runtime_error("unstructured checkpoint missing step");
    input >> key >> time_;
    if (key != "time" || !std::isfinite(time_) || time_ < 0.0) {
        throw std::runtime_error("unstructured checkpoint has invalid time");
    }
    if (checkpoint_v6) {
        std::string collision_model;
        std::uint64_t collision_signature = 0;
        input >> key >> collision_model >> collision_signature;
        const std::string configured_model =
            mcc_model_ ? "null_collision" : "off";
        const std::uint64_t configured_signature =
            mcc_model_ ? mcc_model_->signature() : 0;
        if (key != "collision_model" ||
            collision_model != configured_model ||
            collision_signature != configured_signature) {
            throw std::runtime_error(
                "unstructured checkpoint collision model mismatch");
        }
        std::size_t channel_count = 0;
        input >> key >> collision_totals_.candidates
              >> collision_totals_.null_collisions >> channel_count;
        if (key != "collision_totals" ||
            channel_count !=
                collision_totals_.channel_collisions.size()) {
            throw std::runtime_error(
                "unstructured checkpoint collision totals mismatch");
        }
        for (auto& count : collision_totals_.channel_collisions) {
            input >> count;
        }
        clear_collision_counts(collision_interval_);
    } else if (mcc_model_) {
        throw std::runtime_error(
            "legacy unstructured checkpoint cannot restart MCC");
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
    if (checkpoint_v1) {
        if (!sources_.empty() || !emissions_.empty()) {
            throw std::runtime_error(
                "legacy unstructured checkpoint cannot restart configured "
                "sources or emissions");
        }
    } else {
        std::size_t source_count = 0;
        input >> key >> source_count;
        if (key != "source_count" || source_count != sources_.size()) {
            throw std::runtime_error(
                "unstructured checkpoint source count mismatch");
        }
        std::set<std::string> loaded_sources;
        for (std::size_t i = 0; i < source_count; ++i) {
            std::string name;
            std::string species;
            std::string boundary;
            std::size_t particles_per_step = 0;
            std::size_t start_step = 0;
            std::size_t end_step = 0;
            double normal_velocity = 0.0;
            double tangential_velocity = 0.0;
            double thermal_velocity = 0.0;
            double out_of_plane_velocity = 0.0;
            std::size_t injected_particles = 0;
            input >> key >> std::quoted(name) >> std::quoted(species)
                  >> std::quoted(boundary) >> particles_per_step
                  >> start_step >> end_step >> normal_velocity
                  >> tangential_velocity >> thermal_velocity;
            if (checkpoint_v5 || checkpoint_v6) {
                input >> out_of_plane_velocity;
            }
            input >> injected_particles;
            const auto source = std::find_if(
                sources_.begin(), sources_.end(),
                [&](const BoundarySourceRuntime& candidate) {
                    return candidate.config.name == name;
                });
            if (key != "source" || source == sources_.end() ||
                !loaded_sources.insert(name).second ||
                source->config.species != species ||
                source->config.boundary != boundary ||
                source->config.particles_per_step != particles_per_step ||
                source->config.start_step != start_step ||
                source->config.end_step != end_step ||
                source->config.normal_velocity != normal_velocity ||
                source->config.tangential_velocity != tangential_velocity ||
                source->config.thermal_velocity != thermal_velocity ||
                source->config.out_of_plane_velocity !=
                    out_of_plane_velocity) {
                throw std::runtime_error(
                    "unstructured checkpoint source configuration mismatch");
            }
            source->injected_particles = injected_particles;
        }
    }
    if (!checkpoint_v3 && !checkpoint_v4 && !checkpoint_v5 &&
        !checkpoint_v6) {
        if (!emissions_.empty()) {
            throw std::runtime_error(
                "legacy unstructured checkpoint cannot restart configured "
                "emissions");
        }
        for (auto& [species, boundaries] : impact_flux_) {
            (void)species;
            for (auto& [boundary, flux] : boundaries) {
                (void)boundary;
                flux = {};
            }
        }
    } else {
        std::size_t emission_count = 0;
        input >> key >> emission_count;
        if (key != "emission_count" ||
            emission_count != emissions_.size()) {
            throw std::runtime_error(
                "unstructured checkpoint emission count mismatch");
        }
        std::set<std::string> loaded_emissions;
        for (std::size_t i = 0; i < emission_count; ++i) {
            std::string name;
            std::string boundary;
            std::string incident_species;
            std::string emitted_species;
            double yield = 0.0;
            std::size_t max_particles_per_impact = 0;
            double normal_velocity = 0.0;
            double tangential_velocity = 0.0;
            double thermal_velocity = 0.0;
            double out_of_plane_velocity = 0.0;
            std::size_t emitted_particles = 0;
            input >> key >> std::quoted(name) >> std::quoted(boundary)
                  >> std::quoted(incident_species)
                  >> std::quoted(emitted_species) >> yield
                  >> max_particles_per_impact >> normal_velocity
                  >> tangential_velocity >> thermal_velocity;
            if (checkpoint_v5 || checkpoint_v6) {
                input >> out_of_plane_velocity;
            }
            input >> emitted_particles;
            const auto emission = std::find_if(
                emissions_.begin(), emissions_.end(),
                [&](const SecondaryEmissionRuntime& candidate) {
                    return candidate.config.name == name;
                });
            if (key != "emission" || emission == emissions_.end() ||
                !loaded_emissions.insert(name).second ||
                emission->config.boundary != boundary ||
                emission->config.incident_species != incident_species ||
                emission->config.emitted_species != emitted_species ||
                emission->config.yield != yield ||
                emission->config.max_particles_per_impact !=
                    max_particles_per_impact ||
                emission->config.normal_velocity != normal_velocity ||
                emission->config.tangential_velocity !=
                    tangential_velocity ||
                emission->config.thermal_velocity != thermal_velocity ||
                emission->config.out_of_plane_velocity !=
                    out_of_plane_velocity) {
                throw std::runtime_error(
                    "unstructured checkpoint emission configuration mismatch");
            }
            emission->emitted_particles = emitted_particles;
        }
        std::size_t flux_count = 0;
        input >> key >> flux_count;
        std::size_t expected_flux_count = 0;
        for (const auto& [species, boundaries] : impact_flux_) {
            (void)species;
            expected_flux_count += boundaries.size();
        }
        if (key != "impact_flux_count" ||
            flux_count != expected_flux_count) {
            throw std::runtime_error(
                "unstructured checkpoint impact-flux count mismatch");
        }
        std::set<std::pair<std::string, std::string>> loaded_fluxes;
        for (std::size_t i = 0; i < flux_count; ++i) {
            std::string species;
            std::string boundary;
            UnstructuredBoundaryFlux2D flux;
            input >> key >> std::quoted(species) >> std::quoted(boundary)
                  >> flux.macroparticles
                  >> flux.physical_particles
                  >> flux.charge
                  >> flux.kinetic_energy
                  >> flux.last_step_macroparticles
                  >> flux.last_step_physical_particles
                  >> flux.physical_particle_rate
                  >> flux.physical_particle_flux;
            if (key != "impact_flux" ||
                !impact_flux_.contains(species) ||
                !impact_flux_.at(species).contains(boundary) ||
                !loaded_fluxes.emplace(species, boundary).second ||
                !std::isfinite(flux.physical_particles) ||
                !std::isfinite(flux.charge) ||
                !std::isfinite(flux.kinetic_energy) ||
                !std::isfinite(flux.last_step_physical_particles) ||
                !std::isfinite(flux.physical_particle_rate) ||
                !std::isfinite(flux.physical_particle_flux) ||
                flux.physical_particles < 0.0 ||
                flux.kinetic_energy < 0.0 ||
                flux.last_step_physical_particles < 0.0 ||
                flux.physical_particle_rate < 0.0 ||
                flux.physical_particle_flux < 0.0 ||
                flux.last_step_macroparticles > flux.macroparticles ||
                flux.last_step_physical_particles >
                    flux.physical_particles) {
                throw std::runtime_error(
                    "unstructured checkpoint impact-flux data mismatch");
            }
            impact_flux_.at(species).at(boundary) = flux;
        }
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
            particle_count < species_configs_[species_id].particles ||
            particle_count > config_.max_particles_per_species) {
            throw std::runtime_error("unstructured checkpoint species metadata mismatch");
        }
        auto& particles = species_[species_id].particles();
        particles.resize(particle_count);
        particle_locations_[species_id].assign(
            particle_count, UnstructuredParticleLocation2D{});
        for (auto& particle : particles) {
            int alive = 0;
            input >> particle.position.x >> particle.position.y
                  >> particle.velocity.x >> particle.velocity.y;
            if (checkpoint_v5 || checkpoint_v6) {
                input >> particle.velocity_z;
            } else {
                particle.velocity_z = 0.0;
            }
            input >> particle.velocity_half.x >> particle.velocity_half.y;
            if (checkpoint_v5 || checkpoint_v6) {
                input >> particle.velocity_half_z;
            } else {
                particle.velocity_half_z = 0.0;
            }
            input >> alive;
            particle.alive = alive != 0;
            if (!std::isfinite(particle.position.x) ||
                !std::isfinite(particle.position.y) ||
                !std::isfinite(particle.velocity.x) ||
                !std::isfinite(particle.velocity.y) ||
                !std::isfinite(particle.velocity_z) ||
                !std::isfinite(particle.velocity_half.x) ||
                !std::isfinite(particle.velocity_half.y) ||
                !std::isfinite(particle.velocity_half_z) ||
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
    write_unit_metadata(config_.output_dir, config_.units, 2);
    if (!config_.initial_state_path.empty()) {
        write_external_particle_state_metadata(
            config_.output_dir /
                "initial_state_metadata.txt",
            config_.initial_state_path,
            initial_state_metadata_,
            config_.initial_state_signature);
    }
    std::vector<InitializationSpeciesMoments> initialization_moments;
    initialization_moments.reserve(species_.size());
    for (std::size_t species_id = 0;
         species_id < species_.size(); ++species_id) {
        initialization_moments.push_back(
            summarize_initialization(
                species_[species_id],
                config_.restart_path.empty()
                    ? species_configs_[species_id]
                          .initialization_region
                    : std::string{}));
    }
    write_initialization_report(
        config_.output_dir / "initialization.csv", 2,
        !config_.restart_path.empty()
            ? "restart"
            : (!config_.initial_state_path.empty()
                   ? "external"
                   : "generated"),
        initialization_moments);
    const auto initialization_acceptance =
        assess_initialization_acceptance(
            config_.initialization_acceptance,
            initialization_moments, 3);
    write_initialization_acceptance_report(
        config_.output_dir /
            "initialization_acceptance.csv",
        initialization_acceptance);
    enforce_initialization_acceptance(
        initialization_acceptance);
    std::ofstream diagnostics(config_.output_dir / "scalars.csv");
    if (!diagnostics) throw std::runtime_error("cannot open unstructured diagnostics output");
    write_diagnostics_header(diagnostics);
    std::ofstream collision_output;
    if (mcc_model_) {
        write_collision_metadata(
            config_.output_dir, config_.collisions,
            mcc_model_->signature(),
            mcc_model_->neutral_velocity_stddev(),
            mcc_model_->neutral_speed_limit_sigma());
        collision_output.open(config_.output_dir / "collisions.csv");
        if (!collision_output) {
            throw std::runtime_error(
                "cannot open imported collision diagnostics output");
        }
        write_collision_header(
            collision_output, collision_totals_);
        write_collision_sample(
            collision_output, step_, time_,
            collision_interval_, collision_totals_);
        clear_collision_counts(collision_interval_);
    }

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
        if (mcc_model_) {
            write_collision_sample(
                collision_output, step_, time_,
                collision_interval_, collision_totals_);
            clear_collision_counts(collision_interval_);
        }
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
