#include "pic/Simulation2D.hpp"
#include "pic/Convergence.hpp"
#include "pic/ParticleState.hpp"
#include "pic/Pusher.hpp"
#include "pic/Runtime.hpp"
#include "pic/Units.hpp"
#include "pic/VTKWriter.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <utility>

namespace pic {
namespace {
constexpr const char* kCheckpointMagicV1 = "AuroraPIC-checkpoint-v1";
constexpr const char* kCheckpointMagicV2 = "AuroraPIC-checkpoint-v2";
constexpr const char* kCheckpointMagicV3 = "AuroraPIC-checkpoint-v3";
constexpr const char* kCheckpointMagicV4 = "AuroraPIC-checkpoint-v4";
constexpr const char* kCheckpointMagicV5 = "AuroraPIC-checkpoint-v5";
constexpr const char* kCheckpointMagicV6 = "AuroraPIC-checkpoint-v6";

double wrap_periodic(double value, double length) {
    return std::fmod(std::fmod(value, length) + length, length);
}

ParticleBoundary resolve_particle_boundary(ParticleBoundary configured, Boundary field_boundary) {
    if (configured != ParticleBoundary::Auto) return configured;
    return field_boundary == Boundary::Periodic ? ParticleBoundary::Periodic : ParticleBoundary::Absorbing;
}

void reflect_coordinate(double& coordinate, double& velocity, double length) {
    while (coordinate < 0.0 || coordinate > length) {
        if (coordinate < 0.0) {
            coordinate = -coordinate;
            velocity = -velocity;
        } else {
            coordinate = 2.0 * length - coordinate;
            velocity = -velocity;
        }
    }
}

bool apply_lower_boundary(double& coordinate,
                          double& velocity,
                          double length,
                          ParticleBoundary boundary,
                          std::size_t& absorbed_count) {
    if (coordinate >= 0.0) return true;
    switch (boundary) {
        case ParticleBoundary::Absorbing:
            ++absorbed_count;
            return false;
        case ParticleBoundary::Reflecting:
            reflect_coordinate(coordinate, velocity, length);
            return true;
        case ParticleBoundary::Periodic:
            coordinate = wrap_periodic(coordinate, length);
            return true;
        case ParticleBoundary::Auto:
            throw std::logic_error("unresolved lower 2D particle boundary policy");
    }
    return true;
}

bool apply_upper_boundary(double& coordinate,
                          double& velocity,
                          double length,
                          ParticleBoundary boundary,
                          std::size_t& absorbed_count) {
    if (coordinate <= length) return true;
    switch (boundary) {
        case ParticleBoundary::Absorbing:
            ++absorbed_count;
            return false;
        case ParticleBoundary::Reflecting:
            reflect_coordinate(coordinate, velocity, length);
            return true;
        case ParticleBoundary::Periodic:
            coordinate = wrap_periodic(coordinate, length);
            return true;
        case ParticleBoundary::Auto:
            throw std::logic_error("unresolved upper 2D particle boundary policy");
    }
    return true;
}

std::filesystem::path checkpoint_path_for_step(const Simulation2DConfig& cfg, std::size_t step) {
    if (!cfg.checkpoint_path.empty()) return cfg.checkpoint_path;
    return cfg.output_dir / ("checkpoint_" + std::to_string(step) + ".apc");
}

void ensure_parent_directory(const std::filesystem::path& path) {
    const auto parent = path.parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
}

template <typename T>
void require_stream(T& stream, const std::string& message) {
    if (!stream) throw std::runtime_error(message);
}

bool has_magnetic_field(const Simulation2DConfig& cfg) {
    return cfg.magnetic_field_profile.has_value() ||
           cfg.magnetic_field_x != 0.0 ||
           cfg.magnetic_field_y != 0.0 ||
           cfg.magnetic_field_z != 0.0;
}

Vec3 magnetic_field(
    const Simulation2DConfig& cfg,
    Vec2 position) {
    if (cfg.magnetic_field_profile) {
        return cfg.magnetic_field_profile->evaluate(
            {position.x, position.y, 0.0});
    }
    return {
        cfg.magnetic_field_x,
        cfg.magnetic_field_y,
        cfg.magnetic_field_z};
}

void initialize_particle_pusher(Particle2D& particle, Vec2 electric, double charge_to_mass, const Simulation2DConfig& cfg) {
    if (has_magnetic_field(cfg)) {
        initialize_boris_half_step(
            particle, electric, magnetic_field(cfg, particle.position),
            charge_to_mass, cfg.dt);
    } else {
        initialize_leapfrog_half_step(particle, electric, charge_to_mass, cfg.dt);
    }
}

void kick_particle(Particle2D& particle, Vec2 electric, double charge_to_mass, const Simulation2DConfig& cfg) {
    if (has_magnetic_field(cfg)) {
        kick_boris(
            particle, electric, magnetic_field(cfg, particle.position),
            charge_to_mass, cfg.dt);
    } else {
        kick_leapfrog(particle, electric, charge_to_mass, cfg.dt);
    }
}

void synchronize_particle(Particle2D& particle, Vec2 electric, double charge_to_mass, const Simulation2DConfig& cfg) {
    if (has_magnetic_field(cfg)) {
        synchronize_boris(
            particle, electric, magnetic_field(cfg, particle.position),
            charge_to_mass, cfg.dt);
    } else {
        synchronize_leapfrog(particle, electric, charge_to_mass, cfg.dt);
    }
}

void write_vtk_outputs(const Mesh2D& mesh, const std::filesystem::path& output_dir, std::size_t step, VTKOutputFormat format) {
    const auto stem = output_dir / ("fields_" + std::to_string(step));
    if (format == VTKOutputFormat::Legacy || format == VTKOutputFormat::Both) {
        write_legacy_vtk(mesh, stem.string() + ".vtk");
    }
    if (format == VTKOutputFormat::Xml || format == VTKOutputFormat::Both) {
        write_vtk_xml(mesh, stem.string() + ".vts");
    }
}

double source_profile_integral(
    const VolumetricPairSource2DConfig& source,
    double xmax, double ymax) {
    const double width = xmax - source.x_min;
    const double height = ymax - source.y_min;
    if (source.spatial_profile.density_profile ==
        DensityProfileKind::Uniform) {
        return width * height;
    }
    if (source.spatial_profile.density_profile ==
        DensityProfileKind::Sinusoidal) {
        return width * height /
            (1.0 + std::abs(
                *source.spatial_profile.profile_amplitude));
    }
    const auto gaussian_integral =
        [](double minimum, double maximum,
           double center, double scale) {
            const double denominator =
                std::numbers::sqrt2 * scale;
            return scale *
                std::sqrt(std::numbers::pi / 2.0) *
                (std::erf((maximum - center) / denominator) -
                 std::erf((minimum - center) / denominator));
        };
    return gaussian_integral(
               source.x_min, xmax,
               *source.spatial_profile.profile_center_x,
               *source.spatial_profile.profile_scale_x) *
           gaussian_integral(
               source.y_min, ymax,
               *source.spatial_profile.profile_center_y,
               *source.spatial_profile.profile_scale_y);
}
} // namespace

Simulation2D::Simulation2D(Simulation2DConfig cfg)
    : cfg_(std::move(cfg)),
      mesh_(cfg_.nx, cfg_.ny, cfg_.length_x, cfg_.length_y,
            cfg_.boundary_x.value_or(cfg_.boundary),
            cfg_.boundary_y.value_or(cfg_.boundary),
            cfg_.boundary_config),
      solver_(cfg_.units.permittivity()),
      rng_(cfg_.seed) {
    if (cfg_.checkpoint_output && cfg_.checkpoint_interval == 0) cfg_.checkpoint_interval = cfg_.output_interval;
    if (!cfg_.restart_path.empty() &&
        !cfg_.initial_state_path.empty()) {
        throw std::invalid_argument(
            "2D restart_path and initial_state_path are mutually exclusive");
    }
    if (cfg_.initial_state_signature &&
        cfg_.initial_state_path.empty()) {
        throw std::invalid_argument(
            "2D initial_state_signature requires initial_state_path");
    }
    validate_initialization_acceptance(
        cfg_.initialization_acceptance,
        "2D initialization acceptance config");
    if (!std::isfinite(cfg_.dt) || cfg_.dt <= 0.0) throw std::invalid_argument("2D simulation dt must be positive and finite");
    if (!std::isfinite(cfg_.out_of_plane_depth) ||
        !(cfg_.out_of_plane_depth > 0.0)) {
        throw std::invalid_argument(
            "2D out_of_plane_depth must be positive and finite");
    }
    if (cfg_.output_interval == 0) throw std::invalid_argument("2D output_interval must be positive");
    if (cfg_.particle_output_stride == 0) throw std::invalid_argument("2D particle_output_stride must be positive");
    if (!std::isfinite(cfg_.magnetic_field_x) ||
        !std::isfinite(cfg_.magnetic_field_y) ||
        !std::isfinite(cfg_.magnetic_field_z)) {
        throw std::invalid_argument(
            "2D magnetic_field components must be finite");
    }
    if (cfg_.magnetic_field_profile) {
        if (cfg_.magnetic_field_x != 0.0 ||
            cfg_.magnetic_field_y != 0.0 ||
            cfg_.magnetic_field_z != 0.0) {
            throw std::invalid_argument(
                "2D uniform magnetic_field components and magnetic_field_profile are mutually exclusive");
        }
        cfg_.magnetic_field_profile->validate_domain(
            {0.0, 0.0, 0.0},
            {cfg_.length_x, cfg_.length_y, 0.0},
            "2D simulation");
    }
    if (cfg_.mode == RunMode::SteadyState) {
        if (cfg_.max_steps == 0) throw std::invalid_argument("2D max_steps must be positive for steady-state mode");
        if (cfg_.steady_window == 0) throw std::invalid_argument("2D steady_window must be positive");
        if (!std::isfinite(cfg_.steady_tolerance) || cfg_.steady_tolerance <= 0.0) {
            throw std::invalid_argument("2D steady_tolerance must be positive and finite");
        }
    }
    validate_runtime_policy(cfg_.runtime);
    if (cfg_.max_particles_per_species == 0) {
        throw std::invalid_argument(
            "2D max_particles_per_species must be positive");
    }
    if (cfg_.checkpoint_output && cfg_.checkpoint_interval == 0) {
        throw std::invalid_argument("2D checkpoint_interval must be positive when checkpoint_output is enabled");
    }
    cfg_.particle_boundary_config.left =
        resolve_particle_boundary(
            cfg_.particle_boundary_config.left,
            mesh_.boundary_x());
    cfg_.particle_boundary_config.right =
        resolve_particle_boundary(
            cfg_.particle_boundary_config.right,
            mesh_.boundary_x());
    cfg_.particle_boundary_config.bottom =
        resolve_particle_boundary(
            cfg_.particle_boundary_config.bottom,
            mesh_.boundary_y());
    cfg_.particle_boundary_config.top =
        resolve_particle_boundary(
            cfg_.particle_boundary_config.top,
            mesh_.boundary_y());
    for (const auto& sc : cfg_.species) species_.emplace_back(sc);
    if (species_.empty()) species_.emplace_back(Species2DConfig{});
    for (const auto& species : species_) {
        if (species.config().particles >
            cfg_.max_particles_per_species) {
            throw std::invalid_argument(
                "2D initial species population exceeds max_particles_per_species");
        }
    }
    for (const auto& source : cfg_.sources) {
        if (!std::all_of(
                source.name.begin(), source.name.end(),
                [](unsigned char character) {
                    return std::isalnum(character) ||
                           character == '_' || character == '-';
                })) {
            throw std::invalid_argument(
                "2D source names may contain only letters, digits, '_' and '-'");
        }
        const bool fixed_rate = source.pairs_per_step != 0;
        const bool physical_rate =
            source.represented_pair_rate.has_value();
        const bool volumetric_rate =
            source.peak_volumetric_pair_rate.has_value();
        if (source.name.empty() ||
            static_cast<unsigned>(fixed_rate) +
                    static_cast<unsigned>(physical_rate) +
                    static_cast<unsigned>(volumetric_rate) !=
                1U) {
            throw std::invalid_argument(
                "2D volumetric pair source requires a name and exactly one rate specification");
        }
        if (std::any_of(
                sources_.begin(), sources_.end(),
                [&](const auto& existing) {
                    return existing.config.name == source.name;
                })) {
            throw std::invalid_argument(
                "duplicate 2D source name '" + source.name + "'");
        }
        const double xmax =
            source.x_max < 0.0 ? mesh_.length_x() : source.x_max;
        const double ymax =
            source.y_max < 0.0 ? mesh_.length_y() : source.y_max;
        if (source.x_min < 0.0 || source.y_min < 0.0 ||
            xmax > mesh_.length_x() || ymax > mesh_.length_y() ||
            !(source.x_min < xmax) || !(source.y_min < ymax)) {
            throw std::invalid_argument(
                "2D source '" + source.name +
                "' region must have positive area inside the domain");
        }
        if (source.end_step != 0 &&
            source.end_step <= source.start_step) {
            throw std::invalid_argument(
                "2D source '" + source.name +
                "' end_step must exceed start_step");
        }
        if (!std::isfinite(source.first_thermal_velocity) ||
            source.first_thermal_velocity < 0.0 ||
            !std::isfinite(source.second_thermal_velocity) ||
            source.second_thermal_velocity < 0.0) {
            throw std::invalid_argument(
                "2D source '" + source.name +
                "' thermal velocities must be non-negative and finite");
        }
        const auto finite_vector = [](Vec3 vector) {
            return std::isfinite(vector.x) &&
                   std::isfinite(vector.y) &&
                   std::isfinite(vector.z);
        };
        if (!finite_vector(source.first_drift) ||
            !finite_vector(source.second_drift)) {
            throw std::invalid_argument(
                "2D source '" + source.name +
                "' drift velocities must be finite");
        }
        validate_density_profile(
            source.spatial_profile, 2, 1,
            "2D source '" + source.name + "' spatial profile");
        const auto find_species =
            [&](const std::string& name) -> std::size_t {
                const auto found = std::find_if(
                    species_.begin(), species_.end(),
                    [&](const Species2D& species) {
                        return species.name() == name;
                    });
                if (found == species_.end()) {
                    throw std::invalid_argument(
                        "2D source '" + source.name +
                        "' references unknown species '" + name + "'");
                }
                return static_cast<std::size_t>(
                    std::distance(species_.begin(), found));
            };
        const std::size_t first =
            find_species(source.first_species);
        const std::size_t second =
            find_species(source.second_species);
        const auto& first_species = species_[first];
        const auto& second_species = species_[second];
        const double weight_scale = std::max(
            first_species.weight(), second_species.weight());
        if (std::abs(first_species.weight() -
                     second_species.weight()) >
            1e-12 * weight_scale) {
            throw std::invalid_argument(
                "2D source '" + source.name +
                "' requires equal macro-particle weights");
        }
        if ((source.represented_pair_rate &&
             (!std::isfinite(*source.represented_pair_rate) ||
              !(*source.represented_pair_rate > 0.0))) ||
            (source.peak_volumetric_pair_rate &&
             (!std::isfinite(*source.peak_volumetric_pair_rate) ||
              !(*source.peak_volumetric_pair_rate > 0.0)))) {
                throw std::invalid_argument(
                    "2D source '" + source.name +
                    "' configured rate must be positive and finite");
        }
        const double effective_area =
            source_profile_integral(source, xmax, ymax);
        const double resolved_rate =
            source.represented_pair_rate.value_or(
                source.peak_volumetric_pair_rate.value_or(0.0) *
                effective_area * cfg_.out_of_plane_depth);
        if (physical_rate || volumetric_rate) {
            const double macro_pairs_per_step =
                resolved_rate * cfg_.dt /
                first_species.weight();
            if (!std::isfinite(macro_pairs_per_step) ||
                macro_pairs_per_step >
                    static_cast<double>(
                        cfg_.max_particles_per_species)) {
                throw std::invalid_argument(
                    "2D source '" + source.name +
                    "' physical rate exceeds per-step storage capacity");
            }
        }
        const double charge_scale = std::max(
            std::abs(first_species.charge()),
            std::abs(second_species.charge()));
        if (!((first_species.charge() < 0.0 &&
               second_species.charge() > 0.0) ||
              (first_species.charge() > 0.0 &&
               second_species.charge() < 0.0)) ||
            charge_scale == 0.0 ||
            std::abs(first_species.charge() +
                     second_species.charge()) >
                1e-12 * charge_scale) {
            throw std::invalid_argument(
                "2D source '" + source.name +
                "' requires opposite equal species charges");
        }
        sources_.push_back(
            {source, first, second, resolved_rate, effective_area});
        source_diagnostics_.push_back(
            {source.name, 0, 0.0, 0.0, 0.0});
    }
}

void Simulation2D::initialize() {
    time_ = 0.0;
    step_ = 0;
    boundary_losses_ = {};
    for (auto& diagnostics : source_diagnostics_) {
        diagnostics.macro_pairs_created = 0;
        diagnostics.represented_pairs_created = 0.0;
        diagnostics.fractional_macro_pair_remainder = 0.0;
        diagnostics.injected_kinetic_energy = 0.0;
    }
    if (cfg_.initial_state_path.empty()) {
        for (auto& sp : species_) sp.initialize(mesh_, rng_);
    } else {
        std::vector<ExternalSpeciesExpectation> expected;
        expected.reserve(species_.size());
        for (const auto& species : species_) {
            expected.push_back({
                species.name(),
                species.config().particles});
        }
        for (auto& species : species_) {
            species.particles().resize(
                species.config().particles);
        }
        initial_state_metadata_ =
            load_validated_external_particle_state_bounded(
                cfg_.initial_state_path, 2,
                cfg_.units.system, expected,
                "2D simulation",
                [&](std::size_t species_index,
                    std::size_t record_index,
                    const ExternalParticleRecord& record) {
                    auto& species =
                        species_.at(species_index);
                    const auto& species_config =
                        species.config();
                    const double maximum_x =
                        species_config.init_x_max < 0.0
                            ? mesh_.length_x()
                            : species_config.init_x_max;
                    const double maximum_y =
                        species_config.init_y_max < 0.0
                            ? mesh_.length_y()
                            : species_config.init_y_max;
                    const bool outside =
                        record.position.x <
                            species_config.init_x_min ||
                        record.position.x > maximum_x ||
                        record.position.y <
                            species_config.init_y_min ||
                        record.position.y > maximum_y ||
                        (mesh_.boundary_x() == Boundary::Periodic &&
                         record.position.x == mesh_.length_x()) ||
                        (mesh_.boundary_y() == Boundary::Periodic &&
                         record.position.y == mesh_.length_y());
                    if (outside) {
                        throw std::runtime_error(
                            "external particle for species '" +
                            species.name() +
                            "' lies outside the 2D domain");
                    }
                    auto& particle =
                        species.particles().at(record_index);
                    particle.position = {
                        record.position.x,
                        record.position.y};
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
                },
                cfg_.initial_state_signature);
    }
    deposit_and_solve();
    for (auto& sp : species_) {
        const double qm = sp.charge() / sp.mass();
        auto& particles = sp.particles();
        runtime_parallel_for(std::size_t{0}, particles.size(), cfg_.runtime, [&](std::size_t particle_id) {
            auto& particle = particles[particle_id];
            if (particle.alive) initialize_particle_pusher(particle, interpolate_electric(mesh_, particle.position), qm, cfg_);
        });
    }
    initialized_ = true;
}

void Simulation2D::deposit_and_solve() {
    mesh_.clear_charge();
    for (const auto& sp : species_) {
        sp.deposit_charge(mesh_, cfg_.out_of_plane_depth);
    }
    solver_.solve(mesh_);
}

void Simulation2D::inject_volumetric_pair_sources() {
    const auto available_slots =
        [&](const Species2D& species) {
            const auto dead = static_cast<std::size_t>(
                std::count_if(
                    species.particles().begin(),
                    species.particles().end(),
                    [](const Particle2D& particle) {
                        return !particle.alive;
                    }));
            return dead +
                (cfg_.max_particles_per_species -
                 std::min(cfg_.max_particles_per_species,
                          species.particles().size()));
        };
    const auto insert_particle =
        [](Species2D& species, Particle2D particle) {
            const auto dead = std::find_if(
                species.particles().begin(),
                species.particles().end(),
                [](const Particle2D& candidate) {
                    return !candidate.alive;
                });
            if (dead != species.particles().end()) {
                *dead = particle;
            } else {
                species.particles().push_back(particle);
            }
        };

    for (std::size_t source_id = 0;
         source_id < sources_.size(); ++source_id) {
        const auto& source = sources_[source_id];
        if (step_ < source.config.start_step ||
            (source.config.end_step != 0 &&
             step_ >= source.config.end_step)) {
            continue;
        }
        auto& first = species_[source.first_species];
        auto& second = species_[source.second_species];
        auto& diagnostics = source_diagnostics_[source_id];
        std::size_t pairs_this_step =
            source.config.pairs_per_step;
        double next_fractional_remainder =
            diagnostics.fractional_macro_pair_remainder;
        if (source.config.represented_pair_rate ||
            source.config.peak_volumetric_pair_rate) {
            const double increment =
                source.represented_pair_rate * cfg_.dt /
                first.weight();
            const double accumulated =
                diagnostics.fractional_macro_pair_remainder +
                increment;
            if (!std::isfinite(accumulated) ||
                accumulated >
                    static_cast<double>(
                        std::numeric_limits<std::size_t>::max())) {
                throw std::runtime_error(
                    "2D source '" + source.config.name +
                    "' fractional rate accumulator overflow");
            }
            pairs_this_step =
                static_cast<std::size_t>(std::floor(accumulated));
            next_fractional_remainder =
                accumulated -
                static_cast<double>(pairs_this_step);
        }
        if (available_slots(first) <
                pairs_this_step ||
            available_slots(second) <
                pairs_this_step) {
            throw std::runtime_error(
                "2D source '" + source.config.name +
                "' would exceed max_particles_per_species");
        }
        if (pairs_this_step >
            std::numeric_limits<std::size_t>::max() -
                diagnostics.macro_pairs_created) {
            throw std::runtime_error(
                "2D source '" + source.config.name +
                "' macro-pair counter overflow");
        }
        const double represented_increment =
            static_cast<double>(pairs_this_step) *
            first.weight();
        if (!std::isfinite(represented_increment) ||
            !std::isfinite(
                diagnostics.represented_pairs_created +
                represented_increment)) {
            throw std::runtime_error(
                "2D source '" + source.config.name +
                "' represented-pair counter overflow");
        }
        const double xmax =
            source.config.x_max < 0.0
                ? mesh_.length_x()
                : source.config.x_max;
        const double ymax =
            source.config.y_max < 0.0
                ? mesh_.length_y()
                : source.config.y_max;
        std::uniform_real_distribution<double> x_distribution(
            source.config.x_min, xmax);
        std::uniform_real_distribution<double> y_distribution(
            source.config.y_min, ymax);
        std::uniform_real_distribution<double> unit_distribution(
            0.0, 1.0);
        std::normal_distribution<double> first_vx(
            source.config.first_drift.x,
            source.config.first_thermal_velocity);
        std::normal_distribution<double> first_vy(
            source.config.first_drift.y,
            source.config.first_thermal_velocity);
        std::normal_distribution<double> first_vz(
            source.config.first_drift.z,
            source.config.first_thermal_velocity);
        std::normal_distribution<double> second_vx(
            source.config.second_drift.x,
            source.config.second_thermal_velocity);
        std::normal_distribution<double> second_vy(
            source.config.second_drift.y,
            source.config.second_thermal_velocity);
        std::normal_distribution<double> second_vz(
            source.config.second_drift.z,
            source.config.second_thermal_velocity);

        std::vector<Vec2> positions;
        positions.reserve(pairs_this_step);
        std::size_t profile_attempts = 0;
        while (positions.size() < pairs_this_step) {
            Vec2 position{
                x_distribution(rng_), y_distribution(rng_)};
            if (source.config.spatial_profile.density_profile ==
                    DensityProfileKind::Uniform) {
                positions.push_back(position);
                continue;
            }
            while (true) {
                if (profile_attempts >=
                    source.config.spatial_profile
                        .max_profile_sampling_attempts) {
                    throw std::runtime_error(
                        "2D source '" + source.config.name +
                        "' spatial-profile sampling exceeded max_profile_sampling_attempts");
                }
                ++profile_attempts;
                if (unit_distribution(rng_) <=
                        density_profile_acceptance(
                            source.config.spatial_profile,
                            {position.x, position.y, 0.0},
                            {source.config.x_min,
                             source.config.y_min, 0.0},
                            {xmax, ymax, 1.0})) {
                    positions.push_back(position);
                    break;
                }
                position = {
                    x_distribution(rng_),
                    y_distribution(rng_)};
            }
        }
        double injected_energy_this_step = 0.0;
        for (const Vec2 position : positions) {
            Particle2D first_particle;
            first_particle.position = position;
            first_particle.velocity = {
                first_vx(rng_), first_vy(rng_)};
            first_particle.velocity_z = first_vz(rng_);
            first_particle.velocity_half =
                first_particle.velocity;
            first_particle.velocity_half_z =
                first_particle.velocity_z;
            initialize_particle_pusher(
                first_particle,
                interpolate_electric(mesh_, position),
                first.charge() / first.mass(), cfg_);

            Particle2D second_particle;
            second_particle.position = position;
            second_particle.velocity = {
                second_vx(rng_), second_vy(rng_)};
            second_particle.velocity_z = second_vz(rng_);
            second_particle.velocity_half =
                second_particle.velocity;
            second_particle.velocity_half_z =
                second_particle.velocity_z;
            const double first_speed_squared =
                first_particle.velocity.x *
                    first_particle.velocity.x +
                first_particle.velocity.y *
                    first_particle.velocity.y +
                first_particle.velocity_z *
                    first_particle.velocity_z;
            const double second_speed_squared =
                second_particle.velocity.x *
                    second_particle.velocity.x +
                second_particle.velocity.y *
                    second_particle.velocity.y +
                second_particle.velocity_z *
                    second_particle.velocity_z;
            injected_energy_this_step +=
                0.5 * first.weight() *
                (first.mass() * first_speed_squared +
                 second.mass() * second_speed_squared);
            if (!std::isfinite(injected_energy_this_step) ||
                !std::isfinite(
                    diagnostics.injected_kinetic_energy +
                    injected_energy_this_step)) {
                throw std::runtime_error(
                    "2D source '" + source.config.name +
                    "' injected-energy counter overflow");
            }
            initialize_particle_pusher(
                second_particle,
                interpolate_electric(mesh_, position),
                second.charge() / second.mass(), cfg_);
            insert_particle(first, first_particle);
            insert_particle(second, second_particle);
        }
        diagnostics.macro_pairs_created +=
            pairs_this_step;
        diagnostics.represented_pairs_created +=
            represented_increment;
        diagnostics.fractional_macro_pair_remainder =
            next_fractional_remainder;
        diagnostics.injected_kinetic_energy +=
            injected_energy_this_step;
    }
}

void Simulation2D::apply_particle_boundaries(Particle2D& particle) {
    if (!apply_lower_boundary(particle.position.x, particle.velocity_half.x, mesh_.length_x(),
                              cfg_.particle_boundary_config.left, boundary_losses_.absorbed_left)) {
        particle.alive = false;
        return;
    }
    if (!apply_upper_boundary(particle.position.x, particle.velocity_half.x, mesh_.length_x(),
                              cfg_.particle_boundary_config.right, boundary_losses_.absorbed_right)) {
        particle.alive = false;
        return;
    }
    if (!apply_lower_boundary(particle.position.y, particle.velocity_half.y, mesh_.length_y(),
                              cfg_.particle_boundary_config.bottom, boundary_losses_.absorbed_bottom)) {
        particle.alive = false;
        return;
    }
    if (!apply_upper_boundary(particle.position.y, particle.velocity_half.y, mesh_.length_y(),
                              cfg_.particle_boundary_config.top, boundary_losses_.absorbed_top)) {
        particle.alive = false;
        return;
    }
}

void Simulation2D::step() {
    if (!initialized_) initialize();
    inject_volumetric_pair_sources();
    for (auto& sp : species_) {
        const double qm = sp.charge() / sp.mass();
        for (auto& particle : sp.particles()) {
            if (!particle.alive) continue;
            kick_particle(particle, interpolate_electric(mesh_, particle.position), qm, cfg_);
            drift_leapfrog(particle, cfg_.dt);
            apply_particle_boundaries(particle);
        }
    }

    deposit_and_solve();
    for (auto& sp : species_) {
        const double qm = sp.charge() / sp.mass();
        auto& particles = sp.particles();
        runtime_parallel_for(std::size_t{0}, particles.size(), cfg_.runtime, [&](std::size_t particle_id) {
            auto& particle = particles[particle_id];
            if (particle.alive) synchronize_particle(particle, interpolate_electric(mesh_, particle.position), qm, cfg_);
        });
    }
    ++step_;
    time_ += cfg_.dt;
}

void Simulation2D::save_checkpoint(const std::filesystem::path& path) const {
    ensure_parent_directory(path);
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot open 2D checkpoint for writing: " + path.string());
    out << std::setprecision(17);
    out << kCheckpointMagicV6 << '\n';
    out << "dimension 2\n";
    out << "units " << to_string(cfg_.units.system) << ' '
        << cfg_.units.relative_permittivity << ' '
        << cfg_.units.permittivity() << ' '
        << cfg_.out_of_plane_depth << "\n";
    out << "step " << step_ << "\n";
    out << "time " << time_ << "\n";
    out << "boundary_losses " << boundary_losses_.absorbed_left << ' ' << boundary_losses_.absorbed_right << ' '
        << boundary_losses_.absorbed_bottom << ' ' << boundary_losses_.absorbed_top << "\n";
    out << "species_count " << species_.size() << "\n";
    out << "rng " << rng_ << "\n";
    out << "source_count " << source_diagnostics_.size() << "\n";
    for (std::size_t source_id = 0;
         source_id < sources_.size(); ++source_id) {
        const auto& config = sources_[source_id].config;
        const auto& diagnostics = source_diagnostics_[source_id];
        out << "source " << config.name << ' '
            << config.first_species << ' '
            << config.second_species << ' '
            << config.pairs_per_step << ' '
            << (config.represented_pair_rate ? 1 : 0) << ' '
            << config.represented_pair_rate.value_or(0.0) << ' '
            << (config.peak_volumetric_pair_rate ? 1 : 0) << ' '
            << config.peak_volumetric_pair_rate.value_or(0.0) << ' '
            << config.start_step << ' ' << config.end_step << ' '
            << config.x_min << ' ' << config.x_max << ' '
            << config.y_min << ' ' << config.y_max << ' '
            << config.first_drift.x << ' '
            << config.first_drift.y << ' '
            << config.first_drift.z << ' '
            << config.second_drift.x << ' '
            << config.second_drift.y << ' '
            << config.second_drift.z << ' '
            << config.first_thermal_velocity << ' '
            << config.second_thermal_velocity << ' '
            << to_string(
                   config.spatial_profile.density_profile) << ' ';
        const auto write_optional_double =
            [&](const std::optional<double>& value) {
                out << (value ? 1 : 0) << ' '
                    << value.value_or(0.0) << ' ';
            };
        const auto write_optional_size =
            [&](const std::optional<std::size_t>& value) {
                out << (value ? 1 : 0) << ' '
                    << value.value_or(0) << ' ';
            };
        write_optional_double(
            config.spatial_profile.profile_center_x);
        write_optional_double(
            config.spatial_profile.profile_center_y);
        write_optional_double(
            config.spatial_profile.profile_scale_x);
        write_optional_double(
            config.spatial_profile.profile_scale_y);
        write_optional_double(
            config.spatial_profile.profile_amplitude);
        write_optional_double(
            config.spatial_profile.profile_phase);
        write_optional_size(
            config.spatial_profile.profile_mode_x);
        write_optional_size(
            config.spatial_profile.profile_mode_y);
        out << config.spatial_profile
                   .max_profile_sampling_attempts << ' '
            << diagnostics.macro_pairs_created << ' '
            << diagnostics.represented_pairs_created << ' '
            << diagnostics.fractional_macro_pair_remainder << ' '
            << diagnostics.injected_kinetic_energy << "\n";
    }
    for (std::size_t species_id = 0; species_id < species_.size(); ++species_id) {
        const auto& sp = species_[species_id];
        out << "species " << species_id << ' ' << sp.name() << ' ' << sp.particles().size() << "\n";
        for (const auto& p : sp.particles()) {
            out << p.position.x << ' ' << p.position.y << ' '
                << p.velocity.x << ' ' << p.velocity.y << ' '
                << p.velocity_z << ' '
                << p.velocity_half.x << ' ' << p.velocity_half.y << ' '
                << p.velocity_half_z << ' '
                << (p.alive ? 1 : 0) << "\n";
        }
    }
    require_stream(out, "failed while writing 2D checkpoint: " + path.string());
}

void Simulation2D::load_checkpoint(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open 2D checkpoint for reading: " + path.string());
    std::string magic;
    std::getline(in, magic);
    const bool checkpoint_v1 = magic == kCheckpointMagicV1;
    const bool checkpoint_v2 = magic == kCheckpointMagicV2;
    const bool checkpoint_v3 = magic == kCheckpointMagicV3;
    const bool checkpoint_v4 = magic == kCheckpointMagicV4;
    const bool checkpoint_v5 = magic == kCheckpointMagicV5;
    const bool checkpoint_v6 = magic == kCheckpointMagicV6;
    if (!checkpoint_v1 && !checkpoint_v2 && !checkpoint_v3 &&
        !checkpoint_v4 && !checkpoint_v5 && !checkpoint_v6) {
        throw std::runtime_error("invalid checkpoint magic in: " + path.string());
    }

    std::string key;
    unsigned dimension = 0;
    in >> key >> dimension;
    if (key != "dimension" || dimension != 2) throw std::runtime_error("checkpoint dimension does not match 2D simulation");
    in >> key;
    if (key == "units") {
        std::string unit_system;
        double relative_permittivity = 0.0;
        double permittivity = 0.0;
        in >> unit_system >> relative_permittivity >> permittivity;
        double out_of_plane_depth = 1.0;
        if (checkpoint_v6) in >> out_of_plane_depth;
        if ((!checkpoint_v2 && !checkpoint_v3 && !checkpoint_v4 &&
             !checkpoint_v5 && !checkpoint_v6) ||
            unit_system != to_string(cfg_.units.system) ||
            relative_permittivity != cfg_.units.relative_permittivity ||
            permittivity != cfg_.units.permittivity() ||
            out_of_plane_depth != cfg_.out_of_plane_depth) {
            throw std::runtime_error(
                "checkpoint unit system does not match 2D config");
        }
        in >> key;
    } else if (checkpoint_v2 || checkpoint_v3 || checkpoint_v4 ||
               checkpoint_v5 ||
               checkpoint_v6 ||
               cfg_.units.system != UnitSystem::Normalized ||
               cfg_.units.relative_permittivity != 1.0) {
        throw std::runtime_error(
            "legacy checkpoint without unit metadata requires normalized units");
    }
    in >> step_;
    if (key != "step") throw std::runtime_error("checkpoint missing step");
    in >> key >> time_;
    if (key != "time") throw std::runtime_error("checkpoint missing time");
    in >> key >> boundary_losses_.absorbed_left >> boundary_losses_.absorbed_right
       >> boundary_losses_.absorbed_bottom >> boundary_losses_.absorbed_top;
    if (key != "boundary_losses") throw std::runtime_error("checkpoint missing 2D boundary loss counters");
    std::size_t species_count = 0;
    in >> key >> species_count;
    if (key != "species_count" || species_count != species_.size()) throw std::runtime_error("checkpoint species count does not match 2D config");
    in >> key;
    if (key != "rng") throw std::runtime_error("checkpoint missing rng state");
    in >> rng_;
    if (checkpoint_v4 || checkpoint_v5 || checkpoint_v6) {
        std::size_t source_count = 0;
        in >> key >> source_count;
        if (key != "source_count" ||
            source_count != source_diagnostics_.size()) {
            throw std::runtime_error(
                "checkpoint source count does not match 2D config");
        }
        for (std::size_t source_id = 0;
             source_id < sources_.size(); ++source_id) {
            VolumetricPairSource2DConfig stored;
            auto& diagnostics = source_diagnostics_[source_id];
            in >> key >> stored.name >> stored.first_species >>
                stored.second_species >> stored.pairs_per_step;
            if (checkpoint_v5 || checkpoint_v6) {
                int has_rate = 0;
                double rate = 0.0;
                in >> has_rate >> rate;
                if (has_rate != 0 && has_rate != 1) {
                    throw std::runtime_error(
                        "checkpoint source rate-presence flag is invalid");
                }
                if (has_rate != 0) {
                    stored.represented_pair_rate = rate;
                }
                if (checkpoint_v6) {
                    int has_peak_rate = 0;
                    double peak_rate = 0.0;
                    in >> has_peak_rate >> peak_rate;
                    if (has_peak_rate != 0 &&
                        has_peak_rate != 1) {
                        throw std::runtime_error(
                            "checkpoint source peak-rate-presence flag is invalid");
                    }
                    if (has_peak_rate != 0) {
                        stored.peak_volumetric_pair_rate =
                            peak_rate;
                    }
                }
            }
            in >> stored.start_step >> stored.end_step >>
                stored.x_min >> stored.x_max >>
                stored.y_min >> stored.y_max >>
                stored.first_drift.x >> stored.first_drift.y >>
                stored.first_drift.z >> stored.second_drift.x >>
                stored.second_drift.y >> stored.second_drift.z >>
                stored.first_thermal_velocity >>
                stored.second_thermal_velocity;
            if (checkpoint_v5 || checkpoint_v6) {
                std::string profile;
                in >> profile;
                stored.spatial_profile.density_profile =
                    density_profile_from_string(profile);
                const auto read_optional_double =
                    [&](std::optional<double>& value) {
                        int present = 0;
                        double stored_value = 0.0;
                        in >> present >> stored_value;
                        if (present != 0 && present != 1) {
                            throw std::runtime_error(
                                "checkpoint source optional-value flag is invalid");
                        }
                        if (present != 0) value = stored_value;
                    };
                const auto read_optional_size =
                    [&](std::optional<std::size_t>& value) {
                        int present = 0;
                        std::size_t stored_value = 0;
                        in >> present >> stored_value;
                        if (present != 0 && present != 1) {
                            throw std::runtime_error(
                                "checkpoint source optional-value flag is invalid");
                        }
                        if (present != 0) value = stored_value;
                    };
                read_optional_double(
                    stored.spatial_profile.profile_center_x);
                read_optional_double(
                    stored.spatial_profile.profile_center_y);
                read_optional_double(
                    stored.spatial_profile.profile_scale_x);
                read_optional_double(
                    stored.spatial_profile.profile_scale_y);
                read_optional_double(
                    stored.spatial_profile.profile_amplitude);
                read_optional_double(
                    stored.spatial_profile.profile_phase);
                read_optional_size(
                    stored.spatial_profile.profile_mode_x);
                read_optional_size(
                    stored.spatial_profile.profile_mode_y);
                in >> stored.spatial_profile
                          .max_profile_sampling_attempts;
            }
            in >> diagnostics.macro_pairs_created >>
                diagnostics.represented_pairs_created;
            if (checkpoint_v5 || checkpoint_v6) {
                in >> diagnostics
                          .fractional_macro_pair_remainder >>
                    diagnostics.injected_kinetic_energy;
            } else {
                diagnostics.fractional_macro_pair_remainder =
                    0.0;
                diagnostics.injected_kinetic_energy = 0.0;
            }
            if (!std::isfinite(
                    diagnostics.represented_pairs_created) ||
                diagnostics.represented_pairs_created < 0.0 ||
                !std::isfinite(
                    diagnostics.fractional_macro_pair_remainder) ||
                diagnostics.fractional_macro_pair_remainder < 0.0 ||
                !(diagnostics.fractional_macro_pair_remainder <
                  1.0) ||
                !std::isfinite(
                    diagnostics.injected_kinetic_energy) ||
                diagnostics.injected_kinetic_energy < 0.0) {
                throw std::runtime_error(
                    "checkpoint source diagnostics are invalid");
            }
            const auto& expected = sources_[source_id].config;
            if (key != "source" ||
                stored.name != expected.name ||
                stored.first_species != expected.first_species ||
                stored.second_species != expected.second_species ||
                stored.pairs_per_step != expected.pairs_per_step ||
                stored.represented_pair_rate !=
                    expected.represented_pair_rate ||
                stored.peak_volumetric_pair_rate !=
                    expected.peak_volumetric_pair_rate ||
                stored.start_step != expected.start_step ||
                stored.end_step != expected.end_step ||
                stored.x_min != expected.x_min ||
                stored.x_max != expected.x_max ||
                stored.y_min != expected.y_min ||
                stored.y_max != expected.y_max ||
                stored.first_drift.x != expected.first_drift.x ||
                stored.first_drift.y != expected.first_drift.y ||
                stored.first_drift.z != expected.first_drift.z ||
                stored.second_drift.x != expected.second_drift.x ||
                stored.second_drift.y != expected.second_drift.y ||
                stored.second_drift.z != expected.second_drift.z ||
                stored.first_thermal_velocity !=
                    expected.first_thermal_velocity ||
                stored.second_thermal_velocity !=
                    expected.second_thermal_velocity ||
                stored.spatial_profile.density_profile !=
                    expected.spatial_profile.density_profile ||
                stored.spatial_profile.profile_center_x !=
                    expected.spatial_profile.profile_center_x ||
                stored.spatial_profile.profile_center_y !=
                    expected.spatial_profile.profile_center_y ||
                stored.spatial_profile.profile_scale_x !=
                    expected.spatial_profile.profile_scale_x ||
                stored.spatial_profile.profile_scale_y !=
                    expected.spatial_profile.profile_scale_y ||
                stored.spatial_profile.profile_amplitude !=
                    expected.spatial_profile.profile_amplitude ||
                stored.spatial_profile.profile_phase !=
                    expected.spatial_profile.profile_phase ||
                stored.spatial_profile.profile_mode_x !=
                    expected.spatial_profile.profile_mode_x ||
                stored.spatial_profile.profile_mode_y !=
                    expected.spatial_profile.profile_mode_y ||
                stored.spatial_profile
                        .max_profile_sampling_attempts !=
                    expected.spatial_profile
                        .max_profile_sampling_attempts) {
                throw std::runtime_error(
                    "checkpoint source metadata does not match 2D config");
            }
        }
    } else {
        for (auto& source : source_diagnostics_) {
            source.macro_pairs_created = 0;
            source.represented_pairs_created = 0.0;
            source.fractional_macro_pair_remainder = 0.0;
            source.injected_kinetic_energy = 0.0;
        }
    }

    for (std::size_t species_id = 0; species_id < species_.size(); ++species_id) {
        std::size_t stored_species_id = 0;
        std::string stored_name;
        std::size_t particle_count = 0;
        in >> key >> stored_species_id >> stored_name >> particle_count;
        if (key != "species" || stored_species_id != species_id || stored_name != species_[species_id].name()) {
            throw std::runtime_error("checkpoint species metadata does not match 2D config");
        }
        auto& particles = species_[species_id].particles();
        particles.resize(particle_count);
        for (auto& p : particles) {
            int alive = 0;
            in >> p.position.x >> p.position.y
               >> p.velocity.x >> p.velocity.y;
            if (checkpoint_v3 || checkpoint_v4 ||
                checkpoint_v5 || checkpoint_v6) {
                in >> p.velocity_z;
            } else {
                p.velocity_z = 0.0;
            }
            in >> p.velocity_half.x >> p.velocity_half.y;
            if (checkpoint_v3 || checkpoint_v4 ||
                checkpoint_v5 || checkpoint_v6) {
                in >> p.velocity_half_z;
            } else {
                p.velocity_half_z = 0.0;
            }
            in >> alive;
            p.alive = alive != 0;
        }
    }
    require_stream(in, "failed while reading 2D checkpoint: " + path.string());
    deposit_and_solve();
    initialized_ = true;
}

RunSummary2D Simulation2D::run() {
    if (!cfg_.restart_path.empty()) load_checkpoint(cfg_.restart_path);
    else initialize();

    write_unit_metadata(
        cfg_.output_dir, cfg_.units, 2,
        cfg_.out_of_plane_depth);
    if (!cfg_.initial_state_path.empty()) {
        write_external_particle_state_metadata(
            cfg_.output_dir /
                "initial_state_metadata.txt",
            cfg_.initial_state_path,
            initial_state_metadata_,
            cfg_.initial_state_signature);
    }
    std::vector<InitializationSpeciesMoments> initialization_moments;
    initialization_moments.reserve(species_.size());
    for (const auto& species : species_) {
        initialization_moments.push_back(
            summarize_initialization(species));
    }
    write_initialization_report(
        cfg_.output_dir / "initialization.csv", 2,
        !cfg_.restart_path.empty()
            ? "restart"
            : (!cfg_.initial_state_path.empty()
                   ? "external"
                   : "generated"),
        initialization_moments);
    const auto initialization_acceptance =
        assess_initialization_acceptance(
            cfg_.initialization_acceptance,
            initialization_moments, 3);
    write_initialization_acceptance_report(
        cfg_.output_dir / "initialization_acceptance.csv",
        initialization_acceptance);
    enforce_initialization_acceptance(
        initialization_acceptance);
    Diagnostics2D diag(
        cfg_.output_dir, species_, cfg_.units.permittivity(),
        cfg_.out_of_plane_depth);
    std::ofstream source_output;
    if (!sources_.empty()) {
        source_output.open(cfg_.output_dir / "sources.csv");
        if (!source_output) {
            throw std::runtime_error(
                "cannot open 2D source diagnostics output");
        }
        source_output
            << "step,time,source,macro_pairs_created,"
               "represented_pairs_created,"
               "fractional_macro_pair_remainder,"
               "injected_kinetic_energy,"
               "effective_profile_area,out_of_plane_depth,"
               "configured_peak_volumetric_pair_rate,"
               "resolved_represented_pair_rate\n";
    }
    const auto write_source_sample = [&]() {
        for (std::size_t source_id = 0;
             source_id < sources_.size(); ++source_id) {
            const auto& source = sources_[source_id];
            const auto& diagnostics =
                source_diagnostics_[source_id];
            const double rate =
                (source.config.represented_pair_rate ||
                 source.config.peak_volumetric_pair_rate)
                    ? source.represented_pair_rate
                    :
                    static_cast<double>(
                        source.config.pairs_per_step) *
                    species_[source.first_species].weight() /
                    cfg_.dt;
            source_output << step_ << ',' << std::setprecision(17)
                          << time_ << ',' << diagnostics.name
                          << ',' << diagnostics.macro_pairs_created
                          << ',' << diagnostics.represented_pairs_created
                          << ','
                          << diagnostics
                                 .fractional_macro_pair_remainder
                          << ','
                          << diagnostics.injected_kinetic_energy
                          << ',' << source.effective_profile_area
                          << ',' << cfg_.out_of_plane_depth
                          << ','
                          << source.config
                                 .peak_volumetric_pair_rate
                                 .value_or(0.0)
                          << ',' << rate << '\n';
        }
        if (source_output) source_output.flush();
    };
    diag.write_header();
    auto s0 = diag.sample(step_, time_, mesh_, species_, boundary_losses_);
    diag.write_sample(s0);
    write_source_sample();
    if (cfg_.vtk_output) write_vtk_outputs(mesh_, cfg_.output_dir, step_, cfg_.vtk_format);
    if (cfg_.particle_output) {
        diag.write_particle_sample(step_, species_, cfg_.particle_output_stride, cfg_.particle_sample_count);
    }
    if (cfg_.checkpoint_output) save_checkpoint(checkpoint_path_for_step(cfg_, step_));

    const std::size_t particle_interval = cfg_.particle_output_interval == 0 ? cfg_.output_interval : cfg_.particle_output_interval;
    const std::size_t limit = cfg_.mode == RunMode::SteadyState ? cfg_.max_steps : cfg_.steps;
    RunSummary2D summary;
    summary.final_sample = s0;
    while (step_ < limit) {
        step();
        bool reached_steady = false;
        if (step_ % cfg_.output_interval == 0 || step_ == limit) {
            auto s = diag.sample(step_, time_, mesh_, species_, boundary_losses_);
            diag.write_sample(s);
            write_source_sample();
            summary.final_sample = s;
            if (cfg_.vtk_output) {
                write_vtk_outputs(mesh_, cfg_.output_dir, step_, cfg_.vtk_format);
            }
            reached_steady = cfg_.mode == RunMode::SteadyState &&
                             adjacent_energy_windows_converged(diag.history(), cfg_.steady_window, cfg_.steady_tolerance);
            if (reached_steady) summary.steady_state_reached = true;
        }
        if (cfg_.particle_output && (step_ % particle_interval == 0 || step_ == limit || reached_steady)) {
            diag.write_particle_sample(step_, species_, cfg_.particle_output_stride, cfg_.particle_sample_count);
        }
        if (cfg_.checkpoint_output &&
            (step_ % cfg_.checkpoint_interval == 0 || step_ == limit || reached_steady)) {
            save_checkpoint(checkpoint_path_for_step(cfg_, step_));
        }
        if (reached_steady) break;
    }
    summary.steps_completed = step_;
    summary.final_time = time_;
    if (summary.final_sample.step != step_) summary.final_sample = sample();
    return summary;
}

DiagnosticSample2D Simulation2D::sample() const {
    DiagnosticSample2D s;
    s.step = step_;
    s.time = time_;
    s.boundary_losses = boundary_losses_;
    s.live_particles_by_species.reserve(species_.size());
    for (const auto& sp : species_) {
        s.kinetic_energy += sp.kinetic_energy();
        const auto live = sp.live_count();
        s.live_particles += live;
        s.live_particles_by_species.push_back(live);
    }
    for (std::size_t j = 0; j < mesh_.ny(); ++j) {
        for (std::size_t i = 0; i < mesh_.nx(); ++i) {
            const auto idx = mesh_.index(i, j);
            const double e2 = mesh_.electric_x()[idx] * mesh_.electric_x()[idx] + mesh_.electric_y()[idx] * mesh_.electric_y()[idx];
            const double volume =
                mesh_.node_area(i, j) *
                cfg_.out_of_plane_depth;
            s.field_energy +=
                0.5 * cfg_.units.permittivity() * e2 * volume;
            s.charge_l1 += std::abs(mesh_.rho()[idx]) * volume;
        }
    }
    s.total_energy = s.kinetic_energy + s.field_energy;
    return s;
}
}
