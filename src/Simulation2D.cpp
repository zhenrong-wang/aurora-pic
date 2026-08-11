#include "pic/Simulation2D.hpp"
#include "pic/Convergence.hpp"
#include "pic/ParticleState.hpp"
#include "pic/Pusher.hpp"
#include "pic/Runtime.hpp"
#include "pic/Units.hpp"
#include "pic/VTKWriter.hpp"
#include <algorithm>
#include <array>
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
constexpr const char* kCheckpointMagicV7 = "AuroraPIC-checkpoint-v7";
constexpr const char* kCheckpointMagicV8 = "AuroraPIC-checkpoint-v8";
constexpr const char* kCheckpointMagicV9 = "AuroraPIC-checkpoint-v9";
constexpr const char* kCheckpointMagicV10 = "AuroraPIC-checkpoint-v10";

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

std::size_t boundary_loss(
    const BoundaryLoss2D& losses, BoundarySide2DName side) {
    switch (side) {
        case BoundarySide2DName::Left:
            return losses.absorbed_left;
        case BoundarySide2DName::Right:
            return losses.absorbed_right;
        case BoundarySide2DName::Bottom:
            return losses.absorbed_bottom;
        case BoundarySide2DName::Top:
            return losses.absorbed_top;
    }
    throw std::logic_error("invalid 2D boundary side");
}

BoundarySide2DName boundary_side_from_string(
    const std::string& value) {
    if (value == "left") return BoundarySide2DName::Left;
    if (value == "right") return BoundarySide2DName::Right;
    if (value == "bottom") return BoundarySide2DName::Bottom;
    if (value == "top") return BoundarySide2DName::Top;
    throw std::runtime_error(
        "invalid 2D checkpoint boundary side '" + value + "'");
}

double line_average(
    const Mesh2D& mesh, CoordinateAxis axis, double coordinate) {
    const bool along_x = axis == CoordinateAxis::X;
    const double spacing = along_x ? mesh.dx() : mesh.dy();
    const std::size_t extent = along_x ? mesh.nx() : mesh.ny();
    const Boundary axis_boundary =
        along_x ? mesh.boundary_x() : mesh.boundary_y();
    const double axis_length =
        along_x ? mesh.length_x() : mesh.length_y();
    const double normalized_coordinate =
        axis_boundary == Boundary::Periodic
            ? wrap_periodic(coordinate, axis_length)
            : coordinate;
    const double grid_coordinate =
        normalized_coordinate / spacing;
    std::size_t lower = static_cast<std::size_t>(
        std::floor(grid_coordinate));
    double fraction =
        grid_coordinate - static_cast<double>(lower);
    std::size_t upper = 0;
    if (axis_boundary == Boundary::Periodic) {
        lower %= extent;
        upper = (lower + 1) % extent;
    } else {
        lower = std::min(lower, extent - 2);
        fraction = std::clamp(
            grid_coordinate - static_cast<double>(lower),
            0.0, 1.0);
        upper = lower + 1;
    }
    double weighted_sum = 0.0;
    double weight_sum = 0.0;
    const std::size_t transverse_extent =
        along_x ? mesh.ny() : mesh.nx();
    const Boundary transverse_boundary =
        along_x ? mesh.boundary_y() : mesh.boundary_x();
    for (std::size_t transverse = 0;
         transverse < transverse_extent; ++transverse) {
        const std::size_t first =
            along_x
                ? mesh.index(lower, transverse)
                : mesh.index(transverse, lower);
        const std::size_t second =
            along_x
                ? mesh.index(upper, transverse)
                : mesh.index(transverse, upper);
        const double value =
            (1.0 - fraction) * mesh.phi()[first] +
            fraction * mesh.phi()[second];
        const double weight =
            transverse_boundary == Boundary::Periodic ||
                    (transverse != 0 &&
                     transverse + 1 != transverse_extent)
                ? 1.0
                : 0.5;
        weighted_sum += weight * value;
        weight_sum += weight;
    }
    return weighted_sum / weight_sum;
}
} // namespace

std::string to_string(BoundarySide2DName side) {
    switch (side) {
        case BoundarySide2DName::Left: return "left";
        case BoundarySide2DName::Right: return "right";
        case BoundarySide2DName::Bottom: return "bottom";
        case BoundarySide2DName::Top: return "top";
    }
    return "unknown";
}

std::string to_string(CurrentSourceControlMode mode) {
    switch (mode) {
        case CurrentSourceControlMode::Cumulative:
            return "cumulative";
        case CurrentSourceControlMode::TimestepLocal:
            return "timestep_local";
    }
    return "unknown";
}

std::string to_string(PotentialReferenceCorrection correction) {
    switch (correction) {
        case PotentialReferenceCorrection::Gauge:
            return "gauge";
        case PotentialReferenceCorrection::Affine:
            return "affine";
    }
    return "unknown";
}

Simulation2D::Simulation2D(Simulation2DConfig cfg)
    : cfg_(std::move(cfg)),
      mesh_(cfg_.nx, cfg_.ny, cfg_.length_x, cfg_.length_y,
            cfg_.boundary_x.value_or(cfg_.boundary),
            cfg_.boundary_y.value_or(cfg_.boundary),
            cfg_.boundary_config),
      solver_(cfg_.units.permittivity()),
      rng_(cfg_.seed) {
    if (cfg_.resolved_diagnostics.enabled &&
        cfg_.resolved_diagnostics.interval == 0) {
        cfg_.resolved_diagnostics.interval = cfg_.output_interval;
    }
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
    if (cfg_.resolved_diagnostics.enabled &&
        cfg_.resolved_diagnostics.interval == 0) {
        throw std::invalid_argument(
            "2D resolved diagnostic interval must be positive");
    }
    if (cfg_.resolved_diagnostics.enabled) {
        const auto profile_axis =
            cfg_.resolved_diagnostics.profile_axis;
        const auto mode_axis =
            cfg_.resolved_diagnostics.mode_axis;
        if (profile_axis == CoordinateAxis::Z ||
            mode_axis == CoordinateAxis::Z ||
            profile_axis == mode_axis) {
            throw std::invalid_argument(
                "2D resolved diagnostics require distinct x/y profile and mode axes");
        }
        const Boundary mode_boundary =
            mode_axis == CoordinateAxis::X
                ? mesh_.boundary_x()
                : mesh_.boundary_y();
        const std::size_t mode_nodes =
            mode_axis == CoordinateAxis::X
                ? mesh_.nx()
                : mesh_.ny();
        if (mode_boundary != Boundary::Periodic) {
            throw std::invalid_argument(
                "2D resolved diagnostic mode axis must be periodic");
        }
        if (cfg_.resolved_diagnostics.max_mode >
            mode_nodes / 2) {
            throw std::invalid_argument(
                "2D resolved max_mode exceeds the mode-axis Nyquist limit");
        }
    }
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
    species_boundary_losses_.resize(species_.size());
    if (cfg_.current_regulated_source) {
        const auto& source = *cfg_.current_regulated_source;
        const auto found = std::find_if(
            species_.begin(), species_.end(),
            [&](const Species2D& species) {
                return species.name() == source.species;
            });
        if (found == species_.end()) {
            throw std::invalid_argument(
                "2D current-regulated source references unknown species '" +
                source.species + "'");
        }
        current_regulated_species_ = static_cast<std::size_t>(
            std::distance(species_.begin(), found));
        if (found->charge() == 0.0) {
            throw std::invalid_argument(
                "2D current-regulated source species must be charged");
        }
        const auto boundary_policy =
            [&](BoundarySide2DName side) {
                switch (side) {
                    case BoundarySide2DName::Left:
                        return cfg_.particle_boundary_config.left;
                    case BoundarySide2DName::Right:
                        return cfg_.particle_boundary_config.right;
                    case BoundarySide2DName::Bottom:
                        return cfg_.particle_boundary_config.bottom;
                    case BoundarySide2DName::Top:
                        return cfg_.particle_boundary_config.top;
                }
                return ParticleBoundary::Auto;
            };
        if (boundary_policy(source.monitor_boundary) !=
            ParticleBoundary::Absorbing) {
            throw std::invalid_argument(
                "2D current-regulated source monitor boundary must be absorbing");
        }
        const double normal_length =
            source.emission_boundary == BoundarySide2DName::Left ||
                    source.emission_boundary == BoundarySide2DName::Right
                ? mesh_.length_x()
                : mesh_.length_y();
        if (!std::isfinite(source.emission_inset) ||
            source.emission_inset < 0.0 ||
            !(source.emission_inset < normal_length)) {
            throw std::invalid_argument(
                "2D current-regulated source emission_inset must be finite and inside the domain");
        }
        if (!std::isfinite(source.thermal_velocity) ||
            source.thermal_velocity < 0.0 ||
            !std::isfinite(source.drift.x) ||
            !std::isfinite(source.drift.y) ||
            !std::isfinite(source.drift.z)) {
            throw std::invalid_argument(
                "2D current-regulated source velocities must be finite and thermal_velocity non-negative");
        }
        current_regulated_source_diagnostics_.emplace();
    }
    if (cfg_.potential_reference) {
        const auto& reference = *cfg_.potential_reference;
        if (reference.axis == CoordinateAxis::Z ||
            !std::isfinite(reference.coordinate) ||
            !std::isfinite(reference.target)) {
            throw std::invalid_argument(
                "2D potential reference requires a finite x/y coordinate and target");
        }
        const double length =
            reference.axis == CoordinateAxis::X
                ? mesh_.length_x()
                : mesh_.length_y();
        if (reference.coordinate < 0.0 ||
            reference.coordinate > length) {
            throw std::invalid_argument(
                "2D potential reference coordinate lies outside the domain");
        }
        if (reference.correction ==
                PotentialReferenceCorrection::Affine &&
            reference.coordinate == 0.0) {
            throw std::invalid_argument(
                "2D affine potential reference coordinate must be positive");
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
    std::fill(
        species_boundary_losses_.begin(),
        species_boundary_losses_.end(), BoundaryLoss2D{});
    potential_reference_offset_ = 0.0;
    if (current_regulated_source_diagnostics_) {
        *current_regulated_source_diagnostics_ = {};
    }
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
                cfg_.initial_state_path, 2, 3,
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
    apply_potential_reference();
}

void Simulation2D::apply_potential_reference() {
    potential_reference_offset_ = 0.0;
    if (!cfg_.potential_reference) return;
    const auto& reference = *cfg_.potential_reference;
    const double mean =
        line_average(mesh_, reference.axis, reference.coordinate);
    potential_reference_offset_ = mean - reference.target;
    if (!std::isfinite(potential_reference_offset_)) {
        throw std::runtime_error(
            "2D potential-reference offset is not finite");
    }
    if (reference.correction ==
        PotentialReferenceCorrection::Gauge) {
        for (double& potential : mesh_.phi()) {
            potential -= potential_reference_offset_;
        }
        return;
    }
    const double field_correction =
        potential_reference_offset_ / reference.coordinate;
    for (std::size_t j = 0; j < mesh_.ny(); ++j) {
        for (std::size_t i = 0; i < mesh_.nx(); ++i) {
            const double coordinate =
                reference.axis == CoordinateAxis::X
                    ? static_cast<double>(i) * mesh_.dx()
                    : static_cast<double>(j) * mesh_.dy();
            mesh_.phi()[mesh_.index(i, j)] -=
                coordinate * field_correction;
        }
    }
    auto& corrected_field =
        reference.axis == CoordinateAxis::X
            ? mesh_.electric_x()
            : mesh_.electric_y();
    for (double& field : corrected_field) {
        field += field_correction;
    }
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

void Simulation2D::inject_current_regulated_source() {
    if (!cfg_.current_regulated_source ||
        !current_regulated_species_ ||
        !current_regulated_source_diagnostics_) {
        return;
    }
    const auto& config = *cfg_.current_regulated_source;
    auto& emitted_species = species_[*current_regulated_species_];
    auto& diagnostics =
        *current_regulated_source_diagnostics_;
    double cumulative_charge = 0.0;
    double cumulative_negative_charge = 0.0;
    double cumulative_positive_charge = 0.0;
    for (std::size_t species_id = 0;
         species_id < species_.size(); ++species_id) {
        const double contribution =
            static_cast<double>(boundary_loss(
                species_boundary_losses_[species_id],
                config.monitor_boundary)) *
            species_[species_id].charge() *
            species_[species_id].weight();
        cumulative_charge += contribution;
        if (contribution < 0.0) {
            cumulative_negative_charge += contribution;
        } else {
            cumulative_positive_charge += contribution;
        }
    }
    if (!std::isfinite(cumulative_charge)) {
        throw std::runtime_error(
            "2D current-regulated source monitored charge overflow");
    }
    const double delta_charge =
        cumulative_charge -
        diagnostics.processed_monitored_charge;
    const double delta_negative_charge =
        cumulative_negative_charge -
        diagnostics.processed_monitored_negative_charge;
    const double delta_positive_charge =
        cumulative_positive_charge -
        diagnostics.processed_monitored_positive_charge;
    const double macro_charge =
        emitted_species.charge() * emitted_species.weight();
    const double requested =
        diagnostics.control_macro_remainder +
        delta_charge / macro_charge;
    if (!std::isfinite(requested)) {
        throw std::runtime_error(
            "2D current-regulated source accumulator overflow");
    }
    if (
        diagnostics.control_updates ==
        std::numeric_limits<std::size_t>::max()
    ) {
        throw std::runtime_error(
            "2D current-regulated source update counter overflow");
    }
    ++diagnostics.control_updates;
    if (
        config.control_mode ==
            CurrentSourceControlMode::TimestepLocal &&
        requested < 0.0
    ) {
        const double reverse_demand = -requested;
        const double squared_reverse_demand =
            reverse_demand * reverse_demand;
        if (
            diagnostics.reverse_demand_steps ==
                std::numeric_limits<std::size_t>::max() ||
            !std::isfinite(
                diagnostics.cumulative_reverse_demand_macroparticles +
                reverse_demand) ||
            !std::isfinite(
                diagnostics.squared_reverse_demand_macroparticles +
                squared_reverse_demand) ||
            !std::isfinite(
                diagnostics.distributed_reverse_demand_macroparticles +
                reverse_demand) ||
            !std::isfinite(
                diagnostics.reverse_monitored_negative_charge +
                delta_negative_charge) ||
            !std::isfinite(
                diagnostics.reverse_monitored_positive_charge +
                delta_positive_charge)
        ) {
            throw std::runtime_error(
                "2D current-regulated source reverse-demand "
                "diagnostic overflow");
        }
        ++diagnostics.reverse_demand_steps;
        diagnostics.cumulative_reverse_demand_macroparticles +=
            reverse_demand;
        diagnostics.maximum_reverse_demand_macroparticles =
            std::max(
                diagnostics.maximum_reverse_demand_macroparticles,
                reverse_demand);
        diagnostics.squared_reverse_demand_macroparticles +=
            squared_reverse_demand;
        diagnostics.distributed_reverse_demand_macroparticles +=
            reverse_demand;
        diagnostics.reverse_monitored_negative_charge +=
            delta_negative_charge;
        diagnostics.reverse_monitored_positive_charge +=
            delta_positive_charge;
        auto* distribution_bin =
            reverse_demand <= 1.5
                ? &diagnostics.reverse_one_macro_steps
                : reverse_demand <= 2.5
                    ? &diagnostics.reverse_two_macro_steps
                    : &diagnostics.reverse_multi_macro_steps;
        if (*distribution_bin ==
            std::numeric_limits<std::size_t>::max()) {
            throw std::runtime_error(
                "2D current-regulated source reverse-demand "
                "distribution counter overflow");
        }
        ++(*distribution_bin);
    }
    const double accumulated =
        config.control_mode ==
                CurrentSourceControlMode::TimestepLocal
            ? std::max(0.0, requested)
            : requested;
    if (accumulated >=
        static_cast<double>(
            std::numeric_limits<std::size_t>::max())) {
        throw std::runtime_error(
            "2D current-regulated source per-step count overflow");
    }
    const std::size_t particles_this_step =
        accumulated > 0.0
            ? static_cast<std::size_t>(std::floor(accumulated))
            : 0;
    const auto dead_slots = static_cast<std::size_t>(
        std::count_if(
            emitted_species.particles().begin(),
            emitted_species.particles().end(),
            [](const Particle2D& particle) {
                return !particle.alive;
            }));
    const std::size_t growth_capacity =
        cfg_.max_particles_per_species -
        std::min(
            cfg_.max_particles_per_species,
            emitted_species.particles().size());
    if (particles_this_step > dead_slots &&
        particles_this_step - dead_slots > growth_capacity) {
        throw std::runtime_error(
            "2D current-regulated source would exceed max_particles_per_species");
    }
    if (particles_this_step >
        std::numeric_limits<std::size_t>::max() -
            diagnostics.macro_particles_created) {
        throw std::runtime_error(
            "2D current-regulated source counter overflow");
    }

    std::uniform_real_distribution<double> transverse(0.0, 1.0);
    std::normal_distribution<double> velocity_x(
        config.drift.x, config.thermal_velocity);
    std::normal_distribution<double> velocity_y(
        config.drift.y, config.thermal_velocity);
    std::normal_distribution<double> velocity_z(
        config.drift.z, config.thermal_velocity);
    double injected_energy = 0.0;
    for (std::size_t particle_id = 0;
         particle_id < particles_this_step; ++particle_id) {
        Particle2D particle;
        const double sample = transverse(rng_);
        switch (config.emission_boundary) {
            case BoundarySide2DName::Left:
                particle.position = {
                    config.emission_inset,
                    sample * mesh_.length_y()};
                break;
            case BoundarySide2DName::Right:
                particle.position = {
                    mesh_.length_x() - config.emission_inset,
                    sample * mesh_.length_y()};
                break;
            case BoundarySide2DName::Bottom:
                particle.position = {
                    sample * mesh_.length_x(),
                    config.emission_inset};
                break;
            case BoundarySide2DName::Top:
                particle.position = {
                    sample * mesh_.length_x(),
                    mesh_.length_y() - config.emission_inset};
                break;
        }
        particle.velocity = {
            velocity_x(rng_), velocity_y(rng_)};
        particle.velocity_z = velocity_z(rng_);
        particle.velocity_half = particle.velocity;
        particle.velocity_half_z = particle.velocity_z;
        const double speed_squared =
            particle.velocity.x * particle.velocity.x +
            particle.velocity.y * particle.velocity.y +
            particle.velocity_z * particle.velocity_z;
        injected_energy +=
            0.5 * emitted_species.weight() *
            emitted_species.mass() * speed_squared;
        initialize_particle_pusher(
            particle,
            interpolate_electric(mesh_, particle.position),
            emitted_species.charge() / emitted_species.mass(),
            cfg_);
        const auto dead = std::find_if(
            emitted_species.particles().begin(),
            emitted_species.particles().end(),
            [](const Particle2D& candidate) {
                return !candidate.alive;
            });
        if (dead == emitted_species.particles().end()) {
            emitted_species.particles().push_back(particle);
        } else {
            *dead = particle;
        }
    }
    const double represented_increment =
        static_cast<double>(particles_this_step) *
        emitted_species.weight();
    if (!std::isfinite(injected_energy) ||
        !std::isfinite(
            diagnostics.injected_kinetic_energy +
            injected_energy) ||
        !std::isfinite(
            diagnostics.represented_particles_created +
            represented_increment)) {
        throw std::runtime_error(
            "2D current-regulated source diagnostic overflow");
    }
    diagnostics.macro_particles_created += particles_this_step;
    diagnostics.represented_particles_created +=
        represented_increment;
    diagnostics.control_macro_remainder =
        accumulated - static_cast<double>(particles_this_step);
    diagnostics.processed_monitored_charge =
        cumulative_charge;
    diagnostics.processed_monitored_negative_charge =
        cumulative_negative_charge;
    diagnostics.processed_monitored_positive_charge =
        cumulative_positive_charge;
    diagnostics.injected_kinetic_energy += injected_energy;
}

void Simulation2D::apply_particle_boundaries(
    Particle2D& particle, std::size_t species_id) {
    if (!apply_lower_boundary(particle.position.x, particle.velocity_half.x, mesh_.length_x(),
                              cfg_.particle_boundary_config.left, boundary_losses_.absorbed_left)) {
        ++species_boundary_losses_.at(species_id).absorbed_left;
        particle.alive = false;
        return;
    }
    if (!apply_upper_boundary(particle.position.x, particle.velocity_half.x, mesh_.length_x(),
                              cfg_.particle_boundary_config.right, boundary_losses_.absorbed_right)) {
        ++species_boundary_losses_.at(species_id).absorbed_right;
        particle.alive = false;
        return;
    }
    if (!apply_lower_boundary(particle.position.y, particle.velocity_half.y, mesh_.length_y(),
                              cfg_.particle_boundary_config.bottom, boundary_losses_.absorbed_bottom)) {
        ++species_boundary_losses_.at(species_id).absorbed_bottom;
        particle.alive = false;
        return;
    }
    if (!apply_upper_boundary(particle.position.y, particle.velocity_half.y, mesh_.length_y(),
                              cfg_.particle_boundary_config.top, boundary_losses_.absorbed_top)) {
        ++species_boundary_losses_.at(species_id).absorbed_top;
        particle.alive = false;
        return;
    }
}

void Simulation2D::step() {
    if (!initialized_) initialize();
    inject_current_regulated_source();
    inject_volumetric_pair_sources();
    for (std::size_t species_id = 0;
         species_id < species_.size(); ++species_id) {
        auto& sp = species_[species_id];
        const double qm = sp.charge() / sp.mass();
        for (auto& particle : sp.particles()) {
            if (!particle.alive) continue;
            kick_particle(particle, interpolate_electric(mesh_, particle.position), qm, cfg_);
            drift_leapfrog(particle, cfg_.dt);
            apply_particle_boundaries(particle, species_id);
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
    out << kCheckpointMagicV10 << '\n';
    out << "dimension 2\n";
    out << "units " << to_string(cfg_.units.system) << ' '
        << cfg_.units.relative_permittivity << ' '
        << cfg_.units.permittivity() << ' '
        << cfg_.out_of_plane_depth << "\n";
    out << "step " << step_ << "\n";
    out << "time " << time_ << "\n";
    out << "boundary_losses " << boundary_losses_.absorbed_left << ' ' << boundary_losses_.absorbed_right << ' '
        << boundary_losses_.absorbed_bottom << ' ' << boundary_losses_.absorbed_top << "\n";
    out << "species_boundary_losses "
        << species_boundary_losses_.size() << "\n";
    for (std::size_t species_id = 0;
         species_id < species_boundary_losses_.size();
         ++species_id) {
        const auto& losses = species_boundary_losses_[species_id];
        out << "species_boundary_loss " << species_id << ' '
            << losses.absorbed_left << ' '
            << losses.absorbed_right << ' '
            << losses.absorbed_bottom << ' '
            << losses.absorbed_top << "\n";
    }
    out << "current_regulated_source "
        << (cfg_.current_regulated_source ? 1 : 0);
    if (cfg_.current_regulated_source) {
        const auto& config = *cfg_.current_regulated_source;
        const auto& diagnostics =
            *current_regulated_source_diagnostics_;
        out << ' ' << config.species
            << ' ' << to_string(config.control_mode)
            << ' ' << to_string(config.monitor_boundary)
            << ' ' << to_string(config.emission_boundary)
            << ' ' << config.emission_inset
            << ' ' << config.drift.x
            << ' ' << config.drift.y
            << ' ' << config.drift.z
            << ' ' << config.thermal_velocity
            << ' ' << diagnostics.macro_particles_created
            << ' ' << diagnostics.represented_particles_created
            << ' ' << diagnostics.control_macro_remainder
            << ' ' << diagnostics.processed_monitored_charge
            << ' ' << diagnostics.injected_kinetic_energy
            << ' ' << diagnostics.processed_monitored_negative_charge
            << ' ' << diagnostics.processed_monitored_positive_charge
            << ' ' << diagnostics.control_updates
            << ' ' << diagnostics.reverse_demand_steps
            << ' ' << diagnostics.reverse_diagnostics_start_step
            << ' ' << diagnostics.cumulative_reverse_demand_macroparticles
            << ' ' << diagnostics.maximum_reverse_demand_macroparticles
            << ' ' << diagnostics.reverse_distribution_start_step
            << ' ' << diagnostics.reverse_one_macro_steps
            << ' ' << diagnostics.reverse_two_macro_steps
            << ' ' << diagnostics.reverse_multi_macro_steps
            << ' ' <<
                diagnostics.distributed_reverse_demand_macroparticles
            << ' ' << diagnostics.squared_reverse_demand_macroparticles
            << ' ' << diagnostics.reverse_monitored_negative_charge
            << ' ' << diagnostics.reverse_monitored_positive_charge;
    }
    out << "\n";
    out << "potential_reference "
        << (cfg_.potential_reference ? 1 : 0);
    if (cfg_.potential_reference) {
        out << ' ' << to_string(cfg_.potential_reference->axis)
            << ' ' << cfg_.potential_reference->coordinate
            << ' ' << cfg_.potential_reference->target
            << ' ' << to_string(
                   cfg_.potential_reference->correction);
    }
    out << "\n";
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
    const bool checkpoint_v7 = magic == kCheckpointMagicV7;
    const bool checkpoint_v8 = magic == kCheckpointMagicV8;
    const bool checkpoint_v9 = magic == kCheckpointMagicV9;
    const bool checkpoint_v10 = magic == kCheckpointMagicV10;
    if (!checkpoint_v1 && !checkpoint_v2 && !checkpoint_v3 &&
        !checkpoint_v4 && !checkpoint_v5 && !checkpoint_v6 &&
        !checkpoint_v7 && !checkpoint_v8 && !checkpoint_v9 &&
        !checkpoint_v10) {
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
        if (checkpoint_v6 || checkpoint_v7 || checkpoint_v8 ||
            checkpoint_v9 || checkpoint_v10) {
            in >> out_of_plane_depth;
        }
        if ((!checkpoint_v2 && !checkpoint_v3 && !checkpoint_v4 &&
             !checkpoint_v5 && !checkpoint_v6 && !checkpoint_v7 &&
             !checkpoint_v8 && !checkpoint_v9 && !checkpoint_v10) ||
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
               checkpoint_v7 ||
               checkpoint_v8 ||
               checkpoint_v9 ||
               checkpoint_v10 ||
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
    if (checkpoint_v7 || checkpoint_v8 || checkpoint_v9 ||
        checkpoint_v10) {
        std::size_t loss_species_count = 0;
        in >> key >> loss_species_count;
        if (key != "species_boundary_losses" ||
            loss_species_count != species_.size()) {
            throw std::runtime_error(
                "checkpoint species boundary-loss count does not match 2D config");
        }
        for (std::size_t species_id = 0;
             species_id < loss_species_count; ++species_id) {
            std::size_t stored_species_id = 0;
            auto& losses = species_boundary_losses_[species_id];
            in >> key >> stored_species_id >>
                losses.absorbed_left >> losses.absorbed_right >>
                losses.absorbed_bottom >> losses.absorbed_top;
            if (key != "species_boundary_loss" ||
                stored_species_id != species_id) {
                throw std::runtime_error(
                    "checkpoint species boundary-loss metadata is invalid");
            }
        }
        int has_current_source = 0;
        in >> key >> has_current_source;
        if (key != "current_regulated_source" ||
            (has_current_source != 0 && has_current_source != 1) ||
            (has_current_source != 0) !=
                cfg_.current_regulated_source.has_value()) {
            throw std::runtime_error(
                "checkpoint current-regulated source presence does not match 2D config");
        }
        if (has_current_source != 0) {
            CurrentRegulatedSource2DConfig stored;
            std::string control_mode;
            std::string monitor_boundary;
            std::string emission_boundary;
            auto& diagnostics =
                *current_regulated_source_diagnostics_;
            in >> stored.species;
            if (checkpoint_v8 || checkpoint_v9 || checkpoint_v10) {
                in >> control_mode;
                if (control_mode == "cumulative") {
                    stored.control_mode =
                        CurrentSourceControlMode::Cumulative;
                } else if (control_mode == "timestep_local") {
                    stored.control_mode =
                        CurrentSourceControlMode::TimestepLocal;
                } else {
                    throw std::runtime_error(
                        "invalid current-source control mode in checkpoint");
                }
            }
            in >> monitor_boundary >>
                emission_boundary >> stored.emission_inset >>
                stored.drift.x >> stored.drift.y >>
                stored.drift.z >> stored.thermal_velocity >>
                diagnostics.macro_particles_created >>
                diagnostics.represented_particles_created >>
                diagnostics.control_macro_remainder >>
                diagnostics.processed_monitored_charge >>
                diagnostics.injected_kinetic_energy;
            if (checkpoint_v9 || checkpoint_v10) {
                in >>
                    diagnostics.processed_monitored_negative_charge >>
                    diagnostics.processed_monitored_positive_charge >>
                    diagnostics.control_updates >>
                    diagnostics.reverse_demand_steps >>
                    diagnostics.reverse_diagnostics_start_step >>
                    diagnostics.cumulative_reverse_demand_macroparticles >>
                    diagnostics.maximum_reverse_demand_macroparticles;
                if (checkpoint_v10) {
                    in >>
                        diagnostics.reverse_distribution_start_step >>
                        diagnostics.reverse_one_macro_steps >>
                        diagnostics.reverse_two_macro_steps >>
                        diagnostics.reverse_multi_macro_steps >>
                        diagnostics
                            .distributed_reverse_demand_macroparticles >>
                        diagnostics
                            .squared_reverse_demand_macroparticles >>
                        diagnostics.reverse_monitored_negative_charge >>
                        diagnostics.reverse_monitored_positive_charge;
                } else {
                    diagnostics.reverse_distribution_start_step = step_;
                }
            } else {
                diagnostics.reverse_diagnostics_start_step = step_;
                diagnostics.reverse_distribution_start_step = step_;
                diagnostics.processed_monitored_negative_charge =
                    std::min(
                        0.0,
                        diagnostics.processed_monitored_charge);
                diagnostics.processed_monitored_positive_charge =
                    std::max(
                        0.0,
                        diagnostics.processed_monitored_charge);
            }
            stored.monitor_boundary =
                boundary_side_from_string(monitor_boundary);
            stored.emission_boundary =
                boundary_side_from_string(emission_boundary);
            const auto& expected =
                *cfg_.current_regulated_source;
            if (stored.species != expected.species ||
                stored.control_mode != expected.control_mode ||
                stored.monitor_boundary != expected.monitor_boundary ||
                stored.emission_boundary != expected.emission_boundary ||
                stored.emission_inset != expected.emission_inset ||
                stored.drift.x != expected.drift.x ||
                stored.drift.y != expected.drift.y ||
                stored.drift.z != expected.drift.z ||
                stored.thermal_velocity !=
                    expected.thermal_velocity ||
                !std::isfinite(
                    diagnostics.represented_particles_created) ||
                diagnostics.represented_particles_created < 0.0 ||
                !std::isfinite(
                    diagnostics.control_macro_remainder) ||
                !std::isfinite(
                    diagnostics.processed_monitored_charge) ||
                ((checkpoint_v9 || checkpoint_v10) &&
                 (!std::isfinite(
                      diagnostics.processed_monitored_negative_charge) ||
                  diagnostics.processed_monitored_negative_charge > 0.0 ||
                  !std::isfinite(
                      diagnostics.processed_monitored_positive_charge) ||
                  diagnostics.processed_monitored_positive_charge < 0.0 ||
                  !std::isfinite(
                      diagnostics.processed_monitored_negative_charge +
                      diagnostics.processed_monitored_positive_charge) ||
                  std::abs(
                      diagnostics.processed_monitored_negative_charge +
                      diagnostics.processed_monitored_positive_charge -
                      diagnostics.processed_monitored_charge) >
                      1e-12 * std::max({
                          1e-30,
                          std::abs(
                              diagnostics
                                  .processed_monitored_negative_charge),
                          std::abs(
                              diagnostics
                                  .processed_monitored_positive_charge),
                          std::abs(
                              diagnostics.processed_monitored_charge)}))) ||
                !std::isfinite(
                    diagnostics.injected_kinetic_energy) ||
                diagnostics.injected_kinetic_energy < 0.0 ||
                diagnostics.reverse_demand_steps >
                    diagnostics.control_updates ||
                diagnostics.reverse_diagnostics_start_step > step_ ||
                !std::isfinite(
                    diagnostics.cumulative_reverse_demand_macroparticles) ||
                diagnostics.cumulative_reverse_demand_macroparticles < 0.0 ||
                !std::isfinite(
                    diagnostics.maximum_reverse_demand_macroparticles) ||
                diagnostics.maximum_reverse_demand_macroparticles < 0.0 ||
                diagnostics.maximum_reverse_demand_macroparticles >
                    diagnostics.cumulative_reverse_demand_macroparticles ||
                diagnostics.reverse_distribution_start_step > step_ ||
                diagnostics.reverse_one_macro_steps >
                    diagnostics.reverse_demand_steps ||
                diagnostics.reverse_two_macro_steps >
                    diagnostics.reverse_demand_steps -
                        diagnostics.reverse_one_macro_steps ||
                diagnostics.reverse_multi_macro_steps >
                    diagnostics.reverse_demand_steps -
                        diagnostics.reverse_one_macro_steps -
                        diagnostics.reverse_two_macro_steps ||
                !std::isfinite(
                    diagnostics
                        .distributed_reverse_demand_macroparticles) ||
                diagnostics
                    .distributed_reverse_demand_macroparticles < 0.0 ||
                diagnostics
                    .distributed_reverse_demand_macroparticles >
                        diagnostics
                            .cumulative_reverse_demand_macroparticles ||
                !std::isfinite(
                    diagnostics.squared_reverse_demand_macroparticles) ||
                diagnostics.squared_reverse_demand_macroparticles < 0.0 ||
                !std::isfinite(
                    diagnostics.reverse_monitored_negative_charge) ||
                diagnostics.reverse_monitored_negative_charge > 0.0 ||
                !std::isfinite(
                    diagnostics.reverse_monitored_positive_charge) ||
                diagnostics.reverse_monitored_positive_charge < 0.0) {
                throw std::runtime_error(
                    "checkpoint current-regulated source metadata is invalid or does not match 2D config");
            }
        }
        int has_potential_reference = 0;
        in >> key >> has_potential_reference;
        if (key != "potential_reference" ||
            (has_potential_reference != 0 &&
             has_potential_reference != 1) ||
            (has_potential_reference != 0) !=
                cfg_.potential_reference.has_value()) {
            throw std::runtime_error(
                "checkpoint potential-reference presence does not match 2D config");
        }
        if (has_potential_reference != 0) {
            std::string axis;
            std::string correction;
            PotentialReference2DConfig stored;
            in >> axis >> stored.coordinate >> stored.target;
            stored.axis = parse_coordinate_axis(axis);
            if (checkpoint_v8 || checkpoint_v9 || checkpoint_v10) {
                in >> correction;
                if (correction == "gauge") {
                    stored.correction =
                        PotentialReferenceCorrection::Gauge;
                } else if (correction == "affine") {
                    stored.correction =
                        PotentialReferenceCorrection::Affine;
                } else {
                    throw std::runtime_error(
                        "invalid potential-reference correction in checkpoint");
                }
            }
            const auto& expected = *cfg_.potential_reference;
            if (stored.axis != expected.axis ||
                stored.coordinate != expected.coordinate ||
                stored.target != expected.target ||
                stored.correction != expected.correction) {
                throw std::runtime_error(
                    "checkpoint potential-reference metadata does not match 2D config");
            }
        }
    } else {
        std::fill(
            species_boundary_losses_.begin(),
            species_boundary_losses_.end(),
            BoundaryLoss2D{});
        if (cfg_.current_regulated_source ||
            cfg_.potential_reference) {
            throw std::runtime_error(
                "legacy checkpoint cannot restore 2D current regulation or potential reference state");
        }
    }
    std::size_t species_count = 0;
    in >> key >> species_count;
    if (key != "species_count" || species_count != species_.size()) throw std::runtime_error("checkpoint species count does not match 2D config");
    in >> key;
    if (key != "rng") throw std::runtime_error("checkpoint missing rng state");
    in >> rng_;
    if (checkpoint_v4 || checkpoint_v5 || checkpoint_v6 ||
        checkpoint_v7 || checkpoint_v8 || checkpoint_v9 ||
        checkpoint_v10) {
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
            if (checkpoint_v5 || checkpoint_v6 || checkpoint_v7 ||
                checkpoint_v8 || checkpoint_v9 || checkpoint_v10) {
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
                if (checkpoint_v6 || checkpoint_v7 || checkpoint_v8 ||
                    checkpoint_v9 || checkpoint_v10) {
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
            if (checkpoint_v5 || checkpoint_v6 || checkpoint_v7 ||
                checkpoint_v8 || checkpoint_v9 || checkpoint_v10) {
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
            if (checkpoint_v5 || checkpoint_v6 || checkpoint_v7 ||
                checkpoint_v8 || checkpoint_v9 || checkpoint_v10) {
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
                checkpoint_v5 || checkpoint_v6 ||
                checkpoint_v7 || checkpoint_v8 || checkpoint_v9 ||
                checkpoint_v10) {
                in >> p.velocity_z;
            } else {
                p.velocity_z = 0.0;
            }
            in >> p.velocity_half.x >> p.velocity_half.y;
            if (checkpoint_v3 || checkpoint_v4 ||
                checkpoint_v5 || checkpoint_v6 ||
                checkpoint_v7 || checkpoint_v8 || checkpoint_v9 ||
                checkpoint_v10) {
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
    std::unique_ptr<ResolvedDiagnostics2D> resolved_diagnostics;
    if (cfg_.resolved_diagnostics.enabled) {
        resolved_diagnostics =
            std::make_unique<ResolvedDiagnostics2D>(
                cfg_.output_dir, cfg_.resolved_diagnostics,
                mesh_, species_, cfg_.units.system,
                cfg_.out_of_plane_depth);
    }
    const auto write_resolved_sample = [&]() {
        if (!resolved_diagnostics ||
            step_ < cfg_.resolved_diagnostics.start_step) {
            return;
        }
        const std::size_t offset =
            step_ - cfg_.resolved_diagnostics.start_step;
        if (offset % cfg_.resolved_diagnostics.interval != 0) {
            return;
        }
        (void)resolved_diagnostics->sample(
            step_, time_, mesh_, species_);
    };
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
    std::ofstream current_source_output;
    if (cfg_.current_regulated_source) {
        current_source_output.open(
            cfg_.output_dir / "current_source.csv");
        if (!current_source_output) {
            throw std::runtime_error(
                "cannot open 2D current-source diagnostics output");
        }
        current_source_output
            << "step,time,species,control_mode,monitor_boundary,"
               "emission_boundary,macro_particles_created,"
               "represented_particles_created,"
               "control_macro_remainder,"
               "control_updates,"
               "reverse_diagnostics_start_step,"
               "reverse_demand_steps,"
               "reverse_demand_step_fraction,"
               "cumulative_reverse_demand_macroparticles,"
               "maximum_reverse_demand_macroparticles,"
               "reverse_distribution_start_step,"
               "reverse_distribution_steps,"
               "reverse_one_macro_steps,"
               "reverse_two_macro_steps,"
               "reverse_multi_macro_steps,"
               "distributed_reverse_demand_macroparticles,"
               "mean_reverse_demand_macroparticles,"
               "rms_reverse_demand_macroparticles,"
               "reverse_monitored_negative_charge,"
               "reverse_monitored_positive_charge,"
               "reverse_monitored_net_charge,"
               "cumulative_monitored_negative_charge,"
               "cumulative_monitored_positive_charge,"
               "cumulative_processed_monitored_charge,"
               "cumulative_emitted_charge,"
               "raw_charge_balance_residual,"
               "unserved_reverse_charge,"
               "charge_balance_residual,"
               "injected_kinetic_energy\n";
    }
    const auto write_current_source_sample = [&]() {
        if (!current_source_output.is_open()) return;
        const auto& config = *cfg_.current_regulated_source;
        const auto& diagnostics =
            *current_regulated_source_diagnostics_;
        const auto& emitted =
            species_[*current_regulated_species_];
        const double emitted_charge =
            diagnostics.represented_particles_created *
            emitted.charge();
        const double raw_residual =
            diagnostics.processed_monitored_charge -
            emitted_charge;
        const double control_residual =
            diagnostics.control_macro_remainder *
            emitted.charge() * emitted.weight();
        const std::size_t reverse_distribution_steps =
            diagnostics.reverse_one_macro_steps +
            diagnostics.reverse_two_macro_steps +
            diagnostics.reverse_multi_macro_steps;
        const double mean_reverse_demand =
            reverse_distribution_steps > 0
                ? diagnostics
                    .distributed_reverse_demand_macroparticles /
                    static_cast<double>(
                        reverse_distribution_steps)
                : 0.0;
        const double rms_reverse_demand =
            reverse_distribution_steps > 0
                ? std::sqrt(
                    diagnostics
                        .squared_reverse_demand_macroparticles /
                    static_cast<double>(reverse_distribution_steps))
                : 0.0;
        current_source_output
            << step_ << ',' << std::setprecision(17) << time_
            << ',' << config.species
            << ',' << to_string(config.control_mode)
            << ',' << to_string(config.monitor_boundary)
            << ',' << to_string(config.emission_boundary)
            << ',' << diagnostics.macro_particles_created
            << ',' << diagnostics.represented_particles_created
            << ',' << diagnostics.control_macro_remainder
            << ',' << diagnostics.control_updates
            << ',' << diagnostics.reverse_diagnostics_start_step
            << ',' << diagnostics.reverse_demand_steps
            << ',' << (
                diagnostics.control_updates > 0
                    ? static_cast<double>(
                          diagnostics.reverse_demand_steps) /
                          static_cast<double>(
                              diagnostics.control_updates)
                    : 0.0)
            << ',' <<
                diagnostics.cumulative_reverse_demand_macroparticles
            << ',' <<
                diagnostics.maximum_reverse_demand_macroparticles
            << ',' << diagnostics.reverse_distribution_start_step
            << ',' << reverse_distribution_steps
            << ',' << diagnostics.reverse_one_macro_steps
            << ',' << diagnostics.reverse_two_macro_steps
            << ',' << diagnostics.reverse_multi_macro_steps
            << ',' <<
                diagnostics.distributed_reverse_demand_macroparticles
            << ',' << mean_reverse_demand
            << ',' << rms_reverse_demand
            << ',' << diagnostics.reverse_monitored_negative_charge
            << ',' << diagnostics.reverse_monitored_positive_charge
            << ',' << (
                diagnostics.reverse_monitored_negative_charge +
                diagnostics.reverse_monitored_positive_charge)
            << ',' <<
                diagnostics.processed_monitored_negative_charge
            << ',' <<
                diagnostics.processed_monitored_positive_charge
            << ',' << diagnostics.processed_monitored_charge
            << ',' << emitted_charge
            << ',' << raw_residual
            << ',' << raw_residual - control_residual
            << ',' << control_residual
            << ',' << diagnostics.injected_kinetic_energy
            << '\n';
        current_source_output.flush();
    };
    std::ofstream boundary_flux_output(
        cfg_.output_dir / "boundary_flux.csv");
    if (!boundary_flux_output) {
        throw std::runtime_error(
            "cannot open 2D boundary-flux diagnostics output");
    }
    boundary_flux_output
        << "step,time,window_start_step,window_start_time,"
           "window_duration,species,boundary,"
           "absorbed_macroparticles,"
           "cumulative_absorbed_macroparticles,"
           "represented_particles,represented_charge,"
           "represented_particle_rate,charge_rate\n";
    auto previous_species_boundary_losses =
        species_boundary_losses_;
    std::size_t boundary_flux_start_step = step_;
    double boundary_flux_start_time = time_;
    const auto write_boundary_flux_sample = [&]() {
        const double duration = time_ - boundary_flux_start_time;
        if (!std::isfinite(duration) || duration < 0.0) {
            throw std::runtime_error(
                "2D boundary-flux diagnostic time is invalid");
        }
        constexpr std::array sides{
            BoundarySide2DName::Left,
            BoundarySide2DName::Right,
            BoundarySide2DName::Bottom,
            BoundarySide2DName::Top,
        };
        for (std::size_t species_id = 0;
             species_id < species_.size(); ++species_id) {
            const auto& species = species_[species_id];
            for (const auto side : sides) {
                const std::size_t cumulative = boundary_loss(
                    species_boundary_losses_[species_id], side);
                const std::size_t previous = boundary_loss(
                    previous_species_boundary_losses[species_id],
                    side);
                if (cumulative < previous) {
                    throw std::runtime_error(
                        "2D boundary-flux counter moved backward");
                }
                const std::size_t absorbed = cumulative - previous;
                const double represented =
                    static_cast<double>(absorbed) * species.weight();
                const double charge = represented * species.charge();
                const double represented_rate =
                    duration > 0.0 ? represented / duration : 0.0;
                const double charge_rate =
                    duration > 0.0 ? charge / duration : 0.0;
                if (!std::isfinite(represented) ||
                    !std::isfinite(charge) ||
                    !std::isfinite(represented_rate) ||
                    !std::isfinite(charge_rate)) {
                    throw std::runtime_error(
                        "2D boundary-flux diagnostic overflow");
                }
                boundary_flux_output
                    << step_ << ',' << std::setprecision(17) << time_
                    << ',' << boundary_flux_start_step
                    << ',' << boundary_flux_start_time
                    << ',' << duration
                    << ',' << species.name()
                    << ',' << to_string(side)
                    << ',' << absorbed
                    << ',' << cumulative
                    << ',' << represented
                    << ',' << charge
                    << ',' << represented_rate
                    << ',' << charge_rate
                    << '\n';
            }
        }
        previous_species_boundary_losses = species_boundary_losses_;
        boundary_flux_start_step = step_;
        boundary_flux_start_time = time_;
        boundary_flux_output.flush();
    };
    std::ofstream potential_reference_output;
    if (cfg_.potential_reference) {
        potential_reference_output.open(
            cfg_.output_dir / "potential_reference.csv");
        if (!potential_reference_output) {
            throw std::runtime_error(
                "cannot open 2D potential-reference diagnostics output");
        }
        potential_reference_output
            << "step,time,axis,correction,coordinate,target,"
               "unshifted_line_mean,applied_offset,"
               "corrected_line_mean\n";
    }
    const auto write_potential_reference_sample = [&]() {
        if (!potential_reference_output.is_open()) return;
        const auto& reference = *cfg_.potential_reference;
        potential_reference_output
            << step_ << ',' << std::setprecision(17) << time_
            << ',' << to_string(reference.axis)
            << ',' << to_string(reference.correction)
            << ',' << reference.coordinate
            << ',' << reference.target
            << ',' << reference.target +
                           potential_reference_offset_
            << ',' << potential_reference_offset_
            << ',' << line_average(
                           mesh_, reference.axis,
                           reference.coordinate)
            << '\n';
        potential_reference_output.flush();
    };
    diag.write_header();
    auto s0 = diag.sample(step_, time_, mesh_, species_, boundary_losses_);
    diag.write_sample(s0);
    write_source_sample();
    write_current_source_sample();
    write_boundary_flux_sample();
    write_potential_reference_sample();
    write_resolved_sample();
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
            write_current_source_sample();
            write_boundary_flux_sample();
            write_potential_reference_sample();
            summary.final_sample = s;
            if (cfg_.vtk_output) {
                write_vtk_outputs(mesh_, cfg_.output_dir, step_, cfg_.vtk_format);
            }
            reached_steady = cfg_.mode == RunMode::SteadyState &&
                             adjacent_energy_windows_converged(diag.history(), cfg_.steady_window, cfg_.steady_tolerance);
            if (reached_steady) summary.steady_state_reached = true;
        }
        write_resolved_sample();
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
    if (resolved_diagnostics) {
        resolved_diagnostics->write_time_averages();
    }
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
