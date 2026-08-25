#include "pic/Simulation.hpp"
#include "pic/Convergence.hpp"
#include "pic/ParticleState.hpp"
#include "pic/Pusher.hpp"
#include "pic/Runtime.hpp"
#include "pic/Units.hpp"
#include <algorithm>
#include <bit>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <numbers>
#include <set>
#include <sstream>
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
constexpr const char* kCheckpointMagicV11 = "AuroraPIC-checkpoint-v11";
constexpr const char* kCheckpointMagicV12 = "AuroraPIC-checkpoint-v12";
constexpr const char* kCheckpointMagicV13 = "AuroraPIC-checkpoint-v13";
constexpr const char* kCheckpointMagicV14 = "AuroraPIC-checkpoint-v14";
constexpr const char* kCheckpointMagicV15 = "AuroraPIC-checkpoint-v15";
constexpr const char* kCheckpointMagicV16 = "AuroraPIC-checkpoint-v16";
constexpr const char* kCheckpointMagicV17 = "AuroraPIC-checkpoint-v17";
constexpr const char* kCheckpointMagicV18 = "AuroraPIC-checkpoint-v18";
constexpr const char* kCheckpointMagicV19 = "AuroraPIC-checkpoint-v19";
constexpr const char* kCheckpointMagicV20 = "AuroraPIC-checkpoint-v20";
constexpr const char* kCheckpointMagicV21 = "AuroraPIC-checkpoint-v21";

void validate_runtime_config(const Config& cfg) {
    if (cfg.velocity_dimensions != 1 &&
        cfg.velocity_dimensions != 3) {
        throw std::invalid_argument(
            "simulation velocity_dimensions must be 1 or 3");
    }
    if (!std::isfinite(cfg.dt) || cfg.dt <= 0.0) throw std::invalid_argument("simulation dt must be positive and finite");
    if (cfg.output_interval == 0) throw std::invalid_argument("output_interval must be positive");
    validate_spatial_average_1d(cfg);
    const auto& impact = cfg.wall_impact_spectrum;
    if (!impact.enabled) {
        if (impact.reset_on_restart || impact.energy_bins != 0 ||
            impact.energy_max != 0.0) {
            throw std::invalid_argument(
                "disabled wall_impact_spectrum cannot configure "
                "restart reset, energy bins, or maximum");
        }
    } else if (impact.energy_bins == 0 ||
               impact.energy_bins > 1000000 ||
               !std::isfinite(impact.energy_max) ||
               !(impact.energy_max > 0.0)) {
        throw std::invalid_argument(
            "wall_impact_spectrum requires 1..1000000 bins and "
            "a positive finite energy maximum");
    }
    if (!std::isfinite(cfg.phi_left) || !std::isfinite(cfg.phi_right)) {
        throw std::invalid_argument("Dirichlet boundary potentials must be finite");
    }
    const auto validate_voltage_drive =
        [](const SinusoidalVoltageConfig& drive,
           const std::string& name) {
            if (!std::isfinite(drive.amplitude) ||
                !std::isfinite(drive.frequency) ||
                drive.frequency < 0.0 ||
                !std::isfinite(drive.phase)) {
                throw std::invalid_argument(
                    name + " values must be finite and frequency "
                    "must be non-negative");
            }
            if (drive.amplitude != 0.0 &&
                !(drive.frequency > 0.0)) {
                throw std::invalid_argument(
                    name + " nonzero amplitude requires positive "
                    "frequency");
            }
        };
    validate_voltage_drive(
        cfg.phi_left_drive, "left sinusoidal electrode drive");
    validate_voltage_drive(
        cfg.phi_right_drive, "right sinusoidal electrode drive");
    const bool driven =
        cfg.phi_left_drive.amplitude != 0.0 ||
        cfg.phi_right_drive.amplitude != 0.0;
    if (cfg.boundary != Boundary::Dirichlet && driven) {
        throw std::invalid_argument(
            "sinusoidal electrode drives require a Dirichlet "
            "boundary");
    }
    if (cfg.mode == RunMode::SteadyState && driven) {
        throw std::invalid_argument(
            "sinusoidal electrode drives require transient mode "
            "until cycle-averaged convergence is implemented");
    }
    if (cfg.checkpoint_output && cfg.checkpoint_interval == 0) {
        throw std::invalid_argument("checkpoint_interval must be positive when checkpoint_output is enabled");
    }
    if (cfg.mode == RunMode::SteadyState) {
        if (cfg.max_steps == 0) throw std::invalid_argument("max_steps must be positive for steady-state mode");
        if (cfg.steady_window == 0) throw std::invalid_argument("steady_window must be positive");
        if (!std::isfinite(cfg.steady_tolerance) || cfg.steady_tolerance <= 0.0) {
            throw std::invalid_argument("steady_tolerance must be positive and finite");
        }
    }
    validate_runtime_policy(cfg.runtime);
    if (!std::isfinite(cfg.collisions.frequency) || cfg.collisions.frequency < 0.0) {
        throw std::invalid_argument("collision frequency must be non-negative and finite");
    }
    if (!std::isfinite(cfg.collisions.neutral_temperature_velocity) || cfg.collisions.neutral_temperature_velocity < 0.0) {
        throw std::invalid_argument("neutral_temperature_velocity must be non-negative and finite");
    }
}

void clear_surface_flux_accumulator(
    PhaseSurfaceFluxAccumulator1D& value) {
    value.macro_crossings = 0;
    value.overflow_macro_crossings = 0;
    value.represented_crossings = 0.0;
    value.overflow_represented_crossings = 0.0;
    value.represented_kinetic_energy = 0.0;
    std::fill(value.represented_histogram.begin(),
              value.represented_histogram.end(), 0.0);
}

void add_surface_flux_accumulator(
    PhaseSurfaceFluxAccumulator1D& target,
    const PhaseSurfaceFluxAccumulator1D& source) {
    target.macro_crossings += source.macro_crossings;
    target.overflow_macro_crossings += source.overflow_macro_crossings;
    target.represented_crossings += source.represented_crossings;
    target.overflow_represented_crossings +=
        source.overflow_represented_crossings;
    target.represented_kinetic_energy += source.represented_kinetic_energy;
    if (target.represented_histogram.size() !=
        source.represented_histogram.size()) {
        throw std::logic_error("surface-flux histogram shape differs");
    }
    for (std::size_t bin = 0; bin < target.represented_histogram.size(); ++bin) {
        target.represented_histogram[bin] +=
            source.represented_histogram[bin];
    }
}
std::filesystem::path checkpoint_path_for_step(const Config& cfg, std::size_t step) {
    if (!cfg.checkpoint_path.empty()) return cfg.checkpoint_path;
    return std::filesystem::path(cfg.output_dir) / ("checkpoint_" + std::to_string(step) + ".apc");
}

void ensure_parent_directory(const std::filesystem::path& path) {
    const auto parent = path.parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
}

template <typename T>
void require_stream(T& stream, const std::string& message) {
    if (!stream) throw std::runtime_error(message);
}

void add_collision_statistics(
    CollisionDiagnostics& destination,
    const CollisionStepStatistics& source,
    std::size_t channel_offset = 0) {
    destination.candidates += source.candidates;
    destination.null_collisions += source.null_collisions;
    if (destination.channel_collisions.size() <
        channel_offset + source.channel_collisions.size()) {
        throw std::logic_error("collision diagnostic channel mismatch");
    }
    if (source.channel_projectile_energy_change.size() !=
            source.channel_collisions.size() ||
        destination.channel_energy_change.size() <
            channel_offset + source.channel_projectile_energy_change.size()) {
        throw std::logic_error("collision energy diagnostic channel mismatch");
    }
    for (std::size_t channel = 0;
         channel < source.channel_collisions.size(); ++channel) {
        destination.channel_collisions[
            channel_offset + channel] +=
            source.channel_collisions[channel];
        destination.channel_energy_change[
            channel_offset + channel] +=
            source.channel_projectile_energy_change[channel];
    }
}

void clear_collision_counts(CollisionDiagnostics& diagnostics) {
    diagnostics.candidates = 0;
    diagnostics.null_collisions = 0;
    std::fill(
        diagnostics.channel_collisions.begin(),
        diagnostics.channel_collisions.end(), 0);
    std::fill(
        diagnostics.channel_energy_change.begin(),
        diagnostics.channel_energy_change.end(), 0.0);
}

void write_collision_header(
    std::ofstream& output,
    const CollisionDiagnostics& diagnostics,
    bool si) {
    const char* energy_suffix = si ? "_J_m-2" : "_normalized";
    output << "step,time,candidates,null_collisions";
    for (const auto& name : diagnostics.channel_names) {
        output << ",collisions_" << name;
    }
    for (const auto& name : diagnostics.channel_names) {
        output << ",tracked_kinetic_energy_change_" << name
               << energy_suffix;
    }
    output << ",cumulative_candidates,cumulative_null_collisions";
    for (const auto& name : diagnostics.channel_names) {
        output << ",cumulative_collisions_" << name;
    }
    for (const auto& name : diagnostics.channel_names) {
        output << ",cumulative_tracked_kinetic_energy_change_" << name
               << energy_suffix;
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
    for (const auto energy : interval.channel_energy_change) {
        output << ',' << energy;
    }
    output << ',' << totals.candidates << ','
           << totals.null_collisions;
    for (const auto count : totals.channel_collisions) {
        output << ',' << count;
    }
    for (const auto energy : totals.channel_energy_change) {
        output << ',' << energy;
    }
    output << '\n';
    output.flush();
}

std::string csv_cell(const std::string& value);

void write_boundary_loss_header(
    std::ofstream& output,
    const std::vector<Species>& species,
    bool si) {
    output << "step,time,counter_origin_step";
    for (const auto& item : species) {
        const auto& name = item.name();
        output << ',' << csv_cell("absorbed_left_count_" + name)
               << ',' << csv_cell("absorbed_right_count_" + name)
               << ',' << csv_cell(
                      "absorbed_left_charge_" + name +
                      (si ? "_C_m-2" : "_normalized"))
               << ',' << csv_cell(
                      "absorbed_right_charge_" + name +
                      (si ? "_C_m-2" : "_normalized"))
               << ',' << csv_cell(
                      "absorbed_left_kinetic_energy_" + name +
                      (si ? "_J_m-2" : "_normalized"))
               << ',' << csv_cell(
                      "absorbed_right_kinetic_energy_" + name +
                      (si ? "_J_m-2" : "_normalized"));
    }
    output << '\n';
}

void write_boundary_loss_sample(
    std::ofstream& output,
    std::size_t step,
    double time,
    std::size_t counter_origin_step,
    const std::vector<Species>& species,
    const std::vector<BoundaryLoss1D>& losses) {
    if (species.size() != losses.size()) {
        throw std::logic_error(
            "1D boundary-loss diagnostics do not match species");
    }
    output << step << ',' << std::setprecision(17) << time
           << ',' << counter_origin_step;
    for (std::size_t species_id = 0;
         species_id < species.size(); ++species_id) {
        const auto& item = species[species_id];
        const auto& loss = losses[species_id];
        const double charge_per_macro =
            item.charge() * item.config().weight;
        output << ',' << loss.absorbed_left
               << ',' << loss.absorbed_right
               << ',' << charge_per_macro *
                              static_cast<double>(
                                  loss.absorbed_left)
               << ',' << charge_per_macro *
                              static_cast<double>(
                                  loss.absorbed_right)
               << ',' << loss.kinetic_energy_left
               << ',' << loss.kinetic_energy_right;
    }
    output << '\n';
    output.flush();
}

void clear_wall_impact_side(
    WallImpactSideSpectrum1D& side) {
    side.macro_impacts = 0;
    side.overflow_macro_impacts = 0;
    side.represented_impacts = 0.0;
    side.overflow_represented_impacts = 0.0;
    side.represented_kinetic_energy = 0.0;
    std::fill(
        side.macro_histogram.begin(),
        side.macro_histogram.end(), 0);
    std::fill(
        side.represented_histogram.begin(),
        side.represented_histogram.end(), 0.0);
}

void add_wall_impact_side(
    WallImpactSideSpectrum1D& destination,
    const WallImpactSideSpectrum1D& source) {
    if (destination.macro_histogram.size() !=
            source.macro_histogram.size() ||
        destination.represented_histogram.size() !=
            source.represented_histogram.size()) {
        throw std::logic_error(
            "1D wall-impact histogram shape mismatch");
    }
    destination.macro_impacts += source.macro_impacts;
    destination.overflow_macro_impacts +=
        source.overflow_macro_impacts;
    destination.represented_impacts += source.represented_impacts;
    destination.overflow_represented_impacts +=
        source.overflow_represented_impacts;
    destination.represented_kinetic_energy +=
        source.represented_kinetic_energy;
    for (std::size_t bin = 0;
         bin < destination.macro_histogram.size(); ++bin) {
        destination.macro_histogram[bin] +=
            source.macro_histogram[bin];
        destination.represented_histogram[bin] +=
            source.represented_histogram[bin];
    }
}

void write_power_transfer_header(
    std::ofstream& output,
    const std::vector<Species>& species,
    bool si) {
    output << "step,time,counter_origin_step";
    for (const auto& item : species) {
        output << ',' << csv_cell(
            "electric_work_" + item.name() +
            (si ? "_J_m-2" : "_normalized"));
    }
    output << '\n';
}

void write_power_transfer_sample(
    std::ofstream& output,
    std::size_t step,
    double time,
    std::size_t counter_origin_step,
    const std::vector<Species>& species,
    const std::vector<SpeciesPower1D>& power) {
    if (species.size() != power.size()) {
        throw std::logic_error(
            "1D power-transfer diagnostics do not match species");
    }
    output << step << ',' << std::setprecision(17) << time
           << ',' << counter_origin_step;
    for (const auto& item : power) {
        output << ',' << item.electric_work;
    }
    output << '\n';
    output.flush();
}

std::string csv_cell(const std::string& value) {
    if (value.find_first_of(",\"\r\n") == std::string::npos) {
        return value;
    }
    std::string escaped{"\""};
    for (const char character : value) {
        if (character == '"') escaped.push_back('"');
        escaped.push_back(character);
    }
    escaped.push_back('"');
    return escaped;
}

std::string json_string(const std::string& value) {
    std::string escaped{"\""};
    for (const unsigned char character : value) {
        switch (character) {
            case '"': escaped += "\\\""; break;
            case '\\': escaped += "\\\\"; break;
            case '\b': escaped += "\\b"; break;
            case '\f': escaped += "\\f"; break;
            case '\n': escaped += "\\n"; break;
            case '\r': escaped += "\\r"; break;
            case '\t': escaped += "\\t"; break;
            default:
                if (character < 0x20) {
                    constexpr char hex[] = "0123456789abcdef";
                    escaped += "\\u00";
                    escaped.push_back(hex[character >> 4]);
                    escaped.push_back(hex[character & 0x0f]);
                } else {
                    escaped.push_back(
                        static_cast<char>(character));
                }
        }
    }
    escaped.push_back('"');
    return escaped;
}

} // namespace

Simulation::Simulation(Config cfg)
    : cfg_(std::move(cfg)),
      grid_(cfg_.nx, cfg_.length, cfg_.boundary),
      solver_(cfg_.units.permittivity()),
      rng_(cfg_.seed) {
    if (cfg_.checkpoint_output && cfg_.checkpoint_interval == 0) cfg_.checkpoint_interval = cfg_.output_interval;
    validate_runtime_config(cfg_);
    if (!cfg_.restart_path.empty() &&
        !cfg_.initial_state_path.empty()) {
        throw std::invalid_argument(
            "restart_path and initial_state_path are mutually exclusive");
    }
    if (cfg_.initial_state_signature &&
        cfg_.initial_state_path.empty()) {
        throw std::invalid_argument(
            "initial_state_signature requires initial_state_path");
    }
    validate_initialization_acceptance(
        cfg_.initialization_acceptance,
        "1D initialization acceptance config");
    if (cfg_.max_particles_per_species == 0) {
        throw std::invalid_argument(
            "max_particles_per_species must be positive");
    }
    for (const auto& sc : cfg_.species) {
        if (sc.timestep_multiplier == 0 ||
            !std::isfinite(
                cfg_.dt * static_cast<double>(
                    sc.timestep_multiplier))) {
            throw std::invalid_argument(
                "species '" + sc.name +
                "' timestep_multiplier must produce a positive "
                "finite particle timestep");
        }
        if (sc.particles > cfg_.max_particles_per_species) {
            throw std::invalid_argument(
                "initial particle count exceeds "
                "max_particles_per_species for species '" +
                sc.name + "'");
        }
        if (std::any_of(
                species_.begin(), species_.end(),
                [&](const Species& species) {
                    return species.name() == sc.name;
                })) {
            throw std::invalid_argument(
                "1D species names must be unique: " + sc.name);
        }
        species_.emplace_back(sc, cfg_.velocity_dimensions);
    }
    if (!cfg_.collision_models.empty() &&
        cfg_.collisions.enabled) {
        throw std::invalid_argument(
            "named collision models cannot be combined with the "
            "legacy collision model");
    }
    const auto species_id =
        [&](const std::string& name) {
        const auto target = std::find_if(
            species_.begin(), species_.end(),
            [&](const Species& species) {
                return species.name() == name;
            });
        if (target == species_.end()) {
            throw std::invalid_argument(
                "MCC species does not exist: " + name);
        }
        return static_cast<std::size_t>(
            target - species_.begin());
    };
    const auto add_mcc =
        [&](const std::string& name,
            const CollisionConfig& collision,
            bool qualify_diagnostics) {
        if (collision.model !=
            CollisionModelKind::NullCollision) {
            throw std::invalid_argument(
                "named 1D collision models support only "
                "null_collision");
        }
        MccRuntime runtime;
        runtime.name = name;
        runtime.species_id =
            species_id(collision.species);
        runtime.diagnostic_offset =
            collision_totals_.channel_names.size();
        runtime.model =
            std::make_unique<NullCollisionModel>(
                collision,
                species_[runtime.species_id].mass());
        runtime.ionization_channels.resize(
            collision.channels.size());
        runtime.channel_processes.reserve(
            collision.channels.size());
        const auto& target =
            species_[runtime.species_id].config();
        for (std::size_t channel = 0;
             channel < collision.channels.size(); ++channel) {
            const auto& channel_config =
                collision.channels[channel];
            runtime.channel_processes.push_back(
                channel_config.process);
            if (channel_config.process ==
                CollisionProcessKind::Attachment) {
                throw std::invalid_argument(
                    "1D MCC does not support attachment");
            }
            if (channel_config.process !=
                CollisionProcessKind::Ionization) {
                continue;
            }
            if (cfg_.velocity_dimensions != 3) {
                throw std::invalid_argument(
                    "1D ionization requires "
                    "velocity_dimensions = 3");
            }
            const std::size_t secondary =
                species_id(
                    channel_config.secondary_species);
            const std::size_t ion =
                species_id(channel_config.ion_species);
            const auto& secondary_config =
                species_[secondary].config();
            const auto& ion_config =
                species_[ion].config();
            if (target.charge == 0.0 ||
                target.weight != secondary_config.weight ||
                target.weight != ion_config.weight ||
                target.mass != secondary_config.mass ||
                target.charge != secondary_config.charge ||
                ion_config.charge != -target.charge) {
                throw std::invalid_argument(
                    "1D ionization products violate the "
                    "equal-weight charge-conservation contract");
            }
            runtime.ionization_channels[channel] =
                IonizationChannelRuntime{secondary, ion};
        }
        for (const auto& channel :
             runtime.model->channel_names()) {
            collision_totals_.channel_names.push_back(
                qualify_diagnostics
                    ? name + "." + channel
                    : channel);
        }
        mcc_models_.push_back(std::move(runtime));
    };
    if (cfg_.collisions.enabled) {
        if (cfg_.collisions.model ==
            CollisionModelKind::BGK) {
            legacy_bgk_enabled_ = true;
            collision_totals_.channel_names = {"bgk"};
        } else {
            add_mcc("legacy", cfg_.collisions, false);
        }
    }
    std::set<std::string> model_names;
    std::set<std::size_t> target_species;
    for (const auto& named : cfg_.collision_models) {
        if (named.name.empty() ||
            !model_names.insert(named.name).second) {
            throw std::invalid_argument(
                "named collision model names must be unique");
        }
        if (!std::all_of(
                named.name.begin(), named.name.end(),
                [](unsigned char character) {
                    return std::isalnum(character) ||
                           character == '_' ||
                           character == '-';
                })) {
            throw std::invalid_argument(
                "named collision model name has an invalid "
                "character");
        }
        if (!named.config.enabled) continue;
        const std::size_t target =
            species_id(named.config.species);
        if (!target_species.insert(target).second) {
            throw std::invalid_argument(
                "multiple collision models target species '" +
                named.config.species + "'");
        }
        add_mcc(named.name, named.config, true);
    }
    collision_totals_.channel_collisions.assign(
        collision_totals_.channel_names.size(), 0);
    collision_totals_.channel_energy_change.assign(
        collision_totals_.channel_names.size(), 0.0);
    collision_interval_.channel_names =
        collision_totals_.channel_names;
    collision_interval_.channel_collisions.assign(
        collision_totals_.channel_names.size(), 0);
    collision_interval_.channel_energy_change.assign(
        collision_totals_.channel_names.size(), 0.0);
    species_boundary_losses_.assign(
        species_.size(), BoundaryLoss1D{});
    if (cfg_.wall_impact_spectrum.enabled) {
        wall_impact_spectra_.resize(species_.size());
        wall_impact_chunks_.assign(
            species_.size(),
            std::vector<SpeciesWallImpactSpectrum1D>(
                runtime_info(cfg_.runtime).active_threads));
        const auto initialize_spectrum = [&](
            SpeciesWallImpactSpectrum1D& spectrum) {
            spectrum.left.macro_histogram.assign(
                cfg_.wall_impact_spectrum.energy_bins, 0);
            spectrum.left.represented_histogram.assign(
                cfg_.wall_impact_spectrum.energy_bins, 0.0);
            spectrum.right.macro_histogram.assign(
                cfg_.wall_impact_spectrum.energy_bins, 0);
            spectrum.right.represented_histogram.assign(
                cfg_.wall_impact_spectrum.energy_bins, 0.0);
        };
        for (auto& spectrum : wall_impact_spectra_) {
            initialize_spectrum(spectrum);
        }
        for (auto& species_chunks : wall_impact_chunks_) {
            for (auto& spectrum : species_chunks) {
                initialize_spectrum(spectrum);
            }
        }
    }
    boundary_loss_chunks_.assign(
        species_.size(),
        std::vector<BoundaryLoss1D>(
            runtime_info(cfg_.runtime).active_threads,
            BoundaryLoss1D{}));
    species_power_transfer_.assign(
        species_.size(), SpeciesPower1D{});
    power_transfer_chunks_.assign(
        species_.size(),
        std::vector<SpeciesPower1D>(
            runtime_info(cfg_.runtime).active_threads,
            SpeciesPower1D{}));
    if (cfg_.spatial_average.enabled) {
        spatial_density_sums_.assign(
            species_.size(),
            std::vector<double>(grid_.nx(), 0.0));
        spatial_density_scratch_.assign(
            grid_.nx(), 0.0);
        spatial_kinetic_energy_sums_.assign(
            species_.size(), std::vector<double>(grid_.nx(), 0.0));
        spatial_kinetic_energy_scratch_.assign(grid_.nx(), 0.0);
        spatial_velocity_x_scratch_.assign(grid_.nx(), 0.0);
        spatial_velocity_y_scratch_.assign(grid_.nx(), 0.0);
        spatial_velocity_z_scratch_.assign(grid_.nx(), 0.0);
        spatial_potential_sums_.assign(grid_.nx(), 0.0);
        spatial_electric_sums_.assign(grid_.nx(), 0.0);
        spatial_electric_squared_sums_.assign(grid_.nx(), 0.0);
        spatial_phase_bins_.resize(cfg_.spatial_average.phase_bins);
        for (auto& bin : spatial_phase_bins_) {
            const auto species_nodes = std::vector<std::vector<double>>(
                species_.size(), std::vector<double>(grid_.nx(), 0.0));
            bin.density = species_nodes;
            bin.velocity_x_density = species_nodes;
            bin.velocity_y_density = species_nodes;
            bin.velocity_z_density = species_nodes;
            bin.kinetic_energy_density = species_nodes;
            bin.potential.assign(grid_.nx(), 0.0);
            bin.electric.assign(grid_.nx(), 0.0);
            bin.electric_squared.assign(grid_.nx(), 0.0);
        }
        spatial_collision_energy_sums_.assign(
            collision_totals_.channel_names.size(),
            std::vector<double>(grid_.nx(), 0.0));
        spatial_collision_phase_steps_.assign(
            cfg_.spatial_average.phase_bins, 0);
        spatial_collision_phase_energy_sums_.assign(
            cfg_.spatial_average.phase_bins,
            std::vector<std::vector<double>>(
                collision_totals_.channel_names.size(),
                std::vector<double>(grid_.nx(), 0.0)));
        spatial_collision_event_sums_.assign(
            collision_totals_.channel_names.size(),
            std::vector<double>(grid_.nx(), 0.0));
        spatial_collision_phase_event_sums_.assign(
            cfg_.spatial_average.phase_bins,
            std::vector<std::vector<double>>(
                collision_totals_.channel_names.size(),
                std::vector<double>(grid_.nx(), 0.0)));
    }
    if (cfg_.phase_eedf.enabled) {
        phase_eedf_species_id_ = species_id(cfg_.phase_eedf.species);
        phase_eedf_accumulators_.assign(
            cfg_.spatial_average.phase_bins,
            std::vector<PhaseEedfAccumulator1D>(
                cfg_.phase_eedf.regions.size()));
        for (auto& phase : phase_eedf_accumulators_) {
            for (auto& region : phase) {
                region.histogram.assign(
                    cfg_.phase_eedf.energy_bins, 0.0);
            }
        }
        if (cfg_.phase_eedf.history_enabled) {
            phase_eedf_threshold_crossings_.assign(
                cfg_.spatial_average.phase_bins,
                std::vector<PhaseEedfThresholdCrossingAccumulator1D>(
                    cfg_.phase_eedf.regions.size()));
        }
    }
    if (cfg_.phase_surface_flux.enabled) {
        phase_surface_flux_species_id_ =
            species_id(cfg_.phase_surface_flux.species);
        const auto make_surfaces = [&] {
            auto surfaces = std::vector<std::vector<
                PhaseSurfaceFluxAccumulator1D>>(
                    cfg_.phase_surface_flux.positions.size(),
                    std::vector<PhaseSurfaceFluxAccumulator1D>(2));
            for (auto& surface : surfaces) {
                for (auto& direction : surface) {
                    direction.represented_histogram.assign(
                        cfg_.phase_surface_flux.energy_bins, 0.0);
                }
            }
            return surfaces;
        };
        phase_surface_flux_accumulators_.resize(
            cfg_.spatial_average.phase_bins);
        for (auto& phase : phase_surface_flux_accumulators_) {
            phase = make_surfaces();
        }
        phase_surface_flux_chunks_.resize(
            runtime_info(cfg_.runtime).active_threads);
        for (auto& chunk : phase_surface_flux_chunks_) {
            chunk = make_surfaces();
        }
    }
}

std::uint64_t Simulation::collision_signature() const {
    if (mcc_models_.empty()) {
        if (!legacy_bgk_enabled_) return 0;
        constexpr std::uint64_t offset =
            1469598103934665603ULL;
        constexpr std::uint64_t prime =
            1099511628211ULL;
        std::uint64_t hash = offset;
        const auto append = [&](std::uint64_t value) {
            for (unsigned byte = 0; byte < 8; ++byte) {
                hash ^= static_cast<unsigned char>(
                    value >> (byte * 8));
                hash *= prime;
            }
        };
        append(std::bit_cast<std::uint64_t>(
            cfg_.collisions.frequency));
        append(std::bit_cast<std::uint64_t>(
            cfg_.collisions.neutral_temperature_velocity));
        return hash;
    }
    if (mcc_models_.size() == 1 &&
        cfg_.collision_models.empty()) {
        return mcc_models_.front().model->signature();
    }
    constexpr std::uint64_t offset =
        1469598103934665603ULL;
    constexpr std::uint64_t prime =
        1099511628211ULL;
    std::uint64_t hash = offset;
    const auto append_byte = [&](unsigned char value) {
        hash ^= value;
        hash *= prime;
    };
    const auto append_string = [&](const std::string& value) {
        for (const unsigned char character : value) {
            append_byte(character);
        }
        append_byte(0);
    };
    const auto append_integer = [&](std::uint64_t value) {
        for (unsigned byte = 0; byte < 8; ++byte) {
            append_byte(static_cast<unsigned char>(
                value >> (byte * 8)));
        }
    };
    for (const auto& runtime : mcc_models_) {
        append_string(runtime.name);
        const auto append_species =
            [&](std::size_t species_id) {
            const auto& species =
                species_[species_id].config();
            append_string(species.name);
            append_integer(std::bit_cast<std::uint64_t>(
                species.charge));
            append_integer(std::bit_cast<std::uint64_t>(
                species.mass));
            append_integer(std::bit_cast<std::uint64_t>(
                species.weight));
        };
        append_species(runtime.species_id);
        append_integer(runtime.model->signature());
        for (const auto& channel :
             runtime.ionization_channels) {
            append_integer(channel.has_value() ? 1 : 0);
            if (!channel) continue;
            append_species(
                channel->secondary_species_id);
            append_species(channel->ion_species_id);
        }
    }
    return hash;
}

std::string Simulation::collision_identity() const {
    std::string identity;
    if (!cfg_.collision_models.empty()) {
        identity = "named_null_collision";
    } else {
        identity = to_string(cfg_.collisions.model);
    }
    if (cfg_.collision_velocity_sampling ==
        CollisionVelocitySampling1D::LeapfrogHalfStep) {
        identity += "@leapfrog_half_step";
    }
    return identity;
}

double Simulation::electrode_potential(
    double offset,
    const SinusoidalVoltageConfig& drive,
    double field_time) const {
    return offset + drive.amplitude * std::sin(
        2.0 * std::numbers::pi * drive.frequency * field_time +
        drive.phase);
}

bool Simulation::species_due(std::size_t species_id) const {
    return step_ % species_[species_id].config().timestep_multiplier == 0;
}

double Simulation::species_timestep(std::size_t species_id) const {
    return cfg_.dt * static_cast<double>(
        species_[species_id].config().timestep_multiplier);
}

void Simulation::deposit_and_solve(double field_time) {
    grid_.clear_charge();
    for (const auto& sp : species_) sp.deposit_charge(grid_);
    solver_.solve(
        grid_,
        electrode_potential(
            cfg_.phi_left, cfg_.phi_left_drive, field_time),
        electrode_potential(
            cfg_.phi_right, cfg_.phi_right_drive, field_time));
}

void Simulation::initialize() {
    time_ = 0.0;
    step_ = 0;
    boundary_loss_origin_step_ = 0;
    std::fill(
        species_boundary_losses_.begin(),
        species_boundary_losses_.end(),
        BoundaryLoss1D{});
    wall_impact_origin_step_ = 0;
    for (auto& spectrum : wall_impact_spectra_) {
        spectrum.baseline_loss = {};
        clear_wall_impact_side(spectrum.left);
        clear_wall_impact_side(spectrum.right);
    }
    power_transfer_origin_step_ = 0;
    std::fill(
        species_power_transfer_.begin(),
        species_power_transfer_.end(),
        SpeciesPower1D{});
    if (cfg_.initial_state_path.empty()) {
        for (auto& sp : species_) sp.initialize(grid_, rng_);
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
                cfg_.initial_state_path, 1,
                cfg_.velocity_dimensions,
                cfg_.units.system, expected,
                "1D simulation",
                [&](std::size_t species_index,
                    std::size_t record_index,
                    const ExternalParticleRecord& record) {
                    auto& species =
                        species_.at(species_index);
                    const double minimum =
                        species.config().init_x_min;
                    const double maximum =
                        species.config().init_x_max < 0.0
                            ? grid_.length()
                            : species.config().init_x_max;
                    if (record.position.x < minimum ||
                        record.position.x > maximum ||
                        (grid_.boundary() == Boundary::Periodic &&
                         record.position.x == grid_.length())) {
                        throw std::runtime_error(
                            "external particle for species '" +
                            species.name() +
                            "' lies outside its 1D initialization interval");
                    }
                    auto& particle =
                        species.particles().at(record_index);
                    particle.x = record.position.x;
                    particle.v = record.velocity.x;
                    particle.velocity_y = record.velocity.y;
                    particle.velocity_z = record.velocity.z;
                    particle.v_half = record.velocity.x;
                    particle.alive = true;
                },
                cfg_.initial_state_signature);
    }
    if (cfg_.phase_eedf.history_enabled) {
        phase_eedf_particle_histories_.assign(
            species_[phase_eedf_species_id_].particles().size(), {});
    }
    deposit_and_solve(time_);
    for (std::size_t species_id = 0;
         species_id < species_.size(); ++species_id) {
        auto& sp = species_[species_id];
        const double qm = sp.charge() / sp.mass();
        const double timestep = species_timestep(species_id);
        auto& particles = sp.particles();
        runtime_parallel_for(std::size_t{0}, particles.size(), cfg_.runtime, [&](std::size_t particle_id) {
            auto& p = particles[particle_id];
            if (p.alive) initialize_leapfrog_half_step(p, interpolate_electric(grid_, p.x), qm, timestep);
        });
    }
    initialized_ = true;
}

void Simulation::apply_collisions() {
    const auto collision_velocity = [&](const Particle& particle) {
        return Vec3{
            cfg_.collision_velocity_sampling ==
                    CollisionVelocitySampling1D::LeapfrogHalfStep
                ? particle.v_half
                : particle.v,
            particle.velocity_y,
            particle.velocity_z};
    };
    const auto store_collision_velocity = [this](
        Particle& particle, const Vec3& velocity,
        double charge_to_mass, double timestep) {
        particle.velocity_y = velocity.y;
        particle.velocity_z = velocity.z;
        if (cfg_.collision_velocity_sampling ==
            CollisionVelocitySampling1D::LeapfrogHalfStep) {
            particle.v_half = velocity.x;
            synchronize_leapfrog(
                particle,
                interpolate_electric(grid_, particle.x),
                charge_to_mass, timestep);
        } else {
            particle.v = velocity.x;
            initialize_leapfrog_half_step(
                particle,
                interpolate_electric(grid_, particle.x),
                charge_to_mass, timestep);
        }
    };
    if (legacy_bgk_enabled_ &&
        cfg_.collisions.frequency > 0.0) {
        for (std::size_t species_id = 0;
             species_id < species_.size(); ++species_id) {
            if (!species_due(species_id)) continue;
            auto& species = species_[species_id];
            const double timestep = species_timestep(species_id);
            const double probability =
                1.0 - std::exp(
                    -cfg_.collisions.frequency * timestep);
            std::uniform_real_distribution<double> unit(
                0.0, 1.0);
            std::normal_distribution<double> neutral_velocity(
                0.0,
                cfg_.collisions.neutral_temperature_velocity);
            const double charge_to_mass =
                species.charge() / species.mass();
            auto& particles = species.particles();
            for (std::size_t particle_id = 0;
                 particle_id < particles.size(); ++particle_id) {
                auto& particle = particles[particle_id];
                if (!particle.alive ||
                    unit(rng_) >= probability) {
                    continue;
                }
                Vec3 velocity = collision_velocity(particle);
                const double energy_before = 0.5 * species.mass() *
                    species.weight() *
                    (velocity.x * velocity.x +
                     (cfg_.velocity_dimensions == 3
                          ? velocity.y * velocity.y +
                            velocity.z * velocity.z
                          : 0.0));
                velocity.x = neutral_velocity(rng_);
                if (cfg_.velocity_dimensions == 3) {
                    velocity.y = neutral_velocity(rng_);
                    velocity.z = neutral_velocity(rng_);
                }
                ++collision_totals_.candidates;
                ++collision_totals_.channel_collisions[0];
                ++collision_interval_.candidates;
                ++collision_interval_.channel_collisions[0];
                if (phase_eedf_history_active() &&
                    species_id == phase_eedf_species_id_) {
                    auto& history =
                        phase_eedf_particle_histories_[particle_id];
                    ++history.bgk_collisions;
                }
                deposit_spatial_collision_events(
                    particle.x, 0, species.weight());
                const double energy_after = 0.5 * species.mass() *
                    species.weight() *
                    (velocity.x * velocity.x +
                     (cfg_.velocity_dimensions == 3
                          ? velocity.y * velocity.y +
                            velocity.z * velocity.z
                          : 0.0));
                if (species_id == phase_eedf_species_id_ &&
                    phase_eedf_history_active()) {
                    const double energy_scale =
                        cfg_.units.system == UnitSystem::SI
                            ? ELEMENTARY_CHARGE_SI : 1.0;
                    const double represented_scale =
                        species.weight() * energy_scale;
                    add_phase_eedf_bgk_transition(
                        particle.x,
                        energy_before / represented_scale >=
                            cfg_.phase_eedf.tail_threshold,
                        energy_after / represented_scale >=
                            cfg_.phase_eedf.tail_threshold);
                }
                const double energy_change = energy_after - energy_before;
                collision_totals_.channel_energy_change[0] += energy_change;
                collision_interval_.channel_energy_change[0] += energy_change;
                deposit_spatial_collision_energy(
                    particle.x, 0, energy_change);
                store_collision_velocity(
                    particle, velocity, charge_to_mass, timestep);
            }
        }
    }
    struct IonizationProduct {
        double position{0.0};
        Vec3 secondary_velocity{};
        Vec3 ion_velocity{};
        IonizationChannelRuntime channel{};
    };
    std::vector<IonizationProduct> products;
    for (auto& runtime : mcc_models_) {
        if (!species_due(runtime.species_id)) continue;
        auto& sp = species_[runtime.species_id];
        const double qm = sp.charge() / sp.mass();
        const double timestep = species_timestep(runtime.species_id);
        const std::size_t initial_particle_count =
            sp.particles().size();
        for (std::size_t particle_id = 0;
             particle_id < initial_particle_count;
             ++particle_id) {
            auto& part = sp.particles()[particle_id];
            if (!part.alive) continue;
            auto& statistics =
                runtime.collision_workspace.statistics;
            Vec3 collision_velocity_before = collision_velocity(part);
            Vec3 collision_velocity_after = collision_velocity_before;
            if (cfg_.velocity_dimensions == 3) {
                runtime.model->collide_reusing_storage(
                    collision_velocity_after, timestep, rng_,
                    runtime.collision_workspace);
                store_collision_velocity(
                    part, collision_velocity_after, qm, timestep);
            } else {
                double velocity = collision_velocity_before.x;
                runtime.model->collide_reusing_storage(
                    velocity, timestep, rng_,
                    runtime.collision_workspace);
                collision_velocity_after = Vec3{velocity, 0.0, 0.0};
                store_collision_velocity(
                    part, collision_velocity_after, qm, timestep);
            }
            for (auto& energy_change :
                 statistics.channel_projectile_energy_change) {
                energy_change *= sp.weight();
            }
            runtime.tracked_energy_scratch =
                statistics.channel_projectile_energy_change;
            auto& tracked_energy_change =
                runtime.tracked_energy_scratch;
            add_collision_statistics(
                collision_totals_, statistics,
                runtime.diagnostic_offset);
            add_collision_statistics(
                collision_interval_, statistics,
                runtime.diagnostic_offset);
            bool transition_process_present = false;
            CollisionProcessKind transition_process =
                CollisionProcessKind::Elastic;
            for (std::size_t channel = 0;
                 channel < statistics.channel_collisions.size(); ++channel) {
                if (statistics.channel_collisions[channel] != 0) {
                    if (runtime.species_id == phase_eedf_species_id_ &&
                        channel < runtime.channel_processes.size()) {
                        add_phase_eedf_collision_history(
                            particle_id,
                            runtime.channel_processes[channel],
                            statistics.channel_collisions[channel]);
                        if (!transition_process_present) {
                            transition_process =
                                runtime.channel_processes[channel];
                            transition_process_present = true;
                        }
                    }
                    deposit_spatial_collision_events(
                        part.x, runtime.diagnostic_offset + channel,
                        sp.weight() * static_cast<double>(
                            statistics.channel_collisions[channel]));
                }
            }
            if (transition_process_present &&
                runtime.species_id == phase_eedf_species_id_ &&
                phase_eedf_history_active()) {
                const auto energy_eV = [&](const Vec3& velocity) {
                    const double squared = velocity.x * velocity.x +
                        (cfg_.velocity_dimensions == 3
                            ? velocity.y * velocity.y +
                              velocity.z * velocity.z
                            : 0.0);
                    const double energy_scale =
                        cfg_.units.system == UnitSystem::SI
                            ? ELEMENTARY_CHARGE_SI : 1.0;
                    return 0.5 * sp.mass() * squared / energy_scale;
                };
                add_phase_eedf_collision_transition(
                    part.x, transition_process,
                    energy_eV(collision_velocity_before) >=
                        cfg_.phase_eedf.tail_threshold,
                    energy_eV(collision_velocity_after) >=
                        cfg_.phase_eedf.tail_threshold);
            }
            for (const auto& secondary :
                 statistics.secondaries) {
                if (secondary.channel >=
                        runtime.ionization_channels.size() ||
                    !runtime.ionization_channels[
                        secondary.channel]) {
                    throw std::logic_error(
                        "MCC produced an unmapped 1D "
                        "ionization channel");
                }
                products.push_back({
                    part.x, secondary.velocity,
                    secondary.ion_velocity,
                    *runtime.ionization_channels[
                        secondary.channel]});
                const auto& channel = *runtime.ionization_channels[
                    secondary.channel];
                const auto product_energy = [](const Vec3& velocity) {
                    return velocity.x * velocity.x +
                           velocity.y * velocity.y +
                           velocity.z * velocity.z;
                };
                const double added_energy =
                    0.5 * species_[channel.secondary_species_id].mass() *
                        species_[channel.secondary_species_id].weight() *
                        product_energy(secondary.velocity) +
                    0.5 * species_[channel.ion_species_id].mass() *
                        species_[channel.ion_species_id].weight() *
                        product_energy(secondary.ion_velocity);
                collision_totals_.channel_energy_change[
                    runtime.diagnostic_offset + secondary.channel] +=
                    added_energy;
                collision_interval_.channel_energy_change[
                    runtime.diagnostic_offset + secondary.channel] +=
                    added_energy;
                tracked_energy_change[secondary.channel] +=
                    added_energy;
            }
            if (statistics.primary_removal_channel) {
                throw std::logic_error(
                    "1D MCC produced an unsupported primary "
                    "removal event");
            }
            for (std::size_t channel = 0;
                 channel < tracked_energy_change.size(); ++channel) {
                if (tracked_energy_change[channel] != 0.0) {
                    deposit_spatial_collision_energy(
                        part.x,
                        runtime.diagnostic_offset + channel,
                        tracked_energy_change[channel]);
                }
            }
        }
    }
    std::vector<std::size_t> required_products(
        species_.size(), 0);
    for (const auto& product : products) {
        ++required_products[
            product.channel.secondary_species_id];
        ++required_products[
            product.channel.ion_species_id];
    }
    for (std::size_t species_id = 0;
         species_id < species_.size(); ++species_id) {
        const auto& particles =
            species_[species_id].particles();
        if (particles.size() >
            cfg_.max_particles_per_species) {
            throw std::runtime_error(
                "species storage already exceeds "
                "max_particles_per_species");
        }
        const std::size_t reusable =
            static_cast<std::size_t>(
                std::count_if(
                    particles.begin(), particles.end(),
                    [](const Particle& particle) {
                        return !particle.alive;
                    }));
        const std::size_t available =
            cfg_.max_particles_per_species -
                particles.size() +
            reusable;
        if (required_products[species_id] > available) {
            throw std::runtime_error(
                "reactive collisions exceeded "
                "max_particles_per_species for species '" +
                species_[species_id].name() + "'");
        }
    }
    std::vector<std::size_t> next_reusable_slot(
        species_.size(), 0);
    const auto append_product =
        [&](std::size_t species_id,
            double position,
            Vec3 velocity) {
        auto& product_species =
            species_[species_id];
        auto& particles =
            product_species.particles();
        auto& slot = next_reusable_slot[species_id];
        while (slot < particles.size() &&
               particles[slot].alive) {
            ++slot;
        }
        if (slot == particles.size()) {
            particles.emplace_back();
            if (cfg_.phase_eedf.history_enabled &&
                species_id == phase_eedf_species_id_) {
                phase_eedf_particle_histories_.emplace_back();
            }
        }
        auto& product = particles[slot++];
        product = {};
        product.x = position;
        product.v = velocity.x;
        product.velocity_y = velocity.y;
        product.velocity_z = velocity.z;
        product.alive = true;
        if (cfg_.phase_eedf.history_enabled &&
            species_id == phase_eedf_species_id_) {
            if (slot - 1 >= phase_eedf_particle_histories_.size()) {
                throw std::logic_error(
                    "phase EEDF product history is misaligned");
            }
            auto& history =
                phase_eedf_particle_histories_[slot - 1];
            history = {};
            history.born_during_window =
                phase_eedf_history_active();
            if (phase_eedf_history_active()) {
                const double squared = velocity.x * velocity.x +
                    (cfg_.velocity_dimensions == 3
                        ? velocity.y * velocity.y + velocity.z * velocity.z
                        : 0.0);
                const double energy_scale =
                    cfg_.units.system == UnitSystem::SI
                        ? ELEMENTARY_CHARGE_SI : 1.0;
                add_phase_eedf_birth(
                    position,
                    0.5 * product_species.mass() * squared / energy_scale >=
                        cfg_.phase_eedf.tail_threshold);
            }
        }
        store_collision_velocity(
            product, velocity,
            product_species.charge() /
                product_species.mass(),
            species_timestep(species_id));
    };
    for (const auto& product : products) {
        append_product(
            product.channel.secondary_species_id,
            product.position,
            product.secondary_velocity);
        append_product(
            product.channel.ion_species_id,
            product.position,
            product.ion_velocity);
    }
}

void Simulation::step() {
    if (!initialized_) initialize();
    const std::size_t surface_flux_phase =
        phase_surface_flux_phase(step_ + 1);
    for (std::size_t species_id = 0;
         species_id < species_.size(); ++species_id) {
        auto& sp = species_[species_id];
        const double qm = sp.charge() / sp.mass();
        auto& particles = sp.particles();
        auto& chunk_losses =
            boundary_loss_chunks_[species_id];
        std::fill(
            chunk_losses.begin(), chunk_losses.end(),
            BoundaryLoss1D{});
        if (cfg_.wall_impact_spectrum.enabled) {
            for (auto& spectrum :
                 wall_impact_chunks_[species_id]) {
                clear_wall_impact_side(spectrum.left);
                clear_wall_impact_side(spectrum.right);
            }
        }
        auto& chunk_power =
            power_transfer_chunks_[species_id];
        std::fill(
            chunk_power.begin(), chunk_power.end(),
            SpeciesPower1D{});
        if (!species_due(species_id)) continue;
        const double timestep = species_timestep(species_id);
        if (species_id == phase_eedf_species_id_ &&
            phase_eedf_history_active()) {
            phase_eedf_field_push_origin_energetic_.assign(
                particles.size(), 0);
        }
        runtime_static_chunks(
            std::size_t{0}, particles.size(), cfg_.runtime,
            [&](std::size_t chunk, std::size_t begin,
                std::size_t end) {
                auto& loss = chunk_losses[chunk];
                auto& power = chunk_power[chunk];
                auto* impact = cfg_.wall_impact_spectrum.enabled
                    ? &wall_impact_chunks_[species_id][chunk]
                    : nullptr;
                for (std::size_t particle_id = begin;
                     particle_id < end; ++particle_id) {
                    auto& p = particles[particle_id];
                    if (!p.alive) continue;
                    if (species_id == phase_eedf_species_id_ &&
                        phase_eedf_history_active()) {
                        phase_eedf_field_push_origin_energetic_[particle_id] =
                            phase_eedf_collision_state_energetic(p, sp) ? 1 : 0;
                    }
                    const double old_position = p.x;
                    kick_leapfrog(
                        p, interpolate_electric(grid_, p.x),
                        qm, timestep);
                    drift_leapfrog(p, timestep);
                    if (species_id == phase_surface_flux_species_id_ &&
                        surface_flux_phase <
                            phase_surface_flux_accumulators_.size()) {
                        for (std::size_t surface = 0;
                             surface < cfg_.phase_surface_flux.positions.size();
                             ++surface) {
                            const double position =
                                cfg_.phase_surface_flux.positions[surface];
                            if (old_position < position && p.x >= position) {
                                accumulate_phase_surface_crossing(
                                    chunk, surface, 0, p, sp);
                            } else if (old_position > position &&
                                       p.x <= position) {
                                accumulate_phase_surface_crossing(
                                    chunk, surface, 1, p, sp);
                            }
                        }
                    }
                    if (grid_.boundary() ==
                        Boundary::Periodic) {
                        p.x = std::fmod(
                            std::fmod(
                                p.x, grid_.length()) +
                                grid_.length(),
                            grid_.length());
                    } else if (
                        p.x < 0.0 ||
                        p.x > grid_.length()) {
                        const double old_longitudinal_energy =
                            0.5 * sp.mass() *
                            sp.config().weight * p.v * p.v;
                        const double crossing_longitudinal_energy =
                            0.5 * sp.mass() *
                            sp.config().weight *
                            p.v_half * p.v_half;
                        power.electric_work +=
                            crossing_longitudinal_energy -
                            old_longitudinal_energy;
                        const double speed_squared =
                            p.v_half * p.v_half +
                            p.velocity_y * p.velocity_y +
                            p.velocity_z * p.velocity_z;
                        const double represented_energy =
                            0.5 * sp.mass() *
                            sp.config().weight *
                            speed_squared;
                        const double particle_energy =
                            represented_energy /
                            sp.config().weight;
                        if (p.x < 0.0) {
                            ++loss.absorbed_left;
                            loss.kinetic_energy_left +=
                                represented_energy;
                            if (impact) {
                                accumulate_wall_impact(
                                    impact->left, species_id,
                                    particle_energy,
                                    represented_energy);
                            }
                        } else {
                            ++loss.absorbed_right;
                            loss.kinetic_energy_right +=
                                represented_energy;
                            if (impact) {
                                accumulate_wall_impact(
                                    impact->right, species_id,
                                    particle_energy,
                                    represented_energy);
                            }
                        }
                        p.alive = false;
                    }
                }
            });
        if (species_id == phase_surface_flux_species_id_ &&
            surface_flux_phase < phase_surface_flux_accumulators_.size()) {
            merge_phase_surface_flux_chunks(surface_flux_phase);
        }
        auto& total_loss =
            species_boundary_losses_[species_id];
        for (const auto& loss : chunk_losses) {
            total_loss.absorbed_left +=
                loss.absorbed_left;
            total_loss.absorbed_right +=
                loss.absorbed_right;
            total_loss.kinetic_energy_left +=
                loss.kinetic_energy_left;
            total_loss.kinetic_energy_right +=
                loss.kinetic_energy_right;
        }
        if (cfg_.wall_impact_spectrum.enabled) {
            auto& total = wall_impact_spectra_[species_id];
            for (const auto& chunk :
                 wall_impact_chunks_[species_id]) {
                add_wall_impact_side(total.left, chunk.left);
                add_wall_impact_side(total.right, chunk.right);
            }
        }
    }
    deposit_and_solve(time_ + cfg_.dt);
    for (std::size_t species_id = 0;
         species_id < species_.size(); ++species_id) {
        if (!species_due(species_id)) continue;
        auto& sp = species_[species_id];
        const double qm = sp.charge() / sp.mass();
        const double timestep = species_timestep(species_id);
        auto& particles = sp.particles();
        auto& chunk_power =
            power_transfer_chunks_[species_id];
        runtime_static_chunks(
            std::size_t{0}, particles.size(), cfg_.runtime,
            [&](std::size_t chunk, std::size_t begin,
                std::size_t end) {
                auto& power = chunk_power[chunk];
                for (std::size_t particle_id = begin;
                     particle_id < end; ++particle_id) {
                    auto& p = particles[particle_id];
                    if (!p.alive) continue;
                    const double old_velocity = p.v;
                    synchronize_leapfrog(
                        p, interpolate_electric(grid_, p.x),
                        qm, timestep);
                    power.electric_work +=
                        0.5 * sp.mass() *
                        sp.config().weight *
                        (p.v * p.v -
                         old_velocity * old_velocity);
                }
            });
        auto& total_power =
            species_power_transfer_[species_id];
        for (const auto& power : chunk_power) {
            total_power.electric_work +=
                power.electric_work;
        }
        if (!std::isfinite(total_power.electric_work)) {
            throw std::runtime_error(
                "1D species electric-work diagnostic "
                "became non-finite");
        }
    }
    update_phase_eedf_histories();
    begin_spatial_collision_step();
    if (cfg_.spatial_average.sampling_order ==
        SpatialAverageSamplingOrder1D::PreCollision) {
        accumulate_spatial_average(step_ + 1);
    }
    apply_collisions();
    ++step_;
    time_ += cfg_.dt;
    if (cfg_.spatial_average.sampling_order ==
        SpatialAverageSamplingOrder1D::PostCollision) {
        accumulate_spatial_average(step_);
    }
}

DiagnosticSample Simulation::sample() const {
    DiagnosticSample s;
    s.step = step_;
    s.time = time_;
    for (const auto& sp : species_) {
        s.kinetic_energy += sp.kinetic_energy();
        s.live_particles += sp.live_count();
    }
    for (std::size_t i = 0; i < grid_.nx(); ++i) {
        const double volume = grid_.node_volume(i);
        s.field_energy +=
            0.5 * cfg_.units.permittivity() *
            grid_.electric()[i] * grid_.electric()[i] * volume;
        s.charge_l1 += std::abs(grid_.rho()[i]) * volume;
    }
    s.total_energy = s.kinetic_energy + s.field_energy;
    if (grid_.boundary() == Boundary::Dirichlet) {
        s.phi_left = grid_.phi().front();
        s.phi_right = grid_.phi().back();
    }
    return s;
}

std::size_t Simulation::expected_spatial_average_samples() const {
    if (!cfg_.spatial_average.enabled) return 0;
    return 1 +
        (cfg_.spatial_average.end_step -
         cfg_.spatial_average.start_step) /
            cfg_.spatial_average.interval;
}

void Simulation::deposit_spatial_collision_energy(
    double position,
    std::size_t channel,
    double represented_energy_change) {
    if (!spatial_collision_step_active_) {
        return;
    }
    if (channel >= spatial_collision_energy_sums_.size() ||
        spatial_collision_energy_sums_[channel].size() != grid_.nx() ||
        !std::isfinite(represented_energy_change)) {
        throw std::logic_error(
            "spatial collision-energy storage does not match diagnostics");
    }
    const auto add = [&](std::size_t node, double shape, double volume) {
        const double density = represented_energy_change * shape / volume;
        spatial_collision_energy_sums_[channel][node] += density;
        if (spatial_collision_active_phase_ <
            spatial_collision_phase_energy_sums_.size()) {
            spatial_collision_phase_energy_sums_
                [spatial_collision_active_phase_][channel][node] += density;
        }
    };
    const double dx = grid_.dx();
    if (grid_.boundary() == Boundary::Periodic) {
        const double wrapped = std::fmod(
            std::fmod(position, grid_.length()) + grid_.length(),
            grid_.length());
        const double coordinate = wrapped / dx;
        const auto cell = static_cast<std::size_t>(std::floor(coordinate));
        const double fraction = coordinate - static_cast<double>(cell);
        add(cell % grid_.nx(), 1.0 - fraction, dx);
        add((cell + 1) % grid_.nx(), fraction, dx);
    } else {
        const double bounded = std::clamp(position, 0.0, grid_.length());
        const double coordinate = bounded / dx;
        const auto cell = static_cast<std::size_t>(std::min<double>(
            std::floor(coordinate), grid_.nx() - 2));
        const double fraction = coordinate - static_cast<double>(cell);
        add(cell, 1.0 - fraction, grid_.node_volume(cell));
        add(cell + 1, fraction, grid_.node_volume(cell + 1));
    }
}

void Simulation::deposit_spatial_collision_events(
    double position,
    std::size_t channel,
    double represented_events) {
    if (!spatial_collision_step_active_) return;
    if (channel >= spatial_collision_event_sums_.size() ||
        spatial_collision_event_sums_[channel].size() != grid_.nx() ||
        !std::isfinite(represented_events) || represented_events < 0.0) {
        throw std::logic_error(
            "spatial collision-event storage does not match diagnostics");
    }
    const auto add = [&](std::size_t node, double shape, double volume) {
        const double density = represented_events * shape / volume;
        spatial_collision_event_sums_[channel][node] += density;
        if (spatial_collision_active_phase_ <
            spatial_collision_phase_event_sums_.size()) {
            spatial_collision_phase_event_sums_
                [spatial_collision_active_phase_][channel][node] += density;
        }
    };
    const double dx = grid_.dx();
    if (grid_.boundary() == Boundary::Periodic) {
        const double wrapped = std::fmod(
            std::fmod(position, grid_.length()) + grid_.length(),
            grid_.length());
        const double coordinate = wrapped / dx;
        const auto cell = static_cast<std::size_t>(std::floor(coordinate));
        const double fraction = coordinate - static_cast<double>(cell);
        add(cell % grid_.nx(), 1.0 - fraction, dx);
        add((cell + 1) % grid_.nx(), fraction, dx);
    } else {
        const double bounded = std::clamp(position, 0.0, grid_.length());
        const double coordinate = bounded / dx;
        const auto cell = static_cast<std::size_t>(std::min<double>(
            std::floor(coordinate), grid_.nx() - 2));
        const double fraction = coordinate - static_cast<double>(cell);
        add(cell, 1.0 - fraction, grid_.node_volume(cell));
        add(cell + 1, fraction, grid_.node_volume(cell + 1));
    }
}

void Simulation::begin_spatial_collision_step() {
    const auto& average = cfg_.spatial_average;
    const std::size_t collision_step = step_ + 1;
    spatial_collision_step_active_ = false;
    spatial_collision_active_phase_ =
        spatial_collision_phase_energy_sums_.size();
    if (!average.enabled || spatial_collision_energy_sums_.empty() ||
        collision_step < average.start_step ||
        collision_step > average.end_step) {
        return;
    }
    spatial_collision_step_active_ = true;
    ++spatial_collision_steps_;
    if (!spatial_collision_phase_steps_.empty()) {
        const auto steps_per_cycle = static_cast<std::size_t>(std::llround(
            1.0 / (average.rf_frequency * cfg_.dt)));
        const auto step_in_cycle =
            (collision_step - average.start_step) % steps_per_cycle;
        spatial_collision_active_phase_ = step_in_cycle *
            spatial_collision_phase_steps_.size() / steps_per_cycle;
        ++spatial_collision_phase_steps_[spatial_collision_active_phase_];
    }
}

void Simulation::accumulate_wall_impact(
    WallImpactSideSpectrum1D& accumulator,
    std::size_t species_id,
    double particle_energy,
    double represented_energy) const {
    if (!cfg_.wall_impact_spectrum.enabled) return;
    if (species_id >= species_.size() ||
        !std::isfinite(particle_energy) || particle_energy < 0.0 ||
        !std::isfinite(represented_energy) || represented_energy < 0.0) {
        throw std::logic_error(
            "invalid 1D wall-impact spectrum sample");
    }
    const double energy_scale =
        cfg_.units.system == UnitSystem::SI
            ? ELEMENTARY_CHARGE_SI
            : 1.0;
    const double energy = particle_energy / energy_scale;
    const double represented_count = species_[species_id].weight();
    ++accumulator.macro_impacts;
    accumulator.represented_impacts += represented_count;
    accumulator.represented_kinetic_energy += represented_energy;
    if (energy >= cfg_.wall_impact_spectrum.energy_max) {
        ++accumulator.overflow_macro_impacts;
        accumulator.overflow_represented_impacts += represented_count;
        return;
    }
    const auto bin = std::min(
        cfg_.wall_impact_spectrum.energy_bins - 1,
        static_cast<std::size_t>(
            energy / cfg_.wall_impact_spectrum.energy_max *
            cfg_.wall_impact_spectrum.energy_bins));
    ++accumulator.macro_histogram[bin];
    accumulator.represented_histogram[bin] += represented_count;
}

void Simulation::accumulate_phase_eedf(std::size_t phase) {
    if (!cfg_.phase_eedf.enabled ||
        phase >= phase_eedf_accumulators_.size()) return;
    const auto& species = species_[phase_eedf_species_id_];
    const double energy_scale =
        cfg_.units.system == UnitSystem::SI ? ELEMENTARY_CHARGE_SI : 1.0;
    const double weight = species.weight();
    const auto& particles = species.particles();
    for (std::size_t particle_id = 0;
         particle_id < particles.size(); ++particle_id) {
        const auto& particle = particles[particle_id];
        if (!particle.alive) continue;
        const double velocity_squared =
            particle.v * particle.v +
            (cfg_.velocity_dimensions == 3
                 ? particle.velocity_y * particle.velocity_y +
                       particle.velocity_z * particle.velocity_z
                 : 0.0);
        const double energy =
            0.5 * species.mass() * velocity_squared / energy_scale;
        for (std::size_t region_id = 0;
             region_id < cfg_.phase_eedf.regions.size(); ++region_id) {
            const auto& configured = cfg_.phase_eedf.regions[region_id];
            if (particle.x < configured.x_min ||
                particle.x > configured.x_max) continue;
            auto& accumulator =
                phase_eedf_accumulators_[phase][region_id];
            ++accumulator.macro_observations;
            accumulator.represented_observations += weight;
            accumulator.weighted_energy_sum += weight * energy;
            accumulator.weighted_energy_squared_sum +=
                weight * energy * energy;
            accumulator.weighted_velocity_x_sum += weight * particle.v;
            accumulator.weighted_velocity_y_sum +=
                weight * particle.velocity_y;
            accumulator.weighted_velocity_z_sum +=
                weight * particle.velocity_z;
            accumulator.weighted_velocity_x_squared_sum +=
                weight * particle.v * particle.v;
            accumulator.weighted_velocity_y_squared_sum +=
                weight * particle.velocity_y * particle.velocity_y;
            accumulator.weighted_velocity_z_squared_sum +=
                weight * particle.velocity_z * particle.velocity_z;
            if (energy >= cfg_.phase_eedf.tail_threshold) {
                accumulator.tail_represented_observations += weight;
                if (particle.v >= 0.0) {
                    accumulator.tail_positive_x_represented_observations +=
                        weight;
                } else {
                    accumulator.tail_negative_x_represented_observations +=
                        weight;
                }
                accumulator.tail_weighted_velocity_x_sum +=
                    weight * particle.v;
                accumulator.tail_weighted_velocity_x_squared_sum +=
                    weight * particle.v * particle.v;
                accumulator.tail_weighted_transverse_velocity_squared_sum +=
                    weight * (particle.velocity_y * particle.velocity_y +
                              particle.velocity_z * particle.velocity_z);
                if (cfg_.phase_eedf.history_enabled) {
                    if (particle_id >=
                        phase_eedf_particle_histories_.size()) {
                        throw std::logic_error(
                            "phase EEDF particle history is misaligned");
                    }
                    const auto& history =
                        phase_eedf_particle_histories_[particle_id];
                    accumulator.tail_weighted_age_steps_sum +=
                        weight * static_cast<double>(history.age_steps);
                    accumulator.tail_weighted_energetic_steps_sum +=
                        weight * static_cast<double>(history.energetic_steps);
                    accumulator.tail_weighted_energetic_duty_fraction_sum +=
                        weight * (history.age_steps > 0
                            ? static_cast<double>(history.energetic_steps) /
                                  static_cast<double>(history.age_steps)
                            : 0.0);
                    accumulator
                        .tail_weighted_consecutive_energetic_steps_sum +=
                        weight * static_cast<double>(
                            history.consecutive_energetic_steps);
                    accumulator.tail_weighted_entries_sum +=
                        weight * static_cast<double>(history.tail_entries);
                    accumulator.tail_weighted_elastic_collisions_sum +=
                        weight * static_cast<double>(
                            history.elastic_collisions);
                    accumulator.tail_weighted_excitation_collisions_sum +=
                        weight * static_cast<double>(
                            history.excitation_collisions);
                    accumulator.tail_weighted_ionization_collisions_sum +=
                        weight * static_cast<double>(
                            history.ionization_collisions);
                    accumulator
                        .tail_weighted_charge_exchange_collisions_sum +=
                        weight * static_cast<double>(
                            history.charge_exchange_collisions);
                    accumulator.tail_weighted_bgk_collisions_sum +=
                        weight * static_cast<double>(history.bgk_collisions);
                    if (history.born_during_window) {
                        accumulator
                            .tail_born_during_window_represented_observations +=
                            weight;
                    }
                }
            }
            if (energy >= cfg_.phase_eedf.energy_max) {
                ++accumulator.overflow_macro_observations;
                accumulator.overflow_represented_observations += weight;
            } else {
                const auto bin = std::min(
                    cfg_.phase_eedf.energy_bins - 1,
                    static_cast<std::size_t>(
                        energy / cfg_.phase_eedf.energy_max *
                        cfg_.phase_eedf.energy_bins));
                accumulator.histogram[bin] += weight;
            }
        }
    }
}

bool Simulation::phase_eedf_history_active() const {
    const std::size_t sample_step = step_ + 1;
    return cfg_.phase_eedf.history_enabled &&
        sample_step >= cfg_.spatial_average.start_step &&
        sample_step <= cfg_.spatial_average.end_step;
}

std::size_t Simulation::phase_eedf_history_phase() const {
    if (!phase_eedf_history_active() ||
        phase_eedf_threshold_crossings_.empty()) {
        return phase_eedf_threshold_crossings_.size();
    }
    const auto& average = cfg_.spatial_average;
    const auto steps_per_cycle = static_cast<std::size_t>(std::llround(
        1.0 / (average.rf_frequency * cfg_.dt)));
    const std::size_t sample_step = step_ + 1;
    const std::size_t step_in_cycle =
        (sample_step - average.start_step) % steps_per_cycle;
    return std::min(
        phase_eedf_threshold_crossings_.size() - 1,
        step_in_cycle * phase_eedf_threshold_crossings_.size() /
            steps_per_cycle);
}

bool Simulation::phase_eedf_collision_state_energetic(
    const Particle& particle, const Species& species) const {
    const double longitudinal_velocity =
        cfg_.collision_velocity_sampling ==
                CollisionVelocitySampling1D::LeapfrogHalfStep
            ? particle.v_half
            : particle.v;
    const double velocity_squared =
        longitudinal_velocity * longitudinal_velocity +
        (cfg_.velocity_dimensions == 3
             ? particle.velocity_y * particle.velocity_y +
                   particle.velocity_z * particle.velocity_z
             : 0.0);
    const double energy_scale =
        cfg_.units.system == UnitSystem::SI ? ELEMENTARY_CHARGE_SI : 1.0;
    return 0.5 * species.mass() * velocity_squared / energy_scale >=
        cfg_.phase_eedf.tail_threshold;
}

void Simulation::update_phase_eedf_histories() {
    if (!phase_eedf_history_active() ||
        !species_due(phase_eedf_species_id_)) return;
    const auto& species = species_[phase_eedf_species_id_];
    const auto& particles = species.particles();
    if (phase_eedf_particle_histories_.size() != particles.size()) {
        throw std::logic_error(
            "phase EEDF particle history is misaligned");
    }
    if (phase_eedf_field_push_origin_energetic_.size() != particles.size()) {
        throw std::logic_error(
            "phase EEDF field-push origin state is misaligned");
    }
    const double energy_scale =
        cfg_.units.system == UnitSystem::SI ? ELEMENTARY_CHARGE_SI : 1.0;
    const std::size_t phase = phase_eedf_history_phase();
    if (phase >= phase_eedf_threshold_crossings_.size()) return;
    for (std::size_t particle_id = 0;
         particle_id < particles.size(); ++particle_id) {
        const auto& particle = particles[particle_id];
        if (!particle.alive) continue;
        auto& history = phase_eedf_particle_histories_[particle_id];
        const bool has_previous = history.age_steps > 0;
        const double velocity_squared =
            particle.v * particle.v +
            (cfg_.velocity_dimensions == 3
                 ? particle.velocity_y * particle.velocity_y +
                       particle.velocity_z * particle.velocity_z
                 : 0.0);
        const double energy =
            0.5 * species.mass() * velocity_squared / energy_scale;
        const bool energetic =
            energy >= cfg_.phase_eedf.tail_threshold;
        const bool field_push_energetic =
            phase_eedf_collision_state_energetic(particle, species);
        const bool field_push_origin_energetic =
            phase_eedf_field_push_origin_energetic_[particle_id] != 0;
        for (std::size_t region_id = 0;
             region_id < cfg_.phase_eedf.regions.size(); ++region_id) {
            const auto& configured = cfg_.phase_eedf.regions[region_id];
            if (particle.x < configured.x_min ||
                particle.x > configured.x_max) continue;
            auto& crossing =
                phase_eedf_threshold_crossings_[phase][region_id];
            ++crossing.electron_time_macro_observations;
            if (energetic) {
                ++crossing.energetic_time_macro_observations;
            }
            if (has_previous &&
                !history.energetic_previous_step && energetic) {
                ++crossing.interstep_promotions;
            } else if (has_previous &&
                       history.energetic_previous_step && !energetic) {
                ++crossing.interstep_demotions;
            }
            if (!field_push_origin_energetic && field_push_energetic) {
                ++crossing.field_push_promotions;
            } else if (field_push_origin_energetic &&
                       !field_push_energetic) {
                ++crossing.field_push_demotions;
            }
        }
        if (energetic) {
            ++history.energetic_steps;
            ++history.consecutive_energetic_steps;
            if (!history.energetic_previous_step) {
                ++history.tail_entries;
            }
        } else {
            history.consecutive_energetic_steps = 0;
        }
        history.energetic_previous_step = energetic;
        ++history.age_steps;
    }
}

void Simulation::add_phase_eedf_collision_transition(
    double position, CollisionProcessKind process,
    bool energetic_before, bool energetic_after) {
    if (energetic_before == energetic_after) return;
    const std::size_t phase = phase_eedf_history_phase();
    if (phase >= phase_eedf_threshold_crossings_.size()) return;
    std::size_t process_id = 0;
    switch (process) {
        case CollisionProcessKind::Elastic: process_id = 0; break;
        case CollisionProcessKind::Excitation: process_id = 1; break;
        case CollisionProcessKind::Ionization: process_id = 2; break;
        case CollisionProcessKind::ChargeExchange: process_id = 3; break;
        case CollisionProcessKind::Attachment: process_id = 4; break;
    }
    for (std::size_t region_id = 0;
         region_id < cfg_.phase_eedf.regions.size(); ++region_id) {
        const auto& configured = cfg_.phase_eedf.regions[region_id];
        if (position < configured.x_min || position > configured.x_max) {
            continue;
        }
        auto& crossing =
            phase_eedf_threshold_crossings_[phase][region_id];
        if (energetic_after) ++crossing.collision_promotions[process_id];
        else ++crossing.collision_demotions[process_id];
    }
}

void Simulation::add_phase_eedf_bgk_transition(
    double position, bool energetic_before, bool energetic_after) {
    if (energetic_before == energetic_after) return;
    const std::size_t phase = phase_eedf_history_phase();
    if (phase >= phase_eedf_threshold_crossings_.size()) return;
    for (std::size_t region_id = 0;
         region_id < cfg_.phase_eedf.regions.size(); ++region_id) {
        const auto& configured = cfg_.phase_eedf.regions[region_id];
        if (position < configured.x_min || position > configured.x_max) {
            continue;
        }
        auto& crossing =
            phase_eedf_threshold_crossings_[phase][region_id];
        if (energetic_after) ++crossing.collision_promotions[5];
        else ++crossing.collision_demotions[5];
    }
}

void Simulation::add_phase_eedf_birth(double position, bool energetic) {
    const std::size_t phase = phase_eedf_history_phase();
    if (phase >= phase_eedf_threshold_crossings_.size()) return;
    for (std::size_t region_id = 0;
         region_id < cfg_.phase_eedf.regions.size(); ++region_id) {
        const auto& configured = cfg_.phase_eedf.regions[region_id];
        if (position < configured.x_min || position > configured.x_max) {
            continue;
        }
        auto& crossing =
            phase_eedf_threshold_crossings_[phase][region_id];
        if (energetic) ++crossing.energetic_births;
        else ++crossing.subthreshold_births;
    }
}

void Simulation::add_phase_eedf_collision_history(
    std::size_t particle_id,
    CollisionProcessKind process,
    std::uint64_t count) {
    if (!phase_eedf_history_active() || count == 0) return;
    if (particle_id >= phase_eedf_particle_histories_.size()) {
        throw std::logic_error(
            "phase EEDF collision history is misaligned");
    }
    auto& history = phase_eedf_particle_histories_[particle_id];
    switch (process) {
        case CollisionProcessKind::Elastic:
            history.elastic_collisions += count;
            break;
        case CollisionProcessKind::Excitation:
            history.excitation_collisions += count;
            break;
        case CollisionProcessKind::Ionization:
            history.ionization_collisions += count;
            break;
        case CollisionProcessKind::ChargeExchange:
            history.charge_exchange_collisions += count;
            break;
        case CollisionProcessKind::Attachment:
            break;
    }
}

std::size_t Simulation::phase_surface_flux_phase(
    std::size_t sample_step) const {
    if (!cfg_.phase_surface_flux.enabled) {
        return phase_surface_flux_accumulators_.size();
    }
    const auto& average = cfg_.spatial_average;
    if (sample_step < average.start_step ||
        sample_step > average.end_step ||
        phase_surface_flux_accumulators_.empty()) {
        return phase_surface_flux_accumulators_.size();
    }
    const auto steps_per_cycle = static_cast<std::size_t>(std::llround(
        1.0 / (average.rf_frequency * cfg_.dt)));
    const std::size_t step_in_cycle =
        (sample_step - average.start_step) % steps_per_cycle;
    return std::min(
        phase_surface_flux_accumulators_.size() - 1,
        step_in_cycle * phase_surface_flux_accumulators_.size() /
            steps_per_cycle);
}

void Simulation::accumulate_phase_surface_crossing(
    std::size_t chunk, std::size_t surface, std::size_t direction,
    const Particle& particle, const Species& species) {
    if (chunk >= phase_surface_flux_chunks_.size() ||
        surface >= cfg_.phase_surface_flux.positions.size() ||
        direction >= 2) {
        throw std::logic_error("surface-flux crossing index differs");
    }
    auto& value = phase_surface_flux_chunks_[chunk][surface][direction];
    const double speed_squared =
        particle.v_half * particle.v_half +
        (cfg_.velocity_dimensions == 3
             ? particle.velocity_y * particle.velocity_y +
                   particle.velocity_z * particle.velocity_z
             : 0.0);
    const double particle_energy_joule =
        0.5 * species.mass() * speed_squared;
    const double energy_scale =
        cfg_.units.system == UnitSystem::SI ? ELEMENTARY_CHARGE_SI : 1.0;
    const double diagnostic_energy = particle_energy_joule / energy_scale;
    ++value.macro_crossings;
    value.represented_crossings += species.weight();
    value.represented_kinetic_energy +=
        species.weight() * particle_energy_joule;
    if (diagnostic_energy >= cfg_.phase_surface_flux.energy_max) {
        ++value.overflow_macro_crossings;
        value.overflow_represented_crossings += species.weight();
        return;
    }
    const auto bin = std::min(
        cfg_.phase_surface_flux.energy_bins - 1,
        static_cast<std::size_t>(
            diagnostic_energy / cfg_.phase_surface_flux.energy_max *
            cfg_.phase_surface_flux.energy_bins));
    value.represented_histogram[bin] += species.weight();
}

void Simulation::merge_phase_surface_flux_chunks(std::size_t phase) {
    if (phase >= phase_surface_flux_accumulators_.size()) return;
    for (auto& chunk : phase_surface_flux_chunks_) {
        for (std::size_t surface = 0;
             surface < chunk.size(); ++surface) {
            for (std::size_t direction = 0; direction < 2; ++direction) {
                add_surface_flux_accumulator(
                    phase_surface_flux_accumulators_[phase][surface][direction],
                    chunk[surface][direction]);
                clear_surface_flux_accumulator(chunk[surface][direction]);
            }
        }
    }
}

void Simulation::accumulate_spatial_average(std::size_t sample_step) {
    const auto& average = cfg_.spatial_average;
    if (!average.enabled ||
        sample_step < average.start_step ||
        sample_step > average.end_step ||
        (sample_step - average.start_step) %
                average.interval !=
            0) {
        return;
    }
    if (spatial_density_sums_.size() !=
            species_.size() ||
        spatial_density_scratch_.size() !=
            grid_.nx() ||
        spatial_kinetic_energy_sums_.size() != species_.size() ||
        spatial_kinetic_energy_scratch_.size() != grid_.nx() ||
        spatial_velocity_x_scratch_.size() != grid_.nx() ||
        spatial_velocity_y_scratch_.size() != grid_.nx() ||
        spatial_velocity_z_scratch_.size() != grid_.nx() ||
        spatial_potential_sums_.size() != grid_.nx() ||
        spatial_electric_sums_.size() != grid_.nx() ||
        spatial_electric_squared_sums_.size() != grid_.nx()) {
        throw std::logic_error(
            "spatial-average storage does not match simulation state");
    }
    SpatialPhaseBin1D* phase_bin = nullptr;
    std::size_t phase_id = spatial_phase_bins_.size();
    if (!spatial_phase_bins_.empty()) {
        const auto steps_per_cycle = static_cast<std::size_t>(std::llround(
            1.0 / (average.rf_frequency * cfg_.dt)));
        const auto samples_per_cycle = steps_per_cycle / average.interval;
        const auto sample_in_cycle =
            ((sample_step - average.start_step) / average.interval) %
            samples_per_cycle;
        phase_id =
            sample_in_cycle * spatial_phase_bins_.size() /
            samples_per_cycle;
        phase_bin = &spatial_phase_bins_[phase_id];
    }
    for (std::size_t species_id = 0;
         species_id < species_.size(); ++species_id) {
        species_[species_id].deposit_velocity_moments(
            grid_, spatial_density_scratch_,
            spatial_velocity_x_scratch_, spatial_velocity_y_scratch_,
            spatial_velocity_z_scratch_, spatial_kinetic_energy_scratch_);
        auto& sum = spatial_density_sums_[species_id];
        auto& energy_sum = spatial_kinetic_energy_sums_[species_id];
        for (std::size_t node = 0;
             node < grid_.nx(); ++node) {
            sum[node] += spatial_density_scratch_[node];
            energy_sum[node] += spatial_kinetic_energy_scratch_[node];
            if (phase_bin != nullptr) {
                phase_bin->density[species_id][node] +=
                    spatial_density_scratch_[node];
                phase_bin->velocity_x_density[species_id][node] +=
                    spatial_velocity_x_scratch_[node];
                phase_bin->velocity_y_density[species_id][node] +=
                    spatial_velocity_y_scratch_[node];
                phase_bin->velocity_z_density[species_id][node] +=
                    spatial_velocity_z_scratch_[node];
                phase_bin->kinetic_energy_density[species_id][node] +=
                    spatial_kinetic_energy_scratch_[node];
            }
        }
    }
    for (std::size_t node = 0; node < grid_.nx(); ++node) {
        spatial_potential_sums_[node] += grid_.phi()[node];
        spatial_electric_sums_[node] += grid_.electric()[node];
        spatial_electric_squared_sums_[node] +=
            grid_.electric()[node] * grid_.electric()[node];
        if (phase_bin != nullptr) {
            phase_bin->potential[node] += grid_.phi()[node];
            phase_bin->electric[node] += grid_.electric()[node];
            phase_bin->electric_squared[node] +=
                grid_.electric()[node] * grid_.electric()[node];
        }
    }
    if (phase_bin != nullptr) ++phase_bin->samples;
    if (phase_id < spatial_phase_bins_.size()) {
        accumulate_phase_eedf(phase_id);
    }
    ++spatial_average_samples_;
    ++spatial_moment_samples_;
}

void Simulation::write_wall_impact_spectrum() const {
    if (!cfg_.wall_impact_spectrum.enabled) return;
    if (wall_impact_spectra_.size() != species_.size()) {
        throw std::logic_error(
            "1D wall-impact spectra do not match species");
    }
    const auto output_dir = std::filesystem::path(cfg_.output_dir);
    std::ofstream histogram(
        output_dir / "wall_impact_spectrum.csv");
    std::ofstream summary(
        output_dir / "wall_impact_spectrum_summary.csv");
    if (!histogram || !summary) {
        throw std::runtime_error(
            "cannot open 1D wall-impact spectrum output");
    }
    const bool si = cfg_.units.system == UnitSystem::SI;
    const char* energy_name = si
        ? "impact_energy_eV"
        : "impact_energy_normalized";
    const char* represented_energy_name = si
        ? "represented_kinetic_energy_J_m-2"
        : "represented_kinetic_energy_normalized";
    histogram << "origin_step,species_id,species,electrode,energy_bin,"
              << energy_name
              << ",macro_count,represented_count,probability_density\n"
              << std::setprecision(17);
    summary << "origin_step,species_id,species,electrode,macro_impacts,"
               "represented_impacts,overflow_macro_impacts,"
               "overflow_represented_impacts,overflow_fraction,"
            << represented_energy_name
            << ",boundary_delta_macro_impacts,boundary_delta_"
            << represented_energy_name
            << ",count_closure,energy_closure_residual\n"
            << std::setprecision(17);
    const double bin_width =
        cfg_.wall_impact_spectrum.energy_max /
        static_cast<double>(
            cfg_.wall_impact_spectrum.energy_bins);
    for (std::size_t species_id = 0;
         species_id < species_.size(); ++species_id) {
        const auto& spectrum = wall_impact_spectra_[species_id];
        const auto& loss = species_boundary_losses_[species_id];
        for (const bool left : {true, false}) {
            const auto& side = left
                ? spectrum.left : spectrum.right;
            const auto boundary_count = left
                ? loss.absorbed_left -
                      spectrum.baseline_loss.absorbed_left
                : loss.absorbed_right -
                      spectrum.baseline_loss.absorbed_right;
            const double boundary_energy = left
                ? loss.kinetic_energy_left -
                      spectrum.baseline_loss.kinetic_energy_left
                : loss.kinetic_energy_right -
                      spectrum.baseline_loss.kinetic_energy_right;
            if (boundary_count != side.macro_impacts) {
                throw std::runtime_error(
                    "1D wall-impact count does not close against "
                    "the boundary-loss ledger");
            }
            const double residual =
                side.represented_kinetic_energy - boundary_energy;
            const double tolerance =
                256.0 * std::numeric_limits<double>::epsilon() *
                std::max({std::numeric_limits<double>::min(),
                    std::abs(side.represented_kinetic_energy),
                    std::abs(boundary_energy)});
            if (std::abs(residual) > tolerance) {
                throw std::runtime_error(
                    "1D wall-impact energy does not close against "
                    "the boundary-loss ledger");
            }
            std::uint64_t binned_macro = 0;
            double binned_represented = 0.0;
            for (std::size_t bin = 0;
                 bin < side.macro_histogram.size(); ++bin) {
                binned_macro += side.macro_histogram[bin];
                binned_represented +=
                    side.represented_histogram[bin];
                histogram << wall_impact_origin_step_ << ','
                    << species_id << ','
                    << csv_cell(species_[species_id].name()) << ','
                    << (left ? "left" : "right") << ','
                    << bin << ','
                    << (static_cast<double>(bin) + 0.5) * bin_width
                    << ',' << side.macro_histogram[bin]
                    << ',' << side.represented_histogram[bin]
                    << ','
                    << (side.represented_impacts > 0.0
                            ? side.represented_histogram[bin] /
                                  side.represented_impacts /
                                  bin_width
                            : 0.0)
                    << '\n';
            }
            const double count_residual =
                binned_represented +
                side.overflow_represented_impacts -
                side.represented_impacts;
            const double count_tolerance =
                256.0 * std::numeric_limits<double>::epsilon() *
                std::max(1.0, side.represented_impacts);
            if (binned_macro + side.overflow_macro_impacts !=
                    side.macro_impacts ||
                std::abs(count_residual) > count_tolerance) {
                throw std::runtime_error(
                    "1D wall-impact histogram does not preserve "
                    "its impact count");
            }
            summary << wall_impact_origin_step_ << ','
                << species_id << ','
                << csv_cell(species_[species_id].name()) << ','
                << (left ? "left" : "right") << ','
                << side.macro_impacts << ','
                << side.represented_impacts << ','
                << side.overflow_macro_impacts << ','
                << side.overflow_represented_impacts << ','
                << (side.represented_impacts > 0.0
                        ? side.overflow_represented_impacts /
                              side.represented_impacts
                        : 0.0)
                << ',' << side.represented_kinetic_energy
                << ',' << boundary_count
                << ',' << boundary_energy
                << ",1," << residual << '\n';
        }
    }
    require_stream(
        histogram,
        "failed while writing 1D wall-impact histogram");
    require_stream(
        summary,
        "failed while writing 1D wall-impact summary");
}

void Simulation::write_spatial_average() const {
    if (!cfg_.spatial_average.enabled) return;
    const auto output_dir =
        std::filesystem::path(cfg_.output_dir);
    std::filesystem::create_directories(output_dir);
    const bool si =
        cfg_.units.system == UnitSystem::SI;
    const bool moments_complete =
        spatial_moment_samples_ == spatial_average_samples_;
    std::ofstream profile(
        output_dir / "spatial_average.csv");
    if (!profile) {
        throw std::runtime_error(
            "cannot open 1D spatial-average profile output");
    }
    profile
        << "species_id,species,node,"
        << (si ? "x_m" : "x_normalized") << ','
        << (si
                ? "number_density_mean_m-3"
                : "number_density_mean_normalized")
        << '\n';
    profile << std::setprecision(17);
    if (spatial_average_samples_ > 0) {
        for (std::size_t species_id = 0;
             species_id < species_.size(); ++species_id) {
            for (std::size_t node = 0;
                 node < grid_.nx(); ++node) {
                profile
                    << species_id << ','
                    << csv_cell(species_[species_id].name())
                    << ',' << node << ','
                    << grid_.node_x(node) << ','
                    << spatial_density_sums_[species_id][node] /
                           static_cast<double>(
                               spatial_average_samples_)
                    << '\n';
            }
        }
    }
    require_stream(
        profile,
        "failed while writing 1D spatial-average profile");

    std::ofstream energy_profile(
        output_dir / "spatial_kinetic_energy.csv");
    if (!energy_profile) {
        throw std::runtime_error(
            "cannot open 1D spatial kinetic-energy output");
    }
    energy_profile << "species_id,species,node,"
                   << (si ? "x_m" : "x_normalized") << ','
                   << (si ? "mean_kinetic_energy_eV"
                          : "mean_kinetic_energy_normalized") << ','
                   << (si ? "effective_kinetic_temperature_eV"
                          : "effective_kinetic_temperature_normalized")
                   << '\n' << std::setprecision(17);
    if (spatial_moment_samples_ > 0 && moments_complete) {
        const double energy_scale = si ? ELEMENTARY_CHARGE_SI : 1.0;
        const double dimensions =
            static_cast<double>(cfg_.velocity_dimensions);
        for (std::size_t species_id = 0;
             species_id < species_.size(); ++species_id) {
            for (std::size_t node = 0; node < grid_.nx(); ++node) {
                const double density_sum =
                    spatial_density_sums_[species_id][node];
                const double mean_energy = density_sum > 0.0
                    ? spatial_kinetic_energy_sums_[species_id][node] /
                          density_sum / energy_scale
                    : 0.0;
                energy_profile << species_id << ','
                    << csv_cell(species_[species_id].name()) << ','
                    << node << ',' << grid_.node_x(node) << ','
                    << mean_energy << ','
                    << 2.0 * mean_energy / dimensions << '\n';
            }
        }
    }
    require_stream(energy_profile,
        "failed while writing 1D spatial kinetic-energy output");

    std::ofstream field_profile(
        output_dir / "spatial_field_average.csv");
    if (!field_profile) {
        throw std::runtime_error(
            "cannot open 1D spatial field-average output");
    }
    field_profile << "node," << (si ? "x_m" : "x_normalized")
                  << ',' << (si ? "potential_mean_V"
                                   : "potential_mean_normalized")
                  << ',' << (si ? "electric_field_mean_V_m"
                                   : "electric_field_mean_normalized")
                  << ',' << (si ? "electric_field_rms_V_m"
                                   : "electric_field_rms_normalized")
                  << '\n' << std::setprecision(17);
    if (spatial_moment_samples_ > 0 && moments_complete) {
        const double samples = static_cast<double>(spatial_moment_samples_);
        for (std::size_t node = 0; node < grid_.nx(); ++node) {
            field_profile << node << ',' << grid_.node_x(node) << ','
                << spatial_potential_sums_[node] / samples << ','
                << spatial_electric_sums_[node] / samples << ','
                << std::sqrt(spatial_electric_squared_sums_[node] / samples)
                << '\n';
        }
    }
    require_stream(field_profile,
        "failed while writing 1D spatial field-average output");

    if (!spatial_phase_bins_.empty()) {
        std::ofstream phase_moments(
            output_dir / "spatial_phase_moments.csv");
        if (!phase_moments) {
            throw std::runtime_error(
                "cannot open 1D phase-resolved moment output");
        }
        phase_moments
            << "phase_bin,phase_fraction,samples,species_id,species,node,"
            << (si ? "x_m" : "x_normalized") << ','
            << (si ? "number_density_mean_m-3"
                   : "number_density_mean_normalized")
            << ",mean_velocity_x,mean_velocity_y,mean_velocity_z,"
            << (si ? "mean_kinetic_energy_eV,drift_separated_temperature_eV"
                   : "mean_kinetic_energy_normalized,drift_separated_temperature_normalized")
            << '\n' << std::setprecision(17);
        const double energy_scale = si ? ELEMENTARY_CHARGE_SI : 1.0;
        const double dimensions = static_cast<double>(cfg_.velocity_dimensions);
        for (std::size_t phase = 0; phase < spatial_phase_bins_.size(); ++phase) {
            const auto& bin = spatial_phase_bins_[phase];
            for (std::size_t species_id = 0; species_id < species_.size(); ++species_id) {
                for (std::size_t node = 0; node < grid_.nx(); ++node) {
                    const double density_sum = bin.density[species_id][node];
                    const double ux = density_sum > 0.0
                        ? bin.velocity_x_density[species_id][node] / density_sum : 0.0;
                    const double uy = density_sum > 0.0
                        ? bin.velocity_y_density[species_id][node] / density_sum : 0.0;
                    const double uz = density_sum > 0.0
                        ? bin.velocity_z_density[species_id][node] / density_sum : 0.0;
                    const double mean_energy = density_sum > 0.0
                        ? bin.kinetic_energy_density[species_id][node] /
                              density_sum / energy_scale : 0.0;
                    const double drift_energy = 0.5 * species_[species_id].mass() *
                        (ux * ux + uy * uy + uz * uz) / energy_scale;
                    const double temperature = 2.0 / dimensions *
                        std::max(0.0, mean_energy - drift_energy);
                    phase_moments << phase << ','
                        << (static_cast<double>(phase) + 0.5) /
                               static_cast<double>(spatial_phase_bins_.size())
                        << ',' << bin.samples << ',' << species_id << ','
                        << csv_cell(species_[species_id].name()) << ','
                        << node << ',' << grid_.node_x(node) << ','
                        << (bin.samples > 0 ? density_sum /
                               static_cast<double>(bin.samples) : 0.0)
                        << ',' << ux << ',' << uy << ',' << uz << ','
                        << mean_energy << ',' << temperature << '\n';
                }
            }
        }
        require_stream(phase_moments,
            "failed while writing 1D phase-resolved moment output");

        std::ofstream phase_fields(
            output_dir / "spatial_phase_fields.csv");
        if (!phase_fields) {
            throw std::runtime_error(
                "cannot open 1D phase-resolved field output");
        }
        phase_fields << "phase_bin,phase_fraction,samples,node,"
            << (si ? "x_m,potential_mean_V,electric_field_mean_V_m,electric_field_rms_V_m"
                   : "x_normalized,potential_mean_normalized,electric_field_mean_normalized,electric_field_rms_normalized")
            << '\n' << std::setprecision(17);
        for (std::size_t phase = 0; phase < spatial_phase_bins_.size(); ++phase) {
            const auto& bin = spatial_phase_bins_[phase];
            const double samples = static_cast<double>(bin.samples);
            for (std::size_t node = 0; node < grid_.nx(); ++node) {
                phase_fields << phase << ','
                    << (static_cast<double>(phase) + 0.5) /
                           static_cast<double>(spatial_phase_bins_.size())
                    << ',' << bin.samples << ',' << node << ','
                    << grid_.node_x(node) << ','
                    << (bin.samples > 0 ? bin.potential[node] / samples : 0.0) << ','
                    << (bin.samples > 0 ? bin.electric[node] / samples : 0.0) << ','
                    << (bin.samples > 0 ? std::sqrt(
                           bin.electric_squared[node] / samples) : 0.0) << '\n';
            }
        }
        require_stream(phase_fields,
            "failed while writing 1D phase-resolved field output");
    }

    std::ofstream collision_power(
        output_dir / "spatial_collision_power.csv");
    if (!collision_power) {
        throw std::runtime_error(
            "cannot open 1D spatial collision-power output");
    }
    collision_power
        << "channel_id,channel,node,"
        << (si ? "x_m" : "x_normalized")
        << ",timesteps,"
        << (si ? "duration_s" : "duration_normalized") << ','
        << (si ? "energy_density_sum_J_m-3,mean_power_density_W_m-3"
               : "energy_density_sum_normalized,mean_power_density_normalized")
        << '\n' << std::setprecision(17);
    const double collision_duration =
        static_cast<double>(spatial_collision_steps_) * cfg_.dt;
    for (std::size_t channel = 0;
         channel < spatial_collision_energy_sums_.size(); ++channel) {
        for (std::size_t node = 0; node < grid_.nx(); ++node) {
            const double energy_density =
                spatial_collision_energy_sums_[channel][node];
            collision_power << channel << ','
                << csv_cell(collision_totals_.channel_names[channel]) << ','
                << node << ',' << grid_.node_x(node) << ','
                << spatial_collision_steps_ << ',' << collision_duration << ','
                << energy_density << ','
                << (collision_duration > 0.0
                        ? energy_density / collision_duration
                        : 0.0)
                << '\n';
        }
    }
    require_stream(collision_power,
        "failed while writing 1D spatial collision-power output");

    std::ofstream collision_rate(
        output_dir / "spatial_collision_rate.csv");
    if (!collision_rate) {
        throw std::runtime_error(
            "cannot open 1D spatial collision-rate output");
    }
    collision_rate << "channel_id,channel,node,"
        << (si ? "x_m" : "x_normalized") << ",timesteps,"
        << (si ? "duration_s" : "duration_normalized") << ','
        << (si ? "represented_event_density_sum_m-3,mean_event_rate_m-3_s-1"
               : "represented_event_density_sum_normalized,mean_event_rate_normalized")
        << '\n' << std::setprecision(17);
    for (std::size_t channel = 0;
         channel < spatial_collision_event_sums_.size(); ++channel) {
        for (std::size_t node = 0; node < grid_.nx(); ++node) {
            const double event_density =
                spatial_collision_event_sums_[channel][node];
            collision_rate << channel << ','
                << csv_cell(collision_totals_.channel_names[channel]) << ','
                << node << ',' << grid_.node_x(node) << ','
                << spatial_collision_steps_ << ',' << collision_duration << ','
                << event_density << ','
                << (collision_duration > 0.0
                        ? event_density / collision_duration : 0.0)
                << '\n';
        }
    }
    require_stream(collision_rate,
        "failed while writing 1D spatial collision-rate output");

    if (!spatial_collision_phase_energy_sums_.empty()) {
        std::ofstream phase_collision_power(
            output_dir / "spatial_phase_collision_power.csv");
        if (!phase_collision_power) {
            throw std::runtime_error(
                "cannot open 1D phase-resolved collision-power output");
        }
        phase_collision_power
            << "phase_bin,phase_fraction,timesteps,"
            << (si ? "duration_s" : "duration_normalized")
            << ",channel_id,"
               "channel,node,"
            << (si ? "x_m" : "x_normalized") << ','
            << (si ? "energy_density_sum_J_m-3,mean_power_density_W_m-3"
                   : "energy_density_sum_normalized,mean_power_density_normalized")
            << '\n' << std::setprecision(17);
        for (std::size_t phase = 0;
             phase < spatial_collision_phase_energy_sums_.size(); ++phase) {
            const double duration =
                static_cast<double>(spatial_collision_phase_steps_[phase]) *
                cfg_.dt;
            for (std::size_t channel = 0;
                 channel < spatial_collision_energy_sums_.size(); ++channel) {
                for (std::size_t node = 0; node < grid_.nx(); ++node) {
                    const double energy_density =
                        spatial_collision_phase_energy_sums_[phase]
                            [channel][node];
                    phase_collision_power << phase << ','
                        << (static_cast<double>(phase) + 0.5) /
                               static_cast<double>(
                                   spatial_collision_phase_energy_sums_.size())
                        << ',' << spatial_collision_phase_steps_[phase] << ','
                        << duration << ',' << channel << ','
                        << csv_cell(collision_totals_.channel_names[channel])
                        << ',' << node << ',' << grid_.node_x(node) << ','
                        << energy_density << ','
                        << (duration > 0.0 ? energy_density / duration : 0.0)
                        << '\n';
                }
            }
        }
        require_stream(phase_collision_power,
            "failed while writing 1D phase-resolved collision-power output");

        std::ofstream phase_collision_rate(
            output_dir / "spatial_phase_collision_rate.csv");
        if (!phase_collision_rate) {
            throw std::runtime_error(
                "cannot open 1D phase-resolved collision-rate output");
        }
        phase_collision_rate
            << "phase_bin,phase_fraction,timesteps,"
            << (si ? "duration_s" : "duration_normalized")
            << ",channel_id,channel,node,"
            << (si ? "x_m" : "x_normalized") << ','
            << (si ? "represented_event_density_sum_m-3,mean_event_rate_m-3_s-1"
                   : "represented_event_density_sum_normalized,mean_event_rate_normalized")
            << '\n' << std::setprecision(17);
        for (std::size_t phase = 0;
             phase < spatial_collision_phase_event_sums_.size(); ++phase) {
            const double duration =
                static_cast<double>(spatial_collision_phase_steps_[phase]) *
                cfg_.dt;
            for (std::size_t channel = 0;
                 channel < spatial_collision_event_sums_.size(); ++channel) {
                for (std::size_t node = 0; node < grid_.nx(); ++node) {
                    const double event_density =
                        spatial_collision_phase_event_sums_[phase]
                            [channel][node];
                    phase_collision_rate << phase << ','
                        << (static_cast<double>(phase) + 0.5) /
                               static_cast<double>(
                                   spatial_collision_phase_event_sums_.size())
                        << ',' << spatial_collision_phase_steps_[phase] << ','
                        << duration << ',' << channel << ','
                        << csv_cell(collision_totals_.channel_names[channel])
                        << ',' << node << ',' << grid_.node_x(node) << ','
                        << event_density << ','
                        << (duration > 0.0 ? event_density / duration : 0.0)
                        << '\n';
                }
            }
        }
        require_stream(phase_collision_rate,
            "failed while writing 1D phase-resolved collision-rate output");
    }

    if (cfg_.phase_eedf.enabled) {
        std::ofstream histogram(output_dir / "phase_eedf.csv");
        std::ofstream moments(output_dir / "phase_eedf_moments.csv");
        if (!histogram || !moments) {
            throw std::runtime_error("cannot open 1D phase EEDF output");
        }
        const char* energy_name = si ? "energy_eV" : "energy_normalized";
        histogram << "phase_bin,phase_fraction,region_id,region,x_min,x_max,"
                  << "energy_bin," << energy_name
                  << ",represented_count,probability_density\n"
                  << std::setprecision(17);
        moments << "phase_bin,phase_fraction,region_id,region,x_min,x_max,"
                   "macro_observations,represented_observations,"
                   "overflow_fraction,mean_energy,energy_standard_deviation,"
                   "mean_velocity_x,mean_velocity_y,mean_velocity_z,"
                   "drift_separated_temperature,temperature_x,temperature_y,"
                   "temperature_z,tail_threshold,tail_represented_observations,"
                   "tail_positive_x_fraction,tail_negative_x_fraction,"
                   "tail_directional_population_imbalance,tail_mean_velocity_x,"
                   "tail_longitudinal_energy_fraction,tail_mean_energy,"
                   "history_enabled,tail_mean_age_steps,"
                   "tail_mean_energetic_steps,"
                   "tail_mean_energetic_duty_fraction,"
                   "tail_mean_consecutive_energetic_steps,tail_mean_entries,"
                   "tail_mean_elastic_collisions,"
                   "tail_mean_excitation_collisions,"
                   "tail_mean_ionization_collisions,"
                   "tail_mean_charge_exchange_collisions,"
                   "tail_mean_bgk_collisions,"
                   "tail_born_during_window_fraction\n"
                << std::setprecision(17);
        const double bin_width = cfg_.phase_eedf.energy_max /
            static_cast<double>(cfg_.phase_eedf.energy_bins);
        const double energy_scale = si ? ELEMENTARY_CHARGE_SI : 1.0;
        const double dimensions =
            static_cast<double>(cfg_.velocity_dimensions);
        const auto& target = species_[phase_eedf_species_id_];
        for (std::size_t phase = 0;
             phase < phase_eedf_accumulators_.size(); ++phase) {
            const double phase_fraction =
                (static_cast<double>(phase) + 0.5) /
                static_cast<double>(phase_eedf_accumulators_.size());
            for (std::size_t region_id = 0;
                 region_id < cfg_.phase_eedf.regions.size(); ++region_id) {
                const auto& region = cfg_.phase_eedf.regions[region_id];
                const auto& accumulator =
                    phase_eedf_accumulators_[phase][region_id];
                for (std::size_t bin = 0;
                     bin < accumulator.histogram.size(); ++bin) {
                    histogram << phase << ',' << phase_fraction << ','
                        << region_id << ',' << csv_cell(region.name) << ','
                        << region.x_min << ',' << region.x_max << ',' << bin
                        << ',' << (static_cast<double>(bin) + 0.5) * bin_width
                        << ',' << accumulator.histogram[bin] << ','
                        << (accumulator.represented_observations > 0.0
                                ? accumulator.histogram[bin] /
                                      accumulator.represented_observations /
                                      bin_width
                                : 0.0)
                        << '\n';
                }
                const double count = accumulator.represented_observations;
                const double mean_energy = count > 0.0
                    ? accumulator.weighted_energy_sum / count : 0.0;
                const double energy_variance = count > 0.0
                    ? std::max(0.0,
                        accumulator.weighted_energy_squared_sum / count -
                        mean_energy * mean_energy)
                    : 0.0;
                const double ux = count > 0.0
                    ? accumulator.weighted_velocity_x_sum / count : 0.0;
                const double uy = count > 0.0
                    ? accumulator.weighted_velocity_y_sum / count : 0.0;
                const double uz = count > 0.0
                    ? accumulator.weighted_velocity_z_sum / count : 0.0;
                const double vx2 = count > 0.0
                    ? accumulator.weighted_velocity_x_squared_sum / count : 0.0;
                const double vy2 = count > 0.0
                    ? accumulator.weighted_velocity_y_squared_sum / count : 0.0;
                const double vz2 = count > 0.0
                    ? accumulator.weighted_velocity_z_squared_sum / count : 0.0;
                const double drift_energy = 0.5 * target.mass() *
                    (ux * ux + uy * uy + uz * uz) / energy_scale;
                const double tail_count =
                    accumulator.tail_represented_observations;
                const double tail_velocity_squared =
                    accumulator.tail_weighted_velocity_x_squared_sum +
                    accumulator.tail_weighted_transverse_velocity_squared_sum;
                moments << phase << ',' << phase_fraction << ',' << region_id
                    << ',' << csv_cell(region.name) << ',' << region.x_min
                    << ',' << region.x_max << ','
                    << accumulator.macro_observations << ',' << count << ','
                    << (count > 0.0
                            ? accumulator.overflow_represented_observations /
                                  count
                            : 0.0)
                    << ',' << mean_energy << ',' << std::sqrt(energy_variance)
                    << ',' << ux << ',' << uy << ',' << uz << ','
                    << 2.0 / dimensions *
                           std::max(0.0, mean_energy - drift_energy) << ','
                    << target.mass() / energy_scale *
                           std::max(0.0, vx2 - ux * ux) << ','
                    << target.mass() / energy_scale *
                           std::max(0.0, vy2 - uy * uy) << ','
                    << target.mass() / energy_scale *
                           std::max(0.0, vz2 - uz * uz) << ','
                    << cfg_.phase_eedf.tail_threshold << ',' << tail_count << ','
                    << (tail_count > 0.0
                            ? accumulator.tail_positive_x_represented_observations /
                                  tail_count : 0.0) << ','
                    << (tail_count > 0.0
                            ? accumulator.tail_negative_x_represented_observations /
                                  tail_count : 0.0) << ','
                    << (tail_count > 0.0
                            ? (accumulator.tail_positive_x_represented_observations -
                               accumulator.tail_negative_x_represented_observations) /
                                  tail_count : 0.0) << ','
                    << (tail_count > 0.0
                            ? accumulator.tail_weighted_velocity_x_sum /
                                  tail_count : 0.0) << ','
                    << (tail_velocity_squared > 0.0
                            ? accumulator.tail_weighted_velocity_x_squared_sum /
                                  tail_velocity_squared : 0.0) << ','
                    << (tail_count > 0.0
                            ? 0.5 * target.mass() * tail_velocity_squared /
                                  tail_count / energy_scale : 0.0) << ','
                    << (cfg_.phase_eedf.history_enabled ? 1 : 0) << ','
                    << (tail_count > 0.0
                            ? accumulator.tail_weighted_age_steps_sum /
                                  tail_count : 0.0) << ','
                    << (tail_count > 0.0
                            ? accumulator.tail_weighted_energetic_steps_sum /
                                  tail_count : 0.0) << ','
                    << (tail_count > 0.0
                            ? accumulator
                                  .tail_weighted_energetic_duty_fraction_sum /
                                  tail_count : 0.0) << ','
                    << (tail_count > 0.0
                            ? accumulator
                                  .tail_weighted_consecutive_energetic_steps_sum /
                                  tail_count : 0.0) << ','
                    << (tail_count > 0.0
                            ? accumulator.tail_weighted_entries_sum /
                                  tail_count : 0.0) << ','
                    << (tail_count > 0.0
                            ? accumulator
                                  .tail_weighted_elastic_collisions_sum /
                                  tail_count : 0.0) << ','
                    << (tail_count > 0.0
                            ? accumulator
                                  .tail_weighted_excitation_collisions_sum /
                                  tail_count : 0.0) << ','
                    << (tail_count > 0.0
                            ? accumulator
                                  .tail_weighted_ionization_collisions_sum /
                                  tail_count : 0.0) << ','
                    << (tail_count > 0.0
                            ? accumulator
                                  .tail_weighted_charge_exchange_collisions_sum /
                                  tail_count : 0.0) << ','
                    << (tail_count > 0.0
                            ? accumulator.tail_weighted_bgk_collisions_sum /
                                  tail_count : 0.0) << ','
                    << (tail_count > 0.0
                            ? accumulator
                                  .tail_born_during_window_represented_observations /
                                  tail_count : 0.0)
                    << '\n';
            }
        }
        require_stream(histogram, "failed while writing phase EEDF histogram");
        require_stream(moments, "failed while writing phase EEDF moments");
        if (cfg_.phase_eedf.history_enabled) {
            std::ofstream crossings(
                output_dir / "phase_eedf_threshold_crossings.csv");
            if (!crossings) {
                throw std::runtime_error(
                    "cannot open phase EEDF threshold-crossing output");
            }
            crossings
                << "phase_bin,phase_fraction,region_id,region,x_min,x_max,"
                   "electron_time_macro_observations,"
                   "energetic_time_macro_observations,energetic_fraction,"
                   "interstep_promotions,interstep_demotions,"
                   "interstep_promotions_per_million_electron_steps,"
                   "interstep_demotions_per_million_electron_steps,"
                   "field_push_promotions,field_push_demotions,"
                   "field_push_promotions_per_million_electron_steps,"
                   "field_push_demotions_per_million_electron_steps,"
                   "elastic_collision_promotions,elastic_collision_demotions,"
                   "excitation_collision_promotions,"
                   "excitation_collision_demotions,"
                   "ionization_collision_promotions,"
                   "ionization_collision_demotions,"
                   "charge_exchange_collision_promotions,"
                   "charge_exchange_collision_demotions,"
                   "attachment_collision_promotions,"
                   "attachment_collision_demotions,"
                   "bgk_collision_promotions,bgk_collision_demotions,"
                   "energetic_births,subthreshold_births\n"
                << std::setprecision(17);
            for (std::size_t phase = 0;
                 phase < phase_eedf_threshold_crossings_.size(); ++phase) {
                const double phase_fraction =
                    (static_cast<double>(phase) + 0.5) /
                    static_cast<double>(
                        phase_eedf_threshold_crossings_.size());
                for (std::size_t region_id = 0;
                     region_id < cfg_.phase_eedf.regions.size(); ++region_id) {
                    const auto& region = cfg_.phase_eedf.regions[region_id];
                    const auto& value =
                        phase_eedf_threshold_crossings_[phase][region_id];
                    const double observations = static_cast<double>(
                        value.electron_time_macro_observations);
                    crossings << phase << ',' << phase_fraction << ','
                        << region_id << ',' << csv_cell(region.name) << ','
                        << region.x_min << ',' << region.x_max << ','
                        << value.electron_time_macro_observations << ','
                        << value.energetic_time_macro_observations << ','
                        << (observations > 0.0
                                ? static_cast<double>(
                                      value.energetic_time_macro_observations) /
                                      observations
                                : 0.0) << ','
                        << value.interstep_promotions << ','
                        << value.interstep_demotions << ','
                        << (observations > 0.0
                                ? 1.0e6 * static_cast<double>(
                                      value.interstep_promotions) /
                                      observations
                                : 0.0) << ','
                        << (observations > 0.0
                                ? 1.0e6 * static_cast<double>(
                                      value.interstep_demotions) /
                                      observations
                                : 0.0) << ','
                        << value.field_push_promotions << ','
                        << value.field_push_demotions << ','
                        << (observations > 0.0
                                ? 1.0e6 * static_cast<double>(
                                      value.field_push_promotions) /
                                      observations
                                : 0.0) << ','
                        << (observations > 0.0
                                ? 1.0e6 * static_cast<double>(
                                      value.field_push_demotions) /
                                      observations
                                : 0.0);
                    for (std::size_t process = 0; process < 6; ++process) {
                        crossings << ',' << value.collision_promotions[process]
                                  << ',' << value.collision_demotions[process];
                    }
                    crossings << ',' << value.energetic_births << ','
                              << value.subthreshold_births << '\n';
                }
            }
            require_stream(
                crossings,
                "failed while writing phase EEDF threshold crossings");
        }
    }

    if (cfg_.phase_surface_flux.enabled) {
        std::ofstream histogram(output_dir / "phase_surface_flux.csv");
        std::ofstream summary(output_dir / "phase_surface_flux_summary.csv");
        if (!histogram || !summary) {
            throw std::runtime_error("cannot open 1D phase surface-flux output");
        }
        const char* energy_name = si ? "energy_eV" : "energy_normalized";
        const char* position_name = si ? "position_m" : "position_normalized";
        const char* energy_flux_name =
            si ? "kinetic_energy_flux_W_m-2" :
                 "kinetic_energy_flux_normalized";
        const char* particle_flux_name =
            si ? "represented_particle_flux_m-2_s-1" :
                 "represented_particle_flux_normalized";
        histogram << "phase_bin,phase_fraction,surface_id," << position_name
                  << ",direction,energy_bin," << energy_name <<
                     ",represented_crossings,probability_density\n"
                  << std::setprecision(17);
        summary << "phase_bin,phase_fraction,surface_id," << position_name <<
                   ",direction,"
                   "macro_crossings,overflow_macro_crossings,"
                   "represented_crossings,overflow_fraction," <<
                   particle_flux_name << ',' << energy_flux_name << "\n"
                << std::setprecision(17);
        const double bin_width = cfg_.phase_surface_flux.energy_max /
            static_cast<double>(cfg_.phase_surface_flux.energy_bins);
        const double phase_duration =
            static_cast<double>(cfg_.spatial_average.rf_cycles) /
            cfg_.spatial_average.rf_frequency /
            static_cast<double>(phase_surface_flux_accumulators_.size());
        for (std::size_t phase = 0;
             phase < phase_surface_flux_accumulators_.size(); ++phase) {
            const double phase_fraction =
                (static_cast<double>(phase) + 0.5) /
                static_cast<double>(phase_surface_flux_accumulators_.size());
            for (std::size_t surface = 0;
                 surface < cfg_.phase_surface_flux.positions.size(); ++surface) {
                for (std::size_t direction = 0; direction < 2; ++direction) {
                    const auto& value =
                        phase_surface_flux_accumulators_[phase][surface][direction];
                    const char* direction_name =
                        direction == 0 ? "left_to_right" : "right_to_left";
                    for (std::size_t bin = 0;
                         bin < value.represented_histogram.size(); ++bin) {
                        histogram << phase << ',' << phase_fraction << ','
                            << surface << ','
                            << cfg_.phase_surface_flux.positions[surface] << ','
                            << direction_name << ',' << bin << ','
                            << (static_cast<double>(bin) + 0.5) * bin_width << ','
                            << value.represented_histogram[bin] << ','
                            << (value.represented_crossings > 0.0
                                    ? value.represented_histogram[bin] /
                                          value.represented_crossings /
                                          bin_width
                                    : 0.0)
                            << '\n';
                    }
                    summary << phase << ',' << phase_fraction << ',' << surface
                        << ',' << cfg_.phase_surface_flux.positions[surface]
                        << ',' << direction_name << ',' << value.macro_crossings
                        << ',' << value.overflow_macro_crossings << ','
                        << value.represented_crossings << ','
                        << (value.represented_crossings > 0.0
                                ? value.overflow_represented_crossings /
                                      value.represented_crossings
                                : 0.0)
                        << ',' << value.represented_crossings / phase_duration
                        << ',' << value.represented_kinetic_energy /
                                      phase_duration
                        << '\n';
                }
            }
        }
        require_stream(histogram,
            "failed while writing phase surface-flux histogram");
        require_stream(summary,
            "failed while writing phase surface-flux summary");
    }

    const auto expected =
        expected_spatial_average_samples();
    const bool complete =
        step_ >= cfg_.spatial_average.end_step &&
        spatial_average_samples_ == expected;
    std::ofstream metadata(
        output_dir / "spatial_average_metadata.json");
    if (!metadata) {
        throw std::runtime_error(
            "cannot open 1D spatial-average metadata output");
    }
    metadata << std::setprecision(17)
             << "{\n"
             << "  \"spatial_average_version\": 7,\n"
             << "  \"reset_on_restart\": "
             << (cfg_.spatial_average.reset_on_restart
                     ? "true" : "false")
             << ",\n"
             << "  \"unit_system\": "
             << json_string(to_string(cfg_.units.system))
             << ",\n"
             << "  \"start_step\": "
             << cfg_.spatial_average.start_step << ",\n"
             << "  \"end_step\": "
             << cfg_.spatial_average.end_step << ",\n"
             << "  \"interval\": "
             << cfg_.spatial_average.interval << ",\n"
             << "  \"sampling_order\": "
             << json_string(to_string(
                    cfg_.spatial_average.sampling_order)) << ",\n"
             << "  \"collision_velocity_sampling\": "
             << json_string(to_string(
                    cfg_.collision_velocity_sampling)) << ",\n"
             << "  \"samples\": "
             << spatial_average_samples_ << ",\n"
             << "  \"moment_samples\": "
             << spatial_moment_samples_ << ",\n"
             << "  \"moments_complete\": "
             << (moments_complete ? "true" : "false") << ",\n"
             << "  \"expected_samples\": "
             << expected << ",\n"
             << "  \"final_step\": " << step_ << ",\n"
             << "  \"dt\": " << cfg_.dt << ",\n"
             << "  \"rf_frequency\": "
             << cfg_.spatial_average.rf_frequency << ",\n"
             << "  \"rf_cycles\": "
             << cfg_.spatial_average.rf_cycles << ",\n"
             << "  \"phase_bins\": "
             << cfg_.spatial_average.phase_bins << ",\n"
             << "  \"phase_bin_samples\": [";
    for (std::size_t phase = 0; phase < spatial_phase_bins_.size(); ++phase) {
        if (phase != 0) metadata << ", ";
        metadata << spatial_phase_bins_[phase].samples;
    }
    metadata << "],\n"
             << "  \"collision_timesteps\": "
             << spatial_collision_steps_ << ",\n"
             << "  \"collision_phase_timesteps\": [";
    for (std::size_t phase = 0;
         phase < spatial_collision_phase_steps_.size(); ++phase) {
        if (phase != 0) metadata << ", ";
        metadata << spatial_collision_phase_steps_[phase];
    }
    metadata << "],\n"
             << "  \"phase_eedf_enabled\": "
             << (cfg_.phase_eedf.enabled ? "true" : "false") << ",\n"
             << "  \"phase_eedf_history_enabled\": "
             << (cfg_.phase_eedf.history_enabled ? "true" : "false")
             << ",\n"
             << "  \"phase_eedf_threshold_crossing_enabled\": "
             << (cfg_.phase_eedf.history_enabled ? "true" : "false")
             << ",\n"
             << "  \"phase_eedf_species\": "
             << json_string(cfg_.phase_eedf.species) << ",\n"
             << "  \"phase_eedf_energy_bins\": "
             << cfg_.phase_eedf.energy_bins << ",\n"
             << "  \"phase_eedf_energy_max\": "
             << cfg_.phase_eedf.energy_max << ",\n"
             << "  \"phase_eedf_tail_threshold\": "
             << cfg_.phase_eedf.tail_threshold << ",\n"
             << "  \"phase_surface_flux_enabled\": "
             << (cfg_.phase_surface_flux.enabled ? "true" : "false") << ",\n"
             << "  \"phase_surface_flux_reset_on_restart\": "
             << (cfg_.phase_surface_flux.reset_on_restart ? "true" : "false")
             << ",\n"
             << "  \"phase_surface_flux_species\": "
             << json_string(cfg_.phase_surface_flux.species) << ",\n"
             << "  \"phase_surface_flux_energy_bins\": "
             << cfg_.phase_surface_flux.energy_bins << ",\n"
             << "  \"phase_surface_flux_energy_max\": "
             << cfg_.phase_surface_flux.energy_max << ",\n"
             << "  \"phase_surface_flux_positions\": [";
    for (std::size_t surface = 0;
         surface < cfg_.phase_surface_flux.positions.size(); ++surface) {
        if (surface != 0) metadata << ", ";
        metadata << cfg_.phase_surface_flux.positions[surface];
    }
    metadata << "],\n"
             << "  \"complete\": "
             << (complete ? "true" : "false") << ",\n"
             << "  \"effective_kinetic_temperature_definition\": "
             << json_string("2 * density-weighted mean total kinetic energy / velocity_dimensions; includes directed energy")
             << ",\n"
             << "  \"phase_temperature_definition\": "
             << json_string("2 / velocity_dimensions * (mean kinetic energy - mass * squared mean velocity / 2); clipped at zero")
             << ",\n"
             << "  \"field_statistics\": "
             << json_string("sampled nodal potential mean, electric-field mean, and electric-field RMS; no sheath edge is inferred")
             << ",\n"
             << "  \"species\": [";
    for (std::size_t species_id = 0;
         species_id < species_.size(); ++species_id) {
        if (species_id != 0) metadata << ", ";
        metadata << json_string(species_[species_id].name());
    }
    metadata << "]\n}\n";
    require_stream(
        metadata,
        "failed while writing 1D spatial-average metadata");
}

void Simulation::save_checkpoint(const std::filesystem::path& path) const {
    ensure_parent_directory(path);
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot open checkpoint for writing: " + path.string());
    out << std::setprecision(17);
    out << kCheckpointMagicV21 << '\n';
    out << "dimension 1\n";
    out << "units " << to_string(cfg_.units.system) << ' '
        << cfg_.units.relative_permittivity << ' '
        << cfg_.units.permittivity() << "\n";
    out << "velocity_dimensions "
        << cfg_.velocity_dimensions << "\n";
    out << "species_timestep_multipliers " << species_.size();
    for (const auto& species : species_) {
        out << ' ' << species.config().timestep_multiplier;
    }
    out << "\n";
    const bool collisions_enabled =
        legacy_bgk_enabled_ || !mcc_models_.empty();
    const std::uint64_t configured_collision_signature =
        collision_signature();
    out << "collision_model "
        << collision_identity() << ' '
        << (collisions_enabled ? 1 : 0) << ' '
        << configured_collision_signature << "\n";
    out << "collision_totals " << collision_totals_.candidates
        << ' ' << collision_totals_.null_collisions << ' '
        << collision_totals_.channel_collisions.size();
    for (const auto count : collision_totals_.channel_collisions) {
        out << ' ' << count;
    }
    out << "\n";
    out << "collision_energy_totals "
        << collision_totals_.channel_energy_change.size();
    for (const double energy : collision_totals_.channel_energy_change) {
        out << ' ' << energy;
    }
    out << "\n";
    out << "boundary_loss_origin_step "
        << boundary_loss_origin_step_ << "\n";
    out << "boundary_loss_count "
        << species_boundary_losses_.size() << "\n";
    for (std::size_t species_id = 0;
         species_id < species_boundary_losses_.size();
         ++species_id) {
        const auto& loss =
            species_boundary_losses_[species_id];
        out << "boundary_loss " << species_id << ' '
            << species_[species_id].name() << ' '
            << loss.absorbed_left << ' '
            << loss.absorbed_right << ' '
            << loss.kinetic_energy_left << ' '
            << loss.kinetic_energy_right << "\n";
    }
    out << "wall_impact_spectrum "
        << (cfg_.wall_impact_spectrum.enabled ? 1 : 0) << ' '
        << wall_impact_origin_step_ << ' '
        << cfg_.wall_impact_spectrum.energy_bins << ' '
        << cfg_.wall_impact_spectrum.energy_max << ' '
        << wall_impact_spectra_.size() << "\n";
    for (std::size_t species_id = 0;
         species_id < wall_impact_spectra_.size(); ++species_id) {
        const auto& spectrum = wall_impact_spectra_[species_id];
        out << "wall_impact_species " << species_id << ' '
            << species_[species_id].name() << ' '
            << spectrum.baseline_loss.absorbed_left << ' '
            << spectrum.baseline_loss.absorbed_right << ' '
            << spectrum.baseline_loss.kinetic_energy_left << ' '
            << spectrum.baseline_loss.kinetic_energy_right << "\n";
        for (const bool left : {true, false}) {
            const auto& side = left
                ? spectrum.left : spectrum.right;
            out << "wall_impact_side " << species_id << ' '
                << (left ? "left" : "right") << ' '
                << side.macro_impacts << ' '
                << side.overflow_macro_impacts << ' '
                << side.represented_impacts << ' '
                << side.overflow_represented_impacts << ' '
                << side.represented_kinetic_energy;
            for (std::size_t bin = 0;
                 bin < side.macro_histogram.size(); ++bin) {
                out << ' ' << side.macro_histogram[bin]
                    << ' ' << side.represented_histogram[bin];
            }
            out << "\n";
        }
    }
    out << "power_transfer_origin_step "
        << power_transfer_origin_step_ << "\n";
    out << "power_transfer_count "
        << species_power_transfer_.size() << "\n";
    for (std::size_t species_id = 0;
         species_id < species_power_transfer_.size();
         ++species_id) {
        out << "power_transfer " << species_id << ' '
            << species_[species_id].name() << ' '
            << species_power_transfer_[species_id]
                   .electric_work
            << "\n";
    }
    out << "spatial_average "
        << (cfg_.spatial_average.enabled ? 1 : 0) << ' '
        << to_string(cfg_.spatial_average.sampling_order) << ' '
        << cfg_.spatial_average.interval << ' '
        << cfg_.spatial_average.start_step << ' '
        << cfg_.spatial_average.end_step << ' '
        << cfg_.spatial_average.rf_frequency << ' '
        << cfg_.spatial_average.rf_cycles << ' '
        << spatial_average_samples_ << ' '
        << spatial_density_sums_.size() << ' '
        << grid_.nx() << "\n";
    for (std::size_t species_id = 0;
         species_id < spatial_density_sums_.size();
         ++species_id) {
        out << "spatial_species " << species_id << ' '
            << species_[species_id].name();
        for (const double value :
             spatial_density_sums_[species_id]) {
            out << ' ' << value;
        }
        out << "\n";
    }
    out << "spatial_moments " << spatial_moment_samples_ << ' '
        << spatial_kinetic_energy_sums_.size() << ' '
        << grid_.nx() << "\n";
    for (std::size_t species_id = 0;
         species_id < spatial_kinetic_energy_sums_.size(); ++species_id) {
        out << "spatial_energy " << species_id << ' '
            << species_[species_id].name();
        for (const double value : spatial_kinetic_energy_sums_[species_id]) {
            out << ' ' << value;
        }
        out << "\n";
    }
    out << "spatial_fields";
    for (std::size_t node = 0; node < grid_.nx(); ++node) {
        const bool stored = node < spatial_potential_sums_.size();
        out << ' ' << (stored ? spatial_potential_sums_[node] : 0.0)
            << ' ' << (stored ? spatial_electric_sums_[node] : 0.0)
            << ' ' << (stored ? spatial_electric_squared_sums_[node] : 0.0);
    }
    out << "\n";
    out << "spatial_phase " << spatial_phase_bins_.size() << ' '
        << species_.size() << ' ' << grid_.nx() << "\n";
    for (std::size_t phase = 0; phase < spatial_phase_bins_.size(); ++phase) {
        const auto& bin = spatial_phase_bins_[phase];
        out << "phase_bin " << phase << ' ' << bin.samples << "\n";
        for (std::size_t species_id = 0; species_id < species_.size(); ++species_id) {
            out << "phase_species " << phase << ' ' << species_id << ' '
                << species_[species_id].name();
            for (std::size_t node = 0; node < grid_.nx(); ++node) {
                out << ' ' << bin.density[species_id][node]
                    << ' ' << bin.velocity_x_density[species_id][node]
                    << ' ' << bin.velocity_y_density[species_id][node]
                    << ' ' << bin.velocity_z_density[species_id][node]
                    << ' ' << bin.kinetic_energy_density[species_id][node];
            }
            out << "\n";
        }
        out << "phase_fields " << phase;
        for (std::size_t node = 0; node < grid_.nx(); ++node) {
            out << ' ' << bin.potential[node] << ' ' << bin.electric[node]
                << ' ' << bin.electric_squared[node];
        }
        out << "\n";
    }
    out << "spatial_collision " << spatial_collision_steps_ << ' '
        << spatial_collision_energy_sums_.size() << ' '
        << grid_.nx() << ' '
        << spatial_collision_phase_energy_sums_.size() << "\n";
    for (std::size_t channel = 0;
         channel < spatial_collision_energy_sums_.size(); ++channel) {
        out << "spatial_collision_channel " << channel << ' '
            << collision_totals_.channel_names[channel];
        for (const double value : spatial_collision_energy_sums_[channel]) {
            out << ' ' << value;
        }
        out << "\n";
    }
    out << "spatial_collision_phase_steps";
    for (const auto count : spatial_collision_phase_steps_) {
        out << ' ' << count;
    }
    out << "\n";
    for (std::size_t phase = 0;
         phase < spatial_collision_phase_energy_sums_.size(); ++phase) {
        for (std::size_t channel = 0;
             channel < spatial_collision_energy_sums_.size(); ++channel) {
            out << "spatial_collision_phase_channel " << phase << ' '
                << channel;
            for (const double value :
                 spatial_collision_phase_energy_sums_[phase][channel]) {
                out << ' ' << value;
            }
            out << "\n";
        }
    }
    out << "spatial_collision_events "
        << spatial_collision_event_sums_.size() << ' '
        << grid_.nx() << ' '
        << spatial_collision_phase_event_sums_.size() << "\n";
    for (std::size_t channel = 0;
         channel < spatial_collision_event_sums_.size(); ++channel) {
        out << "spatial_collision_event_channel " << channel << ' '
            << collision_totals_.channel_names[channel];
        for (const double value : spatial_collision_event_sums_[channel]) {
            out << ' ' << value;
        }
        out << "\n";
    }
    for (std::size_t phase = 0;
         phase < spatial_collision_phase_event_sums_.size(); ++phase) {
        for (std::size_t channel = 0;
             channel < spatial_collision_event_sums_.size(); ++channel) {
            out << "spatial_collision_phase_event_channel " << phase << ' '
                << channel;
            for (const double value :
                 spatial_collision_phase_event_sums_[phase][channel]) {
                out << ' ' << value;
            }
            out << "\n";
        }
    }
    out << "phase_eedf " << (cfg_.phase_eedf.enabled ? 1 : 0) << ' '
        << (cfg_.phase_eedf.history_enabled ? 1 : 0) << ' '
        << (cfg_.phase_eedf.species.empty() ? "-" : cfg_.phase_eedf.species)
        << ' ' << cfg_.phase_eedf.energy_bins << ' '
        << cfg_.phase_eedf.energy_max << ' '
        << cfg_.phase_eedf.tail_threshold << ' '
        << phase_eedf_accumulators_.size() << ' '
        << cfg_.phase_eedf.regions.size() << "\n";
    for (std::size_t region = 0;
         region < cfg_.phase_eedf.regions.size(); ++region) {
        const auto& configured = cfg_.phase_eedf.regions[region];
        out << "phase_eedf_region " << region << ' ' << configured.name << ' '
            << configured.x_min << ' ' << configured.x_max << "\n";
    }
    for (std::size_t phase = 0;
         phase < phase_eedf_accumulators_.size(); ++phase) {
        for (std::size_t region = 0;
             region < phase_eedf_accumulators_[phase].size(); ++region) {
            const auto& value = phase_eedf_accumulators_[phase][region];
            out << "phase_eedf_accumulator " << phase << ' ' << region << ' '
                << value.macro_observations << ' '
                << value.overflow_macro_observations << ' '
                << value.represented_observations << ' '
                << value.overflow_represented_observations << ' '
                << value.weighted_energy_sum << ' '
                << value.weighted_energy_squared_sum << ' '
                << value.weighted_velocity_x_sum << ' '
                << value.weighted_velocity_y_sum << ' '
                << value.weighted_velocity_z_sum << ' '
                << value.weighted_velocity_x_squared_sum << ' '
                << value.weighted_velocity_y_squared_sum << ' '
                << value.weighted_velocity_z_squared_sum << ' '
                << value.tail_represented_observations << ' '
                << value.tail_positive_x_represented_observations << ' '
                << value.tail_negative_x_represented_observations << ' '
                << value.tail_weighted_velocity_x_sum << ' '
                << value.tail_weighted_velocity_x_squared_sum << ' '
                << value.tail_weighted_transverse_velocity_squared_sum << ' '
                << value.tail_weighted_age_steps_sum << ' '
                << value.tail_weighted_energetic_steps_sum << ' '
                << value.tail_weighted_energetic_duty_fraction_sum << ' '
                << value.tail_weighted_consecutive_energetic_steps_sum << ' '
                << value.tail_weighted_entries_sum << ' '
                << value.tail_weighted_elastic_collisions_sum << ' '
                << value.tail_weighted_excitation_collisions_sum << ' '
                << value.tail_weighted_ionization_collisions_sum << ' '
                << value.tail_weighted_charge_exchange_collisions_sum << ' '
                << value.tail_weighted_bgk_collisions_sum << ' '
                << value.tail_born_during_window_represented_observations;
            for (const double count : value.histogram) out << ' ' << count;
            out << "\n";
        }
    }
    out << "phase_eedf_threshold_crossings "
        << (cfg_.phase_eedf.history_enabled ? 1 : 0) << ' '
        << phase_eedf_threshold_crossings_.size() << ' '
        << (phase_eedf_threshold_crossings_.empty()
                ? 0 : phase_eedf_threshold_crossings_.front().size())
        << " 6\n";
    for (std::size_t phase = 0;
         phase < phase_eedf_threshold_crossings_.size(); ++phase) {
        for (std::size_t region = 0;
             region < phase_eedf_threshold_crossings_[phase].size();
             ++region) {
            const auto& value =
                phase_eedf_threshold_crossings_[phase][region];
            out << "phase_eedf_threshold_crossing_accumulator "
                << phase << ' ' << region << ' '
                << value.electron_time_macro_observations << ' '
                << value.energetic_time_macro_observations << ' '
                << value.interstep_promotions << ' '
                << value.interstep_demotions << ' '
                << value.field_push_promotions << ' '
                << value.field_push_demotions;
            for (const auto count : value.collision_promotions) {
                out << ' ' << count;
            }
            for (const auto count : value.collision_demotions) {
                out << ' ' << count;
            }
            out << ' ' << value.energetic_births << ' '
                << value.subthreshold_births << '\n';
        }
    }
    out << "phase_surface_flux "
        << (cfg_.phase_surface_flux.enabled ? 1 : 0) << ' '
        << (cfg_.phase_surface_flux.species.empty()
                ? "-" : cfg_.phase_surface_flux.species) << ' '
        << cfg_.phase_surface_flux.energy_bins << ' '
        << cfg_.phase_surface_flux.energy_max << ' '
        << phase_surface_flux_accumulators_.size() << ' '
        << cfg_.phase_surface_flux.positions.size() << "\n";
    for (std::size_t surface = 0;
         surface < cfg_.phase_surface_flux.positions.size(); ++surface) {
        out << "phase_surface_flux_position " << surface << ' '
            << cfg_.phase_surface_flux.positions[surface] << "\n";
    }
    for (std::size_t phase = 0;
         phase < phase_surface_flux_accumulators_.size(); ++phase) {
        for (std::size_t surface = 0;
             surface < phase_surface_flux_accumulators_[phase].size();
             ++surface) {
            for (std::size_t direction = 0; direction < 2; ++direction) {
                const auto& value =
                    phase_surface_flux_accumulators_[phase][surface][direction];
                out << "phase_surface_flux_accumulator " << phase << ' '
                    << surface << ' ' << direction << ' '
                    << value.macro_crossings << ' '
                    << value.overflow_macro_crossings << ' '
                    << value.represented_crossings << ' '
                    << value.overflow_represented_crossings << ' '
                    << value.represented_kinetic_energy;
                for (const double count : value.represented_histogram) {
                    out << ' ' << count;
                }
                out << "\n";
            }
        }
    }
    out << "step " << step_ << "\n";
    out << "time " << time_ << "\n";
    out << "species_count " << species_.size() << "\n";
    out << "rng " << rng_ << "\n";
    for (std::size_t species_id = 0; species_id < species_.size(); ++species_id) {
        const auto& sp = species_[species_id];
        out << "species " << species_id << ' ' << sp.name() << ' ' << sp.particles().size() << "\n";
        for (const auto& p : sp.particles()) {
            out << p.x << ' ' << p.v << ' '
                << p.velocity_y << ' ' << p.velocity_z << ' '
                << p.v_half << ' '
                << (p.alive ? 1 : 0) << "\n";
        }
    }
    out << "phase_eedf_particle_history "
        << (cfg_.phase_eedf.history_enabled ? 1 : 0) << ' '
        << phase_eedf_species_id_ << ' '
        << phase_eedf_particle_histories_.size() << "\n";
    for (std::size_t particle_id = 0;
         particle_id < phase_eedf_particle_histories_.size(); ++particle_id) {
        const auto& history = phase_eedf_particle_histories_[particle_id];
        out << "phase_eedf_particle_history_entry " << particle_id << ' '
            << history.age_steps << ' ' << history.energetic_steps << ' '
            << history.consecutive_energetic_steps << ' '
            << history.tail_entries << ' ' << history.elastic_collisions << ' '
            << history.excitation_collisions << ' '
            << history.ionization_collisions << ' '
            << history.charge_exchange_collisions << ' '
            << history.bgk_collisions << ' '
            << (history.born_during_window ? 1 : 0) << ' '
            << (history.energetic_previous_step ? 1 : 0) << "\n";
    }
    require_stream(out, "failed while writing checkpoint: " + path.string());
}

void Simulation::load_checkpoint(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open checkpoint for reading: " + path.string());
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
    const bool checkpoint_v11 = magic == kCheckpointMagicV11;
    const bool checkpoint_v12 = magic == kCheckpointMagicV12;
    const bool checkpoint_v21 = magic == kCheckpointMagicV21;
    const bool checkpoint_v20 =
        magic == kCheckpointMagicV20 || checkpoint_v21;
    const bool checkpoint_v19 =
        magic == kCheckpointMagicV19 || checkpoint_v20;
    const bool checkpoint_v18 =
        magic == kCheckpointMagicV18 || checkpoint_v19;
    const bool checkpoint_v17 =
        magic == kCheckpointMagicV17 || checkpoint_v18;
    const bool checkpoint_v16 =
        magic == kCheckpointMagicV16 || checkpoint_v17;
    const bool checkpoint_v15 =
        magic == kCheckpointMagicV15 || checkpoint_v16;
    const bool checkpoint_v14 =
        magic == kCheckpointMagicV14 || checkpoint_v15;
    const bool checkpoint_v13 =
        magic == kCheckpointMagicV13 || checkpoint_v14;
    const bool checkpoint_v4_state =
        checkpoint_v4 || checkpoint_v5 ||
        checkpoint_v6 || checkpoint_v7 || checkpoint_v8 || checkpoint_v9 ||
        checkpoint_v10 || checkpoint_v11 || checkpoint_v12 ||
        checkpoint_v13;
    if (!checkpoint_v1 && !checkpoint_v2 &&
        !checkpoint_v3 && !checkpoint_v4 &&
        !checkpoint_v5 && !checkpoint_v6 &&
        !checkpoint_v7 && !checkpoint_v8 && !checkpoint_v9 &&
        !checkpoint_v10 && !checkpoint_v11 && !checkpoint_v12 &&
        !checkpoint_v13) {
        throw std::runtime_error("invalid checkpoint magic in: " + path.string());
    }

    std::string key;
    unsigned dimension = 0;
    in >> key >> dimension;
    if (key != "dimension" || dimension != 1) throw std::runtime_error("checkpoint dimension does not match 1D simulation");
    in >> key;
    if (key == "units") {
        std::string unit_system;
        double relative_permittivity = 0.0;
        double permittivity = 0.0;
        in >> unit_system >> relative_permittivity >> permittivity;
        if ((!checkpoint_v2 && !checkpoint_v3 &&
             !checkpoint_v4 && !checkpoint_v5 &&
             !checkpoint_v6 && !checkpoint_v7 && !checkpoint_v8 &&
             !checkpoint_v9 && !checkpoint_v10 && !checkpoint_v11 &&
             !checkpoint_v12 && !checkpoint_v13) ||
            unit_system != to_string(cfg_.units.system) ||
            relative_permittivity != cfg_.units.relative_permittivity ||
            permittivity != cfg_.units.permittivity()) {
            throw std::runtime_error(
                "checkpoint unit system does not match 1D config");
        }
        in >> key;
    } else if (checkpoint_v2 || checkpoint_v3 ||
               checkpoint_v4 || checkpoint_v5 ||
               checkpoint_v6 || checkpoint_v7 || checkpoint_v8 ||
               checkpoint_v9 || checkpoint_v10 || checkpoint_v11 ||
               checkpoint_v12 || checkpoint_v13 ||
               cfg_.units.system != UnitSystem::Normalized ||
               cfg_.units.relative_permittivity != 1.0) {
        throw std::runtime_error(
            "legacy checkpoint without unit metadata requires normalized units");
    }
    if (key == "velocity_dimensions") {
        std::size_t velocity_dimensions = 0;
        in >> velocity_dimensions;
        if (!checkpoint_v4_state ||
            velocity_dimensions != cfg_.velocity_dimensions) {
            throw std::runtime_error(
                "checkpoint velocity dimensions do not match 1D config");
        }
        in >> key;
    } else if (checkpoint_v4_state ||
               cfg_.velocity_dimensions != 1) {
        throw std::runtime_error(
            "legacy 1D1V checkpoint cannot initialize 1D3V");
    }
    if (key == "species_timestep_multipliers") {
        std::size_t stored_species_count = 0;
        in >> stored_species_count;
        if (!checkpoint_v13 ||
            stored_species_count != species_.size()) {
            throw std::runtime_error(
                "checkpoint species timestep count does not "
                "match 1D config");
        }
        for (std::size_t species_id = 0;
             species_id < stored_species_count; ++species_id) {
            std::size_t stored_multiplier = 0;
            in >> stored_multiplier;
            if (stored_multiplier != species_[species_id]
                    .config().timestep_multiplier) {
                throw std::runtime_error(
                    "checkpoint species timestep multiplier "
                    "does not match 1D config");
            }
        }
        in >> key;
    } else if (checkpoint_v13 ||
               std::any_of(
                   species_.begin(), species_.end(),
                   [](const Species& species) {
                       return species.config()
                                  .timestep_multiplier != 1;
                   })) {
        throw std::runtime_error(
            "legacy checkpoint without species timestep metadata "
            "requires timestep_multiplier = 1");
    }
    if (key == "collision_model") {
        std::string model;
        int enabled = 0;
        std::uint64_t signature = 0;
        in >> model >> enabled >> signature;
        const bool collisions_enabled =
            legacy_bgk_enabled_ || !mcc_models_.empty();
        const std::uint64_t expected_signature =
            collision_signature();
        if ((!checkpoint_v3 && !checkpoint_v4 &&
             !checkpoint_v5 && !checkpoint_v6 &&
             !checkpoint_v7 && !checkpoint_v8 && !checkpoint_v9 &&
             !checkpoint_v10 && !checkpoint_v11 && !checkpoint_v12 &&
             !checkpoint_v13) ||
            model != collision_identity() ||
            enabled != (collisions_enabled ? 1 : 0) ||
            signature != expected_signature) {
            throw std::runtime_error(
                "checkpoint collision model does not match 1D config");
        }
        std::size_t channel_count = 0;
        in >> key >> collision_totals_.candidates
           >> collision_totals_.null_collisions >> channel_count;
        if (key != "collision_totals" ||
            channel_count !=
                collision_totals_.channel_collisions.size()) {
            throw std::runtime_error(
                "checkpoint collision diagnostics do not match 1D config");
        }
        for (auto& count : collision_totals_.channel_collisions) {
            in >> count;
        }
        clear_collision_counts(collision_interval_);
        in >> key;
        if (checkpoint_v10 || checkpoint_v11 || checkpoint_v12 ||
            checkpoint_v13) {
            std::size_t energy_count = 0;
            in >> energy_count;
            if (key != "collision_energy_totals" ||
                energy_count != collision_totals_.channel_energy_change.size()) {
                throw std::runtime_error(
                    "checkpoint collision energy diagnostics do not match config");
            }
            for (auto& energy : collision_totals_.channel_energy_change) {
                in >> energy;
                if (!std::isfinite(energy)) {
                    throw std::runtime_error(
                        "checkpoint collision energy diagnostic is invalid");
                }
            }
            in >> key;
        } else {
            std::fill(collision_totals_.channel_energy_change.begin(),
                      collision_totals_.channel_energy_change.end(), 0.0);
        }
    } else if (checkpoint_v3 || checkpoint_v4 ||
               checkpoint_v5 || checkpoint_v6 ||
               checkpoint_v7 || checkpoint_v8 || checkpoint_v9 ||
               checkpoint_v10 || checkpoint_v11 || checkpoint_v12 ||
               checkpoint_v13 ||
               !mcc_models_.empty()) {
        throw std::runtime_error(
            "legacy checkpoint without MCC metadata cannot restart "
            "null-collision MCC");
    }
    const bool checkpoint_has_boundary_losses =
        checkpoint_v6 || checkpoint_v7 || checkpoint_v8 || checkpoint_v9 ||
        checkpoint_v10 || checkpoint_v11 || checkpoint_v12 ||
        checkpoint_v13;
    const bool legacy_boundary_loss_origin =
        !checkpoint_has_boundary_losses;
    if (checkpoint_has_boundary_losses) {
        in >> boundary_loss_origin_step_;
        if (key != "boundary_loss_origin_step") {
            throw std::runtime_error(
                "checkpoint is missing the 1D boundary-loss "
                "counter origin");
        }
        in >> key;
        std::size_t stored_loss_count = 0;
        in >> stored_loss_count;
        if (key != "boundary_loss_count" ||
            stored_loss_count !=
                species_boundary_losses_.size()) {
            throw std::runtime_error(
                "checkpoint boundary-loss species count does "
                "not match 1D config");
        }
        for (std::size_t species_id = 0;
             species_id < stored_loss_count; ++species_id) {
            std::size_t stored_species_id = 0;
            std::string stored_name;
            auto& loss =
                species_boundary_losses_[species_id];
            in >> key >> stored_species_id >> stored_name >>
                loss.absorbed_left >>
                loss.absorbed_right >>
                loss.kinetic_energy_left >>
                loss.kinetic_energy_right;
            if (key != "boundary_loss" ||
                stored_species_id != species_id ||
                stored_name != species_[species_id].name() ||
                !std::isfinite(loss.kinetic_energy_left) ||
                loss.kinetic_energy_left < 0.0 ||
                !std::isfinite(loss.kinetic_energy_right) ||
                loss.kinetic_energy_right < 0.0) {
                throw std::runtime_error(
                    "checkpoint 1D boundary-loss data are "
                    "invalid");
            }
        }
        in >> key;
    } else {
        std::fill(
            species_boundary_losses_.begin(),
            species_boundary_losses_.end(),
            BoundaryLoss1D{});
    }
    const bool legacy_wall_impact_origin =
        !checkpoint_v14 && cfg_.wall_impact_spectrum.enabled;
    if (checkpoint_v14) {
        int stored_enabled = 0;
        std::size_t stored_origin = 0;
        std::size_t stored_bins = 0;
        double stored_energy_max = 0.0;
        std::size_t stored_species_count = 0;
        in >> stored_enabled >> stored_origin >> stored_bins >>
            stored_energy_max >> stored_species_count;
        const bool enabled = stored_enabled == 1;
        const bool shape_valid =
            (stored_enabled == 0 || stored_enabled == 1) &&
            std::isfinite(stored_energy_max) &&
            stored_energy_max >= 0.0 &&
            stored_species_count ==
                (enabled ? species_.size() : 0) &&
            (!enabled ||
             (stored_bins > 0 && stored_energy_max > 0.0));
        const bool contract_matches =
            stored_enabled ==
                (cfg_.wall_impact_spectrum.enabled ? 1 : 0) &&
            stored_bins ==
                cfg_.wall_impact_spectrum.energy_bins &&
            stored_energy_max ==
                cfg_.wall_impact_spectrum.energy_max &&
            stored_species_count == wall_impact_spectra_.size();
        const bool enable_with_fresh_window =
            !enabled && cfg_.wall_impact_spectrum.enabled &&
            cfg_.wall_impact_spectrum.reset_on_restart &&
            stored_bins == 0 && stored_energy_max == 0.0 &&
            stored_species_count == 0;
        if (key != "wall_impact_spectrum" ||
            !shape_valid ||
            (!contract_matches && !enable_with_fresh_window)) {
            throw std::runtime_error(
                "checkpoint wall-impact spectrum contract is "
                "invalid");
        }
        wall_impact_origin_step_ = stored_origin;
        for (std::size_t species_id = 0;
             species_id < stored_species_count; ++species_id) {
            std::size_t stored_species_id = 0;
            std::string stored_name;
            auto& spectrum = wall_impact_spectra_[species_id];
            in >> key >> stored_species_id >> stored_name >>
                spectrum.baseline_loss.absorbed_left >>
                spectrum.baseline_loss.absorbed_right >>
                spectrum.baseline_loss.kinetic_energy_left >>
                spectrum.baseline_loss.kinetic_energy_right;
            if (key != "wall_impact_species" ||
                stored_species_id != species_id ||
                stored_name != species_[species_id].name() ||
                spectrum.baseline_loss.absorbed_left >
                    species_boundary_losses_[species_id]
                        .absorbed_left ||
                spectrum.baseline_loss.absorbed_right >
                    species_boundary_losses_[species_id]
                        .absorbed_right ||
                !std::isfinite(
                    spectrum.baseline_loss.kinetic_energy_left) ||
                !std::isfinite(
                    spectrum.baseline_loss.kinetic_energy_right) ||
                spectrum.baseline_loss.kinetic_energy_left < 0.0 ||
                spectrum.baseline_loss.kinetic_energy_right < 0.0 ||
                spectrum.baseline_loss.kinetic_energy_left >
                    species_boundary_losses_[species_id]
                        .kinetic_energy_left ||
                spectrum.baseline_loss.kinetic_energy_right >
                    species_boundary_losses_[species_id]
                        .kinetic_energy_right) {
                throw std::runtime_error(
                    "checkpoint wall-impact species baseline is "
                    "invalid");
            }
            for (const bool left : {true, false}) {
                std::size_t side_species_id = 0;
                std::string side_name;
                auto& side = left
                    ? spectrum.left : spectrum.right;
                in >> key >> side_species_id >> side_name >>
                    side.macro_impacts >>
                    side.overflow_macro_impacts >>
                    side.represented_impacts >>
                    side.overflow_represented_impacts >>
                    side.represented_kinetic_energy;
                std::uint64_t binned_macro = 0;
                double binned_represented = 0.0;
                for (std::size_t bin = 0;
                     bin < stored_bins; ++bin) {
                    in >> side.macro_histogram[bin] >>
                        side.represented_histogram[bin];
                    binned_macro += side.macro_histogram[bin];
                    binned_represented +=
                        side.represented_histogram[bin];
                }
                const auto boundary_count = left
                    ? species_boundary_losses_[species_id]
                              .absorbed_left -
                          spectrum.baseline_loss.absorbed_left
                    : species_boundary_losses_[species_id]
                              .absorbed_right -
                          spectrum.baseline_loss.absorbed_right;
                const double boundary_energy = left
                    ? species_boundary_losses_[species_id]
                              .kinetic_energy_left -
                          spectrum.baseline_loss.kinetic_energy_left
                    : species_boundary_losses_[species_id]
                              .kinetic_energy_right -
                          spectrum.baseline_loss.kinetic_energy_right;
                const double energy_tolerance =
                    256.0 * std::numeric_limits<double>::epsilon() *
                    std::max({std::numeric_limits<double>::min(),
                        std::abs(boundary_energy),
                        std::abs(side.represented_kinetic_energy)});
                const double count_tolerance =
                    256.0 * std::numeric_limits<double>::epsilon() *
                    std::max(1.0, side.represented_impacts);
                if (key != "wall_impact_side" ||
                    side_species_id != species_id ||
                    side_name != (left ? "left" : "right") ||
                    side.overflow_macro_impacts >
                        side.macro_impacts ||
                    binned_macro +
                            side.overflow_macro_impacts !=
                        side.macro_impacts ||
                    side.macro_impacts != boundary_count ||
                    !std::isfinite(side.represented_impacts) ||
                    !std::isfinite(
                        side.overflow_represented_impacts) ||
                    !std::isfinite(
                        side.represented_kinetic_energy) ||
                    side.represented_impacts < 0.0 ||
                    side.overflow_represented_impacts < 0.0 ||
                    side.overflow_represented_impacts >
                        side.represented_impacts ||
                    std::abs(
                        binned_represented +
                            side.overflow_represented_impacts -
                            side.represented_impacts) >
                        count_tolerance ||
                    std::abs(
                        side.represented_kinetic_energy -
                            boundary_energy) >
                        energy_tolerance ||
                    std::any_of(
                        side.represented_histogram.begin(),
                        side.represented_histogram.end(),
                        [](double count) {
                            return !std::isfinite(count) ||
                                   count < 0.0;
                        })) {
                    throw std::runtime_error(
                        "checkpoint wall-impact side data are "
                        "invalid");
                }
            }
        }
        in >> key;
    }
    const bool legacy_power_transfer_origin =
        !checkpoint_v7 && !checkpoint_v8 && !checkpoint_v9 &&
        !checkpoint_v10 && !checkpoint_v11 && !checkpoint_v12 &&
        !checkpoint_v13;
    if (checkpoint_v7 || checkpoint_v8 || checkpoint_v9 || checkpoint_v10 ||
        checkpoint_v11 || checkpoint_v12 || checkpoint_v13) {
        in >> power_transfer_origin_step_;
        if (key != "power_transfer_origin_step") {
            throw std::runtime_error(
                "checkpoint is missing the 1D power-transfer "
                "counter origin");
        }
        in >> key;
        std::size_t stored_power_count = 0;
        in >> stored_power_count;
        if (key != "power_transfer_count" ||
            stored_power_count !=
                species_power_transfer_.size()) {
            throw std::runtime_error(
                "checkpoint power-transfer species count does "
                "not match 1D config");
        }
        for (std::size_t species_id = 0;
             species_id < stored_power_count; ++species_id) {
            std::size_t stored_species_id = 0;
            std::string stored_name;
            auto& power =
                species_power_transfer_[species_id];
            in >> key >> stored_species_id >> stored_name >>
                power.electric_work;
            if (key != "power_transfer" ||
                stored_species_id != species_id ||
                stored_name != species_[species_id].name() ||
                !std::isfinite(power.electric_work)) {
                throw std::runtime_error(
                    "checkpoint 1D power-transfer data are "
                    "invalid");
            }
        }
        in >> key;
    } else {
        std::fill(
            species_power_transfer_.begin(),
            species_power_transfer_.end(),
            SpeciesPower1D{});
    }
    if (checkpoint_v5 || checkpoint_v6 ||
        checkpoint_v7 || checkpoint_v8 || checkpoint_v9 || checkpoint_v10 ||
        checkpoint_v11 || checkpoint_v12 || checkpoint_v13) {
        int enabled = 0;
        std::string stored_sampling_order_name = "post_collision";
        std::size_t interval = 0;
        std::size_t start_step = 0;
        std::size_t end_step = 0;
        double rf_frequency = 0.0;
        std::size_t rf_cycles = 0;
        std::size_t stored_samples = 0;
        std::size_t stored_species_count = 0;
        std::size_t stored_nx = 0;
        in >> enabled;
        if (checkpoint_v16) in >> stored_sampling_order_name;
        in >> interval >> start_step >>
            end_step >> rf_frequency >> rf_cycles >>
            stored_samples >> stored_species_count >>
            stored_nx;
        const auto& configured = cfg_.spatial_average;
        SpatialAverageSamplingOrder1D stored_sampling_order{};
        if (stored_sampling_order_name == "post_collision") {
            stored_sampling_order =
                SpatialAverageSamplingOrder1D::PostCollision;
        } else if (stored_sampling_order_name == "pre_collision") {
            stored_sampling_order =
                SpatialAverageSamplingOrder1D::PreCollision;
        } else {
            throw std::runtime_error(
                "checkpoint spatial-average sampling order is invalid");
        }
        const bool reset =
            configured.reset_on_restart;
        const bool stored_enabled = enabled == 1;
        const bool stored_shape_valid =
            (enabled == 0 || enabled == 1) &&
            stored_nx == grid_.nx() &&
            stored_species_count ==
                (stored_enabled ? species_.size() : 0);
        std::size_t stored_expected_samples = 0;
        if (stored_enabled && interval > 0 &&
            start_step > 0 && end_step >= start_step) {
            stored_expected_samples =
                1 + (end_step - start_step) / interval;
        }
        const bool stored_state_valid =
            stored_shape_valid &&
            std::isfinite(rf_frequency) &&
            rf_frequency >= 0.0 &&
            (!stored_enabled ||
             (interval > 0 && start_step > 0 &&
              end_step >= start_step &&
              stored_samples <= stored_expected_samples)) &&
            (stored_enabled || stored_samples == 0);
        const bool configured_contract_matches =
            enabled == (configured.enabled ? 1 : 0) &&
            interval == configured.interval &&
            start_step == configured.start_step &&
            end_step == configured.end_step &&
            rf_frequency == configured.rf_frequency &&
            rf_cycles == configured.rf_cycles &&
            stored_sampling_order == configured.sampling_order &&
            stored_species_count ==
                spatial_density_sums_.size() &&
            stored_samples <=
                expected_spatial_average_samples();
        if (key != "spatial_average" ||
            !stored_state_valid ||
            (!reset && !configured_contract_matches)) {
            throw std::runtime_error(
                "checkpoint spatial-average contract does not "
                "match 1D config");
        }
        spatial_average_samples_ =
            reset ? 0 : stored_samples;
        for (std::size_t species_id = 0;
             species_id < stored_species_count;
             ++species_id) {
            std::size_t stored_species_id = 0;
            std::string stored_name;
            in >> key >> stored_species_id >> stored_name;
            if (key != "spatial_species" ||
                stored_species_id != species_id ||
                stored_name != species_[species_id].name()) {
                throw std::runtime_error(
                    "checkpoint spatial-average species metadata "
                    "does not match 1D config");
            }
            for (std::size_t node = 0;
                 node < stored_nx; ++node) {
                double value = 0.0;
                in >> value;
                if (!std::isfinite(value) || value < 0.0) {
                    throw std::runtime_error(
                        "checkpoint spatial-average sum is invalid");
                }
                if (!reset) {
                    spatial_density_sums_[species_id][node] =
                        value;
                }
            }
        }
        if (reset) {
            for (auto& sum : spatial_density_sums_) {
                std::fill(sum.begin(), sum.end(), 0.0);
            }
        }
        in >> key;
        if (checkpoint_v8 || checkpoint_v9 || checkpoint_v10 ||
            checkpoint_v11 || checkpoint_v12 || checkpoint_v13) {
            std::size_t stored_moment_samples = 0;
            std::size_t stored_moment_species_count = 0;
            std::size_t stored_moment_nx = 0;
            in >> stored_moment_samples >>
                stored_moment_species_count >> stored_moment_nx;
            const bool moment_shape_valid =
                key == "spatial_moments" &&
                stored_moment_nx == grid_.nx() &&
                stored_moment_species_count == stored_species_count &&
                stored_moment_samples == stored_samples &&
                (stored_enabled || stored_moment_samples == 0);
            if (!moment_shape_valid) {
                throw std::runtime_error(
                    "checkpoint spatial-moment contract is invalid");
            }
            spatial_moment_samples_ =
                reset ? 0 : stored_moment_samples;
            for (std::size_t species_id = 0;
                 species_id < stored_moment_species_count; ++species_id) {
                std::size_t stored_species_id = 0;
                std::string stored_name;
                in >> key >> stored_species_id >> stored_name;
                if (key != "spatial_energy" ||
                    stored_species_id != species_id ||
                    stored_name != species_[species_id].name()) {
                    throw std::runtime_error(
                        "checkpoint spatial-energy species metadata is invalid");
                }
                for (std::size_t node = 0; node < stored_moment_nx; ++node) {
                    double value = 0.0;
                    in >> value;
                    if (!std::isfinite(value) || value < 0.0) {
                        throw std::runtime_error(
                            "checkpoint spatial-energy sum is invalid");
                    }
                    if (!reset) {
                        spatial_kinetic_energy_sums_[species_id][node] = value;
                    }
                }
            }
            in >> key;
            if (key != "spatial_fields") {
                throw std::runtime_error(
                    "checkpoint spatial-field data are missing");
            }
            for (std::size_t node = 0; node < stored_moment_nx; ++node) {
                double potential = 0.0;
                double electric = 0.0;
                double electric_squared = 0.0;
                in >> potential >> electric >> electric_squared;
                if (!std::isfinite(potential) ||
                    !std::isfinite(electric) ||
                    !std::isfinite(electric_squared) ||
                    electric_squared < 0.0) {
                    throw std::runtime_error(
                        "checkpoint spatial-field sum is invalid");
                }
                if (!reset && stored_enabled) {
                    spatial_potential_sums_[node] = potential;
                    spatial_electric_sums_[node] = electric;
                    spatial_electric_squared_sums_[node] = electric_squared;
                }
            }
            in >> key;
        } else {
            spatial_moment_samples_ = 0;
        }
        std::size_t stored_phase_count = 0;
        if (checkpoint_v9 || checkpoint_v10 || checkpoint_v11 ||
            checkpoint_v12 || checkpoint_v13) {
            std::size_t stored_phase_species = 0;
            std::size_t stored_phase_nx = 0;
            in >> stored_phase_count >> stored_phase_species >> stored_phase_nx;
            if (key != "spatial_phase" ||
                stored_phase_species != species_.size() ||
                stored_phase_nx != grid_.nx() ||
                (!reset && stored_phase_count != spatial_phase_bins_.size())) {
                throw std::runtime_error(
                    "checkpoint spatial-phase contract does not match config");
            }
            std::size_t stored_samples_per_cycle = 0;
            std::size_t stored_samples_per_bin = 0;
            if (stored_phase_count != 0) {
                if (!(rf_frequency > 0.0) || interval == 0) {
                    throw std::runtime_error(
                        "checkpoint phase bins require a valid RF contract");
                }
                const double stored_steps_per_cycle_value =
                    1.0 / (rf_frequency * cfg_.dt);
                if (!std::isfinite(stored_steps_per_cycle_value) ||
                    stored_steps_per_cycle_value < 1.0) {
                    throw std::runtime_error(
                        "checkpoint phase-bin RF period is invalid");
                }
                const auto stored_steps_per_cycle =
                    static_cast<std::size_t>(std::llround(
                        stored_steps_per_cycle_value));
                stored_samples_per_cycle = stored_steps_per_cycle / interval;
                if (stored_samples_per_cycle % stored_phase_count != 0) {
                    throw std::runtime_error(
                        "checkpoint phase bins do not divide its RF cycle");
                }
                stored_samples_per_bin =
                    stored_samples_per_cycle / stored_phase_count;
            }
            std::size_t phase_sample_total = 0;
            for (std::size_t phase = 0; phase < stored_phase_count; ++phase) {
                std::size_t stored_phase = 0;
                std::size_t stored_phase_samples = 0;
                in >> key >> stored_phase >> stored_phase_samples;
                if (key != "phase_bin" || stored_phase != phase ||
                    stored_phase_samples > stored_samples) {
                    throw std::runtime_error(
                        "checkpoint phase-bin metadata is invalid");
                }
                const std::size_t complete_cycles =
                    stored_samples / stored_samples_per_cycle;
                const std::size_t remainder =
                    stored_samples % stored_samples_per_cycle;
                const std::size_t bin_start = phase * stored_samples_per_bin;
                const std::size_t partial = remainder > bin_start
                    ? std::min(stored_samples_per_bin, remainder - bin_start)
                    : 0;
                if (stored_phase_samples !=
                    complete_cycles * stored_samples_per_bin + partial) {
                    throw std::runtime_error(
                        "checkpoint phase-bin sample count is inconsistent");
                }
                phase_sample_total += stored_phase_samples;
                if (!reset) spatial_phase_bins_[phase].samples = stored_phase_samples;
                for (std::size_t species_id = 0;
                     species_id < stored_phase_species; ++species_id) {
                    std::size_t row_phase = 0;
                    std::size_t stored_species_id = 0;
                    std::string stored_name;
                    in >> key >> row_phase >> stored_species_id >> stored_name;
                    if (key != "phase_species" || row_phase != phase ||
                        stored_species_id != species_id ||
                        stored_name != species_[species_id].name()) {
                        throw std::runtime_error(
                            "checkpoint phase species metadata is invalid");
                    }
                    for (std::size_t node = 0; node < stored_phase_nx; ++node) {
                        double density = 0.0, vx = 0.0, vy = 0.0, vz = 0.0;
                        double energy = 0.0;
                        in >> density >> vx >> vy >> vz >> energy;
                        if (!std::isfinite(density) || density < 0.0 ||
                            !std::isfinite(vx) || !std::isfinite(vy) ||
                            !std::isfinite(vz) || !std::isfinite(energy) ||
                            energy < 0.0) {
                            throw std::runtime_error(
                                "checkpoint phase velocity moments are invalid");
                        }
                        if (!reset) {
                            auto& bin = spatial_phase_bins_[phase];
                            bin.density[species_id][node] = density;
                            bin.velocity_x_density[species_id][node] = vx;
                            bin.velocity_y_density[species_id][node] = vy;
                            bin.velocity_z_density[species_id][node] = vz;
                            bin.kinetic_energy_density[species_id][node] = energy;
                        }
                    }
                }
                std::size_t field_phase = 0;
                in >> key >> field_phase;
                if (key != "phase_fields" || field_phase != phase) {
                    throw std::runtime_error(
                        "checkpoint phase field metadata is invalid");
                }
                for (std::size_t node = 0; node < stored_phase_nx; ++node) {
                    double potential = 0.0, electric = 0.0, squared = 0.0;
                    in >> potential >> electric >> squared;
                    if (!std::isfinite(potential) || !std::isfinite(electric) ||
                        !std::isfinite(squared) || squared < 0.0) {
                        throw std::runtime_error(
                            "checkpoint phase field moments are invalid");
                    }
                    if (!reset) {
                        spatial_phase_bins_[phase].potential[node] = potential;
                        spatial_phase_bins_[phase].electric[node] = electric;
                        spatial_phase_bins_[phase].electric_squared[node] = squared;
                    }
                }
            }
            if (phase_sample_total != (stored_phase_count == 0 ? 0 : stored_samples)) {
                throw std::runtime_error(
                    "checkpoint phase-bin sample counts are incomplete");
            }
            in >> key;
        } else if (!reset && !spatial_phase_bins_.empty()) {
            throw std::runtime_error(
                "legacy checkpoint cannot restore phase-resolved spatial moments");
        }
        std::size_t stored_collision_steps = 0;
        std::size_t stored_collision_channels = 0;
        std::size_t stored_collision_nx = 0;
        std::size_t stored_collision_phases = 0;
        if (checkpoint_v11 || checkpoint_v12 || checkpoint_v13) {
            in >> stored_collision_steps >> stored_collision_channels >>
                stored_collision_nx >> stored_collision_phases;
            const std::size_t window_steps = stored_enabled
                ? end_step - start_step + 1
                : 0;
            if (key != "spatial_collision" ||
                stored_collision_nx != grid_.nx() ||
                stored_collision_steps > window_steps ||
                stored_collision_channels !=
                    (stored_enabled
                         ? collision_totals_.channel_names.size()
                         : 0) ||
                stored_collision_phases != stored_phase_count ||
                (!reset &&
                 (stored_collision_channels !=
                      spatial_collision_energy_sums_.size() ||
                  stored_collision_phases !=
                      spatial_collision_phase_energy_sums_.size()))) {
                throw std::runtime_error(
                    "checkpoint spatial collision-energy contract is invalid");
            }
            spatial_collision_steps_ =
                reset ? 0 : stored_collision_steps;
            for (std::size_t channel = 0;
                 channel < stored_collision_channels; ++channel) {
                std::size_t stored_channel = 0;
                std::string stored_name;
                in >> key >> stored_channel >> stored_name;
                if (key != "spatial_collision_channel" ||
                    stored_channel != channel ||
                    stored_name != collision_totals_.channel_names[channel]) {
                    throw std::runtime_error(
                        "checkpoint spatial collision channel is invalid");
                }
                for (std::size_t node = 0;
                     node < stored_collision_nx; ++node) {
                    double value = 0.0;
                    in >> value;
                    if (!std::isfinite(value)) {
                        throw std::runtime_error(
                            "checkpoint spatial collision energy is invalid");
                    }
                    if (!reset) {
                        spatial_collision_energy_sums_[channel][node] = value;
                    }
                }
            }
            in >> key;
            if (key != "spatial_collision_phase_steps") {
                throw std::runtime_error(
                    "checkpoint collision phase counts are missing");
            }
            std::size_t phase_step_total = 0;
            for (std::size_t phase = 0;
                 phase < stored_collision_phases; ++phase) {
                std::size_t count = 0;
                in >> count;
                if (count > stored_collision_steps) {
                    throw std::runtime_error(
                        "checkpoint collision phase count is invalid");
                }
                phase_step_total += count;
                if (!reset) spatial_collision_phase_steps_[phase] = count;
            }
            if (phase_step_total !=
                (stored_collision_phases == 0
                     ? 0
                     : stored_collision_steps)) {
                throw std::runtime_error(
                    "checkpoint collision phase counts are incomplete");
            }
            for (std::size_t phase = 0;
                 phase < stored_collision_phases; ++phase) {
                for (std::size_t channel = 0;
                     channel < stored_collision_channels; ++channel) {
                    std::size_t stored_phase = 0;
                    std::size_t stored_channel = 0;
                    in >> key >> stored_phase >> stored_channel;
                    if (key != "spatial_collision_phase_channel" ||
                        stored_phase != phase || stored_channel != channel) {
                        throw std::runtime_error(
                            "checkpoint phase collision channel is invalid");
                    }
                    for (std::size_t node = 0;
                         node < stored_collision_nx; ++node) {
                        double value = 0.0;
                        in >> value;
                        if (!std::isfinite(value)) {
                            throw std::runtime_error(
                                "checkpoint phase collision energy is invalid");
                        }
                        if (!reset) {
                            spatial_collision_phase_energy_sums_[phase]
                                [channel][node] = value;
                        }
                    }
                }
            }
            in >> key;
        } else if (!reset && stored_samples != 0 &&
                   !spatial_collision_energy_sums_.empty()) {
            throw std::runtime_error(
                "legacy checkpoint cannot restore spatial collision energy");
        }
        if (checkpoint_v15) {
            std::size_t stored_event_channels = 0;
            std::size_t stored_event_nx = 0;
            std::size_t stored_event_phases = 0;
            in >> stored_event_channels >> stored_event_nx >>
                stored_event_phases;
            if (key != "spatial_collision_events" ||
                stored_event_channels != stored_collision_channels ||
                stored_event_nx != grid_.nx() ||
                stored_event_phases != stored_collision_phases ||
                (!reset &&
                 (stored_event_channels !=
                      spatial_collision_event_sums_.size() ||
                  stored_event_phases !=
                      spatial_collision_phase_event_sums_.size()))) {
                throw std::runtime_error(
                    "checkpoint spatial collision-event contract is invalid");
            }
            for (std::size_t channel = 0;
                 channel < stored_event_channels; ++channel) {
                std::size_t stored_channel = 0;
                std::string stored_name;
                in >> key >> stored_channel >> stored_name;
                if (key != "spatial_collision_event_channel" ||
                    stored_channel != channel ||
                    stored_name != collision_totals_.channel_names[channel]) {
                    throw std::runtime_error(
                        "checkpoint spatial collision-event channel is invalid");
                }
                for (std::size_t node = 0; node < stored_event_nx; ++node) {
                    double value = 0.0;
                    in >> value;
                    if (!std::isfinite(value) || value < 0.0) {
                        throw std::runtime_error(
                            "checkpoint spatial collision-event value is invalid");
                    }
                    if (!reset) {
                        spatial_collision_event_sums_[channel][node] = value;
                    }
                }
            }
            for (std::size_t phase = 0;
                 phase < stored_event_phases; ++phase) {
                for (std::size_t channel = 0;
                     channel < stored_event_channels; ++channel) {
                    std::size_t stored_phase = 0;
                    std::size_t stored_channel = 0;
                    in >> key >> stored_phase >> stored_channel;
                    if (key != "spatial_collision_phase_event_channel" ||
                        stored_phase != phase || stored_channel != channel) {
                        throw std::runtime_error(
                            "checkpoint phase collision-event channel is invalid");
                    }
                    for (std::size_t node = 0; node < stored_event_nx; ++node) {
                        double value = 0.0;
                        in >> value;
                        if (!std::isfinite(value) || value < 0.0) {
                            throw std::runtime_error(
                                "checkpoint phase collision-event value is invalid");
                        }
                        if (!reset) {
                            spatial_collision_phase_event_sums_[phase]
                                [channel][node] = value;
                        }
                    }
                }
            }
            in >> key;
        } else if (!reset && stored_samples != 0 &&
                   !spatial_collision_event_sums_.empty()) {
            throw std::runtime_error(
                "legacy checkpoint cannot restore spatial collision events");
        }
        if (checkpoint_v12 || checkpoint_v13) {
            int stored_eedf_enabled = 0;
            int stored_eedf_history_enabled = 0;
            std::string stored_eedf_species;
            std::size_t stored_eedf_bins = 0;
            double stored_eedf_max = 0.0;
            double stored_eedf_tail_threshold = 0.0;
            std::size_t stored_eedf_phases = 0;
            std::size_t stored_eedf_regions = 0;
            in >> stored_eedf_enabled;
            if (checkpoint_v19) in >> stored_eedf_history_enabled;
            in >> stored_eedf_species >> stored_eedf_bins >> stored_eedf_max;
            if (checkpoint_v18) in >> stored_eedf_tail_threshold;
            in >> stored_eedf_phases >> stored_eedf_regions;
            const bool enabled_shape = stored_eedf_enabled == 1;
            const bool tail_contract_valid =
                !checkpoint_v18 ||
                (std::isfinite(stored_eedf_tail_threshold) &&
                 stored_eedf_tail_threshold >= 0.0 &&
                 (enabled_shape
                      ? stored_eedf_tail_threshold < stored_eedf_max
                      : stored_eedf_tail_threshold == 0.0));
            const bool shape_valid =
                (stored_eedf_enabled == 0 || stored_eedf_enabled == 1) &&
                (stored_eedf_history_enabled == 0 ||
                 stored_eedf_history_enabled == 1) &&
                (stored_eedf_history_enabled == 0 || enabled_shape) &&
                std::isfinite(stored_eedf_max) && stored_eedf_max >= 0.0 &&
                tail_contract_valid &&
                stored_eedf_phases ==
                    (enabled_shape ? stored_phase_count : 0) &&
                (!enabled_shape ||
                 (stored_eedf_bins > 0 && stored_eedf_max > 0.0 &&
                  stored_eedf_regions > 0));
            const bool contract_matches =
                stored_eedf_enabled == (cfg_.phase_eedf.enabled ? 1 : 0) &&
                stored_eedf_history_enabled ==
                    (cfg_.phase_eedf.history_enabled ? 1 : 0) &&
                stored_eedf_species ==
                    (cfg_.phase_eedf.species.empty()
                         ? "-" : cfg_.phase_eedf.species) &&
                stored_eedf_bins == cfg_.phase_eedf.energy_bins &&
                stored_eedf_max == cfg_.phase_eedf.energy_max &&
                stored_eedf_tail_threshold ==
                    cfg_.phase_eedf.tail_threshold &&
                stored_eedf_phases == phase_eedf_accumulators_.size() &&
                stored_eedf_regions == cfg_.phase_eedf.regions.size();
            if (key != "phase_eedf" || !shape_valid ||
                (!reset && (!contract_matches ||
                            (stored_eedf_enabled == 1 && !checkpoint_v18) ||
                            (stored_eedf_history_enabled == 1 &&
                             !checkpoint_v19)))) {
                throw std::runtime_error(
                    "checkpoint phase EEDF contract is invalid");
            }
            std::vector<PhaseEedfRegion1DConfig> stored_regions(
                stored_eedf_regions);
            for (std::size_t region = 0;
                 region < stored_eedf_regions; ++region) {
                std::size_t stored_region = 0;
                in >> key >> stored_region >> stored_regions[region].name >>
                    stored_regions[region].x_min >> stored_regions[region].x_max;
                if (key != "phase_eedf_region" || stored_region != region ||
                    stored_regions[region].name.empty() ||
                    !std::isfinite(stored_regions[region].x_min) ||
                    !std::isfinite(stored_regions[region].x_max) ||
                    stored_regions[region].x_max <=
                        stored_regions[region].x_min ||
                    (!reset &&
                     (stored_regions[region].name !=
                          cfg_.phase_eedf.regions[region].name ||
                      stored_regions[region].x_min !=
                          cfg_.phase_eedf.regions[region].x_min ||
                      stored_regions[region].x_max !=
                          cfg_.phase_eedf.regions[region].x_max))) {
                    throw std::runtime_error(
                        "checkpoint phase EEDF region is invalid");
                }
            }
            for (std::size_t phase = 0;
                 phase < stored_eedf_phases; ++phase) {
                for (std::size_t region = 0;
                     region < stored_eedf_regions; ++region) {
                    std::size_t stored_phase = 0;
                    std::size_t stored_region = 0;
                    PhaseEedfAccumulator1D value;
                    in >> key >> stored_phase >> stored_region >>
                        value.macro_observations >>
                        value.overflow_macro_observations >>
                        value.represented_observations >>
                        value.overflow_represented_observations >>
                        value.weighted_energy_sum >>
                        value.weighted_energy_squared_sum >>
                        value.weighted_velocity_x_sum >>
                        value.weighted_velocity_y_sum >>
                        value.weighted_velocity_z_sum;
                    if (checkpoint_v18) {
                        in >> value.weighted_velocity_x_squared_sum >>
                            value.weighted_velocity_y_squared_sum >>
                            value.weighted_velocity_z_squared_sum >>
                            value.tail_represented_observations >>
                            value.tail_positive_x_represented_observations >>
                            value.tail_negative_x_represented_observations >>
                            value.tail_weighted_velocity_x_sum >>
                            value.tail_weighted_velocity_x_squared_sum >>
                            value.tail_weighted_transverse_velocity_squared_sum;
                    }
                    if (checkpoint_v19) {
                        in >> value.tail_weighted_age_steps_sum >>
                            value.tail_weighted_energetic_steps_sum >>
                            value.tail_weighted_energetic_duty_fraction_sum >>
                            value
                                .tail_weighted_consecutive_energetic_steps_sum >>
                            value.tail_weighted_entries_sum >>
                            value.tail_weighted_elastic_collisions_sum >>
                            value.tail_weighted_excitation_collisions_sum >>
                            value.tail_weighted_ionization_collisions_sum >>
                            value
                                .tail_weighted_charge_exchange_collisions_sum >>
                            value.tail_weighted_bgk_collisions_sum >>
                            value
                                .tail_born_during_window_represented_observations;
                    }
                    value.histogram.resize(stored_eedf_bins);
                    for (auto& count : value.histogram) in >> count;
                    const bool finite = std::isfinite(
                            value.represented_observations) &&
                        std::isfinite(value.overflow_represented_observations) &&
                        std::isfinite(value.weighted_energy_sum) &&
                        std::isfinite(value.weighted_energy_squared_sum) &&
                        std::isfinite(value.weighted_velocity_x_sum) &&
                        std::isfinite(value.weighted_velocity_y_sum) &&
                        std::isfinite(value.weighted_velocity_z_sum) &&
                        std::isfinite(value.weighted_velocity_x_squared_sum) &&
                        std::isfinite(value.weighted_velocity_y_squared_sum) &&
                        std::isfinite(value.weighted_velocity_z_squared_sum) &&
                        std::isfinite(value.tail_represented_observations) &&
                        std::isfinite(
                            value.tail_positive_x_represented_observations) &&
                        std::isfinite(
                            value.tail_negative_x_represented_observations) &&
                        std::isfinite(value.tail_weighted_velocity_x_sum) &&
                        std::isfinite(
                            value.tail_weighted_velocity_x_squared_sum) &&
                        std::isfinite(
                            value.tail_weighted_transverse_velocity_squared_sum) &&
                        std::isfinite(value.tail_weighted_age_steps_sum) &&
                        std::isfinite(value.tail_weighted_energetic_steps_sum) &&
                        std::isfinite(
                            value.tail_weighted_energetic_duty_fraction_sum) &&
                        std::isfinite(value
                            .tail_weighted_consecutive_energetic_steps_sum) &&
                        std::isfinite(value.tail_weighted_entries_sum) &&
                        std::isfinite(
                            value.tail_weighted_elastic_collisions_sum) &&
                        std::isfinite(
                            value.tail_weighted_excitation_collisions_sum) &&
                        std::isfinite(
                            value.tail_weighted_ionization_collisions_sum) &&
                        std::isfinite(value
                            .tail_weighted_charge_exchange_collisions_sum) &&
                        std::isfinite(
                            value.tail_weighted_bgk_collisions_sum) &&
                        std::isfinite(value
                            .tail_born_during_window_represented_observations) &&
                        std::all_of(value.histogram.begin(),
                            value.histogram.end(), [](double count) {
                                return std::isfinite(count) && count >= 0.0;
                            });
                    if (key != "phase_eedf_accumulator" ||
                        stored_phase != phase || stored_region != region ||
                        value.overflow_macro_observations >
                            value.macro_observations || !finite ||
                        value.represented_observations < 0.0 ||
                        value.overflow_represented_observations < 0.0 ||
                        value.overflow_represented_observations >
                            value.represented_observations ||
                        value.weighted_energy_sum < 0.0 ||
                        value.weighted_energy_squared_sum < 0.0 ||
                        value.weighted_velocity_x_squared_sum < 0.0 ||
                        value.weighted_velocity_y_squared_sum < 0.0 ||
                        value.weighted_velocity_z_squared_sum < 0.0 ||
                        value.tail_represented_observations < 0.0 ||
                        value.tail_positive_x_represented_observations < 0.0 ||
                        value.tail_negative_x_represented_observations < 0.0 ||
                        value.tail_positive_x_represented_observations +
                                value.tail_negative_x_represented_observations >
                            value.tail_represented_observations *
                                (1.0 + 1e-12) ||
                        value.tail_weighted_velocity_x_squared_sum < 0.0 ||
                        value.tail_weighted_transverse_velocity_squared_sum <
                            0.0 ||
                        value.tail_weighted_age_steps_sum < 0.0 ||
                        value.tail_weighted_energetic_steps_sum < 0.0 ||
                        value.tail_weighted_energetic_duty_fraction_sum < 0.0 ||
                        value.tail_weighted_energetic_duty_fraction_sum >
                            value.tail_represented_observations *
                                (1.0 + 1e-12) ||
                        value.tail_weighted_consecutive_energetic_steps_sum <
                            0.0 ||
                        value.tail_weighted_entries_sum < 0.0 ||
                        value.tail_weighted_elastic_collisions_sum < 0.0 ||
                        value.tail_weighted_excitation_collisions_sum < 0.0 ||
                        value.tail_weighted_ionization_collisions_sum < 0.0 ||
                        value.tail_weighted_charge_exchange_collisions_sum <
                            0.0 ||
                        value.tail_weighted_bgk_collisions_sum < 0.0 ||
                        value
                                .tail_born_during_window_represented_observations <
                            0.0 ||
                        value
                                .tail_born_during_window_represented_observations >
                            value.tail_represented_observations) {
                        throw std::runtime_error(
                            "checkpoint phase EEDF accumulator is invalid");
                    }
                    if (!reset) {
                        phase_eedf_accumulators_[phase][region] =
                            std::move(value);
                    }
                }
            }
            in >> key;
            if (checkpoint_v20) {
                int stored_crossings_enabled = 0;
                std::size_t stored_crossing_phases = 0;
                std::size_t stored_crossing_regions = 0;
                std::size_t stored_crossing_processes = 0;
                in >> stored_crossings_enabled >> stored_crossing_phases >>
                    stored_crossing_regions >> stored_crossing_processes;
                const bool crossings_enabled =
                    stored_crossings_enabled == 1;
                const bool crossing_shape_valid =
                    (stored_crossings_enabled == 0 ||
                     stored_crossings_enabled == 1) &&
                    stored_crossing_processes == 6 &&
                    stored_crossing_phases ==
                        (crossings_enabled ? stored_phase_count : 0) &&
                    stored_crossing_regions ==
                        (crossings_enabled ? stored_eedf_regions : 0);
                const bool crossing_contract_matches =
                    stored_crossings_enabled ==
                        (cfg_.phase_eedf.history_enabled ? 1 : 0) &&
                    stored_crossing_phases ==
                        phase_eedf_threshold_crossings_.size() &&
                    stored_crossing_regions ==
                        (phase_eedf_threshold_crossings_.empty()
                            ? 0
                            : phase_eedf_threshold_crossings_.front().size());
                if (key != "phase_eedf_threshold_crossings" ||
                    !crossing_shape_valid ||
                    (!reset && !crossing_contract_matches)) {
                    throw std::runtime_error(
                        "checkpoint phase EEDF threshold-crossing contract "
                        "is invalid");
                }
                for (std::size_t phase = 0;
                     phase < stored_crossing_phases; ++phase) {
                    for (std::size_t region = 0;
                         region < stored_crossing_regions; ++region) {
                        std::size_t stored_phase = 0;
                        std::size_t stored_region = 0;
                        PhaseEedfThresholdCrossingAccumulator1D value;
                        in >> key >> stored_phase >> stored_region >>
                            value.electron_time_macro_observations >>
                            value.energetic_time_macro_observations >>
                            value.interstep_promotions >>
                            value.interstep_demotions;
                        if (checkpoint_v21) {
                            in >> value.field_push_promotions >>
                                value.field_push_demotions;
                        }
                        for (auto& count : value.collision_promotions) {
                            in >> count;
                        }
                        for (auto& count : value.collision_demotions) {
                            in >> count;
                        }
                        in >> value.energetic_births >>
                            value.subthreshold_births;
                        const bool counts_valid =
                            value.energetic_time_macro_observations <=
                                value.electron_time_macro_observations &&
                            value.interstep_promotions <=
                                value.electron_time_macro_observations &&
                            value.interstep_demotions <=
                                value.electron_time_macro_observations &&
                            value.field_push_promotions <=
                                value.electron_time_macro_observations &&
                            value.field_push_demotions <=
                                value.electron_time_macro_observations;
                        if (key !=
                                "phase_eedf_threshold_crossing_accumulator" ||
                            stored_phase != phase || stored_region != region ||
                            !counts_valid) {
                            throw std::runtime_error(
                                "checkpoint phase EEDF threshold-crossing "
                                "accumulator is invalid");
                        }
                        if (!reset) {
                            phase_eedf_threshold_crossings_[phase][region] =
                                value;
                        }
                    }
                }
                if (!checkpoint_v21 && !reset && crossings_enabled) {
                    throw std::runtime_error(
                        "legacy checkpoint cannot restore phase EEDF "
                        "field-push threshold crossings");
                }
                in >> key;
            } else if (!reset && cfg_.phase_eedf.history_enabled) {
                throw std::runtime_error(
                    "legacy checkpoint cannot restore phase EEDF "
                    "threshold crossings");
            }
        } else if (!reset && cfg_.phase_eedf.enabled && stored_samples != 0) {
            throw std::runtime_error(
                "legacy checkpoint cannot restore phase EEDF state");
        }
        if (checkpoint_v17) {
            int stored_enabled = 0;
            std::string stored_species;
            std::size_t stored_bins = 0;
            double stored_max = 0.0;
            std::size_t stored_phases = 0;
            std::size_t stored_surfaces = 0;
            in >> stored_enabled >> stored_species >> stored_bins >> stored_max >>
                stored_phases >> stored_surfaces;
            const bool enabled_shape = stored_enabled == 1;
            const bool shape_valid =
                (stored_enabled == 0 || stored_enabled == 1) &&
                std::isfinite(stored_max) && stored_max >= 0.0 &&
                stored_phases == (enabled_shape ? stored_phase_count : 0) &&
                (!enabled_shape || (stored_bins > 0 && stored_max > 0.0 &&
                                    stored_surfaces > 0));
            const auto& configured = cfg_.phase_surface_flux;
            const bool contract_matches =
                stored_enabled == (configured.enabled ? 1 : 0) &&
                stored_species ==
                    (configured.species.empty() ? "-" : configured.species) &&
                stored_bins == configured.energy_bins &&
                stored_max == configured.energy_max &&
                stored_phases == phase_surface_flux_accumulators_.size() &&
                stored_surfaces == configured.positions.size();
            if (key != "phase_surface_flux" || !shape_valid ||
                (!configured.reset_on_restart && !contract_matches)) {
                throw std::runtime_error(
                    "checkpoint phase surface-flux contract is invalid");
            }
            std::vector<double> stored_positions(stored_surfaces);
            for (std::size_t surface = 0; surface < stored_surfaces; ++surface) {
                std::size_t stored_surface = 0;
                in >> key >> stored_surface >> stored_positions[surface];
                if (key != "phase_surface_flux_position" ||
                    stored_surface != surface ||
                    !std::isfinite(stored_positions[surface]) ||
                    (!configured.reset_on_restart &&
                     stored_positions[surface] != configured.positions[surface])) {
                    throw std::runtime_error(
                        "checkpoint phase surface-flux position is invalid");
                }
            }
            for (std::size_t phase = 0; phase < stored_phases; ++phase) {
                for (std::size_t surface = 0;
                     surface < stored_surfaces; ++surface) {
                    for (std::size_t direction = 0; direction < 2; ++direction) {
                        std::size_t stored_phase = 0;
                        std::size_t stored_surface = 0;
                        std::size_t stored_direction = 0;
                        PhaseSurfaceFluxAccumulator1D value;
                        in >> key >> stored_phase >> stored_surface >>
                            stored_direction >> value.macro_crossings >>
                            value.overflow_macro_crossings >>
                            value.represented_crossings >>
                            value.overflow_represented_crossings >>
                            value.represented_kinetic_energy;
                        value.represented_histogram.resize(stored_bins);
                        for (auto& count : value.represented_histogram) {
                            in >> count;
                        }
                        const bool finite =
                            std::isfinite(value.represented_crossings) &&
                            std::isfinite(value.overflow_represented_crossings) &&
                            std::isfinite(value.represented_kinetic_energy) &&
                            std::all_of(
                                value.represented_histogram.begin(),
                                value.represented_histogram.end(),
                                [](double count) {
                                    return std::isfinite(count) && count >= 0.0;
                                });
                        if (key != "phase_surface_flux_accumulator" ||
                            stored_phase != phase ||
                            stored_surface != surface ||
                            stored_direction != direction || !finite ||
                            value.overflow_macro_crossings >
                                value.macro_crossings ||
                            value.represented_crossings < 0.0 ||
                            value.overflow_represented_crossings < 0.0 ||
                            value.overflow_represented_crossings >
                                value.represented_crossings ||
                            value.represented_kinetic_energy < 0.0) {
                            throw std::runtime_error(
                                "checkpoint phase surface-flux accumulator is invalid");
                        }
                        if (!configured.reset_on_restart) {
                            phase_surface_flux_accumulators_[phase]
                                [surface][direction] = std::move(value);
                        }
                    }
                }
            }
            in >> key;
        } else if (cfg_.phase_surface_flux.enabled &&
                   !cfg_.phase_surface_flux.reset_on_restart &&
                   stored_samples != 0) {
            throw std::runtime_error(
                "legacy checkpoint cannot restore phase surface-flux state");
        }
        if (reset) {
            spatial_moment_samples_ = 0;
            spatial_collision_steps_ = 0;
            std::fill(spatial_collision_phase_steps_.begin(),
                      spatial_collision_phase_steps_.end(), 0);
            for (auto& channels : spatial_collision_phase_energy_sums_) {
                for (auto& values : channels) {
                    std::fill(values.begin(), values.end(), 0.0);
                }
            }
            for (auto& values : spatial_collision_energy_sums_) {
                std::fill(values.begin(), values.end(), 0.0);
            }
            for (auto& channels : spatial_collision_phase_event_sums_) {
                for (auto& values : channels) {
                    std::fill(values.begin(), values.end(), 0.0);
                }
            }
            for (auto& values : spatial_collision_event_sums_) {
                std::fill(values.begin(), values.end(), 0.0);
            }
            for (auto& sum : spatial_kinetic_energy_sums_) {
                std::fill(sum.begin(), sum.end(), 0.0);
            }
            std::fill(spatial_potential_sums_.begin(),
                      spatial_potential_sums_.end(), 0.0);
            std::fill(spatial_electric_sums_.begin(),
                      spatial_electric_sums_.end(), 0.0);
            std::fill(spatial_electric_squared_sums_.begin(),
                      spatial_electric_squared_sums_.end(), 0.0);
        }
    } else if (
        cfg_.spatial_average.enabled &&
        !cfg_.spatial_average.reset_on_restart) {
        throw std::runtime_error(
            "legacy checkpoint cannot restore enabled 1D "
            "spatial averaging");
    }
    in >> step_;
    if (key != "step") throw std::runtime_error("checkpoint missing step");
    if (checkpoint_v14 && wall_impact_origin_step_ > step_) {
        throw std::runtime_error(
            "checkpoint wall-impact origin exceeds its step");
    }
    if (checkpoint_v14 &&
        cfg_.wall_impact_spectrum.enabled &&
        cfg_.wall_impact_spectrum.reset_on_restart) {
        wall_impact_origin_step_ = step_;
        for (std::size_t species_id = 0;
             species_id < wall_impact_spectra_.size(); ++species_id) {
            auto& spectrum = wall_impact_spectra_[species_id];
            spectrum.baseline_loss =
                species_boundary_losses_[species_id];
            clear_wall_impact_side(spectrum.left);
            clear_wall_impact_side(spectrum.right);
        }
    }
    if ((checkpoint_v11 || checkpoint_v12 || checkpoint_v13) &&
        cfg_.spatial_average.enabled &&
        !cfg_.spatial_average.reset_on_restart &&
        !spatial_collision_energy_sums_.empty()) {
        const auto& average = cfg_.spatial_average;
        const std::size_t expected_collision_steps =
            step_ < average.start_step
                ? 0
                : std::min(step_, average.end_step) -
                      average.start_step + 1;
        if (spatial_collision_steps_ != expected_collision_steps) {
            throw std::runtime_error(
                "checkpoint spatial collision timestep count is inconsistent");
        }
        if (!spatial_collision_phase_steps_.empty()) {
            const auto steps_per_cycle = static_cast<std::size_t>(std::llround(
                1.0 / (average.rf_frequency * cfg_.dt)));
            const auto steps_per_phase =
                steps_per_cycle / spatial_collision_phase_steps_.size();
            const auto complete_cycles =
                expected_collision_steps / steps_per_cycle;
            const auto remainder =
                expected_collision_steps % steps_per_cycle;
            for (std::size_t phase = 0;
                 phase < spatial_collision_phase_steps_.size(); ++phase) {
                const auto phase_start = phase * steps_per_phase;
                const auto partial = remainder > phase_start
                    ? std::min(steps_per_phase, remainder - phase_start)
                    : 0;
                if (spatial_collision_phase_steps_[phase] !=
                    complete_cycles * steps_per_phase + partial) {
                    throw std::runtime_error(
                        "checkpoint collision phase timestep count is inconsistent");
                }
            }
        }
    }
    if (legacy_boundary_loss_origin) {
        boundary_loss_origin_step_ = step_;
    }
    if (legacy_wall_impact_origin) {
        wall_impact_origin_step_ = step_;
        for (std::size_t species_id = 0;
             species_id < wall_impact_spectra_.size(); ++species_id) {
            wall_impact_spectra_[species_id].baseline_loss =
                species_boundary_losses_[species_id];
        }
    }
    if (legacy_power_transfer_origin) {
        power_transfer_origin_step_ = step_;
    }
    if (cfg_.spatial_average.enabled) {
        const auto& average = cfg_.spatial_average;
        if (average.reset_on_restart &&
            average.start_step <= step_) {
            throw std::runtime_error(
                "restart-reset spatial average must start after "
                "the checkpoint step");
        }
        std::size_t expected_at_step = 0;
        if (step_ >= average.start_step) {
            const std::size_t last =
                std::min(step_, average.end_step);
            expected_at_step =
                1 + (last - average.start_step) /
                        average.interval;
        }
        if (spatial_average_samples_ != expected_at_step) {
            throw std::runtime_error(
                "checkpoint spatial-average sample count does not "
                "match its step");
        }
    }
    in >> key >> time_;
    if (key != "time") throw std::runtime_error("checkpoint missing time");
    std::size_t species_count = 0;
    in >> key >> species_count;
    if (key != "species_count" || species_count != species_.size()) throw std::runtime_error("checkpoint species count does not match config");
    in >> key;
    if (key != "rng") throw std::runtime_error("checkpoint missing rng state");
    in >> rng_;

    for (std::size_t species_id = 0; species_id < species_.size(); ++species_id) {
        std::size_t stored_species_id = 0;
        std::string stored_name;
        std::size_t particle_count = 0;
        in >> key >> stored_species_id >> stored_name >> particle_count;
        if (key != "species" || stored_species_id != species_id || stored_name != species_[species_id].name()) {
            throw std::runtime_error("checkpoint species metadata does not match config");
        }
        if (particle_count >
            cfg_.max_particles_per_species) {
            throw std::runtime_error(
                "checkpoint particle count exceeds "
                "max_particles_per_species");
        }
        auto& particles = species_[species_id].particles();
        particles.resize(particle_count);
        for (auto& p : particles) {
            int alive = 0;
            if (checkpoint_v4_state) {
                in >> p.x >> p.v >>
                    p.velocity_y >> p.velocity_z >>
                    p.v_half >> alive;
            } else {
                in >> p.x >> p.v >> p.v_half >> alive;
                p.velocity_y = 0.0;
                p.velocity_z = 0.0;
            }
            p.alive = alive != 0;
        }
    }
    if (checkpoint_v19) {
        int stored_history_enabled = 0;
        std::size_t stored_history_species = 0;
        std::size_t stored_history_count = 0;
        in >> key >> stored_history_enabled >> stored_history_species >>
            stored_history_count;
        const bool enabled = stored_history_enabled == 1;
        const bool header_valid =
            key == "phase_eedf_particle_history" &&
            (stored_history_enabled == 0 || stored_history_enabled == 1) &&
            (!enabled ||
             (stored_history_species < species_.size() &&
              stored_history_count == species_[stored_history_species]
                  .particles().size() &&
              (cfg_.spatial_average.reset_on_restart ||
               stored_history_species == phase_eedf_species_id_))) &&
            (enabled || stored_history_count == 0) &&
            (cfg_.spatial_average.reset_on_restart ||
             enabled == cfg_.phase_eedf.history_enabled);
        if (!header_valid) {
            throw std::runtime_error(
                "checkpoint phase EEDF particle-history header is invalid");
        }
        std::vector<ParticleHistory1D> stored_histories(
            stored_history_count);
        for (std::size_t particle_id = 0;
             particle_id < stored_history_count; ++particle_id) {
            std::size_t stored_particle_id = 0;
            int born = 0;
            int energetic_previous = 0;
            auto& history = stored_histories[particle_id];
            in >> key >> stored_particle_id >> history.age_steps >>
                history.energetic_steps >>
                history.consecutive_energetic_steps >>
                history.tail_entries >> history.elastic_collisions >>
                history.excitation_collisions >>
                history.ionization_collisions >>
                history.charge_exchange_collisions >>
                history.bgk_collisions >> born >> energetic_previous;
            const bool valid =
                key == "phase_eedf_particle_history_entry" &&
                stored_particle_id == particle_id &&
                (born == 0 || born == 1) &&
                (energetic_previous == 0 || energetic_previous == 1) &&
                history.energetic_steps <= history.age_steps &&
                history.consecutive_energetic_steps <= history.age_steps &&
                history.tail_entries <= history.energetic_steps;
            const bool state_consistent =
                (energetic_previous == 1) ==
                    (history.consecutive_energetic_steps > 0);
            if (!valid || !state_consistent) {
                throw std::runtime_error(
                    "checkpoint phase EEDF particle history is invalid");
            }
            history.born_during_window = born == 1;
            history.energetic_previous_step = energetic_previous == 1;
        }
        if (cfg_.phase_eedf.history_enabled) {
            if (cfg_.spatial_average.reset_on_restart) {
                phase_eedf_particle_histories_.assign(
                    species_[phase_eedf_species_id_].particles().size(), {});
            } else {
                phase_eedf_particle_histories_ =
                    std::move(stored_histories);
            }
        }
    } else if (cfg_.phase_eedf.history_enabled) {
        if (!cfg_.spatial_average.reset_on_restart) {
            throw std::runtime_error(
                "legacy checkpoint cannot restore phase EEDF particle history");
        }
        phase_eedf_particle_histories_.assign(
            species_[phase_eedf_species_id_].particles().size(), {});
    }
    require_stream(in, "failed while reading checkpoint: " + path.string());
    deposit_and_solve(time_);
    initialized_ = true;
}

RunSummary Simulation::run() {
    if (!cfg_.restart_path.empty()) load_checkpoint(cfg_.restart_path);
    else initialize();

    write_unit_metadata(cfg_.output_dir, cfg_.units, 1);
    if (!cfg_.initial_state_path.empty()) {
        write_external_particle_state_metadata(
            std::filesystem::path(cfg_.output_dir) /
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
        std::filesystem::path(cfg_.output_dir) /
            "initialization.csv",
        1,
        !cfg_.restart_path.empty()
            ? "restart"
            : (!cfg_.initial_state_path.empty()
                   ? "external"
                   : "generated"),
        initialization_moments);
    const auto initialization_acceptance =
        assess_initialization_acceptance(
            cfg_.initialization_acceptance,
            initialization_moments,
            cfg_.velocity_dimensions);
    write_initialization_acceptance_report(
        std::filesystem::path(cfg_.output_dir) /
            "initialization_acceptance.csv",
        initialization_acceptance);
    enforce_initialization_acceptance(
        initialization_acceptance);
    Diagnostics diag(
        cfg_.output_dir, species_, cfg_.units.permittivity());
    diag.write_header();
    auto s0 = diag.sample(step_, time_, grid_, species_);
    diag.write_sample(s0);
    diag.write_fields(step_, grid_);
    std::ofstream collision_output;
    const bool collisions_enabled =
        legacy_bgk_enabled_ || !mcc_models_.empty();
    if (collisions_enabled) {
        collision_output.open(
            std::filesystem::path(cfg_.output_dir) /
            "collisions.csv");
        if (!collision_output) {
            throw std::runtime_error(
                "cannot open collision diagnostics output");
        }
        write_collision_header(
            collision_output, collision_totals_,
            cfg_.units.system == UnitSystem::SI);
        write_collision_sample(
            collision_output, step_, time_,
            collision_interval_, collision_totals_);
        clear_collision_counts(collision_interval_);
    }
    std::ofstream boundary_loss_output(
        std::filesystem::path(cfg_.output_dir) /
        "boundary_losses.csv");
    if (!boundary_loss_output) {
        throw std::runtime_error(
            "cannot open 1D boundary-loss diagnostics output");
    }
    write_boundary_loss_header(
        boundary_loss_output, species_,
        cfg_.units.system == UnitSystem::SI);
    write_boundary_loss_sample(
        boundary_loss_output, step_, time_,
        boundary_loss_origin_step_, species_,
        species_boundary_losses_);
    std::ofstream power_transfer_output(
        std::filesystem::path(cfg_.output_dir) /
        "power_transfer.csv");
    if (!power_transfer_output) {
        throw std::runtime_error(
            "cannot open 1D power-transfer diagnostics output");
    }
    write_power_transfer_header(
        power_transfer_output, species_,
        cfg_.units.system == UnitSystem::SI);
    write_power_transfer_sample(
        power_transfer_output, step_, time_,
        power_transfer_origin_step_, species_,
        species_power_transfer_);
    if (cfg_.checkpoint_output) save_checkpoint(checkpoint_path_for_step(cfg_, step_));

    const std::size_t limit = cfg_.mode == RunMode::SteadyState ? cfg_.max_steps : cfg_.steps;
    RunSummary summary;
    summary.final_sample = s0;
    while (step_ < limit) {
        step();
        const bool at_output = step_ % cfg_.output_interval == 0 || step_ == limit;
        bool reached_steady = false;
        if (at_output) {
            auto s = diag.sample(step_, time_, grid_, species_);
            diag.write_sample(s);
            diag.write_fields(step_, grid_);
            if (collisions_enabled) {
                write_collision_sample(
                    collision_output, step_, time_,
                    collision_interval_, collision_totals_);
                clear_collision_counts(collision_interval_);
            }
            write_boundary_loss_sample(
                boundary_loss_output, step_, time_,
                boundary_loss_origin_step_, species_,
                species_boundary_losses_);
            write_power_transfer_sample(
                power_transfer_output, step_, time_,
                power_transfer_origin_step_, species_,
                species_power_transfer_);
            summary.final_sample = s;
            reached_steady = cfg_.mode == RunMode::SteadyState &&
                             adjacent_energy_windows_converged(diag.history(), cfg_.steady_window, cfg_.steady_tolerance);
            if (reached_steady) {
                summary.steady_state_reached = true;
            }
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
    write_wall_impact_spectrum();
    write_spatial_average();
    return summary;
}
}
