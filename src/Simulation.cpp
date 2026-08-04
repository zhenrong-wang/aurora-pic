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

void validate_runtime_config(const Config& cfg) {
    if (cfg.velocity_dimensions != 1 &&
        cfg.velocity_dimensions != 3) {
        throw std::invalid_argument(
            "simulation velocity_dimensions must be 1 or 3");
    }
    if (!std::isfinite(cfg.dt) || cfg.dt <= 0.0) throw std::invalid_argument("simulation dt must be positive and finite");
    if (cfg.output_interval == 0) throw std::invalid_argument("output_interval must be positive");
    validate_spatial_average_1d(cfg);
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
    if (cfg_.velocity_dimensions == 3 &&
        !cfg_.initial_state_path.empty()) {
        throw std::invalid_argument(
            "external particle-state initialization is not yet "
            "available for 1D3V");
    }
    validate_initialization_acceptance(
        cfg_.initialization_acceptance,
        "1D initialization acceptance config");
    if (cfg_.max_particles_per_species == 0) {
        throw std::invalid_argument(
            "max_particles_per_species must be positive");
    }
    for (const auto& sc : cfg_.species) {
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
        const auto& target =
            species_[runtime.species_id].config();
        for (std::size_t channel = 0;
             channel < collision.channels.size(); ++channel) {
            const auto& channel_config =
                collision.channels[channel];
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
    if (!cfg_.collision_models.empty()) {
        return "named_null_collision";
    }
    return to_string(cfg_.collisions.model);
}

double Simulation::electrode_potential(
    double offset,
    const SinusoidalVoltageConfig& drive,
    double field_time) const {
    return offset + drive.amplitude * std::sin(
        2.0 * std::numbers::pi * drive.frequency * field_time +
        drive.phase);
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
                    particle.v_half = record.velocity.x;
                    particle.alive = true;
                },
                cfg_.initial_state_signature);
    }
    deposit_and_solve(time_);
    for (auto& sp : species_) {
        const double qm = sp.charge() / sp.mass();
        auto& particles = sp.particles();
        runtime_parallel_for(std::size_t{0}, particles.size(), cfg_.runtime, [&](std::size_t particle_id) {
            auto& p = particles[particle_id];
            if (p.alive) initialize_leapfrog_half_step(p, interpolate_electric(grid_, p.x), qm, cfg_.dt);
        });
    }
    initialized_ = true;
}

void Simulation::apply_collisions() {
    if (legacy_bgk_enabled_ &&
        cfg_.collisions.frequency > 0.0) {
        const double probability =
            1.0 - std::exp(
                -cfg_.collisions.frequency * cfg_.dt);
        for (auto& species : species_) {
            std::uniform_real_distribution<double> unit(
                0.0, 1.0);
            std::normal_distribution<double> neutral_velocity(
                0.0,
                cfg_.collisions.neutral_temperature_velocity);
            const double charge_to_mass =
                species.charge() / species.mass();
            for (auto& particle : species.particles()) {
                if (!particle.alive ||
                    unit(rng_) >= probability) {
                    continue;
                }
                const double energy_before = 0.5 * species.mass() *
                    species.weight() *
                    (particle.v * particle.v +
                     (cfg_.velocity_dimensions == 3
                          ? particle.velocity_y * particle.velocity_y +
                            particle.velocity_z * particle.velocity_z
                          : 0.0));
                particle.v = neutral_velocity(rng_);
                if (cfg_.velocity_dimensions == 3) {
                    particle.velocity_y =
                        neutral_velocity(rng_);
                    particle.velocity_z =
                        neutral_velocity(rng_);
                }
                ++collision_totals_.candidates;
                ++collision_totals_.channel_collisions[0];
                ++collision_interval_.candidates;
                ++collision_interval_.channel_collisions[0];
                const double energy_after = 0.5 * species.mass() *
                    species.weight() *
                    (particle.v * particle.v +
                     (cfg_.velocity_dimensions == 3
                          ? particle.velocity_y * particle.velocity_y +
                            particle.velocity_z * particle.velocity_z
                          : 0.0));
                collision_totals_.channel_energy_change[0] +=
                    energy_after - energy_before;
                collision_interval_.channel_energy_change[0] +=
                    energy_after - energy_before;
                initialize_leapfrog_half_step(
                    particle,
                    interpolate_electric(
                        grid_, particle.x),
                    charge_to_mass, cfg_.dt);
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
        auto& sp = species_[runtime.species_id];
        const double qm = sp.charge() / sp.mass();
        const std::size_t initial_particle_count =
            sp.particles().size();
        for (std::size_t particle_id = 0;
             particle_id < initial_particle_count;
             ++particle_id) {
            auto& part = sp.particles()[particle_id];
            if (!part.alive) continue;
            CollisionStepStatistics statistics;
            if (cfg_.velocity_dimensions == 3) {
                Vec3 velocity{
                    part.v, part.velocity_y, part.velocity_z};
                statistics =
                    runtime.model->collide(
                        velocity, cfg_.dt, rng_);
                part.v = velocity.x;
                part.velocity_y = velocity.y;
                part.velocity_z = velocity.z;
            } else {
                statistics =
                    runtime.model->collide(
                        part.v, cfg_.dt, rng_);
            }
            for (auto& energy_change :
                 statistics.channel_projectile_energy_change) {
                energy_change *= sp.weight();
            }
            add_collision_statistics(
                collision_totals_, statistics,
                runtime.diagnostic_offset);
            add_collision_statistics(
                collision_interval_, statistics,
                runtime.diagnostic_offset);
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
            }
            if (statistics.primary_removal_channel) {
                throw std::logic_error(
                    "1D MCC produced an unsupported primary "
                    "removal event");
            }
            initialize_leapfrog_half_step(
                part, interpolate_electric(grid_, part.x), qm,
                cfg_.dt);
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
    const auto append_product =
        [&](std::size_t species_id,
            double position,
            Vec3 velocity) {
        auto& product_species =
            species_[species_id];
        auto& particles =
            product_species.particles();
        auto dead = std::find_if(
            particles.begin(), particles.end(),
            [](const Particle& particle) {
                return !particle.alive;
            });
        if (dead == particles.end()) {
            particles.emplace_back();
            dead = std::prev(particles.end());
        }
        *dead = {};
        dead->x = position;
        dead->v = velocity.x;
        dead->velocity_y = velocity.y;
        dead->velocity_z = velocity.z;
        dead->v_half = dead->v;
        dead->alive = true;
        initialize_leapfrog_half_step(
            *dead,
            interpolate_electric(grid_, position),
            product_species.charge() /
                product_species.mass(),
            cfg_.dt);
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
        auto& chunk_power =
            power_transfer_chunks_[species_id];
        std::fill(
            chunk_power.begin(), chunk_power.end(),
            SpeciesPower1D{});
        runtime_static_chunks(
            std::size_t{0}, particles.size(), cfg_.runtime,
            [&](std::size_t chunk, std::size_t begin,
                std::size_t end) {
                auto& loss = chunk_losses[chunk];
                auto& power = chunk_power[chunk];
                for (std::size_t particle_id = begin;
                     particle_id < end; ++particle_id) {
                    auto& p = particles[particle_id];
                    if (!p.alive) continue;
                    kick_leapfrog(
                        p, interpolate_electric(grid_, p.x),
                        qm, cfg_.dt);
                    drift_leapfrog(p, cfg_.dt);
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
                        if (p.x < 0.0) {
                            ++loss.absorbed_left;
                            loss.kinetic_energy_left +=
                                represented_energy;
                        } else {
                            ++loss.absorbed_right;
                            loss.kinetic_energy_right +=
                                represented_energy;
                        }
                        p.alive = false;
                    }
                }
            });
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
    }
    deposit_and_solve(time_ + cfg_.dt);
    for (std::size_t species_id = 0;
         species_id < species_.size(); ++species_id) {
        auto& sp = species_[species_id];
        const double qm = sp.charge() / sp.mass();
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
                        qm, cfg_.dt);
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
    apply_collisions();
    ++step_;
    time_ += cfg_.dt;
    accumulate_spatial_average();
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

void Simulation::accumulate_spatial_average() {
    const auto& average = cfg_.spatial_average;
    if (!average.enabled ||
        step_ < average.start_step ||
        step_ > average.end_step ||
        (step_ - average.start_step) %
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
    if (!spatial_phase_bins_.empty()) {
        const auto steps_per_cycle = static_cast<std::size_t>(std::llround(
            1.0 / (average.rf_frequency * cfg_.dt)));
        const auto samples_per_cycle = steps_per_cycle / average.interval;
        const auto sample_in_cycle =
            ((step_ - average.start_step) / average.interval) %
            samples_per_cycle;
        const auto phase_id =
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
    ++spatial_average_samples_;
    ++spatial_moment_samples_;
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
             << "  \"spatial_average_version\": 3,\n"
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
    out << kCheckpointMagicV10 << '\n';
    out << "dimension 1\n";
    out << "units " << to_string(cfg_.units.system) << ' '
        << cfg_.units.relative_permittivity << ' '
        << cfg_.units.permittivity() << "\n";
    out << "velocity_dimensions "
        << cfg_.velocity_dimensions << "\n";
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
    const bool checkpoint_v4_state =
        checkpoint_v4 || checkpoint_v5 ||
        checkpoint_v6 || checkpoint_v7 || checkpoint_v8 || checkpoint_v9 ||
        checkpoint_v10;
    if (!checkpoint_v1 && !checkpoint_v2 &&
        !checkpoint_v3 && !checkpoint_v4 &&
        !checkpoint_v5 && !checkpoint_v6 &&
        !checkpoint_v7 && !checkpoint_v8 && !checkpoint_v9 &&
        !checkpoint_v10) {
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
             !checkpoint_v9 && !checkpoint_v10) ||
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
               checkpoint_v9 || checkpoint_v10 ||
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
             !checkpoint_v10) ||
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
        if (checkpoint_v10) {
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
               checkpoint_v10 ||
               !mcc_models_.empty()) {
        throw std::runtime_error(
            "legacy checkpoint without MCC metadata cannot restart "
            "null-collision MCC");
    }
    const bool checkpoint_has_boundary_losses =
        checkpoint_v6 || checkpoint_v7 || checkpoint_v8 || checkpoint_v9 ||
        checkpoint_v10;
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
    const bool legacy_power_transfer_origin =
        !checkpoint_v7 && !checkpoint_v8 && !checkpoint_v9 && !checkpoint_v10;
    if (checkpoint_v7 || checkpoint_v8 || checkpoint_v9 || checkpoint_v10) {
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
        checkpoint_v7 || checkpoint_v8 || checkpoint_v9 || checkpoint_v10) {
        int enabled = 0;
        std::size_t interval = 0;
        std::size_t start_step = 0;
        std::size_t end_step = 0;
        double rf_frequency = 0.0;
        std::size_t rf_cycles = 0;
        std::size_t stored_samples = 0;
        std::size_t stored_species_count = 0;
        std::size_t stored_nx = 0;
        in >> enabled >> interval >> start_step >>
            end_step >> rf_frequency >> rf_cycles >>
            stored_samples >> stored_species_count >>
            stored_nx;
        const auto& configured = cfg_.spatial_average;
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
        if (checkpoint_v8 || checkpoint_v9 || checkpoint_v10) {
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
        if (checkpoint_v9 || checkpoint_v10) {
            std::size_t stored_phase_count = 0;
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
        if (reset) {
            spatial_moment_samples_ = 0;
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
    if (legacy_boundary_loss_origin) {
        boundary_loss_origin_step_ = step_;
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
    write_spatial_average();
    return summary;
}
}
