#include "pic/Simulation.hpp"
#include "pic/Convergence.hpp"
#include "pic/ParticleState.hpp"
#include "pic/Pusher.hpp"
#include "pic/Runtime.hpp"
#include "pic/Units.hpp"
#include <algorithm>
#include <bit>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace pic {
namespace {
constexpr const char* kCheckpointMagicV1 = "AuroraPIC-checkpoint-v1";
constexpr const char* kCheckpointMagicV2 = "AuroraPIC-checkpoint-v2";
constexpr const char* kCheckpointMagicV3 = "AuroraPIC-checkpoint-v3";

void validate_runtime_config(const Config& cfg) {
    if (!std::isfinite(cfg.dt) || cfg.dt <= 0.0) throw std::invalid_argument("simulation dt must be positive and finite");
    if (cfg.output_interval == 0) throw std::invalid_argument("output_interval must be positive");
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
    const CollisionStepStatistics& source) {
    destination.candidates += source.candidates;
    destination.null_collisions += source.null_collisions;
    if (destination.channel_collisions.size() !=
        source.channel_collisions.size()) {
        throw std::logic_error("collision diagnostic channel mismatch");
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
        output << ",collisions_" << name;
    }
    output << ",cumulative_candidates,cumulative_null_collisions";
    for (const auto& name : diagnostics.channel_names) {
        output << ",cumulative_collisions_" << name;
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

std::uint64_t collision_signature(
    const Config& config,
    const NullCollisionModel* model) {
    if (!config.collisions.enabled) return 0;
    if (model) return model->signature();
    constexpr std::uint64_t offset = 1469598103934665603ULL;
    constexpr std::uint64_t prime = 1099511628211ULL;
    std::uint64_t hash = offset;
    const auto append = [&](std::uint64_t value) {
        for (unsigned byte = 0; byte < 8; ++byte) {
            hash ^= static_cast<unsigned char>(
                value >> (byte * 8));
            hash *= prime;
        }
    };
    append(std::bit_cast<std::uint64_t>(
        config.collisions.frequency));
    append(std::bit_cast<std::uint64_t>(
        config.collisions.neutral_temperature_velocity));
    return hash;
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
    for (const auto& sc : cfg_.species) species_.emplace_back(sc);
    if (cfg_.collisions.enabled) {
        if (cfg_.collisions.model ==
            CollisionModelKind::NullCollision) {
            const auto target = std::find_if(
                species_.begin(), species_.end(),
                [&](const Species& species) {
                    return species.name() == cfg_.collisions.species;
                });
            if (target == species_.end()) {
                throw std::invalid_argument(
                    "MCC target species does not exist: " +
                    cfg_.collisions.species);
            }
            mcc_species_id_ =
                static_cast<std::size_t>(target - species_.begin());
            mcc_model_ = std::make_unique<NullCollisionModel>(
                cfg_.collisions, target->mass());
            collision_totals_.channel_names =
                mcc_model_->channel_names();
        } else {
            collision_totals_.channel_names = {"bgk"};
        }
        collision_totals_.channel_collisions.assign(
            collision_totals_.channel_names.size(), 0);
        collision_interval_.channel_names =
            collision_totals_.channel_names;
        collision_interval_.channel_collisions.assign(
            collision_totals_.channel_names.size(), 0);
    }
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

void Simulation::apply_collisions(
    Species& sp, std::size_t species_id) {
    if (!cfg_.collisions.enabled) return;
    if (cfg_.collisions.model ==
        CollisionModelKind::NullCollision) {
        if (species_id != mcc_species_id_) return;
        const double qm = sp.charge() / sp.mass();
        for (auto& part : sp.particles()) {
            if (!part.alive) continue;
            const auto statistics =
                mcc_model_->collide(part.v, cfg_.dt, rng_);
            add_collision_statistics(
                collision_totals_, statistics);
            add_collision_statistics(
                collision_interval_, statistics);
            if (statistics.primary_removal_channel) {
                part.alive = false;
                continue;
            }
            initialize_leapfrog_half_step(
                part, interpolate_electric(grid_, part.x), qm,
                cfg_.dt);
        }
        return;
    }
    if (cfg_.collisions.frequency <= 0.0) return;
    const double p = 1.0 - std::exp(-cfg_.collisions.frequency * cfg_.dt);
    const double qm = sp.charge() / sp.mass();
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::normal_distribution<double> nv(0.0, cfg_.collisions.neutral_temperature_velocity);
    for (auto& part : sp.particles()) {
        if (!part.alive || u(rng_) >= p) continue;
        part.v = nv(rng_);
        ++collision_totals_.candidates;
        ++collision_totals_.channel_collisions[0];
        ++collision_interval_.candidates;
        ++collision_interval_.channel_collisions[0];
        initialize_leapfrog_half_step(part, interpolate_electric(grid_, part.x), qm, cfg_.dt);
    }
}

void Simulation::step() {
    if (!initialized_) initialize();
    for (std::size_t species_id = 0;
         species_id < species_.size(); ++species_id) {
        auto& sp = species_[species_id];
        const double qm = sp.charge() / sp.mass();
        auto& particles = sp.particles();
        runtime_parallel_for(std::size_t{0}, particles.size(), cfg_.runtime, [&](std::size_t particle_id) {
            auto& p = particles[particle_id];
            if (!p.alive) return;
            kick_leapfrog(p, interpolate_electric(grid_, p.x), qm, cfg_.dt);
            drift_leapfrog(p, cfg_.dt);
            if (grid_.boundary() == Boundary::Periodic) {
                p.x = std::fmod(std::fmod(p.x, grid_.length()) + grid_.length(), grid_.length());
            } else if (p.x < 0.0 || p.x > grid_.length()) {
                p.alive = false;
            }
        });
    }
    deposit_and_solve(time_ + cfg_.dt);
    for (std::size_t species_id = 0;
         species_id < species_.size(); ++species_id) {
        auto& sp = species_[species_id];
        const double qm = sp.charge() / sp.mass();
        auto& particles = sp.particles();
        runtime_parallel_for(std::size_t{0}, particles.size(), cfg_.runtime, [&](std::size_t particle_id) {
            auto& p = particles[particle_id];
            if (p.alive) synchronize_leapfrog(p, interpolate_electric(grid_, p.x), qm, cfg_.dt);
        });
        apply_collisions(sp, species_id);
    }
    ++step_;
    time_ += cfg_.dt;
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

void Simulation::save_checkpoint(const std::filesystem::path& path) const {
    ensure_parent_directory(path);
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot open checkpoint for writing: " + path.string());
    out << std::setprecision(17);
    out << kCheckpointMagicV3 << '\n';
    out << "dimension 1\n";
    out << "units " << to_string(cfg_.units.system) << ' '
        << cfg_.units.relative_permittivity << ' '
        << cfg_.units.permittivity() << "\n";
    const std::uint64_t configured_collision_signature =
        pic::collision_signature(cfg_, mcc_model_.get());
    out << "collision_model "
        << to_string(cfg_.collisions.model) << ' '
        << (cfg_.collisions.enabled ? 1 : 0) << ' '
        << configured_collision_signature << "\n";
    out << "collision_totals " << collision_totals_.candidates
        << ' ' << collision_totals_.null_collisions << ' '
        << collision_totals_.channel_collisions.size();
    for (const auto count : collision_totals_.channel_collisions) {
        out << ' ' << count;
    }
    out << "\n";
    out << "step " << step_ << "\n";
    out << "time " << time_ << "\n";
    out << "species_count " << species_.size() << "\n";
    out << "rng " << rng_ << "\n";
    for (std::size_t species_id = 0; species_id < species_.size(); ++species_id) {
        const auto& sp = species_[species_id];
        out << "species " << species_id << ' ' << sp.name() << ' ' << sp.particles().size() << "\n";
        for (const auto& p : sp.particles()) {
            out << p.x << ' ' << p.v << ' ' << p.v_half << ' ' << (p.alive ? 1 : 0) << "\n";
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
    if (!checkpoint_v1 && !checkpoint_v2 && !checkpoint_v3) {
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
        if ((!checkpoint_v2 && !checkpoint_v3) ||
            unit_system != to_string(cfg_.units.system) ||
            relative_permittivity != cfg_.units.relative_permittivity ||
            permittivity != cfg_.units.permittivity()) {
            throw std::runtime_error(
                "checkpoint unit system does not match 1D config");
        }
        in >> key;
    } else if (checkpoint_v2 || checkpoint_v3 ||
               cfg_.units.system != UnitSystem::Normalized ||
               cfg_.units.relative_permittivity != 1.0) {
        throw std::runtime_error(
            "legacy checkpoint without unit metadata requires normalized units");
    }
    if (key == "collision_model") {
        std::string model;
        int enabled = 0;
        std::uint64_t signature = 0;
        in >> model >> enabled >> signature;
        const std::uint64_t expected_signature =
            pic::collision_signature(cfg_, mcc_model_.get());
        if (!checkpoint_v3 ||
            model != to_string(cfg_.collisions.model) ||
            enabled != (cfg_.collisions.enabled ? 1 : 0) ||
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
    } else if (checkpoint_v3 ||
               (cfg_.collisions.enabled &&
                cfg_.collisions.model ==
                    CollisionModelKind::NullCollision)) {
        throw std::runtime_error(
            "legacy checkpoint without MCC metadata cannot restart "
            "null-collision MCC");
    }
    in >> step_;
    if (key != "step") throw std::runtime_error("checkpoint missing step");
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
        auto& particles = species_[species_id].particles();
        particles.resize(particle_count);
        for (auto& p : particles) {
            int alive = 0;
            in >> p.x >> p.v >> p.v_half >> alive;
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
            initialization_moments, 1);
    write_initialization_acceptance_report(
        std::filesystem::path(cfg_.output_dir) /
            "initialization_acceptance.csv",
        initialization_acceptance);
    enforce_initialization_acceptance(
        initialization_acceptance);
    Diagnostics diag(cfg_.output_dir, cfg_.units.permittivity());
    diag.write_header();
    auto s0 = diag.sample(step_, time_, grid_, species_);
    diag.write_sample(s0);
    diag.write_fields(step_, grid_);
    std::ofstream collision_output;
    if (cfg_.collisions.enabled) {
        collision_output.open(
            std::filesystem::path(cfg_.output_dir) /
            "collisions.csv");
        if (!collision_output) {
            throw std::runtime_error(
                "cannot open collision diagnostics output");
        }
        write_collision_header(
            collision_output, collision_totals_);
        write_collision_sample(
            collision_output, step_, time_,
            collision_interval_, collision_totals_);
        clear_collision_counts(collision_interval_);
    }
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
            if (cfg_.collisions.enabled) {
                write_collision_sample(
                    collision_output, step_, time_,
                    collision_interval_, collision_totals_);
                clear_collision_counts(collision_interval_);
            }
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
    return summary;
}
}
