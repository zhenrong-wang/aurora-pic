#include "pic/Simulation.hpp"
#include "pic/Pusher.hpp"
#include "pic/Runtime.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace pic {
namespace {
constexpr const char* kCheckpointMagic = "AuroraPIC-checkpoint-v1";

void validate_runtime_config(const Config& cfg) {
    if (cfg.dt <= 0.0) throw std::invalid_argument("simulation dt must be positive");
    if (cfg.output_interval == 0) throw std::invalid_argument("output_interval must be positive");
    if (cfg.checkpoint_output && cfg.checkpoint_interval == 0) {
        throw std::invalid_argument("checkpoint_interval must be positive when checkpoint_output is enabled");
    }
    if (cfg.mode == RunMode::SteadyState) {
        if (cfg.max_steps == 0) throw std::invalid_argument("max_steps must be positive for steady-state mode");
        if (cfg.steady_window == 0) throw std::invalid_argument("steady_window must be positive");
        if (cfg.steady_tolerance <= 0.0) throw std::invalid_argument("steady_tolerance must be positive");
    }
    validate_runtime_policy(cfg.runtime);
    if (cfg.collisions.frequency < 0.0) throw std::invalid_argument("collision frequency must be non-negative");
    if (cfg.collisions.neutral_temperature_velocity < 0.0) {
        throw std::invalid_argument("neutral_temperature_velocity must be non-negative");
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
} // namespace

Simulation::Simulation(Config cfg)
    : cfg_(std::move(cfg)), grid_(cfg_.nx, cfg_.length, cfg_.boundary), rng_(cfg_.seed) {
    if (cfg_.checkpoint_output && cfg_.checkpoint_interval == 0) cfg_.checkpoint_interval = cfg_.output_interval;
    validate_runtime_config(cfg_);
    for (const auto& sc : cfg_.species) species_.emplace_back(sc);
}

void Simulation::deposit_and_solve() {
    grid_.clear_charge();
    for (const auto& sp : species_) sp.deposit_charge(grid_);
    solver_.solve(grid_, cfg_.phi_left, cfg_.phi_right);
}

void Simulation::initialize() {
    time_ = 0.0;
    step_ = 0;
    for (auto& sp : species_) sp.initialize(grid_, rng_);
    deposit_and_solve();
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

void Simulation::apply_collisions(Species& sp) {
    if (!cfg_.collisions.enabled || cfg_.collisions.frequency <= 0.0) return;
    const double p = 1.0 - std::exp(-cfg_.collisions.frequency * cfg_.dt);
    const double qm = sp.charge() / sp.mass();
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::normal_distribution<double> nv(0.0, cfg_.collisions.neutral_temperature_velocity);
    for (auto& part : sp.particles()) {
        if (!part.alive || u(rng_) >= p) continue;
        part.v = nv(rng_);
        initialize_leapfrog_half_step(part, interpolate_electric(grid_, part.x), qm, cfg_.dt);
    }
}

void Simulation::step() {
    if (!initialized_) initialize();
    for (auto& sp : species_) {
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
    deposit_and_solve();
    for (auto& sp : species_) {
        const double qm = sp.charge() / sp.mass();
        auto& particles = sp.particles();
        runtime_parallel_for(std::size_t{0}, particles.size(), cfg_.runtime, [&](std::size_t particle_id) {
            auto& p = particles[particle_id];
            if (p.alive) synchronize_leapfrog(p, interpolate_electric(grid_, p.x), qm, cfg_.dt);
        });
        apply_collisions(sp);
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
        s.field_energy += 0.5 * EPS0 * grid_.electric()[i] * grid_.electric()[i] * volume;
        s.charge_l1 += std::abs(grid_.rho()[i]) * volume;
    }
    s.total_energy = s.kinetic_energy + s.field_energy;
    return s;
}

void Simulation::save_checkpoint(const std::filesystem::path& path) const {
    ensure_parent_directory(path);
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot open checkpoint for writing: " + path.string());
    out << std::setprecision(17);
    out << kCheckpointMagic << '\n';
    out << "dimension 1\n";
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
    if (magic != kCheckpointMagic) throw std::runtime_error("invalid checkpoint magic in: " + path.string());

    std::string key;
    unsigned dimension = 0;
    in >> key >> dimension;
    if (key != "dimension" || dimension != 1) throw std::runtime_error("checkpoint dimension does not match 1D simulation");
    in >> key >> step_;
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
    deposit_and_solve();
    initialized_ = true;
}

bool Simulation::steady_converged(const std::vector<DiagnosticSample>& history) const {
    if (cfg_.mode != RunMode::SteadyState || history.size() < 2 * cfg_.steady_window) return false;

    double current_window_mean = 0.0;
    double previous_window_mean = 0.0;
    for (std::size_t i = history.size() - cfg_.steady_window; i < history.size(); ++i) {
        current_window_mean += history[i].total_energy;
    }
    for (std::size_t i = history.size() - 2 * cfg_.steady_window; i < history.size() - cfg_.steady_window; ++i) {
        previous_window_mean += history[i].total_energy;
    }
    current_window_mean /= static_cast<double>(cfg_.steady_window);
    previous_window_mean /= static_cast<double>(cfg_.steady_window);

    const double rel = std::abs(current_window_mean - previous_window_mean) /
                       std::max(1e-30, std::abs(previous_window_mean));
    return rel < cfg_.steady_tolerance;
}

RunSummary Simulation::run() {
    if (!cfg_.restart_path.empty()) load_checkpoint(cfg_.restart_path);
    else initialize();

    Diagnostics diag(cfg_.output_dir);
    diag.write_header();
    auto s0 = diag.sample(step_, time_, grid_, species_);
    diag.write_sample(s0);
    diag.write_fields(step_, grid_);
    if (cfg_.checkpoint_output) save_checkpoint(checkpoint_path_for_step(cfg_, step_));

    const std::size_t limit = cfg_.mode == RunMode::SteadyState ? cfg_.max_steps : cfg_.steps;
    RunSummary summary;
    summary.final_sample = s0;
    while (step_ < limit) {
        step();
        const bool at_output = step_ % cfg_.output_interval == 0 || step_ == limit;
        if (at_output) {
            auto s = diag.sample(step_, time_, grid_, species_);
            diag.write_sample(s);
            diag.write_fields(step_, grid_);
            summary.final_sample = s;
            if (steady_converged(diag.history())) {
                summary.steady_state_reached = true;
                break;
            }
        }
        if (cfg_.checkpoint_output && (step_ % cfg_.checkpoint_interval == 0 || step_ == limit)) {
            save_checkpoint(checkpoint_path_for_step(cfg_, step_));
        }
    }
    summary.steps_completed = step_;
    summary.final_time = time_;
    if (summary.final_sample.step != step_) summary.final_sample = sample();
    return summary;
}
}
