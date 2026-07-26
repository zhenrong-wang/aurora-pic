#include "pic/Simulation.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace pic {
namespace {
void validate_runtime_config(const Config& cfg) {
    if (cfg.dt <= 0.0) throw std::invalid_argument("simulation dt must be positive");
    if (cfg.output_interval == 0) throw std::invalid_argument("output_interval must be positive");
    if (cfg.mode == RunMode::SteadyState) {
        if (cfg.max_steps == 0) throw std::invalid_argument("max_steps must be positive for steady-state mode");
        if (cfg.steady_window == 0) throw std::invalid_argument("steady_window must be positive");
        if (cfg.steady_tolerance <= 0.0) throw std::invalid_argument("steady_tolerance must be positive");
    }
    if (cfg.collisions.frequency < 0.0) throw std::invalid_argument("collision frequency must be non-negative");
    if (cfg.collisions.neutral_temperature_velocity < 0.0) {
        throw std::invalid_argument("neutral_temperature_velocity must be non-negative");
    }
}
}

Simulation::Simulation(Config cfg)
    : cfg_(std::move(cfg)), grid_(cfg_.nx, cfg_.length, cfg_.boundary), rng_(cfg_.seed) {
    validate_runtime_config(cfg_);
    for (const auto& sc : cfg_.species) species_.emplace_back(sc);
}

void Simulation::initialize() {
    time_ = 0.0;
    step_ = 0;
    for (auto& sp : species_) sp.initialize(grid_, rng_);
    grid_.clear_charge();
    for (const auto& sp : species_) sp.deposit_charge(grid_);
    solver_.solve(grid_, cfg_.phi_left, cfg_.phi_right);
}

void Simulation::apply_collisions(Species& sp) {
    if (!cfg_.collisions.enabled || cfg_.collisions.frequency <= 0.0) return;
    const double p = 1.0 - std::exp(-cfg_.collisions.frequency * cfg_.dt);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::normal_distribution<double> nv(0.0, cfg_.collisions.neutral_temperature_velocity);
    for (auto& part : sp.particles()) if (part.alive && u(rng_) < p) part.v = nv(rng_);
}

void Simulation::step() {
    for (auto& sp : species_) {
        const double qm = sp.charge() / sp.mass();
        for (auto& p : sp.particles()) {
            if (!p.alive) continue;
            const double E = interpolate_electric(grid_, p.x);
            p.v += qm * E * cfg_.dt;
            p.x += p.v * cfg_.dt;
            if (grid_.boundary() == Boundary::Periodic) {
                p.x = std::fmod(std::fmod(p.x, grid_.length()) + grid_.length(), grid_.length());
            } else if (p.x < 0.0 || p.x > grid_.length()) {
                p.alive = false;
            }
        }
        apply_collisions(sp);
    }
    grid_.clear_charge();
    for (const auto& sp : species_) sp.deposit_charge(grid_);
    solver_.solve(grid_, cfg_.phi_left, cfg_.phi_right);
    ++step_;
    time_ += cfg_.dt;
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
    initialize();
    Diagnostics diag(cfg_.output_dir);
    diag.write_header();
    auto s0 = diag.sample(0, 0.0, grid_, species_);
    diag.write_sample(s0);
    diag.write_fields(0, grid_);

    const std::size_t limit = cfg_.mode == RunMode::SteadyState ? cfg_.max_steps : cfg_.steps;
    RunSummary summary;
    summary.final_sample = s0;
    for (std::size_t n = 0; n < limit; ++n) {
        step();
        if (step_ % cfg_.output_interval == 0 || step_ == limit) {
            auto s = diag.sample(step_, time_, grid_, species_);
            diag.write_sample(s);
            diag.write_fields(step_, grid_);
            summary.final_sample = s;
            if (steady_converged(diag.history())) {
                summary.steady_state_reached = true;
                break;
            }
        }
    }
    summary.steps_completed = step_;
    summary.final_time = time_;
    return summary;
}
}
