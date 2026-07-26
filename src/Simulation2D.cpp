#include "pic/Simulation2D.hpp"
#include "pic/Pusher.hpp"
#include "pic/VTKWriter.hpp"
#include <cmath>
#include <stdexcept>
#include <utility>

namespace pic {
namespace {
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
}

Simulation2D::Simulation2D(Simulation2DConfig cfg)
    : cfg_(std::move(cfg)), mesh_(cfg_.nx, cfg_.ny, cfg_.length_x, cfg_.length_y, cfg_.boundary, cfg_.boundary_config), rng_(cfg_.seed) {
    if (cfg_.dt <= 0.0) throw std::invalid_argument("2D simulation dt must be positive");
    if (cfg_.output_interval == 0) throw std::invalid_argument("2D output_interval must be positive");
    if (cfg_.particle_output_stride == 0) throw std::invalid_argument("2D particle_output_stride must be positive");
    cfg_.particle_boundary_config.left = resolve_particle_boundary(cfg_.particle_boundary_config.left, cfg_.boundary);
    cfg_.particle_boundary_config.right = resolve_particle_boundary(cfg_.particle_boundary_config.right, cfg_.boundary);
    cfg_.particle_boundary_config.bottom = resolve_particle_boundary(cfg_.particle_boundary_config.bottom, cfg_.boundary);
    cfg_.particle_boundary_config.top = resolve_particle_boundary(cfg_.particle_boundary_config.top, cfg_.boundary);
    for (const auto& sc : cfg_.species) species_.emplace_back(sc);
    if (species_.empty()) species_.emplace_back(Species2DConfig{});
}

void Simulation2D::initialize() {
    time_ = 0.0;
    step_ = 0;
    boundary_losses_ = {};
    for (auto& sp : species_) sp.initialize(mesh_, rng_);
    deposit_and_solve();
    for (auto& sp : species_) {
        const double qm = sp.charge() / sp.mass();
        for (auto& particle : sp.particles()) {
            if (particle.alive) initialize_leapfrog_half_step(particle, interpolate_electric(mesh_, particle.position), qm, cfg_.dt);
        }
    }
    initialized_ = true;
}

void Simulation2D::deposit_and_solve() {
    mesh_.clear_charge();
    for (const auto& sp : species_) sp.deposit_charge(mesh_);
    solver_.solve(mesh_);
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
    for (auto& sp : species_) {
        const double qm = sp.charge() / sp.mass();
        for (auto& particle : sp.particles()) {
            if (!particle.alive) continue;
            kick_leapfrog(particle, interpolate_electric(mesh_, particle.position), qm, cfg_.dt);
            drift_leapfrog(particle, cfg_.dt);
            apply_particle_boundaries(particle);
        }
    }

    deposit_and_solve();
    for (auto& sp : species_) {
        const double qm = sp.charge() / sp.mass();
        for (auto& particle : sp.particles()) {
            if (particle.alive) synchronize_leapfrog(particle, interpolate_electric(mesh_, particle.position), qm, cfg_.dt);
        }
    }
    ++step_;
    time_ += cfg_.dt;
}
RunSummary2D Simulation2D::run() {
    initialize();
    Diagnostics2D diag(cfg_.output_dir, species_);
    diag.write_header();
    auto s0 = diag.sample(step_, time_, mesh_, species_, boundary_losses_);
    diag.write_sample(s0);
    if (cfg_.vtk_output) write_legacy_vtk(mesh_, cfg_.output_dir / "fields_0.vtk");
    if (cfg_.particle_output) {
        diag.write_particle_sample(step_, species_, cfg_.particle_output_stride, cfg_.particle_sample_count);
    }

    const std::size_t particle_interval = cfg_.particle_output_interval == 0 ? cfg_.output_interval : cfg_.particle_output_interval;
    for (std::size_t n = 0; n < cfg_.steps; ++n) {
        step();
        if (step_ % cfg_.output_interval == 0 || step_ == cfg_.steps) {
            auto s = diag.sample(step_, time_, mesh_, species_, boundary_losses_);
            diag.write_sample(s);
            if (cfg_.vtk_output) {
                write_legacy_vtk(mesh_, cfg_.output_dir / ("fields_" + std::to_string(step_) + ".vtk"));
            }
        }
        if (cfg_.particle_output && (step_ % particle_interval == 0 || step_ == cfg_.steps)) {
            diag.write_particle_sample(step_, species_, cfg_.particle_output_stride, cfg_.particle_sample_count);
        }
    }
    RunSummary2D summary;
    summary.steps_completed = step_;
    summary.final_time = time_;
    summary.final_sample = diag.history().empty() ? sample() : diag.history().back();
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
            const double area = mesh_.node_area(i, j);
            s.field_energy += 0.5 * EPS0 * e2 * area;
            s.charge_l1 += std::abs(mesh_.rho()[idx]) * area;
        }
    }
    s.total_energy = s.kinetic_energy + s.field_energy;
    return s;
}
}
