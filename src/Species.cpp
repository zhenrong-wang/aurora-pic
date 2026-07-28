#include "pic/Species.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace pic {
Species::Species(SpeciesConfig cfg) : cfg_(std::move(cfg)) {
    if (!std::isfinite(cfg_.charge)) throw std::invalid_argument("species charge must be finite");
    if (!std::isfinite(cfg_.mass) || cfg_.mass <= 0.0) throw std::invalid_argument("species mass must be positive and finite");
    if (!std::isfinite(cfg_.weight) || cfg_.weight <= 0.0) throw std::invalid_argument("species weight must be positive and finite");
    if (cfg_.particles == 0) throw std::invalid_argument("species must contain particles");
    if (!std::isfinite(cfg_.drift_velocity)) throw std::invalid_argument("species drift_velocity must be finite");
    validate_particle_initialization(
        cfg_.initialization, 1, cfg_.thermal_velocity, "species");
    validate_density_profile(
        cfg_.initialization, 1, cfg_.particles, "species");
}

void Species::initialize(const Grid& grid, std::mt19937_64& rng) {
    particles_.assign(cfg_.particles, Particle{});
    const double xmin = cfg_.init_x_min;
    const double xmax = cfg_.init_x_max < 0.0 ? grid.length() : cfg_.init_x_max;
    const double thermal_velocity = resolved_thermal_velocity(
        cfg_.initialization, 0, cfg_.thermal_velocity);
    if (cfg_.initialization.density_profile ==
            DensityProfileKind::Uniform &&
        cfg_.initialization.loading == ParticleLoading::Random) {
        std::uniform_real_distribution<double> ux(xmin, xmax);
        std::normal_distribution<double> nv(
            cfg_.drift_velocity, thermal_velocity);
        for (auto& p : particles_) {
            p.x = ux(rng);
            if (grid.boundary() == Boundary::Periodic) {
                p.x = std::fmod(std::fmod(p.x, grid.length()) + grid.length(), grid.length());
            } else {
                p.x = std::clamp(p.x, 0.0, grid.length());
            }
            p.v = nv(rng);
            p.v_half = p.v;
            p.alive = true;
        }
        return;
    }

    const auto velocities = initialize_velocity_component(
        particles_.size(), cfg_.drift_velocity, thermal_velocity,
        cfg_.initialization.loading, rng);
    if (cfg_.initialization.density_profile !=
        DensityProfileKind::Uniform) {
        std::uniform_real_distribution<double> unit(0.0, 1.0);
        const std::array<double, 3> minimum{xmin, 0.0, 0.0};
        const std::array<double, 3> maximum{xmax, 1.0, 1.0};
        std::size_t accepted = 0;
        std::size_t attempts = 0;
        while (accepted < particles_.size()) {
            if (attempts >=
                cfg_.initialization.max_profile_sampling_attempts) {
                throw std::runtime_error(
                    "species density-profile sampling exceeded max_profile_sampling_attempts");
            }
            const std::size_t sequence = attempts++;
            const bool quiet =
                cfg_.initialization.loading ==
                ParticleLoading::QuietStart;
            const double coordinate = quiet
                ? quiet_sequence_coordinate(sequence, 0)
                : unit(rng);
            const double threshold = quiet
                ? quiet_sequence_coordinate(sequence, 1)
                : unit(rng);
            const double x = xmin + (xmax - xmin) * coordinate;
            if (threshold > density_profile_acceptance(
                    cfg_.initialization, {x, 0.0, 0.0},
                    minimum, maximum)) {
                continue;
            }
            auto& particle = particles_[accepted];
            particle.x = x;
            particle.v = velocities[accepted];
            particle.v_half = particle.v;
            particle.alive = true;
            ++accepted;
        }
        return;
    }

    for (std::size_t index = 0; index < particles_.size(); ++index) {
        auto& p = particles_[index];
        p.x = xmin + (xmax - xmin) *
            quiet_unit_coordinate(index, particles_.size(), 0);
        if (grid.boundary() == Boundary::Periodic) {
            p.x = std::fmod(std::fmod(p.x, grid.length()) + grid.length(), grid.length());
        } else {
            p.x = std::clamp(p.x, 0.0, grid.length());
        }
        p.v = velocities[index];
        p.v_half = p.v;
        p.alive = true;
    }
}

void Species::deposit_charge(Grid& grid) const {
    auto& rho = grid.rho();
    const double dx = grid.dx();
    const double qwdx = cfg_.charge * cfg_.weight / dx;
    for (const auto& p : particles_) {
        if (!p.alive) continue;
        if (grid.boundary() == Boundary::Periodic) {
            double xp = std::fmod(std::fmod(p.x, grid.length()) + grid.length(), grid.length());
            double g = xp / dx;
            auto i = static_cast<std::size_t>(std::floor(g));
            double f = g - static_cast<double>(i);
            std::size_t i0 = i % grid.nx();
            std::size_t i1 = (i + 1) % grid.nx();
            rho[i0] += qwdx * (1.0 - f);
            rho[i1] += qwdx * f;
        } else {
            double xp = std::clamp(p.x, 0.0, grid.length());
            double g = xp / dx;
            auto i = static_cast<std::size_t>(std::min<double>(std::floor(g), grid.nx() - 2));
            double f = g - static_cast<double>(i);
            rho[i] += cfg_.charge * cfg_.weight * (1.0 - f) / grid.node_volume(i);
            rho[i + 1] += cfg_.charge * cfg_.weight * f / grid.node_volume(i + 1);
        }
    }
}
double Species::kinetic_energy() const {
    double e = 0.0;
    for (const auto& p : particles_) if (p.alive) e += 0.5 * cfg_.mass * cfg_.weight * p.v * p.v;
    return e;
}

std::size_t Species::live_count() const {
    return static_cast<std::size_t>(std::count_if(particles_.begin(), particles_.end(), [](const Particle& p){ return p.alive; }));
}
}
