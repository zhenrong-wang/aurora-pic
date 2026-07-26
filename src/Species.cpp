#include "pic/Species.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace pic {
Species::Species(SpeciesConfig cfg) : cfg_(std::move(cfg)) {
    if (cfg_.mass <= 0.0) throw std::invalid_argument("species mass must be positive");
    if (cfg_.weight <= 0.0) throw std::invalid_argument("species weight must be positive");
    if (cfg_.particles == 0) throw std::invalid_argument("species must contain particles");
    if (cfg_.thermal_velocity < 0.0) throw std::invalid_argument("species thermal_velocity must be non-negative");
}

void Species::initialize(const Grid& grid, std::mt19937_64& rng) {
    particles_.assign(cfg_.particles, Particle{});
    const double xmin = cfg_.init_x_min;
    const double xmax = cfg_.init_x_max < 0.0 ? grid.length() : cfg_.init_x_max;
    std::uniform_real_distribution<double> ux(xmin, xmax);
    std::normal_distribution<double> nv(cfg_.drift_velocity, cfg_.thermal_velocity);
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
