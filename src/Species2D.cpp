#include "pic/Species2D.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
namespace pic {
namespace {
double wrap_periodic(double value, double length) {
    return std::fmod(std::fmod(value, length) + length, length);
}
}

Species2D::Species2D(Species2DConfig cfg) : cfg_(std::move(cfg)) {
    if (cfg_.mass <= 0.0) throw std::invalid_argument("2D species mass must be positive");
    if (cfg_.weight <= 0.0) throw std::invalid_argument("2D species weight must be positive");
    if (cfg_.particles == 0) throw std::invalid_argument("2D species must contain particles");
    if (cfg_.thermal_velocity < 0.0) throw std::invalid_argument("2D species thermal_velocity must be non-negative");
}

void Species2D::initialize(const Mesh2D& mesh, std::mt19937_64& rng) {
    particles_.assign(cfg_.particles, Particle2D{});
    const double xmin = cfg_.init_x_min;
    const double xmax = cfg_.init_x_max < 0.0 ? mesh.length_x() : cfg_.init_x_max;
    const double ymin = cfg_.init_y_min;
    const double ymax = cfg_.init_y_max < 0.0 ? mesh.length_y() : cfg_.init_y_max;
    std::uniform_real_distribution<double> ux(xmin, xmax);
    std::uniform_real_distribution<double> uy(ymin, ymax);
    std::normal_distribution<double> vx(cfg_.drift_velocity_x, cfg_.thermal_velocity);
    std::normal_distribution<double> vy(cfg_.drift_velocity_y, cfg_.thermal_velocity);

    for (auto& particle : particles_) {
        particle.position.x = ux(rng);
        particle.position.y = uy(rng);
        if (mesh.boundary() == Boundary::Periodic) {
            particle.position.x = wrap_periodic(particle.position.x, mesh.length_x());
            particle.position.y = wrap_periodic(particle.position.y, mesh.length_y());
        } else {
            particle.position.x = std::clamp(particle.position.x, 0.0, mesh.length_x());
            particle.position.y = std::clamp(particle.position.y, 0.0, mesh.length_y());
        }
        particle.velocity.x = vx(rng);
        particle.velocity.y = vy(rng);
        particle.alive = true;
    }
}

void Species2D::deposit_charge(Mesh2D& mesh) const {
    deposit_charge_cic(mesh, particles_, cfg_.charge, cfg_.weight);
}

double Species2D::kinetic_energy() const {
    double energy = 0.0;
    for (const auto& particle : particles_) {
        if (!particle.alive) continue;
        const double v2 = particle.velocity.x * particle.velocity.x + particle.velocity.y * particle.velocity.y;
        energy += 0.5 * cfg_.mass * cfg_.weight * v2;
    }
    return energy;
}

std::size_t Species2D::live_count() const {
    return static_cast<std::size_t>(std::count_if(particles_.begin(), particles_.end(), [](const Particle2D& particle) {
        return particle.alive;
    }));
}
}
