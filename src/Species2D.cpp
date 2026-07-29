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
    if (!std::isfinite(cfg_.charge)) throw std::invalid_argument("2D species charge must be finite");
    if (!std::isfinite(cfg_.mass) || cfg_.mass <= 0.0) throw std::invalid_argument("2D species mass must be positive and finite");
    if (!std::isfinite(cfg_.weight) || cfg_.weight <= 0.0) throw std::invalid_argument("2D species weight must be positive and finite");
    if (cfg_.particles == 0) throw std::invalid_argument("2D species must contain particles");
    if (!std::isfinite(cfg_.drift_velocity_x) ||
        !std::isfinite(cfg_.drift_velocity_y) ||
        !std::isfinite(cfg_.drift_velocity_z)) {
        throw std::invalid_argument("2D species drift velocities must be finite");
    }
    validate_particle_initialization(
        cfg_.initialization, 3, cfg_.thermal_velocity, "2D species");
    validate_density_profile(
        cfg_.initialization, 2, cfg_.particles, "2D species");
}

void Species2D::initialize(const Mesh2D& mesh, std::mt19937_64& rng) {
    particles_.assign(cfg_.particles, Particle2D{});
    const double xmin = cfg_.init_x_min;
    const double xmax = cfg_.init_x_max < 0.0 ? mesh.length_x() : cfg_.init_x_max;
    const double ymin = cfg_.init_y_min;
    const double ymax = cfg_.init_y_max < 0.0 ? mesh.length_y() : cfg_.init_y_max;
    const double thermal_x = resolved_thermal_velocity(
        cfg_.initialization, 0, cfg_.thermal_velocity);
    const double thermal_y = resolved_thermal_velocity(
        cfg_.initialization, 1, cfg_.thermal_velocity);
    const double thermal_z = resolved_thermal_velocity(
        cfg_.initialization, 2, cfg_.thermal_velocity);

    if (cfg_.initialization.density_profile ==
            DensityProfileKind::Uniform &&
        cfg_.initialization.loading == ParticleLoading::Random) {
        std::uniform_real_distribution<double> ux(xmin, xmax);
        std::uniform_real_distribution<double> uy(ymin, ymax);
        std::normal_distribution<double> vx(
            cfg_.drift_velocity_x, thermal_x);
        std::normal_distribution<double> vy(
            cfg_.drift_velocity_y, thermal_y);
        std::normal_distribution<double> vz(
            cfg_.drift_velocity_z, thermal_z);

        for (auto& particle : particles_) {
            particle.position.x = ux(rng);
            particle.position.y = uy(rng);
            if (mesh.boundary_x() == Boundary::Periodic) {
                particle.position.x = wrap_periodic(particle.position.x, mesh.length_x());
            } else {
                particle.position.x = std::clamp(particle.position.x, 0.0, mesh.length_x());
            }
            if (mesh.boundary_y() == Boundary::Periodic) {
                particle.position.y = wrap_periodic(particle.position.y, mesh.length_y());
            } else {
                particle.position.y = std::clamp(particle.position.y, 0.0, mesh.length_y());
            }
            particle.velocity.x = vx(rng);
            particle.velocity.y = vy(rng);
            particle.velocity_z = vz(rng);
            particle.velocity_half = particle.velocity;
            particle.velocity_half_z = particle.velocity_z;
            particle.alive = true;
        }
        return;
    }

    const auto velocity_x = initialize_velocity_component(
        particles_.size(), cfg_.drift_velocity_x, thermal_x,
        cfg_.initialization.loading, rng);
    const auto velocity_y = initialize_velocity_component(
        particles_.size(), cfg_.drift_velocity_y, thermal_y,
        cfg_.initialization.loading, rng);
    const auto velocity_z = initialize_velocity_component(
        particles_.size(), cfg_.drift_velocity_z, thermal_z,
        cfg_.initialization.loading, rng);

    if (cfg_.initialization.density_profile !=
        DensityProfileKind::Uniform) {
        std::uniform_real_distribution<double> unit(0.0, 1.0);
        const std::array<double, 3> minimum{xmin, ymin, 0.0};
        const std::array<double, 3> maximum{xmax, ymax, 1.0};
        std::size_t accepted = 0;
        std::size_t attempts = 0;
        while (accepted < particles_.size()) {
            if (attempts >=
                cfg_.initialization.max_profile_sampling_attempts) {
                throw std::runtime_error(
                    "2D species density-profile sampling exceeded max_profile_sampling_attempts");
            }
            const std::size_t sequence = attempts++;
            const bool quiet =
                cfg_.initialization.loading ==
                ParticleLoading::QuietStart;
            const double x_coordinate = quiet
                ? quiet_sequence_coordinate(sequence, 0)
                : unit(rng);
            const double y_coordinate = quiet
                ? quiet_sequence_coordinate(sequence, 1)
                : unit(rng);
            const double threshold = quiet
                ? quiet_sequence_coordinate(sequence, 2)
                : unit(rng);
            const Vec2 position{
                xmin + (xmax - xmin) * x_coordinate,
                ymin + (ymax - ymin) * y_coordinate};
            if (threshold > density_profile_acceptance(
                    cfg_.initialization,
                    {position.x, position.y, 0.0},
                    minimum, maximum)) {
                continue;
            }
            auto& particle = particles_[accepted];
            particle.position = position;
            particle.velocity.x = velocity_x[accepted];
            particle.velocity.y = velocity_y[accepted];
            particle.velocity_z = velocity_z[accepted];
            particle.velocity_half = particle.velocity;
            particle.velocity_half_z = particle.velocity_z;
            particle.alive = true;
            ++accepted;
        }
        return;
    }

    for (std::size_t index = 0; index < particles_.size(); ++index) {
        auto& particle = particles_[index];
        particle.position.x = xmin + (xmax - xmin) *
            quiet_unit_coordinate(index, particles_.size(), 0);
        particle.position.y = ymin + (ymax - ymin) *
            quiet_unit_coordinate(index, particles_.size(), 1);
        if (mesh.boundary_x() == Boundary::Periodic) {
            particle.position.x = wrap_periodic(particle.position.x, mesh.length_x());
        } else {
            particle.position.x = std::clamp(particle.position.x, 0.0, mesh.length_x());
        }
        if (mesh.boundary_y() == Boundary::Periodic) {
            particle.position.y = wrap_periodic(particle.position.y, mesh.length_y());
        } else {
            particle.position.y = std::clamp(particle.position.y, 0.0, mesh.length_y());
        }
        particle.velocity.x = velocity_x[index];
        particle.velocity.y = velocity_y[index];
        particle.velocity_z = velocity_z[index];
        particle.velocity_half = particle.velocity;
        particle.velocity_half_z = particle.velocity_z;
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
        const double v2 =
            particle.velocity.x * particle.velocity.x +
            particle.velocity.y * particle.velocity.y +
            particle.velocity_z * particle.velocity_z;
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
