#pragma once
#include "pic/Initialization.hpp"
#include "pic/Mesh3D.hpp"
#include <random>
#include <string>
#include <vector>

namespace pic {
struct Species3DConfig {
    std::string name{"electrons"};
    double charge{-1.0};
    double mass{1.0};
    double weight{1.0};
    std::size_t particles{1000};
    double drift_velocity_x{0.0};
    double drift_velocity_y{0.0};
    double drift_velocity_z{0.0};
    double thermal_velocity{0.1};
    double init_x_min{0.0};
    double init_x_max{-1.0}; // negative means full x domain
    double init_y_min{0.0};
    double init_y_max{-1.0}; // negative means full y domain
    double init_z_min{0.0};
    double init_z_max{-1.0}; // negative means full z domain
    ParticleInitializationConfig initialization{};
};

class Species3D {
public:
    explicit Species3D(Species3DConfig cfg);
    const Species3DConfig& config() const { return cfg_; }
    const std::string& name() const { return cfg_.name; }
    double charge() const { return cfg_.charge; }
    double mass() const { return cfg_.mass; }
    double weight() const { return cfg_.weight; }
    std::vector<Particle3D>& particles() { return particles_; }
    const std::vector<Particle3D>& particles() const { return particles_; }
    void initialize(const Mesh3D& mesh, std::mt19937_64& rng);
    void deposit_charge(Mesh3D& mesh) const;
    double kinetic_energy() const;
    std::size_t live_count() const;
private:
    Species3DConfig cfg_;
    std::vector<Particle3D> particles_;
};
}
