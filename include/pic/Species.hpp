#pragma once
#include "pic/Config.hpp"
#include "pic/Grid.hpp"
#include <random>
#include <string>
#include <vector>

namespace pic {
class Species {
public:
    explicit Species(
        SpeciesConfig cfg, std::size_t velocity_dimensions = 1);
    const SpeciesConfig& config() const { return cfg_; }
    const std::string& name() const { return cfg_.name; }
    double charge() const { return cfg_.charge; }
    double mass() const { return cfg_.mass; }
    double weight() const { return cfg_.weight; }
    std::vector<Particle>& particles() { return particles_; }
    const std::vector<Particle>& particles() const { return particles_; }
    void initialize(const Grid& grid, std::mt19937_64& rng);
    void deposit_charge(Grid& grid) const;
    void deposit_number_density(
        const Grid& grid, std::vector<double>& density) const;
    void deposit_kinetic_energy_density(
        const Grid& grid, std::vector<double>& energy_density) const;
    double kinetic_energy() const;
    std::size_t live_count() const;
    std::size_t velocity_dimensions() const {
        return velocity_dimensions_;
    }
private:
    SpeciesConfig cfg_;
    std::size_t velocity_dimensions_{1};
    std::vector<Particle> particles_;
};
}
