#pragma once
#include "pic/Config.hpp"
#include "pic/Grid.hpp"
#include <random>
#include <string>
#include <vector>

namespace pic {
class Species {
public:
    explicit Species(SpeciesConfig cfg);
    const SpeciesConfig& config() const { return cfg_; }
    const std::string& name() const { return cfg_.name; }
    double charge() const { return cfg_.charge; }
    double mass() const { return cfg_.mass; }
    double weight() const { return cfg_.weight; }
    std::vector<Particle>& particles() { return particles_; }
    const std::vector<Particle>& particles() const { return particles_; }
    void initialize(const Grid& grid, std::mt19937_64& rng);
    void deposit_charge(Grid& grid) const;
    double kinetic_energy() const;
    std::size_t live_count() const;
private:
    SpeciesConfig cfg_;
    std::vector<Particle> particles_;
};
}
