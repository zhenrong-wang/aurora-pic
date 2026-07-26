#pragma once
#include "pic/Config.hpp"
#include "pic/Diagnostics.hpp"
#include "pic/FieldSolver.hpp"
#include "pic/Species.hpp"
#include <random>

namespace pic {
struct RunSummary {
    std::size_t steps_completed{0};
    double final_time{0.0};
    bool steady_state_reached{false};
    DiagnosticSample final_sample{};
};

class Simulation {
public:
    explicit Simulation(Config cfg);
    RunSummary run();
    void initialize();
    void step();
    const Grid& grid() const { return grid_; }
    const std::vector<Species>& species() const { return species_; }
private:
    bool steady_converged(const std::vector<DiagnosticSample>& history) const;
    void apply_collisions(Species& sp);
    Config cfg_;
    Grid grid_;
    FieldSolver solver_;
    std::vector<Species> species_;
    std::mt19937_64 rng_;
    double time_{0.0};
    std::size_t step_{0};
};
}
