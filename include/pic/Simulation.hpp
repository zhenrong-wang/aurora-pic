#pragma once
#include "pic/Config.hpp"
#include "pic/Diagnostics.hpp"
#include "pic/FieldSolver.hpp"
#include "pic/Species.hpp"
#include <filesystem>
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
    void save_checkpoint(const std::filesystem::path& path) const;
    void load_checkpoint(const std::filesystem::path& path);
    DiagnosticSample sample() const;
    const Grid& grid() const { return grid_; }
    const std::vector<Species>& species() const { return species_; }
    double time() const { return time_; }
    std::size_t step_count() const { return step_; }
private:
    bool steady_converged(const std::vector<DiagnosticSample>& history) const;
    void apply_collisions(Species& sp);
    void deposit_and_solve();
    Config cfg_;
    Grid grid_;
    FieldSolver solver_;
    std::vector<Species> species_;
    std::mt19937_64 rng_;
    double time_{0.0};
    std::size_t step_{0};
    bool initialized_{false};
};
}
