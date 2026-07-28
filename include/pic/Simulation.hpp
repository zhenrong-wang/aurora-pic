#pragma once
#include "pic/Collision.hpp"
#include "pic/Config.hpp"
#include "pic/Diagnostics.hpp"
#include "pic/FieldSolver.hpp"
#include "pic/Species.hpp"
#include <filesystem>
#include <memory>
#include <random>

namespace pic {
struct RunSummary {
    std::size_t steps_completed{0};
    double final_time{0.0};
    bool steady_state_reached{false};
    DiagnosticSample final_sample{};
};

struct CollisionDiagnostics {
    std::uint64_t candidates{0};
    std::uint64_t null_collisions{0};
    std::vector<std::string> channel_names{};
    std::vector<std::uint64_t> channel_collisions{};
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
    const CollisionDiagnostics& collision_diagnostics() const {
        return collision_totals_;
    }
private:
    void apply_collisions(Species& sp, std::size_t species_id);
    void deposit_and_solve();
    Config cfg_;
    Grid grid_;
    FieldSolver solver_;
    std::vector<Species> species_;
    std::unique_ptr<NullCollisionModel> mcc_model_;
    std::size_t mcc_species_id_{0};
    CollisionDiagnostics collision_totals_{};
    CollisionDiagnostics collision_interval_{};
    std::mt19937_64 rng_;
    double time_{0.0};
    std::size_t step_{0};
    bool initialized_{false};
};
}
