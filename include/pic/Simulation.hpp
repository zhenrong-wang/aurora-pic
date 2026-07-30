#pragma once
#include "pic/Collision.hpp"
#include "pic/Config.hpp"
#include "pic/Diagnostics.hpp"
#include "pic/FieldSolver.hpp"
#include "pic/ParticleState.hpp"
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
    struct IonizationChannelRuntime {
        std::size_t secondary_species_id{0};
        std::size_t ion_species_id{0};
    };
    struct MccRuntime {
        std::string name{};
        std::size_t species_id{0};
        std::size_t diagnostic_offset{0};
        std::unique_ptr<NullCollisionModel> model{};
        std::vector<std::optional<IonizationChannelRuntime>>
            ionization_channels{};
    };
    void apply_collisions();
    void deposit_and_solve(double field_time);
    std::uint64_t collision_signature() const;
    std::string collision_identity() const;
    double electrode_potential(
        double offset,
        const SinusoidalVoltageConfig& drive,
        double field_time) const;
    void accumulate_spatial_average();
    void write_spatial_average() const;
    std::size_t expected_spatial_average_samples() const;
    Config cfg_;
    Grid grid_;
    FieldSolver solver_;
    std::vector<Species> species_;
    std::vector<MccRuntime> mcc_models_;
    bool legacy_bgk_enabled_{false};
    CollisionDiagnostics collision_totals_{};
    CollisionDiagnostics collision_interval_{};
    std::size_t spatial_average_samples_{0};
    std::vector<std::vector<double>>
        spatial_density_sums_{};
    std::vector<double> spatial_density_scratch_{};
    std::mt19937_64 rng_;
    double time_{0.0};
    std::size_t step_{0};
    bool initialized_{false};
    ExternalParticleStateMetadata initial_state_metadata_{};
};
}
