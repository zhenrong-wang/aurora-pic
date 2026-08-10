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

struct BoundaryLoss1D {
    std::size_t absorbed_left{0};
    std::size_t absorbed_right{0};
    double kinetic_energy_left{0.0};
    double kinetic_energy_right{0.0};
};

struct SpeciesPower1D {
    double electric_work{0.0};
};

struct SpatialPhaseBin1D {
    std::size_t samples{0};
    std::vector<std::vector<double>> density{};
    std::vector<std::vector<double>> velocity_x_density{};
    std::vector<std::vector<double>> velocity_y_density{};
    std::vector<std::vector<double>> velocity_z_density{};
    std::vector<std::vector<double>> kinetic_energy_density{};
    std::vector<double> potential{};
    std::vector<double> electric{};
    std::vector<double> electric_squared{};
};

struct PhaseEedfAccumulator1D {
    std::uint64_t macro_observations{0};
    std::uint64_t overflow_macro_observations{0};
    double represented_observations{0.0};
    double overflow_represented_observations{0.0};
    double weighted_energy_sum{0.0};
    double weighted_energy_squared_sum{0.0};
    double weighted_velocity_x_sum{0.0};
    double weighted_velocity_y_sum{0.0};
    double weighted_velocity_z_sum{0.0};
    std::vector<double> histogram{};
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
    const std::vector<BoundaryLoss1D>&
    species_boundary_losses() const {
        return species_boundary_losses_;
    }
    std::size_t boundary_loss_origin_step() const {
        return boundary_loss_origin_step_;
    }
    const std::vector<SpeciesPower1D>&
    species_power_transfer() const {
        return species_power_transfer_;
    }
    std::size_t power_transfer_origin_step() const {
        return power_transfer_origin_step_;
    }
    const std::vector<std::vector<double>>&
    spatial_collision_energy_sums() const {
        return spatial_collision_energy_sums_;
    }
    const std::vector<std::vector<std::vector<double>>>&
    spatial_collision_phase_energy_sums() const {
        return spatial_collision_phase_energy_sums_;
    }
    std::size_t spatial_collision_steps() const {
        return spatial_collision_steps_;
    }
    const std::vector<std::size_t>&
    spatial_collision_phase_steps() const {
        return spatial_collision_phase_steps_;
    }
    const std::vector<std::vector<PhaseEedfAccumulator1D>>&
    phase_eedf_accumulators() const {
        return phase_eedf_accumulators_;
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
    bool species_due(std::size_t species_id) const;
    double species_timestep(std::size_t species_id) const;
    std::uint64_t collision_signature() const;
    std::string collision_identity() const;
    double electrode_potential(
        double offset,
        const SinusoidalVoltageConfig& drive,
        double field_time) const;
    void accumulate_spatial_average();
    void accumulate_phase_eedf(std::size_t phase);
    void deposit_spatial_collision_energy(
        double position,
        std::size_t channel,
        double represented_energy_change);
    void begin_spatial_collision_step();
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
    std::vector<BoundaryLoss1D> species_boundary_losses_{};
    std::vector<std::vector<BoundaryLoss1D>>
        boundary_loss_chunks_{};
    std::size_t boundary_loss_origin_step_{0};
    std::vector<SpeciesPower1D> species_power_transfer_{};
    std::vector<std::vector<SpeciesPower1D>>
        power_transfer_chunks_{};
    std::size_t power_transfer_origin_step_{0};
    std::size_t spatial_average_samples_{0};
    std::vector<std::vector<double>>
        spatial_density_sums_{};
    std::vector<double> spatial_density_scratch_{};
    std::size_t spatial_moment_samples_{0};
    std::vector<std::vector<double>>
        spatial_kinetic_energy_sums_{};
    std::vector<double> spatial_kinetic_energy_scratch_{};
    std::vector<double> spatial_velocity_x_scratch_{};
    std::vector<double> spatial_velocity_y_scratch_{};
    std::vector<double> spatial_velocity_z_scratch_{};
    std::vector<double> spatial_potential_sums_{};
    std::vector<double> spatial_electric_sums_{};
    std::vector<double> spatial_electric_squared_sums_{};
    std::vector<SpatialPhaseBin1D> spatial_phase_bins_{};
    std::size_t spatial_collision_steps_{0};
    std::vector<std::vector<double>>
        spatial_collision_energy_sums_{};
    std::vector<std::size_t> spatial_collision_phase_steps_{};
    std::vector<std::vector<std::vector<double>>>
        spatial_collision_phase_energy_sums_{};
    bool spatial_collision_step_active_{false};
    std::size_t spatial_collision_active_phase_{0};
    std::size_t phase_eedf_species_id_{0};
    std::vector<std::vector<PhaseEedfAccumulator1D>>
        phase_eedf_accumulators_{};
    std::mt19937_64 rng_;
    double time_{0.0};
    std::size_t step_{0};
    bool initialized_{false};
    ExternalParticleStateMetadata initial_state_metadata_{};
};
}
