#pragma once
#include "pic/Collision.hpp"
#include "pic/Config.hpp"
#include "pic/Convergence.hpp"
#include "pic/Diagnostics.hpp"
#include "pic/FieldSolver.hpp"
#include "pic/ParticleState.hpp"
#include "pic/Species.hpp"
#include <array>
#include <filesystem>
#include <memory>
#include <optional>
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

struct WallImpactSideSpectrum1D {
    std::uint64_t macro_impacts{0};
    std::uint64_t overflow_macro_impacts{0};
    double represented_impacts{0.0};
    double overflow_represented_impacts{0.0};
    double represented_kinetic_energy{0.0};
    std::vector<std::uint64_t> macro_histogram{};
    std::vector<double> represented_histogram{};
};

struct SpeciesWallImpactSpectrum1D {
    BoundaryLoss1D baseline_loss{};
    WallImpactSideSpectrum1D left{};
    WallImpactSideSpectrum1D right{};
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
    double weighted_velocity_x_squared_sum{0.0};
    double weighted_velocity_y_squared_sum{0.0};
    double weighted_velocity_z_squared_sum{0.0};
    double tail_represented_observations{0.0};
    double tail_positive_x_represented_observations{0.0};
    double tail_negative_x_represented_observations{0.0};
    double tail_weighted_velocity_x_sum{0.0};
    double tail_weighted_velocity_x_squared_sum{0.0};
    double tail_weighted_transverse_velocity_squared_sum{0.0};
    double tail_weighted_age_steps_sum{0.0};
    double tail_weighted_energetic_steps_sum{0.0};
    double tail_weighted_energetic_duty_fraction_sum{0.0};
    double tail_weighted_consecutive_energetic_steps_sum{0.0};
    double tail_weighted_entries_sum{0.0};
    double tail_weighted_elastic_collisions_sum{0.0};
    double tail_weighted_excitation_collisions_sum{0.0};
    double tail_weighted_ionization_collisions_sum{0.0};
    double tail_weighted_charge_exchange_collisions_sum{0.0};
    double tail_weighted_bgk_collisions_sum{0.0};
    double tail_born_during_window_represented_observations{0.0};
    std::vector<double> histogram{};
};

struct ParticleHistory1D {
    std::uint64_t age_steps{0};
    std::uint64_t energetic_steps{0};
    std::uint64_t consecutive_energetic_steps{0};
    std::uint64_t tail_entries{0};
    std::uint64_t elastic_collisions{0};
    std::uint64_t excitation_collisions{0};
    std::uint64_t ionization_collisions{0};
    std::uint64_t charge_exchange_collisions{0};
    std::uint64_t bgk_collisions{0};
    bool born_during_window{false};
    bool energetic_previous_step{false};
};

struct PhaseEedfThresholdCrossingAccumulator1D {
    std::uint64_t electron_time_macro_observations{0};
    std::uint64_t energetic_time_macro_observations{0};
    std::uint64_t interstep_promotions{0};
    std::uint64_t interstep_demotions{0};
    std::uint64_t field_push_macro_observations{0};
    std::uint64_t field_push_promotions{0};
    std::uint64_t field_push_demotions{0};
    std::uint64_t field_push_promotion_band_observations{0};
    std::uint64_t field_push_promotion_band_promotions{0};
    double field_push_promotion_band_signed_work{0.0};
    double field_push_promotion_band_positive_work{0.0};
    double field_push_promotion_band_negative_work{0.0};
    double field_push_promotion_band_origin_energy{0.0};
    double field_push_promotion_band_origin_longitudinal_energy{0.0};
    double field_push_promotion_band_linear_work{0.0};
    double field_push_promotion_band_positive_linear_work{0.0};
    double field_push_promotion_band_negative_linear_work{0.0};
    double field_push_promotion_band_quadratic_work{0.0};
    std::array<std::uint64_t, 6> collision_promotions{};
    std::array<std::uint64_t, 6> collision_demotions{};
    std::uint64_t energetic_births{0};
    std::uint64_t subthreshold_births{0};
};

struct PhaseSurfaceFluxAccumulator1D {
    std::uint64_t macro_crossings{0};
    std::uint64_t overflow_macro_crossings{0};
    double represented_crossings{0.0};
    double overflow_represented_crossings{0.0};
    double represented_kinetic_energy{0.0};
    std::vector<double> represented_histogram{};
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
    const std::vector<SpeciesWallImpactSpectrum1D>&
    wall_impact_spectra() const {
        return wall_impact_spectra_;
    }
    std::size_t wall_impact_origin_step() const {
        return wall_impact_origin_step_;
    }
    const std::vector<SpeciesPower1D>&
    species_power_transfer() const {
        return species_power_transfer_;
    }
    std::size_t power_transfer_origin_step() const {
        return power_transfer_origin_step_;
    }
    std::vector<BlockConvergenceResult>
    periodic_convergence_results() const {
        return periodic_convergence_
            ? periodic_convergence_->evaluate()
            : std::vector<BlockConvergenceResult>{};
    }
    const std::vector<std::vector<double>>&
    spatial_collision_energy_sums() const {
        return spatial_collision_energy_sums_;
    }
    const std::vector<std::vector<std::vector<double>>>&
    spatial_collision_phase_energy_sums() const {
        return spatial_collision_phase_energy_sums_;
    }
    const std::vector<std::vector<double>>&
    spatial_collision_event_sums() const {
        return spatial_collision_event_sums_;
    }
    const std::vector<std::vector<std::vector<double>>>&
    spatial_collision_phase_event_sums() const {
        return spatial_collision_phase_event_sums_;
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
    const std::vector<ParticleHistory1D>&
    phase_eedf_particle_histories() const {
        return phase_eedf_particle_histories_;
    }
    const auto& phase_eedf_threshold_crossings() const {
        return phase_eedf_threshold_crossings_;
    }
    const auto& phase_surface_flux_accumulators() const {
        return phase_surface_flux_accumulators_;
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
        std::vector<CollisionProcessKind> channel_processes{};
        CollisionWorkspace collision_workspace{};
        std::vector<double> tracked_energy_scratch{};
    };
    void apply_collisions();
    void deposit_and_solve(double field_time);
    void refresh_held_charge(std::size_t species_id);
    bool species_due(std::size_t species_id) const;
    double species_timestep(std::size_t species_id) const;
    std::uint64_t collision_signature() const;
    std::string collision_identity() const;
    double electrode_potential(
        double offset,
        const SinusoidalVoltageConfig& drive,
        double field_time) const;
    void accumulate_spatial_average(std::size_t sample_step);
    void accumulate_phase_eedf(std::size_t phase);
    bool phase_eedf_history_active() const;
    std::size_t phase_eedf_history_phase() const;
    void update_phase_eedf_histories();
    double phase_eedf_collision_state_energy(
        const Particle& particle, const Species& species) const;
    bool phase_eedf_collision_state_energetic(
        const Particle& particle, const Species& species) const;
    void add_phase_eedf_collision_history(
        std::size_t particle_id,
        CollisionProcessKind process,
        std::uint64_t count);
    void add_phase_eedf_collision_transition(
        double position, CollisionProcessKind process,
        bool energetic_before, bool energetic_after);
    void add_phase_eedf_bgk_transition(
        double position, bool energetic_before, bool energetic_after);
    void add_phase_eedf_birth(double position, bool energetic);
    std::size_t phase_surface_flux_phase(std::size_t sample_step) const;
    void accumulate_phase_surface_crossing(
        std::size_t chunk, std::size_t surface, std::size_t direction,
        const Particle& particle, const Species& species);
    void merge_phase_surface_flux_chunks(std::size_t phase);
    void deposit_spatial_collision_energy(
        double position,
        std::size_t channel,
        double represented_energy_change);
    void deposit_spatial_collision_events(
        double position,
        std::size_t channel,
        double represented_events);
    void begin_spatial_collision_step();
    void accumulate_wall_impact(
        WallImpactSideSpectrum1D& accumulator,
        std::size_t species_id,
        double particle_energy,
        double represented_energy) const;
    void write_spatial_average() const;
    void write_wall_impact_spectrum() const;
    void reset_periodic_convergence(std::size_t origin_cycle = 0);
    bool sample_periodic_convergence();
    void write_periodic_convergence() const;
    std::size_t expected_spatial_average_samples() const;
    Config cfg_;
    Grid grid_;
    FieldSolver solver_;
    std::vector<Species> species_;
    std::vector<std::vector<double>> held_charge_density_{};
    std::vector<MccRuntime> mcc_models_;
    bool legacy_bgk_enabled_{false};
    CollisionDiagnostics collision_totals_{};
    CollisionDiagnostics collision_interval_{};
    std::vector<BoundaryLoss1D> species_boundary_losses_{};
    std::vector<std::vector<BoundaryLoss1D>>
        boundary_loss_chunks_{};
    std::size_t boundary_loss_origin_step_{0};
    std::size_t wall_impact_origin_step_{0};
    std::vector<SpeciesWallImpactSpectrum1D>
        wall_impact_spectra_{};
    std::vector<std::vector<SpeciesWallImpactSpectrum1D>>
        wall_impact_chunks_{};
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
    std::vector<std::vector<double>>
        spatial_collision_event_sums_{};
    std::vector<std::vector<std::vector<double>>>
        spatial_collision_phase_event_sums_{};
    bool spatial_collision_step_active_{false};
    std::size_t spatial_collision_active_phase_{0};
    std::size_t phase_eedf_species_id_{0};
    std::vector<std::vector<PhaseEedfAccumulator1D>>
        phase_eedf_accumulators_{};
    std::vector<ParticleHistory1D> phase_eedf_particle_histories_{};
    std::vector<double> phase_eedf_field_push_origin_energy_{};
    std::vector<double>
        phase_eedf_field_push_origin_longitudinal_velocity_{};
    std::vector<std::vector<PhaseEedfThresholdCrossingAccumulator1D>>
        phase_eedf_threshold_crossings_{};
    std::size_t phase_surface_flux_species_id_{0};
    std::vector<std::vector<std::vector<PhaseSurfaceFluxAccumulator1D>>>
        phase_surface_flux_accumulators_{};
    std::vector<std::vector<std::vector<PhaseSurfaceFluxAccumulator1D>>>
        phase_surface_flux_chunks_{};
    std::optional<PeriodicBlockConvergence> periodic_convergence_{};
    std::size_t periodic_convergence_steps_per_cycle_{0};
    bool periodic_convergence_block_closed_{false};
    std::mt19937_64 rng_;
    double time_{0.0};
    std::size_t step_{0};
    bool initialized_{false};
    ExternalParticleStateMetadata initial_state_metadata_{};
};
}
