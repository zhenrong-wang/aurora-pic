#pragma once
#include "pic/Diagnostics.hpp"
#include "pic/FieldSolver.hpp"
#include "pic/ParticleState.hpp"
#include "pic/PrescribedField.hpp"
#include "pic/ResolvedDiagnostics2D.hpp"
#include "pic/Runtime.hpp"
#include "pic/Species2D.hpp"
#include <filesystem>
#include <optional>
#include <random>
#include <string>
#include <vector>

namespace pic {
struct RunSummary2D {
    std::size_t steps_completed{0};
    double final_time{0.0};
    bool steady_state_reached{false};
    DiagnosticSample2D final_sample{};
};

struct ParticleBoundaryConfig2D {
    ParticleBoundary left{ParticleBoundary::Auto};
    ParticleBoundary right{ParticleBoundary::Auto};
    ParticleBoundary bottom{ParticleBoundary::Auto};
    ParticleBoundary top{ParticleBoundary::Auto};
};

struct VolumetricPairSource2DConfig {
    std::string name{};
    std::string first_species{};
    std::string second_species{};
    std::size_t pairs_per_step{0};
    std::optional<double> represented_pair_rate{};
    std::optional<double> peak_volumetric_pair_rate{};
    std::size_t start_step{0};
    std::size_t end_step{0}; // zero keeps the source active indefinitely
    double x_min{0.0};
    double x_max{-1.0}; // negative means the full remaining domain
    double y_min{0.0};
    double y_max{-1.0};
    Vec3 first_drift{};
    Vec3 second_drift{};
    double first_thermal_velocity{0.0};
    double second_thermal_velocity{0.0};
    ParticleInitializationConfig spatial_profile{};
};

struct VolumetricPairSource2DDiagnostics {
    std::string name{};
    std::size_t macro_pairs_created{0};
    double represented_pairs_created{0.0};
    double fractional_macro_pair_remainder{0.0};
    double injected_kinetic_energy{0.0};
};

enum class BoundarySide2DName { Left, Right, Bottom, Top };

std::string to_string(BoundarySide2DName side);

enum class CurrentSourceControlMode {
    Cumulative,
    TimestepLocal
};

std::string to_string(CurrentSourceControlMode mode);

struct CurrentRegulatedSource2DConfig {
    std::string species{};
    BoundarySide2DName monitor_boundary{BoundarySide2DName::Left};
    BoundarySide2DName emission_boundary{BoundarySide2DName::Right};
    double emission_inset{0.0};
    Vec3 drift{};
    double thermal_velocity{0.0};
    CurrentSourceControlMode control_mode{
        CurrentSourceControlMode::Cumulative};
};

struct CurrentRegulatedSource2DDiagnostics {
    std::size_t macro_particles_created{0};
    double represented_particles_created{0.0};
    double control_macro_remainder{0.0};
    double processed_monitored_charge{0.0};
    double processed_monitored_negative_charge{0.0};
    double processed_monitored_positive_charge{0.0};
    double injected_kinetic_energy{0.0};
    std::size_t control_updates{0};
    std::size_t reverse_demand_steps{0};
    std::size_t reverse_diagnostics_start_step{0};
    double cumulative_reverse_demand_macroparticles{0.0};
    double maximum_reverse_demand_macroparticles{0.0};
    std::size_t reverse_distribution_start_step{0};
    std::size_t reverse_one_macro_steps{0};
    std::size_t reverse_two_macro_steps{0};
    std::size_t reverse_multi_macro_steps{0};
    double distributed_reverse_demand_macroparticles{0.0};
    double squared_reverse_demand_macroparticles{0.0};
    double reverse_monitored_negative_charge{0.0};
    double reverse_monitored_positive_charge{0.0};
};

enum class PotentialReferenceCorrection {
    Gauge,
    Affine
};

std::string to_string(PotentialReferenceCorrection correction);

struct PotentialReference2DConfig {
    CoordinateAxis axis{CoordinateAxis::X};
    double coordinate{0.0};
    double target{0.0};
    PotentialReferenceCorrection correction{
        PotentialReferenceCorrection::Gauge};
};

struct Simulation2DConfig {
    UnitSystemConfig units{};
    std::size_t nx{64};
    std::size_t ny{64};
    double length_x{1.0};
    double length_y{1.0};
    double out_of_plane_depth{1.0};
    double dt{0.02};
    std::size_t steps{100};
    RunMode mode{RunMode::Transient};
    double steady_tolerance{1e-6};
    std::size_t steady_window{25};
    std::size_t max_steps{10000};
    Boundary boundary{Boundary::Periodic};
    std::optional<Boundary> boundary_x{};
    std::optional<Boundary> boundary_y{};
    BoundaryConfig2D boundary_config{};
    ParticleBoundaryConfig2D particle_boundary_config{};
    double magnetic_field_z{0.0};
    double magnetic_field_x{0.0};
    double magnetic_field_y{0.0}; // any nonzero uniform B component activates the 2D3V Boris pusher
    std::optional<TabulatedVectorField1D> magnetic_field_profile{};
    unsigned seed{12345};
    bool vtk_output{false};
    VTKOutputFormat vtk_format{VTKOutputFormat::Legacy};
    std::size_t output_interval{10};
    std::filesystem::path output_dir{"output2d"};
    bool particle_output{false};
    std::size_t particle_output_interval{0}; // zero inherits output_interval
    std::size_t particle_output_stride{1};
    std::size_t particle_sample_count{0}; // zero writes all stride-selected particles
    ResolvedDiagnostics2DConfig resolved_diagnostics{};
    bool checkpoint_output{false};
    std::size_t checkpoint_interval{0}; // zero inherits output_interval
    std::filesystem::path checkpoint_path{}; // empty writes output_dir/checkpoint_<step>.apc
    std::filesystem::path restart_path{};
    std::filesystem::path initial_state_path{};
    std::optional<std::uint64_t> initial_state_signature{};
    std::size_t max_particles_per_species{10000000};
    RuntimePolicy runtime{};
    InitializationAcceptanceConfig initialization_acceptance{};
    std::vector<Species2DConfig> species{};
    std::vector<VolumetricPairSource2DConfig> sources{};
    std::optional<CurrentRegulatedSource2DConfig>
        current_regulated_source{};
    std::optional<PotentialReference2DConfig> potential_reference{};
};

Simulation2DConfig load_config_2d(const std::string& path);

class Simulation2D {
public:
    explicit Simulation2D(Simulation2DConfig cfg);
    void initialize();
    void step();
    RunSummary2D run();
    void save_checkpoint(const std::filesystem::path& path) const;
    void load_checkpoint(const std::filesystem::path& path);
    DiagnosticSample2D sample() const;
    const Mesh2D& mesh() const { return mesh_; }
    const ParticleBoundaryConfig2D& particle_boundary_config() const { return cfg_.particle_boundary_config; }
    const BoundaryLoss2D& boundary_losses() const { return boundary_losses_; }
    const std::vector<Species2D>& species() const { return species_; }
    const std::vector<VolumetricPairSource2DDiagnostics>& source_diagnostics() const {
        return source_diagnostics_;
    }
    const std::optional<CurrentRegulatedSource2DDiagnostics>&
    current_regulated_source_diagnostics() const {
        return current_regulated_source_diagnostics_;
    }
    const std::vector<BoundaryLoss2D>&
    species_boundary_losses() const {
        return species_boundary_losses_;
    }
    double potential_reference_offset() const {
        return potential_reference_offset_;
    }
    double time() const { return time_; }
    std::size_t step_count() const { return step_; }
private:
    struct VolumetricPairSourceRuntime {
        VolumetricPairSource2DConfig config{};
        std::size_t first_species{0};
        std::size_t second_species{0};
        double represented_pair_rate{0.0};
        double effective_profile_area{0.0};
    };
    void deposit_and_solve();
    void inject_volumetric_pair_sources();
    void inject_current_regulated_source();
    void apply_potential_reference();
    void apply_particle_boundaries(
        Particle2D& particle, std::size_t species_id);
    Simulation2DConfig cfg_;
    Mesh2D mesh_;
    FieldSolver solver_;
    std::vector<Species2D> species_;
    std::vector<VolumetricPairSourceRuntime> sources_;
    std::vector<VolumetricPairSource2DDiagnostics> source_diagnostics_;
    std::optional<std::size_t> current_regulated_species_{};
    std::optional<CurrentRegulatedSource2DDiagnostics>
        current_regulated_source_diagnostics_{};
    std::vector<BoundaryLoss2D> species_boundary_losses_{};
    double potential_reference_offset_{0.0};
    std::mt19937_64 rng_;
    double time_{0.0};
    std::size_t step_{0};
    BoundaryLoss2D boundary_losses_{};
    bool initialized_{false};
    ExternalParticleStateMetadata initial_state_metadata_{};
};
}
