#pragma once

#include "pic/Collision.hpp"
#include "pic/Runtime.hpp"
#include "pic/Species2D.hpp"
#include "pic/UnstructuredFieldSolver2D.hpp"

#include <filesystem>
#include <iosfwd>
#include <array>
#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <vector>

namespace pic {

struct UnstructuredSpecies2DConfig {
    std::string name{"electrons"};
    double charge{-1.0};
    double mass{1.0};
    double weight{1.0};
    std::size_t particles{1000};
    double drift_velocity_x{0.0};
    double drift_velocity_y{0.0};
    double thermal_velocity{0.1};
    std::optional<Vec2> initialization_minimum;
    std::optional<Vec2> initialization_maximum;
    std::string initialization_region;
    double drift_velocity_z{0.0};
    ParticleInitializationConfig initialization{};
};

struct UnstructuredBoundarySource2DConfig {
    std::string name;
    std::string species;
    std::string boundary;
    std::size_t particles_per_step{0};
    std::size_t start_step{0};
    std::size_t end_step{0}; // zero keeps the source active indefinitely
    double normal_velocity{0.0};
    double tangential_velocity{0.0};
    double thermal_velocity{0.0};
    double out_of_plane_velocity{0.0};
};

struct UnstructuredSecondaryEmission2DConfig {
    std::string name;
    std::string boundary;
    std::string incident_species;
    std::string emitted_species;
    double yield{0.0};
    std::size_t max_particles_per_impact{1000};
    double normal_velocity{0.0};
    double tangential_velocity{0.0};
    double thermal_velocity{0.0};
    double out_of_plane_velocity{0.0};
};

struct UnstructuredBoundaryFlux2D {
    std::size_t macroparticles{0};
    double physical_particles{0.0};
    double charge{0.0};
    double kinetic_energy{0.0};
    std::size_t last_step_macroparticles{0};
    double last_step_physical_particles{0.0};
    double physical_particle_rate{0.0};
    double physical_particle_flux{0.0};
};

struct UnstructuredSimulation2DConfig {
    UnitSystemConfig units{};
    std::filesystem::path mesh_path;
    double dt{0.02};
    std::size_t steps{100};
    RunMode mode{RunMode::Transient};
    double steady_tolerance{1e-6};
    std::size_t steady_window{25};
    std::size_t max_steps{10000};
    std::size_t max_particles_per_species{10000000};
    unsigned seed{12345};
    double magnetic_field_z{0.0};
    double magnetic_field_x{0.0};
    double magnetic_field_y{0.0};
    std::size_t output_interval{10};
    std::filesystem::path output_dir{"output_unstructured_2d"};
    bool vtk_output{false};
    bool particle_output{false};
    std::size_t particle_output_interval{0};
    std::size_t particle_output_stride{1};
    std::size_t particle_sample_count{0};
    bool checkpoint_output{false};
    std::size_t checkpoint_interval{0};
    std::filesystem::path checkpoint_path;
    std::filesystem::path restart_path;
    std::filesystem::path initial_state_path;
    RuntimePolicy runtime{};
    InitializationAcceptanceConfig initialization_acceptance{};
    UnstructuredPoissonOptions2D poisson{};
    std::map<std::string, double> dirichlet_potentials;
    std::map<std::string, double> neumann_normal_derivatives;
    std::map<std::string, ParticleBoundary> particle_boundaries;
    std::vector<UnstructuredSpecies2DConfig> species;
    std::vector<UnstructuredBoundarySource2DConfig> sources;
    std::vector<UnstructuredSecondaryEmission2DConfig> emissions;
    CollisionConfig collisions{};
};

struct UnstructuredDiagnosticSample2D {
    std::size_t step{0};
    double time{0.0};
    double kinetic_energy{0.0};
    double field_energy{0.0};
    double total_energy{0.0};
    double charge_l1{0.0};
    std::size_t live_particles{0};
    std::map<std::string, std::size_t> absorbed_by_label;
    std::map<std::string, std::size_t> injected_by_source;
    std::map<std::string, std::size_t> emitted_by_rule;
    std::map<std::string,
             std::map<std::string, UnstructuredBoundaryFlux2D>> impact_flux;
    UnstructuredPoissonSummary2D poisson{};
};

struct UnstructuredRunSummary2D {
    std::size_t steps_completed{0};
    double final_time{0.0};
    bool steady_state_reached{false};
    UnstructuredDiagnosticSample2D final_sample{};
};

struct UnstructuredTiming2D {
    double particle_seconds{0.0};
    double deposition_seconds{0.0};
    double field_solve_seconds{0.0};
    std::size_t location_cache_hits{0};
    std::size_t location_searches{0};
};

class UnstructuredSimulation2D {
public:
    explicit UnstructuredSimulation2D(UnstructuredSimulation2DConfig config);

    void initialize();
    void step();
    UnstructuredRunSummary2D run();
    void save_checkpoint(const std::filesystem::path& path) const;
    void load_checkpoint(const std::filesystem::path& path);
    UnstructuredDiagnosticSample2D sample() const;

    const UnstructuredMesh2D& mesh() const { return mesh_; }
    const std::vector<Species2D>& species() const { return species_; }
    const std::map<std::string, std::size_t>& absorbed_by_label() const {
        return absorbed_by_label_;
    }
    double time() const { return time_; }
    std::size_t step_count() const { return step_; }
    std::size_t poisson_assembly_count() const {
        return poisson_solver_ ? poisson_solver_->assembly_count() : 0;
    }
    std::size_t poisson_solve_count() const {
        return poisson_solver_ ? poisson_solver_->solve_count() : 0;
    }
    const UnstructuredTiming2D& timing() const { return timing_; }
    const CollisionDiagnostics& collision_diagnostics() const {
        return collision_totals_;
    }

private:
    struct BoundarySegment {
        Vec2 first{};
        Vec2 second{};
        Vec2 inward_normal{};
        std::string label;
    };
    struct SamplingTriangle {
        std::array<Vec2, 3> vertices{};
        double cumulative_area{0.0};
    };
    struct RegionSamplingTriangle {
        std::size_t sampling_triangle_index{0};
        double cumulative_area{0.0};
    };
    struct BoundarySourceRuntime {
        UnstructuredBoundarySource2DConfig config;
        std::size_t species_id{0};
        std::vector<std::size_t> segment_indices;
        std::vector<double> cumulative_lengths;
        std::size_t injected_particles{0};
    };
    struct SecondaryEmissionRuntime {
        UnstructuredSecondaryEmission2DConfig config;
        std::size_t incident_species_id{0};
        std::size_t emitted_species_id{0};
        std::size_t emitted_particles{0};
    };
    struct IonizationChannelRuntime {
        std::size_t secondary_species_id{0};
        std::size_t ion_species_id{0};
    };
    struct AttachmentChannelRuntime {
        std::size_t product_species_id{0};
    };
    struct BoundaryImpact {
        std::size_t species_id{0};
        std::size_t particle_id{0};
        std::size_t segment_id{0};
        Vec2 position{};
        Vec2 incident_velocity{};
        double incident_velocity_z{0.0};
    };

    Vec2 sample_position(
        const UnstructuredSpecies2DConfig& config,
        std::size_t particle_index,
        std::size_t particle_count,
        std::size_t& profile_attempts);
    void inject_boundary_sources();
    void process_boundary_impacts(std::vector<BoundaryImpact> impacts);
    void apply_collisions();
    void deposit_and_solve();
    std::optional<BoundaryImpact> advance_with_boundaries(
        Particle2D& particle, Vec2 previous_position);
    void write_diagnostics_header(std::ofstream& output) const;
    void write_diagnostics_sample(std::ofstream& output,
                                  const UnstructuredDiagnosticSample2D& sample) const;
    void write_particle_sample(std::size_t step) const;
    std::filesystem::path checkpoint_path_for_step(std::size_t step) const;
    std::uint64_t mesh_signature() const;

    UnstructuredSimulation2DConfig config_;
    UnstructuredMesh2D mesh_;
    std::unique_ptr<UnstructuredPoissonSolver2D> poisson_solver_;
    std::vector<Species2D> species_;
    std::vector<UnstructuredSpecies2DConfig> species_configs_;
    std::vector<std::vector<UnstructuredParticleLocation2D>> particle_locations_;
    std::vector<BoundarySegment> boundary_segments_;
    std::vector<SamplingTriangle> sampling_triangles_;
    std::map<std::string, std::vector<RegionSamplingTriangle>>
        region_sampling_triangles_;
    std::map<std::string, std::pair<Vec2, Vec2>>
        region_sampling_bounds_;
    std::vector<BoundarySourceRuntime> sources_;
    std::vector<SecondaryEmissionRuntime> emissions_;
    std::vector<std::optional<IonizationChannelRuntime>>
        ionization_channels_;
    std::vector<std::optional<AttachmentChannelRuntime>>
        attachment_channels_;
    std::unique_ptr<NullCollisionModel> mcc_model_;
    std::size_t mcc_species_id_{0};
    CollisionDiagnostics collision_totals_{};
    CollisionDiagnostics collision_interval_{};
    std::map<std::string, std::size_t> absorbed_by_label_;
    std::map<std::string, double> boundary_lengths_;
    std::map<std::string,
             std::map<std::string, UnstructuredBoundaryFlux2D>> impact_flux_;
    std::mt19937_64 rng_;
    UnstructuredPoissonSummary2D last_poisson_{};
    UnstructuredTiming2D timing_{};
    double time_{0.0};
    std::size_t step_{0};
    bool initialized_{false};
};

bool config_uses_unstructured_mesh_2d(const std::filesystem::path& path);
UnstructuredSimulation2DConfig load_unstructured_config_2d(
    const std::filesystem::path& path);

} // namespace pic
