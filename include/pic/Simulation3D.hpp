#pragma once
#include "pic/Diagnostics.hpp"
#include "pic/FieldSolver.hpp"
#include "pic/ParticleState.hpp"
#include "pic/Runtime.hpp"
#include "pic/Species3D.hpp"
#include <filesystem>
#include <random>
#include <string>
#include <vector>

namespace pic {
struct RunSummary3D {
    std::size_t steps_completed{0};
    double final_time{0.0};
    bool steady_state_reached{false};
    DiagnosticSample3D final_sample{};
};

struct ParticleBoundaryConfig3D {
    ParticleBoundary left{ParticleBoundary::Auto};
    ParticleBoundary right{ParticleBoundary::Auto};
    ParticleBoundary bottom{ParticleBoundary::Auto};
    ParticleBoundary top{ParticleBoundary::Auto};
    ParticleBoundary back{ParticleBoundary::Auto};
    ParticleBoundary front{ParticleBoundary::Auto};
};

struct Simulation3DConfig {
    UnitSystemConfig units{};
    std::size_t nx{32};
    std::size_t ny{32};
    std::size_t nz{32};
    double length_x{1.0};
    double length_y{1.0};
    double length_z{1.0};
    double dt{0.02};
    std::size_t steps{100};
    RunMode mode{RunMode::Transient};
    double steady_tolerance{1e-6};
    std::size_t steady_window{25};
    std::size_t max_steps{10000};
    Boundary boundary{Boundary::Periodic};
    ParticleBoundaryConfig3D particle_boundary_config{};
    Vec3 magnetic_field{}; // uniform B used by the Boris pusher; zero keeps electrostatic leapfrog
    unsigned seed{12345};
    bool vtk_output{false};
    VTKOutputFormat vtk_format{VTKOutputFormat::Legacy};
    std::size_t output_interval{10};
    std::filesystem::path output_dir{"output3d"};
    bool particle_output{false};
    std::size_t particle_output_interval{0}; // zero inherits output_interval
    std::size_t particle_output_stride{1};
    std::size_t particle_sample_count{0}; // zero writes all stride-selected particles
    bool checkpoint_output{false};
    std::size_t checkpoint_interval{0}; // zero inherits output_interval
    std::filesystem::path checkpoint_path{}; // empty writes output_dir/checkpoint_<step>.apc
    std::filesystem::path restart_path{};
    std::filesystem::path initial_state_path{};
    std::optional<std::uint64_t> initial_state_signature{};
    RuntimePolicy runtime{};
    InitializationAcceptanceConfig initialization_acceptance{};
    std::vector<Species3DConfig> species{};
};

Simulation3DConfig load_config_3d(const std::string& path);

class Simulation3D {
public:
    explicit Simulation3D(Simulation3DConfig cfg);
    void initialize();
    void step();
    RunSummary3D run();
    void save_checkpoint(const std::filesystem::path& path) const;
    void load_checkpoint(const std::filesystem::path& path);
    DiagnosticSample3D sample() const;
    const Mesh3D& mesh() const { return mesh_; }
    const ParticleBoundaryConfig3D& particle_boundary_config() const { return cfg_.particle_boundary_config; }
    const BoundaryLoss3D& boundary_losses() const { return boundary_losses_; }
    const std::vector<Species3D>& species() const { return species_; }
    double time() const { return time_; }
    std::size_t step_count() const { return step_; }
private:
    void deposit_and_solve();
    void apply_particle_boundaries(Particle3D& particle);
    Simulation3DConfig cfg_;
    Mesh3D mesh_;
    FieldSolver solver_;
    std::vector<Species3D> species_;
    std::mt19937_64 rng_;
    double time_{0.0};
    std::size_t step_{0};
    BoundaryLoss3D boundary_losses_{};
    bool initialized_{false};
    ExternalParticleStateMetadata initial_state_metadata_{};
};
}
