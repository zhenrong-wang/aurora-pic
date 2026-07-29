#pragma once
#include "pic/Diagnostics.hpp"
#include "pic/FieldSolver.hpp"
#include "pic/ParticleState.hpp"
#include "pic/PrescribedField.hpp"
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

struct Simulation2DConfig {
    UnitSystemConfig units{};
    std::size_t nx{64};
    std::size_t ny{64};
    double length_x{1.0};
    double length_y{1.0};
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
    bool checkpoint_output{false};
    std::size_t checkpoint_interval{0}; // zero inherits output_interval
    std::filesystem::path checkpoint_path{}; // empty writes output_dir/checkpoint_<step>.apc
    std::filesystem::path restart_path{};
    std::filesystem::path initial_state_path{};
    std::optional<std::uint64_t> initial_state_signature{};
    RuntimePolicy runtime{};
    InitializationAcceptanceConfig initialization_acceptance{};
    std::vector<Species2DConfig> species{};
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
    double time() const { return time_; }
    std::size_t step_count() const { return step_; }
private:
    void deposit_and_solve();
    void apply_particle_boundaries(Particle2D& particle);
    Simulation2DConfig cfg_;
    Mesh2D mesh_;
    FieldSolver solver_;
    std::vector<Species2D> species_;
    std::mt19937_64 rng_;
    double time_{0.0};
    std::size_t step_{0};
    BoundaryLoss2D boundary_losses_{};
    bool initialized_{false};
    ExternalParticleStateMetadata initial_state_metadata_{};
};
}
