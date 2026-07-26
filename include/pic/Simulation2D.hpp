#pragma once
#include "pic/Diagnostics.hpp"
#include "pic/FieldSolver.hpp"
#include "pic/Species2D.hpp"
#include <filesystem>
#include <random>
#include <string>
#include <vector>

namespace pic {
struct RunSummary2D {
    std::size_t steps_completed{0};
    double final_time{0.0};
    DiagnosticSample2D final_sample{};
};

struct ParticleBoundaryConfig2D {
    ParticleBoundary left{ParticleBoundary::Auto};
    ParticleBoundary right{ParticleBoundary::Auto};
    ParticleBoundary bottom{ParticleBoundary::Auto};
    ParticleBoundary top{ParticleBoundary::Auto};
};

struct Simulation2DConfig {
    std::size_t nx{64};
    std::size_t ny{64};
    double length_x{1.0};
    double length_y{1.0};
    double dt{0.02};
    std::size_t steps{100};
    Boundary boundary{Boundary::Periodic};
    BoundaryConfig2D boundary_config{};
    ParticleBoundaryConfig2D particle_boundary_config{};
    unsigned seed{12345};
    bool vtk_output{false};
    std::size_t output_interval{10};
    std::filesystem::path output_dir{"output2d"};
    bool particle_output{false};
    std::size_t particle_output_interval{0}; // zero inherits output_interval
    std::size_t particle_output_stride{1};
    std::size_t particle_sample_count{0}; // zero writes all stride-selected particles
    std::vector<Species2DConfig> species{};
};

Simulation2DConfig load_config_2d(const std::string& path);

class Simulation2D {
public:
    explicit Simulation2D(Simulation2DConfig cfg);
    void initialize();
    void step();
    RunSummary2D run();
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
};
}
