#pragma once

#include "pic/Runtime.hpp"
#include "pic/Species2D.hpp"
#include "pic/UnstructuredFieldSolver2D.hpp"

#include <filesystem>
#include <iosfwd>
#include <array>
#include <map>
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
};

struct UnstructuredSimulation2DConfig {
    std::filesystem::path mesh_path;
    double dt{0.02};
    std::size_t steps{100};
    RunMode mode{RunMode::Transient};
    double steady_tolerance{1e-6};
    std::size_t steady_window{25};
    std::size_t max_steps{10000};
    unsigned seed{12345};
    double magnetic_field_z{0.0};
    std::size_t output_interval{10};
    std::filesystem::path output_dir{"output_unstructured_2d"};
    bool vtk_output{false};
    RuntimePolicy runtime{};
    UnstructuredPoissonOptions2D poisson{};
    std::map<std::string, double> dirichlet_potentials;
    std::map<std::string, ParticleBoundary> particle_boundaries;
    std::vector<UnstructuredSpecies2DConfig> species;
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
    UnstructuredPoissonSummary2D poisson{};
};

struct UnstructuredRunSummary2D {
    std::size_t steps_completed{0};
    double final_time{0.0};
    bool steady_state_reached{false};
    UnstructuredDiagnosticSample2D final_sample{};
};

class UnstructuredSimulation2D {
public:
    explicit UnstructuredSimulation2D(UnstructuredSimulation2DConfig config);

    void initialize();
    void step();
    UnstructuredRunSummary2D run();
    UnstructuredDiagnosticSample2D sample() const;

    const UnstructuredMesh2D& mesh() const { return mesh_; }
    const std::vector<Species2D>& species() const { return species_; }
    const std::map<std::string, std::size_t>& absorbed_by_label() const {
        return absorbed_by_label_;
    }
    double time() const { return time_; }
    std::size_t step_count() const { return step_; }

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

    Vec2 sample_position(const UnstructuredSpecies2DConfig& config);
    void deposit_and_solve();
    void advance_with_boundaries(Particle2D& particle, Vec2 previous_position);
    void write_diagnostics_header(std::ofstream& output) const;
    void write_diagnostics_sample(std::ofstream& output,
                                  const UnstructuredDiagnosticSample2D& sample) const;

    UnstructuredSimulation2DConfig config_;
    UnstructuredMesh2D mesh_;
    std::vector<Species2D> species_;
    std::vector<UnstructuredSpecies2DConfig> species_configs_;
    std::vector<BoundarySegment> boundary_segments_;
    std::vector<SamplingTriangle> sampling_triangles_;
    std::map<std::string, std::size_t> absorbed_by_label_;
    std::mt19937_64 rng_;
    UnstructuredPoissonSummary2D last_poisson_{};
    double time_{0.0};
    std::size_t step_{0};
    bool initialized_{false};
};

} // namespace pic
