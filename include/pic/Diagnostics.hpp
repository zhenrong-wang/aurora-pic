#pragma once
#include "pic/Grid.hpp"
#include "pic/Mesh2D.hpp"
#include "pic/Species.hpp"
#include "pic/Species2D.hpp"
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace pic {
struct DiagnosticSample {
    std::size_t step{0};
    double time{0.0};
    double kinetic_energy{0.0};
    double field_energy{0.0};
    double total_energy{0.0};
    double charge_l1{0.0};
    std::size_t live_particles{0};
};

struct BoundaryLoss2D {
    std::size_t absorbed_left{0};
    std::size_t absorbed_right{0};
    std::size_t absorbed_bottom{0};
    std::size_t absorbed_top{0};
};

struct DiagnosticSample2D {
    std::size_t step{0};
    double time{0.0};
    double kinetic_energy{0.0};
    double field_energy{0.0};
    double total_energy{0.0};
    double charge_l1{0.0};
    std::size_t live_particles{0};
    BoundaryLoss2D boundary_losses{};
    std::vector<std::size_t> live_particles_by_species{};
};

class Diagnostics {
public:
    explicit Diagnostics(std::filesystem::path output_dir);
    void write_header();
    DiagnosticSample sample(std::size_t step, double time, const Grid& grid, const std::vector<Species>& species);
    void write_sample(const DiagnosticSample& s);
    void write_fields(std::size_t step, const Grid& grid) const;
    const std::vector<DiagnosticSample>& history() const { return history_; }
private:
    std::filesystem::path output_dir_;
    std::ofstream scalar_file_;
    std::vector<DiagnosticSample> history_;
};

class Diagnostics2D {
public:
    Diagnostics2D(std::filesystem::path output_dir, const std::vector<Species2D>& species);
    void write_header();
    DiagnosticSample2D sample(std::size_t step,
                              double time,
                              const Mesh2D& mesh,
                              const std::vector<Species2D>& species,
                              BoundaryLoss2D boundary_losses = {});
    void write_sample(const DiagnosticSample2D& s);
    void write_particle_sample(std::size_t step,
                               const std::vector<Species2D>& species,
                               std::size_t stride,
                               std::size_t sample_count) const;
    const std::vector<DiagnosticSample2D>& history() const { return history_; }
private:
    std::filesystem::path output_dir_;
    std::ofstream scalar_file_;
    std::vector<std::string> species_names_;
    std::vector<DiagnosticSample2D> history_;
};
}
