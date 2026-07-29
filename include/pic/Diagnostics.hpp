#pragma once
#include "pic/Grid.hpp"
#include "pic/Mesh2D.hpp"
#include "pic/Mesh3D.hpp"
#include "pic/Species.hpp"
#include "pic/Species2D.hpp"
#include "pic/Species3D.hpp"
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
    double phi_left{0.0};
    double phi_right{0.0};
};

struct BoundaryLoss2D {
    std::size_t absorbed_left{0};
    std::size_t absorbed_right{0};
    std::size_t absorbed_bottom{0};
    std::size_t absorbed_top{0};
};

struct BoundaryLoss3D {
    std::size_t absorbed_left{0};
    std::size_t absorbed_right{0};
    std::size_t absorbed_bottom{0};
    std::size_t absorbed_top{0};
    std::size_t absorbed_back{0};
    std::size_t absorbed_front{0};
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

struct DiagnosticSample3D {
    std::size_t step{0};
    double time{0.0};
    double kinetic_energy{0.0};
    double field_energy{0.0};
    double total_energy{0.0};
    double charge_l1{0.0};
    std::size_t live_particles{0};
    BoundaryLoss3D boundary_losses{};
    std::vector<std::size_t> live_particles_by_species{};
};

class Diagnostics {
public:
    explicit Diagnostics(
        std::filesystem::path output_dir,
        double permittivity = EPS0);
    void write_header();
    DiagnosticSample sample(std::size_t step, double time, const Grid& grid, const std::vector<Species>& species);
    void write_sample(const DiagnosticSample& s);
    void write_fields(std::size_t step, const Grid& grid) const;
    const std::vector<DiagnosticSample>& history() const { return history_; }
private:
    std::filesystem::path output_dir_;
    std::ofstream scalar_file_;
    double permittivity_{EPS0};
    std::vector<DiagnosticSample> history_;
};

class Diagnostics2D {
public:
    Diagnostics2D(
        std::filesystem::path output_dir,
        const std::vector<Species2D>& species,
        double permittivity = EPS0,
        double out_of_plane_depth = 1.0);
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
    double permittivity_{EPS0};
    double out_of_plane_depth_{1.0};
    std::vector<DiagnosticSample2D> history_;
};

class Diagnostics3D {
public:
    Diagnostics3D(
        std::filesystem::path output_dir,
        const std::vector<Species3D>& species,
        double permittivity = EPS0);
    void write_header();
    DiagnosticSample3D sample(std::size_t step,
                              double time,
                              const Mesh3D& mesh,
                              const std::vector<Species3D>& species,
                              BoundaryLoss3D boundary_losses = {});
    void write_sample(const DiagnosticSample3D& s);
    void write_particle_sample(std::size_t step,
                               const std::vector<Species3D>& species,
                               std::size_t stride,
                               std::size_t sample_count) const;
    const std::vector<DiagnosticSample3D>& history() const { return history_; }
private:
    std::filesystem::path output_dir_;
    std::ofstream scalar_file_;
    std::vector<std::string> species_names_;
    double permittivity_{EPS0};
    std::vector<DiagnosticSample3D> history_;
};
}
