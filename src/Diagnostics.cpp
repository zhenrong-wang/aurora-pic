#include "pic/Diagnostics.hpp"
#include <cmath>
#include <iomanip>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

namespace pic {
Diagnostics::Diagnostics(
    std::filesystem::path output_dir, double permittivity)
    : output_dir_(std::move(output_dir)),
      permittivity_(permittivity) {
    if (!std::isfinite(permittivity_) || !(permittivity_ > 0.0)) {
        throw std::invalid_argument(
            "diagnostic permittivity must be positive and finite");
    }
    std::filesystem::create_directories(output_dir_);
    scalar_file_.open(output_dir_ / "scalars.csv");
    if (!scalar_file_) throw std::runtime_error("cannot open diagnostics output");
}
void Diagnostics::write_header() { scalar_file_ << "step,time,kinetic_energy,field_energy,total_energy,charge_l1,live_particles\n"; }
DiagnosticSample Diagnostics::sample(std::size_t step, double time, const Grid& grid, const std::vector<Species>& species) {
    DiagnosticSample s; s.step = step; s.time = time;
    for (const auto& sp : species) { s.kinetic_energy += sp.kinetic_energy(); s.live_particles += sp.live_count(); }
    for (std::size_t i = 0; i < grid.nx(); ++i) {
        const double volume = grid.node_volume(i);
        s.field_energy +=
            0.5 * permittivity_ *
            grid.electric()[i] * grid.electric()[i] * volume;
        s.charge_l1 += std::abs(grid.rho()[i]) * volume;
    }
    s.total_energy = s.kinetic_energy + s.field_energy;
    history_.push_back(s);
    return s;
}
void Diagnostics::write_sample(const DiagnosticSample& s) {
    scalar_file_ << s.step << ',' << std::setprecision(17) << s.time << ',' << s.kinetic_energy << ',' << s.field_energy << ',' << s.total_energy << ',' << s.charge_l1 << ',' << s.live_particles << '\n';
    scalar_file_.flush();
}
void Diagnostics::write_fields(std::size_t step, const Grid& grid) const {
    std::ofstream out(output_dir_ / ("fields_" + std::to_string(step) + ".csv"));
    out << "x,rho,phi,E\n";
    for (std::size_t i = 0; i < grid.nx(); ++i) out << grid.node_x(i) << ',' << grid.rho()[i] << ',' << grid.phi()[i] << ',' << grid.electric()[i] << '\n';
}

Diagnostics2D::Diagnostics2D(
    std::filesystem::path output_dir,
    const std::vector<Species2D>& species,
    double permittivity)
    : output_dir_(std::move(output_dir)),
      permittivity_(permittivity) {
    if (!std::isfinite(permittivity_) || !(permittivity_ > 0.0)) {
        throw std::invalid_argument(
            "2D diagnostic permittivity must be positive and finite");
    }
    std::filesystem::create_directories(output_dir_);
    scalar_file_.open(output_dir_ / "scalars.csv");
    if (!scalar_file_) throw std::runtime_error("cannot open 2D diagnostics output");
    species_names_.reserve(species.size());
    for (const auto& sp : species) species_names_.push_back(sp.name());
}

void Diagnostics2D::write_header() {
    scalar_file_ << "step,time,kinetic_energy,field_energy,total_energy,charge_l1,live_particles"
                 << ",absorbed_left,absorbed_right,absorbed_bottom,absorbed_top";
    for (const auto& name : species_names_) scalar_file_ << ",live_particles_" << name;
    scalar_file_ << '\n';
}

DiagnosticSample2D Diagnostics2D::sample(std::size_t step,
                                         double time,
                                         const Mesh2D& mesh,
                                         const std::vector<Species2D>& species,
                                         BoundaryLoss2D boundary_losses) {
    DiagnosticSample2D s;
    s.step = step;
    s.time = time;
    s.boundary_losses = boundary_losses;
    s.live_particles_by_species.reserve(species.size());
    for (const auto& sp : species) {
        s.kinetic_energy += sp.kinetic_energy();
        const auto live = sp.live_count();
        s.live_particles += live;
        s.live_particles_by_species.push_back(live);
    }
    for (std::size_t j = 0; j < mesh.ny(); ++j) {
        for (std::size_t i = 0; i < mesh.nx(); ++i) {
            const auto idx = mesh.index(i, j);
            const double e2 = mesh.electric_x()[idx] * mesh.electric_x()[idx] + mesh.electric_y()[idx] * mesh.electric_y()[idx];
            const double area = mesh.node_area(i, j);
            s.field_energy += 0.5 * permittivity_ * e2 * area;
            s.charge_l1 += std::abs(mesh.rho()[idx]) * area;
        }
    }
    s.total_energy = s.kinetic_energy + s.field_energy;
    history_.push_back(s);
    return s;
}

void Diagnostics2D::write_sample(const DiagnosticSample2D& s) {
    scalar_file_ << s.step << ',' << std::setprecision(17) << s.time << ',' << s.kinetic_energy << ','
                 << s.field_energy << ',' << s.total_energy << ',' << s.charge_l1 << ',' << s.live_particles << ','
                 << s.boundary_losses.absorbed_left << ',' << s.boundary_losses.absorbed_right << ','
                 << s.boundary_losses.absorbed_bottom << ',' << s.boundary_losses.absorbed_top;
    for (const auto live : s.live_particles_by_species) scalar_file_ << ',' << live;
    scalar_file_ << '\n';
    scalar_file_.flush();
}

void Diagnostics2D::write_particle_sample(std::size_t step,
                                          const std::vector<Species2D>& species,
                                          std::size_t stride,
                                          std::size_t sample_count) const {
    if (stride == 0) throw std::invalid_argument("2D particle output stride must be positive");
    std::ofstream out(output_dir_ / ("particles_" + std::to_string(step) + ".csv"));
    if (!out) throw std::runtime_error("cannot open 2D particle output");
    out << "species_id,species,x,y,vx,vy,vz,alive\n";
    std::size_t written = 0;
    for (std::size_t species_id = 0; species_id < species.size(); ++species_id) {
        const auto& sp = species[species_id];
        const auto& particles = sp.particles();
        for (std::size_t particle_id = 0; particle_id < particles.size(); ++particle_id) {
            if (particle_id % stride != 0) continue;
            const auto& p = particles[particle_id];
            out << species_id << ',' << sp.name() << ',' << std::setprecision(17)
                << p.position.x << ',' << p.position.y << ','
                << p.velocity.x << ',' << p.velocity.y << ','
                << p.velocity_z << ','
                << (p.alive ? 1 : 0) << '\n';
            ++written;
            if (sample_count != 0 && written >= sample_count) return;
        }
    }
}

Diagnostics3D::Diagnostics3D(
    std::filesystem::path output_dir,
    const std::vector<Species3D>& species,
    double permittivity)
    : output_dir_(std::move(output_dir)),
      permittivity_(permittivity) {
    if (!std::isfinite(permittivity_) || !(permittivity_ > 0.0)) {
        throw std::invalid_argument(
            "3D diagnostic permittivity must be positive and finite");
    }
    std::filesystem::create_directories(output_dir_);
    scalar_file_.open(output_dir_ / "scalars.csv");
    if (!scalar_file_) throw std::runtime_error("cannot open 3D diagnostics output");
    species_names_.reserve(species.size());
    for (const auto& sp : species) species_names_.push_back(sp.name());
}

void Diagnostics3D::write_header() {
    scalar_file_ << "step,time,kinetic_energy,field_energy,total_energy,charge_l1,live_particles"
                 << ",absorbed_left,absorbed_right,absorbed_bottom,absorbed_top,absorbed_back,absorbed_front";
    for (const auto& name : species_names_) scalar_file_ << ",live_particles_" << name;
    scalar_file_ << '\n';
}

DiagnosticSample3D Diagnostics3D::sample(std::size_t step,
                                         double time,
                                         const Mesh3D& mesh,
                                         const std::vector<Species3D>& species,
                                         BoundaryLoss3D boundary_losses) {
    DiagnosticSample3D s;
    s.step = step;
    s.time = time;
    s.boundary_losses = boundary_losses;
    s.live_particles_by_species.reserve(species.size());
    for (const auto& sp : species) {
        s.kinetic_energy += sp.kinetic_energy();
        const auto live = sp.live_count();
        s.live_particles += live;
        s.live_particles_by_species.push_back(live);
    }
    for (std::size_t k = 0; k < mesh.nz(); ++k) {
        for (std::size_t j = 0; j < mesh.ny(); ++j) {
            for (std::size_t i = 0; i < mesh.nx(); ++i) {
                const auto idx = mesh.index(i, j, k);
                const double e2 = mesh.electric_x()[idx] * mesh.electric_x()[idx]
                                + mesh.electric_y()[idx] * mesh.electric_y()[idx]
                                + mesh.electric_z()[idx] * mesh.electric_z()[idx];
                const double volume = mesh.node_volume(i, j, k);
                s.field_energy +=
                    0.5 * permittivity_ * e2 * volume;
                s.charge_l1 += std::abs(mesh.rho()[idx]) * volume;
            }
        }
    }
    s.total_energy = s.kinetic_energy + s.field_energy;
    history_.push_back(s);
    return s;
}

void Diagnostics3D::write_sample(const DiagnosticSample3D& s) {
    scalar_file_ << s.step << ',' << std::setprecision(17) << s.time << ',' << s.kinetic_energy << ','
                 << s.field_energy << ',' << s.total_energy << ',' << s.charge_l1 << ',' << s.live_particles << ','
                 << s.boundary_losses.absorbed_left << ',' << s.boundary_losses.absorbed_right << ','
                 << s.boundary_losses.absorbed_bottom << ',' << s.boundary_losses.absorbed_top << ','
                 << s.boundary_losses.absorbed_back << ',' << s.boundary_losses.absorbed_front;
    for (const auto live : s.live_particles_by_species) scalar_file_ << ',' << live;
    scalar_file_ << '\n';
    scalar_file_.flush();
}

void Diagnostics3D::write_particle_sample(std::size_t step,
                                          const std::vector<Species3D>& species,
                                          std::size_t stride,
                                          std::size_t sample_count) const {
    if (stride == 0) throw std::invalid_argument("3D particle output stride must be positive");
    std::ofstream out(output_dir_ / ("particles_" + std::to_string(step) + ".csv"));
    if (!out) throw std::runtime_error("cannot open 3D particle output");
    out << "species_id,species,x,y,z,vx,vy,vz,alive\n";
    std::size_t written = 0;
    for (std::size_t species_id = 0; species_id < species.size(); ++species_id) {
        const auto& sp = species[species_id];
        const auto& particles = sp.particles();
        for (std::size_t particle_id = 0; particle_id < particles.size(); ++particle_id) {
            if (particle_id % stride != 0) continue;
            const auto& p = particles[particle_id];
            out << species_id << ',' << sp.name() << ',' << std::setprecision(17)
                << p.position.x << ',' << p.position.y << ',' << p.position.z << ','
                << p.velocity.x << ',' << p.velocity.y << ',' << p.velocity.z << ','
                << (p.alive ? 1 : 0) << '\n';
            ++written;
            if (sample_count != 0 && written >= sample_count) return;
        }
    }
}
}
