#include "pic/Mesh2D.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
namespace pic {
namespace {
std::size_t checked_extent(std::size_t n, const char* name) {
    if (n < 3) throw std::invalid_argument(std::string(name) + " must be at least 3");
    return n;
}

double checked_length(double length, const char* name) {
    if (length <= 0.0) throw std::invalid_argument(std::string(name) + " must be positive");
    return length;
}

double wrap_periodic(double value, double length) {
    return std::fmod(std::fmod(value, length) + length, length);
}
}

Mesh2D::Mesh2D(std::size_t nx, std::size_t ny, double length_x, double length_y,
               Boundary boundary, BoundaryConfig2D boundary_config)
    : nx_(checked_extent(nx, "mesh nx")), ny_(checked_extent(ny, "mesh ny")),
      length_x_(checked_length(length_x, "domain length_x")),
      length_y_(checked_length(length_y, "domain length_y")),
      dx_(length_x_ / static_cast<double>(boundary == Boundary::Periodic ? nx_ : nx_ - 1)),
      dy_(length_y_ / static_cast<double>(boundary == Boundary::Periodic ? ny_ : ny_ - 1)),
      boundary_(boundary), boundary_config_(std::move(boundary_config)),
      rho_(nx_ * ny_, 0.0), phi_(nx_ * ny_, 0.0), electric_x_(nx_ * ny_, 0.0), electric_y_(nx_ * ny_, 0.0) {}

void Mesh2D::clear_charge() { std::fill(rho_.begin(), rho_.end(), 0.0); }

void deposit_charge_cic(Mesh2D& mesh, const std::vector<Particle2D>& particles, double charge, double weight) {
    auto& rho = mesh.rho();
    const double q_over_area = charge * weight / (mesh.dx() * mesh.dy());
    for (const auto& particle : particles) {
        if (!particle.alive) continue;

        const double x = mesh.boundary() == Boundary::Periodic
                             ? wrap_periodic(particle.position.x, mesh.length_x())
                             : std::clamp(particle.position.x, 0.0, mesh.length_x());
        const double y = mesh.boundary() == Boundary::Periodic
                             ? wrap_periodic(particle.position.y, mesh.length_y())
                             : std::clamp(particle.position.y, 0.0, mesh.length_y());
        const double gx = x / mesh.dx();
        const double gy = y / mesh.dy();

        std::size_t i = static_cast<std::size_t>(std::floor(gx));
        std::size_t j = static_cast<std::size_t>(std::floor(gy));
        double fx = gx - static_cast<double>(i);
        double fy = gy - static_cast<double>(j);

        std::size_t i0, i1, j0, j1;
        if (mesh.boundary() == Boundary::Periodic) {
            i0 = i % mesh.nx();
            i1 = (i + 1) % mesh.nx();
            j0 = j % mesh.ny();
            j1 = (j + 1) % mesh.ny();
        } else {
            i = std::min(i, mesh.nx() - 2);
            j = std::min(j, mesh.ny() - 2);
            fx = std::clamp(gx - static_cast<double>(i), 0.0, 1.0);
            fy = std::clamp(gy - static_cast<double>(j), 0.0, 1.0);
            i0 = i;
            i1 = i + 1;
            j0 = j;
            j1 = j + 1;
        }

        rho[mesh.index(i0, j0)] += q_over_area * (1.0 - fx) * (1.0 - fy);
        rho[mesh.index(i1, j0)] += q_over_area * fx * (1.0 - fy);
        rho[mesh.index(i0, j1)] += q_over_area * (1.0 - fx) * fy;
        rho[mesh.index(i1, j1)] += q_over_area * fx * fy;
    }
}
}
