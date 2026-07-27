#include "pic/Mesh2D.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
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
    if (!std::isfinite(length) || length <= 0.0) throw std::invalid_argument(std::string(name) + " must be positive and finite");
    return length;
}

std::size_t checked_node_count(std::size_t nx, std::size_t ny) {
    if (nx > std::numeric_limits<std::size_t>::max() / ny) {
        throw std::invalid_argument("2D mesh node count exceeds supported size");
    }
    return nx * ny;
}

BoundaryConfig2D checked_boundary_config(BoundaryConfig2D boundary_config) {
    const auto check_side = [](const BoundarySide2D& side, const char* name) {
        if (!std::isfinite(side.potential)) {
            throw std::invalid_argument(std::string("2D boundary side '") + name + "' potential must be finite");
        }
    };
    check_side(boundary_config.left, "left");
    check_side(boundary_config.right, "right");
    check_side(boundary_config.bottom, "bottom");
    check_side(boundary_config.top, "top");
    return boundary_config;
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
      boundary_(boundary), boundary_config_(checked_boundary_config(std::move(boundary_config))),
      rho_(checked_node_count(nx_, ny_), 0.0), phi_(rho_.size(), 0.0),
      electric_x_(rho_.size(), 0.0), electric_y_(rho_.size(), 0.0) {}

void Mesh2D::clear_charge() { std::fill(rho_.begin(), rho_.end(), 0.0); }

double Mesh2D::node_area(std::size_t i, std::size_t j) const {
    if (i >= nx_ || j >= ny_) throw std::out_of_range("2D mesh node index out of range");
    if (boundary_ == Boundary::Periodic) return dx_ * dy_;

    const double wx = (i == 0 || i + 1 == nx_) ? 0.5 * dx_ : dx_;
    const double wy = (j == 0 || j + 1 == ny_) ? 0.5 * dy_ : dy_;
    return wx * wy;
}

void deposit_charge_cic(Mesh2D& mesh, const std::vector<Particle2D>& particles, double charge, double weight) {
    auto& rho = mesh.rho();
    const double q = charge * weight;
    const double q_over_area = q / (mesh.dx() * mesh.dy());
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

        if (mesh.boundary() == Boundary::Periodic) {
            rho[mesh.index(i0, j0)] += q_over_area * (1.0 - fx) * (1.0 - fy);
            rho[mesh.index(i1, j0)] += q_over_area * fx * (1.0 - fy);
            rho[mesh.index(i0, j1)] += q_over_area * (1.0 - fx) * fy;
            rho[mesh.index(i1, j1)] += q_over_area * fx * fy;
        } else {
            rho[mesh.index(i0, j0)] += q * (1.0 - fx) * (1.0 - fy) / mesh.node_area(i0, j0);
            rho[mesh.index(i1, j0)] += q * fx * (1.0 - fy) / mesh.node_area(i1, j0);
            rho[mesh.index(i0, j1)] += q * (1.0 - fx) * fy / mesh.node_area(i0, j1);
            rho[mesh.index(i1, j1)] += q * fx * fy / mesh.node_area(i1, j1);
        }
    }
}
}
