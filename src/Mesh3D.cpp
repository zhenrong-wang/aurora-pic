#include "pic/Mesh3D.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

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

std::size_t checked_node_count(std::size_t nx, std::size_t ny, std::size_t nz) {
    if (nx > std::numeric_limits<std::size_t>::max() / ny) {
        throw std::invalid_argument("3D mesh node count exceeds supported size");
    }
    const std::size_t plane = nx * ny;
    if (plane > std::numeric_limits<std::size_t>::max() / nz) {
        throw std::invalid_argument("3D mesh node count exceeds supported size");
    }
    return plane * nz;
}

double wrap_periodic(double value, double length) {
    return std::fmod(std::fmod(value, length) + length, length);
}
}

Mesh3D::Mesh3D(std::size_t nx, std::size_t ny, std::size_t nz,
               double length_x, double length_y, double length_z, Boundary boundary)
    : nx_(checked_extent(nx, "mesh nx")), ny_(checked_extent(ny, "mesh ny")), nz_(checked_extent(nz, "mesh nz")),
      length_x_(checked_length(length_x, "domain length_x")),
      length_y_(checked_length(length_y, "domain length_y")),
      length_z_(checked_length(length_z, "domain length_z")),
      dx_(length_x_ / static_cast<double>(boundary == Boundary::Periodic ? nx_ : nx_ - 1)),
      dy_(length_y_ / static_cast<double>(boundary == Boundary::Periodic ? ny_ : ny_ - 1)),
      dz_(length_z_ / static_cast<double>(boundary == Boundary::Periodic ? nz_ : nz_ - 1)),
      boundary_(boundary),
      rho_(checked_node_count(nx_, ny_, nz_), 0.0), phi_(rho_.size(), 0.0),
      electric_x_(rho_.size(), 0.0), electric_y_(rho_.size(), 0.0), electric_z_(rho_.size(), 0.0) {}

void Mesh3D::clear_charge() { std::fill(rho_.begin(), rho_.end(), 0.0); }

double Mesh3D::node_volume(std::size_t i, std::size_t j, std::size_t k) const {
    if (i >= nx_ || j >= ny_ || k >= nz_) throw std::out_of_range("3D mesh node index out of range");
    if (boundary_ == Boundary::Periodic) return dx_ * dy_ * dz_;

    const double wx = (i == 0 || i + 1 == nx_) ? 0.5 * dx_ : dx_;
    const double wy = (j == 0 || j + 1 == ny_) ? 0.5 * dy_ : dy_;
    const double wz = (k == 0 || k + 1 == nz_) ? 0.5 * dz_ : dz_;
    return wx * wy * wz;
}

void deposit_charge_cic(Mesh3D& mesh, const std::vector<Particle3D>& particles, double charge, double weight) {
    auto& rho = mesh.rho();
    const double q = charge * weight;
    const double q_over_volume = q / (mesh.dx() * mesh.dy() * mesh.dz());

    for (const auto& particle : particles) {
        if (!particle.alive) continue;

        const double x = mesh.boundary() == Boundary::Periodic
                             ? wrap_periodic(particle.position.x, mesh.length_x())
                             : std::clamp(particle.position.x, 0.0, mesh.length_x());
        const double y = mesh.boundary() == Boundary::Periodic
                             ? wrap_periodic(particle.position.y, mesh.length_y())
                             : std::clamp(particle.position.y, 0.0, mesh.length_y());
        const double z = mesh.boundary() == Boundary::Periodic
                             ? wrap_periodic(particle.position.z, mesh.length_z())
                             : std::clamp(particle.position.z, 0.0, mesh.length_z());
        const double gx = x / mesh.dx();
        const double gy = y / mesh.dy();
        const double gz = z / mesh.dz();

        std::size_t i = static_cast<std::size_t>(std::floor(gx));
        std::size_t j = static_cast<std::size_t>(std::floor(gy));
        std::size_t k = static_cast<std::size_t>(std::floor(gz));
        double fx = gx - static_cast<double>(i);
        double fy = gy - static_cast<double>(j);
        double fz = gz - static_cast<double>(k);

        std::array<std::size_t, 2> ii{};
        std::array<std::size_t, 2> jj{};
        std::array<std::size_t, 2> kk{};
        if (mesh.boundary() == Boundary::Periodic) {
            ii = {i % mesh.nx(), (i + 1) % mesh.nx()};
            jj = {j % mesh.ny(), (j + 1) % mesh.ny()};
            kk = {k % mesh.nz(), (k + 1) % mesh.nz()};
        } else {
            i = std::min(i, mesh.nx() - 2);
            j = std::min(j, mesh.ny() - 2);
            k = std::min(k, mesh.nz() - 2);
            fx = std::clamp(gx - static_cast<double>(i), 0.0, 1.0);
            fy = std::clamp(gy - static_cast<double>(j), 0.0, 1.0);
            fz = std::clamp(gz - static_cast<double>(k), 0.0, 1.0);
            ii = {i, i + 1};
            jj = {j, j + 1};
            kk = {k, k + 1};
        }

        const std::array<double, 2> wx{1.0 - fx, fx};
        const std::array<double, 2> wy{1.0 - fy, fy};
        const std::array<double, 2> wz{1.0 - fz, fz};
        for (std::size_t dz_i = 0; dz_i < 2; ++dz_i) {
            for (std::size_t dy_i = 0; dy_i < 2; ++dy_i) {
                for (std::size_t dx_i = 0; dx_i < 2; ++dx_i) {
                    const double shape = wx[dx_i] * wy[dy_i] * wz[dz_i];
                    const auto idx = mesh.index(ii[dx_i], jj[dy_i], kk[dz_i]);
                    if (mesh.boundary() == Boundary::Periodic) {
                        rho[idx] += q_over_volume * shape;
                    } else {
                        rho[idx] += q * shape / mesh.node_volume(ii[dx_i], jj[dy_i], kk[dz_i]);
                    }
                }
            }
        }
    }
}
}
