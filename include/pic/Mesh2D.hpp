#pragma once
#include "pic/Types.hpp"
#include <cstddef>
#include <string>
#include <vector>

namespace pic {
struct BoundarySide2D {
    std::string tag{"wall"};
    double potential{0.0};
};

struct BoundaryConfig2D {
    BoundarySide2D left{"left", 0.0};
    BoundarySide2D right{"right", 0.0};
    BoundarySide2D bottom{"bottom", 0.0};
    BoundarySide2D top{"top", 0.0};
};

class Mesh2D {
public:
    Mesh2D(std::size_t nx, std::size_t ny, double length_x, double length_y,
           Boundary boundary, BoundaryConfig2D boundary_config = {});

    std::size_t nx() const { return nx_; }
    std::size_t ny() const { return ny_; }
    std::size_t size() const { return rho_.size(); }
    double length_x() const { return length_x_; }
    double length_y() const { return length_y_; }
    double dx() const { return dx_; }
    double dy() const { return dy_; }
    Boundary boundary() const { return boundary_; }
    const BoundaryConfig2D& boundary_config() const { return boundary_config_; }

    std::size_t index(std::size_t i, std::size_t j) const { return j * nx_ + i; }
    double node_x(std::size_t i) const { return static_cast<double>(i) * dx_; }
    double node_y(std::size_t j) const { return static_cast<double>(j) * dy_; }

    std::vector<double>& rho() { return rho_; }
    std::vector<double>& phi() { return phi_; }
    std::vector<double>& electric_x() { return electric_x_; }
    std::vector<double>& electric_y() { return electric_y_; }
    const std::vector<double>& rho() const { return rho_; }
    const std::vector<double>& phi() const { return phi_; }
    const std::vector<double>& electric_x() const { return electric_x_; }
    const std::vector<double>& electric_y() const { return electric_y_; }

    void clear_charge();

private:
    std::size_t nx_;
    std::size_t ny_;
    double length_x_;
    double length_y_;
    double dx_;
    double dy_;
    Boundary boundary_;
    BoundaryConfig2D boundary_config_;
    std::vector<double> rho_, phi_, electric_x_, electric_y_;
};

void deposit_charge_cic(Mesh2D& mesh, const std::vector<Particle2D>& particles, double charge, double weight);
}
