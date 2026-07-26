#pragma once
#include "pic/Types.hpp"
#include <cstddef>
#include <vector>

namespace pic {
class Grid {
public:
    Grid(std::size_t nx, double length, Boundary boundary);
    std::size_t nx() const { return nx_; }
    double length() const { return length_; }
    double dx() const { return dx_; }
    Boundary boundary() const { return boundary_; }
    std::vector<double>& rho() { return rho_; }
    std::vector<double>& phi() { return phi_; }
    std::vector<double>& electric() { return electric_; }
    const std::vector<double>& rho() const { return rho_; }
    const std::vector<double>& phi() const { return phi_; }
    const std::vector<double>& electric() const { return electric_; }
    void clear_charge();
    double node_x(std::size_t i) const { return static_cast<double>(i) * dx_; }
    double node_volume(std::size_t i) const;
private:
    std::size_t nx_;
    double length_;
    double dx_;
    Boundary boundary_;
    std::vector<double> rho_, phi_, electric_;
};
}
