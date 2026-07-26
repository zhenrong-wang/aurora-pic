#include "pic/Grid.hpp"
#include <algorithm>
#include <stdexcept>

namespace pic {
namespace {
std::size_t checked_nx(std::size_t nx) {
    if (nx < 3) throw std::invalid_argument("grid nx must be at least 3");
    return nx;
}

double checked_length(double length) {
    if (length <= 0.0) throw std::invalid_argument("domain length must be positive");
    return length;
}
}

Grid::Grid(std::size_t nx, double length, Boundary boundary)
    : nx_(checked_nx(nx)), length_(checked_length(length)),
      dx_(length_ / static_cast<double>(boundary == Boundary::Periodic ? nx_ : nx_ - 1)), boundary_(boundary),
      rho_(nx_, 0.0), phi_(nx_, 0.0), electric_(nx_, 0.0) {}

void Grid::clear_charge() { std::fill(rho_.begin(), rho_.end(), 0.0); }

double Grid::node_volume(std::size_t i) const {
    if (i >= nx_) throw std::out_of_range("grid node index out of range");
    if (boundary_ == Boundary::Periodic) return dx_;
    return (i == 0 || i + 1 == nx_) ? 0.5 * dx_ : dx_;
}
}
