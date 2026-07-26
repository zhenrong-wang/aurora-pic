#pragma once
#include "pic/Grid.hpp"
#include "pic/Mesh2D.hpp"

namespace pic {
class FieldSolver {
public:
    void solve(Grid& grid, double phi_left = 0.0, double phi_right = 0.0) const;
    void solve(Mesh2D& mesh) const;
private:
    void solve_periodic_spectral(Grid& grid) const;
    void solve_periodic_spectral(Mesh2D& mesh) const;
    void solve_dirichlet_tridiagonal(Grid& grid, double phi_left, double phi_right) const;
    void solve_dirichlet_iterative(Mesh2D& mesh) const;
};

double interpolate_electric(const Grid& grid, double x);
Vec2 interpolate_electric(const Mesh2D& mesh, Vec2 position);
}
