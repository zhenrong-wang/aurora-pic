#pragma once
#include "pic/Grid.hpp"
#include "pic/Mesh2D.hpp"
#include "pic/Mesh3D.hpp"

namespace pic {
class FieldSolver {
public:
    explicit FieldSolver(double permittivity = EPS0);
    double permittivity() const { return permittivity_; }
    void solve(Grid& grid, double phi_left = 0.0, double phi_right = 0.0) const;
    void solve(Mesh2D& mesh) const;
    void solve(Mesh3D& mesh) const;
private:
    double permittivity_{EPS0};
    void solve_periodic_spectral(Grid& grid) const;
    void solve_periodic_spectral(Mesh2D& mesh) const;
    void solve_periodic_spectral(Mesh3D& mesh) const;
    void solve_mixed_spectral_tridiagonal(Mesh2D& mesh) const;
    void solve_dirichlet_tridiagonal(Grid& grid, double phi_left, double phi_right) const;
    void solve_dirichlet_iterative(Mesh2D& mesh) const;
    void solve_dirichlet_iterative(Mesh3D& mesh) const;
};

double interpolate_electric(const Grid& grid, double x);
Vec2 interpolate_electric(const Mesh2D& mesh, Vec2 position);
Vec3 interpolate_electric(const Mesh3D& mesh, Vec3 position);
}
