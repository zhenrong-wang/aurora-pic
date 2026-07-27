#pragma once

#include "pic/UnstructuredMesh2D.hpp"

#include <cstddef>
#include <map>
#include <memory>
#include <string>

namespace pic {

struct UnstructuredPoissonOptions2D {
    double relative_tolerance{1e-10};
    double absolute_tolerance{1e-12};
    std::size_t max_iterations{0}; // zero selects a mesh-dependent bound
};

struct UnstructuredPoissonSummary2D {
    std::size_t iterations{0};
    double initial_residual{0.0};
    double final_residual{0.0};
    bool converged{false};
};

class UnstructuredPoissonSolver2D {
public:
    UnstructuredPoissonSolver2D(
        const UnstructuredMesh2D& mesh,
        std::map<std::string, double> dirichlet_potentials,
        UnstructuredPoissonOptions2D options = {});
    ~UnstructuredPoissonSolver2D();
    UnstructuredPoissonSolver2D(UnstructuredPoissonSolver2D&&) noexcept;
    UnstructuredPoissonSolver2D& operator=(UnstructuredPoissonSolver2D&&) noexcept;
    UnstructuredPoissonSolver2D(const UnstructuredPoissonSolver2D&) = delete;
    UnstructuredPoissonSolver2D& operator=(const UnstructuredPoissonSolver2D&) = delete;

    UnstructuredPoissonSummary2D solve(UnstructuredMesh2D& mesh) const;
    std::size_t assembly_count() const;
    std::size_t solve_count() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

UnstructuredPoissonSummary2D solve_unstructured_poisson(
    UnstructuredMesh2D& mesh,
    const std::map<std::string, double>& dirichlet_potentials,
    UnstructuredPoissonOptions2D options = {});

} // namespace pic
