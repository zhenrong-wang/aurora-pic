#pragma once

#include "pic/UnstructuredMesh2D.hpp"

#include <cstddef>
#include <map>
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

UnstructuredPoissonSummary2D solve_unstructured_poisson(
    UnstructuredMesh2D& mesh,
    const std::map<std::string, double>& dirichlet_potentials,
    UnstructuredPoissonOptions2D options = {});

} // namespace pic
