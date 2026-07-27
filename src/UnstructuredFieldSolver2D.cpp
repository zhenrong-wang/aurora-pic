#include "pic/UnstructuredFieldSolver2D.hpp"

#include "pic/Types.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <vector>

namespace pic {
namespace {

using SparseRows = std::vector<std::map<std::size_t, double>>;

struct CompressedSparseMatrix {
    std::vector<std::size_t> row_offsets;
    std::vector<std::size_t> columns;
    std::vector<double> coefficients;
    std::vector<double> diagonal;
};

struct ElementQuadraturePoint2D {
    std::vector<double> weights;
    std::vector<Vec2> gradients;
    double area_weight{0.0};
};

double dot(const std::vector<double>& first, const std::vector<double>& second) {
    if (first.size() != second.size()) throw std::invalid_argument("vector dot-product size mismatch");
    double result = 0.0;
    for (std::size_t i = 0; i < first.size(); ++i) result += first[i] * second[i];
    return result;
}

CompressedSparseMatrix compress(const SparseRows& rows) {
    CompressedSparseMatrix result;
    result.row_offsets.reserve(rows.size() + 1);
    result.diagonal.assign(rows.size(), 0.0);
    result.row_offsets.push_back(0);
    for (std::size_t row = 0; row < rows.size(); ++row) {
        for (const auto& [column, coefficient] : rows[row]) {
            if (coefficient == 0.0) continue;
            if (column >= rows.size() || !std::isfinite(coefficient)) {
                throw std::runtime_error("invalid sparse matrix coefficient");
            }
            result.columns.push_back(column);
            result.coefficients.push_back(coefficient);
            if (column == row) result.diagonal[row] = coefficient;
        }
        result.row_offsets.push_back(result.columns.size());
        if (!(result.diagonal[row] > 0.0) || !std::isfinite(result.diagonal[row])) {
            throw std::runtime_error("unstructured Poisson matrix has a non-positive diagonal");
        }
    }
    return result;
}

std::vector<double> multiply(const CompressedSparseMatrix& matrix,
                             const std::vector<double>& values) {
    if (matrix.diagonal.size() != values.size() ||
        matrix.row_offsets.size() != values.size() + 1) {
        throw std::invalid_argument("sparse matrix-vector size mismatch");
    }
    std::vector<double> result(values.size(), 0.0);
    for (std::size_t row = 0; row < values.size(); ++row) {
        for (std::size_t entry = matrix.row_offsets[row];
             entry < matrix.row_offsets[row + 1]; ++entry) {
            const std::size_t column = matrix.columns[entry];
            if (column >= values.size()) throw std::runtime_error("sparse matrix column out of range");
            result[row] += matrix.coefficients[entry] * values[column];
        }
        if (!std::isfinite(result[row])) throw std::overflow_error("sparse matrix-vector product overflow");
    }
    return result;
}

std::vector<ElementQuadraturePoint2D> triangle_quadrature(
    const ImportedMesh2D& topology, const ImportedCell2D& cell) {
    const Vec2 first = topology.node_by_id(cell.node_ids[0]).position;
    const Vec2 second = topology.node_by_id(cell.node_ids[1]).position;
    const Vec2 third = topology.node_by_id(cell.node_ids[2]).position;
    const double determinant =
        (second.x - first.x) * (third.y - first.y) -
        (second.y - first.y) * (third.x - first.x);
    if (determinant == 0.0 || !std::isfinite(determinant)) {
        throw std::runtime_error("degenerate triangle in unstructured Poisson assembly");
    }
    const double inverse = 1.0 / determinant;
    std::vector<Vec2> gradients{
        {(second.y - third.y) * inverse, (third.x - second.x) * inverse},
        {(third.y - first.y) * inverse, (first.x - third.x) * inverse},
        {(first.y - second.y) * inverse, (second.x - first.x) * inverse},
    };
    return {{{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
             std::move(gradients), 0.5 * std::abs(determinant)}};
}

std::vector<ElementQuadraturePoint2D> quadrilateral_quadrature(
    const ImportedMesh2D& topology, const ImportedCell2D& cell) {
    std::array<Vec2, 4> nodes{};
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        nodes[i] = topology.node_by_id(cell.node_ids[i]).position;
    }

    const double point = 1.0 / std::sqrt(3.0);
    std::vector<ElementQuadraturePoint2D> result;
    result.reserve(4);
    for (const double xi : {-point, point}) {
        for (const double eta : {-point, point}) {
            const std::array<double, 4> weights{
                0.25 * (1.0 - xi) * (1.0 - eta),
                0.25 * (1.0 + xi) * (1.0 - eta),
                0.25 * (1.0 + xi) * (1.0 + eta),
                0.25 * (1.0 - xi) * (1.0 + eta),
            };
            const std::array<Vec2, 4> reference_gradients{
                Vec2{-0.25 * (1.0 - eta), -0.25 * (1.0 - xi)},
                Vec2{0.25 * (1.0 - eta), -0.25 * (1.0 + xi)},
                Vec2{0.25 * (1.0 + eta), 0.25 * (1.0 + xi)},
                Vec2{-0.25 * (1.0 + eta), 0.25 * (1.0 - xi)},
            };
            double dx_dxi = 0.0;
            double dx_deta = 0.0;
            double dy_dxi = 0.0;
            double dy_deta = 0.0;
            for (std::size_t i = 0; i < nodes.size(); ++i) {
                dx_dxi += reference_gradients[i].x * nodes[i].x;
                dx_deta += reference_gradients[i].y * nodes[i].x;
                dy_dxi += reference_gradients[i].x * nodes[i].y;
                dy_deta += reference_gradients[i].y * nodes[i].y;
            }
            const double determinant = dx_dxi * dy_deta - dx_deta * dy_dxi;
            if (determinant == 0.0 || !std::isfinite(determinant)) {
                throw std::runtime_error("degenerate quadrilateral in unstructured Poisson assembly");
            }
            std::vector<Vec2> physical_gradients(4);
            for (std::size_t i = 0; i < physical_gradients.size(); ++i) {
                physical_gradients[i] = {
                    (dy_deta * reference_gradients[i].x -
                     dy_dxi * reference_gradients[i].y) / determinant,
                    (-dx_deta * reference_gradients[i].x +
                     dx_dxi * reference_gradients[i].y) / determinant,
                };
            }
            result.push_back({
                {weights.begin(), weights.end()},
                std::move(physical_gradients),
                std::abs(determinant),
            });
        }
    }
    return result;
}

std::vector<ElementQuadraturePoint2D> element_quadrature(
    const ImportedMesh2D& topology, const ImportedCell2D& cell) {
    return cell.shape == ImportedCellShape2D::Triangle
               ? triangle_quadrature(topology, cell)
               : quadrilateral_quadrature(topology, cell);
}

void validate_options(const UnstructuredPoissonOptions2D& options) {
    if (!std::isfinite(options.relative_tolerance) || options.relative_tolerance < 0.0) {
        throw std::invalid_argument("unstructured Poisson relative tolerance must be finite and non-negative");
    }
    if (!std::isfinite(options.absolute_tolerance) || options.absolute_tolerance < 0.0) {
        throw std::invalid_argument("unstructured Poisson absolute tolerance must be finite and non-negative");
    }
    if (options.relative_tolerance == 0.0 && options.absolute_tolerance == 0.0) {
        throw std::invalid_argument("unstructured Poisson requires a positive convergence tolerance");
    }
}

std::vector<std::optional<double>> map_dirichlet_nodes(
    const UnstructuredMesh2D& mesh,
    const std::map<std::string, double>& dirichlet_potentials) {
    if (dirichlet_potentials.empty()) {
        throw std::invalid_argument("unstructured Poisson requires tagged Dirichlet potentials");
    }
    std::set<std::string> boundary_labels;
    for (const auto& face : mesh.topology().boundary_faces()) boundary_labels.insert(face.label);
    for (const auto& [label, potential] : dirichlet_potentials) {
        if (label.empty()) throw std::invalid_argument("Dirichlet boundary label must not be empty");
        if (!std::isfinite(potential)) {
            throw std::invalid_argument("Dirichlet boundary potential must be finite");
        }
        if (!boundary_labels.contains(label)) {
            throw std::invalid_argument("Dirichlet boundary label not found in imported mesh: " + label);
        }
    }
    for (const auto& label : boundary_labels) {
        if (!dirichlet_potentials.contains(label)) {
            throw std::invalid_argument("missing Dirichlet potential for imported boundary label: " + label);
        }
    }

    std::vector<std::optional<double>> fixed(mesh.size());
    for (const auto& face : mesh.topology().boundary_faces()) {
        const double potential = dirichlet_potentials.at(face.label);
        for (const auto node_id : face.node_ids) {
            const std::size_t index = mesh.node_index(node_id);
            if (fixed[index]) {
                const double scale = std::max({1.0, std::abs(*fixed[index]), std::abs(potential)});
                if (std::abs(*fixed[index] - potential) >
                    128.0 * std::numeric_limits<double>::epsilon() * scale) {
                    throw std::invalid_argument(
                        "conflicting Dirichlet potentials meet at imported node " +
                        std::to_string(node_id));
                }
            } else {
                fixed[index] = potential;
            }
        }
    }
    return fixed;
}

void recover_electric_field(UnstructuredMesh2D& mesh) {
    std::fill(mesh.electric().begin(), mesh.electric().end(), Vec2{});
    std::vector<double> projected_areas(mesh.size(), 0.0);
    for (const auto& cell : mesh.topology().cells()) {
        const auto quadrature = element_quadrature(mesh.topology(), cell);
        for (const auto& point : quadrature) {
            Vec2 electric{};
            for (std::size_t local = 0; local < cell.node_ids.size(); ++local) {
                const double potential = mesh.phi()[mesh.node_index(cell.node_ids[local])];
                electric.x -= potential * point.gradients[local].x;
                electric.y -= potential * point.gradients[local].y;
            }
            for (std::size_t local = 0; local < cell.node_ids.size(); ++local) {
                const std::size_t global = mesh.node_index(cell.node_ids[local]);
                const double projection_weight = point.weights[local] * point.area_weight;
                mesh.electric()[global].x += projection_weight * electric.x;
                mesh.electric()[global].y += projection_weight * electric.y;
                projected_areas[global] += projection_weight;
            }
        }
    }
    for (std::size_t i = 0; i < mesh.size(); ++i) {
        if (!(projected_areas[i] > 0.0) || !std::isfinite(projected_areas[i])) {
            throw std::runtime_error("invalid nodal area during unstructured electric-field recovery");
        }
        mesh.electric()[i].x /= projected_areas[i];
        mesh.electric()[i].y /= projected_areas[i];
        if (!std::isfinite(mesh.electric()[i].x) || !std::isfinite(mesh.electric()[i].y)) {
            throw std::overflow_error("unstructured electric-field recovery overflow");
        }
    }
}

} // namespace

UnstructuredPoissonSummary2D solve_unstructured_poisson(
    UnstructuredMesh2D& mesh,
    const std::map<std::string, double>& dirichlet_potentials,
    UnstructuredPoissonOptions2D options) {
    validate_options(options);
    if (mesh.rho().size() != mesh.size() || mesh.phi().size() != mesh.size() ||
        mesh.electric().size() != mesh.size() ||
        mesh.node_control_areas().size() != mesh.size()) {
        throw std::invalid_argument("unstructured mesh nodal field arrays have inconsistent sizes");
    }
    for (const double density : mesh.rho()) {
        if (!std::isfinite(density)) {
            throw std::invalid_argument("unstructured charge density must be finite");
        }
    }
    for (const double potential : mesh.phi()) {
        if (!std::isfinite(potential)) {
            throw std::invalid_argument("unstructured initial potential must be finite");
        }
    }

    SparseRows matrix(mesh.size());
    for (const auto& cell : mesh.topology().cells()) {
        const auto quadrature = element_quadrature(mesh.topology(), cell);
        for (const auto& point : quadrature) {
            for (std::size_t local_row = 0; local_row < cell.node_ids.size(); ++local_row) {
                const std::size_t row = mesh.node_index(cell.node_ids[local_row]);
                for (std::size_t local_column = 0;
                     local_column < cell.node_ids.size(); ++local_column) {
                    const std::size_t column = mesh.node_index(cell.node_ids[local_column]);
                    const double coefficient =
                        (point.gradients[local_row].x * point.gradients[local_column].x +
                         point.gradients[local_row].y * point.gradients[local_column].y) *
                        point.area_weight;
                    matrix[row][column] += coefficient;
                }
            }
        }
    }

    std::vector<double> right_hand_side(mesh.size(), 0.0);
    for (std::size_t i = 0; i < mesh.size(); ++i) {
        right_hand_side[i] = mesh.rho()[i] * mesh.node_control_areas()[i] / EPS0;
        if (!std::isfinite(right_hand_side[i])) {
            throw std::overflow_error("unstructured Poisson right-hand side overflow");
        }
    }

    const auto fixed = map_dirichlet_nodes(mesh, dirichlet_potentials);
    for (std::size_t row = 0; row < matrix.size(); ++row) {
        if (fixed[row]) continue;
        for (auto it = matrix[row].begin(); it != matrix[row].end();) {
            if (fixed[it->first]) {
                right_hand_side[row] -= it->second * *fixed[it->first];
                it = matrix[row].erase(it);
            } else {
                ++it;
            }
        }
    }
    for (std::size_t row = 0; row < matrix.size(); ++row) {
        if (!fixed[row]) continue;
        matrix[row].clear();
        matrix[row][row] = 1.0;
        right_hand_side[row] = *fixed[row];
        mesh.phi()[row] = *fixed[row];
    }
    for (const double value : right_hand_side) {
        if (!std::isfinite(value)) {
            throw std::overflow_error("unstructured Poisson constrained right-hand side overflow");
        }
    }
    const CompressedSparseMatrix compressed_matrix = compress(matrix);

    std::vector<double> matrix_times_solution = multiply(compressed_matrix, mesh.phi());
    std::vector<double> residual(mesh.size(), 0.0);
    for (std::size_t i = 0; i < mesh.size(); ++i) {
        residual[i] = right_hand_side[i] - matrix_times_solution[i];
    }
    const double residual_squared = dot(residual, residual);
    if (!std::isfinite(residual_squared)) {
        throw std::overflow_error("unstructured Poisson initial residual overflow");
    }
    UnstructuredPoissonSummary2D summary;
    summary.initial_residual = std::sqrt(residual_squared);
    summary.final_residual = summary.initial_residual;
    const double rhs_squared = dot(right_hand_side, right_hand_side);
    if (!std::isfinite(rhs_squared)) {
        throw std::overflow_error("unstructured Poisson right-hand-side norm overflow");
    }
    const double rhs_norm = std::sqrt(rhs_squared);
    const double target =
        std::max(options.absolute_tolerance, options.relative_tolerance * rhs_norm);
    if (summary.final_residual <= target) {
        summary.converged = true;
        recover_electric_field(mesh);
        return summary;
    }

    std::vector<double> preconditioned_residual(mesh.size(), 0.0);
    for (std::size_t i = 0; i < mesh.size(); ++i) {
        preconditioned_residual[i] = residual[i] / compressed_matrix.diagonal[i];
    }
    std::vector<double> direction = preconditioned_residual;
    double residual_preconditioned = dot(residual, preconditioned_residual);
    if (!(residual_preconditioned > 0.0) || !std::isfinite(residual_preconditioned)) {
        throw std::runtime_error("unstructured Poisson preconditioner is not positive definite");
    }

    std::size_t max_iterations = options.max_iterations;
    if (max_iterations == 0) {
        max_iterations = mesh.size() > std::numeric_limits<std::size_t>::max() / 10
                             ? std::numeric_limits<std::size_t>::max()
                             : std::max<std::size_t>(100, 10 * mesh.size());
    }
    for (std::size_t iteration = 0; iteration < max_iterations; ++iteration) {
        const std::vector<double> matrix_times_direction = multiply(compressed_matrix, direction);
        const double curvature = dot(direction, matrix_times_direction);
        if (!(curvature > 0.0) || !std::isfinite(curvature)) {
            throw std::runtime_error("unstructured Poisson matrix is not positive definite");
        }
        const double alpha = residual_preconditioned / curvature;
        if (!std::isfinite(alpha)) {
            throw std::overflow_error("unstructured Poisson step length overflow");
        }
        for (std::size_t i = 0; i < mesh.size(); ++i) {
            mesh.phi()[i] += alpha * direction[i];
            residual[i] -= alpha * matrix_times_direction[i];
            if (!std::isfinite(mesh.phi()[i]) || !std::isfinite(residual[i])) {
                throw std::overflow_error("unstructured Poisson iteration overflow");
            }
        }
        const double next_residual_squared = dot(residual, residual);
        if (!std::isfinite(next_residual_squared)) {
            throw std::overflow_error("unstructured Poisson residual overflow");
        }
        summary.iterations = iteration + 1;
        summary.final_residual = std::sqrt(next_residual_squared);
        if (summary.final_residual <= target) {
            summary.converged = true;
            recover_electric_field(mesh);
            return summary;
        }
        for (std::size_t i = 0; i < mesh.size(); ++i) {
            preconditioned_residual[i] = residual[i] / compressed_matrix.diagonal[i];
        }
        const double next_residual_preconditioned = dot(residual, preconditioned_residual);
        if (!(next_residual_preconditioned > 0.0) ||
            !std::isfinite(next_residual_preconditioned)) {
            throw std::runtime_error("unstructured Poisson preconditioned residual is invalid");
        }
        const double beta = next_residual_preconditioned / residual_preconditioned;
        for (std::size_t i = 0; i < direction.size(); ++i) {
            direction[i] = preconditioned_residual[i] + beta * direction[i];
        }
        residual_preconditioned = next_residual_preconditioned;
    }
    return summary;
}

} // namespace pic
