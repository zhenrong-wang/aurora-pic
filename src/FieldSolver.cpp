#include "pic/FieldSolver.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace pic {
namespace {
using Complex = std::complex<double>;

double wrap_periodic(double value, double length) {
    return std::fmod(std::fmod(value, length) + length, length);
}

bool is_power_of_two(std::size_t n) {
    return n != 0 && (n & (n - 1)) == 0;
}

void dft_unscaled(std::vector<Complex>& values, bool inverse) {
    const std::size_t n = values.size();
    const double sign = inverse ? 1.0 : -1.0;
    const double twopi = 2.0 * std::numbers::pi;
    std::vector<Complex> out(n, Complex{0.0, 0.0});

    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            const double angle = sign * twopi * static_cast<double>(j * k) / static_cast<double>(n);
            out[k] += values[j] * Complex{std::cos(angle), std::sin(angle)};
        }
    }

    values.swap(out);
}

void fft_radix2(std::vector<Complex>& values, bool inverse);

void fft_bluestein_unscaled(
    std::vector<Complex>& values, bool inverse) {
    const std::size_t n = values.size();
    if (n >
        std::numeric_limits<std::size_t>::max() / 2U + 1U) {
        throw std::length_error(
            "spectral transform size exceeds supported range");
    }
    const std::size_t convolution_extent = 2 * n - 1;
    std::size_t padded_size = 1;
    while (padded_size < convolution_extent) {
        if (padded_size >
            std::numeric_limits<std::size_t>::max() / 2U) {
            throw std::length_error(
                "spectral convolution size exceeds supported range");
        }
        padded_size *= 2;
    }

    const double sign = inverse ? 1.0 : -1.0;
    const double pi_over_n =
        std::numbers::pi / static_cast<double>(n);
    std::vector<Complex> signal(
        padded_size, Complex{0.0, 0.0});
    std::vector<Complex> chirp(
        padded_size, Complex{0.0, 0.0});
    for (std::size_t index = 0; index < n; ++index) {
        const double square =
            static_cast<double>(index) *
            static_cast<double>(index);
        const double signal_angle =
            sign * pi_over_n * square;
        const double chirp_angle =
            -sign * pi_over_n * square;
        signal[index] =
            values[index] *
            Complex{std::cos(signal_angle),
                    std::sin(signal_angle)};
        const Complex chirp_value{
            std::cos(chirp_angle),
            std::sin(chirp_angle)};
        chirp[index] = chirp_value;
        if (index != 0) {
            chirp[padded_size - index] = chirp_value;
        }
    }

    fft_radix2(signal, false);
    fft_radix2(chirp, false);
    for (std::size_t index = 0;
         index < padded_size;
         ++index) {
        signal[index] *= chirp[index];
    }
    fft_radix2(signal, true);

    for (std::size_t frequency = 0;
         frequency < n;
         ++frequency) {
        const double square =
            static_cast<double>(frequency) *
            static_cast<double>(frequency);
        const double angle = sign * pi_over_n * square;
        values[frequency] =
            signal[frequency] *
            Complex{std::cos(angle), std::sin(angle)};
    }
}

std::size_t smallest_factor(std::size_t n) {
    for (std::size_t factor = 2;
         factor <= n / factor;
         ++factor) {
        if (n % factor == 0) return factor;
    }
    return n;
}

void fft_mixed_radix_unscaled(
    std::vector<Complex>& values, bool inverse) {
    const std::size_t n = values.size();
    if (n <= 1) return;
    const std::size_t radix = smallest_factor(n);
    if (radix == n) {
        if (n < 32) {
            dft_unscaled(values, inverse);
        } else {
            fft_bluestein_unscaled(values, inverse);
        }
        return;
    }

    const std::size_t subtransform_size = n / radix;
    std::vector<std::vector<Complex>> subtransforms(
        radix,
        std::vector<Complex>(subtransform_size));
    for (std::size_t residue = 0; residue < radix; ++residue) {
        for (std::size_t sample = 0;
             sample < subtransform_size;
             ++sample) {
            subtransforms[residue][sample] =
                values[residue + radix * sample];
        }
        fft_mixed_radix_unscaled(
            subtransforms[residue], inverse);
    }

    const double sign = inverse ? 1.0 : -1.0;
    const double twopi = 2.0 * std::numbers::pi;
    std::vector<Complex> output(n, Complex{0.0, 0.0});
    for (std::size_t frequency = 0;
         frequency < n;
         ++frequency) {
        const std::size_t subfrequency =
            frequency % subtransform_size;
        for (std::size_t residue = 0;
             residue < radix;
             ++residue) {
            const double angle =
                sign * twopi *
                static_cast<double>(residue * frequency) /
                static_cast<double>(n);
            output[frequency] +=
                subtransforms[residue][subfrequency] *
                Complex{std::cos(angle), std::sin(angle)};
        }
    }
    values.swap(output);
}

void fft_radix2(std::vector<Complex>& values, bool inverse) {
    const std::size_t n = values.size();

    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(values[i], values[j]);
    }

    const double twopi = 2.0 * std::numbers::pi;
    for (std::size_t len = 2; len <= n; len <<= 1) {
        const double angle = (inverse ? 1.0 : -1.0) * twopi / static_cast<double>(len);
        const Complex wlen{std::cos(angle), std::sin(angle)};
        for (std::size_t i = 0; i < n; i += len) {
            Complex w{1.0, 0.0};
            for (std::size_t j = 0; j < len / 2; ++j) {
                const Complex u = values[i + j];
                const Complex v = values[i + j + len / 2] * w;
                values[i + j] = u + v;
                values[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }

    if (inverse) {
        const double scale = 1.0 / static_cast<double>(n);
        for (auto& value : values) value *= scale;
    }
}

void transform_1d(std::vector<Complex>& values, bool inverse) {
    if (is_power_of_two(values.size())) {
        fft_radix2(values, inverse);
    } else {
        fft_mixed_radix_unscaled(values, inverse);
        if (inverse) {
            const double scale =
                1.0 / static_cast<double>(values.size());
            for (auto& value : values) value *= scale;
        }
    }
}

void transform_2d(std::vector<Complex>& values, std::size_t nx, std::size_t ny, bool inverse) {
    std::vector<Complex> line;
    line.reserve(std::max(nx, ny));

    line.resize(nx);
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) line[i] = values[j * nx + i];
        transform_1d(line, inverse);
        for (std::size_t i = 0; i < nx; ++i) values[j * nx + i] = line[i];
    }

    line.resize(ny);
    for (std::size_t i = 0; i < nx; ++i) {
        for (std::size_t j = 0; j < ny; ++j) line[j] = values[j * nx + i];
        transform_1d(line, inverse);
        for (std::size_t j = 0; j < ny; ++j) values[j * nx + i] = line[j];
    }
}

void transform_3d(std::vector<Complex>& values, std::size_t nx, std::size_t ny, std::size_t nz, bool inverse) {
    std::vector<Complex> line;
    line.reserve(std::max({nx, ny, nz}));

    line.resize(nx);
    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j < ny; ++j) {
            for (std::size_t i = 0; i < nx; ++i) line[i] = values[(k * ny + j) * nx + i];
            transform_1d(line, inverse);
            for (std::size_t i = 0; i < nx; ++i) values[(k * ny + j) * nx + i] = line[i];
        }
    }

    line.resize(ny);
    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t i = 0; i < nx; ++i) {
            for (std::size_t j = 0; j < ny; ++j) line[j] = values[(k * ny + j) * nx + i];
            transform_1d(line, inverse);
            for (std::size_t j = 0; j < ny; ++j) values[(k * ny + j) * nx + i] = line[j];
        }
    }

    line.resize(nz);
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            for (std::size_t k = 0; k < nz; ++k) line[k] = values[(k * ny + j) * nx + i];
            transform_1d(line, inverse);
            for (std::size_t k = 0; k < nz; ++k) values[(k * ny + j) * nx + i] = line[k];
        }
    }
}

void apply_dirichlet_boundary_potentials(Mesh2D& mesh) {
    auto& phi = mesh.phi();
    const auto& bc = mesh.boundary_config();
    const std::size_t nx = mesh.nx();
    const std::size_t ny = mesh.ny();

    if (mesh.boundary_x() == Boundary::Dirichlet) {
        for (std::size_t j = 0; j < ny; ++j) {
            phi[mesh.index(0, j)] = bc.left.potential;
            phi[mesh.index(nx - 1, j)] = bc.right.potential;
        }
    }
    if (mesh.boundary_y() == Boundary::Dirichlet) {
        for (std::size_t i = 0; i < nx; ++i) {
            phi[mesh.index(i, 0)] = bc.bottom.potential;
            phi[mesh.index(i, ny - 1)] = bc.top.potential;
        }
    }

    if (mesh.fully_dirichlet()) {
        phi[mesh.index(0, 0)] =
            0.5 * (bc.left.potential + bc.bottom.potential);
        phi[mesh.index(nx - 1, 0)] =
            0.5 * (bc.right.potential + bc.bottom.potential);
        phi[mesh.index(0, ny - 1)] =
            0.5 * (bc.left.potential + bc.top.potential);
        phi[mesh.index(nx - 1, ny - 1)] =
            0.5 * (bc.right.potential + bc.top.potential);
    }
}

void compute_electric_field(Mesh2D& mesh) {
    const std::size_t nx = mesh.nx();
    const std::size_t ny = mesh.ny();
    const auto& phi = mesh.phi();
    auto& ex = mesh.electric_x();
    auto& ey = mesh.electric_y();

    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            const auto idx = mesh.index(i, j);
            if (mesh.boundary_x() == Boundary::Periodic) {
                const std::size_t im = (i + nx - 1) % nx;
                const std::size_t ip = (i + 1) % nx;
                ex[idx] =
                    -(phi[mesh.index(ip, j)] -
                      phi[mesh.index(im, j)]) /
                    (2.0 * mesh.dx());
            } else if (i == 0) {
                ex[idx] =
                    -(phi[mesh.index(1, j)] - phi[idx]) / mesh.dx();
            } else if (i + 1 == nx) {
                ex[idx] =
                    -(phi[idx] - phi[mesh.index(i - 1, j)]) / mesh.dx();
            } else {
                ex[idx] =
                    -(phi[mesh.index(i + 1, j)] -
                      phi[mesh.index(i - 1, j)]) /
                    (2.0 * mesh.dx());
            }

            if (mesh.boundary_y() == Boundary::Periodic) {
                const std::size_t jm = (j + ny - 1) % ny;
                const std::size_t jp = (j + 1) % ny;
                ey[idx] =
                    -(phi[mesh.index(i, jp)] -
                      phi[mesh.index(i, jm)]) /
                    (2.0 * mesh.dy());
            } else if (j == 0) {
                ey[idx] =
                    -(phi[mesh.index(i, 1)] - phi[idx]) / mesh.dy();
            } else if (j + 1 == ny) {
                ey[idx] =
                    -(phi[idx] - phi[mesh.index(i, j - 1)]) / mesh.dy();
            } else {
                ey[idx] =
                    -(phi[mesh.index(i, j + 1)] -
                      phi[mesh.index(i, j - 1)]) /
                    (2.0 * mesh.dy());
            }
        }
    }
}

void apply_grounded_dirichlet_boundary(Mesh3D& mesh) {
    auto& phi = mesh.phi();
    const std::size_t nx = mesh.nx();
    const std::size_t ny = mesh.ny();
    const std::size_t nz = mesh.nz();

    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j < ny; ++j) {
            phi[mesh.index(0, j, k)] = 0.0;
            phi[mesh.index(nx - 1, j, k)] = 0.0;
        }
        for (std::size_t i = 0; i < nx; ++i) {
            phi[mesh.index(i, 0, k)] = 0.0;
            phi[mesh.index(i, ny - 1, k)] = 0.0;
        }
    }
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            phi[mesh.index(i, j, 0)] = 0.0;
            phi[mesh.index(i, j, nz - 1)] = 0.0;
        }
    }
}
}

FieldSolver::FieldSolver(double permittivity)
    : permittivity_(permittivity) {
    if (!std::isfinite(permittivity_) || !(permittivity_ > 0.0)) {
        throw std::invalid_argument(
            "field-solver permittivity must be positive and finite");
    }
}

void FieldSolver::solve(Grid& grid, double phi_left, double phi_right) const {
    if (grid.boundary() == Boundary::Periodic) solve_periodic_spectral(grid);
    else solve_dirichlet_tridiagonal(grid, phi_left, phi_right);
}

void FieldSolver::solve(Mesh2D& mesh) const {
    if (mesh.fully_periodic()) {
        solve_periodic_spectral(mesh);
    } else if (mesh.boundary_x() != mesh.boundary_y()) {
        solve_mixed_spectral_tridiagonal(mesh);
    } else {
        solve_dirichlet_iterative(mesh);
    }
}

void FieldSolver::solve(Mesh3D& mesh) const {
    if (mesh.boundary() == Boundary::Periodic) {
        solve_periodic_spectral(mesh);
    } else {
        solve_dirichlet_iterative(mesh);
    }
}

void FieldSolver::solve_periodic_spectral(Grid& grid) const {
    const std::size_t n = grid.nx();
    const double L = grid.length();
    const double dx = grid.dx();
    const double twopi = 2.0 * std::numbers::pi;
    auto& rho = grid.rho();
    auto& phi = grid.phi();
    auto& E = grid.electric();
    double mean = 0.0;
    for (double r : rho) mean += r;
    mean /= static_cast<double>(n);
    std::fill(phi.begin(), phi.end(), 0.0);
    std::fill(E.begin(), E.end(), 0.0);
    for (std::size_t k = 1; k < n; ++k) {
        double kk = (k <= n/2) ? static_cast<double>(k) : -static_cast<double>(n-k);
        double wave = twopi * kk / L;
        std::complex<double> rhok{0.0, 0.0};
        for (std::size_t j = 0; j < n; ++j) {
            double angle = -twopi * static_cast<double>(k*j) / static_cast<double>(n);
            rhok += (rho[j] - mean) * std::complex<double>(std::cos(angle), std::sin(angle));
        }
        rhok /= static_cast<double>(n);
        std::complex<double> phik =
            rhok / (permittivity_ * wave * wave);
        std::complex<double> Ek = -std::complex<double>(0.0, 1.0) * wave * phik;
        for (std::size_t j = 0; j < n; ++j) {
            double angle = twopi * static_cast<double>(k*j) / static_cast<double>(n);
            std::complex<double> basis(std::cos(angle), std::sin(angle));
            phi[j] += (phik * basis).real();
            E[j] += (Ek * basis).real();
        }
    }
    (void)dx;
}

void FieldSolver::solve_periodic_spectral(Mesh2D& mesh) const {
    const std::size_t nx = mesh.nx();
    const std::size_t ny = mesh.ny();
    const double Lx = mesh.length_x();
    const double Ly = mesh.length_y();
    const double twopi = 2.0 * std::numbers::pi;
    auto& rho = mesh.rho();
    auto& phi = mesh.phi();
    auto& Ex = mesh.electric_x();
    auto& Ey = mesh.electric_y();

    double mean = 0.0;
    for (double r : rho) mean += r;
    mean /= static_cast<double>(rho.size());

    std::vector<Complex> rho_hat(nx * ny, Complex{0.0, 0.0});
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            const auto idx = mesh.index(i, j);
            rho_hat[idx] = Complex{rho[idx] - mean, 0.0};
        }
    }
    transform_2d(rho_hat, nx, ny, false);

    std::vector<Complex> phi_hat(nx * ny, Complex{0.0, 0.0});
    std::vector<Complex> ex_hat(nx * ny, Complex{0.0, 0.0});
    std::vector<Complex> ey_hat(nx * ny, Complex{0.0, 0.0});

    for (std::size_t ky_index = 0; ky_index < ny; ++ky_index) {
        const double ky_mode = (ky_index <= ny / 2) ? static_cast<double>(ky_index) : -static_cast<double>(ny - ky_index);
        const double ky = twopi * ky_mode / Ly;
        for (std::size_t kx_index = 0; kx_index < nx; ++kx_index) {
            if (kx_index == 0 && ky_index == 0) continue;
            const double kx_mode = (kx_index <= nx / 2) ? static_cast<double>(kx_index) : -static_cast<double>(nx - kx_index);
            const double kx = twopi * kx_mode / Lx;
            const double k2 = kx * kx + ky * ky;
            const auto idx = mesh.index(kx_index, ky_index);
            const Complex phik =
                rho_hat[idx] / (permittivity_ * k2);
            phi_hat[idx] = phik;
            ex_hat[idx] = -Complex{0.0, 1.0} * kx * phik;
            ey_hat[idx] = -Complex{0.0, 1.0} * ky * phik;
        }
    }

    transform_2d(phi_hat, nx, ny, true);
    transform_2d(ex_hat, nx, ny, true);
    transform_2d(ey_hat, nx, ny, true);

    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            const auto idx = mesh.index(i, j);
            phi[idx] = phi_hat[idx].real();
            Ex[idx] = ex_hat[idx].real();
            Ey[idx] = ey_hat[idx].real();
        }
    }
}

void FieldSolver::solve_mixed_spectral_tridiagonal(Mesh2D& mesh) const {
    const bool periodic_y = mesh.boundary_y() == Boundary::Periodic;
    if (periodic_y ==
        (mesh.boundary_x() == Boundary::Periodic)) {
        throw std::invalid_argument(
            "mixed 2D spectral-tridiagonal solver requires exactly one periodic axis");
    }

    const std::size_t direct_size =
        periodic_y ? mesh.nx() : mesh.ny();
    const std::size_t periodic_size =
        periodic_y ? mesh.ny() : mesh.nx();
    const double direct_spacing =
        periodic_y ? mesh.dx() : mesh.dy();
    const double periodic_spacing =
        periodic_y ? mesh.dy() : mesh.dx();
    const double inverse_direct_spacing_squared =
        1.0 / (direct_spacing * direct_spacing);
    const double inverse_periodic_spacing_squared =
        1.0 / (periodic_spacing * periodic_spacing);
    const auto& boundary = mesh.boundary_config();
    const double lower_potential =
        periodic_y ? boundary.left.potential : boundary.bottom.potential;
    const double upper_potential =
        periodic_y ? boundary.right.potential : boundary.top.potential;
    const std::size_t interior_size = direct_size - 2;

    const auto mesh_index =
        [&mesh, periodic_y](std::size_t direct,
                            std::size_t periodic) {
            return periodic_y
                ? mesh.index(direct, periodic)
                : mesh.index(periodic, direct);
        };
    const auto spectral_index =
        [direct_size](std::size_t mode, std::size_t direct) {
            return mode * direct_size + direct;
        };

    std::vector<Complex> rho_hat(
        direct_size * periodic_size, Complex{0.0, 0.0});
    std::vector<Complex> line(periodic_size);
    for (std::size_t direct = 0; direct < direct_size; ++direct) {
        for (std::size_t periodic = 0;
             periodic < periodic_size;
             ++periodic) {
            line[periodic] =
                mesh.rho()[mesh_index(direct, periodic)];
        }
        transform_1d(line, false);
        for (std::size_t mode = 0; mode < periodic_size; ++mode) {
            rho_hat[spectral_index(mode, direct)] = line[mode];
        }
    }

    std::vector<Complex> phi_hat(
        direct_size * periodic_size, Complex{0.0, 0.0});
    std::vector<Complex> modified_upper(
        interior_size, Complex{0.0, 0.0});
    std::vector<Complex> modified_rhs(
        interior_size, Complex{0.0, 0.0});
    const double off_diagonal =
        -inverse_direct_spacing_squared;
    const double boundary_transform_scale =
        static_cast<double>(periodic_size);

    for (std::size_t mode = 0; mode < periodic_size; ++mode) {
        const double angle =
            std::numbers::pi * static_cast<double>(mode) /
            static_cast<double>(periodic_size);
        const double periodic_eigenvalue =
            4.0 * std::sin(angle) * std::sin(angle) *
            inverse_periodic_spacing_squared;
        const double diagonal =
            2.0 * inverse_direct_spacing_squared +
            periodic_eigenvalue;
        const Complex lower_boundary =
            mode == 0
                ? Complex{lower_potential *
                              boundary_transform_scale,
                          0.0}
                : Complex{0.0, 0.0};
        const Complex upper_boundary =
            mode == 0
                ? Complex{upper_potential *
                              boundary_transform_scale,
                          0.0}
                : Complex{0.0, 0.0};

        phi_hat[spectral_index(mode, 0)] = lower_boundary;
        phi_hat[spectral_index(mode, direct_size - 1)] =
            upper_boundary;

        for (std::size_t interior = 0;
             interior < interior_size;
             ++interior) {
            const std::size_t direct = interior + 1;
            Complex rhs =
                rho_hat[spectral_index(mode, direct)] /
                permittivity_;
            if (interior == 0) {
                rhs += inverse_direct_spacing_squared *
                       lower_boundary;
            }
            if (interior + 1 == interior_size) {
                rhs += inverse_direct_spacing_squared *
                       upper_boundary;
            }

            if (interior == 0) {
                modified_upper[interior] =
                    off_diagonal / diagonal;
                modified_rhs[interior] = rhs / diagonal;
            } else {
                const Complex denominator =
                    diagonal -
                    off_diagonal *
                        modified_upper[interior - 1];
                modified_upper[interior] =
                    interior + 1 == interior_size
                        ? Complex{0.0, 0.0}
                        : off_diagonal / denominator;
                modified_rhs[interior] =
                    (rhs -
                     off_diagonal *
                         modified_rhs[interior - 1]) /
                    denominator;
            }
        }

        phi_hat[spectral_index(mode, direct_size - 2)] =
            modified_rhs[interior_size - 1];
        for (std::size_t interior = interior_size - 1;
             interior-- > 0;) {
            phi_hat[spectral_index(mode, interior + 1)] =
                modified_rhs[interior] -
                modified_upper[interior] *
                    phi_hat[spectral_index(
                        mode, interior + 2)];
        }
    }

    for (std::size_t direct = 0; direct < direct_size; ++direct) {
        for (std::size_t mode = 0; mode < periodic_size; ++mode) {
            line[mode] = phi_hat[spectral_index(mode, direct)];
        }
        transform_1d(line, true);
        for (std::size_t periodic = 0;
             periodic < periodic_size;
             ++periodic) {
            mesh.phi()[mesh_index(direct, periodic)] =
                line[periodic].real();
        }
    }

    compute_electric_field(mesh);
}

void FieldSolver::solve_periodic_spectral(Mesh3D& mesh) const {
    const std::size_t nx = mesh.nx();
    const std::size_t ny = mesh.ny();
    const std::size_t nz = mesh.nz();
    const double Lx = mesh.length_x();
    const double Ly = mesh.length_y();
    const double Lz = mesh.length_z();
    const double twopi = 2.0 * std::numbers::pi;
    auto& rho = mesh.rho();
    auto& phi = mesh.phi();
    auto& Ex = mesh.electric_x();
    auto& Ey = mesh.electric_y();
    auto& Ez = mesh.electric_z();

    double mean = 0.0;
    for (double r : rho) mean += r;
    mean /= static_cast<double>(rho.size());

    std::vector<Complex> rho_hat(nx * ny * nz, Complex{0.0, 0.0});
    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j < ny; ++j) {
            for (std::size_t i = 0; i < nx; ++i) {
                const auto idx = mesh.index(i, j, k);
                rho_hat[idx] = Complex{rho[idx] - mean, 0.0};
            }
        }
    }
    transform_3d(rho_hat, nx, ny, nz, false);

    std::vector<Complex> phi_hat(nx * ny * nz, Complex{0.0, 0.0});
    std::vector<Complex> ex_hat(nx * ny * nz, Complex{0.0, 0.0});
    std::vector<Complex> ey_hat(nx * ny * nz, Complex{0.0, 0.0});
    std::vector<Complex> ez_hat(nx * ny * nz, Complex{0.0, 0.0});

    for (std::size_t kz_index = 0; kz_index < nz; ++kz_index) {
        const double kz_mode = (kz_index <= nz / 2) ? static_cast<double>(kz_index) : -static_cast<double>(nz - kz_index);
        const double kz = twopi * kz_mode / Lz;
        for (std::size_t ky_index = 0; ky_index < ny; ++ky_index) {
            const double ky_mode = (ky_index <= ny / 2) ? static_cast<double>(ky_index) : -static_cast<double>(ny - ky_index);
            const double ky = twopi * ky_mode / Ly;
            for (std::size_t kx_index = 0; kx_index < nx; ++kx_index) {
                if (kx_index == 0 && ky_index == 0 && kz_index == 0) continue;
                const double kx_mode = (kx_index <= nx / 2) ? static_cast<double>(kx_index) : -static_cast<double>(nx - kx_index);
                const double kx = twopi * kx_mode / Lx;
                const double k2 = kx * kx + ky * ky + kz * kz;
                const auto idx = mesh.index(kx_index, ky_index, kz_index);
                const Complex phik =
                    rho_hat[idx] / (permittivity_ * k2);
                phi_hat[idx] = phik;
                ex_hat[idx] = -Complex{0.0, 1.0} * kx * phik;
                ey_hat[idx] = -Complex{0.0, 1.0} * ky * phik;
                ez_hat[idx] = -Complex{0.0, 1.0} * kz * phik;
            }
        }
    }

    transform_3d(phi_hat, nx, ny, nz, true);
    transform_3d(ex_hat, nx, ny, nz, true);
    transform_3d(ey_hat, nx, ny, nz, true);
    transform_3d(ez_hat, nx, ny, nz, true);

    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j < ny; ++j) {
            for (std::size_t i = 0; i < nx; ++i) {
                const auto idx = mesh.index(i, j, k);
                phi[idx] = phi_hat[idx].real();
                Ex[idx] = ex_hat[idx].real();
                Ey[idx] = ey_hat[idx].real();
                Ez[idx] = ez_hat[idx].real();
            }
        }
    }
}

void FieldSolver::solve_dirichlet_tridiagonal(Grid& grid, double phi_left, double phi_right) const {
    const std::size_t n = grid.nx();
    auto& phi = grid.phi();
    auto& E = grid.electric();
    const auto& rho = grid.rho();
    const double dx2 = grid.dx() * grid.dx();
    phi[0] = phi_left;
    phi[n-1] = phi_right;
    const std::size_t m = n - 2;
    std::vector<double> c(m, -1.0), d(m, 0.0);
    for (std::size_t i = 0; i < m; ++i) {
        d[i] = rho[i + 1] * dx2 / permittivity_;
    }
    d[0] += phi_left; d[m-1] += phi_right;
    std::vector<double> cp(m, 0.0), dp(m, 0.0);
    cp[0] = c[0] / 2.0; dp[0] = d[0] / 2.0;
    for (std::size_t i = 1; i < m; ++i) {
        double denom = 2.0 + cp[i-1];
        cp[i] = (i == m-1) ? 0.0 : c[i] / denom;
        dp[i] = (d[i] + dp[i-1]) / denom;
    }
    phi[n-2] = dp[m-1];
    for (std::size_t ii = m - 1; ii-- > 0;) phi[ii+1] = dp[ii] - cp[ii] * phi[ii+2];
    for (std::size_t i = 1; i + 1 < n; ++i) E[i] = -(phi[i+1] - phi[i-1]) / (2.0 * grid.dx());
    // The Dirichlet endpoints represent half control volumes.  The potential
    // difference gives the field at the adjacent half-cell face; integrate
    // Gauss's law over the boundary half-cell to recover the field at the
    // electrode node itself.  This is also the endpoint convention used by
    // the pinned eduPIC reference.
    E[0] = -(phi[1] - phi[0]) / grid.dx() -
        rho[0] * grid.dx() / (2.0 * permittivity_);
    E[n-1] = -(phi[n-1] - phi[n-2]) / grid.dx() +
        rho[n-1] * grid.dx() / (2.0 * permittivity_);
}

void FieldSolver::solve_dirichlet_iterative(Mesh2D& mesh) const {
    const std::size_t nx = mesh.nx();
    const std::size_t ny = mesh.ny();
    auto& phi = mesh.phi();
    const auto& rho = mesh.rho();

    apply_dirichlet_boundary_potentials(mesh);

    const double inv_dx2 = 1.0 / (mesh.dx() * mesh.dx());
    const double inv_dy2 = 1.0 / (mesh.dy() * mesh.dy());
    const double diagonal = 2.0 * (inv_dx2 + inv_dy2);
    constexpr double omega = 1.7;
    constexpr double tolerance = 1.0e-11;
    const std::size_t max_iterations = std::max<std::size_t>(20000, 100 * nx * ny);

    bool converged = false;
    for (std::size_t iter = 0; iter < max_iterations; ++iter) {
        double max_delta = 0.0;
        for (std::size_t j = 0; j < ny; ++j) {
            if (mesh.boundary_y() == Boundary::Dirichlet &&
                (j == 0 || j + 1 == ny)) {
                continue;
            }
            const std::size_t jm =
                mesh.boundary_y() == Boundary::Periodic
                    ? (j + ny - 1) % ny
                    : j - 1;
            const std::size_t jp =
                mesh.boundary_y() == Boundary::Periodic
                    ? (j + 1) % ny
                    : j + 1;
            for (std::size_t i = 0; i < nx; ++i) {
                if (mesh.boundary_x() == Boundary::Dirichlet &&
                    (i == 0 || i + 1 == nx)) {
                    continue;
                }
                const std::size_t im =
                    mesh.boundary_x() == Boundary::Periodic
                        ? (i + nx - 1) % nx
                        : i - 1;
                const std::size_t ip =
                    mesh.boundary_x() == Boundary::Periodic
                        ? (i + 1) % nx
                        : i + 1;
                const auto idx = mesh.index(i, j);
                const double rhs = rho[idx] / permittivity_;
                const double updated =
                    ((phi[mesh.index(im, j)] +
                      phi[mesh.index(ip, j)]) * inv_dx2
                     + (phi[mesh.index(i, jm)] +
                        phi[mesh.index(i, jp)]) * inv_dy2
                                      + rhs) / diagonal;
                const double delta = omega * (updated - phi[idx]);
                phi[idx] += delta;
                max_delta = std::max(max_delta, std::abs(delta));
            }
        }
        if (max_delta < tolerance) {
            converged = true;
            break;
        }
    }
    if (!converged) {
        throw std::runtime_error("2D Dirichlet Poisson solve did not converge");
    }

    compute_electric_field(mesh);
}

void FieldSolver::solve_dirichlet_iterative(Mesh3D& mesh) const {
    const std::size_t nx = mesh.nx();
    const std::size_t ny = mesh.ny();
    const std::size_t nz = mesh.nz();
    auto& phi = mesh.phi();
    auto& Ex = mesh.electric_x();
    auto& Ey = mesh.electric_y();
    auto& Ez = mesh.electric_z();
    const auto& rho = mesh.rho();

    apply_grounded_dirichlet_boundary(mesh);

    const double inv_dx2 = 1.0 / (mesh.dx() * mesh.dx());
    const double inv_dy2 = 1.0 / (mesh.dy() * mesh.dy());
    const double inv_dz2 = 1.0 / (mesh.dz() * mesh.dz());
    const double diagonal = 2.0 * (inv_dx2 + inv_dy2 + inv_dz2);
    constexpr double omega = 1.65;
    constexpr double tolerance = 1.0e-11;
    const std::size_t max_iterations = std::max<std::size_t>(50000, 150 * nx * ny * nz);

    bool converged = false;
    for (std::size_t iter = 0; iter < max_iterations; ++iter) {
        double max_delta = 0.0;
        for (std::size_t k = 1; k + 1 < nz; ++k) {
            for (std::size_t j = 1; j + 1 < ny; ++j) {
                for (std::size_t i = 1; i + 1 < nx; ++i) {
                    const auto idx = mesh.index(i, j, k);
                    const double rhs = rho[idx] / permittivity_;
                    const double updated = ((phi[mesh.index(i - 1, j, k)] + phi[mesh.index(i + 1, j, k)]) * inv_dx2
                                          + (phi[mesh.index(i, j - 1, k)] + phi[mesh.index(i, j + 1, k)]) * inv_dy2
                                          + (phi[mesh.index(i, j, k - 1)] + phi[mesh.index(i, j, k + 1)]) * inv_dz2
                                          + rhs) / diagonal;
                    const double delta = omega * (updated - phi[idx]);
                    phi[idx] += delta;
                    max_delta = std::max(max_delta, std::abs(delta));
                }
            }
        }
        if (max_delta < tolerance) {
            converged = true;
            break;
        }
    }
    if (!converged) {
        throw std::runtime_error("3D Dirichlet Poisson solve did not converge");
    }

    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j < ny; ++j) {
            for (std::size_t i = 0; i < nx; ++i) {
                const auto idx = mesh.index(i, j, k);
                if (i == 0) {
                    Ex[idx] = -(phi[mesh.index(1, j, k)] - phi[idx]) / mesh.dx();
                } else if (i + 1 == nx) {
                    Ex[idx] = -(phi[idx] - phi[mesh.index(i - 1, j, k)]) / mesh.dx();
                } else {
                    Ex[idx] = -(phi[mesh.index(i + 1, j, k)] - phi[mesh.index(i - 1, j, k)]) / (2.0 * mesh.dx());
                }

                if (j == 0) {
                    Ey[idx] = -(phi[mesh.index(i, 1, k)] - phi[idx]) / mesh.dy();
                } else if (j + 1 == ny) {
                    Ey[idx] = -(phi[idx] - phi[mesh.index(i, j - 1, k)]) / mesh.dy();
                } else {
                    Ey[idx] = -(phi[mesh.index(i, j + 1, k)] - phi[mesh.index(i, j - 1, k)]) / (2.0 * mesh.dy());
                }

                if (k == 0) {
                    Ez[idx] = -(phi[mesh.index(i, j, 1)] - phi[idx]) / mesh.dz();
                } else if (k + 1 == nz) {
                    Ez[idx] = -(phi[idx] - phi[mesh.index(i, j, k - 1)]) / mesh.dz();
                } else {
                    Ez[idx] = -(phi[mesh.index(i, j, k + 1)] - phi[mesh.index(i, j, k - 1)]) / (2.0 * mesh.dz());
                }
            }
        }
    }
}

double interpolate_electric(const Grid& grid, double x) {
    const auto& E = grid.electric();
    const double dx = grid.dx();
    if (grid.boundary() == Boundary::Periodic) {
        double xp = wrap_periodic(x, grid.length());
        double g = xp / dx;
        auto i = static_cast<std::size_t>(std::floor(g));
        double f = g - static_cast<double>(i);
        return E[i % grid.nx()] * (1.0 - f) + E[(i + 1) % grid.nx()] * f;
    }
    double xp = std::clamp(x, 0.0, grid.length());
    double g = xp / dx;
    auto i = static_cast<std::size_t>(std::min<double>(std::floor(g), grid.nx() - 2));
    double f = g - static_cast<double>(i);
    return E[i] * (1.0 - f) + E[i+1] * f;
}

Vec2 interpolate_electric(const Mesh2D& mesh, Vec2 position) {
    const double x = mesh.boundary_x() == Boundary::Periodic
                         ? wrap_periodic(position.x, mesh.length_x())
                         : std::clamp(position.x, 0.0, mesh.length_x());
    const double y = mesh.boundary_y() == Boundary::Periodic
                         ? wrap_periodic(position.y, mesh.length_y())
                         : std::clamp(position.y, 0.0, mesh.length_y());
    const double gx = x / mesh.dx();
    const double gy = y / mesh.dy();

    std::size_t i = static_cast<std::size_t>(std::floor(gx));
    std::size_t j = static_cast<std::size_t>(std::floor(gy));
    double fx = gx - static_cast<double>(i);
    double fy = gy - static_cast<double>(j);

    std::size_t i0, i1, j0, j1;
    if (mesh.boundary_x() == Boundary::Periodic) {
        i0 = i % mesh.nx();
        i1 = (i + 1) % mesh.nx();
    } else {
        i = std::min(i, mesh.nx() - 2);
        fx = std::clamp(
            gx - static_cast<double>(i), 0.0, 1.0);
        i0 = i;
        i1 = i + 1;
    }
    if (mesh.boundary_y() == Boundary::Periodic) {
        j0 = j % mesh.ny();
        j1 = (j + 1) % mesh.ny();
    } else {
        j = std::min(j, mesh.ny() - 2);
        fy = std::clamp(gy - static_cast<double>(j), 0.0, 1.0);
        j0 = j;
        j1 = j + 1;
    }

    const double w00 = (1.0 - fx) * (1.0 - fy);
    const double w10 = fx * (1.0 - fy);
    const double w01 = (1.0 - fx) * fy;
    const double w11 = fx * fy;

    const auto& Ex = mesh.electric_x();
    const auto& Ey = mesh.electric_y();
    const auto idx00 = mesh.index(i0, j0);
    const auto idx10 = mesh.index(i1, j0);
    const auto idx01 = mesh.index(i0, j1);
    const auto idx11 = mesh.index(i1, j1);

    return Vec2{Ex[idx00] * w00 + Ex[idx10] * w10 + Ex[idx01] * w01 + Ex[idx11] * w11,
                Ey[idx00] * w00 + Ey[idx10] * w10 + Ey[idx01] * w01 + Ey[idx11] * w11};
}

Vec3 interpolate_electric(const Mesh3D& mesh, Vec3 position) {
    const double x = mesh.boundary() == Boundary::Periodic
                         ? wrap_periodic(position.x, mesh.length_x())
                         : std::clamp(position.x, 0.0, mesh.length_x());
    const double y = mesh.boundary() == Boundary::Periodic
                         ? wrap_periodic(position.y, mesh.length_y())
                         : std::clamp(position.y, 0.0, mesh.length_y());
    const double z = mesh.boundary() == Boundary::Periodic
                         ? wrap_periodic(position.z, mesh.length_z())
                         : std::clamp(position.z, 0.0, mesh.length_z());
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
    const auto& Ex = mesh.electric_x();
    const auto& Ey = mesh.electric_y();
    const auto& Ez = mesh.electric_z();

    Vec3 electric{};
    for (std::size_t dz_i = 0; dz_i < 2; ++dz_i) {
        for (std::size_t dy_i = 0; dy_i < 2; ++dy_i) {
            for (std::size_t dx_i = 0; dx_i < 2; ++dx_i) {
                const double shape = wx[dx_i] * wy[dy_i] * wz[dz_i];
                const auto idx = mesh.index(ii[dx_i], jj[dy_i], kk[dz_i]);
                electric.x += Ex[idx] * shape;
                electric.y += Ey[idx] * shape;
                electric.z += Ez[idx] * shape;
            }
        }
    }
    return electric;
}
}
