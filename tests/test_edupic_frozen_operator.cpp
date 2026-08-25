#include "pic/FieldSolver.hpp"
#include "pic/Pusher.hpp"
#include "pic/Units.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace {
constexpr const char* pinned_edupic_sha256 =
    "7c7679c0f0c98844940ea911bbb7581ec33f818e8d14427c9837ffdcf1ecea41";
constexpr double electron_charge_c = -1.6021766200000001e-19;
constexpr double electron_mass_kg = 9.1093835599999998e-31;
constexpr double dt_s = 1.8436578171091445e-11;
constexpr double threshold_ev = 15.8;

struct NativeField {
    std::vector<double> potential;
    std::vector<double> electric;
};

NativeField native_edupic_poisson(
    const std::vector<double>& rho, double length, double permittivity,
    double phi_left, double phi_right) {
    const std::size_t n = rho.size();
    const double dx = length / static_cast<double>(n - 1);
    std::vector<double> potential(n, 0.0);
    std::vector<double> electric(n, 0.0);
    std::vector<double> w(n, 0.0);
    std::vector<double> g(n, 0.0);
    std::vector<double> f(n, 0.0);
    potential.front() = phi_left;
    potential.back() = phi_right;
    const double alpha = -dx * dx / permittivity;
    for (std::size_t i = 1; i + 1 < n; ++i) f[i] = alpha * rho[i];
    f[1] -= potential.front();
    f[n - 2] -= potential.back();
    w[1] = -0.5;
    g[1] = f[1] / -2.0;
    for (std::size_t i = 2; i + 1 < n; ++i) {
        const double denominator = -2.0 - w[i - 1];
        w[i] = 1.0 / denominator;
        g[i] = (f[i] - g[i - 1]) / denominator;
    }
    potential[n - 2] = g[n - 2];
    for (std::size_t i = n - 2; i-- > 1;) {
        potential[i] = g[i] - w[i] * potential[i + 1];
    }
    for (std::size_t i = 1; i + 1 < n; ++i) {
        electric[i] = (potential[i - 1] - potential[i + 1]) /
            (2.0 * dx);
    }
    electric.front() = (potential.front() - potential[1]) / dx -
        rho.front() * dx / (2.0 * permittivity);
    electric.back() = (potential[n - 2] - potential.back()) / dx +
        rho.back() * dx / (2.0 * permittivity);
    return {std::move(potential), std::move(electric)};
}

double native_interpolate(
    const std::vector<double>& electric, double length, double x) {
    const double dx = length / static_cast<double>(electric.size() - 1);
    const double c0 = x / dx;
    const auto cell = static_cast<std::size_t>(std::min<double>(
        std::floor(c0), electric.size() - 2));
    const double c1 = static_cast<double>(cell) + 1.0 - c0;
    const double c2 = c0 - static_cast<double>(cell);
    return c1 * electric[cell] + c2 * electric[cell + 1];
}

double energy_ev(double vx, double vy, double vz) {
    return 0.5 * electron_mass_kg *
        (vx * vx + vy * vy + vz * vz) /
        pic::ELEMENTARY_CHARGE_SI;
}

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
} // namespace

int main() {
    try {
        constexpr std::size_t nodes = 401;
        constexpr std::size_t particles = 20000;
        constexpr double length_m = 0.025;
        constexpr double phi_left_v = 250.0;
        constexpr double phi_right_v = 0.0;
        const double permittivity = pic::VACUUM_PERMITTIVITY_SI;
        pic::Grid grid(nodes, length_m, pic::Boundary::Dirichlet);
        for (std::size_t i = 0; i < nodes; ++i) {
            const double fraction = static_cast<double>(i) /
                static_cast<double>(nodes - 1);
            grid.rho()[i] = 2.0e-6 * std::sin(
                2.0 * std::numbers::pi * fraction) +
                0.7e-6 * std::cos(6.0 * std::numbers::pi * fraction);
        }
        grid.rho().front() += 5.0e-6;
        grid.rho().back() -= 4.0e-6;
        const auto native = native_edupic_poisson(
            grid.rho(), length_m, permittivity, phi_left_v, phi_right_v);
        pic::FieldSolver solver(permittivity);
        solver.solve(grid, phi_left_v, phi_right_v);

        double maximum_potential_error = 0.0;
        double maximum_nodal_field_error = 0.0;
        for (std::size_t i = 0; i < nodes; ++i) {
            maximum_potential_error = std::max(
                maximum_potential_error,
                std::abs(grid.phi()[i] - native.potential[i]));
            maximum_nodal_field_error = std::max(
                maximum_nodal_field_error,
                std::abs(grid.electric()[i] - native.electric[i]));
        }

        const double threshold_speed = std::sqrt(
            2.0 * threshold_ev * pic::ELEMENTARY_CHARGE_SI /
            electron_mass_kg);
        const double charge_to_mass = electron_charge_c / electron_mass_kg;
        double maximum_interpolated_field_error = 0.0;
        double maximum_velocity_error = 0.0;
        double maximum_position_error = 0.0;
        std::size_t native_promotions = 0;
        std::size_t native_demotions = 0;
        std::size_t aurora_promotions = 0;
        std::size_t aurora_demotions = 0;
        std::size_t wall_adjacent_samples = 0;
        for (std::size_t i = 0; i < particles; ++i) {
            const double fraction =
                (static_cast<double>(i) + 0.375) /
                static_cast<double>(particles);
            const double x = length_m * fraction;
            const double transverse =
                (static_cast<int>(i % 9) - 4) * 12000.0;
            const double vy = transverse;
            const double vz = -0.5 * transverse;
            const double transverse_squared = vy * vy + vz * vz;
            const double threshold_vx = std::sqrt(std::max(
                0.0, threshold_speed * threshold_speed -
                    transverse_squared));
            const double vx = (i % 2 == 0)
                ? -(threshold_vx - 12000.0)
                : +(threshold_vx + 12000.0);
            const double native_field = native_interpolate(
                native.electric, length_m, x);
            const double aurora_field = pic::interpolate_electric(grid, x);
            maximum_interpolated_field_error = std::max(
                maximum_interpolated_field_error,
                std::abs(native_field - aurora_field));
            if (x < grid.dx() || x > length_m - grid.dx()) {
                ++wall_adjacent_samples;
            }

            const double native_vx_after =
                vx + charge_to_mass * native_field * dt_s;
            const double native_x_after = x + native_vx_after * dt_s;
            pic::Particle particle{};
            particle.x = x;
            particle.v_half = vx;
            particle.velocity_y = vy;
            particle.velocity_z = vz;
            pic::kick_leapfrog(
                particle, aurora_field, charge_to_mass, dt_s);
            pic::drift_leapfrog(particle, dt_s);
            maximum_velocity_error = std::max(
                maximum_velocity_error,
                std::abs(particle.v_half - native_vx_after));
            maximum_position_error = std::max(
                maximum_position_error,
                std::abs(particle.x - native_x_after));

            const bool before = energy_ev(vx, vy, vz) >= threshold_ev;
            const bool native_after =
                energy_ev(native_vx_after, vy, vz) >= threshold_ev;
            const bool aurora_after =
                energy_ev(particle.v_half, vy, vz) >= threshold_ev;
            native_promotions += !before && native_after;
            native_demotions += before && !native_after;
            aurora_promotions += !before && aurora_after;
            aurora_demotions += before && !aurora_after;
        }

        require(maximum_potential_error < 2.0e-12,
                "frozen eduPIC/AuroraPIC potential mismatch");
        require(maximum_nodal_field_error < 1.0e-7,
                "frozen eduPIC/AuroraPIC nodal field mismatch");
        require(maximum_interpolated_field_error < 1.0e-7,
                "frozen eduPIC/AuroraPIC interpolated field mismatch");
        require(maximum_velocity_error < 1.0e-8,
                "frozen eduPIC/AuroraPIC kick mismatch");
        require(maximum_position_error < 1.0e-16,
                "frozen eduPIC/AuroraPIC drift mismatch");
        require(native_promotions > 0 && native_demotions > 0,
                "frozen threshold population is not discriminating");
        require(native_promotions == aurora_promotions &&
                    native_demotions == aurora_demotions,
                "frozen eduPIC/AuroraPIC threshold traffic mismatch");
        require(wall_adjacent_samples > 0,
                "frozen state omitted wall-adjacent interpolation");

        std::cout << std::setprecision(17)
                  << "{\n"
                  << "  \"schema_version\": 1,\n"
                  << "  \"scope\": \"frozen_edupic_aurorapic_operator_equivalence\",\n"
                  << "  \"pinned_edupic_source_sha256\": \""
                  << pinned_edupic_sha256 << "\",\n"
                  << "  \"nodes\": " << nodes << ",\n"
                  << "  \"particles\": " << particles << ",\n"
                  << "  \"wall_adjacent_samples\": "
                  << wall_adjacent_samples << ",\n"
                  << "  \"maximum_potential_error_V\": "
                  << maximum_potential_error << ",\n"
                  << "  \"maximum_nodal_field_error_V_m\": "
                  << maximum_nodal_field_error << ",\n"
                  << "  \"maximum_interpolated_field_error_V_m\": "
                  << maximum_interpolated_field_error << ",\n"
                  << "  \"maximum_velocity_error_m_s\": "
                  << maximum_velocity_error << ",\n"
                  << "  \"maximum_position_error_m\": "
                  << maximum_position_error << ",\n"
                  << "  \"native_promotions\": " << native_promotions
                  << ",\n"
                  << "  \"aurorapic_promotions\": " << aurora_promotions
                  << ",\n"
                  << "  \"native_demotions\": " << native_demotions
                  << ",\n"
                  << "  \"aurorapic_demotions\": " << aurora_demotions
                  << ",\n"
                  << "  \"all_equivalence_gates_passed\": true,\n"
                  << "  \"claim_boundary\": \"Deterministic frozen-state formula equivalence; no collisions, self-consistent evolution, or experimental validation.\"\n"
                  << "}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "frozen operator comparison failed: "
                  << error.what() << '\n';
        return 1;
    }
}
