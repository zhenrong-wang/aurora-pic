#include "pic/Config.hpp"
#include "pic/FieldSolver.hpp"
#include "pic/Grid.hpp"
#include "pic/Mesh2D.hpp"
#include "pic/Mesh3D.hpp"
#include "pic/Simulation.hpp"
#include "pic/Simulation2D.hpp"
#include "pic/Simulation3D.hpp"
#include "pic/Pusher.hpp"
#include "pic/VTKWriter.hpp"
#include "pic/Species3D.hpp"
#include <algorithm>
#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void require_near(double actual, double expected, double tolerance, const std::string& message) {
    if (std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(message);
    }
}

std::string read_file_text(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read test file: " + path.string());
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}
std::size_t count_lines(const std::string& text) {
    return static_cast<std::size_t>(std::count(text.begin(), text.end(), '\n'));
}
template <typename Fn>
void require_throws(Fn&& fn, const std::string& message) {
    try {
        fn();
    } catch (const std::exception&) {
        return;
    }
    throw std::runtime_error(message);
}

void require_checkpoint_samples_close(const pic::DiagnosticSample& a,
                                      const pic::DiagnosticSample& b,
                                      const std::string& label) {
    require(a.step == b.step, label + ": step mismatch");
    require_near(a.time, b.time, 1e-15, label + ": time mismatch");
    require(a.live_particles == b.live_particles, label + ": live-particle mismatch");
    require_near(a.kinetic_energy, b.kinetic_energy, 1e-12, label + ": kinetic-energy mismatch");
    require_near(a.field_energy, b.field_energy, 1e-12, label + ": field-energy mismatch");
    require_near(a.total_energy, b.total_energy, 1e-12, label + ": total-energy mismatch");
    require_near(a.charge_l1, b.charge_l1, 1e-12, label + ": charge-l1 mismatch");
}

void require_checkpoint_samples_close(const pic::DiagnosticSample2D& a,
                                      const pic::DiagnosticSample2D& b,
                                      const std::string& label) {
    require(a.step == b.step, label + ": step mismatch");
    require_near(a.time, b.time, 1e-15, label + ": time mismatch");
    require(a.live_particles == b.live_particles, label + ": live-particle mismatch");
    require(a.live_particles_by_species == b.live_particles_by_species, label + ": per-species live-particle mismatch");
    require(a.boundary_losses.absorbed_left == b.boundary_losses.absorbed_left &&
            a.boundary_losses.absorbed_right == b.boundary_losses.absorbed_right &&
            a.boundary_losses.absorbed_bottom == b.boundary_losses.absorbed_bottom &&
            a.boundary_losses.absorbed_top == b.boundary_losses.absorbed_top,
            label + ": boundary-loss mismatch");
    require_near(a.kinetic_energy, b.kinetic_energy, 1e-12, label + ": kinetic-energy mismatch");
    require_near(a.field_energy, b.field_energy, 1e-12, label + ": field-energy mismatch");
    require_near(a.total_energy, b.total_energy, 1e-12, label + ": total-energy mismatch");
    require_near(a.charge_l1, b.charge_l1, 1e-12, label + ": charge-l1 mismatch");
}

void require_checkpoint_samples_close(const pic::DiagnosticSample3D& a,
                                      const pic::DiagnosticSample3D& b,
                                      const std::string& label) {
    require(a.step == b.step, label + ": step mismatch");
    require_near(a.time, b.time, 1e-15, label + ": time mismatch");
    require(a.live_particles == b.live_particles, label + ": live-particle mismatch");
    require(a.live_particles_by_species == b.live_particles_by_species, label + ": per-species live-particle mismatch");
    require(a.boundary_losses.absorbed_left == b.boundary_losses.absorbed_left &&
            a.boundary_losses.absorbed_right == b.boundary_losses.absorbed_right &&
            a.boundary_losses.absorbed_bottom == b.boundary_losses.absorbed_bottom &&
            a.boundary_losses.absorbed_top == b.boundary_losses.absorbed_top &&
            a.boundary_losses.absorbed_back == b.boundary_losses.absorbed_back &&
            a.boundary_losses.absorbed_front == b.boundary_losses.absorbed_front,
            label + ": boundary-loss mismatch");
    require_near(a.kinetic_energy, b.kinetic_energy, 1e-12, label + ": kinetic-energy mismatch");
    require_near(a.field_energy, b.field_energy, 1e-12, label + ": field-energy mismatch");
    require_near(a.total_energy, b.total_energy, 1e-12, label + ": total-energy mismatch");
    require_near(a.charge_l1, b.charge_l1, 1e-12, label + ": charge-l1 mismatch");
}

void require_species_close(const std::vector<pic::Species>& a,
                           const std::vector<pic::Species>& b,
                           const std::string& label) {
    require(a.size() == b.size(), label + ": species-count mismatch");
    for (std::size_t species_id = 0; species_id < a.size(); ++species_id) {
        require(a[species_id].name() == b[species_id].name(), label + ": species-name mismatch");
        const auto& pa = a[species_id].particles();
        const auto& pb = b[species_id].particles();
        require(pa.size() == pb.size(), label + ": particle-count mismatch");
        for (std::size_t i = 0; i < pa.size(); ++i) {
            require_near(pa[i].x, pb[i].x, 1e-12, label + ": particle x mismatch");
            require_near(pa[i].v, pb[i].v, 1e-12, label + ": particle v mismatch");
            require_near(pa[i].v_half, pb[i].v_half, 1e-12, label + ": particle v_half mismatch");
            require(pa[i].alive == pb[i].alive, label + ": particle alive mismatch");
        }
    }
}

void require_species_close(const std::vector<pic::Species2D>& a,
                           const std::vector<pic::Species2D>& b,
                           const std::string& label) {
    require(a.size() == b.size(), label + ": species-count mismatch");
    for (std::size_t species_id = 0; species_id < a.size(); ++species_id) {
        require(a[species_id].name() == b[species_id].name(), label + ": species-name mismatch");
        const auto& pa = a[species_id].particles();
        const auto& pb = b[species_id].particles();
        require(pa.size() == pb.size(), label + ": particle-count mismatch");
        for (std::size_t i = 0; i < pa.size(); ++i) {
            require_near(pa[i].position.x, pb[i].position.x, 1e-12, label + ": particle x mismatch");
            require_near(pa[i].position.y, pb[i].position.y, 1e-12, label + ": particle y mismatch");
            require_near(pa[i].velocity.x, pb[i].velocity.x, 1e-12, label + ": particle vx mismatch");
            require_near(pa[i].velocity.y, pb[i].velocity.y, 1e-12, label + ": particle vy mismatch");
            require_near(pa[i].velocity_half.x, pb[i].velocity_half.x, 1e-12, label + ": particle vx_half mismatch");
            require_near(pa[i].velocity_half.y, pb[i].velocity_half.y, 1e-12, label + ": particle vy_half mismatch");
            require(pa[i].alive == pb[i].alive, label + ": particle alive mismatch");
        }
    }
}

void require_species_close(const std::vector<pic::Species3D>& a,
                           const std::vector<pic::Species3D>& b,
                           const std::string& label) {
    require(a.size() == b.size(), label + ": species-count mismatch");
    for (std::size_t species_id = 0; species_id < a.size(); ++species_id) {
        require(a[species_id].name() == b[species_id].name(), label + ": species-name mismatch");
        const auto& pa = a[species_id].particles();
        const auto& pb = b[species_id].particles();
        require(pa.size() == pb.size(), label + ": particle-count mismatch");
        for (std::size_t i = 0; i < pa.size(); ++i) {
            require_near(pa[i].position.x, pb[i].position.x, 1e-12, label + ": particle x mismatch");
            require_near(pa[i].position.y, pb[i].position.y, 1e-12, label + ": particle y mismatch");
            require_near(pa[i].position.z, pb[i].position.z, 1e-12, label + ": particle z mismatch");
            require_near(pa[i].velocity.x, pb[i].velocity.x, 1e-12, label + ": particle vx mismatch");
            require_near(pa[i].velocity.y, pb[i].velocity.y, 1e-12, label + ": particle vy mismatch");
            require_near(pa[i].velocity.z, pb[i].velocity.z, 1e-12, label + ": particle vz mismatch");
            require_near(pa[i].velocity_half.x, pb[i].velocity_half.x, 1e-12, label + ": particle vx_half mismatch");
            require_near(pa[i].velocity_half.y, pb[i].velocity_half.y, 1e-12, label + ": particle vy_half mismatch");
            require_near(pa[i].velocity_half.z, pb[i].velocity_half.z, 1e-12, label + ": particle vz_half mismatch");
            require(pa[i].alive == pb[i].alive, label + ": particle alive mismatch");
        }
    }
}
double dot(pic::Vec3 a, pic::Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

double norm(pic::Vec3 v) {
    return std::sqrt(dot(v, v));
}

pic::Vec3 scale(pic::Vec3 v, double factor) {
    return pic::Vec3{factor * v.x, factor * v.y, factor * v.z};
}

pic::Vec3 cross(pic::Vec3 a, pic::Vec3 b) {
    return pic::Vec3{a.y * b.z - a.z * b.y,
                     a.z * b.x - a.x * b.z,
                     a.x * b.y - a.y * b.x};
}

pic::Vec3 rotate_about_axis(pic::Vec3 v, pic::Vec3 axis, double angle) {
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    const pic::Vec3 axis_cross_v = cross(axis, v);
    const double axis_dot_v = dot(axis, v);
    return pic::Vec3{v.x * c + axis_cross_v.x * s + axis.x * axis_dot_v * (1.0 - c),
                     v.y * c + axis_cross_v.y * s + axis.y * axis_dot_v * (1.0 - c),
                     v.z * c + axis_cross_v.z * s + axis.z * axis_dot_v * (1.0 - c)};
}

double boris_rotation_angle(double magnetic_magnitude, double charge_to_mass, double dt) {
    return -2.0 * std::atan(0.5 * charge_to_mass * magnetic_magnitude * dt);
}
}

int main() {
    try {
        {
            pic::Grid g(64, 1.0, pic::Boundary::Periodic);
            for (std::size_t i = 0; i < g.nx(); ++i) {
                g.rho()[i] = std::sin(2.0 * std::numbers::pi * g.node_x(i));
            }
            pic::FieldSolver solver;
            solver.solve(g);
            double max_err = 0.0;
            for (std::size_t i = 0; i < g.nx(); ++i) {
                double expected = -std::cos(2.0 * std::numbers::pi * g.node_x(i)) / (2.0 * std::numbers::pi);
                max_err = std::max(max_err, std::abs(g.electric()[i] - expected));
            }
            require(max_err < 1e-12, "periodic Poisson solve exceeded analytic error tolerance");
        }
        {
            constexpr double dt = 0.125;
            constexpr double acceleration = 2.5;
            pic::Particle particle{};
            particle.x = 0.2;
            particle.v = -0.4;
            pic::initialize_leapfrog_half_step(particle, acceleration, 1.0, dt);
            for (std::size_t n = 0; n < 8; ++n) {
                pic::kick_leapfrog(particle, acceleration, 1.0, dt);
                pic::drift_leapfrog(particle, dt);
                pic::synchronize_leapfrog(particle, acceleration, 1.0, dt);
            }
            const double t = 8.0 * dt;
            require_near(particle.x, 0.2 + (-0.4) * t + 0.5 * acceleration * t * t, 1e-14,
                         "1D leapfrog pusher does not match constant-acceleration position");
            require_near(particle.v, -0.4 + acceleration * t, 1e-14,
                         "1D leapfrog pusher does not match constant-acceleration velocity");
        }
        {
            constexpr double dt = 0.1;
            const pic::Vec2 acceleration{1.5, -0.75};
            pic::Particle2D particle{};
            particle.position = pic::Vec2{0.25, 0.75};
            particle.velocity = pic::Vec2{-0.2, 0.3};
            pic::initialize_leapfrog_half_step(particle, acceleration, 1.0, dt);
            for (std::size_t n = 0; n < 6; ++n) {
                pic::kick_leapfrog(particle, acceleration, 1.0, dt);
                pic::drift_leapfrog(particle, dt);
                pic::synchronize_leapfrog(particle, acceleration, 1.0, dt);
            }
            const double t = 6.0 * dt;
            require_near(particle.position.x, 0.25 + (-0.2) * t + 0.5 * acceleration.x * t * t, 1e-14,
                         "2D leapfrog pusher does not match constant-acceleration x position");
            require_near(particle.position.y, 0.75 + 0.3 * t + 0.5 * acceleration.y * t * t, 1e-14,
                         "2D leapfrog pusher does not match constant-acceleration y position");
            require_near(particle.velocity.x, -0.2 + acceleration.x * t, 1e-14,
                         "2D leapfrog pusher does not match constant-acceleration x velocity");
            require_near(particle.velocity.y, 0.3 + acceleration.y * t, 1e-14,
                         "2D leapfrog pusher does not match constant-acceleration y velocity");
        }
        {
            constexpr double dt = 0.04;
            constexpr std::size_t steps = 37;
            constexpr double charge_to_mass = 1.25;
            constexpr double magnetic_z = 1.7;
            const pic::Vec2 initial_velocity{0.8, -0.35};
            pic::Particle2D particle{};
            particle.velocity = initial_velocity;
            const pic::Vec2 electric{0.0, 0.0};
            pic::initialize_boris_half_step(particle, electric, magnetic_z, charge_to_mass, dt);
            for (std::size_t n = 0; n < steps; ++n) {
                pic::kick_boris(particle, electric, magnetic_z, charge_to_mass, dt);
            }
            pic::synchronize_boris(particle, electric, magnetic_z, charge_to_mass, dt);

            const double angle = static_cast<double>(steps) * boris_rotation_angle(std::abs(magnetic_z), charge_to_mass, dt);
            const double expected_x = initial_velocity.x * std::cos(angle) - initial_velocity.y * std::sin(angle);
            const double expected_y = initial_velocity.x * std::sin(angle) + initial_velocity.y * std::cos(angle);
            require_near(particle.velocity.x, expected_x, 1e-13,
                         "2D Boris pusher does not match cyclotron x-velocity rotation");
            require_near(particle.velocity.y, expected_y, 1e-13,
                         "2D Boris pusher does not match cyclotron y-velocity rotation");
            const double initial_speed = std::sqrt(initial_velocity.x * initial_velocity.x + initial_velocity.y * initial_velocity.y);
            const double final_speed = std::sqrt(particle.velocity.x * particle.velocity.x + particle.velocity.y * particle.velocity.y);
            require_near(final_speed, initial_speed, 1e-13, "2D Boris pusher did not conserve perpendicular speed");
        }
        {
            constexpr double dt = 0.03;
            constexpr std::size_t steps = 29;
            constexpr double charge_to_mass = -0.75;
            const pic::Vec3 magnetic{0.4, -0.8, 1.1};
            const pic::Vec3 initial_velocity{0.6, -0.25, 0.9};
            pic::Particle3D particle{};
            particle.velocity = initial_velocity;
            const pic::Vec3 electric{0.0, 0.0, 0.0};
            pic::initialize_boris_half_step(particle, electric, magnetic, charge_to_mass, dt);
            for (std::size_t n = 0; n < steps; ++n) {
                pic::kick_boris(particle, electric, magnetic, charge_to_mass, dt);
            }
            pic::synchronize_boris(particle, electric, magnetic, charge_to_mass, dt);

            const double magnetic_magnitude = norm(magnetic);
            const pic::Vec3 axis = scale(magnetic, 1.0 / magnetic_magnitude);
            const double angle = static_cast<double>(steps) * boris_rotation_angle(magnetic_magnitude, charge_to_mass, dt);
            const pic::Vec3 expected = rotate_about_axis(initial_velocity, axis, angle);
            require_near(particle.velocity.x, expected.x, 1e-13,
                         "3D Boris pusher does not match arbitrary-axis x-velocity rotation");
            require_near(particle.velocity.y, expected.y, 1e-13,
                         "3D Boris pusher does not match arbitrary-axis y-velocity rotation");
            require_near(particle.velocity.z, expected.z, 1e-13,
                         "3D Boris pusher does not match arbitrary-axis z-velocity rotation");
            require_near(norm(particle.velocity), norm(initial_velocity), 1e-13,
                         "3D Boris pusher did not conserve speed");
            require_near(dot(particle.velocity, axis), dot(initial_velocity, axis), 1e-13,
                         "3D Boris pusher did not conserve parallel velocity");
        }
        {
            constexpr double dt = 0.02;
            pic::Particle particle{};
            particle.x = 1.0;
            particle.v = 0.0;
            pic::initialize_leapfrog_half_step(particle, -particle.x, 1.0, dt);
            const double initial_energy = 0.5 * (particle.x * particle.x + particle.v * particle.v);
            double max_energy_error = 0.0;
            for (std::size_t n = 0; n < 5000; ++n) {
                pic::kick_leapfrog(particle, -particle.x, 1.0, dt);
                pic::drift_leapfrog(particle, dt);
                pic::synchronize_leapfrog(particle, -particle.x, 1.0, dt);
                const double energy = 0.5 * (particle.x * particle.x + particle.v * particle.v);
                max_energy_error = std::max(max_energy_error, std::abs(energy - initial_energy));
            }
            require(max_energy_error < 1e-3, "1D leapfrog pusher harmonic oscillator energy is not bounded");
        }
        {
            pic::Mesh2D mesh(16, 20, 1.0, 1.5, pic::Boundary::Periodic);
            const double kx = 2.0 * std::numbers::pi / mesh.length_x();
            const double ky = 4.0 * std::numbers::pi / mesh.length_y();
            const double k2 = kx * kx + ky * ky;
            for (std::size_t j = 0; j < mesh.ny(); ++j) {
                for (std::size_t i = 0; i < mesh.nx(); ++i) {
                    const double x = mesh.node_x(i);
                    const double y = mesh.node_y(j);
                    mesh.rho()[mesh.index(i, j)] = std::sin(kx * x) * std::cos(ky * y);
                }
            }
            pic::FieldSolver solver;
            solver.solve(mesh);
            double max_phi_err = 0.0;
            double max_ex_err = 0.0;
            double max_ey_err = 0.0;
            for (std::size_t j = 0; j < mesh.ny(); ++j) {
                for (std::size_t i = 0; i < mesh.nx(); ++i) {
                    const double x = mesh.node_x(i);
                    const double y = mesh.node_y(j);
                    const auto idx = mesh.index(i, j);
                    const double expected_phi = std::sin(kx * x) * std::cos(ky * y) / k2;
                    const double expected_ex = -kx * std::cos(kx * x) * std::cos(ky * y) / k2;
                    const double expected_ey = ky * std::sin(kx * x) * std::sin(ky * y) / k2;
                    max_phi_err = std::max(max_phi_err, std::abs(mesh.phi()[idx] - expected_phi));
                    max_ex_err = std::max(max_ex_err, std::abs(mesh.electric_x()[idx] - expected_ex));
                    max_ey_err = std::max(max_ey_err, std::abs(mesh.electric_y()[idx] - expected_ey));
                }
            }
            require(max_phi_err < 1e-12, "2D periodic Poisson potential exceeded analytic error tolerance");
            require(max_ex_err < 1e-12, "2D periodic Poisson Ex exceeded analytic error tolerance");
            require(max_ey_err < 1e-12, "2D periodic Poisson Ey exceeded analytic error tolerance");
        }
        {
            pic::Mesh3D mesh(8, 10, 12, 1.0, 1.25, 1.5, pic::Boundary::Periodic);
            const double kx = 2.0 * std::numbers::pi / mesh.length_x();
            const double ky = 4.0 * std::numbers::pi / mesh.length_y();
            const double kz = 6.0 * std::numbers::pi / mesh.length_z();
            const double k2 = kx * kx + ky * ky + kz * kz;
            for (std::size_t k = 0; k < mesh.nz(); ++k) {
                for (std::size_t j = 0; j < mesh.ny(); ++j) {
                    for (std::size_t i = 0; i < mesh.nx(); ++i) {
                        const double x = mesh.node_x(i);
                        const double y = mesh.node_y(j);
                        const double z = mesh.node_z(k);
                        mesh.rho()[mesh.index(i, j, k)] = std::sin(kx * x) * std::cos(ky * y) * std::sin(kz * z);
                    }
                }
            }
            pic::FieldSolver solver;
            solver.solve(mesh);
            double max_phi_err = 0.0;
            double max_ex_err = 0.0;
            double max_ey_err = 0.0;
            double max_ez_err = 0.0;
            for (std::size_t k = 0; k < mesh.nz(); ++k) {
                for (std::size_t j = 0; j < mesh.ny(); ++j) {
                    for (std::size_t i = 0; i < mesh.nx(); ++i) {
                        const double x = mesh.node_x(i);
                        const double y = mesh.node_y(j);
                        const double z = mesh.node_z(k);
                        const auto idx = mesh.index(i, j, k);
                        const double expected_phi = std::sin(kx * x) * std::cos(ky * y) * std::sin(kz * z) / k2;
                        const double expected_ex = -kx * std::cos(kx * x) * std::cos(ky * y) * std::sin(kz * z) / k2;
                        const double expected_ey = ky * std::sin(kx * x) * std::sin(ky * y) * std::sin(kz * z) / k2;
                        const double expected_ez = -kz * std::sin(kx * x) * std::cos(ky * y) * std::cos(kz * z) / k2;
                        max_phi_err = std::max(max_phi_err, std::abs(mesh.phi()[idx] - expected_phi));
                        max_ex_err = std::max(max_ex_err, std::abs(mesh.electric_x()[idx] - expected_ex));
                        max_ey_err = std::max(max_ey_err, std::abs(mesh.electric_y()[idx] - expected_ey));
                        max_ez_err = std::max(max_ez_err, std::abs(mesh.electric_z()[idx] - expected_ez));
                    }
                }
            }
            require(max_phi_err < 1e-12, "3D periodic Poisson potential exceeded analytic error tolerance");
            require(max_ex_err < 1e-12, "3D periodic Poisson Ex exceeded analytic error tolerance");
            require(max_ey_err < 1e-12, "3D periodic Poisson Ey exceeded analytic error tolerance");
            require(max_ez_err < 1e-12, "3D periodic Poisson Ez exceeded analytic error tolerance");
        }
        {
            pic::Mesh2D mesh(18, 14, 1.2, 0.8, pic::Boundary::Dirichlet);
            const double kx = std::numbers::pi / mesh.length_x();
            const double ky = std::numbers::pi / mesh.length_y();
            const double inv_dx2 = 1.0 / (mesh.dx() * mesh.dx());
            const double inv_dy2 = 1.0 / (mesh.dy() * mesh.dy());
            std::vector<double> expected(mesh.size(), 0.0);
            for (std::size_t j = 0; j < mesh.ny(); ++j) {
                for (std::size_t i = 0; i < mesh.nx(); ++i) {
                    const double x = mesh.node_x(i);
                    const double y = mesh.node_y(j);
                    expected[mesh.index(i, j)] = std::sin(kx * x) * std::sin(ky * y);
                }
            }
            for (std::size_t j = 0; j < mesh.ny(); ++j) {
                for (std::size_t i = 0; i < mesh.nx(); ++i) {
                    const auto idx = mesh.index(i, j);
                    if (i == 0 || j == 0 || i + 1 == mesh.nx() || j + 1 == mesh.ny()) {
                        mesh.phi()[idx] = expected[idx];
                    } else {
                        mesh.rho()[idx] = pic::EPS0 * ((2.0 * inv_dx2 + 2.0 * inv_dy2) * expected[idx]
                                      - inv_dx2 * (expected[mesh.index(i - 1, j)] + expected[mesh.index(i + 1, j)])
                                      - inv_dy2 * (expected[mesh.index(i, j - 1)] + expected[mesh.index(i, j + 1)]));
                    }
                }
            }
            pic::FieldSolver solver;
            solver.solve(mesh);
            double max_phi_err = 0.0;
            double max_ex_err = 0.0;
            double max_ey_err = 0.0;
            for (std::size_t j = 0; j < mesh.ny(); ++j) {
                for (std::size_t i = 0; i < mesh.nx(); ++i) {
                    const auto idx = mesh.index(i, j);
                    max_phi_err = std::max(max_phi_err, std::abs(mesh.phi()[idx] - expected[idx]));
                    const double expected_ex = (i == 0)
                        ? -(expected[mesh.index(1, j)] - expected[idx]) / mesh.dx()
                        : (i + 1 == mesh.nx())
                            ? -(expected[idx] - expected[mesh.index(i - 1, j)]) / mesh.dx()
                            : -(expected[mesh.index(i + 1, j)] - expected[mesh.index(i - 1, j)]) / (2.0 * mesh.dx());
                    const double expected_ey = (j == 0)
                        ? -(expected[mesh.index(i, 1)] - expected[idx]) / mesh.dy()
                        : (j + 1 == mesh.ny())
                            ? -(expected[idx] - expected[mesh.index(i, j - 1)]) / mesh.dy()
                            : -(expected[mesh.index(i, j + 1)] - expected[mesh.index(i, j - 1)]) / (2.0 * mesh.dy());
                    max_ex_err = std::max(max_ex_err, std::abs(mesh.electric_x()[idx] - expected_ex));
                    max_ey_err = std::max(max_ey_err, std::abs(mesh.electric_y()[idx] - expected_ey));
                }
            }
            require(max_phi_err < 1e-8, "2D Dirichlet Poisson potential exceeded discrete analytic tolerance");
            require(max_ex_err < 1e-8, "2D Dirichlet Poisson Ex exceeded discrete analytic tolerance");
            require(max_ey_err < 1e-8, "2D Dirichlet Poisson Ey exceeded discrete analytic tolerance");
        }
        {
            pic::Mesh3D mesh(8, 7, 6, 1.1, 0.9, 0.7, pic::Boundary::Dirichlet);
            const double kx = std::numbers::pi / mesh.length_x();
            const double ky = std::numbers::pi / mesh.length_y();
            const double kz = std::numbers::pi / mesh.length_z();
            const double inv_dx2 = 1.0 / (mesh.dx() * mesh.dx());
            const double inv_dy2 = 1.0 / (mesh.dy() * mesh.dy());
            const double inv_dz2 = 1.0 / (mesh.dz() * mesh.dz());
            std::vector<double> expected(mesh.size(), 0.0);
            for (std::size_t k = 0; k < mesh.nz(); ++k) {
                for (std::size_t j = 0; j < mesh.ny(); ++j) {
                    for (std::size_t i = 0; i < mesh.nx(); ++i) {
                        const double x = mesh.node_x(i);
                        const double y = mesh.node_y(j);
                        const double z = mesh.node_z(k);
                        expected[mesh.index(i, j, k)] = std::sin(kx * x) * std::sin(ky * y) * std::sin(kz * z);
                    }
                }
            }
            for (std::size_t k = 0; k < mesh.nz(); ++k) {
                for (std::size_t j = 0; j < mesh.ny(); ++j) {
                    for (std::size_t i = 0; i < mesh.nx(); ++i) {
                        const auto idx = mesh.index(i, j, k);
                        if (i == 0 || j == 0 || k == 0 || i + 1 == mesh.nx() || j + 1 == mesh.ny() || k + 1 == mesh.nz()) {
                            mesh.phi()[idx] = expected[idx];
                        } else {
                            mesh.rho()[idx] = pic::EPS0 * ((2.0 * inv_dx2 + 2.0 * inv_dy2 + 2.0 * inv_dz2) * expected[idx]
                                          - inv_dx2 * (expected[mesh.index(i - 1, j, k)] + expected[mesh.index(i + 1, j, k)])
                                          - inv_dy2 * (expected[mesh.index(i, j - 1, k)] + expected[mesh.index(i, j + 1, k)])
                                          - inv_dz2 * (expected[mesh.index(i, j, k - 1)] + expected[mesh.index(i, j, k + 1)]));
                        }
                    }
                }
            }
            pic::FieldSolver solver;
            solver.solve(mesh);
            double max_phi_err = 0.0;
            double max_ex_err = 0.0;
            double max_ey_err = 0.0;
            double max_ez_err = 0.0;
            for (std::size_t k = 0; k < mesh.nz(); ++k) {
                for (std::size_t j = 0; j < mesh.ny(); ++j) {
                    for (std::size_t i = 0; i < mesh.nx(); ++i) {
                        const auto idx = mesh.index(i, j, k);
                        max_phi_err = std::max(max_phi_err, std::abs(mesh.phi()[idx] - expected[idx]));
                        const double expected_ex = (i == 0)
                            ? -(expected[mesh.index(1, j, k)] - expected[idx]) / mesh.dx()
                            : (i + 1 == mesh.nx())
                                ? -(expected[idx] - expected[mesh.index(i - 1, j, k)]) / mesh.dx()
                                : -(expected[mesh.index(i + 1, j, k)] - expected[mesh.index(i - 1, j, k)]) / (2.0 * mesh.dx());
                        const double expected_ey = (j == 0)
                            ? -(expected[mesh.index(i, 1, k)] - expected[idx]) / mesh.dy()
                            : (j + 1 == mesh.ny())
                                ? -(expected[idx] - expected[mesh.index(i, j - 1, k)]) / mesh.dy()
                                : -(expected[mesh.index(i, j + 1, k)] - expected[mesh.index(i, j - 1, k)]) / (2.0 * mesh.dy());
                        const double expected_ez = (k == 0)
                            ? -(expected[mesh.index(i, j, 1)] - expected[idx]) / mesh.dz()
                            : (k + 1 == mesh.nz())
                                ? -(expected[idx] - expected[mesh.index(i, j, k - 1)]) / mesh.dz()
                                : -(expected[mesh.index(i, j, k + 1)] - expected[mesh.index(i, j, k - 1)]) / (2.0 * mesh.dz());
                        max_ex_err = std::max(max_ex_err, std::abs(mesh.electric_x()[idx] - expected_ex));
                        max_ey_err = std::max(max_ey_err, std::abs(mesh.electric_y()[idx] - expected_ey));
                        max_ez_err = std::max(max_ez_err, std::abs(mesh.electric_z()[idx] - expected_ez));
                    }
                }
            }
            require(max_phi_err < 1e-8, "3D Dirichlet Poisson potential exceeded discrete analytic tolerance");
            require(max_ex_err < 1e-8, "3D Dirichlet Poisson Ex exceeded discrete analytic tolerance");
            require(max_ey_err < 1e-8, "3D Dirichlet Poisson Ey exceeded discrete analytic tolerance");
            require(max_ez_err < 1e-8, "3D Dirichlet Poisson Ez exceeded discrete analytic tolerance");
        }
        {
            pic::BoundaryConfig2D electrodes;
            electrodes.left = {"cathode", -5.0};
            electrodes.right = {"anode", 5.0};
            electrodes.bottom = {"grounded_wall", 0.0};
            electrodes.top = {"biased_wall", 2.0};
            pic::Mesh2D mesh(6, 5, 1.0, 0.8, pic::Boundary::Dirichlet, electrodes);
            require(mesh.boundary_config().left.tag == "cathode", "2D mesh did not retain left boundary tag");
            require(mesh.boundary_config().right.tag == "anode", "2D mesh did not retain right boundary tag");
            require(std::abs(mesh.boundary_config().left.potential + 5.0) < 1e-15, "2D mesh did not retain left boundary potential");
            require(std::abs(mesh.boundary_config().right.potential - 5.0) < 1e-15, "2D mesh did not retain right boundary potential");

            pic::FieldSolver solver;
            solver.solve(mesh);
            for (std::size_t j = 1; j + 1 < mesh.ny(); ++j) {
                require(std::abs(mesh.phi()[mesh.index(0, j)] + 5.0) < 1e-12, "2D Dirichlet left electrode potential was not applied");
                require(std::abs(mesh.phi()[mesh.index(mesh.nx() - 1, j)] - 5.0) < 1e-12, "2D Dirichlet right electrode potential was not applied");
            }
            for (std::size_t i = 1; i + 1 < mesh.nx(); ++i) {
                require(std::abs(mesh.phi()[mesh.index(i, 0)]) < 1e-12, "2D Dirichlet bottom electrode potential was not applied");
                require(std::abs(mesh.phi()[mesh.index(i, mesh.ny() - 1)] - 2.0) < 1e-12, "2D Dirichlet top electrode potential was not applied");
            }
            require(std::abs(mesh.phi()[mesh.index(0, 0)] - (-2.5)) < 1e-12, "2D Dirichlet lower-left corner potential was not averaged");
            require(std::abs(mesh.phi()[mesh.index(mesh.nx() - 1, 0)] - 2.5) < 1e-12, "2D Dirichlet lower-right corner potential was not averaged");
            require(std::abs(mesh.phi()[mesh.index(0, mesh.ny() - 1)] - (-1.5)) < 1e-12, "2D Dirichlet upper-left corner potential was not averaged");
            require(std::abs(mesh.phi()[mesh.index(mesh.nx() - 1, mesh.ny() - 1)] - 3.5) < 1e-12, "2D Dirichlet upper-right corner potential was not averaged");
        }
        {
            pic::Mesh2D mesh(4, 4, 1.0, 1.0, pic::Boundary::Periodic);
            const std::vector<pic::Particle2D> particles{
                pic::Particle2D{pic::Vec2{0.125, 0.125}, pic::Vec2{}, true},
                pic::Particle2D{pic::Vec2{1.125, -0.125}, pic::Vec2{}, true},
                pic::Particle2D{pic::Vec2{0.5, 0.5}, pic::Vec2{}, false},
            };
            pic::deposit_charge_cic(mesh, particles, 2.0, 0.5);
            double total_charge = 0.0;
            for (double rho : mesh.rho()) total_charge += rho * mesh.dx() * mesh.dy();
            require(std::abs(total_charge - 2.0) < 1e-12, "2D CIC deposition did not conserve live-particle charge");
            require(std::abs(mesh.rho()[mesh.index(0, 0)] - 8.0) < 1e-12, "2D CIC wrapped corner weight is wrong");
            require(std::abs(mesh.rho()[mesh.index(1, 0)] - 8.0) < 1e-12, "2D CIC wrapped x weight is wrong");
            require(std::abs(mesh.rho()[mesh.index(0, 1)] - 4.0) < 1e-12, "2D CIC lower-y weight is wrong");
        }
        {
            pic::Mesh3D mesh(4, 4, 4, 1.0, 1.0, 1.0, pic::Boundary::Periodic);
            const std::vector<pic::Particle3D> particles{
                pic::Particle3D{pic::Vec3{0.125, 0.125, 0.125}, pic::Vec3{}, true},
                pic::Particle3D{pic::Vec3{1.125, -0.125, 0.875}, pic::Vec3{}, true},
                pic::Particle3D{pic::Vec3{0.5, 0.5, 0.5}, pic::Vec3{}, false},
            };
            pic::deposit_charge_cic(mesh, particles, 2.0, 0.5);
            double total_charge = 0.0;
            for (double rho : mesh.rho()) total_charge += rho * mesh.dx() * mesh.dy() * mesh.dz();
            require(std::abs(total_charge - 2.0) < 1e-12, "3D CIC deposition did not conserve live-particle charge");
            require(std::abs(mesh.rho()[mesh.index(0, 0, 0)] - 16.0) < 1e-12, "3D CIC wrapped corner weight is wrong");
            require(std::abs(mesh.rho()[mesh.index(1, 0, 0)] - 16.0) < 1e-12, "3D CIC wrapped x weight is wrong");
            require(std::abs(mesh.rho()[mesh.index(0, 1, 0)] - 8.0) < 1e-12, "3D CIC lower-y weight is wrong");
            require(std::abs(mesh.rho()[mesh.index(0, 0, 3)] - 8.0) < 1e-12, "3D CIC wrapped z weight is wrong");
        }
        {
            pic::Mesh3D mesh(5, 4, 3, 1.0, 0.8, 0.6, pic::Boundary::Periodic);
            for (std::size_t k = 0; k < mesh.nz(); ++k) {
                for (std::size_t j = 0; j < mesh.ny(); ++j) {
                    for (std::size_t i = 0; i < mesh.nx(); ++i) {
                        const double x = mesh.node_x(i);
                        const double y = mesh.node_y(j);
                        const double z = mesh.node_z(k);
                        const auto idx = mesh.index(i, j, k);
                        mesh.electric_x()[idx] = 1.0 + 2.0 * x - 3.0 * y + 0.5 * z;
                        mesh.electric_y()[idx] = -0.25 + 0.75 * x + 1.5 * y - 2.0 * z;
                        mesh.electric_z()[idx] = 3.0 - x + 0.25 * y + 4.0 * z;
                    }
                }
            }
            const pic::Vec3 position{0.37, 0.31, 0.26};
            const pic::Vec3 electric = pic::interpolate_electric(mesh, position);
            require(std::abs(electric.x - (1.0 + 2.0 * position.x - 3.0 * position.y + 0.5 * position.z)) < 1e-12,
                    "3D electric interpolation Ex is not trilinear-exact for affine fields");
            require(std::abs(electric.y - (-0.25 + 0.75 * position.x + 1.5 * position.y - 2.0 * position.z)) < 1e-12,
                    "3D electric interpolation Ey is not trilinear-exact for affine fields");
            require(std::abs(electric.z - (3.0 - position.x + 0.25 * position.y + 4.0 * position.z)) < 1e-12,
                    "3D electric interpolation Ez is not trilinear-exact for affine fields");

            const pic::Vec3 wrapped = pic::interpolate_electric(mesh, pic::Vec3{position.x + mesh.length_x(), position.y - mesh.length_y(), position.z});
            require(std::abs(wrapped.x - electric.x) < 1e-12, "3D periodic electric interpolation did not wrap x/y coordinates");
            require(std::abs(wrapped.y - electric.y) < 1e-12, "3D periodic electric interpolation did not wrap x/y coordinates");
            require(std::abs(wrapped.z - electric.z) < 1e-12, "3D periodic electric interpolation did not wrap x/y coordinates");
        }
        {
            pic::Grid periodic_grid(4, 1.0, pic::Boundary::Periodic);
            double periodic_volume = 0.0;
            for (std::size_t i = 0; i < periodic_grid.nx(); ++i) periodic_volume += periodic_grid.node_volume(i);
            require(std::abs(periodic_volume - periodic_grid.length()) < 1e-15, "periodic grid nodal volumes do not sum to domain length");
            require(std::abs(periodic_grid.node_volume(0) - periodic_grid.dx()) < 1e-15, "periodic grid node volume is wrong");

            pic::Grid dirichlet_grid(5, 1.0, pic::Boundary::Dirichlet);
            double dirichlet_volume = 0.0;
            for (std::size_t i = 0; i < dirichlet_grid.nx(); ++i) dirichlet_volume += dirichlet_grid.node_volume(i);
            require(std::abs(dirichlet_volume - dirichlet_grid.length()) < 1e-15, "Dirichlet grid nodal volumes do not sum to domain length");
            require(std::abs(dirichlet_grid.node_volume(0) - 0.5 * dirichlet_grid.dx()) < 1e-15, "Dirichlet grid boundary node volume is wrong");

            pic::Mesh2D periodic_mesh(4, 5, 1.0, 2.0, pic::Boundary::Periodic);
            double periodic_area = 0.0;
            for (std::size_t j = 0; j < periodic_mesh.ny(); ++j) {
                for (std::size_t i = 0; i < periodic_mesh.nx(); ++i) periodic_area += periodic_mesh.node_area(i, j);
            }
            require(std::abs(periodic_area - periodic_mesh.length_x() * periodic_mesh.length_y()) < 1e-15,
                    "periodic mesh nodal areas do not sum to domain area");

            pic::Mesh2D dirichlet_mesh(5, 4, 1.0, 2.0, pic::Boundary::Dirichlet);
            double dirichlet_area = 0.0;
            for (std::size_t j = 0; j < dirichlet_mesh.ny(); ++j) {
                for (std::size_t i = 0; i < dirichlet_mesh.nx(); ++i) dirichlet_area += dirichlet_mesh.node_area(i, j);
            }
            require(std::abs(dirichlet_area - dirichlet_mesh.length_x() * dirichlet_mesh.length_y()) < 1e-15,
                    "Dirichlet mesh nodal areas do not sum to domain area");
            require(std::abs(dirichlet_mesh.node_area(0, 0) - 0.25 * dirichlet_mesh.dx() * dirichlet_mesh.dy()) < 1e-15,
                    "Dirichlet mesh corner node area is wrong");

            pic::Mesh3D periodic_mesh3d(4, 5, 6, 1.0, 2.0, 3.0, pic::Boundary::Periodic);
            double periodic_volume3d = 0.0;
            for (std::size_t k = 0; k < periodic_mesh3d.nz(); ++k) {
                for (std::size_t j = 0; j < periodic_mesh3d.ny(); ++j) {
                    for (std::size_t i = 0; i < periodic_mesh3d.nx(); ++i) periodic_volume3d += periodic_mesh3d.node_volume(i, j, k);
                }
            }
            require(std::abs(periodic_volume3d - periodic_mesh3d.length_x() * periodic_mesh3d.length_y() * periodic_mesh3d.length_z()) < 1e-12,
                    "periodic 3D mesh nodal volumes do not sum to domain volume");

            pic::Mesh3D dirichlet_mesh3d(5, 4, 6, 1.0, 2.0, 3.0, pic::Boundary::Dirichlet);
            double dirichlet_volume3d = 0.0;
            for (std::size_t k = 0; k < dirichlet_mesh3d.nz(); ++k) {
                for (std::size_t j = 0; j < dirichlet_mesh3d.ny(); ++j) {
                    for (std::size_t i = 0; i < dirichlet_mesh3d.nx(); ++i) dirichlet_volume3d += dirichlet_mesh3d.node_volume(i, j, k);
                }
            }
            require(std::abs(dirichlet_volume3d - dirichlet_mesh3d.length_x() * dirichlet_mesh3d.length_y() * dirichlet_mesh3d.length_z()) < 1e-12,
                    "Dirichlet 3D mesh nodal volumes do not sum to domain volume");
            require(std::abs(dirichlet_mesh3d.node_volume(0, 0, 0) - 0.125 * dirichlet_mesh3d.dx() * dirichlet_mesh3d.dy() * dirichlet_mesh3d.dz()) < 1e-15,
                    "Dirichlet 3D mesh corner node volume is wrong");
        }
        {
            pic::Grid grid(5, 1.0, pic::Boundary::Dirichlet);
            pic::SpeciesConfig species_cfg;
            species_cfg.charge = 2.0;
            species_cfg.weight = 0.5;
            species_cfg.particles = 1;
            pic::Species species(species_cfg);
            species.particles() = {pic::Particle{0.0, 0.0, true}};
            species.deposit_charge(grid);
            double deposited_charge = 0.0;
            for (std::size_t i = 0; i < grid.nx(); ++i) deposited_charge += grid.rho()[i] * grid.node_volume(i);
            require(std::abs(deposited_charge - 1.0) < 1e-12, "1D Dirichlet CIC deposition did not conserve boundary-node charge");
            require(std::abs(grid.rho()[0] - 8.0) < 1e-12, "1D Dirichlet boundary-node density did not use half control volume");
        }
        {
            pic::Mesh2D mesh(5, 5, 1.0, 1.0, pic::Boundary::Dirichlet);
            const std::vector<pic::Particle2D> particles{pic::Particle2D{pic::Vec2{0.0, 0.0}, pic::Vec2{}, true}};
            pic::deposit_charge_cic(mesh, particles, 2.0, 0.5);
            double deposited_charge = 0.0;
            for (std::size_t j = 0; j < mesh.ny(); ++j) {
                for (std::size_t i = 0; i < mesh.nx(); ++i) deposited_charge += mesh.rho()[mesh.index(i, j)] * mesh.node_area(i, j);
            }
            require(std::abs(deposited_charge - 1.0) < 1e-12, "2D Dirichlet CIC deposition did not conserve corner-node charge");
            require(std::abs(mesh.rho()[mesh.index(0, 0)] - 64.0) < 1e-12, "2D Dirichlet corner density did not use quarter control volume");

            mesh.electric_x()[mesh.index(0, 0)] = 2.0;
            mesh.rho()[mesh.index(0, 0)] = -64.0;
            pic::Diagnostics2D diagnostics("test_output_quadrature", {});
            auto sample = diagnostics.sample(0, 0.0, mesh, {});
            require(std::abs(sample.charge_l1 - 1.0) < 1e-12, "2D diagnostics charge_l1 did not use nodal control volumes");
            require(std::abs(sample.field_energy - 0.03125) < 1e-12, "2D diagnostics field energy did not use nodal control volumes");
            std::filesystem::remove_all("test_output_quadrature");
        }
        {
            pic::Mesh3D mesh(5, 5, 5, 1.0, 1.0, 1.0, pic::Boundary::Dirichlet);
            const std::vector<pic::Particle3D> particles{pic::Particle3D{pic::Vec3{0.0, 0.0, 0.0}, pic::Vec3{}, true}};
            pic::deposit_charge_cic(mesh, particles, 2.0, 0.5);
            double deposited_charge = 0.0;
            for (std::size_t k = 0; k < mesh.nz(); ++k) {
                for (std::size_t j = 0; j < mesh.ny(); ++j) {
                    for (std::size_t i = 0; i < mesh.nx(); ++i) deposited_charge += mesh.rho()[mesh.index(i, j, k)] * mesh.node_volume(i, j, k);
                }
            }
            require(std::abs(deposited_charge - 1.0) < 1e-12, "3D Dirichlet CIC deposition did not conserve corner-node charge");
            require(std::abs(mesh.rho()[mesh.index(0, 0, 0)] - 512.0) < 1e-12, "3D Dirichlet corner density did not use eighth control volume");
        }
        {
            pic::Mesh3D mesh(5, 5, 5, 1.0, 1.0, 1.0, pic::Boundary::Dirichlet);
            pic::Species3DConfig cfg;
            cfg.mass = 2.0;
            cfg.weight = 0.25;
            cfg.particles = 3;
            cfg.drift_velocity_x = 1.0;
            cfg.drift_velocity_y = -2.0;
            cfg.drift_velocity_z = 0.5;
            cfg.thermal_velocity = 0.0;
            cfg.init_x_min = 0.1;
            cfg.init_x_max = 0.2;
            cfg.init_y_min = 0.3;
            cfg.init_y_max = 0.4;
            cfg.init_z_min = 0.5;
            cfg.init_z_max = 0.6;
            pic::Species3D species(cfg);
            std::mt19937_64 rng(1234);
            species.initialize(mesh, rng);
            require(species.live_count() == 3, "3D species did not initialize all particles as live");
            for (const auto& particle : species.particles()) {
                require(particle.position.x >= 0.1 && particle.position.x <= 0.2, "3D species x initialization interval was not honored");
                require(particle.position.y >= 0.3 && particle.position.y <= 0.4, "3D species y initialization interval was not honored");
                require(particle.position.z >= 0.5 && particle.position.z <= 0.6, "3D species z initialization interval was not honored");
                require_near(particle.velocity.x, 1.0, 1e-15, "3D species vx drift initialization is wrong");
                require_near(particle.velocity.y, -2.0, 1e-15, "3D species vy drift initialization is wrong");
                require_near(particle.velocity.z, 0.5, 1e-15, "3D species vz drift initialization is wrong");
                require_near(particle.velocity_half.x, particle.velocity.x, 1e-15, "3D species half-step vx initialization is wrong");
                require_near(particle.velocity_half.y, particle.velocity.y, 1e-15, "3D species half-step vy initialization is wrong");
                require_near(particle.velocity_half.z, particle.velocity.z, 1e-15, "3D species half-step vz initialization is wrong");
            }
            require_near(species.kinetic_energy(), 3.0 * 0.5 * cfg.mass * cfg.weight * (1.0 + 4.0 + 0.25), 1e-15,
                         "3D species kinetic energy is wrong");
            mesh.clear_charge();
            species.deposit_charge(mesh);
            double deposited_charge = 0.0;
            for (std::size_t k = 0; k < mesh.nz(); ++k) {
                for (std::size_t j = 0; j < mesh.ny(); ++j) {
                    for (std::size_t i = 0; i < mesh.nx(); ++i) deposited_charge += mesh.rho()[mesh.index(i, j, k)] * mesh.node_volume(i, j, k);
                }
            }
            require(std::abs(deposited_charge - cfg.charge * cfg.weight * static_cast<double>(cfg.particles)) < 1e-12,
                    "3D species deposition did not conserve charge");
        }
        {
            const auto output_dir = std::filesystem::path("test_output_vtk");
            std::filesystem::remove_all(output_dir);

            pic::Mesh2D mesh(3, 4, 1.0, 2.0, pic::Boundary::Dirichlet);
            for (std::size_t j = 0; j < mesh.ny(); ++j) {
                for (std::size_t i = 0; i < mesh.nx(); ++i) {
                    const auto idx = mesh.index(i, j);
                    mesh.rho()[idx] = static_cast<double>(idx);
                    mesh.phi()[idx] = 10.0 + static_cast<double>(idx);
                    mesh.electric_x()[idx] = 0.25 * static_cast<double>(i);
                    mesh.electric_y()[idx] = -0.5 * static_cast<double>(j);
                }
            }
            pic::write_legacy_vtk(mesh, output_dir / "manual_fields.vtk", "AuroraPIC test fields");
            const auto vtk = read_file_text(output_dir / "manual_fields.vtk");
            require(vtk.find("# vtk DataFile Version 3.0") != std::string::npos, "VTK header is missing");
            require(vtk.find("DATASET STRUCTURED_GRID") != std::string::npos, "VTK dataset type is missing");
            require(vtk.find("DIMENSIONS 3 4 1") != std::string::npos, "VTK dimensions are wrong");
            require(vtk.find("POINT_DATA 12") != std::string::npos, "VTK point-data count is wrong");
            require(vtk.find("SCALARS rho double 1") != std::string::npos, "VTK rho scalar is missing");
            require(vtk.find("SCALARS phi double 1") != std::string::npos, "VTK phi scalar is missing");
            require(vtk.find("VECTORS electric double") != std::string::npos, "VTK electric vector is missing");
        }
        {
            require_throws([] { pic::Grid(2, 1.0, pic::Boundary::Periodic); }, "grid nx validation did not throw");
            require_throws([] { pic::Grid(32, 0.0, pic::Boundary::Periodic); }, "grid length validation did not throw");
            require_throws([] { pic::Mesh2D(2, 4, 1.0, 1.0, pic::Boundary::Periodic); }, "2D mesh nx validation did not throw");
            require_throws([] { pic::Mesh2D(4, 4, 0.0, 1.0, pic::Boundary::Periodic); }, "2D mesh length validation did not throw");
            require_throws([] { pic::Simulation2DConfig cfg; cfg.output_interval = 0; pic::Simulation2D sim(cfg); }, "2D output_interval validation did not throw");
            require_throws([] { pic::Config cfg; cfg.output_interval = 0; pic::Simulation sim(cfg); }, "1D output_interval validation did not throw");
        }
        {
            const auto config_path = std::filesystem::path("test_density_config.ini");
            {
                std::ofstream out(config_path);
                out << "nx = 16\n"
                    << "length = 2.0\n"
                    << "dt = 0.01\n"
                    << "output_interval = 2\n"
                    << "[species]\n"
                    << "name = density_weighted\n"
                    << "charge = -1\n"
                    << "mass = 1\n"
                    << "density = 5\n"
                    << "particles = 10\n"
                    << "thermal_velocity = 0\n"
                    << "init_x_min = 0.5\n"
                    << "init_x_max = 1.5\n";
            }
            auto cfg = pic::load_config(config_path.string());
            require(cfg.species.size() == 1, "config did not load one species");
            require(std::abs(cfg.species[0].weight - 0.5) < 1e-15, "density-derived macro-particle weight is wrong");
            std::filesystem::remove(config_path);

            const auto config_2d_path = std::filesystem::path("test_2d_config.ini");
            {
                std::ofstream out(config_2d_path);
                out << "dimension = 2\n"
                    << "nx = 8\n"
                    << "ny = 6\n"
                    << "length_x = 2.0\n"
                    << "length_y = 1.5\n"
                    << "dt = 0.01\n"
                    << "steps = 2\n"
                    << "output_interval = 1\n"
                    << "output_dir = test_output_config_2d\n"
                    << "vtk_output = true\n"
                    << "particle_output = true\n"
                    << "particle_output_interval = 3\n"
                    << "particle_output_stride = 2\n"
                    << "particle_sample_count = 7\n"
                    << "particle_boundary = reflecting\n"
                    << "particle_boundary_left = absorbing\n"
                    << "particle_boundary_right = periodic\n"
                    << "particle_boundary_bottom = auto\n"
                    << "particle_boundary_top = reflect\n"
                    << "boundary = dirichlet\n"
                    << "phi_left = -2.5\n"
                    << "phi_right = 3.5\n"
                    << "phi_bottom = 0.25\n"
                    << "phi_top = 1.25\n"
                    << "boundary_left_tag = cathode\n"
                    << "boundary_right_tag = anode\n"
                    << "boundary_bottom_tag = lower_wall\n"
                    << "boundary_top_tag = upper_wall\n"
                    << "magnetic_field_z = 1.75\n"
                    << "[species.electrons]\n"
                    << "charge = -1\n"
                    << "mass = 1\n"
                    << "density = 4\n"
                    << "particles = 12\n"
                    << "thermal_velocity = 0\n"
                    << "init_x_min = 0.5\n"
                    << "init_x_max = 1.5\n"
                    << "init_y_min = 0.25\n"
                    << "init_y_max = 1.25\n";
            }
            require(pic::detect_config_dimension(config_2d_path.string()) == 2, "2D config dimension was not detected");
            auto cfg2 = pic::load_config_2d(config_2d_path.string());
            require(cfg2.nx == 8 && cfg2.ny == 6, "2D config did not load mesh dimensions");
            require(cfg2.vtk_output, "2D config did not load vtk_output");
            require(cfg2.particle_output, "2D config did not load particle_output");
            require(cfg2.particle_output_interval == 3, "2D config did not load particle_output_interval");
            require(cfg2.particle_output_stride == 2, "2D config did not load particle_output_stride");
            require(cfg2.particle_sample_count == 7, "2D config did not load particle_sample_count");
            require(cfg2.particle_boundary_config.left == pic::ParticleBoundary::Absorbing, "2D config did not load left particle boundary");
            require(cfg2.particle_boundary_config.right == pic::ParticleBoundary::Periodic, "2D config did not load right particle boundary");
            require(cfg2.particle_boundary_config.bottom == pic::ParticleBoundary::Auto, "2D config did not load bottom particle boundary");
            require(cfg2.particle_boundary_config.top == pic::ParticleBoundary::Reflecting, "2D config did not load top particle boundary");
            require(cfg2.boundary == pic::Boundary::Dirichlet, "2D config did not load Dirichlet boundary mode");
            require(std::abs(cfg2.boundary_config.left.potential + 2.5) < 1e-15, "2D config did not load phi_left");
            require(std::abs(cfg2.boundary_config.right.potential - 3.5) < 1e-15, "2D config did not load phi_right");
            require(std::abs(cfg2.boundary_config.bottom.potential - 0.25) < 1e-15, "2D config did not load phi_bottom");
            require(std::abs(cfg2.boundary_config.top.potential - 1.25) < 1e-15, "2D config did not load phi_top");
            require(cfg2.boundary_config.left.tag == "cathode", "2D config did not load left boundary tag");
            require(cfg2.boundary_config.right.tag == "anode", "2D config did not load right boundary tag");
            require(cfg2.boundary_config.bottom.tag == "lower_wall", "2D config did not load bottom boundary tag");
            require(cfg2.boundary_config.top.tag == "upper_wall", "2D config did not load top boundary tag");
            require(std::abs(cfg2.magnetic_field_z - 1.75) < 1e-15, "2D config did not load magnetic_field_z");
            require(cfg2.species.size() == 1, "2D config did not load one species");
            require(std::abs(cfg2.species[0].weight - (4.0 / 12.0)) < 1e-15, "2D density-derived macro-particle weight is wrong");
            std::filesystem::remove(config_2d_path);
            require_throws([] {
                const auto path = std::filesystem::path("test_invalid_2d_electrode.ini");
                { std::ofstream out(path); out << "dimension = 2\nphi_left = inf\n"; }
                try { (void)pic::load_config_2d(path.string()); } catch (...) { std::filesystem::remove(path); throw; }
                std::filesystem::remove(path);
            }, "invalid 2D electrode potential validation did not throw");
            require_throws([] {
                const auto path = std::filesystem::path("test_empty_2d_boundary_tag.ini");
                { std::ofstream out(path); out << "dimension = 2\nboundary_left_tag = \n"; }
                try { (void)pic::load_config_2d(path.string()); } catch (...) { std::filesystem::remove(path); throw; }
                std::filesystem::remove(path);
            }, "empty 2D boundary tag validation did not throw");
            require_throws([] {
                const auto path = std::filesystem::path("test_invalid_2d_particle_stride.ini");
                { std::ofstream out(path); out << "dimension = 2\nparticle_output_stride = 0\n"; }
                try { (void)pic::load_config_2d(path.string()); } catch (...) { std::filesystem::remove(path); throw; }
                std::filesystem::remove(path);
            }, "invalid 2D particle output stride validation did not throw");
            require_throws([] {
                const auto path = std::filesystem::path("test_invalid_2d_particle_boundary.ini");
                { std::ofstream out(path); out << "dimension = 2\nparticle_boundary_left = bounce\n"; }
                try { (void)pic::load_config_2d(path.string()); } catch (...) { std::filesystem::remove(path); throw; }
                std::filesystem::remove(path);
            }, "invalid 2D particle boundary validation did not throw");
            require_throws([] {
                const auto path = std::filesystem::path("test_invalid_2d_magnetic_field.ini");
                { std::ofstream out(path); out << "dimension = 2\nmagnetic_field_z = nan\n"; }
                try { (void)pic::load_config_2d(path.string()); } catch (...) { std::filesystem::remove(path); throw; }
                std::filesystem::remove(path);
            }, "invalid 2D magnetic field validation did not throw");
            require_throws([] {
                const auto path = std::filesystem::path("test_invalid_key.ini");
                { std::ofstream out(path); out << "nx = 16\nunknown = 7\n"; }
                try { (void)pic::load_config(path.string()); } catch (...) { std::filesystem::remove(path); throw; }
                std::filesystem::remove(path);
            }, "unknown config key validation did not throw");
            require_throws([] {
                const auto path = std::filesystem::path("test_missing_2d_dimension.ini");
                { std::ofstream out(path); out << "nx = 8\nny = 8\nlength_x = 1\nlength_y = 1\n"; }
                try { (void)pic::load_config_2d(path.string()); } catch (...) { std::filesystem::remove(path); throw; }
                std::filesystem::remove(path);
            }, "missing 2D dimension validation did not throw");
            require_throws([] {
                const auto path = std::filesystem::path("test_invalid_boundary.ini");
                { std::ofstream out(path); out << "boundary = reflective\n"; }
                try { (void)pic::load_config(path.string()); } catch (...) { std::filesystem::remove(path); throw; }
                std::filesystem::remove(path);
            }, "invalid boundary validation did not throw");
            require_throws([] {
                const auto path = std::filesystem::path("test_invalid_bool.ini");
                { std::ofstream out(path); out << "[collisions]\nenabled = maybe\n"; }
                try { (void)pic::load_config(path.string()); } catch (...) { std::filesystem::remove(path); throw; }
                std::filesystem::remove(path);
            }, "invalid boolean validation did not throw");
        {
            const auto config_3d_path = std::filesystem::path("test_config_3d.ini");
            {
                std::ofstream out(config_3d_path);
                out << "dimension = 3\n"
                    << "nx = 6\n"
                    << "ny = 5\n"
                    << "nz = 4\n"
                    << "length_x = 2.0\n"
                    << "length_y = 1.5\n"
                    << "length_z = 1.25\n"
                    << "dt = 0.01\n"
                    << "steps = 2\n"
                    << "output_interval = 1\n"
                    << "output_dir = test_output_config_3d\n"
                    << "vtk_output = true\n"
                    << "particle_output = true\n"
                    << "particle_output_interval = 3\n"
                    << "particle_output_stride = 2\n"
                    << "particle_sample_count = 7\n"
                    << "particle_boundary = reflecting\n"
                    << "particle_boundary_left = absorbing\n"
                    << "particle_boundary_right = periodic\n"
                    << "particle_boundary_bottom = auto\n"
                    << "particle_boundary_top = reflect\n"
                    << "particle_boundary_back = absorbing\n"
                    << "particle_boundary_front = periodic\n"
                    << "boundary = periodic\n"
                    << "magnetic_field_x = 0.25\n"
                    << "magnetic_field_y = -0.5\n"
                    << "magnetic_field_z = 1.75\n"
                    << "[species.electrons]\n"
                    << "charge = -1\n"
                    << "mass = 1\n"
                    << "density = 4\n"
                    << "particles = 12\n"
                    << "thermal_velocity = 0\n"
                    << "init_x_min = 0.5\n"
                    << "init_x_max = 1.5\n"
                    << "init_y_min = 0.25\n"
                    << "init_y_max = 1.25\n"
                    << "init_z_min = 0.125\n"
                    << "init_z_max = 0.625\n";
            }
            require(pic::detect_config_dimension(config_3d_path.string()) == 3, "3D config dimension was not detected");
            auto cfg3 = pic::load_config_3d(config_3d_path.string());
            require(cfg3.nx == 6 && cfg3.ny == 5 && cfg3.nz == 4, "3D config did not load mesh dimensions");
            require(cfg3.vtk_output, "3D config did not load vtk_output");
            require(cfg3.particle_output, "3D config did not load particle_output");
            require(cfg3.particle_output_interval == 3, "3D config did not load particle_output_interval");
            require(cfg3.particle_output_stride == 2, "3D config did not load particle_output_stride");
            require(cfg3.particle_sample_count == 7, "3D config did not load particle_sample_count");
            require(cfg3.particle_boundary_config.left == pic::ParticleBoundary::Absorbing, "3D config did not load left particle boundary");
            require(cfg3.particle_boundary_config.right == pic::ParticleBoundary::Periodic, "3D config did not load right particle boundary");
            require(cfg3.particle_boundary_config.bottom == pic::ParticleBoundary::Auto, "3D config did not load bottom particle boundary");
            require(cfg3.particle_boundary_config.top == pic::ParticleBoundary::Reflecting, "3D config did not load top particle boundary");
            require(cfg3.particle_boundary_config.back == pic::ParticleBoundary::Absorbing, "3D config did not load back particle boundary");
            require(cfg3.particle_boundary_config.front == pic::ParticleBoundary::Periodic, "3D config did not load front particle boundary");
            require(std::abs(cfg3.magnetic_field.x - 0.25) < 1e-15, "3D config did not load magnetic_field_x");
            require(std::abs(cfg3.magnetic_field.y + 0.5) < 1e-15, "3D config did not load magnetic_field_y");
            require(std::abs(cfg3.magnetic_field.z - 1.75) < 1e-15, "3D config did not load magnetic_field_z");
            require(cfg3.species.size() == 1, "3D config did not load one species");
            require(std::abs(cfg3.species[0].weight - (4.0 * 1.0 * 1.0 * 0.5 / 12.0)) < 1e-15,
                    "3D density-derived macro-particle weight is wrong");
            std::filesystem::remove(config_3d_path);
            require_throws([] {
                const auto path = std::filesystem::path("test_missing_3d_dimension.ini");
                { std::ofstream out(path); out << "nx = 4\nny = 4\nnz = 4\nlength_x = 1\nlength_y = 1\nlength_z = 1\n"; }
                try { (void)pic::load_config_3d(path.string()); } catch (...) { std::filesystem::remove(path); throw; }
                std::filesystem::remove(path);
            }, "missing 3D dimension validation did not throw");
            require_throws([] {
                const auto path = std::filesystem::path("test_invalid_3d_magnetic_field.ini");
                { std::ofstream out(path); out << "dimension = 3\nmagnetic_field_y = inf\n"; }
                try { (void)pic::load_config_3d(path.string()); } catch (...) { std::filesystem::remove(path); throw; }
                std::filesystem::remove(path);
            }, "invalid 3D magnetic field validation did not throw");
            require_throws([] {
                const auto path = std::filesystem::path("test_invalid_3d_particle_boundary.ini");
                { std::ofstream out(path); out << "dimension = 3\nparticle_boundary_front = bounce\n"; }
                try { (void)pic::load_config_3d(path.string()); } catch (...) { std::filesystem::remove(path); throw; }
                std::filesystem::remove(path);
            }, "invalid 3D particle boundary validation did not throw");
        }
        }
        {
            const auto output_dir = std::filesystem::path("test_output");
            std::filesystem::remove_all(output_dir);

            pic::Config cfg;
            cfg.nx = 32;
            cfg.length = 1.0;
            cfg.dt = 0.005;
            cfg.steps = 4;
            cfg.output_interval = 2;
            cfg.output_dir = output_dir.string();
            cfg.species = {pic::SpeciesConfig{"e", -1.0, 1.0, 1.0, 200, 1.0, 0.0, 0.01, 0.0, -1.0},
                           pic::SpeciesConfig{"i", 1.0, 1836.0, 1.0, 200, 1.0, 0.0, 0.0, 0.0, -1.0}};
            pic::Simulation sim(cfg);
            auto s = sim.run();
            require(s.steps_completed == cfg.steps, "simulation did not complete requested transient steps");
            require(std::filesystem::exists(output_dir / "scalars.csv"), "simulation did not write scalar diagnostics");
            require(std::filesystem::exists(output_dir / "fields_4.csv"), "simulation did not write final field diagnostics");
        }
        {
            pic::Simulation2DConfig cfg;
            cfg.nx = 16;
            cfg.ny = 12;
            cfg.length_x = 1.0;
            cfg.length_y = 0.75;
            cfg.dt = 0.002;
            cfg.steps = 3;
            cfg.boundary = pic::Boundary::Periodic;
            cfg.seed = 7;
            cfg.vtk_output = true;
            cfg.particle_output = true;
            cfg.output_interval = 2;
            cfg.particle_output_interval = 2;
            cfg.particle_output_stride = 3;
            cfg.particle_sample_count = 5;
            cfg.output_dir = "test_output_2d";
            std::filesystem::remove_all(cfg.output_dir);
            cfg.species = {pic::Species2DConfig{"e2", -1.0, 1.0, 0.01, 128, 0.02, -0.01, 0.001, 0.0, -1.0, 0.0, -1.0},
                           pic::Species2DConfig{"i2", 1.0, 1836.0, 0.01, 128, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, -1.0}};
            pic::Simulation2D sim(cfg);
            auto s = sim.run();
            require(s.steps_completed == cfg.steps, "2D simulation did not complete requested transient steps");
            require(std::abs(s.final_time - cfg.dt * static_cast<double>(cfg.steps)) < 1e-15, "2D simulation final time is wrong");
            require(s.final_sample.live_particles == 256, "2D periodic simulation lost live particles");
            require(s.final_sample.boundary_losses.absorbed_left == 0 && s.final_sample.boundary_losses.absorbed_right == 0 &&
                    s.final_sample.boundary_losses.absorbed_bottom == 0 && s.final_sample.boundary_losses.absorbed_top == 0,
                    "2D periodic simulation reported unexpected boundary losses");

            double total_charge = 0.0;
            for (double rho : sim.mesh().rho()) total_charge += rho * sim.mesh().dx() * sim.mesh().dy();
            require(std::abs(total_charge) < 1e-12, "2D simulation did not conserve net neutral charge");
            require(std::filesystem::exists(cfg.output_dir / "fields_0.vtk"), "2D simulation did not write initial VTK fields");
            require(std::filesystem::exists(cfg.output_dir / "fields_2.vtk"), "2D simulation did not write interval VTK fields");
            require(std::filesystem::exists(cfg.output_dir / "fields_3.vtk"), "2D simulation did not write final VTK fields");
            require(!std::filesystem::exists(cfg.output_dir / "fields_1.vtk"), "2D simulation wrote an unexpected VTK interval");
            require(std::filesystem::exists(cfg.output_dir / "scalars.csv"), "2D simulation did not write scalar diagnostics");
            const auto scalars = read_file_text(cfg.output_dir / "scalars.csv");
            require(scalars.find("step,time,kinetic_energy,field_energy,total_energy,charge_l1,live_particles,absorbed_left,absorbed_right,absorbed_bottom,absorbed_top,live_particles_e2,live_particles_i2\n") == 0,
                    "2D scalar diagnostics header is wrong");
            require(count_lines(scalars) == 4, "2D scalar diagnostics wrote unexpected number of rows");
            require(std::filesystem::exists(cfg.output_dir / "particles_0.csv"), "2D simulation did not write initial particle sample");
            require(std::filesystem::exists(cfg.output_dir / "particles_2.csv"), "2D simulation did not write interval particle sample");
            require(std::filesystem::exists(cfg.output_dir / "particles_3.csv"), "2D simulation did not write final particle sample");
            require(!std::filesystem::exists(cfg.output_dir / "particles_1.csv"), "2D simulation wrote an unexpected particle sample interval");
            const auto particles = read_file_text(cfg.output_dir / "particles_0.csv");
            require(particles.find("species_id,species,x,y,vx,vy,alive\n") == 0, "2D particle diagnostics header is wrong");
            require(count_lines(particles) == 6, "2D particle diagnostics did not honor sample_count");
            const auto vtk = read_file_text(cfg.output_dir / "fields_3.vtk");
            require(vtk.find("DIMENSIONS 16 12 1") != std::string::npos, "2D simulation VTK dimensions are wrong");
        }
        {
            pic::Simulation3DConfig cfg;
            cfg.nx = 8;
            cfg.ny = 6;
            cfg.nz = 5;
            cfg.length_x = 1.0;
            cfg.length_y = 0.75;
            cfg.length_z = 0.5;
            cfg.dt = 0.001;
            cfg.steps = 3;
            cfg.boundary = pic::Boundary::Periodic;
            cfg.seed = 11;
            cfg.vtk_output = true;
            cfg.particle_output = true;
            cfg.output_interval = 2;
            cfg.particle_output_interval = 2;
            cfg.particle_output_stride = 2;
            cfg.particle_sample_count = 5;
            cfg.output_dir = "test_output_3d";
            std::filesystem::remove_all(cfg.output_dir);
            cfg.species = {pic::Species3DConfig{"e3", -1.0, 1.0, 0.01, 64, 0.01, -0.005, 0.002, 0.001, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0},
                           pic::Species3DConfig{"i3", 1.0, 1836.0, 0.01, 64, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0}};
            pic::Simulation3D sim(cfg);
            auto s = sim.run();
            require(s.steps_completed == cfg.steps, "3D simulation did not complete requested transient steps");
            require(std::abs(s.final_time - cfg.dt * static_cast<double>(cfg.steps)) < 1e-15, "3D simulation final time is wrong");
            require(s.final_sample.live_particles == 128, "3D periodic simulation lost live particles");
            require(s.final_sample.boundary_losses.absorbed_left == 0 && s.final_sample.boundary_losses.absorbed_right == 0 &&
                    s.final_sample.boundary_losses.absorbed_bottom == 0 && s.final_sample.boundary_losses.absorbed_top == 0 &&
                    s.final_sample.boundary_losses.absorbed_back == 0 && s.final_sample.boundary_losses.absorbed_front == 0,
                    "3D periodic simulation reported unexpected boundary losses");

            double total_charge = 0.0;
            for (std::size_t k = 0; k < sim.mesh().nz(); ++k) {
                for (std::size_t j = 0; j < sim.mesh().ny(); ++j) {
                    for (std::size_t i = 0; i < sim.mesh().nx(); ++i) {
                        total_charge += sim.mesh().rho()[sim.mesh().index(i, j, k)] * sim.mesh().node_volume(i, j, k);
                    }
                }
            }
            require(std::abs(total_charge) < 1e-12, "3D simulation did not conserve net neutral charge");
            require(std::filesystem::exists(cfg.output_dir / "fields_0.vtk"), "3D simulation did not write initial VTK fields");
            require(std::filesystem::exists(cfg.output_dir / "fields_2.vtk"), "3D simulation did not write interval VTK fields");
            require(std::filesystem::exists(cfg.output_dir / "fields_3.vtk"), "3D simulation did not write final VTK fields");
            require(!std::filesystem::exists(cfg.output_dir / "fields_1.vtk"), "3D simulation wrote an unexpected VTK interval");
            require(std::filesystem::exists(cfg.output_dir / "scalars.csv"), "3D simulation did not write scalar diagnostics");
            const auto scalars3d = read_file_text(cfg.output_dir / "scalars.csv");
            require(scalars3d.find("step,time,kinetic_energy,field_energy,total_energy,charge_l1,live_particles,absorbed_left,absorbed_right,absorbed_bottom,absorbed_top,absorbed_back,absorbed_front,live_particles_e3,live_particles_i3\n") == 0,
                    "3D scalar diagnostics header is wrong");
            require(count_lines(scalars3d) == 4, "3D scalar diagnostics wrote unexpected number of rows");
            require(std::filesystem::exists(cfg.output_dir / "particles_0.csv"), "3D simulation did not write initial particle sample");
            require(std::filesystem::exists(cfg.output_dir / "particles_2.csv"), "3D simulation did not write interval particle sample");
            require(std::filesystem::exists(cfg.output_dir / "particles_3.csv"), "3D simulation did not write final particle sample");
            require(!std::filesystem::exists(cfg.output_dir / "particles_1.csv"), "3D simulation wrote an unexpected particle sample interval");
            const auto particles3d = read_file_text(cfg.output_dir / "particles_0.csv");
            require(particles3d.find("species_id,species,x,y,z,vx,vy,vz,alive\n") == 0, "3D particle diagnostics header is wrong");
            require(count_lines(particles3d) == 6, "3D particle diagnostics did not honor sample_count");
            const auto vtk3d = read_file_text(cfg.output_dir / "fields_3.vtk");
            require(vtk3d.find("DIMENSIONS 8 6 5") != std::string::npos, "3D simulation VTK dimensions are wrong");
        }
        {
            pic::Simulation2DConfig cfg;
            cfg.nx = 8;
            cfg.ny = 8;
            cfg.length_x = 1.0;
            cfg.length_y = 1.0;
            cfg.dt = 0.1;
            cfg.steps = 1;
            cfg.boundary = pic::Boundary::Dirichlet;
            cfg.output_dir = "test_output_2d_absorb";
            cfg.species = {pic::Species2DConfig{"tracer", 0.0, 1.0, 1.0, 4, -10.0, 0.0, 0.0, 0.49, 0.51, 0.45, 0.55}};
            std::filesystem::remove_all(cfg.output_dir);
            pic::Simulation2D sim(cfg);
            auto s = sim.run();
            require(s.final_sample.live_particles == 0, "2D absorbing particle boundary did not remove escaped particles");
            require(s.final_sample.boundary_losses.absorbed_left == 4, "2D absorbing particle boundary did not count left losses");
            require(sim.boundary_losses().absorbed_left == 4, "2D absorbing boundary loss accessor is wrong");
        }
        {
            pic::Simulation2DConfig cfg;
            cfg.nx = 8;
            cfg.ny = 8;
            cfg.length_x = 1.0;
            cfg.length_y = 1.0;
            cfg.dt = 0.1;
            cfg.steps = 1;
            cfg.boundary = pic::Boundary::Dirichlet;
            cfg.output_dir = "test_output_2d_absorb_bottom";
            cfg.species = {pic::Species2DConfig{"tracer", 0.0, 1.0, 1.0, 4, 0.0, -10.0, 0.0, 0.45, 0.55, 0.49, 0.51}};
            std::filesystem::remove_all(cfg.output_dir);
            pic::Simulation2D sim(cfg);
            auto s = sim.run();
            require(s.final_sample.live_particles == 0, "2D bottom absorbing particle boundary did not remove escaped particles");
            require(s.final_sample.boundary_losses.absorbed_bottom == 4, "2D bottom absorbing particle boundary did not count bottom losses");
            require(sim.boundary_losses().absorbed_bottom == 4, "2D bottom absorbing boundary loss accessor is wrong");
        }
        {
            pic::Simulation2DConfig cfg;
            cfg.nx = 8;
            cfg.ny = 8;
            cfg.length_x = 1.0;
            cfg.length_y = 1.0;
            cfg.dt = 0.1;
            cfg.steps = 1;
            cfg.boundary = pic::Boundary::Dirichlet;
            cfg.output_dir = "test_output_2d_reflect";
            cfg.particle_boundary_config.left = pic::ParticleBoundary::Reflecting;
            cfg.species = {pic::Species2DConfig{"tracer", 0.0, 1.0, 1.0, 4, -10.0, 0.0, 0.0, 0.49, 0.51, 0.45, 0.55}};
            std::filesystem::remove_all(cfg.output_dir);
            pic::Simulation2D sim(cfg);
            auto s = sim.run();
            require(s.final_sample.live_particles == 4, "2D reflecting particle boundary lost particles");
            require(s.final_sample.boundary_losses.absorbed_left == 0, "2D reflecting particle boundary reported absorption");
            for (const auto& particle : sim.species().front().particles()) {
                require(particle.position.x >= 0.0 && particle.position.x <= cfg.length_x, "2D reflecting boundary left particle outside domain");
                require(particle.velocity.x > 0.0, "2D reflecting boundary did not reverse x velocity");
            }
        }
        {
            pic::Simulation2DConfig cfg;
            cfg.nx = 8;
            cfg.ny = 8;
            cfg.length_x = 1.0;
            cfg.length_y = 1.0;
            cfg.dt = 0.1;
            cfg.steps = 1;
            cfg.boundary = pic::Boundary::Dirichlet;
            cfg.output_dir = "test_output_2d_side_policies";
            cfg.particle_boundary_config.left = pic::ParticleBoundary::Periodic;
            cfg.particle_boundary_config.right = pic::ParticleBoundary::Periodic;
            cfg.particle_boundary_config.bottom = pic::ParticleBoundary::Absorbing;
            cfg.particle_boundary_config.top = pic::ParticleBoundary::Reflecting;
            cfg.species = {pic::Species2DConfig{"tracer", 0.0, 1.0, 1.0, 4, -10.0, 10.0, 0.0, 0.49, 0.51, 0.49, 0.51}};
            std::filesystem::remove_all(cfg.output_dir);
            pic::Simulation2D sim(cfg);
            auto s = sim.run();
            require(s.final_sample.live_particles == 4, "2D periodic/reflecting side policies lost particles");
            require(s.final_sample.boundary_losses.absorbed_left == 0 && s.final_sample.boundary_losses.absorbed_right == 0 &&
                    s.final_sample.boundary_losses.absorbed_bottom == 0 && s.final_sample.boundary_losses.absorbed_top == 0,
                    "2D periodic/reflecting side policies reported absorption");
            for (const auto& particle : sim.species().front().particles()) {
                require(particle.position.x >= 0.0 && particle.position.x <= cfg.length_x, "2D periodic boundary left particle outside domain");
                require(particle.position.y >= 0.0 && particle.position.y <= cfg.length_y, "2D top reflecting boundary left particle outside domain");
                require(particle.velocity.y < 0.0, "2D top reflecting boundary did not reverse y velocity");
            }
        }
        {
            const auto output_dir = std::filesystem::path("test_output_checkpoint_1d");
            const auto checkpoint_path = output_dir / "manual.apc";
            std::filesystem::remove_all(output_dir);

            pic::Config cfg;
            cfg.nx = 24;
            cfg.length = 1.0;
            cfg.dt = 0.003;
            cfg.steps = 5;
            cfg.output_interval = 5;
            cfg.output_dir = output_dir.string();
            cfg.seed = 99;
            cfg.collisions.enabled = true;
            cfg.collisions.frequency = 200.0;
            cfg.collisions.neutral_temperature_velocity = 0.02;
            cfg.species = {pic::SpeciesConfig{"e", -1.0, 1.0, 0.02, 48, 1.0, 0.03, 0.01, 0.0, -1.0},
                           pic::SpeciesConfig{"i", 1.0, 1836.0, 0.02, 48, 1.0, -0.01, 0.0, 0.0, -1.0}};

            pic::Simulation continuous(cfg);
            continuous.initialize();
            continuous.step();
            continuous.step();
            continuous.save_checkpoint(checkpoint_path);
            require(std::filesystem::exists(checkpoint_path), "1D checkpoint file was not written");
            for (std::size_t n = continuous.step_count(); n < cfg.steps; ++n) continuous.step();

            pic::Simulation restarted(cfg);
            restarted.load_checkpoint(checkpoint_path);
            require(restarted.step_count() == 2, "1D checkpoint did not restore step count");
            require_near(restarted.time(), 2.0 * cfg.dt, 1e-15, "1D checkpoint did not restore time");
            for (std::size_t n = restarted.step_count(); n < cfg.steps; ++n) restarted.step();

            require_checkpoint_samples_close(continuous.sample(), restarted.sample(), "1D checkpoint restart");
            require_species_close(continuous.species(), restarted.species(), "1D checkpoint restart");

            require_throws([&] {
                pic::Simulation2DConfig bad_cfg;
                bad_cfg.species = {pic::Species2DConfig{"e", -1.0, 1.0, 0.02, 48, 0.03, 0.0, 0.01, 0.0, -1.0, 0.0, -1.0},
                                   pic::Species2DConfig{"i", 1.0, 1836.0, 0.02, 48, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, -1.0}};
                pic::Simulation2D bad(bad_cfg);
                bad.load_checkpoint(checkpoint_path);
            }, "loading 1D checkpoint into 2D simulation did not throw");
        }
        {
            const auto output_dir = std::filesystem::path("test_output_checkpoint_2d");
            const auto checkpoint_path = output_dir / "manual.apc";
            std::filesystem::remove_all(output_dir);

            pic::Simulation2DConfig cfg;
            cfg.nx = 10;
            cfg.ny = 8;
            cfg.length_x = 1.0;
            cfg.length_y = 0.8;
            cfg.dt = 0.002;
            cfg.steps = 5;
            cfg.output_interval = 5;
            cfg.output_dir = output_dir;
            cfg.seed = 123;
            cfg.boundary = pic::Boundary::Periodic;
            cfg.species = {pic::Species2DConfig{"e2", -1.0, 1.0, 0.02, 48, 0.03, -0.02, 0.01, 0.0, -1.0, 0.0, -1.0},
                           pic::Species2DConfig{"i2", 1.0, 1836.0, 0.02, 48, -0.01, 0.0, 0.0, 0.0, -1.0, 0.0, -1.0}};

            pic::Simulation2D continuous(cfg);
            continuous.initialize();
            continuous.step();
            continuous.step();
            continuous.save_checkpoint(checkpoint_path);
            require(std::filesystem::exists(checkpoint_path), "2D checkpoint file was not written");
            for (std::size_t n = continuous.step_count(); n < cfg.steps; ++n) continuous.step();

            pic::Simulation2D restarted(cfg);
            restarted.load_checkpoint(checkpoint_path);
            require(restarted.step_count() == 2, "2D checkpoint did not restore step count");
            require_near(restarted.time(), 2.0 * cfg.dt, 1e-15, "2D checkpoint did not restore time");
            for (std::size_t n = restarted.step_count(); n < cfg.steps; ++n) restarted.step();

            require_checkpoint_samples_close(continuous.sample(), restarted.sample(), "2D checkpoint restart");
            require_species_close(continuous.species(), restarted.species(), "2D checkpoint restart");

            auto run_cfg = cfg;
            run_cfg.output_dir = output_dir / "run";
            run_cfg.steps = 3;
            run_cfg.output_interval = 2;
            run_cfg.checkpoint_output = true;
            run_cfg.checkpoint_interval = 2;
            std::filesystem::remove_all(run_cfg.output_dir);
            pic::Simulation2D run_sim(run_cfg);
            auto summary = run_sim.run();
            require(summary.steps_completed == 3, "2D checkpoint run did not complete");
            require(std::filesystem::exists(run_cfg.output_dir / "checkpoint_0.apc"), "2D run did not write initial checkpoint");
            require(std::filesystem::exists(run_cfg.output_dir / "checkpoint_2.apc"), "2D run did not write interval checkpoint");
            require(std::filesystem::exists(run_cfg.output_dir / "checkpoint_3.apc"), "2D run did not write final checkpoint");
        }
        {
            const auto output_dir = std::filesystem::path("test_output_checkpoint_3d");
            const auto checkpoint_path = output_dir / "manual.apc";
            std::filesystem::remove_all(output_dir);

            pic::Simulation3DConfig cfg;
            cfg.nx = 6;
            cfg.ny = 5;
            cfg.nz = 4;
            cfg.length_x = 1.0;
            cfg.length_y = 0.8;
            cfg.length_z = 0.6;
            cfg.dt = 0.0015;
            cfg.steps = 5;
            cfg.output_interval = 5;
            cfg.output_dir = output_dir;
            cfg.seed = 321;
            cfg.boundary = pic::Boundary::Periodic;
            cfg.species = {pic::Species3DConfig{"e3", -1.0, 1.0, 0.02, 36, 0.02, -0.01, 0.005, 0.01, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0},
                           pic::Species3DConfig{"i3", 1.0, 1836.0, 0.02, 36, -0.005, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0}};

            pic::Simulation3D continuous(cfg);
            continuous.initialize();
            continuous.step();
            continuous.step();
            continuous.save_checkpoint(checkpoint_path);
            require(std::filesystem::exists(checkpoint_path), "3D checkpoint file was not written");
            for (std::size_t n = continuous.step_count(); n < cfg.steps; ++n) continuous.step();

            pic::Simulation3D restarted(cfg);
            restarted.load_checkpoint(checkpoint_path);
            require(restarted.step_count() == 2, "3D checkpoint did not restore step count");
            require_near(restarted.time(), 2.0 * cfg.dt, 1e-15, "3D checkpoint did not restore time");
            for (std::size_t n = restarted.step_count(); n < cfg.steps; ++n) restarted.step();

            require_checkpoint_samples_close(continuous.sample(), restarted.sample(), "3D checkpoint restart");
            require_species_close(continuous.species(), restarted.species(), "3D checkpoint restart");

            auto run_cfg = cfg;
            run_cfg.output_dir = output_dir / "run";
            run_cfg.steps = 3;
            run_cfg.output_interval = 2;
            run_cfg.checkpoint_output = true;
            run_cfg.checkpoint_interval = 2;
            std::filesystem::remove_all(run_cfg.output_dir);
            pic::Simulation3D run_sim(run_cfg);
            auto summary = run_sim.run();
            require(summary.steps_completed == 3, "3D checkpoint run did not complete");
            require(std::filesystem::exists(run_cfg.output_dir / "checkpoint_0.apc"), "3D run did not write initial checkpoint");
            require(std::filesystem::exists(run_cfg.output_dir / "checkpoint_2.apc"), "3D run did not write interval checkpoint");
            require(std::filesystem::exists(run_cfg.output_dir / "checkpoint_3.apc"), "3D run did not write final checkpoint");
        }
        {
            const double qm = 1.0;
            const double magnetic_z = 2.0;
            const double dt = 0.1;
            const double angle = -2.0 * std::atan(0.5 * qm * magnetic_z * dt);

            pic::Particle2D p2;
            p2.velocity_half = pic::Vec2{1.0, 0.0};
            pic::kick_boris(p2, pic::Vec2{0.0, 0.0}, magnetic_z, qm, dt);
            require_near(p2.velocity_half.x, std::cos(angle), 1e-15, "2D Boris cyclotron vx mismatch");
            require_near(p2.velocity_half.y, std::sin(angle), 1e-15, "2D Boris cyclotron vy mismatch");
            require_near(p2.velocity_half.x * p2.velocity_half.x + p2.velocity_half.y * p2.velocity_half.y,
                         1.0, 1e-15, "2D Boris cyclotron speed was not conserved");

            pic::Particle3D p3;
            p3.velocity_half = pic::Vec3{1.0, 0.0, 0.25};
            pic::kick_boris(p3, pic::Vec3{0.0, 0.0, 0.0}, pic::Vec3{0.0, 0.0, magnetic_z}, qm, dt);
            require_near(p3.velocity_half.x, std::cos(angle), 1e-15, "3D Boris cyclotron vx mismatch");
            require_near(p3.velocity_half.y, std::sin(angle), 1e-15, "3D Boris cyclotron vy mismatch");
            require_near(p3.velocity_half.z, 0.25, 1e-15, "3D Boris cyclotron parallel velocity changed");
            require_near(p3.velocity_half.x * p3.velocity_half.x + p3.velocity_half.y * p3.velocity_half.y,
                         1.0, 1e-15, "3D Boris cyclotron perpendicular speed was not conserved");
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "aurorapic core test failure: " << e.what() << '\n';
        return 1;
    }
}
