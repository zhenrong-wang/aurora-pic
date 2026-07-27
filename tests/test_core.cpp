#include "pic/Config.hpp"
#include "pic/Diagnostics.hpp"
#include "pic/FieldSolver.hpp"
#include "pic/Grid.hpp"
#include "pic/ImportedMesh2D.hpp"
#include "pic/Mesh2D.hpp"
#include "pic/Mesh3D.hpp"
#include "pic/Simulation.hpp"
#include "pic/Simulation2D.hpp"
#include "pic/Simulation3D.hpp"
#include "pic/Pusher.hpp"
#include "pic/Runtime.hpp"
#include "pic/VTKWriter.hpp"
#include "pic/Species.hpp"
#include "pic/Species2D.hpp"
#include "pic/Species3D.hpp"
#include "pic/UnstructuredFieldSolver2D.hpp"
#include "pic/UnstructuredMesh2D.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
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

template <typename Fn>
void require_throws_contains(Fn&& fn, const std::string& expected, const std::string& message) {
    try {
        fn();
    } catch (const std::exception& exc) {
        const std::string actual = exc.what();
        if (actual.find(expected) != std::string::npos) return;
        throw std::runtime_error(message + ": expected substring '" + expected + "' in '" + actual + "'");
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
            require(max_err < 1e-12, "M1 1D periodic Poisson sine benchmark exceeded electric-field tolerance");
        }
        {
            constexpr std::size_t nx = 16;
            constexpr std::size_t ny = 12;
            constexpr double length_x = 1.0;
            constexpr double length_y = 1.5;
            constexpr double kx = 2.0 * std::numbers::pi * 2.0 / length_x;
            constexpr double ky = 2.0 * std::numbers::pi * 3.0 / length_y;
            constexpr double k2 = kx * kx + ky * ky;
            pic::Mesh2D mesh(nx, ny, length_x, length_y, pic::Boundary::Periodic);
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
                    const double expected_phi = std::sin(kx * x) * std::cos(ky * y) / k2;
                    const double expected_ex = -kx * std::cos(kx * x) * std::cos(ky * y) / k2;
                    const double expected_ey = ky * std::sin(kx * x) * std::sin(ky * y) / k2;
                    const auto idx = mesh.index(i, j);
                    max_phi_err = std::max(max_phi_err, std::abs(mesh.phi()[idx] - expected_phi));
                    max_ex_err = std::max(max_ex_err, std::abs(mesh.electric_x()[idx] - expected_ex));
                    max_ey_err = std::max(max_ey_err, std::abs(mesh.electric_y()[idx] - expected_ey));
                }
            }
            require(max_phi_err < 1e-12, "M1 2D periodic Poisson sine-cosine benchmark exceeded potential tolerance");
            require(max_ex_err < 1e-12, "M1 2D periodic Poisson sine-cosine benchmark exceeded Ex tolerance");
            require(max_ey_err < 1e-12, "M1 2D periodic Poisson sine-cosine benchmark exceeded Ey tolerance");
        }
        {
            constexpr std::size_t nx = 8;
            constexpr std::size_t ny = 10;
            constexpr std::size_t nz = 12;
            constexpr double length_x = 1.0;
            constexpr double length_y = 1.25;
            constexpr double length_z = 1.5;
            constexpr double kx = 2.0 * std::numbers::pi * 1.0 / length_x;
            constexpr double ky = 2.0 * std::numbers::pi * 2.0 / length_y;
            constexpr double kz = 2.0 * std::numbers::pi * 3.0 / length_z;
            constexpr double k2 = kx * kx + ky * ky + kz * kz;
            pic::Mesh3D mesh(nx, ny, nz, length_x, length_y, length_z, pic::Boundary::Periodic);
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
                        const double expected_phi = std::sin(kx * x) * std::cos(ky * y) * std::sin(kz * z) / k2;
                        const double expected_ex = -kx * std::cos(kx * x) * std::cos(ky * y) * std::sin(kz * z) / k2;
                        const double expected_ey = ky * std::sin(kx * x) * std::sin(ky * y) * std::sin(kz * z) / k2;
                        const double expected_ez = -kz * std::sin(kx * x) * std::cos(ky * y) * std::cos(kz * z) / k2;
                        const auto idx = mesh.index(i, j, k);
                        max_phi_err = std::max(max_phi_err, std::abs(mesh.phi()[idx] - expected_phi));
                        max_ex_err = std::max(max_ex_err, std::abs(mesh.electric_x()[idx] - expected_ex));
                        max_ey_err = std::max(max_ey_err, std::abs(mesh.electric_y()[idx] - expected_ey));
                        max_ez_err = std::max(max_ez_err, std::abs(mesh.electric_z()[idx] - expected_ez));
                    }
                }
            }
            require(max_phi_err < 1e-12, "M1 3D periodic Poisson sine-cosine benchmark exceeded potential tolerance");
            require(max_ex_err < 1e-12, "M1 3D periodic Poisson sine-cosine benchmark exceeded Ex tolerance");
            require(max_ey_err < 1e-12, "M1 3D periodic Poisson sine-cosine benchmark exceeded Ey tolerance");
            require(max_ez_err < 1e-12, "M1 3D periodic Poisson sine-cosine benchmark exceeded Ez tolerance");
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
            // M5 regression: Simulation2D config-level uniform-B activation should use the
            // same Boris rotation as the pusher-level cyclotron benchmark while keeping the
            // field solve electrostatic. A tiny macro-particle weight suppresses self-field
            // feedback so the prescribed-B rotation is deterministic end to end.
            constexpr double dt = 0.04;
            constexpr std::size_t steps = 37;
            constexpr double charge_to_mass = 1.25;
            constexpr double magnetic_z = 1.7;
            const pic::Vec2 initial_velocity{0.8, -0.35};

            pic::Simulation2DConfig cfg;
            cfg.nx = 8;
            cfg.ny = 8;
            cfg.length_x = 1.0;
            cfg.length_y = 1.0;
            cfg.dt = dt;
            cfg.steps = steps;
            cfg.boundary = pic::Boundary::Periodic;
            cfg.output_interval = steps;
            cfg.magnetic_field_z = magnetic_z;
            cfg.species = {pic::Species2DConfig{"m5_boris_ion", charge_to_mass, 1.0, 1e-30, 1,
                                                initial_velocity.x, initial_velocity.y, 0.0,
                                                0.25, 0.26, 0.25, 0.26}};

            pic::Simulation2D sim(cfg);
            sim.initialize();
            require_near(sim.species()[0].particles()[0].velocity.x, initial_velocity.x, 1e-12,
                         "M5 2D uniform-B simulation changed initial vx during initialization");
            require_near(sim.species()[0].particles()[0].velocity.y, initial_velocity.y, 1e-12,
                         "M5 2D uniform-B simulation changed initial vy during initialization");
            sim.run();

            const auto& particle = sim.species()[0].particles()[0];
            const double angle = static_cast<double>(steps) * boris_rotation_angle(std::abs(magnetic_z), charge_to_mass, dt);
            const double expected_x = initial_velocity.x * std::cos(angle) - initial_velocity.y * std::sin(angle);
            const double expected_y = initial_velocity.x * std::sin(angle) + initial_velocity.y * std::cos(angle);
            require_near(particle.velocity.x, expected_x, 1e-10,
                         "M5 2D uniform-B Simulation2D run did not apply Boris x-velocity rotation");
            require_near(particle.velocity.y, expected_y, 1e-10,
                         "M5 2D uniform-B Simulation2D run did not apply Boris y-velocity rotation");
            const double initial_speed = std::sqrt(initial_velocity.x * initial_velocity.x + initial_velocity.y * initial_velocity.y);
            const double final_speed = std::sqrt(particle.velocity.x * particle.velocity.x + particle.velocity.y * particle.velocity.y);
            require_near(final_speed, initial_speed, 1e-10,
                         "M5 2D uniform-B Simulation2D run did not conserve magnetic-rotation speed");
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
            // M1 benchmark: exact 2D CIC shape-function deposition for a single interior particle.
            pic::Mesh2D mesh(5, 6, 1.0, 1.5, pic::Boundary::Periodic);
            const pic::Vec2 position{0.31, 0.44};
            const double charge = -1.5;
            const double weight = 2.0;
            const double macro_charge = charge * weight;
            const double gx = position.x / mesh.dx();
            const double gy = position.y / mesh.dy();
            const std::size_t i = static_cast<std::size_t>(std::floor(gx));
            const std::size_t j = static_cast<std::size_t>(std::floor(gy));
            const double fx = gx - static_cast<double>(i);
            const double fy = gy - static_cast<double>(j);
            const std::vector<pic::Particle2D> particles{pic::Particle2D{position, pic::Vec2{}, true}};
            pic::deposit_charge_cic(mesh, particles, charge, weight);

            std::vector<double> expected(mesh.size(), 0.0);
            expected[mesh.index(i, j)] += macro_charge * (1.0 - fx) * (1.0 - fy) / (mesh.dx() * mesh.dy());
            expected[mesh.index(i + 1, j)] += macro_charge * fx * (1.0 - fy) / (mesh.dx() * mesh.dy());
            expected[mesh.index(i, j + 1)] += macro_charge * (1.0 - fx) * fy / (mesh.dx() * mesh.dy());
            expected[mesh.index(i + 1, j + 1)] += macro_charge * fx * fy / (mesh.dx() * mesh.dy());

            double deposited_charge = 0.0;
            for (std::size_t idx = 0; idx < mesh.size(); ++idx) {
                require_near(mesh.rho()[idx], expected[idx], 1e-12,
                             "M1 2D CIC exact shape-function benchmark density mismatch");
                deposited_charge += mesh.rho()[idx] * mesh.dx() * mesh.dy();
            }
            require_near(deposited_charge, macro_charge, 1e-12,
                         "M1 2D CIC exact shape-function benchmark did not conserve charge");
        }
        {
            // M1 benchmark: exact 3D CIC shape-function deposition for a single interior particle.
            pic::Mesh3D mesh(6, 5, 4, 1.2, 1.0, 0.8, pic::Boundary::Periodic);
            const pic::Vec3 position{0.37, 0.46, 0.51};
            const double charge = 1.25;
            const double weight = 0.8;
            const double macro_charge = charge * weight;
            const double gx = position.x / mesh.dx();
            const double gy = position.y / mesh.dy();
            const double gz = position.z / mesh.dz();
            const std::size_t i = static_cast<std::size_t>(std::floor(gx));
            const std::size_t j = static_cast<std::size_t>(std::floor(gy));
            const std::size_t k = static_cast<std::size_t>(std::floor(gz));
            const double fx = gx - static_cast<double>(i);
            const double fy = gy - static_cast<double>(j);
            const double fz = gz - static_cast<double>(k);
            const std::vector<pic::Particle3D> particles{pic::Particle3D{position, pic::Vec3{}, true}};
            pic::deposit_charge_cic(mesh, particles, charge, weight);

            std::vector<double> expected(mesh.size(), 0.0);
            for (std::size_t dz_i = 0; dz_i < 2; ++dz_i) {
                const double wz = dz_i == 0 ? 1.0 - fz : fz;
                for (std::size_t dy_i = 0; dy_i < 2; ++dy_i) {
                    const double wy = dy_i == 0 ? 1.0 - fy : fy;
                    for (std::size_t dx_i = 0; dx_i < 2; ++dx_i) {
                        const double wx = dx_i == 0 ? 1.0 - fx : fx;
                        expected[mesh.index(i + dx_i, j + dy_i, k + dz_i)] +=
                            macro_charge * wx * wy * wz / (mesh.dx() * mesh.dy() * mesh.dz());
                    }
                }
            }

            double deposited_charge = 0.0;
            for (std::size_t idx = 0; idx < mesh.size(); ++idx) {
                require_near(mesh.rho()[idx], expected[idx], 1e-12,
                             "M1 3D CIC exact shape-function benchmark density mismatch");
                deposited_charge += mesh.rho()[idx] * mesh.dx() * mesh.dy() * mesh.dz();
            }
            require_near(deposited_charge, macro_charge, 1e-12,
                         "M1 3D CIC exact shape-function benchmark did not conserve charge");
        }
        {
            // M1 benchmark: 2D bilinear electric interpolation is exact for affine fields.
            pic::Mesh2D mesh(6, 5, 1.2, 1.0, pic::Boundary::Periodic);
            for (std::size_t j = 0; j < mesh.ny(); ++j) {
                for (std::size_t i = 0; i < mesh.nx(); ++i) {
                    const double x = mesh.node_x(i);
                    const double y = mesh.node_y(j);
                    const auto idx = mesh.index(i, j);
                    mesh.electric_x()[idx] = -1.0 + 2.5 * x - 0.75 * y;
                    mesh.electric_y()[idx] = 0.5 - 1.25 * x + 1.5 * y;
                }
            }
            const pic::Vec2 position{0.53, 0.47};
            const pic::Vec2 electric = pic::interpolate_electric(mesh, position);
            require_near(electric.x, -1.0 + 2.5 * position.x - 0.75 * position.y, 1e-12,
                         "M1 2D affine interpolation benchmark Ex mismatch");
            require_near(electric.y, 0.5 - 1.25 * position.x + 1.5 * position.y, 1e-12,
                         "M1 2D affine interpolation benchmark Ey mismatch");

            const pic::Vec2 wrapped = pic::interpolate_electric(mesh, pic::Vec2{position.x + mesh.length_x(), position.y - mesh.length_y()});
            require_near(wrapped.x, electric.x, 1e-12,
                         "M1 2D affine interpolation benchmark did not wrap periodic x/y coordinates");
            require_near(wrapped.y, electric.y, 1e-12,
                         "M1 2D affine interpolation benchmark did not wrap periodic x/y coordinates");
        }
        {
            // M1 benchmark: 2D particle-boundary policies have analytic one-step outcomes.
            pic::Simulation2DConfig cfg;
            cfg.nx = 5;
            cfg.ny = 5;
            cfg.length_x = 1.0;
            cfg.length_y = 1.0;
            cfg.dt = 0.2;
            cfg.steps = 1;
            cfg.boundary = pic::Boundary::Dirichlet;
            cfg.output_interval = 1;
            cfg.particle_boundary_config.left = pic::ParticleBoundary::Absorbing;
            cfg.particle_boundary_config.right = pic::ParticleBoundary::Reflecting;
            cfg.particle_boundary_config.bottom = pic::ParticleBoundary::Periodic;
            cfg.particle_boundary_config.top = pic::ParticleBoundary::Reflecting;
            cfg.species = {
                pic::Species2DConfig{"left_absorb", 0.0, 1.0, 1.0, 1, -0.4, 0.0, 0.0, 0.03, 0.03, 0.50, 0.50},
                pic::Species2DConfig{"right_reflect", 0.0, 1.0, 1.0, 1, 0.5, 0.0, 0.0, 0.94, 0.94, 0.25, 0.25},
                pic::Species2DConfig{"bottom_periodic", 0.0, 1.0, 1.0, 1, 0.0, -0.4, 0.0, 0.20, 0.20, 0.03, 0.03},
                pic::Species2DConfig{"top_reflect", 0.0, 1.0, 1.0, 1, 0.0, 0.6, 0.0, 0.70, 0.70, 0.92, 0.92},
            };

            pic::Simulation2D sim(cfg);
            sim.run();
            const auto& losses = sim.boundary_losses();
            require(losses.absorbed_left == 1, "M1 2D particle-boundary benchmark left absorber count mismatch");
            require(losses.absorbed_right == 0 && losses.absorbed_bottom == 0 && losses.absorbed_top == 0,
                    "M1 2D particle-boundary benchmark reported unexpected non-left absorption");

            const auto& absorbed = sim.species()[0].particles()[0];
            const auto& right = sim.species()[1].particles()[0];
            const auto& bottom = sim.species()[2].particles()[0];
            const auto& top = sim.species()[3].particles()[0];
            require(!absorbed.alive, "M1 2D particle-boundary benchmark did not absorb left-crossing particle");
            require(right.alive && bottom.alive && top.alive,
                    "M1 2D particle-boundary benchmark lost a non-absorbing particle");
            require_near(right.position.x, 0.96, 1e-12,
                         "M1 2D particle-boundary benchmark right reflection position mismatch");
            require_near(right.velocity.x, -0.5, 1e-12,
                         "M1 2D particle-boundary benchmark right reflection velocity mismatch");
            require_near(bottom.position.y, 0.95, 1e-12,
                         "M1 2D particle-boundary benchmark bottom periodic wrap position mismatch");
            require_near(bottom.velocity.y, -0.4, 1e-12,
                         "M1 2D particle-boundary benchmark bottom periodic velocity mismatch");
            require_near(top.position.y, 0.96, 1e-12,
                         "M1 2D particle-boundary benchmark top reflection position mismatch");
            require_near(top.velocity.y, -0.6, 1e-12,
                         "M1 2D particle-boundary benchmark top reflection velocity mismatch");
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
            const auto output_dir = std::filesystem::path("test_output_diagnostics");
            std::filesystem::remove_all(output_dir);

            pic::Mesh2D mesh2d(3, 4, 1.0, 2.0, pic::Boundary::Dirichlet);
            for (std::size_t j = 0; j < mesh2d.ny(); ++j) {
                for (std::size_t i = 0; i < mesh2d.nx(); ++i) {
                    const auto idx = mesh2d.index(i, j);
                    mesh2d.rho()[idx] = static_cast<double>(idx);
                    mesh2d.phi()[idx] = 10.0 + static_cast<double>(idx);
                    mesh2d.electric_x()[idx] = 0.25 * static_cast<double>(i);
                    mesh2d.electric_y()[idx] = 0.5 * static_cast<double>(j);
                }
            }
            pic::write_legacy_vtk(mesh2d, output_dir / "manual_fields_2d.vtk", "AuroraPIC 2D test fields");
            const auto vtk2d = read_file_text(output_dir / "manual_fields_2d.vtk");
            require(vtk2d.find("# vtk DataFile Version 3.0") != std::string::npos, "2D VTK header is missing");
            require(vtk2d.find("DATASET STRUCTURED_GRID") != std::string::npos, "2D VTK dataset type is missing");
            require(vtk2d.find("DIMENSIONS 3 4 1") != std::string::npos, "2D VTK dimensions are wrong");
            require(vtk2d.find("POINTS 12 double") != std::string::npos, "2D VTK point count is wrong");
            require(vtk2d.find("POINT_DATA 12") != std::string::npos, "2D VTK point-data count is wrong");
            require(vtk2d.find("SCALARS rho double 1") != std::string::npos, "2D VTK rho scalar is missing");
            require(vtk2d.find("SCALARS phi double 1") != std::string::npos, "2D VTK phi scalar is missing");
            require(vtk2d.find("VECTORS electric double") != std::string::npos, "2D VTK electric vector is missing");
            require(vtk2d.find("0 0 0\n0.5 0 0\n1 0 0") != std::string::npos, "2D VTK point ordering changed");
            require(vtk2d.find("0 0 0\n0.25 0 0\n0.5 0 0") != std::string::npos, "2D VTK electric vector values are missing");

            pic::write_vtk_xml(mesh2d, output_dir / "manual_fields_2d.vts");
            const auto vts2d = read_file_text(output_dir / "manual_fields_2d.vts");
            require(vts2d.find("<VTKFile type=\"StructuredGrid\"") != std::string::npos, "2D VTK XML root is missing");
            require(vts2d.find("WholeExtent=\"0 2 0 3 0 0\"") != std::string::npos, "2D VTK XML extent is wrong");
            require(vts2d.find("<PointData Scalars=\"rho\" Vectors=\"electric\">") != std::string::npos, "2D VTK XML point data is missing");
            require(vts2d.find("Name=\"rho\" format=\"ascii\"") != std::string::npos, "2D VTK XML rho array is missing");
            require(vts2d.find("Name=\"phi\" format=\"ascii\"") != std::string::npos, "2D VTK XML phi array is missing");
            require(vts2d.find("Name=\"electric\" NumberOfComponents=\"3\" format=\"ascii\"") != std::string::npos, "2D VTK XML electric array is missing");
            require(vts2d.find("Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\"") != std::string::npos, "2D VTK XML points array is missing");
            require(vts2d.find("0 0 0 0.5 0 0 1 0 0") != std::string::npos, "2D VTK XML point ordering changed");

            pic::Mesh3D mesh3d(3, 3, 4, 2.0, 2.0, 3.0, pic::Boundary::Dirichlet);
            for (std::size_t k = 0; k < mesh3d.nz(); ++k) {
                for (std::size_t j = 0; j < mesh3d.ny(); ++j) {
                    for (std::size_t i = 0; i < mesh3d.nx(); ++i) {
                        const auto idx = mesh3d.index(i, j, k);
                        mesh3d.rho()[idx] = -static_cast<double>(idx);
                        mesh3d.phi()[idx] = 20.0 + static_cast<double>(idx);
                        mesh3d.electric_x()[idx] = 0.5 * static_cast<double>(i);
                        mesh3d.electric_y()[idx] = 0.25 * static_cast<double>(j);
                        mesh3d.electric_z()[idx] = -0.125 * static_cast<double>(k);
                    }
                }
            }
            pic::write_legacy_vtk(mesh3d, output_dir / "manual_fields_3d.vtk", "AuroraPIC 3D test fields");
            const auto vtk3d = read_file_text(output_dir / "manual_fields_3d.vtk");
            require(vtk3d.find("DIMENSIONS 3 3 4") != std::string::npos, "3D VTK dimensions are wrong");
            require(vtk3d.find("POINTS 36 double") != std::string::npos, "3D VTK point count is wrong");
            require(vtk3d.find("POINT_DATA 36") != std::string::npos, "3D VTK point-data count is wrong");
            require(vtk3d.find("SCALARS rho double 1") != std::string::npos, "3D VTK rho scalar is missing");
            require(vtk3d.find("SCALARS phi double 1") != std::string::npos, "3D VTK phi scalar is missing");
            require(vtk3d.find("VECTORS electric double") != std::string::npos, "3D VTK electric vector is missing");
            require(vtk3d.find("0 0 0\n1 0 0\n2 0 0\n0 1 0") != std::string::npos, "3D VTK point ordering changed");
            require(vtk3d.find("0 0 -0.125") != std::string::npos, "3D VTK electric z values are missing");

            pic::write_vtk_xml(mesh3d, output_dir / "manual_fields_3d.vts");
            const auto vts3d = read_file_text(output_dir / "manual_fields_3d.vts");
            require(vts3d.find("WholeExtent=\"0 2 0 2 0 3\"") != std::string::npos, "3D VTK XML extent is wrong");
            require(vts3d.find("Name=\"rho\" format=\"ascii\"") != std::string::npos, "3D VTK XML rho array is missing");
            require(vts3d.find("Name=\"phi\" format=\"ascii\"") != std::string::npos, "3D VTK XML phi array is missing");
            require(vts3d.find("Name=\"electric\" NumberOfComponents=\"3\" format=\"ascii\"") != std::string::npos, "3D VTK XML electric array is missing");
            require(vts3d.find("0 0 0 1 0 0 2 0 0 0 1 0") != std::string::npos, "3D VTK XML point ordering changed");
            require(vts3d.find("0 0 -0.125") != std::string::npos, "3D VTK XML electric z values are missing");

            pic::Species2DConfig cfg2d;
            cfg2d.name = "ions";
            cfg2d.particles = 2;
            pic::Species2D species2d(cfg2d);
            species2d.particles() = {
                pic::Particle2D{pic::Vec2{0.25, 0.5}, pic::Vec2{1.5, -0.25}, true, pic::Vec2{}},
                pic::Particle2D{pic::Vec2{0.75, 1.5}, pic::Vec2{-0.5, 0.25}, false, pic::Vec2{}},
            };
            pic::Diagnostics2D diagnostics2d(output_dir / "particles2d", {species2d});
            diagnostics2d.write_particle_sample(7, {species2d}, 1, 2);
            const auto particles2d = read_file_text(output_dir / "particles2d" / "particles_7.csv");
            require(particles2d.find("species_id,species,x,y,vx,vy,alive\n") == 0, "2D particle CSV header is wrong");
            require(particles2d.find("0,ions,0.25,0.5,1.5,-0.25,1\n") != std::string::npos, "2D particle CSV live row is wrong");
            require(particles2d.find("0,ions,0.75,1.5,-0.5,0.25,0\n") != std::string::npos, "2D particle CSV dead row is wrong");
            require(count_lines(particles2d) == 3, "2D particle CSV sample count is wrong");

            pic::Species3DConfig cfg3d;
            cfg3d.name = "electrons";
            cfg3d.particles = 2;
            pic::Species3D species3d(cfg3d);
            species3d.particles() = {
                pic::Particle3D{pic::Vec3{0.25, 0.5, 0.75}, pic::Vec3{1.0, -1.0, 0.5}, true, pic::Vec3{}},
                pic::Particle3D{pic::Vec3{0.75, 1.5, 2.25}, pic::Vec3{-0.25, 0.25, -0.5}, false, pic::Vec3{}},
            };
            pic::Diagnostics3D diagnostics3d(output_dir / "particles3d", {species3d});
            diagnostics3d.write_particle_sample(9, {species3d}, 2, 0);
            const auto particles3d = read_file_text(output_dir / "particles3d" / "particles_9.csv");
            require(particles3d.find("species_id,species,x,y,z,vx,vy,vz,alive\n") == 0, "3D particle CSV header is wrong");
            require(particles3d.find("0,electrons,0.25,0.5,0.75,1,-1,0.5,1\n") != std::string::npos, "3D particle CSV stride row is wrong");
            require(particles3d.find("2.25") == std::string::npos, "3D particle CSV stride was not honored");
            require(count_lines(particles3d) == 2, "3D particle CSV sample count is wrong");
            std::filesystem::remove_all(output_dir);
        }
        {
            require_throws([] { pic::Grid(2, 1.0, pic::Boundary::Periodic); }, "grid nx validation did not throw");
            require_throws([] { pic::Grid(32, 0.0, pic::Boundary::Periodic); }, "grid length validation did not throw");
            require_throws([] { pic::Grid(32, std::numeric_limits<double>::infinity(), pic::Boundary::Periodic); }, "grid non-finite length validation did not throw");
            require_throws([] { pic::Mesh2D(2, 4, 1.0, 1.0, pic::Boundary::Periodic); }, "2D mesh nx validation did not throw");
            require_throws([] { pic::Mesh2D(4, 4, 0.0, 1.0, pic::Boundary::Periodic); }, "2D mesh length validation did not throw");
            require_throws([] { pic::Mesh2D(4, 4, std::numeric_limits<double>::quiet_NaN(), 1.0, pic::Boundary::Periodic); }, "2D mesh non-finite length validation did not throw");
            require_throws([] {
                const auto overflowing_nx = std::numeric_limits<std::size_t>::max() / 3 + 1;
                pic::Mesh2D(overflowing_nx, 3, 1.0, 1.0, pic::Boundary::Periodic);
            }, "2D mesh node-count overflow validation did not throw");
            require_throws([] { pic::Mesh3D(4, 4, 4, 1.0, 1.0, std::numeric_limits<double>::infinity(), pic::Boundary::Periodic); }, "3D mesh non-finite length validation did not throw");
            require_throws([] {
                const auto overflowing_nx = std::numeric_limits<std::size_t>::max() / 9 + 1;
                pic::Mesh3D(overflowing_nx, 3, 3, 1.0, 1.0, 1.0, pic::Boundary::Periodic);
            }, "3D mesh node-count overflow validation did not throw");
            require_throws([] { pic::Simulation2DConfig cfg; cfg.output_interval = 0; pic::Simulation2D sim(cfg); }, "2D output_interval validation did not throw");
            require_throws([] { pic::Config cfg; cfg.output_interval = 0; pic::Simulation sim(cfg); }, "1D output_interval validation did not throw");
            require_throws([] { pic::Config cfg; cfg.dt = std::numeric_limits<double>::quiet_NaN(); pic::Simulation sim(cfg); }, "1D simulation accepted non-finite dt");
            require_throws([] { pic::Config cfg; cfg.phi_left = std::numeric_limits<double>::infinity(); pic::Simulation sim(cfg); }, "1D simulation accepted non-finite boundary potential");
            require_throws([] { pic::Config cfg; cfg.collisions.frequency = std::numeric_limits<double>::infinity(); pic::Simulation sim(cfg); }, "1D simulation accepted non-finite collision frequency");
            require_throws([] { pic::Config cfg; cfg.collisions.neutral_temperature_velocity = std::numeric_limits<double>::quiet_NaN(); pic::Simulation sim(cfg); }, "1D simulation accepted non-finite neutral collision velocity");
            require_throws([] { pic::Config cfg; cfg.mode = pic::RunMode::SteadyState; cfg.steady_tolerance = std::numeric_limits<double>::quiet_NaN(); pic::Simulation sim(cfg); }, "1D simulation accepted non-finite steady tolerance");
            require_throws([] { pic::Simulation2DConfig cfg; cfg.dt = std::numeric_limits<double>::quiet_NaN(); pic::Simulation2D sim(cfg); }, "2D simulation accepted non-finite dt");
            require_throws([] { pic::Simulation2DConfig cfg; cfg.magnetic_field_z = std::numeric_limits<double>::infinity(); pic::Simulation2D sim(cfg); }, "2D simulation accepted non-finite magnetic field");
            require_throws([] { pic::Simulation2DConfig cfg; cfg.boundary_config.left.potential = std::numeric_limits<double>::quiet_NaN(); pic::Simulation2D sim(cfg); }, "2D simulation accepted non-finite boundary potential");
            require_throws([] {
                pic::Simulation2DConfig cfg;
                cfg.mode = pic::RunMode::SteadyState;
                cfg.max_steps = 0;
                pic::Simulation2D sim(cfg);
            }, "2D steady-state simulation accepted zero max_steps");
            require_throws([] { pic::Simulation3DConfig cfg; cfg.dt = std::numeric_limits<double>::quiet_NaN(); pic::Simulation3D sim(cfg); }, "3D simulation accepted non-finite dt");
            require_throws([] { pic::Simulation3DConfig cfg; cfg.magnetic_field.y = std::numeric_limits<double>::infinity(); pic::Simulation3D sim(cfg); }, "3D simulation accepted non-finite magnetic field");
            require_throws([] {
                pic::Simulation3DConfig cfg;
                cfg.mode = pic::RunMode::SteadyState;
                cfg.steady_tolerance = std::numeric_limits<double>::quiet_NaN();
                pic::Simulation3D sim(cfg);
            }, "3D steady-state simulation accepted non-finite tolerance");
            require_throws([] {
                pic::SpeciesConfig cfg;
                cfg.charge = std::numeric_limits<double>::quiet_NaN();
                (void)pic::Species(cfg);
            }, "1D species constructor accepted non-finite charge");
            require_throws([] {
                pic::SpeciesConfig cfg;
                cfg.drift_velocity = std::numeric_limits<double>::infinity();
                (void)pic::Species(cfg);
            }, "1D species constructor accepted non-finite drift velocity");
            require_throws([] {
                pic::Species2DConfig cfg;
                cfg.charge = std::numeric_limits<double>::quiet_NaN();
                (void)pic::Species2D(cfg);
            }, "2D species constructor accepted non-finite charge");
            require_throws([] {
                pic::Species2DConfig cfg;
                cfg.drift_velocity_y = std::numeric_limits<double>::infinity();
                (void)pic::Species2D(cfg);
            }, "2D species constructor accepted non-finite drift velocity");
            require_throws([] {
                pic::Species3DConfig cfg;
                cfg.charge = std::numeric_limits<double>::quiet_NaN();
                (void)pic::Species3D(cfg);
            }, "3D species constructor accepted non-finite charge");
            require_throws([] {
                pic::Species3DConfig cfg;
                cfg.drift_velocity_z = std::numeric_limits<double>::infinity();
                (void)pic::Species3D(cfg);
            }, "3D species constructor accepted non-finite drift velocity");
        }
        {
            // M4 runtime scaling smoke tests: serial remains the deterministic baseline while OpenMP,
            // when compiled in, uses static scheduling and preserves single-rank results.
            pic::RuntimePolicy serial{};
            require(pic::to_string(serial.backend) == "serial", "M4 runtime backend string for serial changed");
            const auto serial_info = pic::runtime_info(serial);
            require(serial_info.backend == pic::RuntimeBackend::Serial, "M4 runtime_info did not preserve serial backend");
            require(serial_info.active_threads == 1, "M4 serial runtime should use one active thread");

            std::vector<std::size_t> squares(16, 0);
            pic::runtime_parallel_for(std::size_t{0}, squares.size(), serial, [&](std::size_t i) {
                squares[i] = i * i;
            });
            require(squares[7] == 49 && squares[15] == 225, "M4 serial runtime_parallel_for produced wrong values");

            require_throws([] { pic::RuntimePolicy p; p.threads = 0; pic::validate_runtime_policy(p); },
                           "M4 runtime_threads=0 validation did not throw");
            require_throws([] { pic::RuntimePolicy p; p.threads = 2; pic::validate_runtime_policy(p); },
                           "M4 serial runtime accepted multiple threads");
            require_throws([] { pic::RuntimePolicy p; p.backend = pic::RuntimeBackend::MPI; pic::validate_runtime_policy(p); },
                           "M4 MPI placeholder runtime did not throw");
            require_throws([] { pic::RuntimePolicy p; p.backend = pic::RuntimeBackend::GPU; pic::validate_runtime_policy(p); },
                           "M4 GPU placeholder runtime did not throw");

#ifdef AURORA_HAVE_OPENMP
            pic::RuntimePolicy openmp{pic::RuntimeBackend::OpenMP, 2};
            const auto openmp_info = pic::runtime_info(openmp);
            require(openmp_info.openmp_compiled, "M4 OpenMP runtime_info did not report compiled support");
            require(openmp_info.active_threads == 2, "M4 OpenMP runtime active thread count mismatch");
            std::vector<std::size_t> openmp_squares(16, 0);
            pic::runtime_parallel_for(std::size_t{0}, openmp_squares.size(), openmp, [&](std::size_t i) {
                openmp_squares[i] = i * i;
            });
            require(openmp_squares == squares, "M4 OpenMP runtime_parallel_for changed deterministic loop output");

            pic::Simulation2DConfig baseline;
            baseline.nx = 9;
            baseline.ny = 8;
            baseline.length_x = 1.0;
            baseline.length_y = 1.0;
            baseline.dt = 0.01;
            baseline.steps = 4;
            baseline.output_interval = 4;
            baseline.output_dir = "test_output_m4_runtime_serial";
            baseline.seed = 2468;
            baseline.species = {pic::Species2DConfig{"runtime_smoke", -1.0, 1.0, 0.01, 8,
                                                       0.03, -0.02, 0.0,
                                                       0.1, 0.9, 0.2, 0.8}};
            std::filesystem::remove_all(baseline.output_dir);
            pic::Simulation2D serial_sim(baseline);
            const auto serial_summary = serial_sim.run();

            auto parallel_cfg = baseline;
            parallel_cfg.runtime = openmp;
            parallel_cfg.output_dir = "test_output_m4_runtime_openmp";
            std::filesystem::remove_all(parallel_cfg.output_dir);
            pic::Simulation2D openmp_sim(parallel_cfg);
            const auto openmp_summary = openmp_sim.run();
            require_checkpoint_samples_close(serial_summary.final_sample, openmp_summary.final_sample,
                                             "M4 2D OpenMP deterministic scaling smoke");
            require_species_close(serial_sim.species(), openmp_sim.species(),
                                  "M4 2D OpenMP deterministic scaling smoke");
            std::filesystem::remove_all(baseline.output_dir);
            std::filesystem::remove_all(parallel_cfg.output_dir);
#endif
        }
        {
            const auto config_path = std::filesystem::path("test_density_config.ini");
            {
                std::ofstream out(config_path);
                out << "config_version = 1\n"
                    << "nx = 16\n"
                    << "length = 2.0\n"
                    << "dt = 0.01\n"
                    << "output_interval = 2\n"
                    << "runtime_backend = single\n"
                    << "runtime_threads = 1\n"
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
            require(cfg.runtime.backend == pic::RuntimeBackend::Serial && cfg.runtime.threads == 1,
                    "M4 1D runtime config aliases were not parsed");
            require(std::abs(cfg.species[0].weight - 0.5) < 1e-15, "density-derived macro-particle weight is wrong");
            std::filesystem::remove(config_path);

            auto require_config_rejects = [](const std::filesystem::path& path,
                                             const std::string& text,
                                             auto loader,
                                             const std::string& message) {
                require_throws([&] {
                    { std::ofstream out(path); out << text; }
                    try { (void)loader(path.string()); } catch (...) { std::filesystem::remove(path); throw; }
                    std::filesystem::remove(path);
                }, message);
            };

            require_config_rejects(
                "test_missing_scale_1d.ini",
                "nx = 16\nlength = 1\ndt = 0.01\n[species]\nname = missing_scale\ncharge = -1\nmass = 1\nparticles = 10\n",
                [](const std::string& path) { return pic::load_config(path); },
                "1D species without weight or density validation did not throw");
            require_config_rejects(
                "test_nonfinite_weight_1d.ini",
                "nx = 16\nlength = 1\ndt = 0.01\n[species]\nname = bad_weight\ncharge = -1\nmass = 1\nweight = inf\nparticles = 10\n",
                [](const std::string& path) { return pic::load_config(path); },
                "1D non-finite species weight validation did not throw");
            require_config_rejects(
                "test_invalid_runtime_threads.ini",
                "nx = 16\nlength = 1\ndt = 0.01\nruntime_threads = 0\n[species]\nname = bad_runtime\ncharge = -1\nmass = 1\nweight = 1\nparticles = 10\n",
                [](const std::string& path) { return pic::load_config(path); },
                "M4 invalid runtime_threads config validation did not throw");
            require_config_rejects(
                "test_invalid_runtime_backend.ini",
                "dimension = 2\nnx = 8\nny = 6\nlength_x = 2\nlength_y = 1\ndt = 0.01\nruntime_backend = distributed\n[species.electrons]\ncharge = -1\nmass = 1\nweight = 1\nparticles = 12\n",
                [](const std::string& path) { return pic::load_config_2d(path); },
                "M4 invalid runtime_backend config validation did not throw");
            require_config_rejects(
                "test_missing_scale_2d.ini",
                "dimension = 2\nnx = 8\nny = 6\nlength_x = 2\nlength_y = 1\ndt = 0.01\n[species.electrons]\ncharge = -1\nmass = 1\nparticles = 12\n",
                [](const std::string& path) { return pic::load_config_2d(path); },
                "2D species without weight or density validation did not throw");
            require_config_rejects(
                "test_nonfinite_density_2d.ini",
                "dimension = 2\nnx = 8\nny = 6\nlength_x = 2\nlength_y = 1\ndt = 0.01\n[species.electrons]\ncharge = -1\nmass = 1\ndensity = nan\nparticles = 12\n",
                [](const std::string& path) { return pic::load_config_2d(path); },
                "2D non-finite species density validation did not throw");
            require_config_rejects(
                "test_missing_scale_3d.ini",
                "dimension = 3\nnx = 6\nny = 5\nnz = 4\nlength_x = 2\nlength_y = 1\nlength_z = 1\ndt = 0.01\n[species.electrons]\ncharge = -1\nmass = 1\nparticles = 12\n",
                [](const std::string& path) { return pic::load_config_3d(path); },
                "3D species without weight or density validation did not throw");
            require_config_rejects(
                "test_nonfinite_weight_3d.ini",
                "dimension = 3\nnx = 6\nny = 5\nnz = 4\nlength_x = 2\nlength_y = 1\nlength_z = 1\ndt = 0.01\n[species.electrons]\ncharge = -1\nmass = 1\nweight = inf\nparticles = 12\n",
                [](const std::string& path) { return pic::load_config_3d(path); },
                "3D non-finite species weight validation did not throw");

            {
                const auto version_path = std::filesystem::path("test_unsupported_config_version.ini");
                { std::ofstream out(version_path); out << "config_version = 2\nnx = 16\nlength = 1\ndt = 0.01\n[species]\nname = bad_version\ncharge = -1\nmass = 1\nweight = 1\nparticles = 10\n"; }
                require_throws_contains([&] {
                    try { (void)pic::load_config(version_path.string()); } catch (...) { std::filesystem::remove(version_path); throw; }
                    std::filesystem::remove(version_path);
                }, "unsupported config_version 2", "M6 unsupported config_version diagnostic did not throw clearly");
            }

            const auto config_2d_path = std::filesystem::path("test_2d_config.ini");
            {
                std::ofstream out(config_2d_path);
                out << "config_version = 1\n"
                    << "dimension = 2\n"
                    << "nx = 8\n"
                    << "ny = 6\n"
                    << "length_x = 2.0\n"
                    << "length_y = 1.5\n"
                    << "dt = 0.01\n"
                    << "steps = 2\n"
                    << "mode = steady_state\n"
                    << "steady_tolerance = 0.001\n"
                    << "steady_window = 3\n"
                    << "max_steps = 7\n"
                    << "output_interval = 1\n"
                    << "output_dir = test_output_config_2d\n"
                    << "runtime_backend = serial\n"
                    << "runtime_threads = 1\n"
                    << "vtk_output = true\n"
                    << "vtk_format = both\n"
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
            require(cfg2.mode == pic::RunMode::SteadyState && cfg2.max_steps == 7,
                    "2D config did not load steady-state mode and step cap");
            require_near(cfg2.steady_tolerance, 0.001, 1e-15, "2D config did not load steady_tolerance");
            require(cfg2.steady_window == 3, "2D config did not load steady_window");
            require(cfg2.vtk_output, "2D config did not load vtk_output");
            require(cfg2.vtk_format == pic::VTKOutputFormat::Both, "2D config did not load vtk_format");
            require(cfg2.runtime.backend == pic::RuntimeBackend::Serial && cfg2.runtime.threads == 1,
                    "M4 2D runtime config was not parsed");
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
                const auto path = std::filesystem::path("test_invalid_1d_species_charge.ini");
                { std::ofstream out(path); out << "[species]\ncharge = nan\nweight = 1\n"; }
                try { (void)pic::load_config(path.string()); } catch (...) { std::filesystem::remove(path); throw; }
                std::filesystem::remove(path);
            }, "invalid 1D species charge validation did not throw");
            require_throws([] {
                const auto path = std::filesystem::path("test_invalid_1d_species_drift.ini");
                { std::ofstream out(path); out << "[species]\ndrift_velocity = inf\nweight = 1\n"; }
                try { (void)pic::load_config(path.string()); } catch (...) { std::filesystem::remove(path); throw; }
                std::filesystem::remove(path);
            }, "invalid 1D species drift velocity validation did not throw");
            require_throws([] {
                const auto path = std::filesystem::path("test_invalid_2d_species_charge.ini");
                { std::ofstream out(path); out << "dimension = 2\n[species]\ncharge = nan\nweight = 1\n"; }
                try { (void)pic::load_config_2d(path.string()); } catch (...) { std::filesystem::remove(path); throw; }
                std::filesystem::remove(path);
            }, "invalid 2D species charge validation did not throw");
            require_throws([] {
                const auto path = std::filesystem::path("test_invalid_2d_species_drift.ini");
                { std::ofstream out(path); out << "dimension = 2\n[species]\ndrift_velocity_x = inf\nweight = 1\n"; }
                try { (void)pic::load_config_2d(path.string()); } catch (...) { std::filesystem::remove(path); throw; }
                std::filesystem::remove(path);
            }, "invalid 2D species drift velocity validation did not throw");
            require_throws([] {
                const auto path = std::filesystem::path("test_invalid_3d_species_charge.ini");
                { std::ofstream out(path); out << "dimension = 3\n[species]\ncharge = nan\nweight = 1\n"; }
                try { (void)pic::load_config_3d(path.string()); } catch (...) { std::filesystem::remove(path); throw; }
                std::filesystem::remove(path);
            }, "invalid 3D species charge validation did not throw");
            require_throws([] {
                const auto path = std::filesystem::path("test_invalid_3d_species_drift.ini");
                { std::ofstream out(path); out << "dimension = 3\n[species]\ndrift_velocity_z = inf\nweight = 1\n"; }
                try { (void)pic::load_config_3d(path.string()); } catch (...) { std::filesystem::remove(path); throw; }
                std::filesystem::remove(path);
            }, "invalid 3D species drift velocity validation did not throw");
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
                out << "config_version = 1\n"
                    << "dimension = 3\n"
                    << "nx = 6\n"
                    << "ny = 5\n"
                    << "nz = 4\n"
                    << "length_x = 2.0\n"
                    << "length_y = 1.5\n"
                    << "length_z = 1.25\n"
                    << "dt = 0.01\n"
                    << "steps = 2\n"
                    << "mode = steady_state\n"
                    << "steady_tolerance = 0.002\n"
                    << "steady_window = 4\n"
                    << "max_steps = 9\n"
                    << "output_interval = 1\n"
                    << "output_dir = test_output_config_3d\n"
                    << "runtime_backend = serial\n"
                    << "runtime_threads = 1\n"
                    << "vtk_output = true\n"
                    << "vtk_format = vts\n"
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
            require(cfg3.mode == pic::RunMode::SteadyState && cfg3.max_steps == 9,
                    "3D config did not load steady-state mode and step cap");
            require_near(cfg3.steady_tolerance, 0.002, 1e-15, "3D config did not load steady_tolerance");
            require(cfg3.steady_window == 4, "3D config did not load steady_window");
            require(cfg3.vtk_output, "3D config did not load vtk_output");
            require(cfg3.vtk_format == pic::VTKOutputFormat::Xml, "3D config did not load vtk_format");
            require(cfg3.runtime.backend == pic::RuntimeBackend::Serial && cfg3.runtime.threads == 1,
                    "M4 3D runtime config was not parsed");
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
            cfg.vtk_format = pic::VTKOutputFormat::Both;
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
            require(std::filesystem::exists(cfg.output_dir / "fields_0.vtk"), "2D simulation did not write initial legacy VTK fields");
            require(std::filesystem::exists(cfg.output_dir / "fields_2.vtk"), "2D simulation did not write interval legacy VTK fields");
            require(std::filesystem::exists(cfg.output_dir / "fields_3.vtk"), "2D simulation did not write final legacy VTK fields");
            require(std::filesystem::exists(cfg.output_dir / "fields_0.vts"), "2D simulation did not write initial VTK XML fields");
            require(std::filesystem::exists(cfg.output_dir / "fields_2.vts"), "2D simulation did not write interval VTK XML fields");
            require(std::filesystem::exists(cfg.output_dir / "fields_3.vts"), "2D simulation did not write final VTK XML fields");
            require(!std::filesystem::exists(cfg.output_dir / "fields_1.vtk"), "2D simulation wrote an unexpected legacy VTK interval");
            require(!std::filesystem::exists(cfg.output_dir / "fields_1.vts"), "2D simulation wrote an unexpected VTK XML interval");
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
            const auto vts = read_file_text(cfg.output_dir / "fields_3.vts");
            require(vts.find("WholeExtent=\"0 15 0 11 0 0\"") != std::string::npos, "2D simulation VTK XML dimensions are wrong");
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
            cfg.vtk_format = pic::VTKOutputFormat::Xml;
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
            require(std::filesystem::exists(cfg.output_dir / "fields_0.vts"), "3D simulation did not write initial VTK XML fields");
            require(std::filesystem::exists(cfg.output_dir / "fields_2.vts"), "3D simulation did not write interval VTK XML fields");
            require(std::filesystem::exists(cfg.output_dir / "fields_3.vts"), "3D simulation did not write final VTK XML fields");
            require(!std::filesystem::exists(cfg.output_dir / "fields_0.vtk"), "3D XML-only simulation wrote legacy VTK fields");
            require(!std::filesystem::exists(cfg.output_dir / "fields_1.vts"), "3D simulation wrote an unexpected VTK XML interval");
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
            const auto vts3d = read_file_text(cfg.output_dir / "fields_3.vts");
            require(vts3d.find("WholeExtent=\"0 7 0 5 0 4\"") != std::string::npos, "3D simulation VTK XML dimensions are wrong");
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
            const auto checkpoint_sample = continuous.sample();
            const auto checkpoint_species = continuous.species();
            for (std::size_t n = continuous.step_count(); n < cfg.steps; ++n) continuous.step();

            pic::Simulation restarted(cfg);
            restarted.load_checkpoint(checkpoint_path);
            require(restarted.step_count() == 2, "1D checkpoint did not restore step count");
            require_near(restarted.time(), 2.0 * cfg.dt, 1e-15, "1D checkpoint did not restore time");
            require_checkpoint_samples_close(checkpoint_sample, restarted.sample(), "M1 1D checkpoint reload determinism");
            require_species_close(checkpoint_species, restarted.species(), "M1 1D checkpoint reload determinism");
            for (std::size_t n = restarted.step_count(); n < cfg.steps; ++n) restarted.step();

            require_checkpoint_samples_close(continuous.sample(), restarted.sample(), "M1 1D checkpoint restart determinism");
            require_species_close(continuous.species(), restarted.species(), "M1 1D checkpoint restart determinism");

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
            const auto checkpoint_sample = continuous.sample();
            const auto checkpoint_species = continuous.species();
            for (std::size_t n = continuous.step_count(); n < cfg.steps; ++n) continuous.step();

            pic::Simulation2D restarted(cfg);
            restarted.load_checkpoint(checkpoint_path);
            require(restarted.step_count() == 2, "2D checkpoint did not restore step count");
            require_near(restarted.time(), 2.0 * cfg.dt, 1e-15, "2D checkpoint did not restore time");
            require_checkpoint_samples_close(checkpoint_sample, restarted.sample(), "M1 2D checkpoint reload determinism");
            require_species_close(checkpoint_species, restarted.species(), "M1 2D checkpoint reload determinism");
            for (std::size_t n = restarted.step_count(); n < cfg.steps; ++n) restarted.step();

            require_checkpoint_samples_close(continuous.sample(), restarted.sample(), "M1 2D checkpoint restart determinism");
            require_species_close(continuous.species(), restarted.species(), "M1 2D checkpoint restart determinism");

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
            const auto checkpoint_sample = continuous.sample();
            const auto checkpoint_species = continuous.species();
            for (std::size_t n = continuous.step_count(); n < cfg.steps; ++n) continuous.step();

            pic::Simulation3D restarted(cfg);
            restarted.load_checkpoint(checkpoint_path);
            require(restarted.step_count() == 2, "3D checkpoint did not restore step count");
            require_near(restarted.time(), 2.0 * cfg.dt, 1e-15, "3D checkpoint did not restore time");
            require_checkpoint_samples_close(checkpoint_sample, restarted.sample(), "M1 3D checkpoint reload determinism");
            require_species_close(checkpoint_species, restarted.species(), "M1 3D checkpoint reload determinism");
            for (std::size_t n = restarted.step_count(); n < cfg.steps; ++n) restarted.step();

            require_checkpoint_samples_close(continuous.sample(), restarted.sample(), "M1 3D checkpoint restart determinism");
            require_species_close(continuous.species(), restarted.species(), "M1 3D checkpoint restart determinism");

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
            pic::Simulation2DConfig zero_cfg;
            zero_cfg.nx = 8;
            zero_cfg.ny = 8;
            zero_cfg.length_x = 1.0;
            zero_cfg.length_y = 1.0;
            zero_cfg.dt = 0.1;
            zero_cfg.steps = 1;
            zero_cfg.boundary = pic::Boundary::Periodic;
            zero_cfg.seed = 19;
            zero_cfg.output_interval = 1;
            zero_cfg.species = {pic::Species2DConfig{"probe2", 1.0, 1.0, 1e-30, 1, 1.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5}};

            auto boris_cfg = zero_cfg;
            boris_cfg.magnetic_field_z = 2.0;

            pic::Simulation2D zero_b(zero_cfg);
            pic::Simulation2D with_b(boris_cfg);
            zero_b.step();
            with_b.step();

            const auto& zero_particle = zero_b.species().front().particles().front();
            const auto& boris_particle = with_b.species().front().particles().front();
            const double angle = boris_rotation_angle(std::abs(boris_cfg.magnetic_field_z), 1.0, boris_cfg.dt);
            require_near(zero_particle.velocity.x, 1.0, 1e-12, "2D zero-B simulation changed x velocity unexpectedly");
            require_near(zero_particle.velocity.y, 0.0, 1e-12, "2D zero-B simulation changed y velocity unexpectedly");
            require_near(boris_particle.velocity.x, std::cos(angle), 1e-12,
                         "Simulation2D did not activate Boris x-velocity rotation");
            require_near(boris_particle.velocity.y, std::sin(angle), 1e-12,
                         "Simulation2D did not activate Boris y-velocity rotation");
            require(std::abs(boris_particle.velocity.y - zero_particle.velocity.y) > 1e-3,
                    "Simulation2D magnetic-field run did not diverge from zero-B velocity evolution");
        }
        {
            pic::Simulation3DConfig zero_cfg;
            zero_cfg.nx = 8;
            zero_cfg.ny = 8;
            zero_cfg.nz = 8;
            zero_cfg.length_x = 1.0;
            zero_cfg.length_y = 1.0;
            zero_cfg.length_z = 1.0;
            zero_cfg.dt = 0.05;
            zero_cfg.steps = 1;
            zero_cfg.boundary = pic::Boundary::Periodic;
            zero_cfg.seed = 23;
            zero_cfg.output_interval = 1;
            const pic::Vec3 initial_velocity{0.6, -0.25, 0.9};
            zero_cfg.species = {pic::Species3DConfig{"probe3", 1.0, 1.0, 1e-30, 1,
                                                      initial_velocity.x, initial_velocity.y, initial_velocity.z,
                                                      0.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5}};

            auto boris_cfg = zero_cfg;
            boris_cfg.magnetic_field = pic::Vec3{0.4, -0.8, 1.1};

            pic::Simulation3D zero_b(zero_cfg);
            pic::Simulation3D with_b(boris_cfg);
            zero_b.step();
            with_b.step();

            const auto& zero_particle = zero_b.species().front().particles().front();
            const auto& boris_particle = with_b.species().front().particles().front();
            const double magnetic_magnitude = norm(boris_cfg.magnetic_field);
            const pic::Vec3 axis = scale(boris_cfg.magnetic_field, 1.0 / magnetic_magnitude);
            const double angle = boris_rotation_angle(magnetic_magnitude, 1.0, boris_cfg.dt);
            const pic::Vec3 expected = rotate_about_axis(initial_velocity, axis, angle);
            require_near(zero_particle.velocity.x, initial_velocity.x, 1e-12, "3D zero-B simulation changed x velocity unexpectedly");
            require_near(zero_particle.velocity.y, initial_velocity.y, 1e-12, "3D zero-B simulation changed y velocity unexpectedly");
            require_near(zero_particle.velocity.z, initial_velocity.z, 1e-12, "3D zero-B simulation changed z velocity unexpectedly");
            require_near(boris_particle.velocity.x, expected.x, 1e-12,
                         "Simulation3D did not activate Boris x-velocity rotation");
            require_near(boris_particle.velocity.y, expected.y, 1e-12,
                         "Simulation3D did not activate Boris y-velocity rotation");
            require_near(boris_particle.velocity.z, expected.z, 1e-12,
                         "Simulation3D did not activate Boris z-velocity rotation");
            require(norm(pic::Vec3{boris_particle.velocity.x - zero_particle.velocity.x,
                                   boris_particle.velocity.y - zero_particle.velocity.y,
                                   boris_particle.velocity.z - zero_particle.velocity.z}) > 1e-3,
                    "Simulation3D magnetic-field run did not diverge from zero-B velocity evolution");
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
        {
            // M1 benchmark: representative 2D periodic neutral-tracer example has analytic drift invariants.
            pic::Simulation2DConfig cfg;
            cfg.nx = 9;
            cfg.ny = 7;
            cfg.length_x = 1.0;
            cfg.length_y = 0.75;
            cfg.dt = 0.25;
            cfg.steps = 3;
            cfg.boundary = pic::Boundary::Periodic;
            cfg.output_interval = 3;
            cfg.output_dir = "test_output_m1_representative_2d";
            cfg.seed = 101;
            cfg.species = {pic::Species2DConfig{"neutral_drift_2d", 0.0, 2.0, 0.5, 1,
                                                  0.35, -0.4, 0.0,
                                                  0.92, 0.92, 0.15, 0.15}};
            std::filesystem::remove_all(cfg.output_dir);

            pic::Simulation2D sim(cfg);
            const auto summary = sim.run();
            const auto& particle = sim.species().front().particles().front();
            require(summary.steps_completed == cfg.steps, "M1 representative 2D drift example did not complete requested steps");
            require_near(summary.final_time, cfg.dt * static_cast<double>(cfg.steps), 1e-15,
                         "M1 representative 2D drift example final time mismatch");
            require(summary.final_sample.live_particles == 1,
                    "M1 representative 2D drift example lost neutral tracer");
            require_near(summary.final_sample.charge_l1, 0.0, 1e-15,
                         "M1 representative 2D drift example accumulated charge density");
            require_near(summary.final_sample.field_energy, 0.0, 1e-15,
                         "M1 representative 2D drift example accumulated field energy");
            require_near(summary.final_sample.kinetic_energy, 0.14125, 1e-14,
                         "M1 representative 2D drift example kinetic-energy invariant mismatch");
            require(particle.alive, "M1 representative 2D drift example particle is not alive");
            require_near(particle.position.x, 0.1825, 1e-12,
                         "M1 representative 2D drift example periodic x position mismatch");
            require_near(particle.position.y, 0.6, 1e-12,
                         "M1 representative 2D drift example periodic y position mismatch");
            require_near(particle.velocity.x, 0.35, 1e-12,
                         "M1 representative 2D drift example vx invariant mismatch");
            require_near(particle.velocity.y, -0.4, 1e-12,
                         "M1 representative 2D drift example vy invariant mismatch");
            std::filesystem::remove_all(cfg.output_dir);
        }
        {
            // M1 benchmark: representative 3D periodic neutral-tracer example has analytic drift invariants.
            pic::Simulation3DConfig cfg;
            cfg.nx = 7;
            cfg.ny = 6;
            cfg.nz = 5;
            cfg.length_x = 1.2;
            cfg.length_y = 1.0;
            cfg.length_z = 0.8;
            cfg.dt = 0.2;
            cfg.steps = 4;
            cfg.boundary = pic::Boundary::Periodic;
            cfg.output_interval = 4;
            cfg.output_dir = "test_output_m1_representative_3d";
            cfg.seed = 202;
            cfg.species = {pic::Species3DConfig{"neutral_drift_3d", 0.0, 1.5, 0.4, 1,
                                                  0.25, -0.35, 0.30, 0.0,
                                                  1.05, 1.05, 0.12, 0.12, 0.72, 0.72}};
            std::filesystem::remove_all(cfg.output_dir);

            pic::Simulation3D sim(cfg);
            const auto summary = sim.run();
            const auto& particle = sim.species().front().particles().front();
            require(summary.steps_completed == cfg.steps, "M1 representative 3D drift example did not complete requested steps");
            require_near(summary.final_time, cfg.dt * static_cast<double>(cfg.steps), 1e-15,
                         "M1 representative 3D drift example final time mismatch");
            require(summary.final_sample.live_particles == 1,
                    "M1 representative 3D drift example lost neutral tracer");
            require_near(summary.final_sample.charge_l1, 0.0, 1e-15,
                         "M1 representative 3D drift example accumulated charge density");
            require_near(summary.final_sample.field_energy, 0.0, 1e-15,
                         "M1 representative 3D drift example accumulated field energy");
            require_near(summary.final_sample.kinetic_energy, 0.0825, 1e-14,
                         "M1 representative 3D drift example kinetic-energy invariant mismatch");
            require(particle.alive, "M1 representative 3D drift example particle is not alive");
            require_near(particle.position.x, 0.05, 1e-12,
                         "M1 representative 3D drift example periodic x position mismatch");
            require_near(particle.position.y, 0.84, 1e-12,
                         "M1 representative 3D drift example periodic y position mismatch");
            require_near(particle.position.z, 0.16, 1e-12,
                         "M1 representative 3D drift example periodic z position mismatch");
            require_near(particle.velocity.x, 0.25, 1e-12,
                         "M1 representative 3D drift example vx invariant mismatch");
            require_near(particle.velocity.y, -0.35, 1e-12,
                         "M1 representative 3D drift example vy invariant mismatch");
            require_near(particle.velocity.z, 0.30, 1e-12,
                         "M1 representative 3D drift example vz invariant mismatch");
            std::filesystem::remove_all(cfg.output_dir);
        }
        {
            const auto output_dir = std::filesystem::path("test_output_steady_1d");
            std::filesystem::remove_all(output_dir);

            pic::Config cfg;
            cfg.nx = 12;
            cfg.dt = 0.01;
            cfg.mode = pic::RunMode::SteadyState;
            cfg.steady_window = 2;
            cfg.steady_tolerance = 1e-12;
            cfg.max_steps = 10;
            cfg.output_interval = 1;
            cfg.output_dir = output_dir.string();
            cfg.checkpoint_output = true;
            cfg.checkpoint_interval = 10;
            cfg.species = {pic::SpeciesConfig{"steady_neutral_1d", 0.0, 1.0, 1.0, 4,
                                               1.0, 0.0, 0.0, 0.0, -1.0}};

            pic::Simulation simulation(cfg);
            const auto summary = simulation.run();
            require(summary.steady_state_reached, "1D steady-state run did not report convergence");
            require(summary.steps_completed == 3, "1D steady-state run did not stop at the first complete energy windows");
            require(std::filesystem::exists(output_dir / "checkpoint_3.apc"),
                    "1D steady-state convergence did not force a final checkpoint");
            std::filesystem::remove_all(output_dir);
        }
        {
            const auto output_dir = std::filesystem::path("test_output_steady_2d");
            std::filesystem::remove_all(output_dir);

            pic::Simulation2DConfig cfg;
            cfg.nx = 6;
            cfg.ny = 5;
            cfg.dt = 0.01;
            cfg.mode = pic::RunMode::SteadyState;
            cfg.steady_window = 2;
            cfg.steady_tolerance = 1e-12;
            cfg.max_steps = 10;
            cfg.output_interval = 1;
            cfg.output_dir = output_dir;
            cfg.checkpoint_output = true;
            cfg.checkpoint_interval = 10;
            cfg.species = {pic::Species2DConfig{"steady_neutral_2d", 0.0, 1.0, 1.0, 4,
                                                 0.0, 0.0, 0.0,
                                                 0.0, -1.0, 0.0, -1.0}};

            pic::Simulation2D simulation(cfg);
            const auto summary = simulation.run();
            require(summary.steady_state_reached, "2D steady-state run did not report convergence");
            require(summary.steps_completed == 3, "2D steady-state run did not stop at the first complete energy windows");
            require_near(summary.final_time, 3.0 * cfg.dt, 1e-15, "2D steady-state final time mismatch");
            require(std::filesystem::exists(output_dir / "checkpoint_3.apc"),
                    "2D steady-state convergence did not force a final checkpoint");
            std::filesystem::remove_all(output_dir);

            auto capped_cfg = cfg;
            capped_cfg.output_dir = "test_output_steady_2d_capped";
            capped_cfg.max_steps = 2;
            capped_cfg.checkpoint_output = false;
            std::filesystem::remove_all(capped_cfg.output_dir);
            pic::Simulation2D capped_simulation(capped_cfg);
            const auto capped_summary = capped_simulation.run();
            require(!capped_summary.steady_state_reached,
                    "2D capped steady-state run incorrectly reported convergence");
            require(capped_summary.steps_completed == capped_cfg.max_steps,
                    "2D capped steady-state run did not stop at max_steps");
            std::filesystem::remove_all(capped_cfg.output_dir);
        }
        {
            const auto output_dir = std::filesystem::path("test_output_steady_3d");
            std::filesystem::remove_all(output_dir);

            pic::Simulation3DConfig cfg;
            cfg.nx = 5;
            cfg.ny = 4;
            cfg.nz = 3;
            cfg.dt = 0.01;
            cfg.mode = pic::RunMode::SteadyState;
            cfg.steady_window = 2;
            cfg.steady_tolerance = 1e-12;
            cfg.max_steps = 10;
            cfg.output_interval = 1;
            cfg.output_dir = output_dir;
            cfg.checkpoint_output = true;
            cfg.checkpoint_interval = 10;
            cfg.species = {pic::Species3DConfig{"steady_neutral_3d", 0.0, 1.0, 1.0, 4,
                                                 0.0, 0.0, 0.0, 0.0,
                                                 0.0, -1.0, 0.0, -1.0, 0.0, -1.0}};

            pic::Simulation3D simulation(cfg);
            const auto summary = simulation.run();
            require(summary.steady_state_reached, "3D steady-state run did not report convergence");
            require(summary.steps_completed == 3, "3D steady-state run did not stop at the first complete energy windows");
            require_near(summary.final_time, 3.0 * cfg.dt, 1e-15, "3D steady-state final time mismatch");
            require(std::filesystem::exists(output_dir / "checkpoint_3.apc"),
                    "3D steady-state convergence did not force a final checkpoint");
            std::filesystem::remove_all(output_dir);
        }
        {
            // M2 benchmark: tagged 2D Gmsh v2 ASCII meshes import into AuroraPIC-owned labels.
            const auto fixture = std::filesystem::path(AURORA_TEST_SOURCE_DIR) / "tests" / "fixtures" / "tagged_square_v2.msh";
            const pic::ImportedMesh2D mesh = pic::load_gmsh2_ascii_mesh2d(fixture);
            require(mesh.nodes().size() == 6, "M2 Gmsh importer node count mismatch");
            require(mesh.cells().size() == 3, "M2 Gmsh importer cell count mismatch");
            require(mesh.boundary_faces().size() == 6, "M2 Gmsh importer boundary face count mismatch");
            require(mesh.physical_names().size() == 5, "M2 Gmsh importer physical-name count mismatch");

            const auto labels = mesh.boundary_labels();
            require(labels == std::vector<std::string>({"electrode", "inlet", "outlet", "wall"}),
                    "M2 Gmsh importer boundary labels were not preserved as internal tags");
            require_near(mesh.min_corner().x, 0.0, 1e-15, "M2 Gmsh importer min x mismatch");
            require_near(mesh.min_corner().y, 0.0, 1e-15, "M2 Gmsh importer min y mismatch");
            require_near(mesh.max_corner().x, 1.0, 1e-15, "M2 Gmsh importer max x mismatch");
            require_near(mesh.max_corner().y, 1.0, 1e-15, "M2 Gmsh importer max y mismatch");
            require_near(mesh.node_by_id(6).position.x, 0.0, 1e-15, "M2 Gmsh importer node lookup x mismatch");
            require_near(mesh.node_by_id(6).position.y, 1.0, 1e-15, "M2 Gmsh importer node lookup y mismatch");
            require(mesh.label_for_physical_tag(1, 3) == "wall", "M2 Gmsh importer boundary tag lookup mismatch");
            require(mesh.label_for_physical_tag(2, 10) == "plasma", "M2 Gmsh importer region tag lookup mismatch");
            require(mesh.label_for_physical_tag(1, 99) == "boundary_physical_99",
                    "M2 Gmsh importer boundary fallback label mismatch");
            require(mesh.label_for_physical_tag(2, 42) == "region_physical_42",
                    "M2 Gmsh importer region fallback label mismatch");

            const auto& first_cell = mesh.cells().front();
            require(first_cell.shape == pic::ImportedCellShape2D::Triangle,
                    "M2 Gmsh importer first cell shape mismatch");
            require(first_cell.label == "plasma", "M2 Gmsh importer cell label mismatch");
            require(first_cell.node_ids == std::vector<std::size_t>({1, 2, 6}),
                    "M2 Gmsh importer triangle connectivity mismatch");
            const auto& quad_cell = mesh.cells().back();
            require(quad_cell.shape == pic::ImportedCellShape2D::Quadrilateral,
                    "M2 Gmsh importer quadrilateral cell shape mismatch");
            require(quad_cell.node_ids == std::vector<std::size_t>({2, 3, 4, 5}),
                    "M2 Gmsh importer quadrilateral connectivity mismatch");
            require(mesh.boundary_faces().front().label == "inlet", "M2 Gmsh importer boundary label mismatch");
            require(mesh.boundary_faces().front().node_ids == std::array<std::size_t, 2>{1, 2},
                    "M2 Gmsh importer boundary connectivity mismatch");
            require(mesh.cell_by_id(9).shape == pic::ImportedCellShape2D::Quadrilateral,
                    "M2 Gmsh importer cell lookup mismatch");
            require(mesh.boundary_face_by_id(3).label == "outlet",
                    "M2 Gmsh importer boundary-face lookup mismatch");
            require_near(mesh.cell_area(7), 0.25, 1e-15, "M2 Gmsh importer first triangle area mismatch");
            require_near(mesh.cell_area(8), 0.25, 1e-15, "M2 Gmsh importer second triangle area mismatch");
            require_near(mesh.cell_area(9), 0.5, 1e-15, "M2 Gmsh importer quadrilateral area mismatch");
            require_near(mesh.total_area(), 1.0, 1e-15, "M2 Gmsh importer total area mismatch");
            require_near(mesh.cell_centroid(7).x, 1.0 / 6.0, 1e-15,
                         "M2 Gmsh importer triangle centroid x mismatch");
            require_near(mesh.cell_centroid(7).y, 1.0 / 3.0, 1e-15,
                         "M2 Gmsh importer triangle centroid y mismatch");
            require_near(mesh.cell_centroid(9).x, 0.75, 1e-15,
                         "M2 Gmsh importer quadrilateral centroid x mismatch");
            require_near(mesh.cell_centroid(9).y, 0.5, 1e-15,
                         "M2 Gmsh importer quadrilateral centroid y mismatch");
            require_near(mesh.boundary_face_length(1), 0.5, 1e-15,
                         "M2 Gmsh importer boundary-face length mismatch");

            const auto triangle_location = mesh.locate_point({1.0 / 6.0, 1.0 / 3.0});
            require(triangle_location.has_value() && triangle_location->cell_id == 7,
                    "imported triangle point location selected the wrong cell");
            require(triangle_location->node_ids == std::vector<std::size_t>({1, 2, 6}),
                    "imported triangle point location returned wrong connectivity");
            for (const double weight : triangle_location->shape_weights) {
                require_near(weight, 1.0 / 3.0, 1e-14,
                             "imported triangle barycentric coordinate mismatch");
            }

            const auto quad_location = mesh.locate_point({0.8, 0.25});
            require(quad_location.has_value() && quad_location->cell_id == 9,
                    "imported quadrilateral point location selected the wrong cell");
            const std::vector<double> expected_quad_weights{0.3, 0.45, 0.15, 0.1};
            for (std::size_t i = 0; i < expected_quad_weights.size(); ++i) {
                require_near(quad_location->shape_weights[i], expected_quad_weights[i], 1e-14,
                             "imported quadrilateral shape coordinate mismatch");
            }
            require(!mesh.locate_point({1.1, 0.5}).has_value(),
                    "imported mesh point location accepted an exterior point");
            require_throws([&]() {
                (void)mesh.locate_point({std::numeric_limits<double>::quiet_NaN(), 0.0});
            }, "imported mesh point location accepted a non-finite query");

            pic::UnstructuredMesh2D computational_mesh(mesh);
            require(computational_mesh.size() == mesh.nodes().size(),
                    "unstructured computational mesh node count mismatch");
            require_near(computational_mesh.node_control_area(1), 1.0 / 12.0, 1e-14,
                         "unstructured node 1 control area mismatch");
            require_near(computational_mesh.node_control_area(2), 7.0 / 24.0, 1e-14,
                         "unstructured node 2 control area mismatch");
            double control_area_sum = 0.0;
            for (const double area : computational_mesh.node_control_areas()) control_area_sum += area;
            require_near(control_area_sum, mesh.total_area(), 1e-14,
                         "unstructured lumped control areas do not recover domain area");

            std::vector<pic::Particle2D> particles(4);
            particles[0].position = {1.0 / 6.0, 1.0 / 3.0};
            particles[1].position = {0.8, 0.25};
            particles[2].position = {1.1, 0.5};
            particles[3].position = {0.25, 0.25};
            particles[3].alive = false;
            const auto deposit = pic::deposit_charge_shape(computational_mesh, particles, -2.0, 0.5);
            require(deposit.deposited_particles == 2,
                    "unstructured charge deposition counted the wrong in-domain particles");
            require(deposit.outside_particles == 1,
                    "unstructured charge deposition did not report an exterior live particle");
            require_near(deposit.deposited_charge, -2.0, 1e-15,
                         "unstructured charge deposition summary mismatch");
            double integrated_charge = 0.0;
            for (std::size_t i = 0; i < computational_mesh.size(); ++i) {
                integrated_charge +=
                    computational_mesh.rho()[i] * computational_mesh.node_control_areas()[i];
            }
            require_near(integrated_charge, deposit.deposited_charge, 1e-14,
                         "unstructured shape deposition did not conserve charge");
            computational_mesh.clear_charge();
            require(std::all_of(computational_mesh.rho().begin(), computational_mesh.rho().end(),
                                [](double value) { return value == 0.0; }),
                    "unstructured clear_charge did not reset all nodal densities");

            for (const auto& node : mesh.nodes()) {
                computational_mesh.electric()[computational_mesh.node_index(node.id)] = {
                    2.0 * node.position.x - 3.0 * node.position.y + 1.0,
                    -node.position.x + 4.0 * node.position.y - 2.0,
                };
            }
            for (const pic::Vec2 point : {pic::Vec2{0.2, 0.2}, pic::Vec2{0.8, 0.25}}) {
                const auto electric = pic::interpolate_electric(computational_mesh, point);
                require(electric.has_value(), "unstructured field interpolation rejected an interior point");
                require_near(electric->x, 2.0 * point.x - 3.0 * point.y + 1.0, 1e-14,
                             "unstructured affine electric interpolation x mismatch");
                require_near(electric->y, -point.x + 4.0 * point.y - 2.0, 1e-14,
                             "unstructured affine electric interpolation y mismatch");
            }
            require(!pic::interpolate_electric(computational_mesh, {1.1, 0.5}).has_value(),
                    "unstructured field interpolation accepted an exterior point");

            {
                pic::ImportedMesh2D distorted;
                distorted.add_node({1, {0.0, 0.0}});
                distorted.add_node({2, {2.0, 0.0}});
                distorted.add_node({3, {1.5, 1.0}});
                distorted.add_node({4, {0.0, 1.0}});
                distorted.add_cell({5, pic::ImportedCellShape2D::Quadrilateral,
                                    {1, 2, 3, 4}, 1, "plasma"});
                distorted.add_boundary_face({6, {1, 2}, 1, "wall"});
                distorted.add_boundary_face({7, {2, 3}, 1, "wall"});
                distorted.add_boundary_face({8, {3, 4}, 1, "wall"});
                distorted.add_boundary_face({9, {4, 1}, 1, "wall"});
                distorted.validate();

                const double xi = 0.2;
                const double eta = -0.4;
                const std::array<double, 4> weights{
                    0.25 * (1.0 - xi) * (1.0 - eta),
                    0.25 * (1.0 + xi) * (1.0 - eta),
                    0.25 * (1.0 + xi) * (1.0 + eta),
                    0.25 * (1.0 - xi) * (1.0 + eta),
                };
                pic::Vec2 mapped{};
                for (std::size_t i = 0; i < weights.size(); ++i) {
                    const auto position = distorted.node_by_id(i + 1).position;
                    mapped.x += weights[i] * position.x;
                    mapped.y += weights[i] * position.y;
                }
                const auto recovered = distorted.cell_coordinates(5, mapped);
                require(recovered.has_value(),
                        "distorted quadrilateral inverse isoparametric mapping did not converge");
                for (std::size_t i = 0; i < weights.size(); ++i) {
                    require_near(recovered->shape_weights[i], weights[i], 1e-13,
                                 "distorted quadrilateral shape coordinate mismatch");
                }
            }

            {
                // Executable imported-geometry benchmark: a symmetric four-triangle
                // domain has one interior degree of freedom with an analytic solution.
                pic::ImportedMesh2D star;
                star.add_node({1, {0.0, 0.0}});
                star.add_node({2, {1.0, 0.0}});
                star.add_node({3, {1.0, 1.0}});
                star.add_node({4, {0.0, 1.0}});
                star.add_node({5, {0.5, 0.5}});
                star.add_cell({10, pic::ImportedCellShape2D::Triangle, {1, 2, 5}, 1, "plasma"});
                star.add_cell({11, pic::ImportedCellShape2D::Triangle, {2, 3, 5}, 1, "plasma"});
                star.add_cell({12, pic::ImportedCellShape2D::Triangle, {3, 4, 5}, 1, "plasma"});
                star.add_cell({13, pic::ImportedCellShape2D::Triangle, {4, 1, 5}, 1, "plasma"});
                star.add_boundary_face({20, {1, 2}, 1, "ground"});
                star.add_boundary_face({21, {2, 3}, 1, "ground"});
                star.add_boundary_face({22, {3, 4}, 1, "ground"});
                star.add_boundary_face({23, {4, 1}, 1, "ground"});

                pic::UnstructuredMesh2D field_mesh(star);
                const auto constant_summary =
                    pic::solve_unstructured_poisson(field_mesh, {{"ground", 2.0}});
                require(constant_summary.converged,
                        "unstructured Poisson solve did not converge for constant Dirichlet data");
                for (const double potential : field_mesh.phi()) {
                    require_near(potential, 2.0, 1e-12,
                                 "unstructured Poisson solve did not preserve constant potential");
                }
                for (const auto electric : field_mesh.electric()) {
                    require_near(electric.x, 0.0, 1e-12,
                                 "constant unstructured potential produced a nonzero Ex");
                    require_near(electric.y, 0.0, 1e-12,
                                 "constant unstructured potential produced a nonzero Ey");
                }

                field_mesh.clear_charge();
                std::fill(field_mesh.phi().begin(), field_mesh.phi().end(), 0.0);
                const std::size_t center = field_mesh.node_index(5);
                field_mesh.rho()[center] = 1.0 / field_mesh.node_control_area(5);
                const auto source_summary =
                    pic::solve_unstructured_poisson(field_mesh, {{"ground", 0.0}});
                require(source_summary.converged,
                        "unstructured Poisson solve did not converge for the symmetric source");
                require(source_summary.iterations == 1,
                        "one-degree-of-freedom unstructured Poisson solve did not converge in one iteration");
                require_near(field_mesh.phi()[center], 0.25, 1e-12,
                             "unstructured Poisson analytic center potential mismatch");
                require_near(field_mesh.electric()[center].x, 0.0, 1e-12,
                             "symmetric unstructured source produced nonzero center Ex");
                require_near(field_mesh.electric()[center].y, 0.0, 1e-12,
                             "symmetric unstructured source produced nonzero center Ey");

                require_throws(
                    [&]() {
                        (void)pic::solve_unstructured_poisson(field_mesh, {});
                    },
                    "unstructured Poisson solve accepted missing boundary potentials");
                require_throws(
                    [&]() {
                        (void)pic::solve_unstructured_poisson(
                            field_mesh, {{"ground", 0.0}, {"unknown", 0.0}});
                    },
                    "unstructured Poisson solve accepted an unknown boundary label");
                require_throws(
                    [&]() {
                        (void)pic::solve_unstructured_poisson(
                            field_mesh, {{"ground", std::numeric_limits<double>::infinity()}});
                    },
                    "unstructured Poisson solve accepted a non-finite boundary potential");
            }

            {
                // Four Q1 elements provide an independent analytic check of the
                // quadrilateral stiffness/Jacobian assembly.
                pic::ImportedMesh2D quad_patch;
                quad_patch.add_node({1, {0.0, 0.0}});
                quad_patch.add_node({2, {0.5, 0.0}});
                quad_patch.add_node({3, {1.0, 0.0}});
                quad_patch.add_node({4, {0.0, 0.5}});
                quad_patch.add_node({5, {0.5, 0.5}});
                quad_patch.add_node({6, {1.0, 0.5}});
                quad_patch.add_node({7, {0.0, 1.0}});
                quad_patch.add_node({8, {0.5, 1.0}});
                quad_patch.add_node({9, {1.0, 1.0}});
                quad_patch.add_cell({10, pic::ImportedCellShape2D::Quadrilateral,
                                     {1, 2, 5, 4}, 1, "plasma"});
                quad_patch.add_cell({11, pic::ImportedCellShape2D::Quadrilateral,
                                     {2, 3, 6, 5}, 1, "plasma"});
                quad_patch.add_cell({12, pic::ImportedCellShape2D::Quadrilateral,
                                     {4, 5, 8, 7}, 1, "plasma"});
                quad_patch.add_cell({13, pic::ImportedCellShape2D::Quadrilateral,
                                     {5, 6, 9, 8}, 1, "plasma"});
                quad_patch.add_boundary_face({20, {1, 2}, 1, "ground"});
                quad_patch.add_boundary_face({21, {2, 3}, 1, "ground"});
                quad_patch.add_boundary_face({22, {3, 6}, 1, "ground"});
                quad_patch.add_boundary_face({23, {6, 9}, 1, "ground"});
                quad_patch.add_boundary_face({24, {9, 8}, 1, "ground"});
                quad_patch.add_boundary_face({25, {8, 7}, 1, "ground"});
                quad_patch.add_boundary_face({26, {7, 4}, 1, "ground"});
                quad_patch.add_boundary_face({27, {4, 1}, 1, "ground"});

                pic::UnstructuredMesh2D quad_field_mesh(quad_patch);
                const std::size_t center = quad_field_mesh.node_index(5);
                quad_field_mesh.rho()[center] =
                    1.0 / quad_field_mesh.node_control_area(5);
                const auto summary =
                    pic::solve_unstructured_poisson(quad_field_mesh, {{"ground", 0.0}});
                require(summary.converged && summary.iterations == 1,
                        "one-degree-of-freedom Q1 Poisson solve did not converge in one iteration");
                require_near(quad_field_mesh.phi()[center], 3.0 / 8.0, 1e-12,
                             "unstructured Q1 analytic center potential mismatch");
                require_near(quad_field_mesh.electric()[center].x, 0.0, 1e-12,
                             "symmetric Q1 source produced nonzero center Ex");
                require_near(quad_field_mesh.electric()[center].y, 0.0, 1e-12,
                             "symmetric Q1 source produced nonzero center Ey");
            }

            {
                pic::UnstructuredMesh2D mixed_field_mesh(mesh);
                const auto summary = pic::solve_unstructured_poisson(
                    mixed_field_mesh,
                    {{"electrode", 1.25}, {"inlet", 1.25}, {"outlet", 1.25}, {"wall", 1.25}});
                require(summary.converged,
                        "mixed triangle/quadrilateral unstructured Poisson solve did not converge");
                for (const double potential : mixed_field_mesh.phi()) {
                    require_near(potential, 1.25, 1e-12,
                                 "mixed-cell unstructured solve did not preserve constant potential");
                }
                for (const auto electric : mixed_field_mesh.electric()) {
                    require_near(electric.x, 0.0, 1e-12,
                                 "mixed-cell constant potential produced nonzero Ex");
                    require_near(electric.y, 0.0, 1e-12,
                                 "mixed-cell constant potential produced nonzero Ey");
                }
                require_throws(
                    [&]() {
                        (void)pic::solve_unstructured_poisson(
                            mixed_field_mesh,
                            {{"electrode", 1.0}, {"inlet", 0.0},
                             {"outlet", 0.0}, {"wall", 0.0}});
                    },
                    "unstructured Poisson solve accepted conflicting corner potentials");
            }

            pic::Gmsh2ImportLimits small_limits;
            small_limits.max_nodes = 5;
            require_throws([&]() { (void)pic::load_gmsh2_ascii_mesh2d(fixture, small_limits); },
                           "M2 Gmsh importer ignored its configured node limit");

            const auto bad_path = std::filesystem::path("test_output_bad_nonplanar.msh");
            {
                std::ofstream bad(bad_path);
                bad << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n"
                    << "$Nodes\n1\n1 0 0 0.25\n$EndNodes\n"
                    << "$Elements\n0\n$EndElements\n";
            }
            require_throws([&]() { (void)pic::load_gmsh2_ascii_mesh2d(bad_path); },
                           "M2 Gmsh importer accepted non-planar 2D node coordinates");
            std::filesystem::remove(bad_path);

            require_throws(
                []() {
                    pic::ImportedMesh2D invalid;
                    invalid.add_node({1, {std::numeric_limits<double>::infinity(), 0.0}});
                },
                "M2 imported mesh accepted a non-finite node");

            require_throws(
                []() {
                    pic::ImportedMesh2D invalid;
                    invalid.add_node({1, {0.0, 0.0}});
                    invalid.add_node({2, {1.0, 0.0}});
                    invalid.add_node({3, {0.0, 1.0}});
                    invalid.add_cell({1, pic::ImportedCellShape2D::Triangle, {1, 2, 3}, 1, "plasma"});
                    invalid.add_boundary_face({1, {1, 2}, 1, "wall"});
                },
                "M2 imported mesh accepted a duplicate global element id");

            require_throws(
                []() {
                    pic::ImportedMesh2D invalid;
                    invalid.add_node({1, {0.0, 0.0}});
                    invalid.add_node({2, {1.0, 0.0}});
                    invalid.add_node({3, {2.0, 0.0}});
                    invalid.add_cell({4, pic::ImportedCellShape2D::Triangle, {1, 2, 3}, 1, "plasma"});
                    invalid.add_boundary_face({5, {1, 2}, 1, "wall"});
                    invalid.add_boundary_face({6, {2, 3}, 1, "wall"});
                    invalid.add_boundary_face({7, {3, 1}, 1, "wall"});
                    invalid.validate();
                },
                "M2 imported mesh accepted a degenerate cell");

            require_throws(
                []() {
                    pic::ImportedMesh2D invalid;
                    invalid.add_node({1, {0.0, 0.0}});
                    invalid.add_node({2, {1.0, 0.0}});
                    invalid.add_node({3, {0.25, 0.25}});
                    invalid.add_node({4, {0.0, 1.0}});
                    invalid.add_cell({5, pic::ImportedCellShape2D::Quadrilateral, {1, 2, 3, 4}, 1, "plasma"});
                    invalid.add_boundary_face({6, {1, 2}, 1, "wall"});
                    invalid.add_boundary_face({7, {2, 3}, 1, "wall"});
                    invalid.add_boundary_face({8, {3, 4}, 1, "wall"});
                    invalid.add_boundary_face({9, {4, 1}, 1, "wall"});
                    invalid.validate();
                },
                "M2 imported mesh accepted a concave quadrilateral");

            require_throws(
                []() {
                    pic::ImportedMesh2D invalid;
                    invalid.add_node({1, {0.0, 0.0}});
                    invalid.add_node({2, {1.0, 0.0}});
                    invalid.add_node({3, {0.0, 1.0}});
                    invalid.add_cell({4, pic::ImportedCellShape2D::Triangle, {1, 2, 3}, 1, "plasma"});
                    invalid.add_boundary_face({5, {1, 2}, 1, "wall"});
                    invalid.add_boundary_face({6, {2, 3}, 1, "wall"});
                    invalid.validate();
                },
                "M2 imported mesh accepted an incomplete tagged boundary");

            require_throws(
                []() {
                    pic::ImportedMesh2D invalid;
                    invalid.add_node({1, {0.0, 0.0}});
                    invalid.add_node({2, {1.0, 0.0}});
                    invalid.add_node({3, {0.0, 1.0}});
                    invalid.add_node({4, {0.0, -1.0}});
                    invalid.add_node({5, {0.5, 1.0}});
                    invalid.add_cell({6, pic::ImportedCellShape2D::Triangle, {1, 2, 3}, 1, "plasma"});
                    invalid.add_cell({7, pic::ImportedCellShape2D::Triangle, {2, 1, 4}, 1, "plasma"});
                    invalid.add_cell({8, pic::ImportedCellShape2D::Triangle, {1, 2, 5}, 1, "plasma"});
                    invalid.add_boundary_face({9, {1, 3}, 1, "wall"});
                    invalid.validate();
                },
                "M2 imported mesh accepted a non-manifold cell edge");
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "aurorapic core test failure: " << e.what() << '\n';
        return 1;
    }
}
