#include "pic/Collision.hpp"
#include "pic/Config.hpp"
#include "pic/Diagnostics.hpp"
#include "pic/FieldSolver.hpp"
#include "pic/GasDataset.hpp"
#include "pic/Grid.hpp"
#include "pic/ImportedMesh2D.hpp"
#include "pic/Mesh2D.hpp"
#include "pic/Mesh3D.hpp"
#include "pic/PrescribedField.hpp"
#include "pic/Simulation.hpp"
#include "pic/Simulation2D.hpp"
#include "pic/Simulation3D.hpp"
#include "pic/Pusher.hpp"
#include "pic/Runtime.hpp"
#include "pic/VTKWriter.hpp"
#include "pic/Species.hpp"
#include "pic/Species2D.hpp"
#include "pic/Species3D.hpp"
#include "pic/Swarm.hpp"
#include "pic/UnstructuredFieldSolver2D.hpp"
#include "pic/UnstructuredMesh2D.hpp"
#include "pic/UnstructuredSimulation2D.hpp"
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
#include <set>
#include <sstream>
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
    require(a.live_particles_by_species == b.live_particles_by_species,
            label + ": per-species live-particle mismatch");
    require_near(a.kinetic_energy, b.kinetic_energy, 1e-12, label + ": kinetic-energy mismatch");
    require_near(a.field_energy, b.field_energy, 1e-12, label + ": field-energy mismatch");
    require_near(a.total_energy, b.total_energy, 1e-12, label + ": total-energy mismatch");
    require_near(a.charge_l1, b.charge_l1, 1e-12, label + ": charge-l1 mismatch");
    require_near(a.phi_left, b.phi_left, 1e-12, label + ": left-potential mismatch");
    require_near(a.phi_right, b.phi_right, 1e-12, label + ": right-potential mismatch");
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
            require_near(
                pa[i].velocity_y, pb[i].velocity_y, 1e-12,
                label + ": particle vy mismatch");
            require_near(
                pa[i].velocity_z, pb[i].velocity_z, 1e-12,
                label + ": particle vz mismatch");
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
            require_near(pa[i].velocity_z, pb[i].velocity_z, 1e-12, label + ": particle vz mismatch");
            require_near(pa[i].velocity_half.x, pb[i].velocity_half.x, 1e-12, label + ": particle vx_half mismatch");
            require_near(pa[i].velocity_half.y, pb[i].velocity_half.y, 1e-12, label + ": particle vy_half mismatch");
            require_near(pa[i].velocity_half_z, pb[i].velocity_half_z, 1e-12, label + ": particle vz_half mismatch");
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
            pic::Grid normalized(32, 1.0, pic::Boundary::Periodic);
            pic::Grid doubled_permittivity(
                32, 1.0, pic::Boundary::Periodic);
            for (std::size_t i = 0; i < normalized.nx(); ++i) {
                const double density =
                    std::sin(
                        2.0 * std::numbers::pi *
                        normalized.node_x(i));
                normalized.rho()[i] = density;
                doubled_permittivity.rho()[i] = density;
            }
            pic::FieldSolver normalized_solver;
            pic::FieldSolver scaled_solver(2.0);
            normalized_solver.solve(normalized);
            scaled_solver.solve(doubled_permittivity);
            for (std::size_t i = 0; i < normalized.nx(); ++i) {
                require_near(
                    doubled_permittivity.electric()[i],
                    0.5 * normalized.electric()[i], 1e-13,
                    "Poisson field did not scale inversely with permittivity");
            }
            pic::UnitSystemConfig si;
            si.system = pic::UnitSystem::SI;
            si.relative_permittivity = 2.5;
            require_near(
                si.permittivity(),
                2.5 * pic::VACUUM_PERMITTIVITY_SI, 1e-26,
                "SI unit contract computed the wrong permittivity");
            require_throws(
                []() { (void)pic::FieldSolver(0.0); },
                "field solver accepted zero permittivity");
        }
        {
            const auto table_path =
                std::filesystem::path("test_cross_section_table.dat");
            {
                std::ofstream table(table_path);
                table << "# energy cross_section\n"
                      << "0 0\n"
                      << "1 0.5\n"
                      << "2 1\n";
            }
            pic::CrossSectionTable table(table_path);
            require_near(
                table.evaluate(0.5), 0.25, 1e-15,
                "cross-section interpolation is incorrect");
            require_near(
                table.evaluate(3.0), 1.0, 1e-15,
                "cross-section upper endpoint behavior changed");
            pic::CrossSectionTable scaled_table(
                table_path, 2.0, 3.0);
            require_near(
                scaled_table.evaluate(1.0), 0.75, 1e-15,
                "cross-section column scaling is incorrect");
            const auto malformed_path =
                std::filesystem::path(
                    "test_cross_section_table_malformed.dat");
            {
                std::ofstream malformed(malformed_path);
                malformed << "0 0.1\n0 0.2\n";
            }
            require_throws(
                [&] {
                    (void)pic::CrossSectionTable(malformed_path);
                },
                "cross-section table accepted duplicate energies");
            std::filesystem::remove(malformed_path);

            pic::CollisionConfig elastic_config;
            elastic_config.enabled = true;
            elastic_config.model =
                pic::CollisionModelKind::NullCollision;
            elastic_config.neutral_density = 1.0;
            elastic_config.species = "test";
            elastic_config.max_frequency = 2.0;
            elastic_config.channels = {
                pic::CollisionChannelConfig{
                    "elastic",
                    pic::CollisionProcessKind::Elastic,
                    table_path, 0.0, 1.0, 1.0}};
            pic::NullCollisionModel elastic(elastic_config, 1.0);
            std::mt19937_64 rng(8128);
            std::uint64_t candidates = 0;
            std::uint64_t accepted = 0;
            for (std::size_t sample = 0; sample < 10000; ++sample) {
                double velocity = std::sqrt(2.0);
                const auto stats =
                    elastic.collide(velocity, 0.1, rng);
                candidates += stats.candidates;
                accepted += stats.channel_collisions[0];
                require_near(
                    std::abs(velocity), std::sqrt(2.0), 1e-14,
                    "elastic MCC did not preserve particle energy");
            }
            require(
                candidates > 1800 && candidates < 2200,
                "null-collision candidate rate left its statistical envelope");
            require(
                accepted > 620 && accepted < 800,
                "elastic collision acceptance left its statistical envelope");
            std::mt19937_64 vector_rng(8128);
            std::uint64_t vector_accepted = 0;
            double mean_direction_x = 0.0;
            double mean_direction_y = 0.0;
            double mean_direction_z = 0.0;
            for (std::size_t sample = 0; sample < 10000; ++sample) {
                pic::Vec3 velocity{std::sqrt(2.0), 0.0, 0.0};
                const auto stats =
                    elastic.collide(velocity, 0.1, vector_rng);
                vector_accepted += stats.channel_collisions[0];
                require_near(
                    norm(velocity), std::sqrt(2.0), 1e-14,
                    "2D3V elastic MCC did not preserve particle energy");
                if (stats.channel_collisions[0] != 0) {
                    mean_direction_x += velocity.x / std::sqrt(2.0);
                    mean_direction_y += velocity.y / std::sqrt(2.0);
                    mean_direction_z += velocity.z / std::sqrt(2.0);
                }
            }
            require(
                vector_accepted > 620 && vector_accepted < 800,
                "2D3V elastic collision acceptance left its statistical envelope");
            require(
                std::abs(mean_direction_x) <
                        0.1 * static_cast<double>(vector_accepted) &&
                    std::abs(mean_direction_y) <
                        0.1 * static_cast<double>(vector_accepted) &&
                    std::abs(mean_direction_z) <
                        0.1 * static_cast<double>(vector_accepted),
                "2D3V elastic MCC scattering is not statistically isotropic");

            const auto mean_cosine_path =
                std::filesystem::path(
                    "test_mean_cosine_table.dat");
            {
                std::ofstream angular(mean_cosine_path);
                angular << "0 0.8\n2 0.8\n";
            }
            pic::MeanCosineTable mean_cosines(mean_cosine_path);
            require_near(
                mean_cosines.evaluate(1.0), 0.8, 1e-15,
                "mean-cosine interpolation is incorrect");
            auto angular_config = elastic_config;
            angular_config.channels.front().angular_scattering =
                pic::AngularScatteringKind::HenyeyGreenstein;
            angular_config.channels.front().mean_cosine_file =
                mean_cosine_path;
            pic::NullCollisionModel angular_elastic(
                angular_config, 1.0);
            std::mt19937_64 angular_rng(161803);
            std::uint64_t single_angular_collisions = 0;
            double angular_cosine_sum = 0.0;
            for (std::size_t sample = 0; sample < 20000; ++sample) {
                pic::Vec3 velocity{std::sqrt(2.0), 0.0, 0.0};
                const auto stats = angular_elastic.collide(
                    velocity, 0.1, angular_rng);
                require_near(
                    norm(velocity), std::sqrt(2.0), 1e-14,
                    "anisotropic elastic MCC changed particle energy");
                if (stats.channel_collisions[0] != 1) continue;
                angular_cosine_sum +=
                    velocity.x / std::sqrt(2.0);
                ++single_angular_collisions;
            }
            require(
                single_angular_collisions > 1200,
                "anisotropic elastic MCC sampled too few collisions");
            require_near(
                angular_cosine_sum /
                    static_cast<double>(single_angular_collisions),
                0.8, 0.04,
                "Henyey-Greenstein scattering did not preserve "
                "the configured mean cosine");
            require_throws_contains(
                [&] {
                    double velocity = std::sqrt(2.0);
                    (void)angular_elastic.collide(
                        velocity, 0.1, angular_rng);
                },
                "3V collision interface",
                "scalar MCC accepted anisotropic scattering");
            const auto original_angular_signature =
                angular_elastic.signature();
            {
                std::ofstream changed(mean_cosine_path);
                changed << "0 0.7\n2 0.7\n";
            }
            pic::NullCollisionModel changed_angular(
                angular_config, 1.0);
            require(
                changed_angular.signature() !=
                    original_angular_signature,
                "MCC signature ignored changed angular data");
            {
                std::ofstream uncovered(mean_cosine_path);
                uncovered << "0 0.7\n1 0.7\n";
            }
            require_throws_contains(
                [&] {
                    (void)pic::NullCollisionModel(
                        angular_config, 1.0);
                },
                "must cover the cross-section table",
                "MCC accepted incomplete angular energy coverage");
            const auto invalid_mean_cosine_path =
                std::filesystem::path(
                    "test_invalid_mean_cosine_table.dat");
            {
                std::ofstream invalid(invalid_mean_cosine_path);
                invalid << "0 0\n1 1\n";
            }
            require_throws_contains(
                [&] {
                    (void)pic::MeanCosineTable(
                        invalid_mean_cosine_path);
                },
                "strictly between -1 and 1",
                "mean-cosine table accepted a singular endpoint");
            std::filesystem::remove(invalid_mean_cosine_path);
            std::filesystem::remove(mean_cosine_path);

            auto finite_mass_config = elastic_config;
            finite_mass_config.neutral_mass = 3.0;
            pic::NullCollisionModel finite_mass_elastic(
                finite_mass_config, 1.0);
            std::mt19937_64 recoil_rng(57721);
            bool observed_recoil = false;
            for (std::size_t attempt = 0;
                 attempt < 10000 && !observed_recoil; ++attempt) {
                const pic::Vec3 initial{std::sqrt(2.0), 0.0, 0.0};
                pic::Vec3 projectile = initial;
                const auto recoil_stats =
                    finite_mass_elastic.collide(
                        projectile, 0.01, recoil_rng);
                if (recoil_stats.channel_collisions[0] != 1) continue;
                const pic::Vec3 neutral{
                    (initial.x - projectile.x) / 3.0,
                    (initial.y - projectile.y) / 3.0,
                    (initial.z - projectile.z) / 3.0};
                require_near(
                    projectile.x + 3.0 * neutral.x,
                    initial.x, 1e-14,
                    "finite-mass elastic collision lost x momentum");
                require_near(
                    projectile.y + 3.0 * neutral.y,
                    initial.y, 1e-14,
                    "finite-mass elastic collision lost y momentum");
                require_near(
                    projectile.z + 3.0 * neutral.z,
                    initial.z, 1e-14,
                    "finite-mass elastic collision lost z momentum");
                require_near(
                    0.5 * norm(projectile) * norm(projectile) +
                        1.5 * norm(neutral) * norm(neutral),
                    1.0, 1e-13,
                    "finite-mass elastic collision lost total kinetic energy");
                observed_recoil = true;
            }
            require(
                observed_recoil,
                "finite-mass elastic collision was not sampled");

            auto center_of_mass_config = elastic_config;
            center_of_mass_config.neutral_mass = 1.0;
            center_of_mass_config.channels.front().energy_frame =
                pic::CollisionEnergyFrame::CenterOfMass;
            pic::NullCollisionModel center_of_mass_elastic(
                center_of_mass_config, 1.0);
            std::mt19937_64 center_of_mass_rng(141421);
            std::uint64_t center_of_mass_accepted = 0;
            for (std::size_t sample = 0; sample < 10000; ++sample) {
                pic::Vec3 projectile{std::sqrt(2.0), 0.0, 0.0};
                const auto stats = center_of_mass_elastic.collide(
                    projectile, 0.1, center_of_mass_rng);
                center_of_mass_accepted +=
                    stats.channel_collisions[0];
            }
            require(
                center_of_mass_accepted > 280 &&
                    center_of_mass_accepted < 440,
                "center-of-mass energy did not select the expected "
                "cross section");
            auto projectile_frame_config = center_of_mass_config;
            projectile_frame_config.channels.front().energy_frame =
                pic::CollisionEnergyFrame::Projectile;
            pic::NullCollisionModel projectile_frame_elastic(
                projectile_frame_config, 1.0);
            require(
                center_of_mass_elastic.signature() !=
                    projectile_frame_elastic.signature(),
                "MCC signature ignored the collision energy frame");

            auto backward_config = center_of_mass_config;
            backward_config.channels.front().angular_scattering =
                pic::AngularScatteringKind::Backward;
            pic::NullCollisionModel backward_elastic(
                backward_config, 1.0);
            std::mt19937_64 backward_rng(173205);
            bool observed_backward = false;
            for (std::size_t attempt = 0;
                 attempt < 10000 && !observed_backward; ++attempt) {
                pic::Vec3 projectile{std::sqrt(2.0), 0.0, 0.0};
                const auto stats = backward_elastic.collide(
                    projectile, 0.1, backward_rng);
                if (stats.channel_collisions[0] != 1) continue;
                require_near(
                    norm(projectile), 0.0, 1e-14,
                    "equal-mass backward scattering did not exchange "
                    "the projectile and neutral velocities");
                observed_backward = true;
            }
            require(
                observed_backward,
                "backward elastic collision was not sampled");
            require_throws_contains(
                [&] {
                    double projectile = std::sqrt(2.0);
                    (void)backward_elastic.collide(
                        projectile, 0.1, backward_rng);
                },
                "3V collision interface",
                "scalar MCC accepted backward scattering");

            auto excitation_config = elastic_config;
            excitation_config.max_frequency = 4.0;
            excitation_config.channels = {
                pic::CollisionChannelConfig{
                    "excitation",
                    pic::CollisionProcessKind::Excitation,
                    table_path, 0.5, 1.0, 2.0}};
            pic::NullCollisionModel excitation(
                excitation_config, 1.0);
            double excited_velocity = 2.0;
            const double initial_energy =
                0.5 * excited_velocity * excited_velocity;
            const auto excitation_stats =
                excitation.collide(excited_velocity, 2.0, rng);
            require(
                excitation_stats.channel_collisions[0] > 0,
                "excitation MCC did not accept a collision");
            const double final_energy =
                0.5 * excited_velocity * excited_velocity;
            require_near(
                initial_energy - final_energy,
                0.5 * static_cast<double>(
                          excitation_stats.channel_collisions[0]),
                1e-13,
                "excitation MCC removed the wrong threshold energy");
            pic::Vec3 excited_vector{2.0, 0.0, 0.0};
            const double initial_vector_energy =
                0.5 * norm(excited_vector) * norm(excited_vector);
            const auto vector_excitation_stats =
                excitation.collide(excited_vector, 2.0, vector_rng);
            require(
                vector_excitation_stats.channel_collisions[0] > 0,
                "2D3V excitation MCC did not accept a collision");
            require_near(
                initial_vector_energy -
                    0.5 * norm(excited_vector) * norm(excited_vector),
                0.5 * static_cast<double>(
                          vector_excitation_stats.channel_collisions[0]),
                1e-13,
                "2D3V excitation MCC removed the wrong threshold energy");

            auto ionization_config = elastic_config;
            ionization_config.max_frequency = 10.0;
            ionization_config.max_candidates_per_particle = 256;
            ionization_config.channels = {
                pic::CollisionChannelConfig{
                    "ionization",
                    pic::CollisionProcessKind::Ionization,
                    table_path, 0.5, 1.0, 1.0,
                    "electrons", "ions"}};
            pic::NullCollisionModel ionization(
                ionization_config, 1.0);
            pic::Vec3 ionizing_velocity{4.0, 0.0, 0.0};
            const double ionizing_initial_energy =
                0.5 * norm(ionizing_velocity) *
                norm(ionizing_velocity);
            std::mt19937_64 ionization_rng(314159);
            const auto ionization_stats =
                ionization.collide(
                    ionizing_velocity, 10.0, ionization_rng);
            require(
                ionization_stats.channel_collisions[0] == 4 &&
                    ionization_stats.secondaries.size() == 4,
                "2D3V ionization did not create one secondary per event");
            double product_energy =
                0.5 * norm(ionizing_velocity) *
                norm(ionizing_velocity);
            for (const auto& secondary :
                 ionization_stats.secondaries) {
                require(
                    secondary.channel == 0,
                    "ionization secondary channel mapping changed");
                product_energy +=
                    0.5 * norm(secondary.velocity) *
                    norm(secondary.velocity);
            }
            require_near(
                product_energy +
                    0.5 * static_cast<double>(
                              ionization_stats.secondaries.size()),
                ionizing_initial_energy, 1e-12,
                "ionization energy partition is not conservative");
            auto thermal_ionization_config =
                ionization_config;
            thermal_ionization_config.gas_data_units =
                pic::UnitSystem::SI;
            thermal_ionization_config.neutral_mass = 1.0;
            thermal_ionization_config.neutral_temperature =
                1.0 / 1.380649e-23;
            thermal_ionization_config.max_frequency = 20.0;
            pic::NullCollisionModel thermal_ionization(
                thermal_ionization_config, 1.0);
            std::mt19937_64 thermal_ionization_rng(223607);
            bool thermal_ionized = false;
            for (std::size_t attempt = 0;
                 attempt < 10000 && !thermal_ionized; ++attempt) {
                pic::Vec3 electron{4.0, 0.0, 0.0};
                const auto stats = thermal_ionization.collide(
                    electron, 0.01, thermal_ionization_rng);
                if (stats.secondaries.empty()) continue;
                require(
                    norm(stats.secondaries.front().ion_velocity) >
                            0.0 &&
                        norm(stats.secondaries.front().ion_velocity) <=
                            8.0,
                    "thermal ionization product did not inherit "
                    "the target-neutral velocity");
                thermal_ionized = true;
            }
            require(
                thermal_ionized,
                "thermal ionization was not sampled");
            require_throws(
                [&] {
                    double velocity = 4.0;
                    (void)ionization.collide(
                        velocity, 0.1, ionization_rng);
                },
                "scalar MCC accepted an ionization channel");

            auto attachment_config = elastic_config;
            pic::CollisionChannelConfig attachment_channel;
            attachment_channel.name = "attachment";
            attachment_channel.process =
                pic::CollisionProcessKind::Attachment;
            attachment_channel.cross_section_file = table_path;
            attachment_channel.attachment_species = "negative_ions";
            attachment_config.channels = {attachment_channel};
            pic::NullCollisionModel attachment(
                attachment_config, 1.0);
            std::mt19937_64 attachment_rng(271828);
            bool attached = false;
            for (std::size_t attempt = 0;
                 attempt < 10000 && !attached; ++attempt) {
                pic::Vec3 electron{std::sqrt(2.0), 0.0, 0.0};
                const auto stats = attachment.collide(
                    electron, 0.1, attachment_rng);
                if (!stats.primary_removal_channel) continue;
                require(
                    *stats.primary_removal_channel == 0 &&
                        stats.channel_collisions[0] == 1 &&
                        stats.secondaries.empty() &&
                        norm(electron) == 0.0,
                    "attachment did not consume exactly one primary");
                attached = true;
            }
            require(
                attached,
                "electron attachment was not sampled");
            auto thermal_attachment_config =
                attachment_config;
            thermal_attachment_config.gas_data_units =
                pic::UnitSystem::SI;
            thermal_attachment_config.neutral_mass = 1.0;
            thermal_attachment_config.neutral_temperature =
                1.0 / 1.380649e-23;
            thermal_attachment_config.max_frequency = 20.0;
            pic::NullCollisionModel thermal_attachment(
                thermal_attachment_config, 1.0);
            std::mt19937_64 thermal_attachment_rng(244949);
            bool thermal_attached = false;
            for (std::size_t attempt = 0;
                 attempt < 10000 && !thermal_attached; ++attempt) {
                pic::Vec3 electron{2.0, 0.0, 0.0};
                const auto stats = thermal_attachment.collide(
                    electron, 0.01, thermal_attachment_rng);
                if (!stats.primary_removal_channel) continue;
                require(
                    stats.primary_removal_product_velocity &&
                        norm(*stats.primary_removal_product_velocity) >
                            0.0 &&
                        norm(*stats.primary_removal_product_velocity) <=
                            8.0,
                    "thermal attachment product did not inherit "
                    "the target-neutral velocity");
                thermal_attached = true;
            }
            require(
                thermal_attached,
                "thermal attachment was not sampled");
            auto unmapped_attachment = attachment_config;
            unmapped_attachment.channels.front()
                .attachment_species.clear();
            require_throws_contains(
                [&] {
                    (void)pic::NullCollisionModel(
                        unmapped_attachment, 1.0);
                },
                "requires an attachment species",
                "attachment accepted a missing product mapping");

            auto charge_exchange_config = elastic_config;
            charge_exchange_config.neutral_mass = 40.0;
            charge_exchange_config.channels = {
                pic::CollisionChannelConfig{
                    "charge_exchange",
                    pic::CollisionProcessKind::ChargeExchange,
                    table_path, 0.0, 1.0, 1.0}};
            pic::NullCollisionModel charge_exchange(
                charge_exchange_config, 40.0);
            std::mt19937_64 charge_exchange_rng(141421);
            bool exchanged = false;
            for (std::size_t attempt = 0;
                 attempt < 10000 && !exchanged; ++attempt) {
                pic::Vec3 ion{1.0, 0.0, 0.0};
                const auto stats = charge_exchange.collide(
                    ion, 0.01, charge_exchange_rng);
                if (stats.channel_collisions[0] == 0) continue;
                require(
                    stats.channel_collisions[0] == 1 &&
                        norm(ion) == 0.0 &&
                        stats.secondaries.empty(),
                    "resonant charge exchange did not create a slow ion");
                exchanged = true;
            }
            require(
                exchanged,
                "resonant charge exchange was not sampled");
            auto thermal_exchange_config =
                charge_exchange_config;
            thermal_exchange_config.gas_data_units =
                pic::UnitSystem::SI;
            thermal_exchange_config.neutral_mass = 1.0;
            thermal_exchange_config.neutral_temperature =
                1.0 / 1.380649e-23;
            thermal_exchange_config.max_frequency = 20.0;
            pic::NullCollisionModel thermal_exchange(
                thermal_exchange_config, 1.0);
            require_near(
                thermal_exchange.neutral_velocity_stddev(),
                1.0, 1e-15,
                "SI neutral temperature produced the wrong "
                "Maxwellian component speed");
            require(
                thermal_exchange.signature() !=
                    charge_exchange.signature(),
                "MCC signature ignored thermal-neutral kinematics");
            std::mt19937_64 thermal_exchange_rng(173205);
            bool thermal_exchanged = false;
            for (std::size_t attempt = 0;
                 attempt < 10000 && !thermal_exchanged; ++attempt) {
                pic::Vec3 ion{2.0, 0.0, 0.0};
                const auto stats = thermal_exchange.collide(
                    ion, 0.01, thermal_exchange_rng);
                if (stats.channel_collisions[0] == 0) continue;
                require(
                    norm(ion) > 0.0 && norm(ion) <= 8.0 &&
                        stats.secondaries.empty(),
                    "thermal resonant charge exchange did not map "
                    "the sampled neutral velocity onto the ion");
                thermal_exchanged = true;
            }
            require(
                thermal_exchanged,
                "thermal resonant charge exchange was not sampled");
            auto unsafe_thermal_exchange =
                thermal_exchange_config;
            unsafe_thermal_exchange.max_frequency = 1.0;
            pic::NullCollisionModel unsafe_thermal(
                unsafe_thermal_exchange, 1.0);
            require_throws_contains(
                [&] {
                    pic::Vec3 ion{2.0, 0.0, 0.0};
                    (void)unsafe_thermal.collide(
                        ion, 0.01, thermal_exchange_rng);
                },
                "thermal-neutral collision-frequency bound",
                "thermal MCC accepted an unsafe null-collision "
                "majorant");
            auto mismatched_exchange_config =
                charge_exchange_config;
            require_throws_contains(
                [&] {
                    (void)pic::NullCollisionModel(
                        mismatched_exchange_config, 39.0);
                },
                "projectile mass equal to neutral_mass",
                "resonant charge exchange accepted mismatched masses");

            auto unsafe_config = elastic_config;
            unsafe_config.max_frequency = 0.1;
            pic::NullCollisionModel unsafe(unsafe_config, 1.0);
            require_throws(
                [&] {
                    double velocity = std::sqrt(2.0);
                    (void)unsafe.collide(velocity, 0.1, rng);
                },
                "MCC accepted a total rate above max_frequency");
            std::filesystem::remove(table_path);
        }
        {
            const auto table_path =
                std::filesystem::absolute(
                    "test_mcc_checkpoint_table.dat");
            const auto output_dir =
                std::filesystem::path("test_output_mcc_checkpoint");
            const auto checkpoint_path =
                output_dir / "manual.apc";
            std::filesystem::remove_all(output_dir);
            {
                std::ofstream table(table_path);
                table << "0 0.2\n10 0.2\n";
            }

            pic::Config config;
            config.nx = 16;
            config.length = 1.0;
            config.dt = 0.02;
            config.steps = 5;
            config.output_interval = 5;
            config.output_dir = output_dir.string();
            config.seed = 441;
            config.collisions.enabled = true;
            config.collisions.model =
                pic::CollisionModelKind::NullCollision;
            config.collisions.species = "tracer";
            config.collisions.neutral_density = 2.0;
            config.collisions.max_frequency = 1.0;
            config.collisions.channels = {
                pic::CollisionChannelConfig{
                    "elastic",
                    pic::CollisionProcessKind::Elastic,
                    table_path, 0.0, 1.0, 1.0}};
            config.species = {
                pic::SpeciesConfig{
                    "tracer", 0.0, 1.0, 1.0, 128, 1.0,
                    1.0, 0.05, 0.0, -1.0}};

            pic::Simulation continuous(config);
            continuous.initialize();
            continuous.step();
            continuous.step();
            continuous.save_checkpoint(checkpoint_path);
            for (std::size_t step = continuous.step_count();
                 step < config.steps; ++step) {
                continuous.step();
            }

            pic::Simulation restarted(config);
            restarted.load_checkpoint(checkpoint_path);
            for (std::size_t step = restarted.step_count();
                 step < config.steps; ++step) {
                restarted.step();
            }
            require_species_close(
                continuous.species(), restarted.species(),
                "MCC checkpoint restart determinism");
            const auto& expected =
                continuous.collision_diagnostics();
            const auto& actual =
                restarted.collision_diagnostics();
            require(
                expected.candidates == actual.candidates &&
                    expected.null_collisions ==
                        actual.null_collisions &&
                    expected.channel_collisions ==
                        actual.channel_collisions,
                "MCC checkpoint did not preserve collision diagnostics");

            {
                std::ofstream changed_table(table_path);
                changed_table << "0 0.1\n10 0.1\n";
            }
            require_throws(
                [&] {
                    pic::Simulation changed(config);
                    changed.load_checkpoint(checkpoint_path);
                },
                "MCC checkpoint accepted changed cross-section data");
            std::filesystem::remove_all(output_dir);
            std::filesystem::remove(table_path);
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
            const pic::TabulatedVectorField1D profile(
                pic::CoordinateAxis::X,
                {0.0, 0.5, 1.0},
                {{0.0, 1.0, 2.0},
                 {1.0, 3.0, 4.0},
                 {2.0, 5.0, 6.0}});
            const auto interpolated =
                profile.evaluate({0.25, 9.0, -4.0});
            require_near(
                interpolated.x, 0.5, 1e-15,
                "tabulated vector field x interpolation mismatch");
            require_near(
                interpolated.y, 2.0, 1e-15,
                "tabulated vector field y interpolation mismatch");
            require_near(
                interpolated.z, 3.0, 1e-15,
                "tabulated vector field z interpolation mismatch");
            profile.validate_domain(
                {0.1, 0.0, 0.0}, {0.9, 1.0, 1.0},
                "test profile");
            require(
                pic::parse_coordinate_axis("Y") ==
                    pic::CoordinateAxis::Y,
                "coordinate-axis parser is not case insensitive");
            require_throws(
                [&]() { (void)profile.evaluate({-0.01, 0.0, 0.0}); },
                "tabulated vector field extrapolated below its range");
            require_throws(
                [&]() {
                    profile.validate_domain(
                        {-0.1, 0.0, 0.0}, {0.9, 1.0, 1.0},
                        "undersized profile");
                },
                "tabulated vector field accepted incomplete domain coverage");
            require_throws(
                []() {
                    (void)pic::TabulatedVectorField1D(
                        pic::CoordinateAxis::X,
                        {0.0, 0.5, 0.5},
                        {{}, {}, {}});
                },
                "tabulated vector field accepted duplicate coordinates");

            const auto profile_path =
                std::filesystem::path("test_magnetic_profile.dat");
            const auto config_path =
                std::filesystem::path("test_magnetic_profile.cfg");
            {
                std::ofstream output(profile_path);
                output << "# x Bx By Bz\n"
                       << "0 0 0 1\n"
                       << "0.5 0 0 2\n"
                       << "1 0 0 3\n";
            }
            {
                std::ofstream output(config_path);
                output
                    << "config_version = 1\n"
                    << "dimension = 2\n"
                    << "nx = 8\nny = 8\n"
                    << "length_x = 1\nlength_y = 1\n"
                    << "boundary = dirichlet\n"
                    << "boundary_x = dirichlet\n"
                    << "boundary_y = periodic\n"
                    << "magnetic_field_profile_file = "
                    << profile_path.string() << "\n"
                    << "magnetic_field_profile_axis = x\n";
            }
            const auto loaded = pic::load_config_2d(config_path.string());
            require(
                loaded.magnetic_field_profile.has_value(),
                "2D config did not load magnetic field profile");
            require(
                loaded.boundary_x == pic::Boundary::Dirichlet &&
                    loaded.boundary_y == pic::Boundary::Periodic,
                "2D config did not load per-axis boundary topology");
            require_near(
                loaded.magnetic_field_profile
                    ->evaluate({0.75, 0.0, 0.0}).z,
                2.5, 1e-15,
                "2D config magnetic field profile interpolation mismatch");
            pic::Simulation2D mixed_topology_simulation(loaded);
            require(
                mixed_topology_simulation.particle_boundary_config().left ==
                        pic::ParticleBoundary::Absorbing &&
                    mixed_topology_simulation.particle_boundary_config().right ==
                        pic::ParticleBoundary::Absorbing &&
                    mixed_topology_simulation.particle_boundary_config().bottom ==
                        pic::ParticleBoundary::Periodic &&
                    mixed_topology_simulation.particle_boundary_config().top ==
                        pic::ParticleBoundary::Periodic,
                "2D mixed field topology did not resolve default particle boundaries per axis");

            auto conflicting = loaded;
            conflicting.magnetic_field_z = 1.0;
            require_throws(
                [&]() { pic::Simulation2D simulation(conflicting); },
                "2D simulation accepted uniform and profiled magnetic fields together");

            pic::Simulation3DConfig config_3d;
            config_3d.magnetic_field_profile =
                pic::TabulatedVectorField1D(
                    pic::CoordinateAxis::Y,
                    {0.0, 1.0},
                    {{0.0, 0.0, 1.0}, {0.0, 0.0, 2.0}});
            pic::Simulation3D simulation_3d(config_3d);
            (void)simulation_3d;

            std::filesystem::remove(profile_path);
            std::filesystem::remove(config_path);
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
            pic::Particle2D particle_2d{};
            particle_2d.velocity = pic::Vec2{initial_velocity.x, initial_velocity.y};
            particle_2d.velocity_z = initial_velocity.z;
            pic::Particle3D particle{};
            particle.velocity = initial_velocity;
            const pic::Vec2 electric_2d{0.0, 0.0};
            const pic::Vec3 electric{0.0, 0.0, 0.0};
            pic::initialize_boris_half_step(
                particle_2d, electric_2d, magnetic, charge_to_mass, dt);
            pic::initialize_boris_half_step(particle, electric, magnetic, charge_to_mass, dt);
            for (std::size_t n = 0; n < steps; ++n) {
                pic::kick_boris(
                    particle_2d, electric_2d, magnetic, charge_to_mass, dt);
                pic::kick_boris(particle, electric, magnetic, charge_to_mass, dt);
            }
            pic::synchronize_boris(
                particle_2d, electric_2d, magnetic, charge_to_mass, dt);
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
            require_near(particle_2d.velocity.x, particle.velocity.x, 1e-13,
                         "2D3V Boris pusher disagrees with 3D arbitrary-axis x velocity");
            require_near(particle_2d.velocity.y, particle.velocity.y, 1e-13,
                         "2D3V Boris pusher disagrees with 3D arbitrary-axis y velocity");
            require_near(particle_2d.velocity_z, particle.velocity.z, 1e-13,
                         "2D3V Boris pusher disagrees with 3D arbitrary-axis z velocity");
            require_near(
                norm(pic::Vec3{particle_2d.velocity.x, particle_2d.velocity.y,
                               particle_2d.velocity_z}),
                norm(initial_velocity), 1e-13,
                "2D3V Boris pusher did not conserve speed");
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
            pic::BoundaryConfig2D electrodes;
            electrodes.left.potential = 200.0;
            electrodes.right.potential = 0.0;
            pic::Mesh2D mesh(
                17, 8, 0.025, 0.0128,
                pic::Boundary::Dirichlet,
                pic::Boundary::Periodic,
                electrodes);
            require(
                mesh.boundary_x() == pic::Boundary::Dirichlet &&
                    mesh.boundary_y() == pic::Boundary::Periodic,
                "mixed 2D mesh did not retain per-axis topology");
            require_near(
                mesh.dx(), mesh.length_x() /
                    static_cast<double>(mesh.nx() - 1),
                1e-15,
                "mixed 2D mesh used periodic spacing on Dirichlet x");
            require_near(
                mesh.dy(), mesh.length_y() /
                    static_cast<double>(mesh.ny()),
                1e-15,
                "mixed 2D mesh used endpoint spacing on periodic y");
            require_throws(
                [&]() { (void)mesh.boundary(); },
                "mixed 2D mesh global boundary compatibility accessor did not reject ambiguity");

            pic::FieldSolver solver;
            solver.solve(mesh);
            double maximum_phi_error = 0.0;
            double maximum_ex_error = 0.0;
            double maximum_ey = 0.0;
            for (std::size_t j = 0; j < mesh.ny(); ++j) {
                for (std::size_t i = 0; i < mesh.nx(); ++i) {
                    const auto index = mesh.index(i, j);
                    const double expected_phi =
                        200.0 * (1.0 - mesh.node_x(i) /
                                          mesh.length_x());
                    maximum_phi_error = std::max(
                        maximum_phi_error,
                        std::abs(mesh.phi()[index] - expected_phi));
                    maximum_ex_error = std::max(
                        maximum_ex_error,
                        std::abs(mesh.electric_x()[index] -
                                 200.0 / mesh.length_x()));
                    maximum_ey = std::max(
                        maximum_ey,
                        std::abs(mesh.electric_y()[index]));
                }
            }
            require(
                maximum_phi_error < 1e-8,
                "mixed Dirichlet-x/periodic-y Poisson potential exceeded analytic tolerance");
            require(
                maximum_ex_error < 1e-6,
                "mixed Dirichlet-x/periodic-y Poisson Ex exceeded analytic tolerance");
            require(
                maximum_ey < 1e-6,
                "mixed Dirichlet-x/periodic-y Poisson produced spurious Ey");

            const std::vector<pic::Particle2D> particles{
                pic::Particle2D{
                    pic::Vec2{0.5 * mesh.length_x(),
                              1.25 * mesh.length_y()},
                    pic::Vec2{}, true}};
            mesh.clear_charge();
            pic::deposit_charge_cic(mesh, particles, 2.0, 0.5);
            double represented_charge = 0.0;
            for (std::size_t j = 0; j < mesh.ny(); ++j) {
                for (std::size_t i = 0; i < mesh.nx(); ++i) {
                    represented_charge +=
                        mesh.rho()[mesh.index(i, j)] *
                        mesh.node_area(i, j);
                }
            }
            require_near(
                represented_charge, 1.0, 1e-12,
                "mixed-topology CIC deposition did not conserve charge");
        }
        {
            pic::BoundaryConfig2D electrodes;
            electrodes.left.potential = 200.0;
            electrodes.right.potential = 0.0;
            pic::Mesh2D mesh(
                501, 256, 0.025, 0.0128,
                pic::Boundary::Dirichlet,
                pic::Boundary::Periodic,
                electrodes);
            pic::FieldSolver solver;
            solver.solve(mesh);
            double maximum_phi_error = 0.0;
            for (std::size_t j = 0; j < mesh.ny(); ++j) {
                for (std::size_t i = 0; i < mesh.nx(); ++i) {
                    const double expected =
                        200.0 * (1.0 - mesh.node_x(i) /
                                          mesh.length_x());
                    maximum_phi_error = std::max(
                        maximum_phi_error,
                        std::abs(
                            mesh.phi()[mesh.index(i, j)] -
                            expected));
                }
            }
            require(
                maximum_phi_error < 1e-10,
                "500-by-256-cell LANDMARK mixed spectral-tridiagonal vacuum solve exceeded tolerance");
        }
        {
            pic::Mesh2D mesh(
                18, 37, 1.2, 0.8,
                pic::Boundary::Dirichlet,
                pic::Boundary::Periodic);
            const double kx =
                std::numbers::pi / mesh.length_x();
            const double ky =
                2.0 * std::numbers::pi / mesh.length_y();
            const double inv_dx2 =
                1.0 / (mesh.dx() * mesh.dx());
            const double inv_dy2 =
                1.0 / (mesh.dy() * mesh.dy());
            std::vector<double> expected(mesh.size(), 0.0);
            for (std::size_t j = 0; j < mesh.ny(); ++j) {
                for (std::size_t i = 0; i < mesh.nx(); ++i) {
                    expected[mesh.index(i, j)] =
                        std::sin(kx * mesh.node_x(i)) *
                        std::cos(ky * mesh.node_y(j));
                }
            }
            for (std::size_t j = 0; j < mesh.ny(); ++j) {
                const std::size_t jm =
                    (j + mesh.ny() - 1) % mesh.ny();
                const std::size_t jp =
                    (j + 1) % mesh.ny();
                for (std::size_t i = 1;
                     i + 1 < mesh.nx(); ++i) {
                    const auto index = mesh.index(i, j);
                    mesh.rho()[index] = pic::EPS0 *
                        ((2.0 * inv_dx2 + 2.0 * inv_dy2) *
                             expected[index]
                         - inv_dx2 *
                             (expected[mesh.index(i - 1, j)] +
                              expected[mesh.index(i + 1, j)])
                         - inv_dy2 *
                             (expected[mesh.index(i, jm)] +
                              expected[mesh.index(i, jp)]));
                }
            }

            pic::FieldSolver solver;
            solver.solve(mesh);
            double maximum_phi_error = 0.0;
            double maximum_ex_error = 0.0;
            double maximum_ey_error = 0.0;
            for (std::size_t j = 0; j < mesh.ny(); ++j) {
                const std::size_t jm =
                    (j + mesh.ny() - 1) % mesh.ny();
                const std::size_t jp =
                    (j + 1) % mesh.ny();
                for (std::size_t i = 0; i < mesh.nx(); ++i) {
                    const auto index = mesh.index(i, j);
                    maximum_phi_error = std::max(
                        maximum_phi_error,
                        std::abs(mesh.phi()[index] -
                                 expected[index]));
                    const double expected_ex =
                        i == 0
                            ? -(expected[mesh.index(1, j)] -
                                expected[index]) / mesh.dx()
                            : i + 1 == mesh.nx()
                                  ? -(expected[index] -
                                      expected[mesh.index(i - 1, j)]) /
                                        mesh.dx()
                                  : -(expected[mesh.index(i + 1, j)] -
                                      expected[mesh.index(i - 1, j)]) /
                                        (2.0 * mesh.dx());
                    const double expected_ey =
                        -(expected[mesh.index(i, jp)] -
                          expected[mesh.index(i, jm)]) /
                        (2.0 * mesh.dy());
                    maximum_ex_error = std::max(
                        maximum_ex_error,
                        std::abs(mesh.electric_x()[index] -
                                 expected_ex));
                    maximum_ey_error = std::max(
                        maximum_ey_error,
                        std::abs(mesh.electric_y()[index] -
                                 expected_ey));
                }
            }
            require(
                maximum_phi_error < 1e-11,
                "mixed spectral-tridiagonal Poisson potential exceeded tolerance");
            require(
                maximum_ex_error < 1e-10,
                "mixed spectral-tridiagonal Poisson Ex exceeded tolerance");
            require(
                maximum_ey_error < 1e-10,
                "mixed spectral-tridiagonal Poisson Ey exceeded tolerance");
        }
        {
            pic::BoundaryConfig2D electrodes;
            electrodes.bottom.potential = -1.0;
            electrodes.top.potential = 2.0;
            pic::Mesh2D mesh(
                10, 13, 1.0, 0.6,
                pic::Boundary::Periodic,
                pic::Boundary::Dirichlet,
                electrodes);
            const double kx =
                4.0 * std::numbers::pi / mesh.length_x();
            const double ky =
                std::numbers::pi / mesh.length_y();
            const double inv_dx2 =
                1.0 / (mesh.dx() * mesh.dx());
            const double inv_dy2 =
                1.0 / (mesh.dy() * mesh.dy());
            std::vector<double> expected(mesh.size(), 0.0);
            for (std::size_t j = 0; j < mesh.ny(); ++j) {
                const double linear =
                    -1.0 + 3.0 * mesh.node_y(j) /
                               mesh.length_y();
                for (std::size_t i = 0; i < mesh.nx(); ++i) {
                    expected[mesh.index(i, j)] =
                        linear +
                        0.4 * std::cos(kx * mesh.node_x(i)) *
                            std::sin(ky * mesh.node_y(j));
                }
            }
            for (std::size_t j = 1;
                 j + 1 < mesh.ny();
                 ++j) {
                for (std::size_t i = 0; i < mesh.nx(); ++i) {
                    const std::size_t im =
                        (i + mesh.nx() - 1) % mesh.nx();
                    const std::size_t ip =
                        (i + 1) % mesh.nx();
                    const auto index = mesh.index(i, j);
                    mesh.rho()[index] = pic::EPS0 *
                        ((2.0 * inv_dx2 + 2.0 * inv_dy2) *
                             expected[index]
                         - inv_dx2 *
                             (expected[mesh.index(im, j)] +
                              expected[mesh.index(ip, j)])
                         - inv_dy2 *
                             (expected[mesh.index(i, j - 1)] +
                              expected[mesh.index(i, j + 1)]));
                }
            }
            pic::FieldSolver solver;
            solver.solve(mesh);
            for (std::size_t j = 0; j < mesh.ny(); ++j) {
                for (std::size_t i = 0; i < mesh.nx(); ++i) {
                    const auto index = mesh.index(i, j);
                    const std::size_t im =
                        (i + mesh.nx() - 1) % mesh.nx();
                    const std::size_t ip =
                        (i + 1) % mesh.nx();
                    const double expected_ex =
                        -(expected[mesh.index(ip, j)] -
                          expected[mesh.index(im, j)]) /
                        (2.0 * mesh.dx());
                    const double expected_ey =
                        j == 0
                            ? -(expected[mesh.index(i, 1)] -
                                expected[index]) / mesh.dy()
                            : j + 1 == mesh.ny()
                                  ? -(expected[index] -
                                      expected[mesh.index(i, j - 1)]) /
                                        mesh.dy()
                                  : -(expected[mesh.index(i, j + 1)] -
                                      expected[mesh.index(i, j - 1)]) /
                                        (2.0 * mesh.dy());
                    require_near(
                        mesh.phi()[index], expected[index], 1e-11,
                        "mixed periodic-x/Dirichlet-y Poisson potential mismatch");
                    require_near(
                        mesh.electric_x()[index], expected_ex, 1e-10,
                        "mixed periodic-x/Dirichlet-y Poisson Ex mismatch");
                    require_near(
                        mesh.electric_y()[index], expected_ey, 1e-10,
                        "mixed periodic-x/Dirichlet-y Poisson Ey mismatch");
                }
            }
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
            pic::Mesh2D mesh(
                5, 4, 1.0, 2.0,
                pic::Boundary::Dirichlet,
                pic::Boundary::Periodic);
            for (std::size_t j = 0; j < mesh.ny(); ++j) {
                const double phase =
                    2.0 * std::numbers::pi *
                    mesh.node_y(j) / mesh.length_y();
                for (std::size_t i = 0; i < mesh.nx(); ++i) {
                    const auto index = mesh.index(i, j);
                    mesh.phi()[index] = 2.0 * mesh.node_x(i);
                    mesh.rho()[index] =
                        5.0 + 2.0 * std::cos(phase);
                    mesh.electric_x()[index] =
                        3.0 * std::sin(phase);
                    mesh.electric_y()[index] =
                        -4.0 * std::cos(phase);
                }
            }
            pic::Species2DConfig species_config;
            species_config.name = "diagnostic_species";
            species_config.charge = -2.0;
            species_config.mass =
                pic::ELEMENTARY_CHARGE_SI;
            species_config.weight = 2.0;
            species_config.particles = 2;
            pic::Species2D species(species_config);
            pic::Particle2D first;
            first.position = {0.5, 0.25};
            first.velocity = {1.0, 2.0};
            first.velocity_z = -1.0;
            pic::Particle2D second;
            second.position = {0.5, 1.25};
            second.velocity = {3.0, 4.0};
            second.velocity_z = 1.0;
            species.particles() = {first, second};

            pic::ResolvedDiagnostics2DConfig config;
            config.enabled = true;
            config.max_mode = 1;
            const auto snapshot =
                pic::compute_resolved_diagnostics_2d(
                    7, 0.5, mesh, {species}, config,
                    pic::UnitSystem::SI);
            require(
                snapshot.fields.size() == 5 &&
                    snapshot.species.size() == 5,
                "resolved 2D diagnostics emitted the wrong profile shape");
            const auto& field = snapshot.fields[2];
            require_near(
                field.coordinate, 0.5, 1e-15,
                "resolved 2D profile coordinate is wrong");
            require_near(
                field.potential, 1.0, 1e-14,
                "resolved 2D transverse potential average is wrong");
            require_near(
                field.charge_density, 5.0, 1e-14,
                "resolved 2D transverse charge average is wrong");
            require_near(
                field.electric_x, 0.0, 1e-14,
                "resolved 2D transverse Ex average is wrong");
            require_near(
                field.electric_y, 0.0, 1e-14,
                "resolved 2D transverse Ey average is wrong");
            const auto& moments = snapshot.species[2];
            require_near(
                moments.number_density, 8.0, 1e-14,
                "resolved 2D number density is wrong");
            require_near(
                moments.mean_velocity_x, 2.0, 1e-14,
                "resolved 2D mean vx is wrong");
            require_near(
                moments.mean_velocity_y, 3.0, 1e-14,
                "resolved 2D mean vy is wrong");
            require_near(
                moments.mean_velocity_z, 0.0, 1e-14,
                "resolved 2D mean vz is wrong");
            require_near(
                moments.thermal_speed_x, 1.0, 1e-14,
                "resolved 2D x thermal speed is wrong");
            require_near(
                moments.thermal_speed_y, 1.0, 1e-14,
                "resolved 2D y thermal speed is wrong");
            require_near(
                moments.thermal_speed_z, 1.0, 1e-14,
                "resolved 2D z thermal speed is wrong");
            require_near(
                moments.temperature_ev, 1.0, 1e-14,
                "resolved 2D SI temperature is wrong");
            require_near(
                moments.current_density_x, -32.0, 1e-13,
                "resolved 2D x current density is wrong");
            require_near(
                moments.current_density_y, -48.0, 1e-13,
                "resolved 2D y current density is wrong");
            require_near(
                moments.current_density_z, 0.0, 1e-14,
                "resolved 2D z current density is wrong");
            const auto find_mode =
                [&](const std::string& quantity) {
                    return std::find_if(
                        snapshot.modes.begin(),
                        snapshot.modes.end(),
                        [&](const auto& mode) {
                            return mode.mode == 1 &&
                                mode.quantity == quantity &&
                                mode.species.empty();
                        });
                };
            const auto rho_mode = find_mode("charge_density");
            const auto ex_mode = find_mode("electric_x");
            const auto ey_mode = find_mode("electric_y");
            require(
                rho_mode != snapshot.modes.end() &&
                    ex_mode != snapshot.modes.end() &&
                    ey_mode != snapshot.modes.end(),
                "resolved 2D field modes are missing");
            require_near(
                rho_mode->real, 1.0, 1e-14,
                "resolved 2D cosine mode real coefficient is wrong");
            require_near(
                rho_mode->amplitude, 2.0, 1e-14,
                "resolved 2D cosine mode amplitude is wrong");
            require_near(
                ex_mode->imaginary, -1.5, 1e-14,
                "resolved 2D sine mode imaginary coefficient is wrong");
            require_near(
                ex_mode->amplitude, 3.0, 1e-14,
                "resolved 2D sine mode amplitude is wrong");
            require_near(
                ey_mode->real, -2.0, 1e-14,
                "resolved 2D signed cosine coefficient is wrong");
            require_near(
                ey_mode->amplitude, 4.0, 1e-14,
                "resolved 2D signed cosine amplitude is wrong");

            const auto output_dir = std::filesystem::path(
                "test_output_resolved_diagnostics");
            std::filesystem::remove_all(output_dir);
            pic::ResolvedDiagnostics2D diagnostics(
                output_dir, config, mesh, {species},
                pic::UnitSystem::SI);
            (void)diagnostics.sample(
                0, 0.0, mesh, {species});
            for (double& potential : mesh.phi()) {
                potential *= 3.0;
            }
            (void)diagnostics.sample(
                2, 2.0, mesh, {species});
            diagnostics.write_time_averages();
            require(
                diagnostics.sample_count() == 2,
                "resolved 2D time-average sample count is wrong");
            require(
                std::filesystem::exists(
                    output_dir /
                    "resolved_field_profiles.csv") &&
                    std::filesystem::exists(
                        output_dir /
                        "resolved_species_profiles.csv") &&
                    std::filesystem::exists(
                        output_dir / "resolved_modes.csv") &&
                    std::filesystem::exists(
                        output_dir /
                        "resolved_field_time_average.csv") &&
                    std::filesystem::exists(
                        output_dir /
                        "resolved_species_time_average.csv"),
                "resolved 2D diagnostic outputs are incomplete");
            require(
                read_file_text(
                    output_dir /
                    "resolved_field_time_average.csv")
                        .find("0,2,2,2,x,0.5,2,") !=
                    std::string::npos,
                "resolved 2D trapezoidal field average is wrong");
            std::filesystem::remove_all(output_dir);

            auto nonperiodic = config;
            require_throws_contains(
                [&] {
                    pic::Mesh2D invalid_mesh(
                        5, 4, 1.0, 2.0,
                        pic::Boundary::Dirichlet);
                    (void)pic::compute_resolved_diagnostics_2d(
                        0, 0.0, invalid_mesh, {},
                        nonperiodic,
                        pic::UnitSystem::Normalized);
                },
                "mode axis must be periodic",
                "resolved 2D modes accepted a nonperiodic axis");
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
            pic::Mesh2D mesh(5, 5, 1.0, 1.0, pic::Boundary::Dirichlet);
            pic::Species2DConfig cfg;
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
            pic::Species2D species(cfg);
            std::mt19937_64 rng(1234);
            species.initialize(mesh, rng);
            require(species.live_count() == 3,
                    "2D3V species did not initialize all particles as live");
            for (const auto& particle : species.particles()) {
                require(particle.position.x >= 0.1 && particle.position.x <= 0.2,
                        "2D3V species x initialization interval was not honored");
                require(particle.position.y >= 0.3 && particle.position.y <= 0.4,
                        "2D3V species y initialization interval was not honored");
                require_near(particle.velocity.x, 1.0, 1e-15,
                             "2D3V species vx drift initialization is wrong");
                require_near(particle.velocity.y, -2.0, 1e-15,
                             "2D3V species vy drift initialization is wrong");
                require_near(particle.velocity_z, 0.5, 1e-15,
                             "2D3V species vz drift initialization is wrong");
                require_near(particle.velocity_half.x, particle.velocity.x,
                             1e-15,
                             "2D3V species half-step vx initialization is wrong");
                require_near(particle.velocity_half.y, particle.velocity.y,
                             1e-15,
                             "2D3V species half-step vy initialization is wrong");
                require_near(particle.velocity_half_z, particle.velocity_z,
                             1e-15,
                             "2D3V species half-step vz initialization is wrong");
            }
            require_near(
                species.kinetic_energy(),
                3.0 * 0.5 * cfg.mass * cfg.weight * (1.0 + 4.0 + 0.25),
                1e-15, "2D3V species kinetic energy is wrong");
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
            require(particles2d.find("species_id,species,x,y,vx,vy,vz,alive\n") == 0, "2D particle CSV header is wrong");
            require(particles2d.find("0,ions,0.25,0.5,1.5,-0.25,0,1\n") != std::string::npos, "2D particle CSV live row is wrong");
            require(particles2d.find("0,ions,0.75,1.5,-0.5,0.25,0,0\n") != std::string::npos, "2D particle CSV dead row is wrong");
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
                    << "units = si\n"
                    << "relative_permittivity = 2.5\n"
                    << "nx = 16\n"
                    << "velocity_dimensions = 3\n"
                    << "length = 2.0\n"
                    << "dt = 0.01\n"
                    << "output_interval = 2\n"
                    << "boundary = dirichlet\n"
                    << "phi_left = -2\n"
                    << "phi_left_amplitude = 3\n"
                    << "phi_left_frequency = 4\n"
                    << "phi_left_phase = 0.25\n"
                    << "phi_right = 5\n"
                    << "phi_right_amplitude = 6\n"
                    << "phi_right_frequency = 7\n"
                    << "phi_right_phase = -0.5\n"
                    << "runtime_backend = single\n"
                    << "runtime_threads = 1\n"
                    << "[species]\n"
                    << "name = density_weighted\n"
                    << "charge = -1\n"
                    << "mass = 1\n"
                    << "density = 5\n"
                    << "particles = 10\n"
                    << "thermal_velocity = 0\n"
                    << "drift_velocity_y = 0.25\n"
                    << "drift_velocity_z = -0.5\n"
                    << "thermal_velocity_y = 0.1\n"
                    << "thermal_velocity_z = 0.2\n"
                    << "init_x_min = 0.5\n"
                    << "init_x_max = 1.5\n";
            }
            auto cfg = pic::load_config(config_path.string());
            require(cfg.species.size() == 1, "config did not load one species");
            require(
                cfg.velocity_dimensions == 3,
                "1D config did not load velocity_dimensions");
            require(cfg.runtime.backend == pic::RuntimeBackend::Serial && cfg.runtime.threads == 1,
                    "M4 1D runtime config aliases were not parsed");
            require(
                cfg.units.system == pic::UnitSystem::SI &&
                    std::abs(
                        cfg.units.permittivity() -
                        2.5 * pic::VACUUM_PERMITTIVITY_SI) <
                        1e-26,
                "SI unit-system config was not parsed");
            require_near(
                cfg.phi_left, -2.0, 1e-15,
                "1D config did not load left voltage offset");
            require_near(
                cfg.phi_left_drive.amplitude, 3.0, 1e-15,
                "1D config did not load left voltage amplitude");
            require_near(
                cfg.phi_left_drive.frequency, 4.0, 1e-15,
                "1D config did not load left voltage frequency");
            require_near(
                cfg.phi_left_drive.phase, 0.25, 1e-15,
                "1D config did not load left voltage phase");
            require_near(
                cfg.phi_right, 5.0, 1e-15,
                "1D config did not load right voltage offset");
            require_near(
                cfg.phi_right_drive.amplitude, 6.0, 1e-15,
                "1D config did not load right voltage amplitude");
            require_near(
                cfg.phi_right_drive.frequency, 7.0, 1e-15,
                "1D config did not load right voltage frequency");
            require_near(
                cfg.phi_right_drive.phase, -0.5, 1e-15,
                "1D config did not load right voltage phase");
            require_near(
                cfg.species[0].drift_velocity_y, 0.25, 1e-15,
                "1D3V config did not load drift_velocity_y");
            require_near(
                cfg.species[0].drift_velocity_z, -0.5, 1e-15,
                "1D3V config did not load drift_velocity_z");
            require_near(
                *cfg.species[0].initialization.thermal_velocity_y,
                0.1, 1e-15,
                "1D3V config did not load thermal_velocity_y");
            require_near(
                *cfg.species[0].initialization.thermal_velocity_z,
                0.2, 1e-15,
                "1D3V config did not load thermal_velocity_z");
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
                "test_invalid_units.ini",
                "units = cgs\nnx = 16\nlength = 1\ndt = 0.01\n[species]\nname = bad_units\ncharge = -1\nmass = 1\nweight = 1\nparticles = 10\n",
                [](const std::string& path) { return pic::load_config(path); },
                "invalid unit system validation did not throw");
            require_config_rejects(
                "test_invalid_relative_permittivity.ini",
                "units = si\nrelative_permittivity = 0\nnx = 16\nlength = 1\ndt = 0.01\n[species]\nname = bad_permittivity\ncharge = -1\nmass = 1\nweight = 1\nparticles = 10\n",
                [](const std::string& path) { return pic::load_config(path); },
                "invalid relative permittivity validation did not throw");
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
                "test_periodic_voltage_drive.ini",
                "boundary = periodic\nphi_right_amplitude = 1\nphi_right_frequency = 1\n",
                [](const std::string& path) { return pic::load_config(path); },
                "1D periodic config accepted a sinusoidal electrode drive");
            require_config_rejects(
                "test_invalid_velocity_dimensions.ini",
                "velocity_dimensions = 2\n",
                [](const std::string& path) { return pic::load_config(path); },
                "1D config accepted velocity_dimensions other than 1 or 3");
            require_config_rejects(
                "test_mixed_collision_schemas.ini",
                "[collisions]\nenabled = false\n"
                "[collisions.electron_mcc]\nenabled = false\n",
                [](const std::string& path) { return pic::load_config(path); },
                "1D config mixed legacy and named collision schemas");
            require_config_rejects(
                "test_inactive_transverse_drift.ini",
                "[species]\nweight = 1\ndrift_velocity_y = 1\n",
                [](const std::string& path) { return pic::load_config(path); },
                "1D1V config accepted active transverse drift");
            require_config_rejects(
                "test_inactive_transverse_thermal.ini",
                "[species]\nweight = 1\nthermal_velocity_z = 1\n",
                [](const std::string& path) { return pic::load_config(path); },
                "1D1V config accepted active transverse thermal velocity");
            require_config_rejects(
                "test_invalid_voltage_frequency.ini",
                "boundary = dirichlet\nphi_right_amplitude = 1\nphi_right_frequency = -1\n",
                [](const std::string& path) { return pic::load_config(path); },
                "1D config accepted a negative electrode frequency");
            require_config_rejects(
                "test_zero_voltage_frequency.ini",
                "boundary = dirichlet\nphi_right_amplitude = 1\nphi_right_frequency = 0\n",
                [](const std::string& path) { return pic::load_config(path); },
                "1D config accepted a driven electrode without positive frequency");
            require_config_rejects(
                "test_steady_voltage_drive.ini",
                "boundary = dirichlet\nmode = steady_state\nphi_right_amplitude = 1\nphi_right_frequency = 1\n",
                [](const std::string& path) { return pic::load_config(path); },
                "1D config accepted an RF drive with instantaneous steady-state convergence");
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
                    << "magnetic_field_x = 0.25\n"
                    << "magnetic_field_y = -0.5\n"
                    << "magnetic_field_z = 1.75\n"
                    << "[species.electrons]\n"
                    << "charge = -1\n"
                    << "mass = 1\n"
                    << "density = 4\n"
                    << "particles = 12\n"
                    << "drift_velocity_z = 0.75\n"
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
            require_near(cfg2.magnetic_field_x, 0.25, 1e-15,
                         "2D3V config did not load magnetic_field_x");
            require_near(cfg2.magnetic_field_y, -0.5, 1e-15,
                         "2D3V config did not load magnetic_field_y");
            require_near(cfg2.magnetic_field_z, 1.75, 1e-15,
                         "2D3V config did not load magnetic_field_z");
            require(cfg2.species.size() == 1, "2D config did not load one species");
            require(std::abs(cfg2.species[0].weight - (4.0 / 12.0)) < 1e-15, "2D density-derived macro-particle weight is wrong");
            require_near(cfg2.species[0].drift_velocity_z, 0.75, 1e-15,
                         "2D3V config did not load drift_velocity_z");
            std::filesystem::remove(config_2d_path);
            require_throws([] {
                const auto path = std::filesystem::path("test_invalid_2d_electrode.ini");
                { std::ofstream out(path); out << "dimension = 2\nphi_left = inf\n"; }
                try { (void)pic::load_config_2d(path.string()); } catch (...) { std::filesystem::remove(path); throw; }
                std::filesystem::remove(path);
            }, "invalid 2D electrode potential validation did not throw");
            require_throws([] {
                const auto path = std::filesystem::path(
                    "test_periodic_axis_potential.ini");
                {
                    std::ofstream out(path);
                    out << "dimension = 2\n"
                        << "boundary = dirichlet\n"
                        << "boundary_x = periodic\n"
                        << "phi_left = 1\n";
                }
                try {
                    (void)pic::load_config_2d(path.string());
                } catch (...) {
                    std::filesystem::remove(path);
                    throw;
                }
                std::filesystem::remove(path);
            }, "2D config accepted an electrode potential on a periodic axis");
            require_throws([] {
                const auto path = std::filesystem::path(
                    "test_invalid_2d_axis_boundary.ini");
                {
                    std::ofstream out(path);
                    out << "dimension = 2\n"
                        << "boundary_y = open\n";
                }
                try {
                    (void)pic::load_config_2d(path.string());
                } catch (...) {
                    std::filesystem::remove(path);
                    throw;
                }
                std::filesystem::remove(path);
            }, "2D config accepted an invalid per-axis boundary value");
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
                { std::ofstream out(path); out << "dimension = 2\nmagnetic_field_x = nan\n"; }
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
            const auto output_dir =
                std::filesystem::path(
                    "test_output_1d_boundary_losses");
            const auto checkpoint_path =
                output_dir / "checkpoint_1.apc";
            std::filesystem::remove_all(output_dir);

            pic::Config cfg;
            cfg.nx = 3;
            cfg.length = 1.0;
            cfg.dt = 0.2;
            cfg.steps = 1;
            cfg.output_interval = 1;
            cfg.output_dir = output_dir.string();
            cfg.boundary = pic::Boundary::Dirichlet;
            cfg.checkpoint_output = true;
            cfg.checkpoint_interval = 1;
            cfg.species = {
                pic::SpeciesConfig{
                    "left", 0.0, 2.0, 3.0, 1, 1.0,
                    -10.0, 0.0, 0.0, -1.0},
                pic::SpeciesConfig{
                    "right", 0.0, 2.0, 3.0, 1, 1.0,
                    10.0, 0.0, 0.0, -1.0}};

            pic::Simulation simulation(cfg);
            const auto summary = simulation.run();
            const auto& losses =
                simulation.species_boundary_losses();
            require(
                summary.final_sample.live_particles == 0 &&
                    losses.size() == 2 &&
                    losses[0].absorbed_left == 1 &&
                    losses[0].absorbed_right == 0 &&
                    losses[1].absorbed_left == 0 &&
                    losses[1].absorbed_right == 1,
                "1D species-resolved wall counts are wrong");
            require_near(
                losses[0].kinetic_energy_left, 300.0,
                1e-12,
                "1D left-wall represented impact energy is wrong");
            require_near(
                losses[1].kinetic_energy_right, 300.0,
                1e-12,
                "1D right-wall represented impact energy is wrong");
            const auto& power =
                simulation.species_power_transfer();
            require_near(
                power[0].electric_work, 0.0, 1e-14,
                "neutral left-wall particle acquired electric work");
            require_near(
                power[1].electric_work, 0.0, 1e-14,
                "neutral right-wall particle acquired electric work");

            const auto diagnostics = read_file_text(
                output_dir / "boundary_losses.csv");
            require(
                diagnostics.find(
                    "absorbed_left_count_left,"
                    "absorbed_right_count_left,"
                    "absorbed_left_charge_left_normalized,"
                    "absorbed_right_charge_left_normalized,"
                    "absorbed_left_kinetic_energy_left_normalized,"
                    "absorbed_right_kinetic_energy_left_normalized")
                    != std::string::npos &&
                    count_lines(diagnostics) == 3,
                "1D wall-loss CSV contract is wrong");
            const auto power_diagnostics = read_file_text(
                output_dir / "power_transfer.csv");
            require(
                power_diagnostics.find(
                    "step,time,counter_origin_step,"
                    "electric_work_left_normalized,"
                    "electric_work_right_normalized\n") == 0 &&
                    count_lines(power_diagnostics) == 3,
                "1D power-transfer CSV contract is wrong");

            pic::Simulation restarted(cfg);
            restarted.load_checkpoint(checkpoint_path);
            const auto& restarted_losses =
                restarted.species_boundary_losses();
            require(
                restarted_losses[0].absorbed_left == 1 &&
                    restarted_losses[1].absorbed_right == 1 &&
                    restarted.boundary_loss_origin_step() == 0,
                "1D checkpoint lost species-resolved wall counts");
            require_near(
                restarted_losses[0].kinetic_energy_left,
                300.0, 1e-12,
                "1D checkpoint lost left-wall impact energy");
            require_near(
                restarted_losses[1].kinetic_energy_right,
                300.0, 1e-12,
                "1D checkpoint lost right-wall impact energy");
            require(
                read_file_text(checkpoint_path).find(
                    "AuroraPIC-checkpoint-v7\n") == 0,
                "1D power-transfer checkpoint did not use v7");

            const auto legacy_v6_path =
                output_dir / "legacy_v6.apc";
            {
                std::istringstream input(
                    read_file_text(checkpoint_path));
                std::ofstream legacy(legacy_v6_path);
                std::string line;
                bool first = true;
                while (std::getline(input, line)) {
                    if (first) {
                        legacy
                            << "AuroraPIC-checkpoint-v6\n";
                        first = false;
                    } else if (
                        line.starts_with(
                            "power_transfer")) {
                        continue;
                    } else {
                        legacy << line << '\n';
                    }
                }
            }
            pic::Simulation legacy_v6_restarted(cfg);
            legacy_v6_restarted.load_checkpoint(
                legacy_v6_path);
            require(
                legacy_v6_restarted
                        .boundary_loss_origin_step() == 0 &&
                    legacy_v6_restarted
                            .species_boundary_losses()[0]
                            .absorbed_left == 1 &&
                    legacy_v6_restarted
                            .power_transfer_origin_step() == 1,
                "legacy v6 restart did not preserve wall data or "
                "expose its power-counter origin");

            const auto legacy_path =
                output_dir / "legacy_v5.apc";
            {
                std::istringstream input(
                    read_file_text(checkpoint_path));
                std::ofstream legacy(legacy_path);
                std::string line;
                bool first = true;
                while (std::getline(input, line)) {
                    if (first) {
                        legacy
                            << "AuroraPIC-checkpoint-v5\n";
                        first = false;
                    } else if (
                        line.starts_with(
                            "boundary_loss") ||
                        line.starts_with(
                            "power_transfer")) {
                        continue;
                    } else {
                        legacy << line << '\n';
                    }
                }
            }
            pic::Simulation legacy_restarted(cfg);
            legacy_restarted.load_checkpoint(legacy_path);
            require(
                legacy_restarted.boundary_loss_origin_step() == 1 &&
                    legacy_restarted
                            .species_boundary_losses()[0]
                            .absorbed_left == 0 &&
                    legacy_restarted
                            .species_boundary_losses()[1]
                            .absorbed_right == 0,
                "legacy v5 restart did not expose its wall-counter "
                "origin");
            std::filesystem::remove_all(output_dir);
        }
        {
            const auto output_dir =
                std::filesystem::path(
                    "test_output_1d_power_transfer");
            const auto checkpoint_path =
                output_dir / "power.apc";
            std::filesystem::remove_all(output_dir);

            pic::Config cfg;
            cfg.nx = 17;
            cfg.length = 1.0;
            cfg.dt = 1e-4;
            cfg.steps = 1;
            cfg.output_interval = 1;
            cfg.output_dir = output_dir.string();
            cfg.boundary = pic::Boundary::Dirichlet;
            cfg.species = {
                pic::SpeciesConfig{
                    "electrons", -1.0, 1.0, 1.0, 32, 1.0,
                    0.0, 0.0, 0.1, 0.9}};

            pic::Simulation continuous(cfg);
            continuous.initialize();
            const double initial_energy =
                continuous.sample().kinetic_energy;
            continuous.step();
            const double energy_change =
                continuous.sample().kinetic_energy -
                initial_energy;
            const double electric_work =
                continuous.species_power_transfer()[0]
                    .electric_work;
            require(
                std::abs(electric_work) > 1e-12,
                "1D electric-work diagnostic exercised a zero-work "
                "case");
            require_near(
                electric_work, energy_change, 1e-12,
                "collisionless survivor electric work does not "
                "equal its kinetic-energy change");
            continuous.save_checkpoint(checkpoint_path);

            pic::Simulation restarted(cfg);
            restarted.load_checkpoint(checkpoint_path);
            require(
                restarted.power_transfer_origin_step() == 0,
                "1D checkpoint lost the power-transfer origin");
            require_near(
                restarted.species_power_transfer()[0]
                    .electric_work,
                electric_work, 1e-14,
                "1D checkpoint lost cumulative electric work");
            continuous.step();
            restarted.step();
            require_near(
                restarted.species_power_transfer()[0]
                    .electric_work,
                continuous.species_power_transfer()[0]
                    .electric_work,
                1e-14,
                "1D power transfer is not restart deterministic");
            std::filesystem::remove_all(output_dir);
        }
        {
            const auto output_dir =
                std::filesystem::path("test_output_1d_rf_drive");
            const auto checkpoint_path =
                output_dir / "rf_checkpoint.apc";
            std::filesystem::remove_all(output_dir);

            pic::Config cfg;
            cfg.nx = 9;
            cfg.length = 1.0;
            cfg.dt = 0.125;
            cfg.steps = 2;
            cfg.output_interval = 1;
            cfg.output_dir = output_dir.string();
            cfg.boundary = pic::Boundary::Dirichlet;
            cfg.phi_left = -0.25;
            cfg.phi_left_drive = {
                0.5, 1.0, std::numbers::pi / 2.0};
            cfg.phi_right = 0.5;
            cfg.phi_right_drive = {2.0, 1.0, 0.0};
            cfg.species = {
                pic::SpeciesConfig{
                    "neutral", 0.0, 1.0, 1.0, 4, 1.0,
                    0.0, 0.0, 0.0, -1.0}};

            pic::Simulation sim(cfg);
            sim.initialize();
            require_near(
                sim.grid().phi().front(), 0.25, 1e-14,
                "1D RF drive applied the wrong initial left potential");
            require_near(
                sim.grid().phi().back(), 0.5, 1e-14,
                "1D RF drive was not zero-phase at time zero");
            sim.step();
            const double expected_left =
                -0.25 + 0.5 * std::sin(3.0 * std::numbers::pi / 4.0);
            const double expected_right =
                0.5 + 2.0 * std::sin(std::numbers::pi / 4.0);
            require_near(
                sim.grid().phi().front(), expected_left, 1e-14,
                "1D RF drive used the wrong new-time left potential");
            require_near(
                sim.grid().phi().back(), expected_right, 1e-14,
                "1D RF drive used the wrong new-time right potential");
            sim.save_checkpoint(checkpoint_path);

            pic::Simulation restarted(cfg);
            restarted.load_checkpoint(checkpoint_path);
            require_near(
                restarted.grid().phi().front(), expected_left, 1e-14,
                "1D RF restart did not restore the left drive phase");
            require_near(
                restarted.grid().phi().back(), expected_right, 1e-14,
                "1D RF restart did not restore the right drive phase");

            std::filesystem::remove_all(output_dir);
            pic::Simulation output_sim(cfg);
            output_sim.run();
            const auto scalars =
                read_file_text(output_dir / "scalars.csv");
            require(
                scalars.find(
                    "step,time,kinetic_energy,field_energy,"
                    "total_energy,charge_l1,live_particles,"
                    "phi_left,phi_right,"
                    "live_particles_neutral\n") == 0,
                "1D RF scalar diagnostics omitted electrode potentials");
        }
        {
            const auto output_dir =
                std::filesystem::path(
                    "test_output_1d3v_acceptance");
            std::filesystem::remove_all(output_dir);
            pic::Config cfg;
            cfg.velocity_dimensions = 3;
            cfg.steps = 0;
            cfg.output_dir = output_dir.string();
            cfg.initialization_acceptance
                .max_relative_current_imbalance = 0.0;
            cfg.species = {
                pic::SpeciesConfig{
                    "electrons", -1.0, 1.0, 1.0, 4, 1.0,
                    0.0, 0.0, 0.0, -1.0, {}, 1.0, 0.0},
                pic::SpeciesConfig{
                    "ions", 1.0, 1.0, 1.0, 4, 1.0,
                    0.0, 0.0, 0.0, -1.0, {}, 0.0, 0.0}};
            require_throws(
                [&] {
                    pic::Simulation simulation(cfg);
                    (void)simulation.run();
                },
                "1D3V initialization acceptance ignored transverse current");
            std::filesystem::remove_all(output_dir);
        }
        {
            const auto table_path =
                std::filesystem::path(
                    "test_mcc_1d3v_elastic.dat");
            const auto output_dir =
                std::filesystem::path(
                    "test_output_mcc_1d3v");
            const auto checkpoint_path =
                output_dir / "manual.apc";
            std::filesystem::remove_all(output_dir);
            {
                std::ofstream table(table_path);
                table << "0 0.2\n10 0.2\n";
            }

            pic::Config cfg;
            cfg.velocity_dimensions = 3;
            cfg.nx = 16;
            cfg.length = 1.0;
            cfg.dt = 0.1;
            cfg.steps = 8;
            cfg.output_interval = 8;
            cfg.output_dir = output_dir.string();
            cfg.seed = 1307;
            cfg.collisions.enabled = true;
            cfg.collisions.model =
                pic::CollisionModelKind::NullCollision;
            cfg.collisions.species = "electrons";
            cfg.collisions.neutral_density = 2.0;
            cfg.collisions.max_frequency = 1.0;
            cfg.collisions.channels = {
                pic::CollisionChannelConfig{
                    "elastic",
                    pic::CollisionProcessKind::Elastic,
                    table_path, 0.0, 1.0, 1.0}};
            cfg.species = {
                pic::SpeciesConfig{
                    "electrons", 0.0, 1.0, 1.0, 128, 1.0,
                    1.0, 0.0, 0.0, -1.0}};

            pic::Simulation continuous(cfg);
            continuous.initialize();
            const double initial_energy =
                continuous.sample().kinetic_energy;
            for (std::size_t step = 0; step < 4; ++step) {
                continuous.step();
            }
            continuous.save_checkpoint(checkpoint_path);
            for (std::size_t step = 4; step < cfg.steps; ++step) {
                continuous.step();
            }
            const bool populated_transverse =
                std::any_of(
                    continuous.species().front().particles().begin(),
                    continuous.species().front().particles().end(),
                    [](const pic::Particle& particle) {
                        return std::abs(particle.velocity_y) > 1e-12 ||
                               std::abs(particle.velocity_z) > 1e-12;
                    });
            require(
                populated_transverse,
                "1D3V elastic MCC did not scatter into transverse velocity");
            require_near(
                continuous.sample().kinetic_energy,
                initial_energy, 1e-10,
                "1D3V isotropic elastic MCC changed particle energy");

            pic::Simulation restarted(cfg);
            restarted.load_checkpoint(checkpoint_path);
            for (std::size_t step = 4; step < cfg.steps; ++step) {
                restarted.step();
            }
            require_species_close(
                continuous.species(), restarted.species(),
                "1D3V MCC checkpoint restart determinism");
            require(
                continuous.collision_diagnostics().candidates ==
                    restarted.collision_diagnostics().candidates &&
                continuous.collision_diagnostics().channel_collisions ==
                    restarted.collision_diagnostics().channel_collisions,
                "1D3V MCC restart lost collision diagnostics");
            require(
                read_file_text(checkpoint_path).find(
                    "AuroraPIC-checkpoint-v7\n") == 0,
                "1D3V checkpoint did not use the power-transfer-aware "
                "spatial-average velocity format");
            std::filesystem::remove_all(output_dir);
            std::filesystem::remove(table_path);
        }
        {
            const auto parsed =
                pic::load_config(
                    (std::filesystem::path(
                         AURORA_TEST_SOURCE_DIR) /
                     "examples" /
                     "mcc_ionization_1d.cfg")
                        .string());
            require(
                parsed.velocity_dimensions == 3 &&
                    parsed.max_particles_per_species == 256 &&
                    parsed.collision_models.size() == 2 &&
                    parsed.collision_models[0].name ==
                        "electron_mcc" &&
                    parsed.collision_models[0].config.channels.size() ==
                        1 &&
                    parsed.collision_models[0].config.channels[0]
                            .process ==
                        pic::CollisionProcessKind::Ionization &&
                    parsed.collision_models[0].config.channels[0]
                            .secondary_species ==
                        "electrons" &&
                    parsed.collision_models[0].config.channels[0]
                            .ion_species ==
                        "ions" &&
                    parsed.collision_models[1].name == "ion_mcc",
                "named 1D MCC configuration did not preserve model "
                "and ionization product mappings");

            const auto turner_config_path =
                std::filesystem::path(
                    "test_turner_ion_collision.cfg");
            const auto turner_table_path =
                std::filesystem::path("ion_backward.dat");
            const auto turner_output_dir =
                std::filesystem::path(
                    "test_output_turner_collision");
            std::filesystem::remove_all(turner_output_dir);
            {
                std::ofstream table(turner_table_path);
                table << "0 1e-19\n10 1e-19\n";
            }
            {
                std::ofstream config(turner_config_path);
                config
                    << "config_version = 1\n"
                    << "units = si\n"
                    << "dimension = 1\n"
                    << "velocity_dimensions = 3\n"
                    << "nx = 3\n"
                    << "length = 0.067\n"
                    << "dt = 1e-12\n"
                    << "steps = 1\n"
                    << "boundary = dirichlet\n"
                    << "output_interval = 1\n"
                    << "output_dir = "
                    << turner_output_dir.string() << "\n"
                    << "[collisions.ion_mcc]\n"
                    << "model = null_collision\n"
                    << "species = ions\n"
                    << "neutral_density = 9.64e20\n"
                    << "neutral_mass = 6.67e-27\n"
                    << "neutral_temperature = 300\n"
                    << "max_frequency = 1e8\n"
                    << "[collisions.ion_mcc.channel.backward]\n"
                    << "type = elastic\n"
                    << "cross_section_file = ion_backward.dat\n"
                    << "energy_scale = 1.602176634e-19\n"
                    << "angular_model = backward\n"
                    << "energy_frame = center_of_mass\n"
                    << "[species.ions]\n"
                    << "charge = 1.602176634e-19\n"
                    << "mass = 6.67e-27\n"
                    << "density = 1e14\n"
                    << "particles = 2\n";
            }
            const auto turner_config =
                pic::load_config(turner_config_path.string());
            require(
                turner_config.collision_models.size() == 1 &&
                    turner_config.collision_models[0].config
                            .neutral_mass == 6.67e-27 &&
                    turner_config.collision_models[0].config
                            .neutral_temperature == 300.0 &&
                    turner_config.collision_models[0].config
                            .gas_data_units ==
                        pic::UnitSystem::SI &&
                    turner_config.collision_models[0].config.channels[0]
                            .angular_scattering ==
                        pic::AngularScatteringKind::Backward &&
                    turner_config.collision_models[0].config.channels[0]
                            .energy_frame ==
                        pic::CollisionEnergyFrame::CenterOfMass,
                "Turner ion-collision frame and scattering contract "
                "did not survive 1D config parsing");
            pic::Simulation turner_collision_smoke(turner_config);
            turner_collision_smoke.initialize();
            turner_collision_smoke.step();
            std::filesystem::remove(turner_config_path);
            std::filesystem::remove(turner_table_path);
            std::filesystem::remove_all(turner_output_dir);

            const auto table_path =
                std::filesystem::path(
                    "test_mcc_1d_multi_reactive.dat");
            const auto output_dir =
                std::filesystem::path(
                    "test_output_mcc_1d_multi");
            const auto checkpoint_path =
                output_dir / "manual.apc";
            std::filesystem::remove_all(output_dir);
            {
                std::ofstream table(table_path);
                table << "0 0.5\n1000 0.5\n";
            }
            pic::SpeciesConfig electrons;
            electrons.name = "electrons";
            electrons.charge = -1.0;
            electrons.mass = 1.0;
            electrons.weight = 1.0;
            electrons.particles = 64;
            electrons.drift_velocity = 4.0;
            electrons.thermal_velocity = 0.0;
            electrons.initialization.loading =
                pic::ParticleLoading::QuietStart;
            pic::SpeciesConfig ions = electrons;
            ions.name = "ions";
            ions.charge = 1.0;
            ions.mass = 40.0;

            pic::CollisionConfig electron_mcc;
            electron_mcc.enabled = true;
            electron_mcc.model =
                pic::CollisionModelKind::NullCollision;
            electron_mcc.species = "electrons";
            electron_mcc.neutral_density = 1.0;
            electron_mcc.max_frequency = 3.0;
            electron_mcc.channels = {
                pic::CollisionChannelConfig{
                    "ionization",
                    pic::CollisionProcessKind::Ionization,
                    table_path, 1.0, 1.0, 1.0,
                    "electrons", "ions"}};
            pic::CollisionConfig ion_mcc;
            ion_mcc.enabled = true;
            ion_mcc.model =
                pic::CollisionModelKind::NullCollision;
            ion_mcc.species = "ions";
            ion_mcc.neutral_density = 1.0;
            ion_mcc.max_frequency = 3.0;
            ion_mcc.channels = {
                pic::CollisionChannelConfig{
                    "elastic",
                    pic::CollisionProcessKind::Elastic,
                    table_path, 0.0, 1.0, 1.0}};

            pic::Config cfg;
            cfg.velocity_dimensions = 3;
            cfg.nx = 16;
            cfg.dt = 0.1;
            cfg.steps = 4;
            cfg.output_interval = 4;
            cfg.output_dir = output_dir.string();
            cfg.seed = 271828;
            cfg.max_particles_per_species = 256;
            cfg.species = {electrons, ions};
            cfg.collision_models = {
                {"electron_mcc", electron_mcc},
                {"ion_mcc", ion_mcc}};

            auto duplicate_target = cfg;
            duplicate_target.collision_models[1].config.species =
                "electrons";
            require_throws(
                [&] {
                    pic::Simulation invalid(duplicate_target);
                },
                "multiple named MCC models accepted the same target");
            auto invalid_products = cfg;
            invalid_products.species[1].weight = 2.0;
            require_throws(
                [&] {
                    pic::Simulation invalid(invalid_products);
                },
                "1D ionization accepted unequal product macro weights");
            auto inactive_velocity = cfg;
            inactive_velocity.velocity_dimensions = 1;
            require_throws(
                [&] {
                    pic::Simulation invalid(inactive_velocity);
                },
                "1D1V accepted an ionization collision channel");

            pic::Simulation continuous(cfg);
            continuous.initialize();
            const double initial_kinetic =
                continuous.sample().kinetic_energy;
            continuous.step();
            const auto first_collision =
                continuous.collision_diagnostics();
            require(
                first_collision.channel_names ==
                    std::vector<std::string>{
                        "electron_mcc.ionization",
                        "ion_mcc.elastic"} &&
                    first_collision.channel_collisions[0] > 0 &&
                    first_collision.channel_collisions[1] > 0,
                "simultaneous 1D electron and ion MCC did not "
                "exercise both targets");
            const std::size_t ionizations =
                static_cast<std::size_t>(
                    first_collision.channel_collisions[0]);
            require(
                continuous.species()[0].live_count() ==
                        64 + ionizations &&
                    continuous.species()[1].live_count() ==
                        64 + ionizations,
                "1D ionization did not create paired products");
            require_near(
                continuous.sample().kinetic_energy,
                initial_kinetic -
                    static_cast<double>(ionizations),
                1e-10,
                "1D ionization did not remove exactly one threshold "
                "energy per macro-event");

            continuous.step();
            continuous.save_checkpoint(checkpoint_path);
            for (std::size_t step = 2; step < cfg.steps; ++step) {
                continuous.step();
            }
            pic::Simulation restarted(cfg);
            restarted.load_checkpoint(checkpoint_path);
            for (std::size_t step = 2; step < cfg.steps; ++step) {
                restarted.step();
            }
            require_species_close(
                continuous.species(), restarted.species(),
                "multi-model 1D ionization checkpoint restart");
            require(
                continuous.collision_diagnostics().channel_collisions ==
                    restarted.collision_diagnostics()
                        .channel_collisions,
                "multi-model 1D restart lost collision diagnostics");

            auto changed = cfg;
            changed.collision_models[0].config.neutral_density =
                2.0;
            require_throws(
                [&] {
                    pic::Simulation incompatible(changed);
                    incompatible.load_checkpoint(checkpoint_path);
                },
                "multi-model checkpoint accepted changed MCC physics");

            auto bounded = cfg;
            bounded.max_particles_per_species = 64;
            pic::Simulation bounded_simulation(bounded);
            bounded_simulation.initialize();
            require_throws_contains(
                [&] { bounded_simulation.step(); },
                "max_particles_per_species",
                "1D ionization ignored product storage capacity");
            require(
                bounded_simulation.species()[0].particles().size() ==
                        64 &&
                    bounded_simulation.species()[1].particles().size() ==
                        64,
                "1D ionization capacity failure created partial "
                "products");

            std::filesystem::remove_all(output_dir);
            std::filesystem::remove(table_path);
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
            const auto boundary_flux = read_file_text(
                cfg.output_dir / "boundary_flux.csv");
            require(
                boundary_flux.find(
                    "step,time,window_start_step,window_start_time,"
                    "window_duration,species,boundary,"
                    "absorbed_macroparticles,"
                    "cumulative_absorbed_macroparticles,"
                    "represented_particles,represented_charge,"
                    "represented_particle_rate,charge_rate\n") == 0,
                "2D boundary-flux diagnostics header is wrong");
            require(
                count_lines(boundary_flux) == 25,
                "2D boundary-flux diagnostics wrote unexpected rows");
            require(std::filesystem::exists(cfg.output_dir / "particles_0.csv"), "2D simulation did not write initial particle sample");
            require(std::filesystem::exists(cfg.output_dir / "particles_2.csv"), "2D simulation did not write interval particle sample");
            require(std::filesystem::exists(cfg.output_dir / "particles_3.csv"), "2D simulation did not write final particle sample");
            require(!std::filesystem::exists(cfg.output_dir / "particles_1.csv"), "2D simulation wrote an unexpected particle sample interval");
            const auto particles = read_file_text(cfg.output_dir / "particles_0.csv");
            require(particles.find("species_id,species,x,y,vx,vy,vz,alive\n") == 0, "2D particle diagnostics header is wrong");
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
            const auto output_dir = std::filesystem::path(
                "test_output_2d_current_regulation");
            const auto checkpoint_path =
                output_dir / "regulated.apc";
            std::filesystem::remove_all(output_dir);

            pic::Species2DConfig electrons;
            electrons.name = "electrons";
            electrons.charge = -1e-12;
            electrons.mass = 1.0;
            electrons.weight = 2.0;
            electrons.particles = 3;
            electrons.drift_velocity_x = -10.0;
            electrons.thermal_velocity = 0.0;
            electrons.init_x_min = 0.49;
            electrons.init_x_max = 0.51;
            electrons.init_y_min = 0.45;
            electrons.init_y_max = 0.55;
            electrons.initialization.loading =
                pic::ParticleLoading::QuietStart;
            auto ions = electrons;
            ions.name = "ions";
            ions.charge = 1e-12;
            ions.weight = 4.0;
            ions.particles = 1;

            pic::Simulation2DConfig cfg;
            cfg.nx = 6;
            cfg.ny = 6;
            cfg.length_x = 1.0;
            cfg.length_y = 1.0;
            cfg.dt = 0.1;
            cfg.steps = 2;
            cfg.boundary = pic::Boundary::Dirichlet;
            cfg.boundary_config.left.potential = 5.0;
            cfg.output_interval = 1;
            cfg.output_dir = output_dir;
            cfg.seed = 777;
            cfg.max_particles_per_species = 5;
            cfg.species = {electrons, ions};
            cfg.current_regulated_source =
                pic::CurrentRegulatedSource2DConfig{
                    "electrons",
                    pic::BoundarySide2DName::Left,
                    pic::BoundarySide2DName::Right,
                    0.1, {}, 0.0};
            cfg.potential_reference =
                pic::PotentialReference2DConfig{
                    pic::CoordinateAxis::X, 0.8, 2.0};

            pic::Simulation2D continuous(cfg);
            continuous.initialize();
            double referenced_mean = 0.0;
            double reference_weight = 0.0;
            for (std::size_t j = 0; j < cfg.ny; ++j) {
                const double weight =
                    (j == 0 || j + 1 == cfg.ny) ? 0.5 : 1.0;
                referenced_mean += weight *
                    continuous.mesh().phi()[
                        continuous.mesh().index(4, j)];
                reference_weight += weight;
            }
            require_near(
                referenced_mean / reference_weight,
                2.0, 1e-12,
                "2D line-average potential reference missed its target");
            require(
                std::isfinite(
                    continuous.potential_reference_offset()),
                "2D potential-reference offset is not finite");
            continuous.step();
            require(
                continuous.species_boundary_losses()[0]
                        .absorbed_left == 3 &&
                    continuous.species_boundary_losses()[1]
                        .absorbed_left == 1,
                "2D species-resolved boundary loss accounting is wrong");
            continuous.save_checkpoint(checkpoint_path);
            require(
                read_file_text(checkpoint_path).find(
                    "AuroraPIC-checkpoint-v10\n") == 0,
                "2D current regulation checkpoint did not use v10");
            const auto legacy_checkpoint_path =
                output_dir / "regulated-v8.apc";
            const auto v9_checkpoint_path =
                output_dir / "regulated-v9.apc";
            const auto write_legacy_checkpoint = [&](
                const std::filesystem::path& legacy_path,
                const std::string& magic,
                std::size_t removed_fields
            ) {
                std::istringstream input(
                    read_file_text(checkpoint_path));
                std::ofstream legacy(legacy_path);
                require(
                    static_cast<bool>(legacy),
                    "cannot create synthetic legacy checkpoint");
                std::string line;
                while (std::getline(input, line)) {
                    if (line == "AuroraPIC-checkpoint-v10") {
                        legacy << magic << '\n';
                        continue;
                    }
                    if (line.rfind(
                            "current_regulated_source 1 ", 0) == 0) {
                        std::istringstream fields(line);
                        std::vector<std::string> tokens{
                            std::istream_iterator<std::string>(fields),
                            std::istream_iterator<std::string>()};
                        require(
                            tokens.size() > removed_fields,
                            "synthetic legacy controller line is incomplete");
                        tokens.resize(
                            tokens.size() - removed_fields);
                        for (std::size_t index = 0;
                             index < tokens.size(); ++index) {
                            if (index != 0) legacy << ' ';
                            legacy << tokens[index];
                        }
                        legacy << '\n';
                        continue;
                    }
                    legacy << line << '\n';
                }
            };
            write_legacy_checkpoint(
                legacy_checkpoint_path,
                "AuroraPIC-checkpoint-v8", 15);
            write_legacy_checkpoint(
                v9_checkpoint_path,
                "AuroraPIC-checkpoint-v9", 8);
            continuous.step();
            const auto& regulated =
                *continuous
                     .current_regulated_source_diagnostics();
            require(
                regulated.macro_particles_created == 1 &&
                    continuous.species()[0].live_count() == 1 &&
                    regulated.control_updates == 2 &&
                    regulated.reverse_diagnostics_start_step == 0,
                "2D current regulator did not emit the charge-balanced macro-particle");
            require_near(
                regulated.represented_particles_created,
                2.0, 1e-14,
                "2D current regulator represented-particle accounting is wrong");
            require_near(
                regulated.processed_monitored_charge,
                -2e-12, 1e-24,
                "2D current regulator monitored the wrong represented charge");

            pic::Simulation2D restarted(cfg);
            restarted.load_checkpoint(checkpoint_path);
            restarted.step();
            require_species_close(
                continuous.species(), restarted.species(),
                "2D current-regulation checkpoint determinism");
            require(
                restarted.species_boundary_losses()[0]
                        .absorbed_left == 3 &&
                    restarted
                        .current_regulated_source_diagnostics()
                        ->macro_particles_created == 1 &&
                    restarted
                        .current_regulated_source_diagnostics()
                        ->control_updates ==
                    regulated.control_updates &&
                    restarted
                        .current_regulated_source_diagnostics()
                        ->reverse_demand_steps ==
                    regulated.reverse_demand_steps &&
                    restarted
                        .current_regulated_source_diagnostics()
                        ->reverse_one_macro_steps ==
                    regulated.reverse_one_macro_steps &&
                    restarted
                        .current_regulated_source_diagnostics()
                        ->squared_reverse_demand_macroparticles ==
                    regulated.squared_reverse_demand_macroparticles,
                "2D current-regulation checkpoint lost controller state");
            pic::Simulation2D v9_restarted(cfg);
            v9_restarted.load_checkpoint(v9_checkpoint_path);
            require(
                v9_restarted
                        .current_regulated_source_diagnostics()
                        ->control_updates == 1 &&
                    v9_restarted
                        .current_regulated_source_diagnostics()
                        ->reverse_distribution_start_step == 1,
                "2D v9 restart did not mark partial reverse distribution");
            v9_restarted.step();
            require_species_close(
                continuous.species(),
                v9_restarted.species(),
                "2D v9 current-regulation checkpoint compatibility");
            pic::Simulation2D legacy_restarted(cfg);
            legacy_restarted.load_checkpoint(
                legacy_checkpoint_path);
            require(
                legacy_restarted
                        .current_regulated_source_diagnostics()
                        ->reverse_diagnostics_start_step == 1 &&
                    legacy_restarted
                        .current_regulated_source_diagnostics()
                        ->control_updates == 0,
                "2D v8 restart did not mark partial reverse diagnostics");
            legacy_restarted.step();
            require_species_close(
                continuous.species(),
                legacy_restarted.species(),
                "2D v8 current-regulation checkpoint compatibility");
            require(
                legacy_restarted
                        .current_regulated_source_diagnostics()
                        ->control_updates == 1,
                "2D v8 restart did not begin reverse-demand accounting");
            auto changed_controller = cfg;
            changed_controller.current_regulated_source
                ->emission_inset = 0.2;
            require_throws_contains(
                [&] {
                    pic::Simulation2D incompatible(
                        changed_controller);
                    incompatible.load_checkpoint(
                        checkpoint_path);
                },
                "current-regulated source metadata",
                "2D checkpoint accepted changed current-controller physics");
            changed_controller = cfg;
            changed_controller.current_regulated_source
                ->control_mode =
                    pic::CurrentSourceControlMode::TimestepLocal;
            require_throws_contains(
                [&] {
                    pic::Simulation2D incompatible(
                        changed_controller);
                    incompatible.load_checkpoint(
                        checkpoint_path);
                },
                "current-regulated source metadata",
                "2D checkpoint accepted changed current-controller mode");
            auto changed_reference = cfg;
            changed_reference.potential_reference->target = 1.0;
            require_throws_contains(
                [&] {
                    pic::Simulation2D incompatible(
                        changed_reference);
                    incompatible.load_checkpoint(
                        checkpoint_path);
                },
                "potential-reference metadata",
                "2D checkpoint accepted changed potential-reference physics");

            auto reverse_electrons = electrons;
            reverse_electrons.particles = 1;
            auto reverse_ions = ions;
            reverse_ions.particles = 1;
            auto reverse_cfg = cfg;
            reverse_cfg.species = {
                reverse_electrons, reverse_ions};
            reverse_cfg.potential_reference.reset();
            reverse_cfg.current_regulated_source
                ->control_mode =
                    pic::CurrentSourceControlMode::Cumulative;
            pic::Simulation2D cumulative_reverse(reverse_cfg);
            cumulative_reverse.initialize();
            cumulative_reverse.step();
            cumulative_reverse.step();
            require_near(
                cumulative_reverse
                    .current_regulated_source_diagnostics()
                    ->control_macro_remainder,
                -1.0, 1e-14,
                "2D cumulative current control lost signed debt");

            reverse_cfg.current_regulated_source
                ->control_mode =
                    pic::CurrentSourceControlMode::TimestepLocal;
            pic::Simulation2D timestep_local_reverse(reverse_cfg);
            timestep_local_reverse.initialize();
            timestep_local_reverse.step();
            timestep_local_reverse.step();
            require(
                timestep_local_reverse
                        .current_regulated_source_diagnostics()
                        ->macro_particles_created == 0,
                "2D timestep-local current control emitted reverse demand");
            require_near(
                timestep_local_reverse
                    .current_regulated_source_diagnostics()
                    ->control_macro_remainder,
                0.0, 1e-14,
                "2D timestep-local current control retained signed debt");
            require(
                timestep_local_reverse
                        .current_regulated_source_diagnostics()
                        ->reverse_demand_steps > 0 &&
                    timestep_local_reverse
                        .current_regulated_source_diagnostics()
                        ->cumulative_reverse_demand_macroparticles > 0.0 &&
                    timestep_local_reverse
                        .current_regulated_source_diagnostics()
                        ->maximum_reverse_demand_macroparticles > 0.0 &&
                    timestep_local_reverse
                        .current_regulated_source_diagnostics()
                        ->reverse_one_macro_steps > 0 &&
                    timestep_local_reverse
                        .current_regulated_source_diagnostics()
                        ->distributed_reverse_demand_macroparticles > 0.0 &&
                    timestep_local_reverse
                        .current_regulated_source_diagnostics()
                        ->squared_reverse_demand_macroparticles > 0.0 &&
                    timestep_local_reverse
                        .current_regulated_source_diagnostics()
                        ->reverse_monitored_positive_charge > 0.0,
                "2D timestep-local current control did not audit reverse demand");

            auto affine_cfg = cfg;
            affine_cfg.boundary_x = pic::Boundary::Dirichlet;
            affine_cfg.boundary_y = pic::Boundary::Periodic;
            affine_cfg.current_regulated_source.reset();
            affine_cfg.potential_reference.reset();
            affine_cfg.species = {reverse_electrons};
            pic::Simulation2D unreferenced(affine_cfg);
            unreferenced.initialize();
            double unreferenced_mean = 0.0;
            for (std::size_t j = 0; j < affine_cfg.ny; ++j) {
                unreferenced_mean +=
                    unreferenced.mesh().phi()[
                        unreferenced.mesh().index(4, j)];
            }
            unreferenced_mean /= static_cast<double>(affine_cfg.ny);
            const double affine_offset = unreferenced_mean - 2.0;
            affine_cfg.potential_reference =
                pic::PotentialReference2DConfig{
                    pic::CoordinateAxis::X, 0.8, 2.0,
                    pic::PotentialReferenceCorrection::Affine};
            pic::Simulation2D affine(affine_cfg);
            affine.initialize();
            double corrected_mean = 0.0;
            for (std::size_t j = 0; j < affine_cfg.ny; ++j) {
                const auto left = affine.mesh().index(0, j);
                const auto reference = affine.mesh().index(4, j);
                corrected_mean += affine.mesh().phi()[reference];
                require_near(
                    affine.mesh().phi()[left],
                    unreferenced.mesh().phi()[left], 1e-12,
                    "2D affine correction changed the preserved electrode");
                require_near(
                    affine.mesh().electric_x()[reference] -
                        unreferenced.mesh().electric_x()[reference],
                    affine_offset / 0.8, 1e-12,
                    "2D affine correction applied the wrong electric field");
            }
            corrected_mean /= static_cast<double>(affine_cfg.ny);
            require_near(
                corrected_mean, 2.0, 1e-12,
                "2D affine potential correction missed its target");
            std::filesystem::remove_all(output_dir);
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
                auto bad_units_cfg = cfg;
                bad_units_cfg.units.system = pic::UnitSystem::SI;
                pic::Simulation bad_units(bad_units_cfg);
                bad_units.load_checkpoint(checkpoint_path);
            }, "1D checkpoint accepted a different unit system");

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
            cfg.magnetic_field_x = 0.07;
            cfg.magnetic_field_y = -0.03;
            cfg.magnetic_field_z = 0.11;
            cfg.species[0].drift_velocity_z = 0.04;
            cfg.species[1].drift_velocity_z = -0.02;

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
            const auto output_dir =
                std::filesystem::path("test_output_2d_pair_source");
            const auto checkpoint_path =
                output_dir / "pair_source.apc";
            std::filesystem::remove_all(output_dir);

            pic::Species2DConfig electrons;
            electrons.name = "electrons";
            electrons.charge = -1.0;
            electrons.mass = 1.0;
            electrons.weight = 3.0;
            electrons.particles = 4;
            electrons.thermal_velocity = 0.0;
            electrons.initialization.loading =
                pic::ParticleLoading::QuietStart;
            pic::Species2DConfig ions = electrons;
            ions.name = "ions";
            ions.charge = 1.0;
            ions.mass = 10.0;

            pic::VolumetricPairSource2DConfig source;
            source.name = "pair_seed";
            source.first_species = "electrons";
            source.second_species = "ions";
            source.pairs_per_step = 2;
            source.start_step = 1;
            source.end_step = 3;
            source.x_min = 0.2;
            source.x_max = 0.4;
            source.y_min = 0.6;
            source.y_max = 0.8;
            source.first_drift = {0.0, 0.0, 0.03};
            source.second_drift = {0.0, 0.0, -0.03};

            pic::Simulation2DConfig cfg;
            cfg.nx = 8;
            cfg.ny = 8;
            cfg.dt = 0.01;
            cfg.steps = 4;
            cfg.output_interval = 1;
            cfg.output_dir = output_dir;
            cfg.seed = 2027;
            cfg.max_particles_per_species = 12;
            cfg.species = {electrons, ions};
            cfg.sources = {source};

            pic::Simulation2D continuous(cfg);
            continuous.initialize();
            continuous.step();
            require(
                continuous.species()[0].live_count() == 4,
                "2D pair source ignored start_step");
            continuous.step();
            require(
                continuous.species()[0].live_count() == 6 &&
                    continuous.species()[1].live_count() == 6,
                "2D pair source did not create paired species");
            require(
                continuous.source_diagnostics()[0]
                        .macro_pairs_created == 2,
                "2D pair source macro-pair accounting is wrong");
            require_near(
                continuous.source_diagnostics()[0]
                    .represented_pairs_created,
                6.0, 1e-14,
                "2D pair source represented-pair accounting is wrong");
            require_near(
                continuous.sample().charge_l1, 0.0, 1e-14,
                "2D colocated pair source introduced charge");
            for (std::size_t particle = 4;
                 particle < 6; ++particle) {
                const auto& first =
                    continuous.species()[0].particles()[particle];
                const auto& second =
                    continuous.species()[1].particles()[particle];
                require_near(
                    first.position.x, second.position.x, 1e-14,
                    "2D pair source did not share pair position x");
                require_near(
                    first.position.y, second.position.y, 1e-14,
                    "2D pair source did not share pair position y");
                require(
                    first.position.x >= source.x_min &&
                        first.position.x <= source.x_max &&
                        first.position.y >= source.y_min &&
                        first.position.y <= source.y_max,
                    "2D pair source sampled outside its region");
            }

            continuous.save_checkpoint(checkpoint_path);
            require(
                read_file_text(checkpoint_path).find(
                    "AuroraPIC-checkpoint-v10\n") == 0,
                "2D pair source checkpoint did not use v10");
            for (std::size_t step = continuous.step_count();
                 step < cfg.steps; ++step) {
                continuous.step();
            }
            require(
                continuous.source_diagnostics()[0]
                        .macro_pairs_created == 4 &&
                    continuous.species()[0].live_count() == 8 &&
                    continuous.species()[1].live_count() == 8,
                "2D pair source schedule created the wrong total");

            pic::Simulation2D restarted(cfg);
            restarted.load_checkpoint(checkpoint_path);
            for (std::size_t step = restarted.step_count();
                 step < cfg.steps; ++step) {
                restarted.step();
            }
            require_species_close(
                continuous.species(), restarted.species(),
                "2D pair source checkpoint determinism");
            require(
                restarted.source_diagnostics()[0]
                        .macro_pairs_created == 4 &&
                    restarted.source_diagnostics()[0]
                        .represented_pairs_created == 12.0,
                "2D pair source checkpoint lost source counters");
            auto changed_source = cfg;
            changed_source.sources[0].pairs_per_step = 1;
            require_throws_contains(
                [&] {
                    pic::Simulation2D incompatible(changed_source);
                    incompatible.load_checkpoint(checkpoint_path);
                },
                "source metadata",
                "2D pair-source checkpoint accepted changed source physics");

            auto bounded = cfg;
            bounded.max_particles_per_species = 5;
            pic::Simulation2D bounded_simulation(bounded);
            bounded_simulation.initialize();
            bounded_simulation.step();
            require_throws_contains(
                [&] { bounded_simulation.step(); },
                "max_particles_per_species",
                "2D pair source ignored its storage bound");
            require(
                bounded_simulation.species()[0].particles().size() == 4 &&
                    bounded_simulation.species()[1].particles().size() == 4,
                "2D pair source capacity failure created a partial pair");

            auto unequal_weight = cfg;
            unequal_weight.species[1].weight = 2.0;
            require_throws_contains(
                [&] {
                    pic::Simulation2D invalid(unequal_weight);
                },
                "equal macro-particle weights",
                "2D pair source accepted unequal macro weights");

            std::filesystem::remove_all(output_dir);
        }
        {
            const auto output_dir = std::filesystem::path(
                "test_output_2d_physical_rate_source");
            const auto checkpoint_path =
                output_dir / "physical_rate.apc";
            std::filesystem::remove_all(output_dir);

            pic::Species2DConfig first;
            first.name = "negative";
            first.charge = -1.0;
            first.mass = 1.0;
            first.weight = 2.0;
            first.particles = 4;
            first.thermal_velocity = 0.0;
            first.initialization.loading =
                pic::ParticleLoading::QuietStart;
            pic::Species2DConfig second = first;
            second.name = "positive";
            second.charge = 1.0;
            second.mass = 4.0;

            pic::Species2DConfig deposition_config = first;
            deposition_config.particles = 1;
            pic::Species2D deposition_species(deposition_config);
            deposition_species.particles() = {
                pic::Particle2D{{0.35, 0.45}, {}, true, {}}};
            pic::Mesh2D unit_depth_mesh(
                8, 8, 1.0, 1.0, pic::Boundary::Periodic);
            pic::Mesh2D double_depth_mesh(
                8, 8, 1.0, 1.0, pic::Boundary::Periodic);
            deposition_species.deposit_charge(
                unit_depth_mesh, 1.0);
            deposition_species.deposit_charge(
                double_depth_mesh, 2.0);
            double unit_depth_charge = 0.0;
            double double_depth_charge = 0.0;
            for (std::size_t j = 0; j < 8; ++j) {
                for (std::size_t i = 0; i < 8; ++i) {
                    const auto index =
                        unit_depth_mesh.index(i, j);
                    require_near(
                        unit_depth_mesh.rho()[index],
                        2.0 * double_depth_mesh.rho()[index],
                        1e-14,
                        "2D extrusion depth did not scale volume charge density");
                    unit_depth_charge +=
                        unit_depth_mesh.rho()[index] *
                        unit_depth_mesh.node_area(i, j);
                    double_depth_charge +=
                        double_depth_mesh.rho()[index] *
                        double_depth_mesh.node_area(i, j) * 2.0;
                }
            }
            require_near(
                unit_depth_charge, -2.0, 1e-14,
                "unit-depth deposition lost represented charge");
            require_near(
                double_depth_charge, -2.0, 1e-14,
                "explicit-depth deposition lost represented charge");

            pic::VolumetricPairSource2DConfig source;
            source.name = "physical_rate";
            source.first_species = first.name;
            source.second_species = second.name;
            source.represented_pair_rate = 5.0;
            source.first_drift.z = 2.0;
            source.second_drift.z = 3.0;

            pic::Simulation2DConfig cfg;
            cfg.nx = 8;
            cfg.ny = 8;
            cfg.dt = 0.1;
            cfg.steps = 10;
            cfg.output_interval = 10;
            cfg.output_dir = output_dir;
            cfg.max_particles_per_species = 16;
            cfg.species = {first, second};
            cfg.sources = {source};

            pic::Simulation2D continuous(cfg);
            continuous.initialize();
            for (std::size_t step = 0; step < 5; ++step) {
                continuous.step();
            }
            require(
                continuous.source_diagnostics()[0]
                        .macro_pairs_created == 1,
                "2D physical-rate source rounded instead of accumulating");
            require_near(
                continuous.source_diagnostics()[0]
                    .fractional_macro_pair_remainder,
                0.25, 1e-14,
                "2D physical-rate source has the wrong fractional remainder");
            continuous.save_checkpoint(checkpoint_path);
            for (std::size_t step = 5; step < 10; ++step) {
                continuous.step();
            }
            require(
                continuous.source_diagnostics()[0]
                        .macro_pairs_created == 2,
                "2D physical-rate source produced the wrong ten-step count");
            require_near(
                continuous.source_diagnostics()[0]
                    .fractional_macro_pair_remainder,
                0.5, 1e-14,
                "2D physical-rate source lost its ten-step remainder");
            require_near(
                continuous.source_diagnostics()[0]
                    .represented_pairs_created,
                4.0, 1e-14,
                "2D physical-rate represented-pair accounting is wrong");
            require_near(
                continuous.source_diagnostics()[0]
                    .injected_kinetic_energy,
                80.0, 1e-12,
                "2D source-energy accounting is wrong");

            pic::Simulation2D restarted(cfg);
            restarted.load_checkpoint(checkpoint_path);
            for (std::size_t step = restarted.step_count();
                 step < 10; ++step) {
                restarted.step();
            }
            require_species_close(
                continuous.species(), restarted.species(),
                "2D physical-rate source restart");
            require_near(
                restarted.source_diagnostics()[0]
                    .fractional_macro_pair_remainder,
                0.5, 1e-14,
                "2D source restart lost its fractional accumulator");
            require_near(
                restarted.source_diagnostics()[0]
                    .injected_kinetic_energy,
                80.0, 1e-12,
                "2D source restart lost energy accounting");

            auto refined_cfg = cfg;
            refined_cfg.dt = 0.05;
            refined_cfg.steps = 20;
            pic::Simulation2D refined(refined_cfg);
            refined.initialize();
            for (std::size_t step = 0; step < 20; ++step) {
                refined.step();
            }
            require(
                refined.source_diagnostics()[0]
                        .macro_pairs_created == 2,
                "2D physical source changed production under timestep refinement");
            require_near(
                refined.source_diagnostics()[0]
                    .fractional_macro_pair_remainder,
                0.5, 1e-14,
                "2D physical source remainder changed under timestep refinement");
            require_near(
                refined.source_diagnostics()[0]
                    .represented_pairs_created,
                continuous.source_diagnostics()[0]
                    .represented_pairs_created,
                1e-14,
                "2D physical represented rate is not timestep invariant");

            auto volumetric_cfg = cfg;
            volumetric_cfg.out_of_plane_depth = 2.0;
            volumetric_cfg.sources[0].represented_pair_rate.reset();
            volumetric_cfg.sources[0].peak_volumetric_pair_rate =
                5.0;
            volumetric_cfg.sources[0].spatial_profile
                .density_profile =
                pic::DensityProfileKind::Sinusoidal;
            volumetric_cfg.sources[0].spatial_profile
                .profile_amplitude = -1.0;
            volumetric_cfg.sources[0].spatial_profile
                .profile_mode_x = 1;
            pic::Simulation2D volumetric(volumetric_cfg);
            volumetric.initialize();
            for (std::size_t step = 0; step < 10; ++step) {
                volumetric.step();
            }
            require(
                volumetric.source_diagnostics()[0]
                        .macro_pairs_created == 2,
                "2D peak volumetric source used the wrong profile integral");
            require_near(
                volumetric.source_diagnostics()[0]
                    .fractional_macro_pair_remainder,
                0.5, 1e-14,
                "2D peak volumetric conversion has the wrong remainder");

            auto both_rates = cfg;
            both_rates.sources[0].pairs_per_step = 1;
            require_throws_contains(
                [&] { pic::Simulation2D invalid(both_rates); },
                "exactly one rate",
                "2D source accepted two rate specifications");
            auto no_rate = cfg;
            no_rate.sources[0].represented_pair_rate.reset();
            require_throws_contains(
                [&] { pic::Simulation2D invalid(no_rate); },
                "exactly one rate",
                "2D source accepted no rate specification");
            std::filesystem::remove_all(output_dir);
        }
        {
            pic::Species2DConfig negative;
            negative.name = "negative";
            negative.charge = -1.0;
            negative.mass = 1.0;
            negative.weight = 1.0;
            negative.particles = 4;
            negative.thermal_velocity = 0.0;
            negative.initialization.loading =
                pic::ParticleLoading::QuietStart;
            auto positive = negative;
            positive.name = "positive";
            positive.charge = 1.0;

            pic::VolumetricPairSource2DConfig source;
            source.name = "sinusoidal_profile";
            source.first_species = negative.name;
            source.second_species = positive.name;
            source.pairs_per_step = 2000;
            source.x_min = 0.2;
            source.x_max = 0.8;
            source.spatial_profile.density_profile =
                pic::DensityProfileKind::Sinusoidal;
            source.spatial_profile.profile_amplitude = -1.0;
            source.spatial_profile.profile_mode_x = 1;
            source.spatial_profile
                .max_profile_sampling_attempts = 10000;

            pic::Simulation2DConfig cfg;
            cfg.nx = 8;
            cfg.ny = 8;
            cfg.dt = 0.01;
            cfg.max_particles_per_species = 2100;
            cfg.species = {negative, positive};
            cfg.sources = {source};
            cfg.seed = 90817;

            pic::Simulation2D simulation(cfg);
            simulation.initialize();
            simulation.step();
            double cosine_moment = 0.0;
            for (std::size_t particle = 4;
                 particle < simulation.species()[0]
                                .particles().size();
                 ++particle) {
                const auto& first_particle =
                    simulation.species()[0].particles()[particle];
                const auto& second_particle =
                    simulation.species()[1].particles()[particle];
                require_near(
                    first_particle.position.x,
                    second_particle.position.x, 1e-14,
                    "profiled source separated a colocated pair");
                const double normalized =
                    (first_particle.position.x - source.x_min) /
                    (source.x_max - source.x_min);
                cosine_moment += std::cos(
                    2.0 * std::numbers::pi * normalized);
            }
            cosine_moment /= 2000.0;
            require(
                cosine_moment > -0.56 &&
                    cosine_moment < -0.44,
                "2D sinusoidal source does not follow its normalized profile");
            require_near(
                simulation.sample().charge_l1, 0.0, 1e-12,
                "profiled pair source introduced net deposited charge");
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

            const auto imported_example =
                std::filesystem::path(AURORA_TEST_SOURCE_DIR) /
                "examples" / "imported_plasma_2d.cfg";
            require(pic::config_uses_unstructured_mesh_2d(imported_example),
                    "imported 2D config detection did not select the unstructured runtime");
            const auto imported_config =
                pic::load_unstructured_config_2d(imported_example);
            require(imported_config.mesh_path.filename() == "imported_square_v2.msh" &&
                        imported_config.species.size() == 2 &&
                        imported_config.dirichlet_potentials.size() == 1 &&
                        imported_config.dirichlet_potentials.contains("inlet") &&
                        imported_config.neumann_normal_derivatives.size() == 3 &&
                        imported_config.neumann_normal_derivatives.contains("outlet") &&
                        imported_config.sources.size() == 2 &&
                        imported_config.sources.front().species == "electrons" &&
                        imported_config.emissions.size() == 1 &&
                        imported_config.emissions.front().incident_species ==
                            "ions" &&
                        imported_config.particle_boundaries.at("electrode") ==
                        pic::ParticleBoundary::Reflecting,
                    "imported 2D example config did not load expected runtime settings");
            auto profiled_imported_config = imported_config;
            profiled_imported_config.magnetic_field_profile =
                pic::TabulatedVectorField1D(
                    pic::CoordinateAxis::X,
                    {0.0, 0.5, 1.0},
                    {{0.0, 0.0, 0.5},
                     {0.0, 0.0, 1.0},
                     {0.0, 0.0, 1.5}});
            pic::UnstructuredSimulation2D profiled_imported_simulation(
                profiled_imported_config);
            (void)profiled_imported_simulation;
            const auto imported_mcc_example =
                std::filesystem::path(AURORA_TEST_SOURCE_DIR) /
                "examples" / "imported_mcc_2d.cfg";
            const auto imported_mcc_config =
                pic::load_unstructured_config_2d(imported_mcc_example);
            require(
                imported_mcc_config.collisions.enabled &&
                    imported_mcc_config.collisions.model ==
                        pic::CollisionModelKind::NullCollision &&
                    imported_mcc_config.collisions.species == "tracers" &&
                    imported_mcc_config.collisions.gas_name ==
                        "synthetic_validation_gas" &&
                    imported_mcc_config.collisions.neutral_mass == 40.0 &&
                    imported_mcc_config.collisions.neutral_temperature ==
                        300.0 &&
                    imported_mcc_config.collisions.channels.size() == 1 &&
                    imported_mcc_config.collisions.channels.front()
                            .cross_section_file.filename() ==
                        "mcc_2d3v_elastic.dat",
                "imported 2D MCC config did not preserve gas and channel metadata");
            const auto imported_ionization_example =
                std::filesystem::path(AURORA_TEST_SOURCE_DIR) /
                "examples" / "imported_ionization_2d.cfg";
            const auto imported_ionization_config =
                pic::load_unstructured_config_2d(
                    imported_ionization_example);
            require(
                imported_ionization_config.collisions.gas_data_file
                        .filename() == "synthetic_ionization.gas" &&
                    imported_ionization_config.collisions.gas_data_version ==
                        2 &&
                    imported_ionization_config.collisions.gas_data_units ==
                        pic::UnitSystem::Normalized &&
                    imported_ionization_config.collisions.dataset_id ==
                        "aurorapic.synthetic.ionization" &&
                    imported_ionization_config.collisions.dataset_version ==
                        "1" &&
                    imported_ionization_config.collisions.retrieved ==
                        "2026-07-28" &&
                    !imported_ionization_config.collisions.citation.empty() &&
                    !imported_ionization_config.collisions.license.empty() &&
                imported_ionization_config.collisions.channels.size() == 1 &&
                    imported_ionization_config.collisions.channels.front()
                            .process ==
                        pic::CollisionProcessKind::Ionization &&
                    imported_ionization_config.collisions.channels.front()
                            .secondary_species == "electrons" &&
                    imported_ionization_config.collisions.channels.front()
                        .ion_species == "ions",
                "imported ionization config did not preserve product species");
            const auto imported_attachment_example =
                std::filesystem::path(AURORA_TEST_SOURCE_DIR) /
                "examples" / "imported_attachment_2d.cfg";
            const auto imported_attachment_config =
                pic::load_unstructured_config_2d(
                    imported_attachment_example);
            require(
                imported_attachment_config.collisions.dataset_id ==
                        "aurorapic.synthetic.attachment" &&
                    imported_attachment_config.collisions.channels
                            .size() == 1 &&
                    imported_attachment_config.collisions.channels
                            .front().process ==
                        pic::CollisionProcessKind::Attachment &&
                    imported_attachment_config.collisions.channels
                            .front().attachment_species ==
                        "negative_ions",
                "imported attachment config lost its product mapping");
            const auto imported_charge_exchange_example =
                std::filesystem::path(AURORA_TEST_SOURCE_DIR) /
                "examples" / "imported_charge_exchange_2d.cfg";
            const auto imported_charge_exchange_config =
                pic::load_unstructured_config_2d(
                    imported_charge_exchange_example);
            require(
                imported_charge_exchange_config.collisions.dataset_id ==
                        "aurorapic.synthetic.charge_exchange" &&
                    imported_charge_exchange_config.collisions
                            .neutral_mass == 40.0 &&
                    imported_charge_exchange_config.collisions.channels
                            .size() == 1 &&
                    imported_charge_exchange_config.collisions.channels
                            .front().process ==
                        pic::CollisionProcessKind::ChargeExchange,
                "imported charge-exchange gas dataset did not load");
            {
                const auto gas_dataset_path =
                    std::filesystem::path(AURORA_TEST_SOURCE_DIR) /
                    "examples" / "synthetic_ionization.gas";
                const auto gas_dataset =
                    pic::load_gas_dataset(gas_dataset_path);
                require(
                    gas_dataset.format_version == 2 &&
                        gas_dataset.unit_system ==
                            pic::UnitSystem::Normalized &&
                        gas_dataset.gas_name ==
                            "synthetic_validation_gas" &&
                        gas_dataset.neutral_mass == 40.0 &&
                        gas_dataset.channels.size() == 1 &&
                        gas_dataset.channels.front().process ==
                            pic::CollisionProcessKind::Ionization,
                    "standalone gas dataset loader lost manifest data");
                const auto angular_gas_path =
                    std::filesystem::path(AURORA_TEST_SOURCE_DIR) /
                    "examples" / "synthetic_swarm.gas";
                const auto angular_gas =
                    pic::load_gas_dataset(angular_gas_path);
                require(
                    angular_gas.channels.size() == 3 &&
                        angular_gas.channels.front()
                                .angular_scattering ==
                            pic::AngularScatteringKind::
                                HenyeyGreenstein &&
                        angular_gas.channels.front()
                                .mean_cosine_file.filename() ==
                            "synthetic_swarm_mean_cosine.dat",
                    "gas dataset loader lost angular scattering data");

                const auto turner_table_path =
                    std::filesystem::path(
                        "test_turner_ion_cross_section.dat");
                const auto turner_gas_path =
                    std::filesystem::path(
                        "test_turner_ion.gas");
                {
                    std::ofstream table(turner_table_path);
                    table << "0 1e-19\n10 1e-19\n";
                    std::ofstream manifest(turner_gas_path);
                    manifest
                        << "gas_data_version = 2\n"
                        << "units = si\n"
                        << "gas = helium\n"
                        << "neutral_mass = 6.67e-27\n"
                        << "dataset_id = test.turner.he-ion\n"
                        << "dataset_version = 1\n"
                        << "data_provenance = unit test\n"
                        << "citation = Turner et al. 2013\n"
                        << "retrieved = 2026-07-30\n"
                        << "license = test data\n"
                        << "[collision.isotropic]\n"
                        << "type = elastic\n"
                        << "cross_section_file = "
                        << turner_table_path.string() << "\n"
                        << "energy_scale = 1.602176634e-19\n"
                        << "energy_frame = center_of_mass\n"
                        << "[collision.backward]\n"
                        << "type = elastic\n"
                        << "cross_section_file = "
                        << turner_table_path.string() << "\n"
                        << "energy_scale = 1.602176634e-19\n"
                        << "energy_frame = center_of_mass\n"
                        << "angular_model = backward\n";
                }
                const auto turner_gas =
                    pic::load_gas_dataset(turner_gas_path);
                require(
                    turner_gas.channels.size() == 2 &&
                        turner_gas.channels[0].energy_frame ==
                            pic::CollisionEnergyFrame::CenterOfMass &&
                        turner_gas.channels[1].angular_scattering ==
                            pic::AngularScatteringKind::Backward,
                    "gas dataset loader lost Turner ion energy-frame "
                    "or backward-scattering data");
                const auto turner_config_path =
                    std::filesystem::path(
                        "test_turner_named_gas_config.cfg");
                {
                    std::ofstream config(turner_config_path);
                    config
                        << "config_version = 1\n"
                        << "units = si\n"
                        << "dimension = 1\n"
                        << "velocity_dimensions = 3\n"
                        << "nx = 3\n"
                        << "length = 0.067\n"
                        << "dt = 1e-12\n"
                        << "steps = 1\n"
                        << "output_interval = 1\n"
                        << "boundary = dirichlet\n"
                        << "mode = transient\n"
                        << "[collisions.ion_mcc]\n"
                        << "model = null_collision\n"
                        << "species = ions\n"
                        << "neutral_density = 9.64e20\n"
                        << "neutral_temperature = 300\n"
                        << "max_frequency = 1e9\n"
                        << "gas_data_file = "
                        << turner_gas_path.string() << "\n"
                        << "[species.ions]\n"
                        << "charge = 1.602176634e-19\n"
                        << "mass = 6.67e-27\n"
                        << "weight = 1\n"
                        << "particles = 2\n"
                        << "thermal_velocity = 0\n";
                }
                const auto turner_config =
                    pic::load_config(turner_config_path.string());
                require(
                    turner_config.collision_models.size() == 1 &&
                        turner_config.collision_models[0]
                                .config.dataset_id ==
                            "test.turner.he-ion" &&
                        turner_config.collision_models[0]
                                .config.channels.size() == 2 &&
                        turner_config.collision_models[0]
                                .config.gas_data_file.filename() ==
                            turner_gas_path.filename(),
                    "1D named collision model lost its gas dataset "
                    "identity or channels");
                std::filesystem::remove(turner_config_path);

                const auto reactive_config_path =
                    std::filesystem::path(
                        "test_named_reactive_gas_config.cfg");
                {
                    std::ofstream config(reactive_config_path);
                    config
                        << "config_version = 1\n"
                        << "units = normalized\n"
                        << "dimension = 1\n"
                        << "velocity_dimensions = 3\n"
                        << "nx = 3\n"
                        << "length = 1\n"
                        << "dt = 0.01\n"
                        << "steps = 1\n"
                        << "output_interval = 1\n"
                        << "boundary = periodic\n"
                        << "mode = transient\n"
                        << "[collisions.electron_mcc]\n"
                        << "model = null_collision\n"
                        << "species = electrons\n"
                        << "neutral_density = 1\n"
                        << "neutral_temperature = 0\n"
                        << "max_frequency = 1\n"
                        << "gas_data_file = "
                        << (std::filesystem::path(
                                AURORA_TEST_SOURCE_DIR) /
                            "examples" /
                            "synthetic_ionization.gas").string()
                        << "\n"
                        << "[collisions.electron_mcc.channel."
                           "synthetic_ionization]\n"
                        << "secondary_species = electrons\n"
                        << "ion_species = ions\n"
                        << "[species.electrons]\n"
                        << "charge = -1\n"
                        << "mass = 1\n"
                        << "weight = 1\n"
                        << "particles = 2\n"
                        << "[species.ions]\n"
                        << "charge = 1\n"
                        << "mass = 40\n"
                        << "weight = 1\n"
                        << "particles = 2\n";
                }
                const auto reactive_config =
                    pic::load_config(reactive_config_path.string());
                require(
                    reactive_config.collision_models.size() == 1 &&
                        reactive_config.collision_models[0]
                                .config.channels.size() == 1 &&
                        reactive_config.collision_models[0]
                                .config.channels[0].secondary_species ==
                            "electrons" &&
                        reactive_config.collision_models[0]
                                .config.channels[0].ion_species ==
                            "ions",
                    "1D gas dataset ionization product mapping was lost");
                std::filesystem::remove(reactive_config_path);
                std::filesystem::remove(turner_gas_path);
                std::filesystem::remove(turner_table_path);

                const auto invalid_dataset =
                    std::filesystem::path(
                        "test_invalid_retrieval_date.gas");
                std::string invalid_text =
                    read_file_text(gas_dataset_path);
                const auto date =
                    invalid_text.find("retrieved = 2026-07-28");
                require(
                    date != std::string::npos,
                    "gas dataset fixture retrieval date is missing");
                invalid_text.replace(
                    date, std::string("retrieved = 2026-07-28").size(),
                    "retrieved = 2026-02-30");
                {
                    std::ofstream invalid(invalid_dataset);
                    invalid << invalid_text;
                }
                require_throws_contains(
                    [&] {
                        try {
                            (void)pic::load_gas_dataset(
                                invalid_dataset);
                        } catch (...) {
                            std::filesystem::remove(
                                invalid_dataset);
                            throw;
                        }
                        std::filesystem::remove(invalid_dataset);
                    },
                    "valid YYYY-MM-DD",
                    "gas dataset accepted an invalid retrieval date");

                const auto invalid_config =
                    std::filesystem::path(
                        "test_gas_dataset_override.cfg");
                std::string invalid_config_text =
                    read_file_text(imported_ionization_example);
                const auto manifest_reference =
                    invalid_config_text.find(
                        "gas_data_file = synthetic_ionization.gas");
                require(
                    manifest_reference != std::string::npos,
                    "external gas example manifest reference is missing");
                invalid_config_text.replace(
                    manifest_reference,
                    std::string(
                        "gas_data_file = synthetic_ionization.gas").size(),
                    "gas_data_file = " + gas_dataset_path.string());
                const auto product_mapping =
                    invalid_config_text.find(
                        "secondary_species = electrons");
                require(
                    product_mapping != std::string::npos,
                    "external gas example product mapping is missing");
                invalid_config_text.insert(
                    product_mapping,
                    "threshold_energy = 2.0\n");
                {
                    std::ofstream invalid(invalid_config);
                    invalid << invalid_config_text;
                }
                require_throws_contains(
                    [&] {
                        try {
                            (void)pic::load_unstructured_config_2d(
                                invalid_config);
                        } catch (...) {
                            std::filesystem::remove(
                                invalid_config);
                            throw;
                        }
                        std::filesystem::remove(invalid_config);
                    },
                    "physics cannot be overridden",
                    "simulation config overrode packaged gas physics");

                const auto si_manifest =
                    std::filesystem::absolute(
                        "test_si_gas_dataset.gas");
                const auto unit_mismatch_config =
                    std::filesystem::path(
                        "test_gas_dataset_unit_mismatch.cfg");
                std::string si_manifest_text =
                    read_file_text(gas_dataset_path);
                const auto units = si_manifest_text.find(
                    "units = normalized");
                require(
                    units != std::string::npos,
                    "gas dataset fixture unit contract is missing");
                si_manifest_text.replace(
                    units, std::string("units = normalized").size(),
                    "units = si");
                const auto table = si_manifest_text.find(
                    "cross_section_file = mcc_2d3v_ionization.dat");
                require(
                    table != std::string::npos,
                    "gas dataset fixture table reference is missing");
                si_manifest_text.replace(
                    table,
                    std::string(
                        "cross_section_file = "
                        "mcc_2d3v_ionization.dat").size(),
                    "cross_section_file = " +
                        (gas_dataset_path.parent_path() /
                         "mcc_2d3v_ionization.dat").string());
                {
                    std::ofstream manifest(si_manifest);
                    manifest << si_manifest_text;
                }
                std::string mismatch_text =
                    read_file_text(imported_ionization_example);
                const auto mismatch_reference =
                    mismatch_text.find(
                        "gas_data_file = synthetic_ionization.gas");
                require(
                    mismatch_reference != std::string::npos,
                    "gas dataset fixture reference is missing");
                mismatch_text.replace(
                    mismatch_reference,
                    std::string(
                        "gas_data_file = synthetic_ionization.gas").size(),
                    "gas_data_file = " + si_manifest.string());
                {
                    std::ofstream mismatch(unit_mismatch_config);
                    mismatch << mismatch_text;
                }
                require_throws_contains(
                    [&] {
                        try {
                            (void)pic::load_unstructured_config_2d(
                                unit_mismatch_config);
                        } catch (...) {
                            std::filesystem::remove(si_manifest);
                            std::filesystem::remove(
                                unit_mismatch_config);
                            throw;
                        }
                        std::filesystem::remove(si_manifest);
                        std::filesystem::remove(
                            unit_mismatch_config);
                    },
                    "do not match simulation units",
                    "simulation accepted mismatched gas dataset units");

                const auto swarm_table =
                    std::filesystem::absolute(
                        "test_swarm_zero_elastic.dat");
                const auto swarm_manifest =
                    std::filesystem::absolute(
                        "test_swarm_zero_elastic.gas");
                const auto swarm_config_path =
                    std::filesystem::absolute(
                        "test_swarm_zero_elastic.cfg");
                const auto swarm_output =
                    std::filesystem::absolute(
                        "test_swarm_zero_elastic.csv");
                const auto swarm_ionization_table =
                    std::filesystem::absolute(
                        "test_swarm_ionization.dat");
                const auto swarm_attachment_table =
                    std::filesystem::absolute(
                        "test_swarm_attachment.dat");
                const auto swarm_branch_manifest =
                    std::filesystem::absolute(
                        "test_swarm_branching.gas");
                const auto swarm_branch_config_path =
                    std::filesystem::absolute(
                        "test_swarm_branching.cfg");
                const auto swarm_branch_output =
                    std::filesystem::absolute(
                        "test_swarm_branching.csv");
                const auto swarm_spatial_output =
                    std::filesystem::absolute(
                        "test_swarm_spatial.csv");
                {
                    std::ofstream table_output(swarm_table);
                    table_output << "0 0\n50 0\n";
                    std::ofstream manifest_output(swarm_manifest);
                    manifest_output
                        << "gas_data_version = 2\n"
                        << "units = si\n"
                        << "gas = synthetic_swarm\n"
                        << "neutral_mass = 6.6335209e-26\n"
                        << "dataset_id = aurorapic.synthetic.swarm\n"
                        << "dataset_version = 1\n"
                        << "data_provenance = AuroraPIC test\n"
                        << "citation = Synthetic fixture\n"
                        << "retrieved = 2026-07-28\n"
                        << "license = Synthetic test data\n\n"
                        << "[collision.elastic]\n"
                        << "type = elastic\n"
                        << "cross_section_file = "
                        << swarm_table.string() << '\n'
                        << "energy_scale = 1.602176634e-19\n";
                    std::ofstream config_output(swarm_config_path);
                    config_output
                        << "swarm_config_version = 1\n"
                        << "gas_data_file = "
                        << swarm_manifest.string() << '\n'
                        << "neutral_density = 1e20\n"
                        << "neutral_temperature = 300\n"
                        << "reduced_fields_td = 1\n"
                        << "max_frequency = 1e6\n"
                        << "timestep = 1e-8\n"
                        << "steps = 6\n"
                        << "burn_in_steps = 2\n"
                        << "particles = 4\n"
                        << "uncertainty_blocks = 2\n"
                        << "initial_mean_energy_ev = 0\n"
                        << "max_energy_ev = 1\n"
                        << "seed = 17\n"
                        << "output_file = "
                        << swarm_output.string() << '\n';
                    std::ofstream ionization_output(
                        swarm_ionization_table);
                    ionization_output
                        << "0 0\n0.01 1e-18\n50 1e-18\n";
                    std::ofstream attachment_output(
                        swarm_attachment_table);
                    attachment_output
                        << "0 0\n0.01 1e-19\n50 1e-19\n";
                    std::ofstream branch_manifest_output(
                        swarm_branch_manifest);
                    branch_manifest_output
                        << "gas_data_version = 2\n"
                        << "units = si\n"
                        << "gas = synthetic_branching_swarm\n"
                        << "neutral_mass = 6.6335209e-26\n"
                        << "dataset_id = aurorapic.synthetic.branching\n"
                        << "dataset_version = 1\n"
                        << "data_provenance = AuroraPIC test\n"
                        << "citation = Synthetic fixture\n"
                        << "retrieved = 2026-07-28\n"
                        << "license = Synthetic test data\n\n"
                        << "[collision.elastic]\n"
                        << "type = elastic\n"
                        << "cross_section_file = "
                        << swarm_table.string() << '\n'
                        << "energy_scale = 1.602176634e-19\n\n"
                        << "[collision.ionization]\n"
                        << "type = ionization\n"
                        << "cross_section_file = "
                        << swarm_ionization_table.string() << '\n'
                        << "energy_scale = 1.602176634e-19\n"
                        << "threshold_energy = 1.602176634e-21\n\n"
                        << "[collision.attachment]\n"
                        << "type = attachment\n"
                        << "cross_section_file = "
                        << swarm_attachment_table.string() << '\n'
                        << "energy_scale = 1.602176634e-19\n";
                    std::ofstream branch_config_output(
                        swarm_branch_config_path);
                    branch_config_output
                        << "swarm_config_version = 1\n"
                        << "gas_data_file = "
                        << swarm_branch_manifest.string() << '\n'
                        << "neutral_density = 1e20\n"
                        << "reduced_fields_td = 1000\n"
                        << "max_frequency = 5e8\n"
                        << "timestep = 2e-10\n"
                        << "steps = 300\n"
                        << "burn_in_steps = 60\n"
                        << "particles = 64\n"
                        << "population_model = branching_resampled\n"
                        << "population_limit = 256\n"
                        << "uncertainty_blocks = 8\n"
                        << "initial_mean_energy_ev = 1\n"
                        << "max_energy_ev = 50\n"
                        << "seed = 29\n"
                        << "output_file = "
                        << swarm_branch_output.string() << '\n'
                        << "spatial_histories = 128\n"
                        << "spatial_length_m = 0.01\n"
                        << "spatial_bins = 8\n"
                        << "spatial_fit_begin_bin = 1\n"
                        << "spatial_fit_end_bin = 7\n"
                        << "spatial_max_steps = 2000\n"
                        << "spatial_work_item_limit = 1000000\n"
                        << "spatial_min_r_squared = 0.5\n"
                        << "spatial_profile_file = "
                        << swarm_spatial_output.string() << '\n';
                }
                try {
                    const auto swarm_config =
                        pic::load_swarm_benchmark_config(
                            swarm_config_path);
                    const auto swarm_results =
                        pic::run_swarm_benchmark(swarm_config);
                    require(
                        swarm_results.size() == 1 &&
                            swarm_results.front().channels.size() == 1 &&
                            swarm_results.front().channels.front()
                                    .collisions == 0 &&
                            swarm_results.front()
                                    .neutral_velocity_stddev_m_s > 200.0 &&
                            swarm_results.front()
                                    .neutral_velocity_stddev_m_s < 300.0 &&
                            swarm_results.front()
                                    .neutral_speed_limit_sigma == 8.0,
                        "zero-cross-section swarm produced collisions");
                    constexpr double electron_mass =
                        9.1093837139e-31;
                    constexpr double elementary_charge =
                        1.602176634e-19;
                    const double acceleration_step =
                        elementary_charge * 0.1 / electron_mass *
                        swarm_config.timestep;
                    const double expected_drift =
                        4.5 * acceleration_step;
                    require_near(
                        swarm_results.front()
                            .electron_drift_velocity_m_s,
                        expected_drift,
                        expected_drift * 1e-13,
                        "swarm uniform-field drift is incorrect");
                    const double expected_mean_energy =
                        0.5 * electron_mass *
                        acceleration_step * acceleration_step *
                        21.5 / elementary_charge;
                    require_near(
                        swarm_results.front().mean_energy_ev,
                        expected_mean_energy,
                        expected_mean_energy * 1e-13,
                        "swarm uniform-field mean energy is incorrect");
                    require_near(
                        swarm_results.front()
                            .longitudinal_diffusion_m2_s,
                        0.0, 1e-30,
                        "identical swarm particles diffused");
                    const auto swarm_dataset =
                        pic::load_gas_dataset(swarm_manifest);
                    pic::write_swarm_benchmark_csv(
                        swarm_output, swarm_config,
                        swarm_dataset, swarm_results);
                    const std::string csv =
                        read_file_text(swarm_output);
                    const auto csv_header_end = csv.find('\n');
                    const auto csv_row_end =
                        csv.find('\n', csv_header_end + 1);
                    require(
                        csv.find(
                            "fixed_population_no_avalanche") !=
                                std::string::npos &&
                            csv.find(
                                "electron_drift_velocity_m_s") !=
                                std::string::npos &&
                            csv_header_end != std::string::npos &&
                            csv_row_end != std::string::npos &&
                            std::count(
                                csv.begin(),
                                csv.begin() +
                                    static_cast<std::ptrdiff_t>(
                                        csv_header_end),
                                ',') ==
                                std::count(
                                    csv.begin() +
                                        static_cast<std::ptrdiff_t>(
                                            csv_header_end + 1),
                                    csv.begin() +
                                        static_cast<std::ptrdiff_t>(
                                            csv_row_end),
                                    ','),
                        "swarm CSV omitted method metadata");

                    auto invalid_swarm = swarm_config;
                    invalid_swarm.max_energy_ev = 51.0;
                    require_throws_contains(
                        [&] {
                            (void)pic::run_swarm_benchmark(
                                invalid_swarm);
                        },
                        "table coverage",
                        "swarm accepted an uncovered energy range");
                    invalid_swarm = swarm_config;
                    invalid_swarm.work_item_limit = 1;
                    require_throws_contains(
                        [&] {
                            (void)pic::run_swarm_benchmark(
                                invalid_swarm);
                        },
                        "work_item_limit",
                        "swarm ignored its conservative work limit");

                    const auto branch_config =
                        pic::load_swarm_benchmark_config(
                            swarm_branch_config_path);
                    require(
                        branch_config.population_model ==
                                pic::SwarmPopulationModel::
                                    BranchingResampled &&
                            branch_config.population_limit == 256,
                        "branching swarm controls were not parsed");
                    const auto branch_results =
                        pic::run_swarm_benchmark(branch_config);
                    require(
                        branch_results.size() == 1,
                        "branching swarm returned the wrong field count");
                    const auto& branch = branch_results.front();
                    require(
                        !branch.diffusion_available &&
                            branch.final_computational_particles == 64 &&
                            branch.final_total_electron_weight >
                                branch.initial_total_electron_weight &&
                            branch.temporal_growth_rate_s > 0.0 &&
                            branch.ionization_rate_s >
                                branch.attachment_rate_s &&
                            branch.attachment_rate_s > 0.0 &&
                            branch.net_creation_rate_s > 0.0 &&
                            branch.spatial_townsend_available &&
                            branch.spatial_flux_townsend_1_m > 0.0 &&
                            branch.spatial_flux_profile.size() == 8 &&
                            branch.spatial_histories_completed == 128 &&
                            branch.spatial_flux_fit_r_squared >= 0.0 &&
                            branch.spatial_flux_fit_r_squared <= 1.0 &&
                            branch
                                    .spatial_maximum_active_particles <=
                                256 &&
                            branch.spatial_particle_updates <= 1000000,
                        "branching swarm did not preserve its bounded "
                        "ensemble while multiplying electron weight");
                    require(
                        branch.townsend_available &&
                            branch.electron_drift_velocity_m_s > 0.0 &&
                            branch
                                    .rate_balance_effective_townsend_1_m >
                                0.0,
                        "branching swarm did not expose its effective "
                        "Townsend diagnostics");
                    require_near(
                        branch.rate_balance_effective_townsend_1_m,
                        branch.net_creation_rate_s /
                            branch.electron_drift_velocity_m_s,
                        std::abs(
                            branch
                                .rate_balance_effective_townsend_1_m) *
                            1e-14,
                        "rate-balance effective Townsend coefficient "
                        "is inconsistent");
                    const auto branch_dataset =
                        pic::load_gas_dataset(swarm_branch_manifest);
                    pic::write_swarm_benchmark_csv(
                        swarm_branch_output, branch_config,
                        branch_dataset, branch_results);
                    pic::write_swarm_spatial_profile_csv(
                        swarm_spatial_output, branch_config,
                        branch_dataset, branch_results);
                    const std::string branch_csv =
                        read_file_text(swarm_branch_output);
                    const auto branch_header_end =
                        branch_csv.find('\n');
                    const auto branch_row_end =
                        branch_csv.find(
                            '\n', branch_header_end + 1);
                    require(
                        branch_csv.find("branching_resampled") !=
                                std::string::npos &&
                            branch_csv.find(
                                "temporal_growth_rate_s") !=
                                std::string::npos &&
                            branch_csv.find(
                                "diffusion_available") !=
                                std::string::npos &&
                            branch_csv.find(
                                "net_creation_rate_s") !=
                                std::string::npos &&
                            branch_csv.find(
                                "total_attachment_rate_s") !=
                                std::string::npos &&
                            branch_header_end != std::string::npos &&
                            branch_row_end != std::string::npos &&
                            std::count(
                                branch_csv.begin(),
                                branch_csv.begin() +
                                    static_cast<std::ptrdiff_t>(
                                        branch_header_end),
                                ',') ==
                                std::count(
                                    branch_csv.begin() +
                                        static_cast<std::ptrdiff_t>(
                                            branch_header_end + 1),
                                    branch_csv.begin() +
                                        static_cast<std::ptrdiff_t>(
                                            branch_row_end),
                                    ','),
                        "branching swarm CSV omitted population "
                        "diagnostics");
                    const std::string spatial_csv =
                        read_file_text(swarm_spatial_output);
                    require(
                        count_lines(spatial_csv) == 9 &&
                            spatial_csv.find(
                                "net_crossings_per_injected_electron") !=
                                std::string::npos &&
                            spatial_csv.find(",yes,") !=
                                std::string::npos,
                        "spatial Townsend profile CSV is incomplete");
                    auto invalid_branch = branch_config;
                    invalid_branch.population_limit =
                        invalid_branch.particles;
                    require_throws_contains(
                        [&] {
                            (void)pic::run_swarm_benchmark(
                                invalid_branch);
                        },
                        "population_limit must exceed particles",
                        "branching swarm accepted an unsafe population "
                        "limit");
                    invalid_branch = branch_config;
                    invalid_branch.spatial_work_item_limit = 1;
                    require_throws_contains(
                        [&] {
                            (void)pic::run_swarm_benchmark(
                                invalid_branch);
                        },
                        "spatial_work_item_limit",
                        "spatial Townsend experiment ignored its "
                        "work limit");
                    invalid_branch = branch_config;
                    invalid_branch.population_model =
                        pic::SwarmPopulationModel::
                            FixedPopulationNoAvalanche;
                    invalid_branch.population_limit = 0;
                    invalid_branch.spatial_histories = 0;
                    invalid_branch.spatial_length_m = 0.0;
                    invalid_branch.spatial_bins = 0;
                    invalid_branch.spatial_fit_begin_bin = 0;
                    invalid_branch.spatial_fit_end_bin = 0;
                    invalid_branch.spatial_max_steps = 0;
                    invalid_branch.spatial_min_r_squared = 0.0;
                    require_throws_contains(
                        [&] {
                            (void)pic::run_swarm_benchmark(
                                invalid_branch);
                        },
                        "does not support attachment",
                        "fixed-population swarm accepted attachment");
                } catch (...) {
                    std::filesystem::remove(swarm_table);
                    std::filesystem::remove(swarm_manifest);
                    std::filesystem::remove(swarm_config_path);
                    std::filesystem::remove(swarm_output);
                    std::filesystem::remove(swarm_ionization_table);
                    std::filesystem::remove(swarm_attachment_table);
                    std::filesystem::remove(swarm_branch_manifest);
                    std::filesystem::remove(swarm_branch_config_path);
                    std::filesystem::remove(swarm_branch_output);
                    std::filesystem::remove(swarm_spatial_output);
                    throw;
                }
                std::filesystem::remove(swarm_table);
                std::filesystem::remove(swarm_manifest);
                std::filesystem::remove(swarm_config_path);
                std::filesystem::remove(swarm_output);
                std::filesystem::remove(swarm_ionization_table);
                std::filesystem::remove(swarm_attachment_table);
                std::filesystem::remove(swarm_branch_manifest);
                std::filesystem::remove(swarm_branch_config_path);
                std::filesystem::remove(swarm_branch_output);
                std::filesystem::remove(swarm_spatial_output);
            }
            {
                const auto invalid_path =
                    std::filesystem::path(
                        "test_imported_mcc_missing_provenance.cfg");
                std::string invalid_text =
                    read_file_text(imported_mcc_example);
                const auto provenance =
                    invalid_text.find("data_provenance =");
                require(
                    provenance != std::string::npos,
                    "imported MCC fixture is missing provenance");
                invalid_text.erase(
                    provenance,
                    invalid_text.find('\n', provenance) - provenance + 1);
                {
                    std::ofstream invalid(invalid_path);
                    invalid << invalid_text;
                }
                require_throws(
                    [&]() {
                        try {
                            (void)pic::load_unstructured_config_2d(
                                invalid_path);
                        } catch (...) {
                            std::filesystem::remove(invalid_path);
                            throw;
                        }
                        std::filesystem::remove(invalid_path);
                    },
                    "imported MCC parser accepted missing data provenance");
            }

            {
                const auto output_dir =
                    std::filesystem::path(
                        "test_output_unstructured_mcc");
                const auto continued_output_dir =
                    std::filesystem::path(
                        "test_output_unstructured_mcc_continued");
                const auto changed_table =
                    std::filesystem::absolute(
                        "test_unstructured_mcc_changed.dat");
                std::filesystem::remove_all(output_dir);
                std::filesystem::remove_all(continued_output_dir);

                auto config = imported_mcc_config;
                config.output_dir = output_dir;
                pic::UnstructuredSimulation2D continuous(config);
                const auto continuous_summary = continuous.run();
                const auto& continuous_collisions =
                    continuous.collision_diagnostics();
                require(
                    continuous_summary.steps_completed == 6 &&
                        continuous_collisions.candidates == 44 &&
                        continuous_collisions.null_collisions == 20 &&
                        continuous_collisions.channel_collisions ==
                            std::vector<std::uint64_t>{24},
                    "imported 2D3V MCC deterministic collision envelope changed");
                require_near(
                    continuous_summary.final_sample.kinetic_energy,
                    31.560389541145053, 1e-12,
                    "imported finite-mass elastic energy envelope changed");
                const auto collision_csv =
                    read_file_text(output_dir / "collisions.csv");
                require(
                    collision_csv.find("cumulative_synthetic_elastic") !=
                            std::string::npos &&
                        count_lines(collision_csv) == 8,
                    "imported MCC diagnostics are incomplete");

                auto continued_config = config;
                continued_config.restart_path =
                    output_dir / "checkpoint_3.apc";
                continued_config.output_dir = continued_output_dir;
                continued_config.checkpoint_output = false;
                pic::UnstructuredSimulation2D continued(
                    continued_config);
                const auto continued_summary = continued.run();
                require_species_close(
                    continuous.species(), continued.species(),
                    "imported MCC checkpoint restart determinism");
                require(
                    continued.collision_diagnostics().candidates ==
                            continuous_collisions.candidates &&
                        continued.collision_diagnostics().null_collisions ==
                            continuous_collisions.null_collisions &&
                        continued.collision_diagnostics()
                                .channel_collisions ==
                            continuous_collisions.channel_collisions,
                    "imported MCC checkpoint did not preserve collision state");
                require_near(
                    continued_summary.final_sample.total_energy,
                    continuous_summary.final_sample.total_energy, 1e-13,
                    "imported MCC checkpoint changed total energy");

                {
                    std::ofstream changed(changed_table);
                    changed << "0 0.5\n10 0.5\n";
                }
                auto changed_config = continued_config;
                changed_config.collisions.channels.front()
                    .cross_section_file = changed_table;
                require_throws_contains(
                    [&]() {
                        pic::UnstructuredSimulation2D changed(
                            changed_config);
                        (void)changed.run();
                    },
                    "collision model mismatch",
                    "imported MCC checkpoint accepted changed cross sections");

                std::filesystem::remove_all(output_dir);
                std::filesystem::remove_all(continued_output_dir);
                std::filesystem::remove(changed_table);
            }

            {
                const auto output_dir =
                    std::filesystem::path(
                        "test_output_unstructured_ionization");
                const auto continued_output_dir =
                    std::filesystem::path(
                        "test_output_unstructured_ionization_continued");
                std::filesystem::remove_all(output_dir);
                std::filesystem::remove_all(continued_output_dir);

                auto config = imported_ionization_config;
                config.output_dir = output_dir;
                pic::UnstructuredSimulation2D continuous(config);
                const auto summary = continuous.run();
                const auto& collisions =
                    continuous.collision_diagnostics();
                require(
                    collisions.candidates == 5 &&
                        collisions.null_collisions == 2 &&
                        collisions.channel_collisions ==
                            std::vector<std::uint64_t>{3},
                    "imported ionization deterministic collision envelope changed");
                require(
                    continuous.species()[0].live_count() == 19 &&
                        continuous.species()[1].live_count() == 19 &&
                        summary.final_sample.live_particles == 38,
                    "imported ionization did not create paired products");
                const double created_charge =
                    static_cast<double>(
                        continuous.species()[0].live_count() - 16) *
                        continuous.species()[0].charge() *
                        continuous.species()[0].weight() +
                    static_cast<double>(
                        continuous.species()[1].live_count() - 16) *
                        continuous.species()[1].charge() *
                        continuous.species()[1].weight();
                require_near(
                    created_charge, 0.0, 0.0,
                    "imported ionization products do not conserve charge");
                require_near(
                    summary.final_sample.total_energy, 125.0, 1e-11,
                    "imported ionization removed the wrong threshold energy");
                const auto collision_metadata =
                    read_file_text(output_dir / "collision_data.txt");
                require(
                    collision_metadata.find(
                        "dataset_id \"aurorapic.synthetic.ionization\"") !=
                            std::string::npos &&
                        collision_metadata.find(
                            "retrieved \"2026-07-28\"") !=
                            std::string::npos &&
                        collision_metadata.find(
                            "channel \"synthetic_ionization\" "
                            "\"ionization\"") != std::string::npos,
                    "imported gas dataset metadata output is incomplete");

                auto continued_config = config;
                continued_config.restart_path =
                    output_dir / "checkpoint_2.apc";
                continued_config.output_dir =
                    continued_output_dir;
                continued_config.checkpoint_output = false;
                pic::UnstructuredSimulation2D continued(
                    continued_config);
                (void)continued.run();
                require_species_close(
                    continuous.species(), continued.species(),
                    "imported ionization checkpoint restart determinism");
                require(
                    continued.collision_diagnostics().candidates ==
                            collisions.candidates &&
                        continued.collision_diagnostics().null_collisions ==
                            collisions.null_collisions &&
                        continued.collision_diagnostics()
                                .channel_collisions ==
                        collisions.channel_collisions,
                    "imported ionization checkpoint lost collision state");

                auto changed_dataset_config = continued_config;
                changed_dataset_config.collisions.dataset_version = "2";
                require_throws_contains(
                    [&] {
                        pic::UnstructuredSimulation2D changed(
                            changed_dataset_config);
                        (void)changed.run();
                    },
                    "collision model mismatch",
                    "ionization checkpoint accepted changed gas metadata");

                auto bounded_config = config;
                bounded_config.output_dir =
                    "test_output_unstructured_ionization_bounded";
                bounded_config.checkpoint_output = false;
                bounded_config.max_particles_per_species = 16;
                pic::UnstructuredSimulation2D bounded(
                    bounded_config);
                require_throws_contains(
                    [&] { (void)bounded.run(); },
                    "max_particles_per_species",
                    "ionization ignored the product species capacity");
                require(
                    bounded.species()[0].live_count() == 16 &&
                        bounded.species()[1].live_count() == 16,
                    "ionization capacity failure created partial products");

                std::filesystem::remove_all(output_dir);
                std::filesystem::remove_all(
                    continued_output_dir);
                std::filesystem::remove_all(
                    bounded_config.output_dir);
            }

            {
                const auto output_dir =
                    std::filesystem::path(
                        "test_output_unstructured_attachment");
                std::filesystem::remove_all(output_dir);

                auto config = imported_attachment_config;
                config.output_dir = output_dir;
                pic::UnstructuredSimulation2D simulation(config);
                const auto summary = simulation.run();
                const auto& collisions =
                    simulation.collision_diagnostics();
                const std::size_t remaining_electrons =
                    simulation.species()[0].live_count();
                const std::size_t negative_ions =
                    simulation.species()[1].live_count();
                require(
                    collisions.candidates == 6 &&
                        collisions.null_collisions == 2 &&
                        collisions.channel_collisions ==
                            std::vector<std::uint64_t>{4} &&
                        remaining_electrons == 12 &&
                        negative_ions == 20 &&
                        summary.final_sample.live_particles == 32,
                    "imported attachment did not replace consumed "
                    "electrons with negative ions");
                const double final_charge =
                    -static_cast<double>(
                        remaining_electrons + negative_ions);
                require_near(
                    final_charge, -32.0, 0.0,
                    "imported attachment did not conserve charge");

                auto invalid_product = config;
                invalid_product.species[1].charge = 1.0;
                require_throws_contains(
                    [&] {
                        pic::UnstructuredSimulation2D invalid(
                            invalid_product);
                    },
                    "distinct heavier product with target charge",
                    "attachment accepted a positive-ion product");

                auto bounded_config = config;
                bounded_config.output_dir =
                    "test_output_unstructured_attachment_bounded";
                bounded_config.checkpoint_output = false;
                bounded_config.max_particles_per_species = 16;
                pic::UnstructuredSimulation2D bounded(
                    bounded_config);
                require_throws_contains(
                    [&] { (void)bounded.run(); },
                    "max_particles_per_species",
                    "attachment ignored product species capacity");
                require(
                    bounded.species()[0].live_count() == 16 &&
                        bounded.species()[1].live_count() == 16,
                    "attachment capacity failure partially consumed "
                    "electrons");

                std::filesystem::remove_all(output_dir);
                std::filesystem::remove_all(
                    bounded_config.output_dir);
            }

            {
                const auto output_dir =
                    std::filesystem::path(
                        "test_output_unstructured_charge_exchange");
                const auto continued_output_dir =
                    std::filesystem::path(
                        "test_output_unstructured_charge_exchange_continued");
                std::filesystem::remove_all(output_dir);
                std::filesystem::remove_all(continued_output_dir);

                auto config = imported_charge_exchange_config;
                config.output_dir = output_dir;
                pic::UnstructuredSimulation2D continuous(config);
                const auto summary = continuous.run();
                const auto& collisions =
                    continuous.collision_diagnostics();
                require(
                    collisions.candidates == 15 &&
                        collisions.null_collisions == 8 &&
                        collisions.channel_collisions ==
                            std::vector<std::uint64_t>{7},
                    "imported charge-exchange deterministic envelope changed");
                require(
                    summary.final_sample.live_particles == 64,
                    "charge exchange changed particle count");
                require_near(
                    summary.final_sample.total_energy, 500.0, 1e-10,
                    "charge exchange removed the wrong ion energy");
                const auto metadata =
                    read_file_text(output_dir / "collision_data.txt");
                require(
                    metadata.find(
                        "dataset_id "
                        "\"aurorapic.synthetic.charge_exchange\"") !=
                            std::string::npos &&
                        metadata.find(
                            "\"charge_exchange\"") !=
                            std::string::npos,
                    "charge-exchange metadata output is incomplete");

                auto continued_config = config;
                continued_config.restart_path =
                    output_dir / "checkpoint_3.apc";
                continued_config.output_dir =
                    continued_output_dir;
                continued_config.checkpoint_output = false;
                pic::UnstructuredSimulation2D continued(
                    continued_config);
                (void)continued.run();
                require_species_close(
                    continuous.species(), continued.species(),
                    "charge-exchange checkpoint restart determinism");
                require(
                    continued.collision_diagnostics().candidates ==
                            collisions.candidates &&
                        continued.collision_diagnostics().null_collisions ==
                            collisions.null_collisions &&
                        continued.collision_diagnostics()
                                .channel_collisions ==
                            collisions.channel_collisions,
                    "charge-exchange checkpoint lost collision state");

                std::filesystem::remove_all(output_dir);
                std::filesystem::remove_all(
                    continued_output_dir);
            }

            {
                const auto real_case_mesh =
                    std::filesystem::path(AURORA_TEST_SOURCE_DIR) /
                    "examples" / "biased_probe_2d.msh";
                const auto chamber =
                    pic::load_gmsh2_ascii_mesh2d(real_case_mesh);
                const auto quality = chamber.quality();
                const std::set<std::string> expected_labels{
                    "inlet", "outlet", "probe", "wall"};
                const auto labels = chamber.boundary_labels();
                require(
                    chamber.nodes().size() == 725 &&
                        chamber.cells().size() == 1342 &&
                        chamber.boundary_faces().size() == 108,
                    "real-case Gmsh mesh artifact changed unexpectedly");
                require(
                    std::set<std::string>(labels.begin(), labels.end()) ==
                        expected_labels,
                    "real-case Gmsh physical groups were not preserved");
                require(
                    chamber.total_area() > 0.0090 &&
                        chamber.total_area() < 0.0093,
                    "real-case chamber area is outside its geometry envelope");
                require(
                    quality.minimum_corner_angle_degrees > 30.0 &&
                        quality.maximum_cell_edge_ratio < 2.0 &&
                        quality.minimum_cell_area > 0.0,
                    "real-case mesh quality is below its acceptance envelope");
                require(
                    !chamber.locate_point({0.075, 0.0}).has_value() &&
                        chamber.locate_point({0.02, 0.0}).has_value(),
                    "real-case internal probe hole topology is wrong");
                const auto real_case_config =
                    pic::load_unstructured_config_2d(
                        std::filesystem::path(AURORA_TEST_SOURCE_DIR) /
                        "examples" / "biased_probe_2d.cfg");
                require(
                    real_case_config.emissions.size() == 1 &&
                        real_case_config.sources.size() == 2 &&
                        real_case_config.dirichlet_potentials.at("probe") ==
                            10.0 &&
                        real_case_config.particle_boundaries.at("probe") ==
                            pic::ParticleBoundary::Absorbing,
                    "real-case simulation config lost its probe physics");
            }

            const auto bad_imported_config =
                std::filesystem::path("test_output_bad_imported_config.cfg");
            {
                std::ofstream bad(bad_imported_config);
                bad << "config_version = 1\n"
                    << "dimension = 2\n"
                    << "mesh = imported\n"
                    << "mesh_file = " << fixture.string() << "\n"
                    << "unknown_runtime_key = true\n";
            }
            require_throws(
                [&]() {
                    (void)pic::load_unstructured_config_2d(
                        bad_imported_config);
                },
                "imported 2D config loader accepted an unknown key");
            std::filesystem::remove(bad_imported_config);

            const auto all_neumann_config =
                std::filesystem::path("test_output_all_neumann_imported_config.cfg");
            {
                std::ofstream bad(all_neumann_config);
                bad << "config_version = 1\n"
                    << "dimension = 2\n"
                    << "mesh = imported\n"
                    << "mesh_file = " << fixture.string() << "\n";
                for (const std::string label :
                     {"electrode", "inlet", "outlet", "wall"}) {
                    bad << "[boundary." << label << "]\n"
                        << "field = neumann\n"
                        << "normal_derivative = 0\n"
                        << "particle = reflecting\n";
                }
            }
            require_throws(
                [&]() {
                    (void)pic::load_unstructured_config_2d(
                        all_neumann_config);
                },
                "imported 2D config loader accepted an all-Neumann field gauge");
            std::filesystem::remove(all_neumann_config);

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
            {
                pic::UnstructuredMesh2D cached_mesh(mesh);
                std::vector<pic::UnstructuredParticleLocation2D> locations(
                    particles.size());
                const auto initial_cached = pic::deposit_charge_shape(
                    cached_mesh, particles, -2.0, 0.5,
                    pic::RuntimePolicy{}, locations);
                require(initial_cached.location_cache_hits == 0 &&
                            initial_cached.location_searches == 3,
                        "initial unstructured deposition did not populate location caches");
                cached_mesh.clear_charge();
                const auto repeated_cached = pic::deposit_charge_shape(
                    cached_mesh, particles, -2.0, 0.5,
                    pic::RuntimePolicy{}, locations);
                require(repeated_cached.location_cache_hits == 2 &&
                            repeated_cached.location_searches == 1,
                        "repeated unstructured deposition did not reuse valid locations");
                particles[0].position = {0.8, 0.25};
                cached_mesh.clear_charge();
                const auto crossed_cached = pic::deposit_charge_shape(
                    cached_mesh, particles, -2.0, 0.5,
                    pic::RuntimePolicy{}, locations);
                require(crossed_cached.location_cache_hits == 1 &&
                            crossed_cached.location_searches == 2,
                        "unstructured location cache did not fall back after a cell crossing");
                require_throws(
                    [&]() {
                        std::vector<pic::UnstructuredParticleLocation2D> wrong_size;
                        (void)pic::deposit_charge_shape(
                            cached_mesh, particles, -2.0, 0.5,
                            pic::RuntimePolicy{}, wrong_size);
                    },
                    "unstructured deposition accepted a mismatched location cache");
                particles[0].position = {1.0 / 6.0, 1.0 / 3.0};
            }
#ifdef AURORA_HAVE_OPENMP
            {
                std::vector<pic::Particle2D> scaling_particles(128);
                for (std::size_t i = 0; i < scaling_particles.size(); ++i) {
                    scaling_particles[i].position =
                        i % 2 == 0 ? pic::Vec2{0.2, 0.2} : pic::Vec2{0.8, 0.25};
                }
                pic::UnstructuredMesh2D serial_deposit_mesh(mesh);
                pic::UnstructuredMesh2D parallel_deposit_mesh(mesh);
                const auto serial_deposit = pic::deposit_charge_shape(
                    serial_deposit_mesh, scaling_particles, -1.0, 0.125,
                    pic::RuntimePolicy{});
                const auto parallel_deposit = pic::deposit_charge_shape(
                    parallel_deposit_mesh, scaling_particles, -1.0, 0.125,
                    pic::RuntimePolicy{pic::RuntimeBackend::OpenMP, 2});
                require(serial_deposit.deposited_particles ==
                            parallel_deposit.deposited_particles &&
                            serial_deposit.outside_particles ==
                                parallel_deposit.outside_particles,
                        "parallel unstructured deposition changed particle accounting");
                require_near(serial_deposit.deposited_charge,
                             parallel_deposit.deposited_charge, 1e-14,
                             "parallel unstructured deposition changed total charge");
                for (std::size_t i = 0; i < serial_deposit_mesh.size(); ++i) {
                    require_near(serial_deposit_mesh.rho()[i],
                                 parallel_deposit_mesh.rho()[i], 1e-12,
                                 "parallel unstructured deposition changed nodal charge");
                }
            }
#endif
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
            pic::UnstructuredParticleLocation2D interpolation_location;
            bool interpolation_hit = true;
            require(
                pic::interpolate_electric(
                    computational_mesh, {0.2, 0.2},
                    interpolation_location, &interpolation_hit).has_value() &&
                    !interpolation_hit,
                "initial unstructured interpolation incorrectly reported a cache hit");
            require(
                pic::interpolate_electric(
                    computational_mesh, {0.21, 0.2},
                    interpolation_location, &interpolation_hit).has_value() &&
                    interpolation_hit,
                "same-cell unstructured interpolation missed its cached location");
            require(
                pic::interpolate_electric(
                    computational_mesh, {0.8, 0.25},
                    interpolation_location, &interpolation_hit).has_value() &&
                    !interpolation_hit,
                "cross-cell unstructured interpolation did not fall back to point location");
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
                pic::UnstructuredPoissonSolver2D reusable_solver(
                    field_mesh, {{"ground", 0.0}});
                const auto source_summary = reusable_solver.solve(field_mesh);
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
                const auto reference_phi = field_mesh.phi();
                const auto reference_electric = field_mesh.electric();
                std::fill(field_mesh.phi().begin(), field_mesh.phi().end(), 3.0);
                const auto repeated_summary = reusable_solver.solve(field_mesh);
                require(repeated_summary.converged,
                        "reused unstructured Poisson solve did not converge");
                require(reusable_solver.assembly_count() == 1 &&
                            reusable_solver.solve_count() == 2,
                        "unstructured Poisson solver rebuilt or miscounted its cached operator");
                for (std::size_t i = 0; i < field_mesh.size(); ++i) {
                    require_near(field_mesh.phi()[i], reference_phi[i], 1e-12,
                                 "reused unstructured Poisson potential changed");
                    require_near(field_mesh.electric()[i].x, reference_electric[i].x, 1e-12,
                                 "reused unstructured Poisson Ex changed");
                    require_near(field_mesh.electric()[i].y, reference_electric[i].y, 1e-12,
                                 "reused unstructured Poisson Ey changed");
                }
                pic::UnstructuredMesh2D distinct_field_mesh(star);
                require_throws(
                    [&]() { (void)reusable_solver.solve(distinct_field_mesh); },
                    "unstructured Poisson solver accepted a different mesh topology");

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
                // Mixed weak-form boundary benchmark: phi=x satisfies Laplace's
                // equation with phi=0 at x=0, dphi/dn=1 at x=1, and homogeneous
                // Neumann data on the horizontal sides.
                pic::ImportedMesh2D mixed_boundary_patch;
                mixed_boundary_patch.add_node({1, {0.0, 0.0}});
                mixed_boundary_patch.add_node({2, {1.0, 0.0}});
                mixed_boundary_patch.add_node({3, {1.0, 1.0}});
                mixed_boundary_patch.add_node({4, {0.0, 1.0}});
                mixed_boundary_patch.add_cell({
                    10, pic::ImportedCellShape2D::Quadrilateral,
                    {1, 2, 3, 4}, 1, "plasma"});
                mixed_boundary_patch.add_boundary_face({20, {1, 2}, 1, "bottom"});
                mixed_boundary_patch.add_boundary_face({21, {2, 3}, 1, "right"});
                mixed_boundary_patch.add_boundary_face({22, {3, 4}, 1, "top"});
                mixed_boundary_patch.add_boundary_face({23, {4, 1}, 1, "left"});

                pic::UnstructuredMesh2D mixed_boundary_mesh(mixed_boundary_patch);
                const auto summary = pic::solve_unstructured_poisson(
                    mixed_boundary_mesh,
                    {{"left", 0.0}},
                    {{"bottom", 0.0}, {"right", 1.0}, {"top", 0.0}});
                require(summary.converged,
                        "mixed Dirichlet/Neumann Poisson solve did not converge");
                for (const auto& node : mixed_boundary_patch.nodes()) {
                    const std::size_t index = mixed_boundary_mesh.node_index(node.id);
                    require_near(mixed_boundary_mesh.phi()[index], node.position.x, 1e-12,
                                 "mixed-boundary linear potential mismatch");
                    require_near(mixed_boundary_mesh.electric()[index].x, -1.0, 1e-12,
                                 "mixed-boundary linear Ex mismatch");
                    require_near(mixed_boundary_mesh.electric()[index].y, 0.0, 1e-12,
                                 "mixed-boundary linear Ey mismatch");
                }
                require_throws(
                    [&]() {
                        (void)pic::solve_unstructured_poisson(
                            mixed_boundary_mesh, {},
                            {{"bottom", 0.0}, {"right", 1.0},
                             {"top", 0.0}, {"left", -1.0}});
                    },
                    "unstructured Poisson solve accepted an all-Neumann gauge");
                require_throws(
                    [&]() {
                        (void)pic::solve_unstructured_poisson(
                            mixed_boundary_mesh, {{"left", 0.0}},
                            {{"left", 0.0}, {"bottom", 0.0},
                             {"right", 1.0}, {"top", 0.0}});
                    },
                    "unstructured Poisson solve accepted overlapping field conditions");
                require_throws(
                    [&]() {
                        (void)pic::solve_unstructured_poisson(
                            mixed_boundary_mesh, {{"left", 0.0}},
                            {{"bottom", 0.0}, {"right", 1.0}});
                    },
                    "unstructured Poisson solve accepted a missing field condition");
                require_throws(
                    [&]() {
                        (void)pic::solve_unstructured_poisson(
                            mixed_boundary_mesh, {{"left", 0.0}},
                            {{"bottom", 0.0},
                             {"right", std::numeric_limits<double>::infinity()},
                             {"top", 0.0}});
                    },
                    "unstructured Poisson solve accepted non-finite Neumann data");
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

            {
                const auto output_dir =
                    std::filesystem::path("test_output_unstructured_reflecting");
                std::filesystem::remove_all(output_dir);
                pic::UnstructuredSimulation2DConfig config;
                config.mesh_path = fixture;
                config.dt = 1.0;
                config.steps = 1;
                config.output_interval = 1;
                config.output_dir = output_dir;
                config.vtk_output = true;
                config.particle_output = true;
                config.particle_output_interval = 1;
                config.particle_output_stride = 2;
                config.particle_sample_count = 4;
                config.checkpoint_output = true;
                config.checkpoint_interval = 1;
                config.dirichlet_potentials = {
                    {"electrode", 0.0}, {"inlet", 0.0},
                    {"outlet", 0.0}, {"wall", 0.0},
                };
                config.particle_boundaries = {
                    {"electrode", pic::ParticleBoundary::Reflecting},
                    {"inlet", pic::ParticleBoundary::Reflecting},
                    {"outlet", pic::ParticleBoundary::Reflecting},
                    {"wall", pic::ParticleBoundary::Reflecting},
                };
                pic::UnstructuredSpecies2DConfig neutral;
                neutral.name = "neutral";
                neutral.charge = 0.0;
                neutral.mass = 1.0;
                neutral.weight = 1.0;
                neutral.particles = 16;
                neutral.drift_velocity_x = 2.0;
                neutral.drift_velocity_z = 0.7;
                neutral.thermal_velocity = 0.0;
                config.species = {neutral};

                pic::UnstructuredSimulation2D simulation(config);
                const auto summary = simulation.run();
                require(summary.steps_completed == 1 && summary.final_time == config.dt,
                        "unstructured runtime did not complete its configured transient step");
                require(simulation.poisson_assembly_count() == 1 &&
                            simulation.poisson_solve_count() == 2,
                        "unstructured runtime did not reuse its Poisson operator");
                require(simulation.timing().location_cache_hits > 0 &&
                            simulation.timing().location_searches >= neutral.particles,
                        "unstructured runtime did not report location-cache reuse");
                require(summary.final_sample.live_particles == neutral.particles,
                        "reflecting imported boundaries lost particles");
                for (const auto& [label, count] : summary.final_sample.absorbed_by_label) {
                    (void)label;
                    require(count == 0, "reflecting imported boundary reported absorbed particles");
                }
                for (const auto& particle : simulation.species().front().particles()) {
                    require(particle.alive && simulation.mesh().locate_point(particle.position).has_value(),
                            "reflected particle did not remain inside imported geometry");
                    require_near(particle.velocity_z, neutral.drift_velocity_z,
                                 1e-14,
                                 "imported reflection changed out-of-plane velocity");
                }
                require(std::filesystem::exists(output_dir / "scalars.csv"),
                        "unstructured runtime did not write scalar diagnostics");
                require(std::filesystem::exists(output_dir / "fields_0.vtu") &&
                            std::filesystem::exists(output_dir / "fields_1.vtu"),
                        "unstructured runtime did not write requested VTU snapshots");
                require(std::filesystem::exists(output_dir / "particles_1.csv"),
                        "unstructured runtime did not write requested particle samples");
                require(std::filesystem::exists(output_dir / "checkpoint_1.apc"),
                        "unstructured runtime did not write requested checkpoint");
                const std::string vtu = read_file_text(output_dir / "fields_1.vtu");
                require(vtu.find("type=\"UnstructuredGrid\"") != std::string::npos,
                        "unstructured VTU output has the wrong dataset type");
                require(vtu.find("Name=\"connectivity\"") != std::string::npos &&
                            vtu.find("Name=\"types\"") != std::string::npos,
                        "unstructured VTU output is missing cell topology arrays");
                const std::string scalars = read_file_text(output_dir / "scalars.csv");
                require(scalars.find("poisson_final_residual") != std::string::npos &&
                            scalars.find("particle_seconds") != std::string::npos &&
                            scalars.find("deposition_seconds") != std::string::npos &&
                            scalars.find("field_solve_seconds") != std::string::npos &&
                            count_lines(scalars) == 3,
                        "unstructured scalar diagnostics structure mismatch");

#ifdef AURORA_HAVE_OPENMP
                auto parallel_config = config;
                parallel_config.runtime = {
                    pic::RuntimeBackend::OpenMP, 2};
                parallel_config.output_dir =
                    "test_output_unstructured_parallel";
                parallel_config.vtk_output = false;
                parallel_config.particle_output = false;
                parallel_config.checkpoint_output = false;
                std::filesystem::remove_all(parallel_config.output_dir);
                pic::UnstructuredSimulation2D parallel_simulation(parallel_config);
                const auto parallel_summary = parallel_simulation.run();
                require_near(
                    parallel_summary.final_sample.total_energy,
                    summary.final_sample.total_energy, 1e-13,
                    "OpenMP imported runtime changed total energy");
                require(parallel_summary.final_sample.absorbed_by_label ==
                            summary.final_sample.absorbed_by_label,
                        "OpenMP imported runtime changed boundary loss accounting");
                const auto& parallel_particles =
                    parallel_simulation.species().front().particles();
                const auto& serial_particles =
                    simulation.species().front().particles();
                require(parallel_particles.size() == serial_particles.size(),
                        "OpenMP imported runtime changed particle count");
                for (std::size_t i = 0; i < serial_particles.size(); ++i) {
                    require_near(parallel_particles[i].position.x,
                                 serial_particles[i].position.x, 1e-13,
                                 "OpenMP imported runtime changed particle x");
                    require_near(parallel_particles[i].position.y,
                                 serial_particles[i].position.y, 1e-13,
                                 "OpenMP imported runtime changed particle y");
                    require_near(parallel_particles[i].velocity_half.x,
                                 serial_particles[i].velocity_half.x, 1e-13,
                                 "OpenMP imported runtime changed half-step vx");
                    require_near(parallel_particles[i].velocity_half.y,
                                 serial_particles[i].velocity_half.y, 1e-13,
                                 "OpenMP imported runtime changed half-step vy");
                }
                std::filesystem::remove_all(parallel_config.output_dir);
#endif

                auto continued_config = config;
                continued_config.steps = 2;
                continued_config.restart_path = output_dir / "checkpoint_1.apc";
                continued_config.output_dir =
                    "test_output_unstructured_restart_continued";
                continued_config.vtk_output = false;
                continued_config.particle_output = false;
                continued_config.checkpoint_output = false;
                std::filesystem::remove_all(continued_config.output_dir);
                pic::UnstructuredSimulation2D continued(continued_config);
                const auto continued_summary = continued.run();

                auto reference_config = config;
                reference_config.steps = 2;
                reference_config.output_dir =
                    "test_output_unstructured_restart_reference";
                reference_config.vtk_output = false;
                reference_config.particle_output = false;
                reference_config.checkpoint_output = false;
                std::filesystem::remove_all(reference_config.output_dir);
                pic::UnstructuredSimulation2D reference(reference_config);
                const auto reference_summary = reference.run();
                require_near(continued_summary.final_sample.total_energy,
                             reference_summary.final_sample.total_energy, 1e-13,
                             "unstructured checkpoint restart changed total energy");
                const auto& continued_particles =
                    continued.species().front().particles();
                const auto& reference_particles =
                    reference.species().front().particles();
                require(continued_particles.size() == reference_particles.size(),
                        "unstructured checkpoint restart changed particle count");
                for (std::size_t i = 0; i < continued_particles.size(); ++i) {
                    require_near(continued_particles[i].position.x,
                                 reference_particles[i].position.x, 1e-13,
                                 "unstructured checkpoint restart changed particle x");
                    require_near(continued_particles[i].position.y,
                                 reference_particles[i].position.y, 1e-13,
                                 "unstructured checkpoint restart changed particle y");
                    require_near(continued_particles[i].velocity_half.x,
                                 reference_particles[i].velocity_half.x, 1e-13,
                                 "unstructured checkpoint restart changed half-step vx");
                    require_near(continued_particles[i].velocity_half.y,
                                 reference_particles[i].velocity_half.y, 1e-13,
                                 "unstructured checkpoint restart changed half-step vy");
                }
                std::filesystem::remove_all(output_dir);
                std::filesystem::remove_all(continued_config.output_dir);
                std::filesystem::remove_all(reference_config.output_dir);
            }

            {
                const auto output_dir =
                    std::filesystem::path("test_output_unstructured_absorbing");
                std::filesystem::remove_all(output_dir);
                pic::UnstructuredSimulation2DConfig config;
                config.mesh_path = fixture;
                config.dt = 1.0;
                config.steps = 1;
                config.output_interval = 1;
                config.output_dir = output_dir;
                config.dirichlet_potentials = {
                    {"electrode", 0.0}, {"inlet", 0.0},
                    {"outlet", 0.0}, {"wall", 0.0},
                };
                config.particle_boundaries = {
                    {"electrode", pic::ParticleBoundary::Absorbing},
                    {"inlet", pic::ParticleBoundary::Absorbing},
                    {"outlet", pic::ParticleBoundary::Absorbing},
                    {"wall", pic::ParticleBoundary::Absorbing},
                };
                pic::UnstructuredSpecies2DConfig beam;
                beam.name = "beam";
                beam.charge = 0.0;
                beam.particles = 16;
                beam.drift_velocity_x = 2.0;
                beam.thermal_velocity = 0.0;
                config.species = {beam};

                pic::UnstructuredSimulation2D simulation(config);
                const auto summary = simulation.run();
                require(summary.final_sample.live_particles == 0,
                        "absorbing imported boundary did not remove outgoing particles");
                require(summary.final_sample.absorbed_by_label.at("outlet") == beam.particles,
                        "absorbing imported boundary attributed losses to the wrong physical label");
                std::filesystem::remove_all(output_dir);
            }

            {
                const auto output_dir =
                    std::filesystem::path("test_output_unstructured_source");
                std::filesystem::remove_all(output_dir);
                pic::UnstructuredSimulation2DConfig config;
                config.mesh_path = fixture;
                config.dt = 0.1;
                config.steps = 3;
                config.output_interval = 1;
                config.output_dir = output_dir;
                config.checkpoint_output = true;
                config.checkpoint_interval = 1;
                config.max_particles_per_species = 32;
                config.dirichlet_potentials = {
                    {"electrode", 0.0}, {"inlet", 0.0},
                    {"outlet", 0.0}, {"wall", 0.0},
                };
                config.particle_boundaries = {
                    {"electrode", pic::ParticleBoundary::Reflecting},
                    {"inlet", pic::ParticleBoundary::Reflecting},
                    {"outlet", pic::ParticleBoundary::Reflecting},
                    {"wall", pic::ParticleBoundary::Reflecting},
                };
                pic::UnstructuredSpecies2DConfig tracer;
                tracer.name = "tracer";
                tracer.charge = 0.0;
                tracer.particles = 0;
                tracer.thermal_velocity = 0.0;
                config.species = {tracer};
                pic::UnstructuredBoundarySource2DConfig source;
                source.name = "tracer_inlet";
                source.species = "tracer";
                source.boundary = "inlet";
                source.particles_per_step = 3;
                source.start_step = 0;
                source.end_step = 2;
                source.normal_velocity = 0.0;
                source.thermal_velocity = 0.0;
                source.out_of_plane_velocity = 0.6;
                config.sources = {source};

                pic::UnstructuredSimulation2D simulation(config);
                const auto summary = simulation.run();
                require(summary.final_sample.live_particles == 6 &&
                            simulation.species().front().particles().size() == 6,
                        "imported boundary source injected the wrong particle count");
                require(summary.final_sample.injected_by_source.at("tracer_inlet") == 6,
                        "imported boundary source diagnostic count mismatch");
                for (const auto& particle :
                     simulation.species().front().particles()) {
                    require(particle.alive &&
                                simulation.mesh().locate_point(
                                    particle.position).has_value(),
                            "imported boundary source generated an exterior particle");
                    require_near(
                        particle.velocity_z, source.out_of_plane_velocity,
                        1e-14,
                        "imported boundary source did not initialize out-of-plane velocity");
                }
                const std::string scalars =
                    read_file_text(output_dir / "scalars.csv");
                require(scalars.find("injected_tracer_inlet") !=
                            std::string::npos,
                        "imported boundary source diagnostic column is missing");

                auto continued_config = config;
                continued_config.restart_path =
                    output_dir / "checkpoint_1.apc";
                continued_config.output_dir =
                    "test_output_unstructured_source_continued";
                continued_config.checkpoint_output = false;
                std::filesystem::remove_all(continued_config.output_dir);
                pic::UnstructuredSimulation2D continued(continued_config);
                const auto continued_summary = continued.run();
                require(
                    continued_summary.final_sample.injected_by_source.at(
                        "tracer_inlet") == 6,
                    "source checkpoint restart changed cumulative injection");
                const auto& expected_particles =
                    simulation.species().front().particles();
                const auto& continued_particles =
                    continued.species().front().particles();
                require(continued_particles.size() == expected_particles.size(),
                        "source checkpoint restart changed particle storage");
                for (std::size_t i = 0; i < expected_particles.size(); ++i) {
                    require_near(
                        continued_particles[i].position.x,
                        expected_particles[i].position.x, 1e-13,
                        "source checkpoint restart changed particle x");
                    require_near(
                        continued_particles[i].position.y,
                        expected_particles[i].position.y, 1e-13,
                        "source checkpoint restart changed particle y");
                    require_near(
                        continued_particles[i].velocity_z,
                        expected_particles[i].velocity_z, 1e-13,
                        "source checkpoint restart changed particle vz");
                    require_near(
                        continued_particles[i].velocity_half_z,
                        expected_particles[i].velocity_half_z, 1e-13,
                        "source checkpoint restart changed half-step particle vz");
                }
                auto mismatched_units = continued_config;
                mismatched_units.units.system = pic::UnitSystem::SI;
                mismatched_units.output_dir =
                    "test_output_unstructured_source_bad_units";
                std::filesystem::remove_all(mismatched_units.output_dir);
                require_throws(
                    [&]() {
                        pic::UnstructuredSimulation2D mismatched(
                            mismatched_units);
                        (void)mismatched.run();
                    },
                    "unstructured checkpoint accepted a different unit system");

                auto bounded_config = config;
                bounded_config.max_particles_per_species = 2;
                bounded_config.steps = 1;
                bounded_config.checkpoint_output = false;
                bounded_config.output_dir =
                    "test_output_unstructured_source_bounded";
                std::filesystem::remove_all(bounded_config.output_dir);
                require_throws(
                    [&]() {
                        pic::UnstructuredSimulation2D bounded(bounded_config);
                        (void)bounded.run();
                    },
                    "imported boundary source exceeded its particle storage bound");

                auto recycling_config = config;
                recycling_config.steps = 2;
                recycling_config.max_particles_per_species = 3;
                recycling_config.checkpoint_output = false;
                recycling_config.output_dir =
                    "test_output_unstructured_source_recycling";
                recycling_config.sources.front().end_step = 0;
                recycling_config.sources.front().normal_velocity = 20.0;
                for (auto& [label, policy] :
                     recycling_config.particle_boundaries) {
                    (void)label;
                    policy = pic::ParticleBoundary::Absorbing;
                }
                std::filesystem::remove_all(recycling_config.output_dir);
                pic::UnstructuredSimulation2D recycling(recycling_config);
                const auto recycling_summary = recycling.run();
                require(
                    recycling_summary.final_sample.injected_by_source.at(
                        "tracer_inlet") == 6 &&
                        recycling.species().front().particles().size() == 3 &&
                        recycling_summary.final_sample.live_particles == 0,
                    "imported boundary source did not reuse absorbed particle slots");
                std::size_t recycled_absorbed = 0;
                for (const auto& [label, count] :
                     recycling_summary.final_sample.absorbed_by_label) {
                    (void)label;
                    recycled_absorbed += count;
                }
                require(recycled_absorbed == 6,
                        "injected particles were not absorbed after transit");

                std::filesystem::remove_all(output_dir);
                std::filesystem::remove_all(continued_config.output_dir);
                std::filesystem::remove_all(mismatched_units.output_dir);
                std::filesystem::remove_all(bounded_config.output_dir);
                std::filesystem::remove_all(recycling_config.output_dir);
            }

            {
                const auto output_dir =
                    std::filesystem::path("test_output_unstructured_emission");
                std::filesystem::remove_all(output_dir);
                pic::UnstructuredSimulation2DConfig config;
                config.mesh_path = fixture;
                config.dt = 1.0;
                config.steps = 2;
                config.output_interval = 1;
                config.output_dir = output_dir;
                config.checkpoint_output = true;
                config.checkpoint_interval = 1;
                config.max_particles_per_species = 64;
                config.dirichlet_potentials = {
                    {"electrode", 0.0}, {"inlet", 0.0},
                    {"outlet", 0.0}, {"wall", 0.0},
                };
                config.particle_boundaries = {
                    {"electrode", pic::ParticleBoundary::Absorbing},
                    {"inlet", pic::ParticleBoundary::Absorbing},
                    {"outlet", pic::ParticleBoundary::Absorbing},
                    {"wall", pic::ParticleBoundary::Absorbing},
                };
                pic::UnstructuredSpecies2DConfig beam;
                beam.name = "beam";
                beam.charge = 1.0;
                beam.mass = 2.0;
                beam.weight = 2.0;
                beam.particles = 8;
                beam.drift_velocity_x = 2.0;
                beam.thermal_velocity = 0.0;
                pic::UnstructuredSpecies2DConfig secondary;
                secondary.name = "secondary";
                secondary.charge = -1.0;
                secondary.mass = 1.0;
                secondary.weight = 1.0;
                secondary.particles = 0;
                secondary.thermal_velocity = 0.0;
                config.species = {beam, secondary};
                pic::UnstructuredSecondaryEmission2DConfig emission;
                emission.name = "wall_secondaries";
                emission.boundary = "outlet";
                emission.incident_species = "beam";
                emission.emitted_species = "secondary";
                emission.yield = 1.0;
                emission.max_particles_per_impact = 2;
                emission.normal_velocity = 0.25;
                emission.out_of_plane_velocity = 0.4;
                config.emissions = {emission};

                pic::UnstructuredSimulation2D simulation(config);
                const auto summary = simulation.run();
                require(
                    summary.final_sample.absorbed_by_label.at("outlet") == 8 &&
                        summary.final_sample.emitted_by_rule.at(
                            "wall_secondaries") == 16 &&
                        simulation.species()[1].live_count() == 16,
                    "secondary emission count or weight conversion is wrong");
                const auto& flux =
                    summary.final_sample.impact_flux.at("beam").at("outlet");
                require(
                    flux.macroparticles == 8 &&
                        flux.physical_particles == 16.0 &&
                        flux.charge == 16.0 &&
                        flux.kinetic_energy == 64.0,
                    "species-resolved boundary flux accumulation is wrong");
                require(
                    flux.last_step_macroparticles == 0 &&
                        flux.physical_particle_rate == 0.0 &&
                        flux.physical_particle_flux == 0.0,
                    "boundary flux last-step rate did not reset");
                for (const auto& particle :
                     simulation.species()[1].particles()) {
                    require_near(
                        particle.velocity_z,
                        emission.out_of_plane_velocity, 1e-14,
                        "secondary emission did not initialize out-of-plane velocity");
                }
                const std::string scalars =
                    read_file_text(output_dir / "scalars.csv");
                require(
                    scalars.find("emitted_wall_secondaries") !=
                            std::string::npos &&
                        scalars.find("impact_flux_beam@outlet") !=
                            std::string::npos,
                    "secondary emission or boundary flux diagnostics are missing");

                auto continued_config = config;
                continued_config.restart_path =
                    output_dir / "checkpoint_1.apc";
                continued_config.output_dir =
                    "test_output_unstructured_emission_continued";
                continued_config.checkpoint_output = false;
                std::filesystem::remove_all(continued_config.output_dir);
                pic::UnstructuredSimulation2D continued(continued_config);
                const auto continued_summary = continued.run();
                require(
                    continued_summary.final_sample.emitted_by_rule ==
                            summary.final_sample.emitted_by_rule &&
                        continued_summary.final_sample.impact_flux.at(
                            "beam").at("outlet").physical_particles ==
                            flux.physical_particles,
                    "emission checkpoint restart changed accumulated diagnostics");
                const auto& expected_secondaries =
                    simulation.species()[1].particles();
                const auto& continued_secondaries =
                    continued.species()[1].particles();
                require(
                    continued_secondaries.size() ==
                        expected_secondaries.size(),
                    "emission checkpoint restart changed particle storage");
                for (std::size_t i = 0;
                     i < expected_secondaries.size(); ++i) {
                    require_near(
                        continued_secondaries[i].position.x,
                        expected_secondaries[i].position.x, 1e-13,
                        "emission checkpoint restart changed particle x");
                    require_near(
                        continued_secondaries[i].position.y,
                        expected_secondaries[i].position.y, 1e-13,
                        "emission checkpoint restart changed particle y");
                    require_near(
                        continued_secondaries[i].velocity_z,
                        expected_secondaries[i].velocity_z, 1e-13,
                        "emission checkpoint restart changed particle vz");
                    require_near(
                        continued_secondaries[i].velocity_half_z,
                        expected_secondaries[i].velocity_half_z, 1e-13,
                        "emission checkpoint restart changed half-step particle vz");
                }

#ifdef AURORA_HAS_OPENMP
                auto parallel_config = config;
                parallel_config.runtime.backend =
                    pic::RuntimeBackend::OpenMP;
                parallel_config.runtime.threads = 2;
                parallel_config.output_dir =
                    "test_output_unstructured_emission_parallel";
                parallel_config.checkpoint_output = false;
                std::filesystem::remove_all(parallel_config.output_dir);
                pic::UnstructuredSimulation2D parallel(parallel_config);
                const auto parallel_summary = parallel.run();
                require(
                    parallel_summary.final_sample.emitted_by_rule ==
                            summary.final_sample.emitted_by_rule &&
                        parallel_summary.final_sample.impact_flux.at(
                            "beam").at("outlet").physical_particles ==
                            flux.physical_particles,
                    "parallel emission changed deterministic accounting");
                const auto& parallel_secondaries =
                    parallel.species()[1].particles();
                require(
                    parallel_secondaries.size() ==
                        expected_secondaries.size(),
                    "parallel emission changed particle storage");
                for (std::size_t i = 0;
                     i < expected_secondaries.size(); ++i) {
                    require_near(
                        parallel_secondaries[i].position.x,
                        expected_secondaries[i].position.x, 1e-13,
                        "parallel emission changed particle x");
                    require_near(
                        parallel_secondaries[i].position.y,
                        expected_secondaries[i].position.y, 1e-13,
                        "parallel emission changed particle y");
                }
                std::filesystem::remove_all(parallel_config.output_dir);
#endif

                auto invalid_config = config;
                invalid_config.emissions.front().max_particles_per_impact = 1;
                require_throws(
                    [&]() {
                        pic::UnstructuredSimulation2D invalid(invalid_config);
                    },
                    "emission accepted an unsafe macro-particle yield");
                std::filesystem::remove_all(output_dir);
                std::filesystem::remove_all(continued_config.output_dir);
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
