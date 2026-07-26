#include "pic/Config.hpp"
#include "pic/FieldSolver.hpp"
#include "pic/Grid.hpp"
#include "pic/Mesh2D.hpp"
#include "pic/Simulation.hpp"
#include "pic/Simulation2D.hpp"
#include "pic/VTKWriter.hpp"
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
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "aurorapic core test failure: " << e.what() << '\n';
        return 1;
    }
}
