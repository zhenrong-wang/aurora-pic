#include "pic/Config.hpp"
#include "pic/Grid.hpp"
#include "pic/Initialization.hpp"
#include "pic/Mesh2D.hpp"
#include "pic/Mesh3D.hpp"
#include "pic/Species.hpp"
#include "pic/Species2D.hpp"
#include "pic/Species3D.hpp"
#include "pic/Simulation2D.hpp"
#include "pic/UnstructuredSimulation2D.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void require_near(
    double actual, double expected, double tolerance,
    const std::string& message) {
    if (std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(
            message + ": expected " + std::to_string(expected) +
            ", got " + std::to_string(actual));
    }
}

template <typename Function>
void require_throws(Function&& function, const std::string& message) {
    try {
        function();
    } catch (const std::exception&) {
        return;
    }
    throw std::runtime_error(message);
}

double mean(const std::vector<double>& values) {
    return std::accumulate(values.begin(), values.end(), 0.0) /
           static_cast<double>(values.size());
}

double variance(const std::vector<double>& values, double center) {
    double sum = 0.0;
    for (const double value : values) {
        const double delta = value - center;
        sum += delta * delta;
    }
    return sum / static_cast<double>(values.size());
}

} // namespace

int main() {
    try {
        require(
            pic::particle_loading_from_string("quiet_start") ==
                pic::ParticleLoading::QuietStart,
            "quiet_start loading parser failed");
        require_throws(
            [] { (void)pic::particle_loading_from_string("silent"); },
            "unknown loading model was accepted");

        {
            pic::SpeciesConfig config;
            config.particles = 8;
            config.drift_velocity = 1.25;
            config.thermal_velocity = 0.4;
            config.initialization.loading =
                pic::ParticleLoading::QuietStart;
            pic::Species species(config);
            pic::Grid grid(9, 1.0, pic::Boundary::Dirichlet);
            std::mt19937_64 rng(1234);
            species.initialize(grid, rng);

            std::vector<double> positions;
            std::vector<double> velocities;
            for (const auto& particle : species.particles()) {
                positions.push_back(particle.x);
                velocities.push_back(particle.v);
            }
            std::sort(positions.begin(), positions.end());
            for (std::size_t index = 0; index < positions.size(); ++index) {
                require_near(
                    positions[index],
                    (static_cast<double>(index) + 0.5) /
                        static_cast<double>(positions.size()),
                    1e-15,
                    "1D quiet-start position is not centered in its stratum");
            }
            require_near(
                mean(velocities), config.drift_velocity, 1e-15,
                "1D antithetic loading did not preserve mean drift");
            require(
                variance(velocities, config.drift_velocity) > 0.0,
                "1D quiet-start thermal spread vanished");
        }

        {
            pic::Species2DConfig config;
            config.particles = 10;
            config.drift_velocity_x = 0.25;
            config.drift_velocity_y = -0.5;
            config.drift_velocity_z = 0.75;
            config.thermal_velocity = 9.0;
            config.initialization.loading =
                pic::ParticleLoading::QuietStart;
            config.initialization.thermal_velocity_x = 0.2;
            config.initialization.thermal_velocity_y = 0.4;
            config.initialization.thermal_velocity_z = 0.0;
            pic::Species2D species(config);
            pic::Mesh2D mesh(
                7, 7, 1.0, 1.0, pic::Boundary::Dirichlet);
            std::mt19937_64 rng(5678);
            species.initialize(mesh, rng);

            std::vector<double> x;
            std::vector<double> y;
            std::vector<double> z;
            for (const auto& particle : species.particles()) {
                require(
                    particle.position.x > 0.0 &&
                        particle.position.x < 1.0 &&
                        particle.position.y > 0.0 &&
                        particle.position.y < 1.0,
                    "2D quiet-start particle lies on or outside a wall");
                x.push_back(particle.velocity.x);
                y.push_back(particle.velocity.y);
                z.push_back(particle.velocity_z);
            }
            require_near(
                mean(x), config.drift_velocity_x, 1e-15,
                "2D quiet-start x drift is wrong");
            require_near(
                mean(y), config.drift_velocity_y, 1e-15,
                "2D quiet-start y drift is wrong");
            require_near(
                mean(z), config.drift_velocity_z, 1e-15,
                "2D quiet-start z drift is wrong");
            require(
                variance(x, config.drift_velocity_x) > 0.0 &&
                    variance(y, config.drift_velocity_y) > 0.0,
                "2D anisotropic thermal spread vanished");
            require_near(
                variance(z, config.drift_velocity_z), 0.0, 1e-30,
                "2D per-axis zero thermal spread was not honored");
        }

        {
            pic::Species3DConfig config;
            config.particles = 6;
            config.initialization.loading =
                pic::ParticleLoading::QuietStart;
            config.initialization.thermal_velocity_x = 0.1;
            config.initialization.thermal_velocity_y = 0.0;
            config.initialization.thermal_velocity_z = 0.3;
            pic::Species3D species(config);
            pic::Mesh3D mesh(
                5, 5, 5, 1.0, 1.0, 1.0,
                pic::Boundary::Dirichlet);
            std::mt19937_64 rng(9012);
            species.initialize(mesh, rng);

            std::vector<double> x;
            std::vector<double> y;
            std::vector<double> z;
            for (const auto& particle : species.particles()) {
                x.push_back(particle.velocity.x);
                y.push_back(particle.velocity.y);
                z.push_back(particle.velocity.z);
            }
            require_near(mean(x), 0.0, 1e-15, "3D x drift is wrong");
            require_near(mean(y), 0.0, 1e-15, "3D y drift is wrong");
            require_near(mean(z), 0.0, 1e-15, "3D z drift is wrong");
            require(
                variance(x, 0.0) > 0.0 &&
                    variance(z, 0.0) > 0.0,
                "3D anisotropic thermal spread vanished");
            require_near(
                variance(y, 0.0), 0.0, 1e-30,
                "3D zero y thermal spread was not honored");
        }

        {
            const auto config_path =
                std::filesystem::path("test_initialization_config.cfg");
            {
                std::ofstream output(config_path);
                output
                    << "config_version = 1\n"
                    << "dimension = 2\n"
                    << "nx = 5\nny = 5\n"
                    << "length_x = 1\nlength_y = 1\n"
                    << "dt = 0.01\n"
                    << "[species.electrons]\n"
                    << "charge = -1\nmass = 1\nweight = 1\n"
                    << "particles = 8\n"
                    << "initialization_version = 1\n"
                    << "loading = quiet_start\n"
                    << "thermal_velocity = 2\n"
                    << "thermal_velocity_x = 0.1\n"
                    << "thermal_velocity_y = 0.2\n"
                    << "thermal_velocity_z = 0.3\n";
            }
            const auto config =
                pic::load_config_2d(config_path.string());
            std::filesystem::remove(config_path);
            require(
                config.species.front().initialization.loading ==
                    pic::ParticleLoading::QuietStart,
                "2D config did not parse quiet-start loading");
            require_near(
                *config.species.front().initialization.thermal_velocity_z,
                0.3, 1e-15,
                "2D config did not parse per-axis thermal velocity");
        }

        {
            pic::UnstructuredSimulation2DConfig config;
            config.mesh_path =
                std::filesystem::path(AURORA_TEST_SOURCE_DIR) /
                "tests" / "fixtures" / "tagged_square_v2.msh";
            config.dirichlet_potentials = {
                {"electrode", 0.0},
                {"inlet", 0.0},
                {"outlet", 0.0},
                {"wall", 0.0},
            };
            config.particle_boundaries = {
                {"electrode", pic::ParticleBoundary::Reflecting},
                {"inlet", pic::ParticleBoundary::Reflecting},
                {"outlet", pic::ParticleBoundary::Reflecting},
                {"wall", pic::ParticleBoundary::Reflecting},
            };
            pic::UnstructuredSpecies2DConfig species;
            species.name = "quiet";
            species.charge = 0.0;
            species.particles = 8;
            species.drift_velocity_x = 0.5;
            species.drift_velocity_y = -0.25;
            species.initialization.loading =
                pic::ParticleLoading::QuietStart;
            species.initialization.thermal_velocity_x = 0.2;
            species.initialization.thermal_velocity_y = 0.1;
            species.initialization.thermal_velocity_z = 0.0;
            config.species = {species};

            pic::UnstructuredSimulation2D simulation(config);
            simulation.initialize();
            std::vector<double> x;
            std::vector<double> y;
            for (const auto& particle :
                 simulation.species().front().particles()) {
                require(
                    simulation.mesh()
                        .locate_point(particle.position)
                        .has_value(),
                    "imported quiet-start particle lies outside geometry");
                x.push_back(particle.velocity.x);
                y.push_back(particle.velocity.y);
            }
            require_near(
                mean(x), species.drift_velocity_x, 1e-15,
                "imported quiet-start x drift is wrong");
            require_near(
                mean(y), species.drift_velocity_y, 1e-15,
                "imported quiet-start y drift is wrong");

            species.initialization_minimum = pic::Vec2{0.1, 0.1};
            species.initialization_maximum = pic::Vec2{0.9, 0.9};
            config.species = {species};
            require_throws(
                [&] {
                    pic::UnstructuredSimulation2D invalid(config);
                },
                "bounded imported quiet-start loading was silently accepted");
        }

        {
            pic::ParticleInitializationConfig invalid;
            invalid.version = 2;
            require_throws(
                [&] {
                    pic::validate_particle_initialization(
                        invalid, 3, 0.1, "test species");
                },
                "unsupported initialization version was accepted");
        }

        std::cout << "Initialization tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Initialization test failure: "
                  << error.what() << '\n';
        return 1;
    }
}
