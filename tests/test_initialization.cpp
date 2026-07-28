#include "pic/Config.hpp"
#include "pic/Grid.hpp"
#include "pic/Initialization.hpp"
#include "pic/Mesh2D.hpp"
#include "pic/Mesh3D.hpp"
#include "pic/Species.hpp"
#include "pic/Species2D.hpp"
#include "pic/Species3D.hpp"
#include "pic/Simulation.hpp"
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

std::string read_text(const std::filesystem::path& path) {
    std::ifstream input(path);
    return std::string(
        std::istreambuf_iterator<char>(input),
        std::istreambuf_iterator<char>());
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
            const auto moments =
                pic::summarize_initialization(species);
            require(
                moments.macroparticles == config.particles,
                "2D initialization report particle count is wrong");
            require_near(
                moments.mean_velocity_x, config.drift_velocity_x,
                1e-15,
                "2D initialization report mean velocity is wrong");
            require_near(
                moments.thermal_velocity_z, 0.0, 1e-15,
                "2D initialization report thermal velocity is wrong");
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
            const auto mesh_path =
                std::filesystem::path(AURORA_TEST_SOURCE_DIR) /
                "tests" / "fixtures" / "tagged_regions_v2.msh";
            const auto config_path =
                std::filesystem::path(
                    "test_imported_region_config.cfg");
            {
                std::ofstream output(config_path);
                output
                    << "config_version = 1\n"
                    << "dimension = 2\n"
                    << "mesh = imported\n"
                    << "mesh_file = " << mesh_path.string() << '\n'
                    << "dt = 0.01\nsteps = 0\n"
                    << "[boundary.inlet]\n"
                    << "field = dirichlet\npotential = 0\n"
                    << "particle = reflecting\n"
                    << "[boundary.outlet]\n"
                    << "field = neumann\nnormal_derivative = 0\n"
                    << "particle = reflecting\n"
                    << "[boundary.wall]\n"
                    << "field = neumann\nnormal_derivative = 0\n"
                    << "particle = reflecting\n"
                    << "[boundary.electrode]\n"
                    << "field = neumann\nnormal_derivative = 0\n"
                    << "particle = reflecting\n"
                    << "[species.quiet]\n"
                    << "charge = 0\nmass = 1\nweight = 1\n"
                    << "particles = 8\n"
                    << "initialization_region = region_a\n"
                    << "loading = quiet_start\n";
            }
            const auto parsed =
                pic::load_unstructured_config_2d(config_path);
            std::filesystem::remove(config_path);
            require(
                parsed.species.size() == 1 &&
                    parsed.species.front().initialization_region ==
                        "region_a",
                "imported config did not parse initialization_region");
        }

        {
            pic::UnstructuredSimulation2DConfig config;
            config.mesh_path =
                std::filesystem::path(AURORA_TEST_SOURCE_DIR) /
                "tests" / "fixtures" / "tagged_regions_v2.msh";
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
            species.initialization_region = "region_a";
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
                    simulation.mesh().locate_point(
                        particle.position).has_value(),
                    "imported quiet-start particle lies outside geometry");
                const auto location =
                    simulation.mesh().locate_point(
                        particle.position);
                require(
                    simulation.mesh().topology()
                            .cell_by_id(location->cell_id).label ==
                        species.initialization_region,
                    "imported particle was loaded outside its named physical region");
                x.push_back(particle.velocity.x);
                y.push_back(particle.velocity.y);
            }
            require_near(
                mean(x), species.drift_velocity_x, 1e-15,
                "imported quiet-start x drift is wrong");
            require_near(
                mean(y), species.drift_velocity_y, 1e-15,
                "imported quiet-start y drift is wrong");
            const auto moments = pic::summarize_initialization(
                simulation.species().front(),
                species.initialization_region);
            require(
                moments.region == "region_a" &&
                    moments.macroparticles == species.particles,
                "imported initialization summary lost its physical region");

            species.initialization_minimum = pic::Vec2{0.1, 0.1};
            species.initialization_maximum = pic::Vec2{0.9, 0.9};
            config.species = {species};
            require_throws(
                [&] {
                    pic::UnstructuredSimulation2D invalid(config);
                },
                "bounded imported quiet-start loading was silently accepted");

            species.initialization_minimum.reset();
            species.initialization_maximum.reset();
            species.initialization_region = "missing_region";
            config.species = {species};
            require_throws(
                [&] {
                    pic::UnstructuredSimulation2D invalid(config);
                },
                "unknown imported initialization region was accepted");
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

        {
            pic::SpeciesConfig config;
            config.particles = 2;
            config.drift_velocity = 1.0e12;
            config.thermal_velocity = 0.0;
            pic::Species species(config);
            species.particles() = {
                pic::Particle{0.25, 1.0e12 - 1.0, true},
                pic::Particle{0.75, 1.0e12 + 1.0, true},
            };
            const auto moments =
                pic::summarize_initialization(species);
            require_near(
                moments.mean_velocity_x, 1.0e12, 1e-3,
                "stable initialization mean is wrong");
            require_near(
                moments.thermal_velocity_x, 1.0, 1e-12,
                "initialization report lost a small thermal spread on a large drift");
        }

        {
            const auto output_directory =
                std::filesystem::path(
                    "test_output_initialization_report");
            std::filesystem::remove_all(output_directory);
            pic::Config config;
            config.steps = 0;
            config.output_interval = 1;
            config.output_dir = output_directory.string();
            config.checkpoint_output = true;
            config.checkpoint_interval = 1;
            config.species = {pic::SpeciesConfig{}};
            config.species.front().particles = 4;
            config.species.front().weight = 2.0;
            config.species.front().charge = -3.0;
            config.species.front().thermal_velocity = 0.0;
            config.species.front().drift_velocity = 0.5;
            pic::Simulation simulation(config);
            (void)simulation.run();
            const auto report =
                output_directory / "initialization.csv";
            require(
                std::filesystem::exists(report),
                "simulation run did not write initialization.csv");
            const std::string contents = read_text(report);
            require(
                contents.find(
                    "initialization_version,state_source,dimension") !=
                    std::string::npos &&
                    contents.find("\"generated\",1,\"electrons\"") !=
                    std::string::npos &&
                    contents.find(",-24,") != std::string::npos,
                "initialization.csv is missing schema or represented charge");

            pic::Config restart_config = config;
            const auto restart_output =
                std::filesystem::path(
                    "test_output_initialization_restart_report");
            std::filesystem::remove_all(restart_output);
            restart_config.output_dir = restart_output.string();
            restart_config.restart_path =
                (output_directory / "checkpoint_0.apc").string();
            pic::Simulation restarted(restart_config);
            (void)restarted.run();
            const std::string restart_contents = read_text(
                restart_output / "initialization.csv");
            require(
                restart_contents.find(
                    "\"restart\",1,\"electrons\",\"restart\",\"\"") !=
                    std::string::npos,
                "restart initialization report did not identify restored state");
            std::filesystem::remove_all(output_directory);
            std::filesystem::remove_all(restart_output);
        }

        std::cout << "Initialization tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Initialization test failure: "
                  << error.what() << '\n';
        return 1;
    }
}
