#include "pic/Config.hpp"
#include "pic/Grid.hpp"
#include "pic/Initialization.hpp"
#include "pic/Mesh2D.hpp"
#include "pic/Mesh3D.hpp"
#include "pic/ParticleState.hpp"
#include "pic/Species.hpp"
#include "pic/Species2D.hpp"
#include "pic/Species3D.hpp"
#include "pic/Simulation.hpp"
#include "pic/Simulation2D.hpp"
#include "pic/Simulation3D.hpp"
#include "pic/UnstructuredSimulation2D.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <numbers>
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

void write_particle_state(
    const std::filesystem::path& path,
    std::size_t dimension,
    std::size_t particle_count,
    const std::string& records,
    const std::string& units = "normalized") {
    std::ofstream output(path);
    output
        << "AuroraPIC-particle-state-v1\n"
        << "dimension " << dimension << '\n'
        << "units " << units << '\n'
        << "weighting species_constant\n"
        << "velocity_staggering time_centered\n"
        << "particle_count " << particle_count << '\n'
        << "records\n"
        << records
        << "end\n";
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
            pic::SpeciesConfig config;
            config.particles = 2000;
            config.thermal_velocity = 0.0;
            config.initialization.loading =
                pic::ParticleLoading::QuietStart;
            config.initialization.density_profile =
                pic::DensityProfileKind::Gaussian;
            config.initialization.profile_center_x = 0.72;
            config.initialization.profile_scale_x = 0.08;
            config.initialization.max_profile_sampling_attempts =
                100000;
            pic::Species species(config);
            pic::Grid grid(33, 1.0, pic::Boundary::Dirichlet);
            std::mt19937_64 rng(101);
            species.initialize(grid, rng);
            double mean_position = 0.0;
            for (const auto& particle : species.particles()) {
                mean_position += particle.x;
            }
            mean_position /= static_cast<double>(
                species.particles().size());
            require_near(
                mean_position, 0.72, 0.01,
                "1D Gaussian profile realized the wrong center");
            const auto moments =
                pic::summarize_initialization(species);
            require_near(
                moments.mean_position_x, mean_position, 1e-15,
                "Gaussian initialization audit position mean is wrong");
            require(
                moments.position_stddev_x > 0.06 &&
                    moments.position_stddev_x < 0.10,
                "Gaussian initialization audit position spread is wrong");
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
            pic::Species2DConfig config;
            config.particles = 2000;
            config.thermal_velocity = 0.0;
            config.initialization.loading =
                pic::ParticleLoading::QuietStart;
            config.initialization.density_profile =
                pic::DensityProfileKind::Sinusoidal;
            config.initialization.profile_amplitude = 1.0;
            config.initialization.profile_mode_x = 1;
            config.initialization.max_profile_sampling_attempts =
                10000;
            pic::Species2D species(config);
            pic::Mesh2D mesh(
                17, 17, 1.0, 1.0,
                pic::Boundary::Dirichlet);
            std::mt19937_64 rng(202);
            species.initialize(mesh, rng);
            double cosine_moment = 0.0;
            for (const auto& particle : species.particles()) {
                cosine_moment += std::cos(
                    2.0 * std::numbers::pi *
                    particle.position.x);
            }
            cosine_moment /= static_cast<double>(
                species.particles().size());
            require(
                cosine_moment > 0.4 &&
                    cosine_moment < 0.6,
                "2D sinusoidal profile realized the wrong spatial moment");
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
            pic::Species3DConfig config;
            config.particles = 1200;
            config.thermal_velocity = 0.0;
            config.initialization.density_profile =
                pic::DensityProfileKind::Gaussian;
            config.initialization.profile_center_x = 0.25;
            config.initialization.profile_center_y = 0.50;
            config.initialization.profile_center_z = 0.75;
            config.initialization.profile_scale_x = 0.08;
            config.initialization.profile_scale_y = 0.10;
            config.initialization.profile_scale_z = 0.08;
            config.initialization.max_profile_sampling_attempts =
                200000;
            pic::Species3D species(config);
            pic::Mesh3D mesh(
                9, 9, 9, 1.0, 1.0, 1.0,
                pic::Boundary::Dirichlet);
            std::mt19937_64 rng(303);
            species.initialize(mesh, rng);
            pic::Vec3 mean_position{};
            for (const auto& particle : species.particles()) {
                mean_position.x += particle.position.x;
                mean_position.y += particle.position.y;
                mean_position.z += particle.position.z;
            }
            const double inverse = 1.0 /
                static_cast<double>(species.particles().size());
            mean_position.x *= inverse;
            mean_position.y *= inverse;
            mean_position.z *= inverse;
            require_near(
                mean_position.x, 0.25, 0.02,
                "3D Gaussian profile x center is wrong");
            require_near(
                mean_position.y, 0.50, 0.02,
                "3D Gaussian profile y center is wrong");
            require_near(
                mean_position.z, 0.75, 0.02,
                "3D Gaussian profile z center is wrong");
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
                    << "initial_state_path = initial_particles.aps\n"
                    << "initialization_max_relative_charge_imbalance = 0.01\n"
                    << "initialization_max_relative_current_imbalance = 0.02\n"
                    << "initialization_max_relative_pair_imbalance = 0.03\n"
                    << "initialization_charge_pairs = electrons:ions\n"
                    << "[species.electrons]\n"
                    << "charge = -1\nmass = 1\nweight = 1\n"
                    << "particles = 8\n"
                    << "initialization_version = 1\n"
                    << "loading = quiet_start\n"
                    << "thermal_velocity = 2\n"
                    << "thermal_velocity_x = 0.1\n"
                    << "thermal_velocity_y = 0.2\n"
                    << "thermal_velocity_z = 0.3\n";
                output
                    << "density_profile = sinusoidal\n"
                    << "profile_amplitude = 0.5\n"
                    << "profile_mode_x = 2\n"
                    << "max_profile_sampling_attempts = 1000\n"
                    << "[species.ions]\n"
                    << "charge = 1\nmass = 10\nweight = 1\n"
                    << "particles = 8\nthermal_velocity = 0\n";
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
            require(
                config.species.front().initialization.density_profile ==
                        pic::DensityProfileKind::Sinusoidal &&
                    config.species.front().initialization
                            .profile_mode_x.value_or(0) == 2,
                "2D config did not parse analytic density profile");
            require(
                config.initialization_acceptance
                        .max_relative_charge_imbalance.value_or(-1.0) ==
                    0.01 &&
                    config.initialization_acceptance.charge_pairs.size() ==
                        1 &&
                    config.initialization_acceptance.charge_pairs.front()
                            .second_species == "ions" &&
                    config.initial_state_path.filename() ==
                        "initial_particles.aps",
                "2D config did not parse initialization acceptance gates");
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
                    << "initial_state_path = imported_particles.aps\n"
                    << "initialization_max_relative_pair_imbalance = 0.05\n"
                    << "initialization_charge_pairs = quiet:ions\n"
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
                    << "loading = quiet_start\n"
                    << "density_profile = sinusoidal\n"
                    << "profile_amplitude = 0.5\n"
                    << "profile_mode_x = 1\n"
                    << "[species.ions]\n"
                    << "charge = 1\nmass = 10\nweight = 1\n"
                    << "particles = 8\n"
                    << "initialization_region = region_a\n"
                    << "loading = quiet_start\n";
            }
            const auto parsed =
                pic::load_unstructured_config_2d(config_path);
            std::filesystem::remove(config_path);
            require(
                parsed.species.size() == 2 &&
                    parsed.species.front().initialization_region ==
                        "region_a" &&
                    parsed.species.front().initialization
                            .density_profile ==
                        pic::DensityProfileKind::Sinusoidal &&
                    parsed.initialization_acceptance.charge_pairs.size() ==
                        1 &&
                    parsed.initial_state_path.filename() ==
                        "imported_particles.aps",
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
            species.initialization.density_profile =
                pic::DensityProfileKind::Sinusoidal;
            species.initialization.profile_amplitude = 0.5;
            species.initialization.profile_mode_x = 1;
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
                    moments.macroparticles == species.particles &&
                    moments.density_profile == "sinusoidal",
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

            pic::ParticleInitializationConfig invalid_uniform;
            invalid_uniform.profile_amplitude = 0.5;
            require_throws(
                [&] {
                    pic::validate_density_profile(
                        invalid_uniform, 1, 10,
                        "test species");
                },
                "uniform profile silently accepted profile parameters");

            pic::ParticleInitializationConfig invalid_gaussian;
            invalid_gaussian.density_profile =
                pic::DensityProfileKind::Gaussian;
            invalid_gaussian.profile_center_x = 0.5;
            require_throws(
                [&] {
                    pic::validate_density_profile(
                        invalid_gaussian, 1, 10,
                        "test species");
                },
                "Gaussian profile without scale was accepted");

            pic::ParticleInitializationConfig invalid_sinusoidal;
            invalid_sinusoidal.density_profile =
                pic::DensityProfileKind::Sinusoidal;
            invalid_sinusoidal.profile_amplitude = 1.1;
            invalid_sinusoidal.profile_mode_x = 1;
            require_throws(
                [&] {
                    pic::validate_density_profile(
                        invalid_sinusoidal, 1, 10,
                        "test species");
                },
                "sinusoidal profile accepted amplitude above one");
        }

        {
            pic::SpeciesConfig config;
            config.particles = 1;
            config.thermal_velocity = 0.0;
            config.initialization.density_profile =
                pic::DensityProfileKind::Gaussian;
            config.initialization.profile_center_x = 100.0;
            config.initialization.profile_scale_x = 0.1;
            config.initialization.max_profile_sampling_attempts = 10;
            pic::Species species(config);
            pic::Grid grid(5, 1.0, pic::Boundary::Dirichlet);
            std::mt19937_64 rng(404);
            require_throws(
                [&] { species.initialize(grid, rng); },
                "exhausted density-profile work budget did not fail");
        }

        {
            const auto state_path =
                std::filesystem::path(
                    "test_external_particle_state.aps");
            write_particle_state(
                state_path, 1, 2,
                "particle electrons 0.25 0 0 1.5 0 0\n"
                "particle electrons 0.75 0 0 -0.5 0 0\n");
            const auto state =
                pic::load_external_particle_state(
                    state_path, 2);
            require(
                state.spatial_dimension == 1 &&
                    state.particle_count == 2 &&
                    state.species.at("electrons").size() == 2,
                "external particle-state loader lost metadata or records");
            pic::validate_external_particle_state(
                state, 1, pic::UnitSystem::Normalized,
                {{"electrons", 2}}, "test");
            require_throws(
                [&] {
                    (void)pic::load_external_particle_state(
                        state_path, 1);
                },
                "external particle-state count limit was ignored");
            require_throws(
                [&] {
                    pic::validate_external_particle_state(
                        state, 2,
                        pic::UnitSystem::Normalized,
                        {{"electrons", 2}}, "test");
                },
                "external particle-state dimension mismatch was accepted");
            require_throws(
                [&] {
                    pic::validate_external_particle_state(
                        state, 1, pic::UnitSystem::SI,
                        {{"electrons", 2}}, "test");
                },
                "external particle-state unit mismatch was accepted");
            require_throws(
                [&] {
                    pic::validate_external_particle_state(
                        state, 1,
                        pic::UnitSystem::Normalized,
                        {{"ions", 2}}, "test");
                },
                "external particle-state species mismatch was accepted");

            const auto output_directory =
                std::filesystem::path(
                    "test_output_external_state_1d");
            std::filesystem::remove_all(output_directory);
            pic::Config config;
            config.steps = 0;
            config.output_dir = output_directory.string();
            config.initial_state_path = state_path;
            config.species = {pic::SpeciesConfig{}};
            config.species.front().particles = 2;
            config.species.front().thermal_velocity = 0.0;
            pic::Simulation simulation(config);
            (void)simulation.run();
            const auto& particles =
                simulation.species().front().particles();
            require_near(
                particles[0].x, 0.25, 1e-15,
                "1D external position changed");
            require_near(
                particles[0].v, 1.5, 1e-15,
                "1D external velocity changed");
            require(
                read_text(
                    output_directory /
                    "initialization.csv")
                        .find(
                            "\"external\",1,\"electrons\",\"external\",\"external\"") !=
                    std::string::npos,
                "1D external state was not identified in its audit");
            std::filesystem::remove_all(output_directory);
            pic::Config bounded_config = config;
            bounded_config.output_dir =
                "test_output_external_state_bounds";
            bounded_config.species.front().init_x_min =
                0.3;
            std::filesystem::remove_all(
                bounded_config.output_dir);
            pic::Simulation bounded_simulation(
                bounded_config);
            require_throws(
                [&] {
                    (void)bounded_simulation.run();
                },
                "external state escaped its structured initialization interval");
            std::filesystem::remove_all(
                bounded_config.output_dir);

            write_particle_state(
                state_path, 2, 2,
                "particle electrons 0.2 0.3 0 1 2 3\n"
                "particle electrons 0.8 0.7 0 -1 -2 -3\n");
            pic::Simulation2DConfig config_2d;
            pic::Species2DConfig species_2d;
            species_2d.particles = 2;
            species_2d.thermal_velocity = 0.0;
            config_2d.species = {species_2d};
            config_2d.initial_state_path = state_path;
            pic::Simulation2D simulation_2d(config_2d);
            simulation_2d.initialize();
            require_near(
                simulation_2d.species().front()
                    .particles()[0].position.y,
                0.3, 1e-15,
                "2D external position changed");
            require_near(
                simulation_2d.species().front()
                    .particles()[0].velocity_z,
                3.0, 1e-15,
                "2D3V external velocity changed");

            write_particle_state(
                state_path, 3, 1,
                "particle electrons 0.2 0.3 0.4 1 2 3\n");
            pic::Simulation3DConfig config_3d;
            pic::Species3DConfig species_3d;
            species_3d.particles = 1;
            species_3d.thermal_velocity = 0.0;
            config_3d.species = {species_3d};
            config_3d.initial_state_path = state_path;
            pic::Simulation3D simulation_3d(config_3d);
            simulation_3d.initialize();
            require_near(
                simulation_3d.species().front()
                    .particles()[0].position.z,
                0.4, 1e-15,
                "3D external position changed");
            require_near(
                simulation_3d.species().front()
                    .particles()[0].velocity.z,
                3.0, 1e-15,
                "3D external velocity changed");

            write_particle_state(
                state_path, 2, 1,
                "particle imported 0.1 0.1 0 0.4 -0.2 0.3\n");
            pic::UnstructuredSimulation2DConfig imported_config;
            imported_config.mesh_path =
                std::filesystem::path(
                    AURORA_TEST_SOURCE_DIR) /
                "tests" / "fixtures" /
                "tagged_regions_v2.msh";
            imported_config.dirichlet_potentials = {
                {"electrode", 0.0},
                {"inlet", 0.0},
                {"outlet", 0.0},
                {"wall", 0.0},
            };
            imported_config.particle_boundaries = {
                {"electrode",
                 pic::ParticleBoundary::Reflecting},
                {"inlet",
                 pic::ParticleBoundary::Reflecting},
                {"outlet",
                 pic::ParticleBoundary::Reflecting},
                {"wall",
                 pic::ParticleBoundary::Reflecting},
            };
            pic::UnstructuredSpecies2DConfig imported_species;
            imported_species.name = "imported";
            imported_species.charge = 0.0;
            imported_species.particles = 1;
            imported_species.initialization_region =
                "region_a";
            imported_config.species = {imported_species};
            imported_config.initial_state_path =
                state_path;
            pic::UnstructuredSimulation2D imported(
                imported_config);
            imported.initialize();
            require_near(
                imported.species().front()
                    .particles()[0].velocity_z,
                0.3, 1e-15,
                "imported-geometry external velocity changed");

            write_particle_state(
                state_path, 1, 1,
                "particle electrons 0.5 0.1 0 0 0 0\n");
            require_throws(
                [&] {
                    (void)pic::load_external_particle_state(
                        state_path, 1);
                },
                "1D external state accepted an inactive component");
            std::filesystem::remove(state_path);

            pic::Config conflicting;
            conflicting.restart_path = "restart.apc";
            conflicting.initial_state_path = "initial.aps";
            require_throws(
                [&] {
                    pic::Simulation invalid(conflicting);
                },
                "restart and external initial state were accepted together");
        }

        {
            pic::InitializationSpeciesMoments electrons;
            electrons.species = "electrons";
            electrons.represented_charge = -10.0;
            electrons.mean_velocity_x = 2.0;
            pic::InitializationSpeciesMoments ions;
            ions.species = "ions";
            ions.represented_charge = 10.0;
            ions.mean_velocity_x = 2.0;

            pic::InitializationAcceptanceConfig acceptance;
            acceptance.max_relative_charge_imbalance = 0.0;
            acceptance.max_relative_current_imbalance = 0.0;
            acceptance.max_relative_pair_imbalance = 0.0;
            acceptance.charge_pairs.push_back(
                {"electrons", "ions"});
            const auto passed =
                pic::assess_initialization_acceptance(
                    acceptance, {electrons, ions}, 1);
            require(
                passed.enabled && passed.passed &&
                    passed.metrics.size() == 3,
                "balanced charge/current initialization did not pass");

            ions.mean_velocity_x = 0.0;
            const auto current_failure =
                pic::assess_initialization_acceptance(
                    acceptance, {electrons, ions}, 1);
            require(
                !current_failure.passed &&
                    current_failure.metrics[1]
                            .relative_residual == 1.0,
                "net-current imbalance was not detected");
            require_throws(
                [&] {
                    pic::enforce_initialization_acceptance(
                        current_failure);
                },
                "failed initialization acceptance was not enforced");

            ions.mean_velocity_x = 2.0;
            ions.represented_charge = 8.0;
            const auto charge_failure =
                pic::assess_initialization_acceptance(
                    acceptance, {electrons, ions}, 1);
            require(
                !charge_failure.passed &&
                    !charge_failure.metrics.front().passed &&
                    !charge_failure.metrics.back().passed,
                "charge and pair imbalance were not detected");
            ions.represented_charge = -10.0;
            const auto same_sign_failure =
                pic::assess_initialization_acceptance(
                    acceptance, {electrons, ions}, 1);
            require(
                !same_sign_failure.metrics.back().passed,
                "same-sign charge pair was accepted");
            acceptance.charge_pairs.front().second_species =
                "missing";
            require_throws(
                [&] {
                    (void)pic::assess_initialization_acceptance(
                        acceptance, {electrons, ions}, 1);
                },
                "unknown charge-pair species was accepted");
            acceptance.charge_pairs.front().second_species =
                "ions";

            pic::InitializationAcceptanceConfig invalid;
            invalid.max_relative_pair_imbalance = 0.1;
            require_throws(
                [&] {
                    pic::validate_initialization_acceptance(
                        invalid, "test");
                },
                "pair tolerance without named pairs was accepted");
            invalid = {};
            invalid.max_relative_charge_imbalance = 1.1;
            require_throws(
                [&] {
                    pic::validate_initialization_acceptance(
                        invalid, "test");
                },
                "out-of-range initialization tolerance was accepted");

            const auto acceptance_path =
                std::filesystem::path(
                    "test_initialization_acceptance.csv");
            pic::write_initialization_acceptance_report(
                acceptance_path, charge_failure);
            const std::string acceptance_contents =
                read_text(acceptance_path);
            std::filesystem::remove(acceptance_path);
            require(
                acceptance_contents.find(
                    "\"net_charge\"") != std::string::npos &&
                    acceptance_contents.find(
                        "\"charge_pair:electrons:ions\"") !=
                        std::string::npos,
                "initialization acceptance report lost gate metrics");
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
            require(
                read_text(
                    output_directory /
                    "initialization_acceptance.csv")
                        .find("\"acceptance gates disabled\"") !=
                    std::string::npos,
                "simulation run did not audit disabled acceptance gates");
            const std::string contents = read_text(report);
            require(
                contents.find(
                    "initialization_version,state_source,dimension") !=
                    std::string::npos &&
                    contents.find(
                        "\"generated\",1,\"electrons\",\"random\",\"uniform\"") !=
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
                    "\"restart\",1,\"electrons\",\"restart\",\"restart\",\"\"") !=
                    std::string::npos,
                "restart initialization report did not identify restored state");
            std::filesystem::remove_all(output_directory);
            std::filesystem::remove_all(restart_output);
        }

        {
            const auto output_directory =
                std::filesystem::path(
                    "test_output_initialization_gate_failure");
            std::filesystem::remove_all(output_directory);
            pic::Config config;
            config.steps = 0;
            config.output_dir = output_directory.string();
            config.species = {pic::SpeciesConfig{}};
            config.species.front().particles = 4;
            config.species.front().thermal_velocity = 0.0;
            config.initialization_acceptance
                .max_relative_charge_imbalance = 0.01;
            pic::Simulation simulation(config);
            require_throws(
                [&] { (void)simulation.run(); },
                "non-neutral simulation passed its initialization gate");
            const auto report =
                output_directory /
                "initialization_acceptance.csv";
            require(
                std::filesystem::exists(report) &&
                    read_text(report).find(
                        "1,0,\"net_charge\"") !=
                        std::string::npos,
                "failed initialization gate did not leave an audit report");
            std::filesystem::remove_all(output_directory);
        }

        std::cout << "Initialization tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Initialization test failure: "
                  << error.what() << '\n';
        return 1;
    }
}
