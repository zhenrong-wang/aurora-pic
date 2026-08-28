#include "pic/Simulation.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

pic::Config stationary_config(const std::filesystem::path& output) {
    pic::Config cfg;
    cfg.nx = 16;
    cfg.length = 1.0;
    cfg.dt = 0.25;
    cfg.output_interval = 4;
    cfg.output_dir = output.string();
    cfg.boundary = pic::Boundary::Periodic;
    cfg.mode = pic::RunMode::SteadyState;
    cfg.max_steps = 64;
    cfg.runtime.backend = pic::RuntimeBackend::Serial;
    cfg.periodic_convergence.enabled = true;
    cfg.periodic_convergence.rf_frequency = 1.0;
    cfg.periodic_convergence.cycles_per_block = 2;
    cfg.periodic_convergence.minimum_blocks = 4;
    cfg.periodic_convergence.minimum_effective_blocks = 4.0;
    cfg.periodic_convergence.maximum_absolute_projected_fractional_drift =
        1e-12;
    cfg.periodic_convergence.maximum_absolute_split_half_fractional_change =
        1e-12;
    cfg.periodic_convergence.maximum_relative_standard_error = 1e-12;
    pic::SpeciesConfig species;
    species.name = "neutral";
    species.charge = 0.0;
    species.mass = 1.0;
    species.weight = 2.0;
    species.particles = 8;
    species.drift_velocity = 0.125;
    species.thermal_velocity = 0.0;
    cfg.species = {species};
    return cfg;
}

void test_periodic_steady_state_terminates_on_complete_blocks() {
    const auto output = std::filesystem::path(
        "test_output_periodic_convergence");
    std::filesystem::remove_all(output);
    pic::Simulation simulation(stationary_config(output));
    const auto summary = simulation.run();
    require(summary.steady_state_reached,
            "periodic stationary case should converge");
    require(summary.steps_completed == 32,
            "termination must occur after four complete two-cycle blocks");
    const auto results = simulation.periodic_convergence_results();
    require(results.size() == 2,
            "population and phase-energy observables are required");
    require(results[0].converged() && results[1].converged(),
            "all periodic observables should pass");
    require(std::filesystem::exists(
                output / "periodic_convergence_blocks.csv") &&
            std::filesystem::exists(
                output / "periodic_convergence_status.csv"),
            "periodic convergence diagnostics should be written");
}

void test_checkpoint_split_preserves_periodic_history() {
    const auto root = std::filesystem::path(
        "test_output_periodic_convergence_checkpoint");
    std::filesystem::remove_all(root);
    std::filesystem::create_directories(root);
    const auto checkpoint = root / "split.apc";

    auto prefix_cfg = stationary_config(root / "prefix");
    pic::Simulation prefix(prefix_cfg);
    prefix.initialize();
    for (std::size_t step = 0; step < 20; ++step) prefix.step();
    prefix.save_checkpoint(checkpoint);
    const auto prefix_results = prefix.periodic_convergence_results();
    require(prefix_results.front().blocks == 2,
            "prefix should contain two complete blocks and one partial block");

    auto restart_cfg = stationary_config(root / "restart");
    restart_cfg.restart_path = checkpoint.string();
    pic::Simulation restarted(restart_cfg);
    const auto restarted_summary = restarted.run();

    auto continuous_cfg = stationary_config(root / "continuous");
    pic::Simulation continuous(continuous_cfg);
    const auto continuous_summary = continuous.run();
    require(restarted_summary.steps_completed ==
                continuous_summary.steps_completed &&
            restarted_summary.steady_state_reached &&
            continuous_summary.steady_state_reached,
            "checkpoint split and continuous termination must agree");
    const auto restarted_results = restarted.periodic_convergence_results();
    const auto continuous_results = continuous.periodic_convergence_results();
    require(restarted_results.size() == continuous_results.size(),
            "checkpoint result shape must agree");
    for (std::size_t i = 0; i < restarted_results.size(); ++i) {
        require(restarted_results[i].blocks == continuous_results[i].blocks &&
                restarted_results[i].mean == continuous_results[i].mean &&
                restarted_results[i].effective_blocks ==
                    continuous_results[i].effective_blocks &&
                restarted_results[i].classification ==
                    continuous_results[i].classification,
                "checkpoint convergence result must be exact");
    }
}

void test_driven_steady_state_requires_matching_periodic_contract() {
    auto cfg = stationary_config("test_output_periodic_driven_validation");
    cfg.boundary = pic::Boundary::Dirichlet;
    cfg.phi_left_drive.amplitude = 1.0;
    cfg.phi_left_drive.frequency = 1.0;
    pic::Simulation accepted(cfg);
    (void)accepted;

    cfg.periodic_convergence.rf_frequency = 0.8;
    bool mismatch_rejected = false;
    try {
        pic::Simulation rejected(cfg);
    } catch (const std::exception&) {
        mismatch_rejected = true;
    }
    require(mismatch_rejected,
            "drive and periodic convergence frequencies must match");

    cfg.periodic_convergence.enabled = false;
    cfg.periodic_convergence.rf_frequency = 0.0;
    cfg.periodic_convergence.cycles_per_block = 0;
    bool missing_rejected = false;
    try {
        pic::Simulation rejected(cfg);
    } catch (const std::exception&) {
        missing_rejected = true;
    }
    require(missing_rejected,
            "driven steady state must require periodic convergence");
}

void test_periodic_configuration_keys_load() {
    const auto root = std::filesystem::path(
        "test_output_periodic_convergence_config");
    std::filesystem::remove_all(root);
    std::filesystem::create_directories(root);
    const auto path = root / "case.cfg";
    std::ofstream stream(path);
    stream <<
        "config_version = 1\n"
        "dimension = 1\n"
        "nx = 16\n"
        "length = 1\n"
        "dt = 0.25\n"
        "output_interval = 4\n"
        "mode = steady_state\n"
        "max_steps = 64\n"
        "boundary = dirichlet\n"
        "phi_left_amplitude = 1\n"
        "phi_left_frequency = 1\n"
        "periodic_convergence = true\n"
        "periodic_convergence_rf_frequency = 1\n"
        "periodic_convergence_cycles_per_block = 3\n"
        "periodic_convergence_minimum_blocks = 5\n"
        "periodic_convergence_minimum_effective_blocks = 4.5\n"
        "periodic_convergence_maximum_absolute_projected_fractional_drift = 0.02\n"
        "periodic_convergence_maximum_absolute_split_half_fractional_change = 0.03\n"
        "periodic_convergence_maximum_relative_standard_error = 0.04\n"
        "output_dir = test_output_periodic_convergence_config/run\n"
        "[species.neutral]\n"
        "charge = 0\n"
        "mass = 1\n"
        "weight = 1\n"
        "particles = 8\n"
        "thermal_velocity = 0\n";
    stream.close();
    const auto cfg = pic::load_config(path.string());
    require(cfg.periodic_convergence.enabled &&
            cfg.periodic_convergence.rf_frequency == 1.0 &&
            cfg.periodic_convergence.cycles_per_block == 3 &&
            cfg.periodic_convergence.minimum_blocks == 5 &&
            cfg.periodic_convergence.minimum_effective_blocks == 4.5 &&
            cfg.periodic_convergence
                    .maximum_absolute_projected_fractional_drift == 0.02 &&
            cfg.periodic_convergence
                    .maximum_absolute_split_half_fractional_change == 0.03 &&
            cfg.periodic_convergence.maximum_relative_standard_error == 0.04,
            "periodic convergence configuration keys should round-trip");
}

} // namespace

int main() {
    try {
        test_periodic_steady_state_terminates_on_complete_blocks();
        test_checkpoint_split_preserves_periodic_history();
        test_driven_steady_state_requires_matching_periodic_contract();
        test_periodic_configuration_keys_load();
        std::cout << "periodic convergence runtime tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "periodic convergence runtime test failed: "
                  << error.what() << '\n';
        return 1;
    }
}
