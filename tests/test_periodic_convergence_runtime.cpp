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

void write_v24_without_periodic_state(
    const std::filesystem::path& source,
    const std::filesystem::path& destination) {
    std::ifstream input(source);
    std::ofstream output(destination);
    require(static_cast<bool>(input) && static_cast<bool>(output),
            "cannot open synthetic v24 checkpoint streams");
    std::string line;
    bool first = true;
    while (std::getline(input, line)) {
        if (first) {
            output << "AuroraPIC-checkpoint-v24\n";
            first = false;
        } else if (!line.starts_with("periodic_convergence")) {
            output << line << '\n';
        }
    }
    output.flush();
    require(input.eof() && static_cast<bool>(output),
            "failed to create synthetic v24 checkpoint");
}

void test_phase_aligned_legacy_restart_can_begin_fresh_epoch() {
    const auto root = std::filesystem::path(
        "test_output_periodic_convergence_legacy_reset");
    std::filesystem::remove_all(root);
    std::filesystem::create_directories(root);
    const auto v25 = root / "disabled-v25.apc";
    const auto v24 = root / "legacy-v24.apc";
    const auto misaligned_v25 = root / "misaligned-v25.apc";
    const auto misaligned_v24 = root / "misaligned-v24.apc";

    auto source_cfg = stationary_config(root / "source");
    source_cfg.mode = pic::RunMode::Transient;
    source_cfg.periodic_convergence.enabled = false;
    source_cfg.periodic_convergence.rf_frequency = 0.0;
    source_cfg.periodic_convergence.cycles_per_block = 0;
    pic::Simulation source(source_cfg);
    source.initialize();
    for (std::size_t step = 0; step < 8; ++step) source.step();
    source.save_checkpoint(v25);
    write_v24_without_periodic_state(v25, v24);

    auto strict_cfg = stationary_config(root / "strict");
    bool strict_rejected = false;
    try {
        pic::Simulation strict(strict_cfg);
        strict.load_checkpoint(v24);
    } catch (const std::exception&) {
        strict_rejected = true;
    }
    require(strict_rejected,
            "legacy checkpoint must not invent convergence history by default");

    auto reset_cfg = stationary_config(root / "reset");
    reset_cfg.periodic_convergence.reset_on_restart = true;
    reset_cfg.restart_path = v24.string();
    pic::Simulation reset(reset_cfg);
    const auto summary = reset.run();
    require(summary.steady_state_reached && summary.steps_completed == 40,
            "fresh epoch should begin at legacy cycle two and require eight "
            "new complete cycles");
    const auto results = reset.periodic_convergence_results();
    require(results.front().blocks == 4 && results.front().converged(),
            "fresh legacy convergence epoch should contain only new blocks");

    pic::Simulation misaligned(source_cfg);
    misaligned.initialize();
    for (std::size_t step = 0; step < 6; ++step) misaligned.step();
    misaligned.save_checkpoint(misaligned_v25);
    write_v24_without_periodic_state(misaligned_v25, misaligned_v24);
    auto misaligned_cfg = stationary_config(root / "misaligned");
    misaligned_cfg.periodic_convergence.reset_on_restart = true;
    bool misaligned_rejected = false;
    try {
        pic::Simulation rejected(misaligned_cfg);
        rejected.load_checkpoint(misaligned_v24);
    } catch (const std::exception&) {
        misaligned_rejected = true;
    }
    require(misaligned_rejected,
            "restart reset must reject a checkpoint between RF phases");
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
        "periodic_convergence_reset_on_restart = true\n"
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
            cfg.periodic_convergence.reset_on_restart &&
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
        test_phase_aligned_legacy_restart_can_begin_fresh_epoch();
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
