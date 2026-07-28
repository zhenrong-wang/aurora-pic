#pragma once

#include "pic/GasDataset.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace pic {

struct SwarmBenchmarkConfig {
    std::filesystem::path gas_data_file{};
    double neutral_density{0.0};
    std::vector<double> reduced_fields_td{};
    double max_frequency{0.0};
    double timestep{0.0};
    std::size_t steps{0};
    std::size_t burn_in_steps{0};
    std::size_t particles{0};
    std::size_t uncertainty_blocks{10};
    std::uint64_t work_item_limit{100000000};
    double initial_mean_energy_ev{0.0};
    double max_energy_ev{0.0};
    std::uint64_t seed{12345};
    std::filesystem::path output_file{"swarm.csv"};
};

struct SwarmChannelResult {
    std::string name{};
    std::uint64_t collisions{0};
    double rate_per_electron_s{0.0};
    double poisson_standard_error_s{0.0};
};

struct SwarmBenchmarkResult {
    double reduced_field_td{0.0};
    double electric_field_v_m{0.0};
    double mean_velocity_x_m_s{0.0};
    double mean_velocity_x_standard_error_m_s{0.0};
    double electron_drift_velocity_m_s{0.0};
    double reduced_mobility_1_v_m_s{0.0};
    double mean_energy_ev{0.0};
    double mean_energy_standard_error_ev{0.0};
    double longitudinal_diffusion_m2_s{0.0};
    double transverse_diffusion_m2_s{0.0};
    double maximum_observed_energy_ev{0.0};
    std::uint64_t collision_candidates{0};
    std::uint64_t null_collisions{0};
    std::vector<SwarmChannelResult> channels{};
};

SwarmBenchmarkConfig load_swarm_benchmark_config(
    const std::filesystem::path& path);

std::vector<SwarmBenchmarkResult> run_swarm_benchmark(
    const SwarmBenchmarkConfig& config);

void write_swarm_benchmark_csv(
    const std::filesystem::path& path,
    const SwarmBenchmarkConfig& config,
    const GasDataset& dataset,
    const std::vector<SwarmBenchmarkResult>& results);

} // namespace pic
