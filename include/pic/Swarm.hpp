#pragma once

#include "pic/GasDataset.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace pic {

enum class SwarmPopulationModel {
    FixedPopulationNoAvalanche,
    BranchingResampled
};

inline std::string to_string(SwarmPopulationModel model) {
    switch (model) {
        case SwarmPopulationModel::FixedPopulationNoAvalanche:
            return "fixed_population_no_avalanche";
        case SwarmPopulationModel::BranchingResampled:
            return "branching_resampled";
    }
    return "unknown";
}

struct SwarmBenchmarkConfig {
    std::filesystem::path gas_data_file{};
    double neutral_density{0.0};
    double neutral_temperature{0.0};
    std::vector<double> reduced_fields_td{};
    double max_frequency{0.0};
    double timestep{0.0};
    std::size_t steps{0};
    std::size_t burn_in_steps{0};
    std::size_t particles{0};
    SwarmPopulationModel population_model{
        SwarmPopulationModel::FixedPopulationNoAvalanche};
    std::size_t population_limit{0};
    std::size_t uncertainty_blocks{10};
    std::uint64_t work_item_limit{100000000};
    double initial_mean_energy_ev{0.0};
    double max_energy_ev{0.0};
    std::uint64_t seed{12345};
    std::filesystem::path output_file{"swarm.csv"};
    std::size_t spatial_histories{0};
    double spatial_length_m{0.0};
    std::size_t spatial_bins{0};
    std::size_t spatial_fit_begin_bin{0};
    std::size_t spatial_fit_end_bin{0};
    std::size_t spatial_max_steps{0};
    std::uint64_t spatial_work_item_limit{100000000};
    double spatial_min_r_squared{0.0};
    std::filesystem::path spatial_profile_file{
        "swarm-spatial-profile.csv"};
};

struct SwarmChannelResult {
    std::string name{};
    std::uint64_t collisions{0};
    double rate_per_electron_s{0.0};
    double poisson_standard_error_s{0.0};
};

struct SwarmSpatialFluxPoint {
    double distance_m{0.0};
    double net_crossings_per_injected_electron{0.0};
    double standard_error{0.0};
};

struct SwarmBenchmarkResult {
    std::uint64_t collision_model_signature{0};
    double neutral_velocity_stddev_m_s{0.0};
    double neutral_speed_limit_sigma{0.0};
    double reduced_field_td{0.0};
    double electric_field_v_m{0.0};
    double mean_velocity_x_m_s{0.0};
    double mean_velocity_x_standard_error_m_s{0.0};
    double electron_drift_velocity_m_s{0.0};
    double electron_drift_velocity_first_half_m_s{0.0};
    double electron_drift_velocity_second_half_m_s{0.0};
    double reduced_mobility_1_v_m_s{0.0};
    double mean_energy_ev{0.0};
    double mean_energy_standard_error_ev{0.0};
    double mean_energy_first_half_ev{0.0};
    double mean_energy_second_half_ev{0.0};
    double longitudinal_diffusion_m2_s{0.0};
    double transverse_diffusion_m2_s{0.0};
    bool diffusion_available{true};
    double maximum_observed_energy_ev{0.0};
    double initial_total_electron_weight{0.0};
    double final_total_electron_weight{0.0};
    std::size_t final_computational_particles{0};
    double temporal_growth_rate_s{0.0};
    double temporal_growth_rate_standard_error_s{0.0};
    double ionization_rate_s{0.0};
    double ionization_rate_standard_error_s{0.0};
    double attachment_rate_s{0.0};
    double attachment_rate_standard_error_s{0.0};
    double net_creation_rate_s{0.0};
    double net_creation_rate_standard_error_s{0.0};
    bool townsend_available{false};
    double growth_over_flux_drift_townsend_1_m{0.0};
    double growth_over_flux_drift_townsend_standard_error_1_m{0.0};
    double rate_balance_effective_townsend_1_m{0.0};
    double rate_balance_effective_townsend_standard_error_1_m{0.0};
    bool spatial_townsend_available{false};
    double spatial_flux_townsend_1_m{0.0};
    double spatial_flux_townsend_standard_error_1_m{0.0};
    double spatial_flux_fit_r_squared{0.0};
    std::size_t spatial_histories_completed{0};
    std::size_t spatial_maximum_active_particles{0};
    std::uint64_t spatial_particle_updates{0};
    std::vector<SwarmSpatialFluxPoint> spatial_flux_profile{};
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

void write_swarm_spatial_profile_csv(
    const std::filesystem::path& path,
    const SwarmBenchmarkConfig& config,
    const GasDataset& dataset,
    const std::vector<SwarmBenchmarkResult>& results);

} // namespace pic
