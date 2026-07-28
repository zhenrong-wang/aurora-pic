#pragma once

#include <array>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <random>
#include <string>
#include <vector>

namespace pic {

class Species;
class Species2D;
class Species3D;

enum class ParticleLoading {
    Random,
    QuietStart,
};

enum class DensityProfileKind {
    Uniform,
    Gaussian,
    Sinusoidal,
};

struct ParticleInitializationConfig {
    std::size_t version{1};
    ParticleLoading loading{ParticleLoading::Random};
    std::optional<double> thermal_velocity_x{};
    std::optional<double> thermal_velocity_y{};
    std::optional<double> thermal_velocity_z{};
    DensityProfileKind density_profile{DensityProfileKind::Uniform};
    std::optional<double> profile_center_x{};
    std::optional<double> profile_center_y{};
    std::optional<double> profile_center_z{};
    std::optional<double> profile_scale_x{};
    std::optional<double> profile_scale_y{};
    std::optional<double> profile_scale_z{};
    std::optional<double> profile_amplitude{};
    std::optional<double> profile_phase{};
    std::optional<std::size_t> profile_mode_x{};
    std::optional<std::size_t> profile_mode_y{};
    std::optional<std::size_t> profile_mode_z{};
    std::size_t max_profile_sampling_attempts{1000000};
};

struct InitializationSpeciesMoments {
    std::size_t initialization_version{1};
    std::string species;
    std::string loading;
    std::string density_profile;
    std::string region;
    std::size_t macroparticles{0};
    double macro_weight{0.0};
    double physical_particles{0.0};
    double represented_charge{0.0};
    double mean_position_x{0.0};
    double mean_position_y{0.0};
    double mean_position_z{0.0};
    double position_stddev_x{0.0};
    double position_stddev_y{0.0};
    double position_stddev_z{0.0};
    double mean_velocity_x{0.0};
    double mean_velocity_y{0.0};
    double mean_velocity_z{0.0};
    double thermal_velocity_x{0.0};
    double thermal_velocity_y{0.0};
    double thermal_velocity_z{0.0};
};

std::string to_string(ParticleLoading loading);
ParticleLoading particle_loading_from_string(const std::string& value);
std::string to_string(DensityProfileKind profile);
DensityProfileKind density_profile_from_string(
    const std::string& value);

void validate_particle_initialization(
    const ParticleInitializationConfig& config,
    std::size_t velocity_dimensions,
    double fallback_thermal_velocity,
    const std::string& context);

void validate_density_profile(
    const ParticleInitializationConfig& config,
    std::size_t spatial_dimensions,
    std::size_t particle_count,
    const std::string& context);

double resolved_thermal_velocity(
    const ParticleInitializationConfig& config,
    std::size_t component,
    double fallback_thermal_velocity);

double quiet_unit_coordinate(
    std::size_t particle_index,
    std::size_t particle_count,
    std::size_t component);

double quiet_sequence_coordinate(
    std::size_t sequence_index,
    std::size_t component);

double density_profile_acceptance(
    const ParticleInitializationConfig& config,
    const std::array<double, 3>& position,
    const std::array<double, 3>& minimum,
    const std::array<double, 3>& maximum);

std::vector<double> initialize_velocity_component(
    std::size_t particle_count,
    double drift_velocity,
    double thermal_velocity,
    ParticleLoading loading,
    std::mt19937_64& rng);

InitializationSpeciesMoments summarize_initialization(
    const Species& species, const std::string& region = {});
InitializationSpeciesMoments summarize_initialization(
    const Species2D& species, const std::string& region = {});
InitializationSpeciesMoments summarize_initialization(
    const Species3D& species, const std::string& region = {});

void write_initialization_report(
    const std::filesystem::path& path,
    std::size_t spatial_dimension,
    const std::string& state_source,
    const std::vector<InitializationSpeciesMoments>& moments);

} // namespace pic
