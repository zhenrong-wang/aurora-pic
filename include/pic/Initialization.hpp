#pragma once

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

struct ParticleInitializationConfig {
    std::size_t version{1};
    ParticleLoading loading{ParticleLoading::Random};
    std::optional<double> thermal_velocity_x{};
    std::optional<double> thermal_velocity_y{};
    std::optional<double> thermal_velocity_z{};
};

struct InitializationSpeciesMoments {
    std::size_t initialization_version{1};
    std::string species;
    std::string loading;
    std::string region;
    std::size_t macroparticles{0};
    double macro_weight{0.0};
    double physical_particles{0.0};
    double represented_charge{0.0};
    double mean_velocity_x{0.0};
    double mean_velocity_y{0.0};
    double mean_velocity_z{0.0};
    double thermal_velocity_x{0.0};
    double thermal_velocity_y{0.0};
    double thermal_velocity_z{0.0};
};

std::string to_string(ParticleLoading loading);
ParticleLoading particle_loading_from_string(const std::string& value);

void validate_particle_initialization(
    const ParticleInitializationConfig& config,
    std::size_t velocity_dimensions,
    double fallback_thermal_velocity,
    const std::string& context);

double resolved_thermal_velocity(
    const ParticleInitializationConfig& config,
    std::size_t component,
    double fallback_thermal_velocity);

double quiet_unit_coordinate(
    std::size_t particle_index,
    std::size_t particle_count,
    std::size_t component);

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
