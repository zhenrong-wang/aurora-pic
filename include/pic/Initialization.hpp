#pragma once

#include <cstddef>
#include <optional>
#include <random>
#include <string>
#include <vector>

namespace pic {

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

} // namespace pic
