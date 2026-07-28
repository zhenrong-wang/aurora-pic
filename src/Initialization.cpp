#include "pic/Initialization.hpp"

#include <cmath>
#include <numeric>
#include <stdexcept>

namespace pic {
namespace {

std::size_t quiet_multiplier(std::size_t count, std::size_t component) {
    if (count <= 1 || component == 0) return 1;
    std::size_t candidate = 2 * component + 1;
    while (std::gcd(candidate, count) != 1) ++candidate;
    return candidate;
}

void validate_thermal_velocity(
    const std::optional<double>& value,
    const std::string& name,
    const std::string& context) {
    if (value && (!std::isfinite(*value) || *value < 0.0)) {
        throw std::invalid_argument(
            context + " " + name + " must be non-negative and finite");
    }
}

} // namespace

std::string to_string(ParticleLoading loading) {
    switch (loading) {
        case ParticleLoading::Random: return "random";
        case ParticleLoading::QuietStart: return "quiet_start";
    }
    return "unknown";
}

ParticleLoading particle_loading_from_string(const std::string& value) {
    if (value == "random") return ParticleLoading::Random;
    if (value == "quiet_start" || value == "quiet") {
        return ParticleLoading::QuietStart;
    }
    throw std::invalid_argument(
        "particle loading must be random or quiet_start");
}

void validate_particle_initialization(
    const ParticleInitializationConfig& config,
    std::size_t velocity_dimensions,
    double fallback_thermal_velocity,
    const std::string& context) {
    if (config.version != 1) {
        throw std::invalid_argument(
            context + " declares unsupported initialization_version " +
            std::to_string(config.version) +
            "; this AuroraPIC build supports initialization_version = 1");
    }
    if (velocity_dimensions == 0 || velocity_dimensions > 3) {
        throw std::invalid_argument(
            context + " has an invalid velocity dimension");
    }
    if (!std::isfinite(fallback_thermal_velocity) ||
        fallback_thermal_velocity < 0.0) {
        throw std::invalid_argument(
            context + " thermal_velocity must be non-negative and finite");
    }
    validate_thermal_velocity(
        config.thermal_velocity_x, "thermal_velocity_x", context);
    validate_thermal_velocity(
        config.thermal_velocity_y, "thermal_velocity_y", context);
    validate_thermal_velocity(
        config.thermal_velocity_z, "thermal_velocity_z", context);
    if (velocity_dimensions < 2 && config.thermal_velocity_y) {
        throw std::invalid_argument(
            context + " thermal_velocity_y is unavailable in 1D1V");
    }
    if (velocity_dimensions < 3 && config.thermal_velocity_z) {
        throw std::invalid_argument(
            context + " thermal_velocity_z is unavailable for this velocity model");
    }
}

double resolved_thermal_velocity(
    const ParticleInitializationConfig& config,
    std::size_t component,
    double fallback_thermal_velocity) {
    const std::optional<double>* configured = nullptr;
    switch (component) {
        case 0: configured = &config.thermal_velocity_x; break;
        case 1: configured = &config.thermal_velocity_y; break;
        case 2: configured = &config.thermal_velocity_z; break;
        default:
            throw std::invalid_argument(
                "thermal-velocity component must be 0, 1, or 2");
    }
    return configured->value_or(fallback_thermal_velocity);
}

double quiet_unit_coordinate(
    std::size_t particle_index,
    std::size_t particle_count,
    std::size_t component) {
    if (particle_count == 0 || particle_index >= particle_count) {
        throw std::invalid_argument(
            "quiet-start coordinate index is outside the population");
    }
    const std::size_t multiplier =
        quiet_multiplier(particle_count, component);
    const std::size_t offset =
        component == 0 ? 0 : (particle_count / (component + 1));
    const std::size_t stratum =
        (multiplier * particle_index + offset) % particle_count;
    return (static_cast<double>(stratum) + 0.5) /
           static_cast<double>(particle_count);
}

std::vector<double> initialize_velocity_component(
    std::size_t particle_count,
    double drift_velocity,
    double thermal_velocity,
    ParticleLoading loading,
    std::mt19937_64& rng) {
    if (particle_count == 0) return {};
    if (!std::isfinite(drift_velocity)) {
        throw std::invalid_argument(
            "initial drift velocity must be finite");
    }
    if (!std::isfinite(thermal_velocity) || thermal_velocity < 0.0) {
        throw std::invalid_argument(
            "initial thermal velocity must be non-negative and finite");
    }

    std::vector<double> velocities(particle_count, drift_velocity);
    if (thermal_velocity == 0.0) return velocities;

    std::normal_distribution<double> normal(0.0, thermal_velocity);
    if (loading == ParticleLoading::Random) {
        for (double& velocity : velocities) {
            velocity += normal(rng);
        }
        return velocities;
    }

    const std::size_t pairs = particle_count / 2;
    for (std::size_t pair = 0; pair < pairs; ++pair) {
        const double perturbation = normal(rng);
        velocities[2 * pair] += perturbation;
        velocities[2 * pair + 1] -= perturbation;
    }
    return velocities;
}

} // namespace pic
