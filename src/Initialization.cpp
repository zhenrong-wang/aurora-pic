#include "pic/Initialization.hpp"
#include "pic/Species.hpp"
#include "pic/Species2D.hpp"
#include "pic/Species3D.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
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

std::string csv_quote(const std::string& value) {
    std::string result{"\""};
    for (const char character : value) {
        if (character == '"') result += '"';
        result += character;
    }
    result += '"';
    return result;
}

struct VelocityMomentAccumulator {
    std::size_t count{0};
    double mean_x{0.0};
    double mean_y{0.0};
    double mean_z{0.0};
    double moment_x{0.0};
    double moment_y{0.0};
    double moment_z{0.0};

    void add(double x, double y, double z) {
        ++count;
        const double inverse = 1.0 / static_cast<double>(count);
        const double delta_x = x - mean_x;
        const double delta_y = y - mean_y;
        const double delta_z = z - mean_z;
        mean_x += delta_x * inverse;
        mean_y += delta_y * inverse;
        mean_z += delta_z * inverse;
        moment_x += delta_x * (x - mean_x);
        moment_y += delta_y * (y - mean_y);
        moment_z += delta_z * (z - mean_z);
    }
};

InitializationSpeciesMoments finish_moments(
    std::size_t initialization_version,
    std::string species,
    std::string loading,
    std::string region,
    double charge,
    double weight,
    const VelocityMomentAccumulator& accumulator) {
    InitializationSpeciesMoments result;
    result.initialization_version = initialization_version;
    result.species = std::move(species);
    result.loading = std::move(loading);
    result.region = std::move(region);
    result.macroparticles = accumulator.count;
    result.macro_weight = weight;
    result.physical_particles =
        weight * static_cast<double>(accumulator.count);
    result.represented_charge = charge * result.physical_particles;
    if (!std::isfinite(result.physical_particles) ||
        !std::isfinite(result.represented_charge) ||
        !std::isfinite(accumulator.mean_x) ||
        !std::isfinite(accumulator.mean_y) ||
        !std::isfinite(accumulator.mean_z) ||
        !std::isfinite(accumulator.moment_x) ||
        !std::isfinite(accumulator.moment_y) ||
        !std::isfinite(accumulator.moment_z)) {
        throw std::overflow_error(
            "initialization moment accumulation overflow for species '" +
            result.species + "'");
    }
    if (accumulator.count == 0) return result;

    const double inverse =
        1.0 / static_cast<double>(accumulator.count);
    result.mean_velocity_x = accumulator.mean_x;
    result.mean_velocity_y = accumulator.mean_y;
    result.mean_velocity_z = accumulator.mean_z;
    result.thermal_velocity_x = std::sqrt(std::max(
        0.0, accumulator.moment_x * inverse));
    result.thermal_velocity_y = std::sqrt(std::max(
        0.0, accumulator.moment_y * inverse));
    result.thermal_velocity_z = std::sqrt(std::max(
        0.0, accumulator.moment_z * inverse));
    if (!std::isfinite(result.mean_velocity_x) ||
        !std::isfinite(result.mean_velocity_y) ||
        !std::isfinite(result.mean_velocity_z) ||
        !std::isfinite(result.thermal_velocity_x) ||
        !std::isfinite(result.thermal_velocity_y) ||
        !std::isfinite(result.thermal_velocity_z)) {
        throw std::overflow_error(
            "initialization moment reduction overflow for species '" +
            result.species + "'");
    }
    return result;
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

InitializationSpeciesMoments summarize_initialization(
    const Species& species, const std::string& region) {
    VelocityMomentAccumulator accumulator;
    for (const auto& particle : species.particles()) {
        if (!particle.alive) continue;
        accumulator.add(particle.v, 0.0, 0.0);
    }
    return finish_moments(
        species.config().initialization.version, species.name(),
        to_string(species.config().initialization.loading), region,
        species.charge(), species.weight(), accumulator);
}

InitializationSpeciesMoments summarize_initialization(
    const Species2D& species, const std::string& region) {
    VelocityMomentAccumulator accumulator;
    for (const auto& particle : species.particles()) {
        if (!particle.alive) continue;
        accumulator.add(
            particle.velocity.x, particle.velocity.y,
            particle.velocity_z);
    }
    return finish_moments(
        species.config().initialization.version, species.name(),
        to_string(species.config().initialization.loading), region,
        species.charge(), species.weight(), accumulator);
}

InitializationSpeciesMoments summarize_initialization(
    const Species3D& species, const std::string& region) {
    VelocityMomentAccumulator accumulator;
    for (const auto& particle : species.particles()) {
        if (!particle.alive) continue;
        accumulator.add(
            particle.velocity.x, particle.velocity.y,
            particle.velocity.z);
    }
    return finish_moments(
        species.config().initialization.version, species.name(),
        to_string(species.config().initialization.loading), region,
        species.charge(), species.weight(), accumulator);
}

void write_initialization_report(
    const std::filesystem::path& path,
    std::size_t spatial_dimension,
    const std::string& state_source,
    const std::vector<InitializationSpeciesMoments>& moments) {
    if (spatial_dimension == 0 || spatial_dimension > 3) {
        throw std::invalid_argument(
            "initialization report dimension must be 1, 2, or 3");
    }
    if (state_source != "generated" && state_source != "restart") {
        throw std::invalid_argument(
            "initialization report state_source must be generated or restart");
    }
    const auto parent = path.parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error(
            "cannot write initialization report: " + path.string());
    }
    output
        << "initialization_version,state_source,dimension,species,loading,region,"
        << "macroparticles,macro_weight,physical_particles,represented_charge,"
        << "mean_velocity_x,mean_velocity_y,mean_velocity_z,"
        << "thermal_velocity_x,thermal_velocity_y,thermal_velocity_z\n"
        << std::setprecision(17);
    for (const auto& value : moments) {
        output
            << value.initialization_version << ','
            << csv_quote(state_source) << ','
            << spatial_dimension << ','
            << csv_quote(value.species) << ','
            << csv_quote(
                   state_source == "restart" ? "restart" : value.loading)
            << ',' << csv_quote(value.region) << ','
            << value.macroparticles << ','
            << value.macro_weight << ','
            << value.physical_particles << ','
            << value.represented_charge << ','
            << value.mean_velocity_x << ','
            << value.mean_velocity_y << ','
            << value.mean_velocity_z << ','
            << value.thermal_velocity_x << ','
            << value.thermal_velocity_y << ','
            << value.thermal_velocity_z << '\n';
    }
    if (!output) {
        throw std::runtime_error(
            "failed while writing initialization report: " +
            path.string());
    }
}

} // namespace pic
