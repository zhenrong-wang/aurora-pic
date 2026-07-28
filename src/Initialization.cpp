#include "pic/Initialization.hpp"
#include "pic/Species.hpp"
#include "pic/Species2D.hpp"
#include "pic/Species3D.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <numbers>
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

void validate_optional_finite(
    const std::optional<double>& value,
    const std::string& name,
    const std::string& context) {
    if (value && !std::isfinite(*value)) {
        throw std::invalid_argument(
            context + " " + name + " must be finite");
    }
}

double radical_inverse(std::size_t index, std::size_t base) {
    double inverse_base = 1.0 / static_cast<double>(base);
    double factor = inverse_base;
    double result = 0.0;
    while (index != 0) {
        result += factor * static_cast<double>(index % base);
        index /= base;
        factor *= inverse_base;
    }
    return result;
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
    std::string density_profile,
    std::string region,
    double charge,
    double weight,
    const VelocityMomentAccumulator& position_accumulator,
    const VelocityMomentAccumulator& velocity_accumulator) {
    InitializationSpeciesMoments result;
    result.initialization_version = initialization_version;
    result.species = std::move(species);
    result.loading = std::move(loading);
    result.density_profile = std::move(density_profile);
    result.region = std::move(region);
    if (position_accumulator.count !=
        velocity_accumulator.count) {
        throw std::logic_error(
            "initialization position/velocity moment counts differ");
    }
    result.macroparticles = velocity_accumulator.count;
    result.macro_weight = weight;
    result.physical_particles =
        weight * static_cast<double>(velocity_accumulator.count);
    result.represented_charge = charge * result.physical_particles;
    if (!std::isfinite(result.physical_particles) ||
        !std::isfinite(result.represented_charge) ||
        !std::isfinite(position_accumulator.mean_x) ||
        !std::isfinite(position_accumulator.mean_y) ||
        !std::isfinite(position_accumulator.mean_z) ||
        !std::isfinite(position_accumulator.moment_x) ||
        !std::isfinite(position_accumulator.moment_y) ||
        !std::isfinite(position_accumulator.moment_z) ||
        !std::isfinite(velocity_accumulator.mean_x) ||
        !std::isfinite(velocity_accumulator.mean_y) ||
        !std::isfinite(velocity_accumulator.mean_z) ||
        !std::isfinite(velocity_accumulator.moment_x) ||
        !std::isfinite(velocity_accumulator.moment_y) ||
        !std::isfinite(velocity_accumulator.moment_z)) {
        throw std::overflow_error(
            "initialization moment accumulation overflow for species '" +
            result.species + "'");
    }
    if (velocity_accumulator.count == 0) return result;

    const double inverse =
        1.0 / static_cast<double>(velocity_accumulator.count);
    result.mean_position_x = position_accumulator.mean_x;
    result.mean_position_y = position_accumulator.mean_y;
    result.mean_position_z = position_accumulator.mean_z;
    result.position_stddev_x = std::sqrt(std::max(
        0.0, position_accumulator.moment_x * inverse));
    result.position_stddev_y = std::sqrt(std::max(
        0.0, position_accumulator.moment_y * inverse));
    result.position_stddev_z = std::sqrt(std::max(
        0.0, position_accumulator.moment_z * inverse));
    result.mean_velocity_x = velocity_accumulator.mean_x;
    result.mean_velocity_y = velocity_accumulator.mean_y;
    result.mean_velocity_z = velocity_accumulator.mean_z;
    result.thermal_velocity_x = std::sqrt(std::max(
        0.0, velocity_accumulator.moment_x * inverse));
    result.thermal_velocity_y = std::sqrt(std::max(
        0.0, velocity_accumulator.moment_y * inverse));
    result.thermal_velocity_z = std::sqrt(std::max(
        0.0, velocity_accumulator.moment_z * inverse));
    if (!std::isfinite(result.mean_position_x) ||
        !std::isfinite(result.mean_position_y) ||
        !std::isfinite(result.mean_position_z) ||
        !std::isfinite(result.position_stddev_x) ||
        !std::isfinite(result.position_stddev_y) ||
        !std::isfinite(result.position_stddev_z) ||
        !std::isfinite(result.mean_velocity_x) ||
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

std::string to_string(DensityProfileKind profile) {
    switch (profile) {
        case DensityProfileKind::Uniform: return "uniform";
        case DensityProfileKind::Gaussian: return "gaussian";
        case DensityProfileKind::Sinusoidal: return "sinusoidal";
    }
    return "unknown";
}

DensityProfileKind density_profile_from_string(
    const std::string& value) {
    if (value == "uniform") return DensityProfileKind::Uniform;
    if (value == "gaussian") return DensityProfileKind::Gaussian;
    if (value == "sinusoidal" || value == "sinusoid") {
        return DensityProfileKind::Sinusoidal;
    }
    throw std::invalid_argument(
        "density_profile must be uniform, gaussian, or sinusoidal");
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

void validate_density_profile(
    const ParticleInitializationConfig& config,
    std::size_t spatial_dimensions,
    std::size_t particle_count,
    const std::string& context) {
    if (spatial_dimensions == 0 || spatial_dimensions > 3) {
        throw std::invalid_argument(
            context + " has an invalid spatial dimension");
    }
    if (config.max_profile_sampling_attempts == 0) {
        throw std::invalid_argument(
            context +
            " max_profile_sampling_attempts must be positive");
    }

    const std::array<std::optional<double>, 3> centers{
        config.profile_center_x, config.profile_center_y,
        config.profile_center_z};
    const std::array<std::optional<double>, 3> scales{
        config.profile_scale_x, config.profile_scale_y,
        config.profile_scale_z};
    const std::array<std::optional<std::size_t>, 3> modes{
        config.profile_mode_x, config.profile_mode_y,
        config.profile_mode_z};
    for (std::size_t component = 0; component < 3; ++component) {
        validate_optional_finite(
            centers[component],
            "profile_center_" +
                std::string(1, "xyz"[component]),
            context);
        validate_optional_finite(
            scales[component],
            "profile_scale_" +
                std::string(1, "xyz"[component]),
            context);
        if (component >= spatial_dimensions &&
            (centers[component] || scales[component] ||
             modes[component])) {
            throw std::invalid_argument(
                context +
                " density profile configures an inactive spatial component");
        }
    }
    validate_optional_finite(
        config.profile_amplitude, "profile_amplitude", context);
    validate_optional_finite(
        config.profile_phase, "profile_phase", context);

    if (config.density_profile == DensityProfileKind::Uniform) {
        if (config.profile_center_x || config.profile_center_y ||
            config.profile_center_z || config.profile_scale_x ||
            config.profile_scale_y || config.profile_scale_z ||
            config.profile_amplitude || config.profile_phase ||
            config.profile_mode_x || config.profile_mode_y ||
            config.profile_mode_z) {
            throw std::invalid_argument(
                context +
                " uniform density_profile does not accept profile parameters");
        }
        return;
    }
    if (config.max_profile_sampling_attempts < particle_count) {
        throw std::invalid_argument(
            context +
            " max_profile_sampling_attempts must be at least particles");
    }

    if (config.density_profile == DensityProfileKind::Gaussian) {
        if (config.profile_amplitude || config.profile_phase ||
            config.profile_mode_x || config.profile_mode_y ||
            config.profile_mode_z) {
            throw std::invalid_argument(
                context +
                " gaussian density_profile does not accept sinusoidal parameters");
        }
        for (std::size_t component = 0;
             component < spatial_dimensions; ++component) {
            if (!centers[component] || !scales[component] ||
                !(*scales[component] > 0.0)) {
                throw std::invalid_argument(
                    context +
                    " gaussian density_profile requires finite centers and positive scales for every active axis");
            }
        }
        return;
    }

    if (config.profile_center_x || config.profile_center_y ||
        config.profile_center_z || config.profile_scale_x ||
        config.profile_scale_y || config.profile_scale_z) {
        throw std::invalid_argument(
            context +
            " sinusoidal density_profile does not accept gaussian parameters");
    }
    if (!config.profile_amplitude ||
        std::abs(*config.profile_amplitude) > 1.0) {
        throw std::invalid_argument(
            context +
            " sinusoidal density_profile requires profile_amplitude with absolute value at most one");
    }
    bool has_mode = false;
    for (std::size_t component = 0;
         component < spatial_dimensions; ++component) {
        has_mode = has_mode ||
                   modes[component].value_or(0) != 0;
    }
    if (!has_mode) {
        throw std::invalid_argument(
            context +
            " sinusoidal density_profile requires a nonzero active profile mode");
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

double quiet_sequence_coordinate(
    std::size_t sequence_index,
    std::size_t component) {
    static constexpr std::array<std::size_t, 7> primes{
        2, 3, 5, 7, 11, 13, 17};
    if (component >= primes.size()) {
        throw std::invalid_argument(
            "quiet sequence component exceeds supported dimensions");
    }
    if (sequence_index == std::numeric_limits<std::size_t>::max()) {
        throw std::overflow_error("quiet sequence index overflow");
    }
    return radical_inverse(
        sequence_index + 1, primes[component]);
}

double density_profile_acceptance(
    const ParticleInitializationConfig& config,
    const std::array<double, 3>& position,
    const std::array<double, 3>& minimum,
    const std::array<double, 3>& maximum) {
    if (config.density_profile == DensityProfileKind::Uniform) {
        return 1.0;
    }
    if (config.density_profile == DensityProfileKind::Gaussian) {
        const std::array<std::optional<double>, 3> centers{
            config.profile_center_x, config.profile_center_y,
            config.profile_center_z};
        const std::array<std::optional<double>, 3> scales{
            config.profile_scale_x, config.profile_scale_y,
            config.profile_scale_z};
        double exponent = 0.0;
        for (std::size_t component = 0; component < 3; ++component) {
            if (centers[component].has_value() !=
                scales[component].has_value()) {
                throw std::invalid_argument(
                    "gaussian density profile has incomplete center/scale parameters");
            }
            if (!centers[component]) continue;
            if (!(*scales[component] > 0.0) ||
                !std::isfinite(*scales[component])) {
                throw std::invalid_argument(
                    "gaussian density profile scale must be positive and finite");
            }
            const double normalized =
                (position[component] - *centers[component]) /
                *scales[component];
            exponent += normalized * normalized;
        }
        return std::exp(-0.5 * exponent);
    }

    const std::array<std::optional<std::size_t>, 3> modes{
        config.profile_mode_x, config.profile_mode_y,
        config.profile_mode_z};
    double argument = config.profile_phase.value_or(0.0);
    for (std::size_t component = 0; component < 3; ++component) {
        const std::size_t mode = modes[component].value_or(0);
        if (mode == 0) continue;
        const double length = maximum[component] - minimum[component];
        if (!(length > 0.0) || !std::isfinite(length)) {
            throw std::invalid_argument(
                "density profile envelope must have positive finite active lengths");
        }
        const double normalized =
            (position[component] - minimum[component]) / length;
        argument += 2.0 * std::numbers::pi *
                    static_cast<double>(mode) * normalized;
    }
    if (!config.profile_amplitude ||
        !std::isfinite(*config.profile_amplitude) ||
        std::abs(*config.profile_amplitude) > 1.0) {
        throw std::invalid_argument(
            "sinusoidal density profile amplitude is missing or invalid");
    }
    const double amplitude = *config.profile_amplitude;
    const double result =
        (1.0 + amplitude * std::cos(argument)) /
        (1.0 + std::abs(amplitude));
    return std::clamp(result, 0.0, 1.0);
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
    VelocityMomentAccumulator position_accumulator;
    VelocityMomentAccumulator velocity_accumulator;
    for (const auto& particle : species.particles()) {
        if (!particle.alive) continue;
        position_accumulator.add(particle.x, 0.0, 0.0);
        velocity_accumulator.add(particle.v, 0.0, 0.0);
    }
    return finish_moments(
        species.config().initialization.version, species.name(),
        to_string(species.config().initialization.loading),
        to_string(species.config().initialization.density_profile),
        region,
        species.charge(), species.weight(),
        position_accumulator, velocity_accumulator);
}

InitializationSpeciesMoments summarize_initialization(
    const Species2D& species, const std::string& region) {
    VelocityMomentAccumulator position_accumulator;
    VelocityMomentAccumulator velocity_accumulator;
    for (const auto& particle : species.particles()) {
        if (!particle.alive) continue;
        position_accumulator.add(
            particle.position.x, particle.position.y, 0.0);
        velocity_accumulator.add(
            particle.velocity.x, particle.velocity.y,
            particle.velocity_z);
    }
    return finish_moments(
        species.config().initialization.version, species.name(),
        to_string(species.config().initialization.loading),
        to_string(species.config().initialization.density_profile),
        region,
        species.charge(), species.weight(),
        position_accumulator, velocity_accumulator);
}

InitializationSpeciesMoments summarize_initialization(
    const Species3D& species, const std::string& region) {
    VelocityMomentAccumulator position_accumulator;
    VelocityMomentAccumulator velocity_accumulator;
    for (const auto& particle : species.particles()) {
        if (!particle.alive) continue;
        position_accumulator.add(
            particle.position.x, particle.position.y,
            particle.position.z);
        velocity_accumulator.add(
            particle.velocity.x, particle.velocity.y,
            particle.velocity.z);
    }
    return finish_moments(
        species.config().initialization.version, species.name(),
        to_string(species.config().initialization.loading),
        to_string(species.config().initialization.density_profile),
        region,
        species.charge(), species.weight(),
        position_accumulator, velocity_accumulator);
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
    if (state_source != "generated" &&
        state_source != "restart" &&
        state_source != "external") {
        throw std::invalid_argument(
            "initialization report state_source must be generated, external, or restart");
    }
    const auto parent = path.parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error(
            "cannot write initialization report: " + path.string());
    }
    output
        << "initialization_version,state_source,dimension,species,loading,density_profile,region,"
        << "macroparticles,macro_weight,physical_particles,represented_charge,"
        << "mean_position_x,mean_position_y,mean_position_z,"
        << "position_stddev_x,position_stddev_y,position_stddev_z,"
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
                   state_source == "generated"
                       ? value.loading
                       : state_source)
            << ',' << csv_quote(
                   state_source == "generated"
                       ? value.density_profile
                       : state_source)
            << ',' << csv_quote(value.region) << ','
            << value.macroparticles << ','
            << value.macro_weight << ','
            << value.physical_particles << ','
            << value.represented_charge << ','
            << value.mean_position_x << ','
            << value.mean_position_y << ','
            << value.mean_position_z << ','
            << value.position_stddev_x << ','
            << value.position_stddev_y << ','
            << value.position_stddev_z << ','
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

void validate_initialization_acceptance(
    const InitializationAcceptanceConfig& config,
    const std::string& context) {
    const auto validate_tolerance =
        [&](const std::optional<double>& value,
            const std::string& name) {
            if (value &&
                (!std::isfinite(*value) || *value < 0.0 ||
                 *value > 1.0)) {
                throw std::invalid_argument(
                    context + " " + name +
                    " must be finite and in [0, 1]");
            }
        };
    validate_tolerance(
        config.max_relative_charge_imbalance,
        "max_relative_charge_imbalance");
    validate_tolerance(
        config.max_relative_current_imbalance,
        "max_relative_current_imbalance");
    validate_tolerance(
        config.max_relative_pair_imbalance,
        "max_relative_pair_imbalance");
    if (config.charge_pairs.empty() &&
        config.max_relative_pair_imbalance) {
        throw std::invalid_argument(
            context +
            " max_relative_pair_imbalance requires charge_pairs");
    }
    if (!config.charge_pairs.empty() &&
        !config.max_relative_pair_imbalance) {
        throw std::invalid_argument(
            context +
            " charge_pairs require max_relative_pair_imbalance");
    }
    std::map<std::string, std::string> paired_species;
    for (const auto& pair : config.charge_pairs) {
        if (pair.first_species.empty() ||
            pair.second_species.empty()) {
            throw std::invalid_argument(
                context + " charge-pair species names cannot be empty");
        }
        if (pair.first_species == pair.second_species) {
            throw std::invalid_argument(
                context +
                " charge-pair species names must differ");
        }
        for (const auto& [first, second] :
             {std::pair{pair.first_species, pair.second_species},
              std::pair{pair.second_species, pair.first_species}}) {
            if (!paired_species.emplace(first, second).second) {
                throw std::invalid_argument(
                    context + " species '" + first +
                    "' occurs in more than one charge pair");
            }
        }
    }
}

InitializationAcceptanceSummary assess_initialization_acceptance(
    const InitializationAcceptanceConfig& config,
    const std::vector<InitializationSpeciesMoments>& moments,
    std::size_t velocity_dimensions) {
    validate_initialization_acceptance(
        config, "initialization acceptance");
    if (velocity_dimensions == 0 || velocity_dimensions > 3) {
        throw std::invalid_argument(
            "initialization acceptance velocity dimension must be 1, 2, or 3");
    }
    InitializationAcceptanceSummary summary;
    summary.enabled =
        config.max_relative_charge_imbalance.has_value() ||
        config.max_relative_current_imbalance.has_value() ||
        !config.charge_pairs.empty();
    if (summary.enabled && moments.empty()) {
        throw std::invalid_argument(
            "enabled initialization acceptance gates require at least one species");
    }
    std::map<std::string, const InitializationSpeciesMoments*>
        species_by_name;
    for (const auto& species : moments) {
        if (species.species.empty() ||
            !species_by_name.emplace(
                species.species, &species).second) {
            throw std::invalid_argument(
                "initialization acceptance requires unique non-empty species names");
        }
        if (!std::isfinite(species.represented_charge) ||
            !std::isfinite(species.mean_velocity_x) ||
            !std::isfinite(species.mean_velocity_y) ||
            !std::isfinite(species.mean_velocity_z)) {
            throw std::invalid_argument(
                "initialization acceptance received non-finite moments for species '" +
                species.species + "'");
        }
    }
    const auto residual = [](double value, double scale) {
        return scale > 0.0 ? value / scale
                           : (value == 0.0 ? 0.0 : 1.0);
    };
    if (config.max_relative_charge_imbalance) {
        double net_charge = 0.0;
        double charge_scale = 0.0;
        for (const auto& species : moments) {
            net_charge += species.represented_charge;
            charge_scale += std::abs(species.represented_charge);
        }
        if (!std::isfinite(net_charge) ||
            !std::isfinite(charge_scale)) {
            throw std::overflow_error(
                "initialization net-charge acceptance accumulation overflowed");
        }
        const double value = std::abs(net_charge);
        const double relative = residual(value, charge_scale);
        const bool passed =
            relative <= *config.max_relative_charge_imbalance;
        summary.metrics.push_back({
            "net_charge", net_charge, charge_scale, relative,
            *config.max_relative_charge_imbalance, passed,
            "value is signed represented charge; residual uses its magnitude"});
        summary.passed = summary.passed && passed;
    }
    if (config.max_relative_current_imbalance) {
        std::array<double, 3> net_current{};
        double current_scale = 0.0;
        for (const auto& species : moments) {
            const std::array<double, 3> velocity{
                species.mean_velocity_x,
                species.mean_velocity_y,
                species.mean_velocity_z};
            double speed = 0.0;
            for (std::size_t component = 0;
                 component < velocity_dimensions; ++component) {
                net_current[component] +=
                    species.represented_charge *
                    velocity[component];
                speed = std::hypot(
                    speed, velocity[component]);
            }
            current_scale +=
                std::abs(species.represented_charge) *
                speed;
        }
        if (!std::isfinite(current_scale) ||
            !std::all_of(
                net_current.begin(), net_current.end(),
                [](double value) {
                    return std::isfinite(value);
                })) {
            throw std::overflow_error(
                "initialization net-current acceptance accumulation overflowed");
        }
        double value = 0.0;
        for (std::size_t component = 0;
             component < velocity_dimensions; ++component) {
            value = std::hypot(
                value, net_current[component]);
        }
        const double relative = residual(value, current_scale);
        const bool passed =
            relative <= *config.max_relative_current_imbalance;
        summary.metrics.push_back({
            "net_current", value, current_scale, relative,
            *config.max_relative_current_imbalance, passed,
            "charge-weighted mean-velocity norm"});
        summary.passed = summary.passed && passed;
    }
    for (const auto& pair : config.charge_pairs) {
        const auto first = species_by_name.find(pair.first_species);
        const auto second = species_by_name.find(pair.second_species);
        if (first == species_by_name.end() ||
            second == species_by_name.end()) {
            throw std::invalid_argument(
                "initialization charge pair references unknown species '" +
                (first == species_by_name.end()
                     ? pair.first_species
                     : pair.second_species) +
                "'");
        }
        const double first_charge =
            first->second->represented_charge;
        const double second_charge =
            second->second->represented_charge;
        const double value =
            std::abs(std::abs(first_charge) -
                     std::abs(second_charge));
        const double scale =
            std::max(std::abs(first_charge),
                     std::abs(second_charge));
        const double relative = residual(value, scale);
        const bool opposite =
            first_charge != 0.0 && second_charge != 0.0 &&
            std::signbit(first_charge) !=
                std::signbit(second_charge);
        const bool passed =
            opposite &&
            relative <=
                *config.max_relative_pair_imbalance;
        summary.metrics.push_back({
            "charge_pair:" + pair.first_species + ":" +
                pair.second_species,
            value, scale, relative,
            *config.max_relative_pair_imbalance, passed,
            opposite
                ? "opposite-sign represented-charge magnitude balance"
                : "species do not have opposite represented-charge signs"});
        summary.passed = summary.passed && passed;
    }
    return summary;
}

void write_initialization_acceptance_report(
    const std::filesystem::path& path,
    const InitializationAcceptanceSummary& summary) {
    const auto parent = path.parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error(
            "cannot write initialization acceptance report: " +
            path.string());
    }
    output
        << "enabled,overall_passed,metric,value,scale,"
        << "relative_residual,tolerance,passed,details\n"
        << std::setprecision(17);
    if (summary.metrics.empty()) {
        output << (summary.enabled ? 1 : 0) << ','
               << (summary.passed ? 1 : 0)
               << ",\"none\",0,0,0,0,1,"
               << csv_quote(summary.enabled
                                ? "no acceptance metrics"
                                : "acceptance gates disabled")
               << '\n';
    } else {
        for (const auto& metric : summary.metrics) {
            output << (summary.enabled ? 1 : 0) << ','
                   << (summary.passed ? 1 : 0) << ','
                   << csv_quote(metric.metric) << ','
                   << metric.value << ','
                   << metric.scale << ','
                   << metric.relative_residual << ','
                   << metric.tolerance << ','
                   << (metric.passed ? 1 : 0) << ','
                   << csv_quote(metric.details) << '\n';
        }
    }
    if (!output) {
        throw std::runtime_error(
            "failed while writing initialization acceptance report: " +
            path.string());
    }
}

void enforce_initialization_acceptance(
    const InitializationAcceptanceSummary& summary) {
    if (summary.passed) return;
    std::string failed_metrics;
    for (const auto& metric : summary.metrics) {
        if (metric.passed) continue;
        if (!failed_metrics.empty()) failed_metrics += ", ";
        failed_metrics += metric.metric;
    }
    throw std::runtime_error(
        "initialization acceptance gates failed: " +
        failed_metrics);
}

} // namespace pic
