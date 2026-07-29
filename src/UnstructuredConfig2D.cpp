#include "pic/UnstructuredSimulation2D.hpp"
#include "pic/GasDataset.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace pic {
namespace {

using Values = std::map<std::string, std::string>;

struct NamedBlock {
    std::string name;
    Values values;
};

struct ParsedConfig {
    Values global;
    Values collisions;
    std::vector<NamedBlock> boundaries;
    std::vector<NamedBlock> species;
    std::vector<NamedBlock> sources;
    std::vector<NamedBlock> emissions;
    std::vector<NamedBlock> collision_channels;
};

std::string trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char character) {
                       return static_cast<char>(std::tolower(character));
                   });
    return value;
}

[[noreturn]] void config_error(std::size_t line, const std::string& message) {
    throw std::runtime_error("unstructured 2D config line " +
                             std::to_string(line) + ": " + message);
}

void assign(Values& values, std::string key, std::string value,
            std::size_t line) {
    key = lower(trim(std::move(key)));
    value = trim(std::move(value));
    if (key.empty() || value.empty()) config_error(line, "empty key or value");
    if (!values.emplace(key, value).second) {
        config_error(line, "duplicate key '" + key + "'");
    }
}

ParsedConfig parse(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open config: " + path.string());

    static const std::set<std::string> global_keys{
        "config_version", "dimension", "mesh", "mesh_file", "dt", "steps",
        "units", "relative_permittivity",
        "mode", "steady_tolerance", "steady_window", "max_steps",
        "max_particles_per_species", "seed",
        "magnetic_field_x", "magnetic_field_y", "magnetic_field_z",
        "output_interval", "output_dir", "vtk_output",
        "particle_output", "particle_output_interval", "particle_output_stride",
        "particle_sample_count", "checkpoint_output", "checkpoint_interval",
        "checkpoint_path", "restart_path", "initial_state_path",
        "initial_state_signature",
        "runtime_backend", "runtime_threads",
        "poisson_relative_tolerance", "poisson_absolute_tolerance",
        "poisson_max_iterations",
        "initialization_max_relative_charge_imbalance",
        "initialization_max_relative_current_imbalance",
        "initialization_max_relative_pair_imbalance",
        "initialization_charge_pairs",
    };
    static const std::set<std::string> boundary_keys{
        "field", "potential", "normal_derivative", "particle"};
    static const std::set<std::string> species_keys{
        "charge", "mass", "weight", "particles", "drift_velocity_x",
        "drift_velocity_y", "drift_velocity_z", "thermal_velocity",
        "thermal_velocity_x", "thermal_velocity_y", "thermal_velocity_z",
        "initialization_version", "loading",
        "initialization_region",
        "density_profile", "profile_center_x", "profile_center_y",
        "profile_scale_x", "profile_scale_y", "profile_amplitude",
        "profile_phase", "profile_mode_x", "profile_mode_y",
        "max_profile_sampling_attempts",
        "init_x_min", "init_x_max",
        "init_y_min", "init_y_max",
    };
    static const std::set<std::string> source_keys{
        "species", "boundary", "particles_per_step", "start_step", "end_step",
        "normal_velocity", "tangential_velocity", "thermal_velocity",
        "out_of_plane_velocity",
    };
    static const std::set<std::string> emission_keys{
        "boundary", "incident_species", "emitted_species", "yield",
        "max_particles_per_impact", "normal_velocity", "tangential_velocity",
        "thermal_velocity", "out_of_plane_velocity",
    };
    static const std::set<std::string> collision_keys{
        "enabled", "model", "species", "gas", "neutral_density",
        "neutral_mass", "neutral_temperature", "data_provenance",
        "gas_data_file", "max_frequency", "max_candidates_per_particle",
    };
    static const std::set<std::string> collision_channel_keys{
        "type", "cross_section_file", "threshold_energy",
        "energy_scale", "cross_section_scale",
        "angular_model", "mean_cosine_file",
        "mean_cosine_energy_scale",
        "secondary_species", "ion_species", "attachment_species",
    };

    ParsedConfig result;
    Values* current = &result.global;
    const std::set<std::string>* allowed = &global_keys;
    std::set<std::string> boundary_names;
    std::set<std::string> species_names;
    std::set<std::string> source_names;
    std::set<std::string> emission_names;
    std::set<std::string> collision_channel_names;
    bool have_collisions = false;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        const auto comment = line.find_first_of("#;");
        if (comment != std::string::npos) line.resize(comment);
        line = trim(line);
        if (line.empty()) continue;
        if (line.front() == '[' && line.back() == ']') {
            const std::string section = trim(line.substr(1, line.size() - 2));
            const std::string lowered = lower(section);
            constexpr const char* boundary_prefix = "boundary.";
            constexpr const char* species_prefix = "species.";
            constexpr const char* source_prefix = "source.";
            constexpr const char* emission_prefix = "emission.";
            constexpr const char* collision_prefix = "collision.";
            if (lowered == "collisions") {
                if (have_collisions) {
                    config_error(
                        line_number, "duplicate collisions section");
                }
                have_collisions = true;
                current = &result.collisions;
                allowed = &collision_keys;
            } else if (lowered.rfind(boundary_prefix, 0) == 0) {
                const std::string name =
                    trim(section.substr(std::char_traits<char>::length(boundary_prefix)));
                if (name.empty() || !boundary_names.insert(name).second) {
                    config_error(line_number, "empty or duplicate boundary section");
                }
                result.boundaries.push_back({name, {}});
                current = &result.boundaries.back().values;
                allowed = &boundary_keys;
            } else if (lowered.rfind(species_prefix, 0) == 0) {
                const std::string name =
                    trim(section.substr(std::char_traits<char>::length(species_prefix)));
                if (name.empty() || !species_names.insert(name).second) {
                    config_error(line_number, "empty or duplicate species section");
                }
                result.species.push_back({name, {}});
                current = &result.species.back().values;
                allowed = &species_keys;
            } else if (lowered.rfind(source_prefix, 0) == 0) {
                const std::string name =
                    trim(section.substr(std::char_traits<char>::length(source_prefix)));
                if (name.empty() || !source_names.insert(name).second) {
                    config_error(line_number, "empty or duplicate source section");
                }
                result.sources.push_back({name, {}});
                current = &result.sources.back().values;
                allowed = &source_keys;
            } else if (lowered.rfind(emission_prefix, 0) == 0) {
                const std::string name =
                    trim(section.substr(
                        std::char_traits<char>::length(emission_prefix)));
                if (name.empty() || !emission_names.insert(name).second) {
                    config_error(
                        line_number, "empty or duplicate emission section");
                }
                result.emissions.push_back({name, {}});
                current = &result.emissions.back().values;
                allowed = &emission_keys;
            } else if (lowered.rfind(collision_prefix, 0) == 0) {
                const std::string name =
                    trim(section.substr(
                        std::char_traits<char>::length(collision_prefix)));
                if (name.empty() ||
                    !collision_channel_names.insert(name).second) {
                    config_error(
                        line_number,
                        "empty or duplicate collision channel section");
                }
                result.collision_channels.push_back({name, {}});
                current = &result.collision_channels.back().values;
                allowed = &collision_channel_keys;
            } else {
                config_error(line_number, "unknown section '" + section + "'");
            }
            continue;
        }
        const auto separator = line.find('=');
        if (separator == std::string::npos) {
            config_error(line_number, "expected key = value");
        }
        const std::string key = lower(trim(line.substr(0, separator)));
        if (!allowed->contains(key)) {
            config_error(line_number, "unknown key '" + key + "'");
        }
        assign(*current, key, line.substr(separator + 1), line_number);
    }
    return result;
}

const std::string& required(const Values& values, const std::string& key,
                            const std::string& context) {
    const auto it = values.find(key);
    if (it == values.end()) {
        throw std::runtime_error("unstructured 2D config missing " + context +
                                 " key '" + key + "'");
    }
    return it->second;
}

template <typename T>
T number(const Values& values, const std::string& key, T fallback) {
    const auto it = values.find(key);
    if (it == values.end()) return fallback;
    std::size_t consumed = 0;
    try {
        if constexpr (std::is_floating_point_v<T>) {
            const double parsed = std::stod(it->second, &consumed);
            if (consumed != it->second.size() || !std::isfinite(parsed)) {
                throw std::invalid_argument("non-finite or trailing");
            }
            return static_cast<T>(parsed);
        } else {
            if (it->second.empty() || it->second.front() == '-') {
                throw std::invalid_argument("negative");
            }
            const unsigned long long parsed = std::stoull(it->second, &consumed);
            if (consumed != it->second.size() ||
                parsed > std::numeric_limits<T>::max()) {
                throw std::invalid_argument("range");
            }
            return static_cast<T>(parsed);
        }
    } catch (const std::exception&) {
        throw std::runtime_error("invalid numeric value for unstructured key '" +
                                 key + "': " + it->second);
    }
}

bool boolean(const Values& values, const std::string& key, bool fallback) {
    const auto it = values.find(key);
    if (it == values.end()) return fallback;
    const std::string value = lower(it->second);
    if (value == "true" || value == "yes" || value == "1") return true;
    if (value == "false" || value == "no" || value == "0") return false;
    throw std::runtime_error("invalid boolean value for unstructured key '" +
                             key + "': " + it->second);
}

std::optional<std::uint64_t> optional_uint64(
    const Values& values, const std::string& key) {
    const auto found = values.find(key);
    if (found == values.end()) return {};
    if (found->second.empty() || found->second.front() == '-') {
        throw std::runtime_error(
            "invalid unsigned 64-bit value for unstructured key '" +
            key + "'");
    }
    try {
        std::size_t consumed = 0;
        const int base =
            found->second.size() > 2 &&
                    found->second[0] == '0' &&
                    (found->second[1] == 'x' ||
                     found->second[1] == 'X')
                ? 16
                : 10;
        const auto value = std::stoull(
            found->second, &consumed, base);
        if (consumed != found->second.size()) {
            throw std::invalid_argument("trailing");
        }
        return static_cast<std::uint64_t>(value);
    } catch (const std::exception&) {
        throw std::runtime_error(
            "invalid unsigned 64-bit value for unstructured key '" +
            key + "': " + found->second);
    }
}

InitializationAcceptanceConfig parse_initialization_acceptance(
    const Values& values) {
    InitializationAcceptanceConfig result;
    const auto optional_tolerance =
        [&](const std::string& key,
            std::optional<double>& destination) {
            if (values.contains(key)) {
                destination = number<double>(values, key, 0.0);
            }
        };
    optional_tolerance(
        "initialization_max_relative_charge_imbalance",
        result.max_relative_charge_imbalance);
    optional_tolerance(
        "initialization_max_relative_current_imbalance",
        result.max_relative_current_imbalance);
    optional_tolerance(
        "initialization_max_relative_pair_imbalance",
        result.max_relative_pair_imbalance);
    if (values.contains("initialization_charge_pairs")) {
        const std::string pair_list = trim(
            values.at("initialization_charge_pairs"));
        if (pair_list.empty() || pair_list.back() == ',') {
            throw std::runtime_error(
                "initialization_charge_pairs cannot be empty or end with a comma");
        }
        std::stringstream pairs(pair_list);
        std::string entry;
        while (std::getline(pairs, entry, ',')) {
            entry = trim(entry);
            const auto separator = entry.find(':');
            if (entry.empty() ||
                separator == std::string::npos ||
                entry.find(':', separator + 1) !=
                    std::string::npos) {
                throw std::runtime_error(
                    "initialization_charge_pairs must be a comma-separated list of first:second species names");
            }
            result.charge_pairs.push_back({
                trim(entry.substr(0, separator)),
                trim(entry.substr(separator + 1))});
        }
    }
    validate_initialization_acceptance(
        result, "unstructured initialization acceptance config");
    return result;
}

std::filesystem::path resolved_path(const std::filesystem::path& config_path,
                                    const std::string& value) {
    std::filesystem::path path(value);
    if (path.is_absolute()) return path;
    const auto relative_to_config = config_path.parent_path() / path;
    if (std::filesystem::exists(relative_to_config)) return relative_to_config;
    return path;
}

} // namespace

bool config_uses_unstructured_mesh_2d(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open config: " + path.string());
    std::string line;
    bool in_global = true;
    while (std::getline(input, line)) {
        const auto comment = line.find_first_of("#;");
        if (comment != std::string::npos) line.resize(comment);
        line = trim(line);
        if (line.empty()) continue;
        if (line.front() == '[' && line.back() == ']') {
            in_global = false;
            continue;
        }
        if (!in_global) continue;
        const auto separator = line.find('=');
        if (separator == std::string::npos) continue;
        if (lower(trim(line.substr(0, separator))) != "mesh") continue;
        const std::string value = lower(trim(line.substr(separator + 1)));
        return value == "imported" || value == "unstructured";
    }
    return false;
}

UnstructuredSimulation2DConfig load_unstructured_config_2d(
    const std::filesystem::path& path) {
    const ParsedConfig parsed = parse(path);
    const auto& global = parsed.global;
    if (number<std::size_t>(global, "config_version", 1) != 1) {
        throw std::runtime_error("unstructured 2D config supports config_version = 1");
    }
    if (number<std::size_t>(global, "dimension", 0) != 2) {
        throw std::runtime_error("unstructured 2D config requires dimension = 2");
    }
    const std::string mesh_kind = lower(required(global, "mesh", "global"));
    if (mesh_kind != "imported" && mesh_kind != "unstructured") {
        throw std::runtime_error("unstructured 2D config requires mesh = imported");
    }

    UnstructuredSimulation2DConfig result;
    if (global.contains("units")) {
        const std::string units = lower(global.at("units"));
        if (units == "normalized" || units == "normalised") {
            result.units.system = UnitSystem::Normalized;
        } else if (units == "si") {
            result.units.system = UnitSystem::SI;
        } else {
            throw std::runtime_error(
                "invalid unstructured units value: " + units);
        }
    }
    result.units.relative_permittivity = number<double>(
        global, "relative_permittivity",
        result.units.relative_permittivity);
    result.poisson.permittivity = result.units.permittivity();
    result.mesh_path =
        resolved_path(path, required(global, "mesh_file", "global"));
    result.dt = number<double>(global, "dt", result.dt);
    result.steps = number<std::size_t>(global, "steps", result.steps);
    result.steady_tolerance =
        number<double>(global, "steady_tolerance", result.steady_tolerance);
    result.steady_window =
        number<std::size_t>(global, "steady_window", result.steady_window);
    result.max_steps =
        number<std::size_t>(global, "max_steps", result.max_steps);
    result.max_particles_per_species = number<std::size_t>(
        global, "max_particles_per_species",
        result.max_particles_per_species);
    result.seed = number<unsigned>(global, "seed", result.seed);
    result.magnetic_field_x =
        number<double>(
            global, "magnetic_field_x",
            result.magnetic_field_x);
    result.magnetic_field_y =
        number<double>(
            global, "magnetic_field_y",
            result.magnetic_field_y);
    result.magnetic_field_z =
        number<double>(global, "magnetic_field_z", result.magnetic_field_z);
    result.output_interval =
        number<std::size_t>(global, "output_interval", result.output_interval);
    if (global.contains("output_dir")) result.output_dir = global.at("output_dir");
    result.vtk_output = boolean(global, "vtk_output", result.vtk_output);
    result.particle_output =
        boolean(global, "particle_output", result.particle_output);
    result.particle_output_interval = number<std::size_t>(
        global, "particle_output_interval", result.particle_output_interval);
    result.particle_output_stride = number<std::size_t>(
        global, "particle_output_stride", result.particle_output_stride);
    result.particle_sample_count = number<std::size_t>(
        global, "particle_sample_count", result.particle_sample_count);
    result.checkpoint_output =
        boolean(global, "checkpoint_output", result.checkpoint_output);
    result.checkpoint_interval = number<std::size_t>(
        global, "checkpoint_interval", result.checkpoint_interval);
    if (global.contains("checkpoint_path")) {
        result.checkpoint_path = global.at("checkpoint_path");
    }
    if (global.contains("restart_path")) {
        result.restart_path = global.at("restart_path");
    }
    if (global.contains("initial_state_path")) {
        result.initial_state_path =
            global.at("initial_state_path");
        if (result.initial_state_path.is_relative()) {
            result.initial_state_path =
                path.parent_path() /
                result.initial_state_path;
        }
    }
    result.initial_state_signature = optional_uint64(
        global, "initial_state_signature");
    if (!result.restart_path.empty() &&
        !result.initial_state_path.empty()) {
        throw std::runtime_error(
            "unstructured restart_path and initial_state_path are mutually exclusive");
    }
    if (result.initial_state_signature &&
        result.initial_state_path.empty()) {
        throw std::runtime_error(
            "unstructured initial_state_signature requires initial_state_path");
    }
    result.poisson.relative_tolerance = number<double>(
        global, "poisson_relative_tolerance", result.poisson.relative_tolerance);
    result.poisson.absolute_tolerance = number<double>(
        global, "poisson_absolute_tolerance", result.poisson.absolute_tolerance);
    result.poisson.max_iterations = number<std::size_t>(
        global, "poisson_max_iterations", result.poisson.max_iterations);

    if (global.contains("mode")) {
        const std::string value = lower(global.at("mode"));
        if (value == "transient") result.mode = RunMode::Transient;
        else if (value == "steady" || value == "steady_state") {
            result.mode = RunMode::SteadyState;
        } else {
            throw std::runtime_error("invalid unstructured run mode: " + value);
        }
    }
    if (global.contains("runtime_backend")) {
        const std::string value = lower(global.at("runtime_backend"));
        if (value == "serial") result.runtime.backend = RuntimeBackend::Serial;
        else if (value == "openmp") result.runtime.backend = RuntimeBackend::OpenMP;
        else throw std::runtime_error("invalid unstructured runtime backend: " + value);
    }
    result.runtime.threads =
        number<std::size_t>(global, "runtime_threads", result.runtime.threads);
    result.initialization_acceptance =
        parse_initialization_acceptance(global);

    if (!parsed.collisions.empty()) {
        result.collisions.enabled =
            boolean(parsed.collisions, "enabled", true);
        const std::string model = lower(
            parsed.collisions.contains("model")
                ? parsed.collisions.at("model")
                : "null_collision");
        if (model != "null_collision" && model != "null-collision" &&
            model != "mcc") {
            throw std::runtime_error(
                "imported 2D collisions support only model = null_collision");
        }
        result.collisions.model = CollisionModelKind::NullCollision;
        if (result.collisions.enabled) {
            result.collisions.gas_data_units =
                result.units.system;
            result.collisions.species = required(
                parsed.collisions, "species", "collisions");
            (void)required(
                parsed.collisions, "neutral_density", "collisions");
            (void)required(
                parsed.collisions, "neutral_temperature", "collisions");
            (void)required(
                parsed.collisions, "max_frequency", "collisions");
            if (parsed.collisions.contains("gas_data_file")) {
                if (parsed.collisions.contains("gas") ||
                    parsed.collisions.contains("neutral_mass") ||
                    parsed.collisions.contains("data_provenance")) {
                    throw std::runtime_error(
                        "gas_data_file cannot be combined with inline gas, "
                        "neutral_mass, or data_provenance");
                }
                result.collisions.gas_data_file = resolved_path(
                    path, parsed.collisions.at("gas_data_file"));
                const auto dataset = load_gas_dataset(
                    result.collisions.gas_data_file);
                if (dataset.unit_system != result.units.system) {
                    throw std::runtime_error(
                        "gas dataset units '" +
                        to_string(dataset.unit_system) +
                        "' do not match simulation units '" +
                        to_string(result.units.system) + "'");
                }
                result.collisions.gas_data_version =
                    dataset.format_version;
                result.collisions.gas_data_units =
                    dataset.unit_system;
                result.collisions.gas_name = dataset.gas_name;
                result.collisions.neutral_mass =
                    dataset.neutral_mass;
                result.collisions.data_provenance =
                    dataset.data_provenance;
                result.collisions.dataset_id = dataset.dataset_id;
                result.collisions.dataset_version =
                    dataset.dataset_version;
                result.collisions.citation = dataset.citation;
                result.collisions.retrieved = dataset.retrieved;
                result.collisions.license = dataset.license;
                result.collisions.channels = dataset.channels;
            } else {
                result.collisions.gas_name = required(
                    parsed.collisions, "gas", "collisions");
                result.collisions.data_provenance = required(
                    parsed.collisions, "data_provenance", "collisions");
                (void)required(
                    parsed.collisions, "neutral_mass", "collisions");
            }
        }
        result.collisions.neutral_density = number<double>(
            parsed.collisions, "neutral_density", 0.0);
        if (result.collisions.gas_data_file.empty()) {
            result.collisions.neutral_mass = number<double>(
                parsed.collisions, "neutral_mass", 0.0);
        }
        result.collisions.neutral_temperature = number<double>(
            parsed.collisions, "neutral_temperature", 0.0);
        result.collisions.max_frequency = number<double>(
            parsed.collisions, "max_frequency", 0.0);
        result.collisions.max_candidates_per_particle =
            number<std::size_t>(
                parsed.collisions, "max_candidates_per_particle",
                result.collisions.max_candidates_per_particle);
    }

    for (const auto& boundary : parsed.boundaries) {
        const std::string field = boundary.values.contains("field")
                                      ? lower(boundary.values.at("field"))
                                      : "dirichlet";
        if (field == "dirichlet") {
            (void)required(boundary.values, "potential",
                           "Dirichlet boundary '" + boundary.name + "'");
            if (boundary.values.contains("normal_derivative")) {
                throw std::runtime_error(
                    "Dirichlet boundary '" + boundary.name +
                    "' must not define normal_derivative");
            }
            result.dirichlet_potentials.emplace(
                boundary.name,
                number<double>(boundary.values, "potential",
                               std::numeric_limits<double>::quiet_NaN()));
        } else if (field == "neumann") {
            (void)required(boundary.values, "normal_derivative",
                           "Neumann boundary '" + boundary.name + "'");
            if (boundary.values.contains("potential")) {
                throw std::runtime_error(
                    "Neumann boundary '" + boundary.name +
                    "' must not define potential");
            }
            result.neumann_normal_derivatives.emplace(
                boundary.name,
                number<double>(boundary.values, "normal_derivative",
                               std::numeric_limits<double>::quiet_NaN()));
        } else {
            throw std::runtime_error(
                "boundary '" + boundary.name +
                "' field type must be dirichlet or neumann");
        }
        const std::string particle =
            lower(required(boundary.values, "particle",
                           "boundary '" + boundary.name + "'"));
        if (particle == "absorbing") {
            result.particle_boundaries.emplace(
                boundary.name, ParticleBoundary::Absorbing);
        } else if (particle == "reflecting") {
            result.particle_boundaries.emplace(
                boundary.name, ParticleBoundary::Reflecting);
        } else {
            throw std::runtime_error(
                "boundary '" + boundary.name +
                "' particle policy must be absorbing or reflecting");
        }
    }
    if (parsed.boundaries.empty()) {
        throw std::runtime_error("unstructured 2D config requires boundary sections");
    }
    if (result.dirichlet_potentials.empty()) {
        throw std::runtime_error(
            "unstructured 2D config requires at least one Dirichlet field boundary");
    }

    for (const auto& species : parsed.species) {
        UnstructuredSpecies2DConfig value;
        value.name = species.name;
        value.charge = number<double>(species.values, "charge", value.charge);
        value.mass = number<double>(species.values, "mass", value.mass);
        value.weight = number<double>(species.values, "weight", value.weight);
        value.particles =
            number<std::size_t>(species.values, "particles", value.particles);
        value.drift_velocity_x = number<double>(
            species.values, "drift_velocity_x", value.drift_velocity_x);
        value.drift_velocity_y = number<double>(
            species.values, "drift_velocity_y", value.drift_velocity_y);
        value.drift_velocity_z = number<double>(
            species.values, "drift_velocity_z",
            value.drift_velocity_z);
        value.thermal_velocity = number<double>(
            species.values, "thermal_velocity", value.thermal_velocity);
        value.initialization.version = number<std::size_t>(
            species.values, "initialization_version",
            value.initialization.version);
        if (species.values.contains("loading")) {
            value.initialization.loading =
                particle_loading_from_string(
                    required(
                        species.values, "loading",
                        "species '" + species.name + "'"));
        }
        if (species.values.contains("thermal_velocity_x")) {
            value.initialization.thermal_velocity_x = number<double>(
                species.values, "thermal_velocity_x", 0.0);
        }
        if (species.values.contains("thermal_velocity_y")) {
            value.initialization.thermal_velocity_y = number<double>(
                species.values, "thermal_velocity_y", 0.0);
        }
        if (species.values.contains("thermal_velocity_z")) {
            value.initialization.thermal_velocity_z = number<double>(
                species.values, "thermal_velocity_z", 0.0);
        }
        if (species.values.contains("initialization_region")) {
            value.initialization_region = required(
                species.values, "initialization_region",
                "species '" + species.name + "'");
        }
        if (species.values.contains("density_profile")) {
            value.initialization.density_profile =
                density_profile_from_string(lower(required(
                    species.values, "density_profile",
                    "species '" + species.name + "'")));
        }
        if (species.values.contains("profile_center_x")) {
            value.initialization.profile_center_x = number<double>(
                species.values, "profile_center_x", 0.0);
        }
        if (species.values.contains("profile_center_y")) {
            value.initialization.profile_center_y = number<double>(
                species.values, "profile_center_y", 0.0);
        }
        if (species.values.contains("profile_scale_x")) {
            value.initialization.profile_scale_x = number<double>(
                species.values, "profile_scale_x", 0.0);
        }
        if (species.values.contains("profile_scale_y")) {
            value.initialization.profile_scale_y = number<double>(
                species.values, "profile_scale_y", 0.0);
        }
        if (species.values.contains("profile_amplitude")) {
            value.initialization.profile_amplitude = number<double>(
                species.values, "profile_amplitude", 0.0);
        }
        if (species.values.contains("profile_phase")) {
            value.initialization.profile_phase = number<double>(
                species.values, "profile_phase", 0.0);
        }
        if (species.values.contains("profile_mode_x")) {
            value.initialization.profile_mode_x =
                number<std::size_t>(
                    species.values, "profile_mode_x", 0);
        }
        if (species.values.contains("profile_mode_y")) {
            value.initialization.profile_mode_y =
                number<std::size_t>(
                    species.values, "profile_mode_y", 0);
        }
        value.initialization.max_profile_sampling_attempts =
            number<std::size_t>(
                species.values, "max_profile_sampling_attempts",
                value.initialization.max_profile_sampling_attempts);
        const bool any_bounds =
            species.values.contains("init_x_min") ||
            species.values.contains("init_x_max") ||
            species.values.contains("init_y_min") ||
            species.values.contains("init_y_max");
        if (any_bounds) {
            value.initialization_minimum = Vec2{
                number<double>(species.values, "init_x_min",
                               std::numeric_limits<double>::quiet_NaN()),
                number<double>(species.values, "init_y_min",
                               std::numeric_limits<double>::quiet_NaN()),
            };
            value.initialization_maximum = Vec2{
                number<double>(species.values, "init_x_max",
                               std::numeric_limits<double>::quiet_NaN()),
                number<double>(species.values, "init_y_max",
                               std::numeric_limits<double>::quiet_NaN()),
            };
        }
        result.species.push_back(value);
    }
    if (parsed.species.empty()) {
        throw std::runtime_error("unstructured 2D config requires species sections");
    }

    std::set<std::string> configured_species;
    for (const auto& species : result.species) {
        configured_species.insert(species.name);
    }
    for (const auto& source : parsed.sources) {
        UnstructuredBoundarySource2DConfig value;
        value.name = source.name;
        value.species = required(
            source.values, "species", "source '" + source.name + "'");
        value.boundary = required(
            source.values, "boundary", "source '" + source.name + "'");
        (void)required(
            source.values, "particles_per_step",
            "source '" + source.name + "'");
        value.particles_per_step = number<std::size_t>(
            source.values, "particles_per_step", 0);
        value.start_step = number<std::size_t>(
            source.values, "start_step", value.start_step);
        value.end_step = number<std::size_t>(
            source.values, "end_step", value.end_step);
        value.normal_velocity = number<double>(
            source.values, "normal_velocity", value.normal_velocity);
        value.tangential_velocity = number<double>(
            source.values, "tangential_velocity", value.tangential_velocity);
        value.thermal_velocity = number<double>(
            source.values, "thermal_velocity", value.thermal_velocity);
        value.out_of_plane_velocity = number<double>(
            source.values, "out_of_plane_velocity",
            value.out_of_plane_velocity);
        if (!configured_species.contains(value.species)) {
            throw std::runtime_error(
                "source '" + source.name +
                "' references unknown species '" + value.species + "'");
        }
        result.sources.push_back(std::move(value));
    }
    for (const auto& emission : parsed.emissions) {
        UnstructuredSecondaryEmission2DConfig value;
        value.name = emission.name;
        value.boundary = required(
            emission.values, "boundary", "emission '" + emission.name + "'");
        value.incident_species = required(
            emission.values, "incident_species",
            "emission '" + emission.name + "'");
        value.emitted_species = required(
            emission.values, "emitted_species",
            "emission '" + emission.name + "'");
        (void)required(
            emission.values, "yield", "emission '" + emission.name + "'");
        value.yield = number<double>(emission.values, "yield", 0.0);
        value.max_particles_per_impact = number<std::size_t>(
            emission.values, "max_particles_per_impact",
            value.max_particles_per_impact);
        value.normal_velocity = number<double>(
            emission.values, "normal_velocity", value.normal_velocity);
        value.tangential_velocity = number<double>(
            emission.values, "tangential_velocity",
            value.tangential_velocity);
        value.thermal_velocity = number<double>(
            emission.values, "thermal_velocity", value.thermal_velocity);
        value.out_of_plane_velocity = number<double>(
            emission.values, "out_of_plane_velocity",
            value.out_of_plane_velocity);
        if (!configured_species.contains(value.incident_species) ||
            !configured_species.contains(value.emitted_species)) {
            throw std::runtime_error(
                "emission '" + emission.name +
                "' references an unknown species");
        }
        result.emissions.push_back(std::move(value));
    }
    for (const auto& channel : parsed.collision_channels) {
        if (!result.collisions.gas_data_file.empty()) {
            const auto configured = std::find_if(
                result.collisions.channels.begin(),
                result.collisions.channels.end(),
                [&](const auto& value) {
                    return value.name == channel.name;
                });
            if (configured == result.collisions.channels.end()) {
                throw std::runtime_error(
                    "gas dataset has no collision channel named '" +
                    channel.name + "'");
            }
            if (configured->process !=
                    CollisionProcessKind::Ionization &&
                configured->process !=
                    CollisionProcessKind::Attachment) {
                throw std::runtime_error(
                    "only reactive product species may be mapped for "
                    "gas dataset channel '" + channel.name + "'");
            }
            for (const auto& [key, unused] : channel.values) {
                (void)unused;
                const bool ionization_key =
                    configured->process ==
                        CollisionProcessKind::Ionization &&
                    (key == "secondary_species" ||
                     key == "ion_species");
                const bool attachment_key =
                    configured->process ==
                        CollisionProcessKind::Attachment &&
                    key == "attachment_species";
                if (!ionization_key && !attachment_key) {
                    throw std::runtime_error(
                        "gas dataset channel '" + channel.name +
                        "' physics cannot be overridden by a simulation");
                }
            }
            if (configured->process ==
                CollisionProcessKind::Ionization) {
                configured->secondary_species = required(
                    channel.values, "secondary_species",
                    "ionization channel '" + channel.name + "'");
                configured->ion_species = required(
                    channel.values, "ion_species",
                    "ionization channel '" + channel.name + "'");
            } else {
                configured->attachment_species = required(
                    channel.values, "attachment_species",
                    "attachment channel '" + channel.name + "'");
            }
            continue;
        }
        CollisionChannelConfig value;
        value.name = channel.name;
        const std::string type = lower(required(
            channel.values, "type",
            "collision channel '" + channel.name + "'"));
        if (type == "elastic") {
            value.process = CollisionProcessKind::Elastic;
        } else if (type == "excitation") {
            value.process = CollisionProcessKind::Excitation;
        } else if (type == "ionization") {
            value.process = CollisionProcessKind::Ionization;
        } else if (type == "attachment") {
            value.process = CollisionProcessKind::Attachment;
        } else if (type == "charge_exchange" ||
                   type == "charge-exchange") {
            value.process =
                CollisionProcessKind::ChargeExchange;
        } else {
            throw std::runtime_error(
                "collision channel '" + channel.name +
                "' type must be elastic, excitation, ionization, "
                "attachment, or charge_exchange");
        }
        value.cross_section_file = resolved_path(
            path, required(
                channel.values, "cross_section_file",
                "collision channel '" + channel.name + "'"));
        value.threshold_energy = number<double>(
            channel.values, "threshold_energy",
            value.threshold_energy);
        value.energy_scale = number<double>(
            channel.values, "energy_scale", value.energy_scale);
        value.cross_section_scale = number<double>(
            channel.values, "cross_section_scale",
            value.cross_section_scale);
        const std::string angular_model = lower(
            channel.values.contains("angular_model")
                ? channel.values.at("angular_model")
                : "isotropic");
        if (angular_model == "isotropic") {
            value.angular_scattering =
                AngularScatteringKind::Isotropic;
        } else if (
            angular_model == "henyey_greenstein" ||
            angular_model == "henyey-greenstein") {
            value.angular_scattering =
                AngularScatteringKind::HenyeyGreenstein;
        } else {
            throw std::runtime_error(
                "collision channel '" + channel.name +
                "' angular_model must be isotropic or "
                "henyey_greenstein");
        }
        if (channel.values.contains("mean_cosine_file")) {
            value.mean_cosine_file = resolved_path(
                path, channel.values.at("mean_cosine_file"));
        }
        value.mean_cosine_energy_scale = number<double>(
            channel.values, "mean_cosine_energy_scale",
            value.mean_cosine_energy_scale);
        if (value.angular_scattering !=
                AngularScatteringKind::Isotropic &&
            value.process != CollisionProcessKind::Elastic) {
            throw std::runtime_error(
                "collision channel '" + channel.name +
                "' anisotropic scattering is valid only for elastic "
                "channels");
        }
        if (value.angular_scattering ==
            AngularScatteringKind::HenyeyGreenstein) {
            if (value.mean_cosine_file.empty()) {
                throw std::runtime_error(
                    "collision channel '" + channel.name +
                    "' Henyey-Greenstein scattering requires "
                    "mean_cosine_file");
            }
        } else if (!value.mean_cosine_file.empty() ||
                   value.mean_cosine_energy_scale != 1.0) {
            throw std::runtime_error(
                "collision channel '" + channel.name +
                "' mean-cosine data requires angular_model = "
                "henyey_greenstein");
        }
        if (value.process == CollisionProcessKind::Ionization) {
            if (channel.values.contains("attachment_species")) {
                throw std::runtime_error(
                    "ionization channel '" + channel.name +
                    "' does not accept attachment_species");
            }
            value.secondary_species = required(
                channel.values, "secondary_species",
                "ionization channel '" + channel.name + "'");
            value.ion_species = required(
                channel.values, "ion_species",
                "ionization channel '" + channel.name + "'");
        } else if (value.process ==
                   CollisionProcessKind::Attachment) {
            if (channel.values.contains("secondary_species") ||
                channel.values.contains("ion_species")) {
                throw std::runtime_error(
                    "attachment channel '" + channel.name +
                    "' accepts only attachment_species");
            }
            value.attachment_species = required(
                channel.values, "attachment_species",
                "attachment channel '" + channel.name + "'");
        } else if (channel.values.contains("secondary_species") ||
                   channel.values.contains("ion_species") ||
                   channel.values.contains("attachment_species")) {
            throw std::runtime_error(
                "collision channel '" + channel.name +
                "' product species do not match its process");
        }
        result.collisions.channels.push_back(std::move(value));
    }
    if (result.collisions.enabled) {
        if (parsed.collisions.empty()) {
            throw std::runtime_error(
                "enabled imported collisions require a [collisions] section");
        }
        if (!configured_species.contains(result.collisions.species)) {
            throw std::runtime_error(
                "collisions reference unknown species '" +
                result.collisions.species + "'");
        }
        if (result.collisions.gas_name.empty() ||
            result.collisions.data_provenance.empty()) {
            throw std::runtime_error(
                "imported MCC requires non-empty gas and data_provenance");
        }
        if (result.collisions.channels.empty()) {
            throw std::runtime_error(
                "imported MCC requires collision channel sections");
        }
        for (const auto& channel : result.collisions.channels) {
            if (channel.process == CollisionProcessKind::Ionization &&
                (!configured_species.contains(
                     channel.secondary_species) ||
                 !configured_species.contains(channel.ion_species))) {
                throw std::runtime_error(
                    "ionization channel '" + channel.name +
                    "' references unknown product species");
            }
            if (channel.process == CollisionProcessKind::Attachment &&
                !configured_species.contains(
                    channel.attachment_species)) {
                throw std::runtime_error(
                    "attachment channel '" + channel.name +
                    "' references an unknown product species");
            }
        }
        if (!result.collisions.gas_data_file.empty()) {
            for (const auto& channel : result.collisions.channels) {
                if (channel.process == CollisionProcessKind::Ionization &&
                    (channel.secondary_species.empty() ||
                     channel.ion_species.empty())) {
                    throw std::runtime_error(
                        "gas dataset ionization channel '" + channel.name +
                        "' requires a matching [collision." +
                        channel.name + "] product mapping");
                }
                if (channel.process == CollisionProcessKind::Attachment &&
                    channel.attachment_species.empty()) {
                    throw std::runtime_error(
                        "gas dataset attachment channel '" + channel.name +
                        "' requires a matching [collision." +
                        channel.name + "] product mapping");
                }
            }
        }
    } else if (!parsed.collision_channels.empty()) {
        throw std::runtime_error(
            "collision channels require enabled imported collisions");
    }
    return result;
}

} // namespace pic
