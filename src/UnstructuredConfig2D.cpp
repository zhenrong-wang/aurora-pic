#include "pic/UnstructuredSimulation2D.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <set>
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
    std::vector<NamedBlock> boundaries;
    std::vector<NamedBlock> species;
    std::vector<NamedBlock> sources;
    std::vector<NamedBlock> emissions;
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
        "mode", "steady_tolerance", "steady_window", "max_steps",
        "max_particles_per_species", "seed",
        "magnetic_field_z", "output_interval", "output_dir", "vtk_output",
        "particle_output", "particle_output_interval", "particle_output_stride",
        "particle_sample_count", "checkpoint_output", "checkpoint_interval",
        "checkpoint_path", "restart_path", "runtime_backend", "runtime_threads",
        "poisson_relative_tolerance", "poisson_absolute_tolerance",
        "poisson_max_iterations",
    };
    static const std::set<std::string> boundary_keys{
        "field", "potential", "normal_derivative", "particle"};
    static const std::set<std::string> species_keys{
        "charge", "mass", "weight", "particles", "drift_velocity_x",
        "drift_velocity_y", "thermal_velocity", "init_x_min", "init_x_max",
        "init_y_min", "init_y_max",
    };
    static const std::set<std::string> source_keys{
        "species", "boundary", "particles_per_step", "start_step", "end_step",
        "normal_velocity", "tangential_velocity", "thermal_velocity",
    };
    static const std::set<std::string> emission_keys{
        "boundary", "incident_species", "emitted_species", "yield",
        "max_particles_per_impact", "normal_velocity", "tangential_velocity",
        "thermal_velocity",
    };

    ParsedConfig result;
    Values* current = &result.global;
    const std::set<std::string>* allowed = &global_keys;
    std::set<std::string> boundary_names;
    std::set<std::string> species_names;
    std::set<std::string> source_names;
    std::set<std::string> emission_names;
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
            if (lowered.rfind(boundary_prefix, 0) == 0) {
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
        value.thermal_velocity = number<double>(
            species.values, "thermal_velocity", value.thermal_velocity);
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
        if (!configured_species.contains(value.incident_species) ||
            !configured_species.contains(value.emitted_species)) {
            throw std::runtime_error(
                "emission '" + emission.name +
                "' references an unknown species");
        }
        result.emissions.push_back(std::move(value));
    }
    return result;
}

} // namespace pic
