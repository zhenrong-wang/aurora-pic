#include "pic/Config.hpp"
#include "pic/Simulation2D.hpp"
#include "pic/Simulation3D.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace pic {
namespace {
using KeyValue = std::unordered_map<std::string, std::string>;

std::string trim(std::string s) {
    auto not_space = [](unsigned char c){ return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
    return s;
}

std::string lower(std::string s) {
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

std::runtime_error config_error(std::size_t line_number, const std::string& message) {
    std::ostringstream out;
    out << "config line " << line_number << ": " << message;
    return std::runtime_error(out.str());
}

bool starts_with(const std::string& s, const std::string& prefix) {
    return s.rfind(prefix, 0) == 0;
}

void ensure_key_allowed(const std::unordered_set<std::string>& allowed,
                        const std::string& key,
                        const std::string& section,
                        std::size_t line_number) {
    if (!allowed.count(key)) throw config_error(line_number, "unknown key '" + key + "' in [" + section + "]");
}

void assign_key(KeyValue& kv, const std::string& key, const std::string& value, std::size_t line_number) {
    if (key.empty()) throw config_error(line_number, "empty key");
    if (value.empty()) throw config_error(line_number, "empty value for key '" + key + "'");
    if (!kv.emplace(key, value).second) throw config_error(line_number, "duplicate key '" + key + "'");
}

void require_full_parse(const std::string& value, std::size_t parsed, const std::string& key) {
    if (parsed != value.size()) throw std::runtime_error("invalid value for '" + key + "': trailing characters in '" + value + "'");
}

template<class T> T as(const KeyValue& kv, const std::string& key, T def);

template<> double as<double>(const KeyValue& kv, const std::string& key, double def) {
    auto it = kv.find(key);
    if (it == kv.end()) return def;
    try {
        std::size_t parsed = 0;
        const double value = std::stod(it->second, &parsed);
        require_full_parse(it->second, parsed, key);
        if (!std::isfinite(value)) throw std::runtime_error("not finite");
        return value;
    } catch (const std::exception&) {
        throw std::runtime_error("invalid floating-point value for '" + key + "': '" + it->second + "'");
    }
}

template<> std::size_t as<std::size_t>(const KeyValue& kv, const std::string& key, std::size_t def) {
    auto it = kv.find(key);
    if (it == kv.end()) return def;
    try {
        std::size_t parsed = 0;
        const long long signed_value = std::stoll(it->second, &parsed);
        require_full_parse(it->second, parsed, key);
        if (signed_value < 0) throw std::runtime_error("negative");
        return static_cast<std::size_t>(signed_value);
    } catch (const std::exception&) {
        throw std::runtime_error("invalid non-negative integer value for '" + key + "': '" + it->second + "'");
    }
}

template<> unsigned as<unsigned>(const KeyValue& kv, const std::string& key, unsigned def) {
    const std::size_t value = as<std::size_t>(kv, key, static_cast<std::size_t>(def));
    if (value > static_cast<std::size_t>(std::numeric_limits<unsigned>::max())) {
        throw std::runtime_error("integer value for '" + key + "' exceeds unsigned range");
    }
    return static_cast<unsigned>(value);
}

template<> std::string as<std::string>(const KeyValue& kv, const std::string& key, std::string def) {
    auto it = kv.find(key);
    return it == kv.end() ? def : it->second;
}

bool parse_bool(const KeyValue& kv, const std::string& key, bool def) {
    auto it = kv.find(key);
    if (it == kv.end()) return def;
    const auto x = lower(trim(it->second));
    if (x == "1" || x == "true" || x == "yes" || x == "on") return true;
    if (x == "0" || x == "false" || x == "no" || x == "off") return false;
    throw std::runtime_error("invalid boolean value for '" + key + "': '" + it->second + "'");
}

Boundary parse_boundary(const KeyValue& kv, Boundary def) {
    const auto value = lower(as<std::string>(kv, "boundary", to_string(def)));
    if (value == "periodic") return Boundary::Periodic;
    if (value == "dirichlet" || value == "absorbing") return Boundary::Dirichlet;
    throw std::runtime_error("invalid boundary value: '" + value + "'");
}

RunMode parse_mode(const KeyValue& kv, RunMode def) {
    const auto value = lower(as<std::string>(kv, "mode", to_string(def)));
    if (value == "transient") return RunMode::Transient;
    if (value == "steady" || value == "steady_state") return RunMode::SteadyState;
    throw std::runtime_error("invalid mode value: '" + value + "'");
}

ParticleBoundary parse_particle_boundary(const KeyValue& kv, const std::string& key, ParticleBoundary def) {
    const auto value = lower(as<std::string>(kv, key, to_string(def)));
    if (value == "auto") return ParticleBoundary::Auto;
    if (value == "absorbing" || value == "absorb") return ParticleBoundary::Absorbing;
    if (value == "reflecting" || value == "reflective" || value == "reflect") return ParticleBoundary::Reflecting;
    if (value == "periodic" || value == "wrap") return ParticleBoundary::Periodic;
    throw std::runtime_error("invalid particle boundary value for '" + key + "': '" + value + "'");
}

void validate_positive(double value, const std::string& name) {
    if (!(value > 0.0)) throw std::runtime_error(name + " must be positive");
}
void validate_non_negative(double value, const std::string& name) {
    if (value < 0.0) throw std::runtime_error(name + " must be non-negative");
}

void validate_config(const Config& cfg) {
    if (cfg.nx < 3) throw std::runtime_error("nx must be at least 3");
    validate_positive(cfg.length, "length");
    validate_positive(cfg.dt, "dt");
    if (cfg.output_interval == 0) throw std::runtime_error("output_interval must be positive");
    validate_positive(cfg.steady_tolerance, "steady_tolerance");
    if (cfg.steady_window == 0) throw std::runtime_error("steady_window must be positive");
    if (cfg.mode == RunMode::SteadyState && cfg.max_steps == 0) throw std::runtime_error("max_steps must be positive for steady-state mode");
    validate_non_negative(cfg.collisions.frequency, "collision frequency");
    validate_non_negative(cfg.collisions.neutral_temperature_velocity, "neutral_temperature_velocity");
    for (const auto& s : cfg.species) {
        if (s.name.empty()) throw std::runtime_error("species name must not be empty");
        validate_positive(s.mass, "species '" + s.name + "' mass");
        validate_positive(s.weight, "species '" + s.name + "' weight");
        validate_positive(s.density, "species '" + s.name + "' density");
        if (s.particles == 0) throw std::runtime_error("species '" + s.name + "' particles must be positive");
        validate_non_negative(s.thermal_velocity, "species '" + s.name + "' thermal_velocity");
        if (s.init_x_min < 0.0) throw std::runtime_error("species '" + s.name + "' init_x_min must be non-negative");
        const double xmax = s.init_x_max < 0.0 ? cfg.length : s.init_x_max;
        if (xmax > cfg.length) throw std::runtime_error("species '" + s.name + "' init_x_max exceeds domain length");
        if (!(s.init_x_min < xmax)) throw std::runtime_error("species '" + s.name + "' initialization interval must have positive width");
    }
}

void validate_boundary_side(const BoundarySide2D& side, const std::string& name) {
    if (trim(side.tag).empty()) throw std::runtime_error("2D boundary side '" + name + "' tag must not be empty");
    if (!std::isfinite(side.potential)) throw std::runtime_error("2D boundary side '" + name + "' potential must be finite");
}

void validate_config_2d(const Simulation2DConfig& cfg) {
    if (cfg.nx < 3) throw std::runtime_error("2D nx must be at least 3");
    if (cfg.ny < 3) throw std::runtime_error("2D ny must be at least 3");
    validate_positive(cfg.length_x, "length_x");
    validate_positive(cfg.length_y, "length_y");
    validate_positive(cfg.dt, "dt");
    if (cfg.output_interval == 0) throw std::runtime_error("output_interval must be positive");
    if (cfg.particle_output_stride == 0) throw std::runtime_error("particle_output_stride must be positive");
    if (!std::isfinite(cfg.magnetic_field_z)) throw std::runtime_error("magnetic_field_z must be finite");
    validate_boundary_side(cfg.boundary_config.left, "left");
    validate_boundary_side(cfg.boundary_config.right, "right");
    validate_boundary_side(cfg.boundary_config.bottom, "bottom");
    validate_boundary_side(cfg.boundary_config.top, "top");
    for (const auto& s : cfg.species) {
        if (s.name.empty()) throw std::runtime_error("2D species name must not be empty");
        validate_positive(s.mass, "2D species '" + s.name + "' mass");
        validate_positive(s.weight, "2D species '" + s.name + "' weight");
        if (s.particles == 0) throw std::runtime_error("2D species '" + s.name + "' particles must be positive");
        validate_non_negative(s.thermal_velocity, "2D species '" + s.name + "' thermal_velocity");
        if (s.init_x_min < 0.0) throw std::runtime_error("2D species '" + s.name + "' init_x_min must be non-negative");
        if (s.init_y_min < 0.0) throw std::runtime_error("2D species '" + s.name + "' init_y_min must be non-negative");
        const double xmax = s.init_x_max < 0.0 ? cfg.length_x : s.init_x_max;
        const double ymax = s.init_y_max < 0.0 ? cfg.length_y : s.init_y_max;
        if (xmax > cfg.length_x) throw std::runtime_error("2D species '" + s.name + "' init_x_max exceeds domain length_x");
        if (ymax > cfg.length_y) throw std::runtime_error("2D species '" + s.name + "' init_y_max exceeds domain length_y");
        if (!(s.init_x_min < xmax)) throw std::runtime_error("2D species '" + s.name + "' x initialization interval must have positive width");
        if (!(s.init_y_min < ymax)) throw std::runtime_error("2D species '" + s.name + "' y initialization interval must have positive width");
    }
}

void validate_config_3d(const Simulation3DConfig& cfg) {
    if (cfg.nx < 3) throw std::runtime_error("3D nx must be at least 3");
    if (cfg.ny < 3) throw std::runtime_error("3D ny must be at least 3");
    if (cfg.nz < 3) throw std::runtime_error("3D nz must be at least 3");
    validate_positive(cfg.length_x, "length_x");
    validate_positive(cfg.length_y, "length_y");
    validate_positive(cfg.length_z, "length_z");
    validate_positive(cfg.dt, "dt");
    if (cfg.output_interval == 0) throw std::runtime_error("output_interval must be positive");
    if (cfg.particle_output_stride == 0) throw std::runtime_error("particle_output_stride must be positive");
    if (!std::isfinite(cfg.magnetic_field.x) || !std::isfinite(cfg.magnetic_field.y) || !std::isfinite(cfg.magnetic_field.z)) {
        throw std::runtime_error("magnetic_field components must be finite");
    }
    for (const auto& s : cfg.species) {
        if (s.name.empty()) throw std::runtime_error("3D species name must not be empty");
        validate_positive(s.mass, "3D species '" + s.name + "' mass");
        validate_positive(s.weight, "3D species '" + s.name + "' weight");
        if (s.particles == 0) throw std::runtime_error("3D species '" + s.name + "' particles must be positive");
        validate_non_negative(s.thermal_velocity, "3D species '" + s.name + "' thermal_velocity");
        if (s.init_x_min < 0.0) throw std::runtime_error("3D species '" + s.name + "' init_x_min must be non-negative");
        if (s.init_y_min < 0.0) throw std::runtime_error("3D species '" + s.name + "' init_y_min must be non-negative");
        if (s.init_z_min < 0.0) throw std::runtime_error("3D species '" + s.name + "' init_z_min must be non-negative");
        const double xmax = s.init_x_max < 0.0 ? cfg.length_x : s.init_x_max;
        const double ymax = s.init_y_max < 0.0 ? cfg.length_y : s.init_y_max;
        const double zmax = s.init_z_max < 0.0 ? cfg.length_z : s.init_z_max;
        if (xmax > cfg.length_x) throw std::runtime_error("3D species '" + s.name + "' init_x_max exceeds domain length_x");
        if (ymax > cfg.length_y) throw std::runtime_error("3D species '" + s.name + "' init_y_max exceeds domain length_y");
        if (zmax > cfg.length_z) throw std::runtime_error("3D species '" + s.name + "' init_z_max exceeds domain length_z");
        if (!(s.init_x_min < xmax)) throw std::runtime_error("3D species '" + s.name + "' x initialization interval must have positive width");
        if (!(s.init_y_min < ymax)) throw std::runtime_error("3D species '" + s.name + "' y initialization interval must have positive width");
        if (!(s.init_z_min < zmax)) throw std::runtime_error("3D species '" + s.name + "' z initialization interval must have positive width");
    }
}

struct ParsedBlocks {
    KeyValue global;
    KeyValue collisions;
    std::vector<KeyValue> species_blocks;
};

ParsedBlocks parse_config_blocks(const std::string& path,
                                 const std::unordered_set<std::string>& global_keys,
                                 const std::unordered_set<std::string>& species_keys,
                                 const std::unordered_set<std::string>* collision_keys,
                                 const std::string& loader_name) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open config: " + path);

    ParsedBlocks blocks;
    std::string section = "global", line;
    KeyValue* current = &blocks.global;
    const std::unordered_set<std::string>* allowed = &global_keys;
    std::size_t line_number = 0;
    while (std::getline(in, line)) {
        ++line_number;
        auto comment = line.find_first_of("#;");
        if (comment != std::string::npos) line = line.substr(0, comment);
        line = trim(line);
        if (line.empty()) continue;
        if (line.front() == '[' && line.back() == ']') {
            section = lower(trim(line.substr(1, line.size() - 2)));
            if (section.empty()) throw config_error(line_number, "empty section name");
            if (collision_keys && section == "collisions") {
                current = &blocks.collisions;
                allowed = collision_keys;
            } else if (section == "species" || starts_with(section, "species.")) {
                blocks.species_blocks.emplace_back();
                current = &blocks.species_blocks.back();
                allowed = &species_keys;
                if (starts_with(section, "species.")) {
                    const std::string species_name = trim(section.substr(std::string("species.").size()));
                    if (species_name.empty()) throw config_error(line_number, "empty species section suffix");
                    assign_key(*current, "name", species_name, line_number);
                }
            } else {
                throw config_error(line_number, "unknown section '" + section + "' in " + loader_name + " config");
            }
            continue;
        }
        auto eq = line.find('=');
        if (eq == std::string::npos) throw config_error(line_number, "invalid config line: " + line);
        const std::string key = lower(trim(line.substr(0, eq)));
        ensure_key_allowed(*allowed, key, section, line_number);
        assign_key(*current, key, trim(line.substr(eq + 1)), line_number);
    }
    return blocks;
}

} // namespace

Config load_config(const std::string& path) {
    static const std::unordered_set<std::string> global_keys{
        "nx", "length", "dt", "steps", "output_interval", "output_dir", "seed",
        "phi_left", "phi_right", "steady_tolerance", "steady_window", "max_steps",
        "boundary", "mode", "dimension", "checkpoint_output", "checkpoint_interval",
        "checkpoint_path", "restart_path"
    };
    static const std::unordered_set<std::string> collision_keys{
        "enabled", "frequency", "neutral_temperature_velocity"
    };
    static const std::unordered_set<std::string> species_keys{
        "name", "charge", "mass", "weight", "particles", "density", "drift_velocity",
        "thermal_velocity", "init_x_min", "init_x_max"
    };

    auto blocks = parse_config_blocks(path, global_keys, species_keys, &collision_keys, "1D");
    const auto& global = blocks.global;
    const auto& collision = blocks.collisions;

    const auto dimension = as<std::size_t>(global, "dimension", 1);
    if (dimension != 1) throw std::runtime_error("1D config loader requires dimension = 1 or no dimension key");

    Config cfg;
    cfg.nx = as<std::size_t>(global, "nx", cfg.nx);
    cfg.length = as<double>(global, "length", cfg.length);
    cfg.dt = as<double>(global, "dt", cfg.dt);
    cfg.steps = as<std::size_t>(global, "steps", cfg.steps);
    cfg.output_interval = as<std::size_t>(global, "output_interval", cfg.output_interval);
    cfg.output_dir = as<std::string>(global, "output_dir", cfg.output_dir);
    cfg.seed = as<unsigned>(global, "seed", cfg.seed);
    cfg.phi_left = as<double>(global, "phi_left", cfg.phi_left);
    cfg.phi_right = as<double>(global, "phi_right", cfg.phi_right);
    cfg.steady_tolerance = as<double>(global, "steady_tolerance", cfg.steady_tolerance);
    cfg.steady_window = as<std::size_t>(global, "steady_window", cfg.steady_window);
    cfg.max_steps = as<std::size_t>(global, "max_steps", cfg.max_steps);
    cfg.boundary = parse_boundary(global, cfg.boundary);
    cfg.mode = parse_mode(global, cfg.mode);
    cfg.collisions.enabled = parse_bool(collision, "enabled", cfg.collisions.enabled);
    cfg.collisions.frequency = as<double>(collision, "frequency", cfg.collisions.frequency);
    cfg.collisions.neutral_temperature_velocity = as<double>(collision, "neutral_temperature_velocity", cfg.collisions.neutral_temperature_velocity);
    cfg.checkpoint_output = parse_bool(global, "checkpoint_output", cfg.checkpoint_output);
    cfg.checkpoint_interval = as<std::size_t>(global, "checkpoint_interval", cfg.checkpoint_interval);
    cfg.checkpoint_path = as<std::string>(global, "checkpoint_path", cfg.checkpoint_path);
    cfg.restart_path = as<std::string>(global, "restart_path", cfg.restart_path);

    cfg.species.clear();
    for (const auto& block : blocks.species_blocks) {
        SpeciesConfig s;
        s.name = as<std::string>(block, "name", s.name);
        s.charge = as<double>(block, "charge", s.charge);
        s.mass = as<double>(block, "mass", s.mass);
        s.particles = as<std::size_t>(block, "particles", s.particles);
        s.density = as<double>(block, "density", s.density);
        s.drift_velocity = as<double>(block, "drift_velocity", s.drift_velocity);
        s.thermal_velocity = as<double>(block, "thermal_velocity", s.thermal_velocity);
        s.init_x_min = as<double>(block, "init_x_min", s.init_x_min);
        s.init_x_max = as<double>(block, "init_x_max", s.init_x_max);
        if (block.count("weight")) {
            s.weight = as<double>(block, "weight", s.weight);
        } else if (block.count("density")) {
            const double xmax = s.init_x_max < 0.0 ? cfg.length : s.init_x_max;
            s.weight = s.density * (xmax - s.init_x_min) / static_cast<double>(s.particles);
        }
        cfg.species.push_back(s);
    }
    if (cfg.species.empty()) cfg.species.push_back(SpeciesConfig{});
    if (cfg.checkpoint_output && cfg.checkpoint_interval == 0) cfg.checkpoint_interval = cfg.output_interval;
    validate_config(cfg);
    return cfg;
}

unsigned detect_config_dimension(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open config: " + path);
    std::string section = "global", line;
    std::size_t line_number = 0;
    while (std::getline(in, line)) {
        ++line_number;
        auto comment = line.find_first_of("#;");
        if (comment != std::string::npos) line = line.substr(0, comment);
        line = trim(line);
        if (line.empty()) continue;
        if (line.front() == '[' && line.back() == ']') {
            section = lower(trim(line.substr(1, line.size() - 2)));
            continue;
        }
        if (section != "global") continue;
        auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        const std::string key = lower(trim(line.substr(0, eq)));
        if (key != "dimension") continue;
        KeyValue kv;
        assign_key(kv, key, trim(line.substr(eq + 1)), line_number);
        const auto dimension = as<std::size_t>(kv, "dimension", 1);
        if (dimension == 1 || dimension == 2 || dimension == 3) return static_cast<unsigned>(dimension);
        throw config_error(line_number, "dimension must be 1, 2, or 3");
    }
    return 1;
}

Simulation2DConfig load_config_2d(const std::string& path) {
    static const std::unordered_set<std::string> global_keys{
        "dimension", "nx", "ny", "length_x", "length_y", "dt", "steps",
        "output_interval", "output_dir", "seed", "boundary", "vtk_output",
        "particle_output", "particle_output_interval", "particle_output_stride", "particle_sample_count",
        "checkpoint_output", "checkpoint_interval", "checkpoint_path", "restart_path",
        "magnetic_field_z",
        "particle_boundary", "particle_boundary_left", "particle_boundary_right",
        "particle_boundary_bottom", "particle_boundary_top",
        "phi_left", "phi_right", "phi_bottom", "phi_top",
        "boundary_left_tag", "boundary_right_tag", "boundary_bottom_tag", "boundary_top_tag"
    };
    static const std::unordered_set<std::string> species_keys{
        "name", "charge", "mass", "weight", "density", "particles", "drift_velocity_x",
        "drift_velocity_y", "thermal_velocity", "init_x_min", "init_x_max",
        "init_y_min", "init_y_max"
    };

    auto blocks = parse_config_blocks(path, global_keys, species_keys, nullptr, "2D");
    const auto& global = blocks.global;

    if (!global.count("dimension")) throw std::runtime_error("2D config loader requires dimension = 2");
    const auto dimension = as<std::size_t>(global, "dimension", 2);
    if (dimension != 2) throw std::runtime_error("2D config loader requires dimension = 2");

    Simulation2DConfig cfg;
    cfg.nx = as<std::size_t>(global, "nx", cfg.nx);
    cfg.ny = as<std::size_t>(global, "ny", cfg.ny);
    cfg.length_x = as<double>(global, "length_x", cfg.length_x);
    cfg.length_y = as<double>(global, "length_y", cfg.length_y);
    cfg.dt = as<double>(global, "dt", cfg.dt);
    cfg.steps = as<std::size_t>(global, "steps", cfg.steps);
    cfg.output_interval = as<std::size_t>(global, "output_interval", cfg.output_interval);
    cfg.output_dir = as<std::string>(global, "output_dir", cfg.output_dir.string());
    cfg.seed = as<unsigned>(global, "seed", cfg.seed);
    cfg.boundary = parse_boundary(global, cfg.boundary);
    cfg.vtk_output = parse_bool(global, "vtk_output", cfg.vtk_output);
    cfg.particle_output = parse_bool(global, "particle_output", cfg.particle_output);
    cfg.particle_output_interval = as<std::size_t>(global, "particle_output_interval", cfg.particle_output_interval);
    cfg.particle_output_stride = as<std::size_t>(global, "particle_output_stride", cfg.particle_output_stride);
    cfg.particle_sample_count = as<std::size_t>(global, "particle_sample_count", cfg.particle_sample_count);
    cfg.checkpoint_output = parse_bool(global, "checkpoint_output", cfg.checkpoint_output);
    cfg.checkpoint_interval = as<std::size_t>(global, "checkpoint_interval", cfg.checkpoint_interval);
    cfg.checkpoint_path = as<std::string>(global, "checkpoint_path", cfg.checkpoint_path.string());
    cfg.restart_path = as<std::string>(global, "restart_path", cfg.restart_path.string());
    cfg.magnetic_field_z = as<double>(global, "magnetic_field_z", cfg.magnetic_field_z);
    const ParticleBoundary default_particle_boundary = parse_particle_boundary(global, "particle_boundary", ParticleBoundary::Auto);
    cfg.particle_boundary_config.left = parse_particle_boundary(global, "particle_boundary_left", default_particle_boundary);
    cfg.particle_boundary_config.right = parse_particle_boundary(global, "particle_boundary_right", default_particle_boundary);
    cfg.particle_boundary_config.bottom = parse_particle_boundary(global, "particle_boundary_bottom", default_particle_boundary);
    cfg.particle_boundary_config.top = parse_particle_boundary(global, "particle_boundary_top", default_particle_boundary);
    cfg.boundary_config.left.potential = as<double>(global, "phi_left", cfg.boundary_config.left.potential);
    cfg.boundary_config.right.potential = as<double>(global, "phi_right", cfg.boundary_config.right.potential);
    cfg.boundary_config.bottom.potential = as<double>(global, "phi_bottom", cfg.boundary_config.bottom.potential);
    cfg.boundary_config.top.potential = as<double>(global, "phi_top", cfg.boundary_config.top.potential);
    cfg.boundary_config.left.tag = as<std::string>(global, "boundary_left_tag", cfg.boundary_config.left.tag);
    cfg.boundary_config.right.tag = as<std::string>(global, "boundary_right_tag", cfg.boundary_config.right.tag);
    cfg.boundary_config.bottom.tag = as<std::string>(global, "boundary_bottom_tag", cfg.boundary_config.bottom.tag);
    cfg.boundary_config.top.tag = as<std::string>(global, "boundary_top_tag", cfg.boundary_config.top.tag);

    cfg.species.clear();
    for (const auto& block : blocks.species_blocks) {
        Species2DConfig s;
        s.name = as<std::string>(block, "name", s.name);
        s.charge = as<double>(block, "charge", s.charge);
        s.mass = as<double>(block, "mass", s.mass);
        s.particles = as<std::size_t>(block, "particles", s.particles);
        s.drift_velocity_x = as<double>(block, "drift_velocity_x", s.drift_velocity_x);
        s.drift_velocity_y = as<double>(block, "drift_velocity_y", s.drift_velocity_y);
        s.thermal_velocity = as<double>(block, "thermal_velocity", s.thermal_velocity);
        s.init_x_min = as<double>(block, "init_x_min", s.init_x_min);
        s.init_x_max = as<double>(block, "init_x_max", s.init_x_max);
        s.init_y_min = as<double>(block, "init_y_min", s.init_y_min);
        s.init_y_max = as<double>(block, "init_y_max", s.init_y_max);
        if (block.count("density")) {
            const double density = as<double>(block, "density", 1.0);
            validate_positive(density, "2D species '" + s.name + "' density");
            if (!block.count("weight")) {
                const double xmax = s.init_x_max < 0.0 ? cfg.length_x : s.init_x_max;
                const double ymax = s.init_y_max < 0.0 ? cfg.length_y : s.init_y_max;
                s.weight = density * (xmax - s.init_x_min) * (ymax - s.init_y_min) / static_cast<double>(s.particles);
            }
        }
        if (block.count("weight")) s.weight = as<double>(block, "weight", s.weight);
        cfg.species.push_back(s);
    }
    if (cfg.species.empty()) cfg.species.push_back(Species2DConfig{});
    if (cfg.checkpoint_output && cfg.checkpoint_interval == 0) cfg.checkpoint_interval = cfg.output_interval;
    validate_config_2d(cfg);
    return cfg;
}

Simulation3DConfig load_config_3d(const std::string& path) {
    static const std::unordered_set<std::string> global_keys{
        "dimension", "nx", "ny", "nz", "length_x", "length_y", "length_z", "dt", "steps",
        "output_interval", "output_dir", "seed", "boundary", "vtk_output",
        "particle_output", "particle_output_interval", "particle_output_stride", "particle_sample_count",
        "checkpoint_output", "checkpoint_interval", "checkpoint_path", "restart_path",
        "magnetic_field_x", "magnetic_field_y", "magnetic_field_z",
        "particle_boundary", "particle_boundary_left", "particle_boundary_right",
        "particle_boundary_bottom", "particle_boundary_top", "particle_boundary_back", "particle_boundary_front"
    };
    static const std::unordered_set<std::string> species_keys{
        "name", "charge", "mass", "weight", "density", "particles", "drift_velocity_x",
        "drift_velocity_y", "drift_velocity_z", "thermal_velocity", "init_x_min", "init_x_max",
        "init_y_min", "init_y_max", "init_z_min", "init_z_max"
    };

    auto blocks = parse_config_blocks(path, global_keys, species_keys, nullptr, "3D");
    const auto& global = blocks.global;

    if (!global.count("dimension")) throw std::runtime_error("3D config loader requires dimension = 3");
    const auto dimension = as<std::size_t>(global, "dimension", 3);
    if (dimension != 3) throw std::runtime_error("3D config loader requires dimension = 3");

    Simulation3DConfig cfg;
    cfg.nx = as<std::size_t>(global, "nx", cfg.nx);
    cfg.ny = as<std::size_t>(global, "ny", cfg.ny);
    cfg.nz = as<std::size_t>(global, "nz", cfg.nz);
    cfg.length_x = as<double>(global, "length_x", cfg.length_x);
    cfg.length_y = as<double>(global, "length_y", cfg.length_y);
    cfg.length_z = as<double>(global, "length_z", cfg.length_z);
    cfg.dt = as<double>(global, "dt", cfg.dt);
    cfg.steps = as<std::size_t>(global, "steps", cfg.steps);
    cfg.output_interval = as<std::size_t>(global, "output_interval", cfg.output_interval);
    cfg.output_dir = as<std::string>(global, "output_dir", cfg.output_dir.string());
    cfg.seed = as<unsigned>(global, "seed", cfg.seed);
    cfg.boundary = parse_boundary(global, cfg.boundary);
    cfg.vtk_output = parse_bool(global, "vtk_output", cfg.vtk_output);
    cfg.particle_output = parse_bool(global, "particle_output", cfg.particle_output);
    cfg.particle_output_interval = as<std::size_t>(global, "particle_output_interval", cfg.particle_output_interval);
    cfg.particle_output_stride = as<std::size_t>(global, "particle_output_stride", cfg.particle_output_stride);
    cfg.particle_sample_count = as<std::size_t>(global, "particle_sample_count", cfg.particle_sample_count);
    cfg.checkpoint_output = parse_bool(global, "checkpoint_output", cfg.checkpoint_output);
    cfg.checkpoint_interval = as<std::size_t>(global, "checkpoint_interval", cfg.checkpoint_interval);
    cfg.checkpoint_path = as<std::string>(global, "checkpoint_path", cfg.checkpoint_path.string());
    cfg.restart_path = as<std::string>(global, "restart_path", cfg.restart_path.string());
    cfg.magnetic_field.x = as<double>(global, "magnetic_field_x", cfg.magnetic_field.x);
    cfg.magnetic_field.y = as<double>(global, "magnetic_field_y", cfg.magnetic_field.y);
    cfg.magnetic_field.z = as<double>(global, "magnetic_field_z", cfg.magnetic_field.z);
    const ParticleBoundary default_particle_boundary = parse_particle_boundary(global, "particle_boundary", ParticleBoundary::Auto);
    cfg.particle_boundary_config.left = parse_particle_boundary(global, "particle_boundary_left", default_particle_boundary);
    cfg.particle_boundary_config.right = parse_particle_boundary(global, "particle_boundary_right", default_particle_boundary);
    cfg.particle_boundary_config.bottom = parse_particle_boundary(global, "particle_boundary_bottom", default_particle_boundary);
    cfg.particle_boundary_config.top = parse_particle_boundary(global, "particle_boundary_top", default_particle_boundary);
    cfg.particle_boundary_config.back = parse_particle_boundary(global, "particle_boundary_back", default_particle_boundary);
    cfg.particle_boundary_config.front = parse_particle_boundary(global, "particle_boundary_front", default_particle_boundary);

    cfg.species.clear();
    for (const auto& block : blocks.species_blocks) {
        Species3DConfig s;
        s.name = as<std::string>(block, "name", s.name);
        s.charge = as<double>(block, "charge", s.charge);
        s.mass = as<double>(block, "mass", s.mass);
        s.particles = as<std::size_t>(block, "particles", s.particles);
        s.drift_velocity_x = as<double>(block, "drift_velocity_x", s.drift_velocity_x);
        s.drift_velocity_y = as<double>(block, "drift_velocity_y", s.drift_velocity_y);
        s.drift_velocity_z = as<double>(block, "drift_velocity_z", s.drift_velocity_z);
        s.thermal_velocity = as<double>(block, "thermal_velocity", s.thermal_velocity);
        s.init_x_min = as<double>(block, "init_x_min", s.init_x_min);
        s.init_x_max = as<double>(block, "init_x_max", s.init_x_max);
        s.init_y_min = as<double>(block, "init_y_min", s.init_y_min);
        s.init_y_max = as<double>(block, "init_y_max", s.init_y_max);
        s.init_z_min = as<double>(block, "init_z_min", s.init_z_min);
        s.init_z_max = as<double>(block, "init_z_max", s.init_z_max);
        if (block.count("density")) {
            const double density = as<double>(block, "density", 1.0);
            validate_positive(density, "3D species '" + s.name + "' density");
            if (!block.count("weight")) {
                const double xmax = s.init_x_max < 0.0 ? cfg.length_x : s.init_x_max;
                const double ymax = s.init_y_max < 0.0 ? cfg.length_y : s.init_y_max;
                const double zmax = s.init_z_max < 0.0 ? cfg.length_z : s.init_z_max;
                s.weight = density * (xmax - s.init_x_min) * (ymax - s.init_y_min) * (zmax - s.init_z_min)
                         / static_cast<double>(s.particles);
            }
        }
        if (block.count("weight")) s.weight = as<double>(block, "weight", s.weight);
        cfg.species.push_back(s);
    }
    if (cfg.species.empty()) cfg.species.push_back(Species3DConfig{});
    if (cfg.checkpoint_output && cfg.checkpoint_interval == 0) cfg.checkpoint_interval = cfg.output_interval;
    validate_config_3d(cfg);
    return cfg;
}
} // namespace pic
