#include "pic/Config.hpp"
#include "pic/GasDataset.hpp"
#include "pic/Runtime.hpp"
#include "pic/Simulation2D.hpp"
#include "pic/Simulation3D.hpp"
#include "pic/Units.hpp"
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

std::filesystem::path resolved_input_path(
    const std::filesystem::path& config_path,
    const std::string& value) {
    std::filesystem::path result(value);
    if (result.is_absolute()) return result;
    return config_path.parent_path() / result;
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

SpatialAverageSamplingOrder1D parse_spatial_average_sampling_order(
    const KeyValue& values,
    SpatialAverageSamplingOrder1D default_order) {
    const auto value = lower(trim(as<std::string>(
        values, "spatial_average_sampling_order",
        to_string(default_order))));
    if (value == "post_collision") {
        return SpatialAverageSamplingOrder1D::PostCollision;
    }
    if (value == "pre_collision") {
        return SpatialAverageSamplingOrder1D::PreCollision;
    }
    throw std::runtime_error(
        "invalid spatial_average_sampling_order value: '" + value + "'");
}

CollisionVelocitySampling1D parse_collision_velocity_sampling(
    const KeyValue& values,
    CollisionVelocitySampling1D default_sampling) {
    const auto value = lower(trim(as<std::string>(
        values, "collision_velocity_sampling",
        to_string(default_sampling))));
    if (value == "time_centered") {
        return CollisionVelocitySampling1D::TimeCentered;
    }
    if (value == "leapfrog_half_step") {
        return CollisionVelocitySampling1D::LeapfrogHalfStep;
    }
    throw std::runtime_error(
        "invalid collision_velocity_sampling value: '" + value + "'");
}

std::vector<PhaseEedfRegion1DConfig> parse_phase_eedf_regions(
    const KeyValue& values) {
    const auto found = values.find("phase_eedf_regions");
    if (found == values.end()) return {};
    std::vector<PhaseEedfRegion1DConfig> regions;
    std::istringstream entries(found->second);
    std::string entry;
    while (std::getline(entries, entry, ',')) {
        std::istringstream fields(trim(entry));
        PhaseEedfRegion1DConfig region;
        std::string minimum;
        std::string maximum;
        if (!std::getline(fields, region.name, ':') ||
            !std::getline(fields, minimum, ':') ||
            !std::getline(fields, maximum, ':') ||
            fields.rdbuf()->in_avail() != 0) {
            throw std::runtime_error(
                "phase_eedf_regions entries must use name:x_min:x_max");
        }
        region.name = trim(region.name);
        KeyValue coordinates{{"minimum", trim(minimum)},
                             {"maximum", trim(maximum)}};
        region.x_min = as<double>(coordinates, "minimum", 0.0);
        region.x_max = as<double>(coordinates, "maximum", 0.0);
        regions.push_back(std::move(region));
    }
    return regions;
}

std::vector<double> parse_phase_surface_flux_positions(
    const KeyValue& values) {
    const auto found = values.find("phase_surface_flux_positions");
    if (found == values.end()) return {};
    std::vector<double> positions;
    std::istringstream entries(found->second);
    std::string entry;
    while (std::getline(entries, entry, ',')) {
        KeyValue value{{"position", trim(entry)}};
        positions.push_back(as<double>(value, "position", 0.0));
    }
    return positions;
}

std::optional<std::uint64_t> parse_optional_uint64(
    const KeyValue& values, const std::string& key) {
    const auto found = values.find(key);
    if (found == values.end()) return {};
    if (found->second.empty() || found->second.front() == '-') {
        throw std::runtime_error(
            "invalid unsigned 64-bit value for '" + key + "'");
    }
    try {
        std::size_t parsed = 0;
        const int base =
            found->second.size() > 2 &&
                    found->second[0] == '0' &&
                    (found->second[1] == 'x' ||
                     found->second[1] == 'X')
                ? 16
                : 10;
        const auto value = std::stoull(
            found->second, &parsed, base);
        require_full_parse(found->second, parsed, key);
        return static_cast<std::uint64_t>(value);
    } catch (const std::exception&) {
        throw std::runtime_error(
            "invalid unsigned 64-bit value for '" + key +
            "': '" + found->second + "'");
    }
}

Boundary parse_boundary_key(
    const KeyValue& kv, const std::string& key, Boundary def) {
    const auto value =
        lower(as<std::string>(kv, key, to_string(def)));
    if (value == "periodic") return Boundary::Periodic;
    if (value == "dirichlet" || value == "absorbing") return Boundary::Dirichlet;
    throw std::runtime_error(
        "invalid boundary value for '" + key + "': '" + value + "'");
}

Boundary parse_boundary(const KeyValue& kv, Boundary def) {
    return parse_boundary_key(kv, "boundary", def);
}

VTKOutputFormat parse_vtk_output_format(const KeyValue& kv, VTKOutputFormat def) {
    const auto value = lower(trim(as<std::string>(kv, "vtk_format", to_string(def))));
    if (value == "legacy" || value == "vtk") return VTKOutputFormat::Legacy;
    if (value == "xml" || value == "vts") return VTKOutputFormat::Xml;
    if (value == "both" || value == "all") return VTKOutputFormat::Both;
    throw std::runtime_error("invalid vtk_format value: '" + value + "'");
}

RuntimeBackend parse_runtime_backend(const KeyValue& kv, RuntimeBackend def) {
    const auto value = lower(trim(as<std::string>(kv, "runtime_backend", to_string(def))));
    if (value == "serial" || value == "none" || value == "single") return RuntimeBackend::Serial;
    if (value == "openmp" || value == "omp") return RuntimeBackend::OpenMP;
    if (value == "mpi") return RuntimeBackend::MPI;
    if (value == "gpu" || value == "cuda" || value == "accelerator") return RuntimeBackend::GPU;
    throw std::runtime_error("invalid runtime_backend value: '" + value + "'");
}

RuntimePolicy parse_runtime_policy(const KeyValue& kv, const RuntimePolicy& def) {
    RuntimePolicy policy = def;
    policy.backend = parse_runtime_backend(kv, policy.backend);
    policy.threads = as<std::size_t>(kv, "runtime_threads", policy.threads);
    validate_runtime_policy(policy);
    return policy;
}

std::size_t parse_config_version(const KeyValue& kv, const std::string& loader_name) {
    const std::size_t version = as<std::size_t>(kv, "config_version", 1);
    if (version != 1) {
        throw std::runtime_error(loader_name + " config declares unsupported config_version "
                                 + std::to_string(version)
                                 + "; this AuroraPIC build supports config_version = 1");
    }
    return version;
}

RunMode parse_mode(const KeyValue& kv, RunMode def) {
    const auto value = lower(as<std::string>(kv, "mode", to_string(def)));
    if (value == "transient") return RunMode::Transient;
    if (value == "steady" || value == "steady_state") return RunMode::SteadyState;
    throw std::runtime_error("invalid mode value: '" + value + "'");
}

void parse_density_profile(
    const KeyValue& values,
    ParticleInitializationConfig& initialization,
    std::size_t spatial_dimensions) {
    if (values.count("density_profile")) {
        initialization.density_profile =
            density_profile_from_string(lower(trim(
                as<std::string>(
                    values, "density_profile", "uniform"))));
    }
    const auto optional_double =
        [&](const std::string& key,
            std::optional<double>& destination) {
            if (values.count(key)) {
                destination = as<double>(values, key, 0.0);
            }
        };
    const auto optional_size =
        [&](const std::string& key,
            std::optional<std::size_t>& destination) {
            if (values.count(key)) {
                destination =
                    as<std::size_t>(values, key, 0);
            }
        };
    optional_double(
        "profile_center_x", initialization.profile_center_x);
    optional_double(
        "profile_scale_x", initialization.profile_scale_x);
    optional_size(
        "profile_mode_x", initialization.profile_mode_x);
    if (spatial_dimensions >= 2) {
        optional_double(
            "profile_center_y",
            initialization.profile_center_y);
        optional_double(
            "profile_scale_y", initialization.profile_scale_y);
        optional_size(
            "profile_mode_y", initialization.profile_mode_y);
    }
    if (spatial_dimensions >= 3) {
        optional_double(
            "profile_center_z",
            initialization.profile_center_z);
        optional_double(
            "profile_scale_z", initialization.profile_scale_z);
        optional_size(
            "profile_mode_z", initialization.profile_mode_z);
    }
    optional_double(
        "profile_amplitude", initialization.profile_amplitude);
    optional_double(
        "profile_phase", initialization.profile_phase);
    initialization.max_profile_sampling_attempts =
        as<std::size_t>(
            values, "max_profile_sampling_attempts",
            initialization.max_profile_sampling_attempts);
}

InitializationAcceptanceConfig parse_initialization_acceptance(
    const KeyValue& values) {
    InitializationAcceptanceConfig result;
    const auto optional_tolerance =
        [&](const std::string& key,
            std::optional<double>& destination) {
            if (values.count(key)) {
                destination = as<double>(values, key, 0.0);
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
    if (values.count("initialization_charge_pairs")) {
        const std::string pair_list = trim(
            as<std::string>(
                values, "initialization_charge_pairs", ""));
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
        result, "initialization acceptance config");
    return result;
}

CollisionModelKind parse_collision_model(
    const KeyValue& kv, CollisionModelKind def) {
    const auto value = lower(trim(as<std::string>(
        kv, "model", to_string(def))));
    if (value == "bgk") return CollisionModelKind::BGK;
    if (value == "null_collision" || value == "null-collision" ||
        value == "mcc") {
        return CollisionModelKind::NullCollision;
    }
    throw std::runtime_error(
        "invalid collision model: '" + value +
        "'; expected bgk or null_collision");
}

CollisionProcessKind parse_collision_process(
    const KeyValue& kv, CollisionProcessKind def) {
    const auto value = lower(trim(as<std::string>(
        kv, "type", to_string(def))));
    if (value == "elastic") return CollisionProcessKind::Elastic;
    if (value == "excitation") return CollisionProcessKind::Excitation;
    if (value == "ionization") {
        return CollisionProcessKind::Ionization;
    }
    if (value == "charge_exchange" || value == "charge-exchange") {
        return CollisionProcessKind::ChargeExchange;
    }
    throw std::runtime_error(
        "invalid collision channel type: '" + value +
        "'; expected elastic, excitation, ionization, or "
        "charge_exchange");
}

AngularScatteringKind parse_angular_scattering(
    const KeyValue& kv, AngularScatteringKind def) {
    const auto value = lower(trim(as<std::string>(
        kv, "angular_model", to_string(def))));
    if (value == "isotropic") {
        return AngularScatteringKind::Isotropic;
    }
    if (value == "backward") {
        return AngularScatteringKind::Backward;
    }
    if (value == "henyey_greenstein" ||
        value == "henyey-greenstein") {
        return AngularScatteringKind::HenyeyGreenstein;
    }
    throw std::runtime_error(
        "invalid angular_model: '" + value +
        "'; expected isotropic, backward, or henyey_greenstein");
}

CollisionEnergyFrame parse_collision_energy_frame(
    const KeyValue& kv, CollisionEnergyFrame def) {
    const auto value = lower(trim(as<std::string>(
        kv, "energy_frame", to_string(def))));
    if (value == "projectile") {
        return CollisionEnergyFrame::Projectile;
    }
    if (value == "center_of_mass" || value == "center-of-mass" ||
        value == "centre_of_mass" || value == "centre-of-mass") {
        return CollisionEnergyFrame::CenterOfMass;
    }
    throw std::runtime_error(
        "invalid energy_frame: '" + value +
        "'; expected projectile or center_of_mass");
}

CrossSectionInterpolationKind parse_cross_section_interpolation(
    const KeyValue& kv, CrossSectionInterpolationKind def) {
    const auto value = lower(trim(as<std::string>(
        kv, "cross_section_interpolation", to_string(def))));
    if (value == "linear") {
        return CrossSectionInterpolationKind::Linear;
    }
    if (value == "lower_bin" || value == "lower-bin" ||
        value == "step") {
        return CrossSectionInterpolationKind::LowerBin;
    }
    throw std::runtime_error(
        "invalid cross_section_interpolation: '" + value +
        "'; expected linear or lower_bin");
}

InelasticTransformKind parse_inelastic_transform(
    const KeyValue& kv, InelasticTransformKind def) {
    const auto value = lower(trim(as<std::string>(
        kv, "inelastic_transform", to_string(def))));
    if (value == "heavy_target" || value == "heavy-target") {
        return InelasticTransformKind::HeavyTarget;
    }
    if (value == "finite_mass_center_of_mass" ||
        value == "finite-mass-center-of-mass" ||
        value == "finite_mass_com") {
        return InelasticTransformKind::FiniteMassCenterOfMass;
    }
    throw std::runtime_error(
        "invalid inelastic_transform: '" + value +
        "'; expected heavy_target or finite_mass_center_of_mass");
}

IonizationKinematicsKind parse_ionization_kinematics(
    const KeyValue& kv, IonizationKinematicsKind def) {
    const auto value = lower(trim(as<std::string>(
        kv, "ionization_kinematics", to_string(def))));
    if (value == "equal_energy_isotropic" || value == "equal-isotropic") {
        return IonizationKinematicsKind::EqualEnergyIsotropic;
    }
    if (value == "opal_beaty_peterson" || value == "opal") {
        return IonizationKinematicsKind::OpalBeatyPeterson;
    }
    throw std::runtime_error(
        "invalid ionization_kinematics: '" + value +
        "'; expected equal_energy_isotropic or opal_beaty_peterson");
}

UnitSystemConfig parse_units(
    const KeyValue& kv, const UnitSystemConfig& defaults) {
    UnitSystemConfig result = defaults;
    const std::string value =
        lower(trim(as<std::string>(
            kv, "units", to_string(result.system))));
    if (value == "normalized" || value == "normalised") {
        result.system = UnitSystem::Normalized;
    } else if (value == "si") {
        result.system = UnitSystem::SI;
    } else {
        throw std::runtime_error(
            "invalid units value: '" + value +
            "'; expected normalized or si");
    }
    result.relative_permittivity = as<double>(
        kv, "relative_permittivity",
        result.relative_permittivity);
    (void)result.permittivity();
    return result;
}

std::optional<TabulatedVectorField1D>
parse_magnetic_field_profile(
    const KeyValue& values,
    const std::filesystem::path& config_path) {
    const bool has_file =
        values.count("magnetic_field_profile_file") != 0;
    const bool has_axis =
        values.count("magnetic_field_profile_axis") != 0;
    if (!has_file && has_axis) {
        throw std::runtime_error(
            "magnetic_field_profile_axis requires magnetic_field_profile_file");
    }
    if (!has_file) return std::nullopt;
    if (!has_axis) {
        throw std::runtime_error(
            "magnetic_field_profile_file requires magnetic_field_profile_axis");
    }
    const auto file = resolved_input_path(
        config_path,
        as<std::string>(
            values, "magnetic_field_profile_file", ""));
    const auto axis = parse_coordinate_axis(
        as<std::string>(
            values, "magnetic_field_profile_axis", ""));
    return load_tabulated_vector_field_1d(file, axis);
}

ParticleBoundary parse_particle_boundary(const KeyValue& kv, const std::string& key, ParticleBoundary def) {
    const auto value = lower(as<std::string>(kv, key, to_string(def)));
    if (value == "auto") return ParticleBoundary::Auto;
    if (value == "absorbing" || value == "absorb") return ParticleBoundary::Absorbing;
    if (value == "reflecting" || value == "reflective" || value == "reflect") return ParticleBoundary::Reflecting;
    if (value == "periodic" || value == "wrap") return ParticleBoundary::Periodic;
    throw std::runtime_error("invalid particle boundary value for '" + key + "': '" + value + "'");
}

BoundarySide2DName parse_boundary_side_2d(
    const KeyValue& values, const std::string& key,
    BoundarySide2DName default_value) {
    const auto value = lower(trim(as<std::string>(
        values, key, to_string(default_value))));
    if (value == "left") return BoundarySide2DName::Left;
    if (value == "right") return BoundarySide2DName::Right;
    if (value == "bottom") return BoundarySide2DName::Bottom;
    if (value == "top") return BoundarySide2DName::Top;
    throw std::runtime_error(
        "invalid 2D boundary side for '" + key + "': '" +
        value + "'");
}

CurrentSourceControlMode parse_current_source_control_mode(
    const KeyValue& values, const std::string& key,
    CurrentSourceControlMode default_value) {
    const auto value = lower(trim(as<std::string>(
        values, key, to_string(default_value))));
    if (value == "cumulative") {
        return CurrentSourceControlMode::Cumulative;
    }
    if (value == "timestep_local" || value == "timestep-local") {
        return CurrentSourceControlMode::TimestepLocal;
    }
    throw std::runtime_error(
        "invalid current-source control mode for '" + key + "': '" +
        value + "'");
}

PotentialReferenceCorrection parse_potential_reference_correction(
    const KeyValue& values, const std::string& key,
    PotentialReferenceCorrection default_value) {
    const auto value = lower(trim(as<std::string>(
        values, key, to_string(default_value))));
    if (value == "gauge") {
        return PotentialReferenceCorrection::Gauge;
    }
    if (value == "affine") {
        return PotentialReferenceCorrection::Affine;
    }
    throw std::runtime_error(
        "invalid potential-reference correction for '" + key + "': '" +
        value + "'");
}

void validate_positive(double value, const std::string& name) {
    if (!std::isfinite(value) || !(value > 0.0)) throw std::runtime_error(name + " must be positive and finite");
}
void validate_non_negative(double value, const std::string& name) {
    if (!std::isfinite(value) || value < 0.0) throw std::runtime_error(name + " must be non-negative and finite");
}

void require_species_scale_source(const KeyValue& block, const std::string& species_name, const std::string& dimension_label) {
    if (!block.count("weight") && !block.count("density")) {
        const std::string prefix = dimension_label.empty() ? "species" : dimension_label + " species";
        throw std::runtime_error(prefix + " '" + species_name + "' must specify positive weight or density");
    }
}

void validate_config(const Config& cfg) {
    (void)cfg.units.permittivity();
    if (!cfg.restart_path.empty() &&
        !cfg.initial_state_path.empty()) {
        throw std::runtime_error(
            "restart_path and initial_state_path are mutually exclusive");
    }
    if (cfg.initial_state_signature &&
        cfg.initial_state_path.empty()) {
        throw std::runtime_error(
            "initial_state_signature requires initial_state_path");
    }
    if (cfg.nx < 3) throw std::runtime_error("nx must be at least 3");
    if (cfg.velocity_dimensions != 1 &&
        cfg.velocity_dimensions != 3) {
        throw std::runtime_error(
            "velocity_dimensions must be 1 or 3");
    }
    validate_positive(cfg.length, "length");
    validate_positive(cfg.dt, "dt");
    if (!std::isfinite(cfg.phi_left) ||
        !std::isfinite(cfg.phi_right)) {
        throw std::runtime_error(
            "1D electrode potential offsets must be finite");
    }
    const auto validate_voltage_drive =
        [&](const SinusoidalVoltageConfig& drive,
            const std::string& name) {
            if (!std::isfinite(drive.amplitude)) {
                throw std::runtime_error(
                    name + " amplitude must be finite");
            }
            validate_non_negative(
                drive.frequency, name + " frequency");
            if (!std::isfinite(drive.phase)) {
                throw std::runtime_error(
                    name + " phase must be finite");
            }
            if (drive.amplitude != 0.0 &&
                !(drive.frequency > 0.0)) {
                throw std::runtime_error(
                    name + " nonzero amplitude requires positive "
                    "frequency");
            }
        };
    validate_voltage_drive(
        cfg.phi_left_drive, "phi_left sinusoidal drive");
    validate_voltage_drive(
        cfg.phi_right_drive, "phi_right sinusoidal drive");
    const bool driven =
        cfg.phi_left_drive.amplitude != 0.0 ||
        cfg.phi_right_drive.amplitude != 0.0;
    if (cfg.boundary != Boundary::Dirichlet && driven) {
        throw std::runtime_error(
            "sinusoidal electrode drives require boundary = "
            "dirichlet");
    }
    if (cfg.mode == RunMode::SteadyState && driven) {
        throw std::runtime_error(
            "sinusoidal electrode drives require mode = transient "
            "until cycle-averaged convergence is implemented");
    }
    if (cfg.output_interval == 0) throw std::runtime_error("output_interval must be positive");
    validate_spatial_average_1d(cfg);
    validate_positive(cfg.steady_tolerance, "steady_tolerance");
    if (cfg.steady_window == 0) throw std::runtime_error("steady_window must be positive");
    if (cfg.mode == RunMode::SteadyState && cfg.max_steps == 0) throw std::runtime_error("max_steps must be positive for steady-state mode");
    if (cfg.max_particles_per_species == 0) {
        throw std::runtime_error(
            "max_particles_per_species must be positive");
    }
    if (!cfg.collision_models.empty() &&
        cfg.collisions.enabled) {
        throw std::runtime_error(
            "named collision models cannot be combined with the "
            "legacy collision model");
    }
    const auto species_config =
        [&](const std::string& name) -> const SpeciesConfig* {
        const auto found = std::find_if(
            cfg.species.begin(), cfg.species.end(),
            [&](const SpeciesConfig& species) {
                return species.name == name;
            });
        return found == cfg.species.end() ? nullptr : &*found;
    };
    const auto validate_collision =
        [&](const CollisionConfig& collision,
            const std::string& context,
            bool named) {
        validate_non_negative(
            collision.frequency, context + " frequency");
        validate_non_negative(
            collision.neutral_temperature_velocity,
            context + " neutral_temperature_velocity");
        validate_non_negative(
            collision.neutral_density,
            context + " neutral_density");
        validate_non_negative(
            collision.max_frequency,
            context + " max_frequency");
        if (collision.max_candidates_per_particle == 0) {
            throw std::runtime_error(
                context +
                " max_candidates_per_particle must be positive");
        }
        if (!collision.enabled) return;
        if (named &&
            collision.model !=
                CollisionModelKind::NullCollision) {
            throw std::runtime_error(
                context +
                " named models support only null_collision");
        }
        if (collision.model !=
            CollisionModelKind::NullCollision) {
            return;
        }
        validate_positive(
            collision.neutral_density,
            context + " neutral_density");
        validate_positive(
            collision.max_frequency,
            context + " max_frequency");
        if (collision.species.empty()) {
            throw std::runtime_error(
                context + " requires a target species");
        }
        const auto* target =
            species_config(collision.species);
        if (!target) {
            throw std::runtime_error(
                context + " target species does not exist: " +
                collision.species);
        }
        if (collision.channels.empty()) {
            throw std::runtime_error(
                context + " requires at least one channel");
        }
        std::unordered_set<std::string> channel_names;
        for (const auto& channel : collision.channels) {
            const std::string channel_context =
                context + " channel '" + channel.name + "'";
            if (channel.name.empty()) {
                throw std::runtime_error(
                    context +
                    " channel name must not be empty");
            }
            if (!std::all_of(
                    channel.name.begin(), channel.name.end(),
                    [](unsigned char character) {
                        return std::isalnum(character) ||
                               character == '_' ||
                               character == '-' ||
                               character == '.';
                    })) {
                throw std::runtime_error(
                    context +
                    " channel names contain an invalid character");
            }
            if (!channel_names.insert(channel.name).second) {
                throw std::runtime_error(
                    context + " duplicate channel name: " +
                    channel.name);
            }
            if (channel.cross_section_file.empty()) {
                throw std::runtime_error(
                    channel_context +
                    " requires cross_section_file");
            }
            validate_non_negative(
                channel.threshold_energy,
                channel_context + " threshold_energy");
            validate_positive(
                channel.energy_scale,
                channel_context + " energy_scale");
            validate_positive(
                channel.cross_section_scale,
                channel_context + " cross_section_scale");
            if (channel.angular_scattering !=
                    AngularScatteringKind::Isotropic &&
                cfg.velocity_dimensions != 3) {
                throw std::runtime_error(
                    channel_context +
                    " anisotropic scattering requires "
                    "velocity_dimensions = 3");
            }
            if (channel.angular_scattering ==
                AngularScatteringKind::HenyeyGreenstein) {
                if (channel.mean_cosine_file.empty()) {
                    throw std::runtime_error(
                        channel_context +
                        " Henyey-Greenstein scattering requires "
                        "mean_cosine_file");
                }
            } else if (!channel.mean_cosine_file.empty() ||
                       channel.mean_cosine_energy_scale != 1.0) {
                throw std::runtime_error(
                    channel_context +
                    " mean-cosine data requires angular_model = "
                    "henyey_greenstein");
            }
            if (channel.process ==
                CollisionProcessKind::Attachment) {
                throw std::runtime_error(
                    channel_context +
                    " process is not supported in 1D");
            }
            if (channel.process ==
                    CollisionProcessKind::ChargeExchange &&
                cfg.velocity_dimensions != 3) {
                throw std::runtime_error(
                    channel_context +
                    " charge exchange requires "
                    "velocity_dimensions = 3");
            }
            if (channel.energy_frame ==
                    CollisionEnergyFrame::CenterOfMass &&
                !(collision.neutral_mass > 0.0)) {
                throw std::runtime_error(
                    channel_context +
                    " center_of_mass energy requires positive "
                    "neutral_mass");
            }
            if (channel.energy_frame ==
                    CollisionEnergyFrame::CenterOfMass &&
                channel.process != CollisionProcessKind::Elastic &&
                channel.process !=
                    CollisionProcessKind::ChargeExchange) {
                throw std::runtime_error(
                    channel_context +
                    " center_of_mass energy is supported only for "
                    "elastic and charge-exchange channels");
            }
            if (channel.process ==
                CollisionProcessKind::Elastic) {
                if (channel.threshold_energy != 0.0) {
                    throw std::runtime_error(
                        channel_context +
                        " elastic threshold_energy must be zero");
                }
            } else if (
                channel.process ==
                    CollisionProcessKind::Excitation ||
                channel.process ==
                    CollisionProcessKind::Ionization) {
                if (!(channel.threshold_energy > 0.0)) {
                    throw std::runtime_error(
                        channel_context +
                        " inelastic threshold_energy must be positive");
                }
            }
            if (channel.process ==
                CollisionProcessKind::Ionization) {
                if (channel.ionization_kinematics ==
                        IonizationKinematicsKind::OpalBeatyPeterson) {
                    validate_positive(
                        channel.ionization_ejected_energy_scale,
                        channel_context +
                            " ionization_ejected_energy_scale");
                } else if (
                    channel.ionization_ejected_energy_scale != 0.0) {
                    throw std::runtime_error(
                        channel_context +
                        " ionization_ejected_energy_scale requires "
                        "ionization_kinematics = opal_beaty_peterson");
                }
                if (cfg.velocity_dimensions != 3) {
                    throw std::runtime_error(
                        channel_context +
                        " ionization requires velocity_dimensions = 3");
                }
                const auto* secondary =
                    species_config(channel.secondary_species);
                const auto* ion =
                    species_config(channel.ion_species);
                if (!secondary || !ion) {
                    throw std::runtime_error(
                        channel_context +
                        " references an unknown product species");
                }
                if (target->charge == 0.0 ||
                    target->weight != secondary->weight ||
                    target->weight != ion->weight ||
                    target->mass != secondary->mass ||
                    target->charge != secondary->charge ||
                    ion->charge != -target->charge) {
                    throw std::runtime_error(
                        channel_context +
                        " requires a charged target, equal macro "
                        "weights, a secondary with target mass/charge, "
                        "and an oppositely charged ion");
                }
            } else if (
                !channel.secondary_species.empty() ||
                !channel.ion_species.empty() ||
                channel.ionization_kinematics !=
                    IonizationKinematicsKind::EqualEnergyIsotropic ||
                channel.ionization_ejected_energy_scale != 0.0) {
                throw std::runtime_error(
                    channel_context +
                    " product species and ionization kinematics are "
                    "valid only for ionization");
            }
            const bool inelastic =
                channel.process == CollisionProcessKind::Excitation ||
                channel.process == CollisionProcessKind::Ionization;
            if (channel.inelastic_transform !=
                    InelasticTransformKind::HeavyTarget &&
                (!inelastic || !(collision.neutral_mass > 0.0))) {
                throw std::runtime_error(
                    channel_context +
                    " finite-mass inelastic transform requires an "
                    "excitation or ionization channel and positive "
                    "neutral_mass");
            }
        }
    };
    validate_collision(
        cfg.collisions, "legacy collision model", false);
    std::unordered_set<std::string> model_names;
    std::unordered_set<std::string> target_names;
    for (const auto& named : cfg.collision_models) {
        if (named.name.empty() ||
            !model_names.insert(named.name).second) {
            throw std::runtime_error(
                "named collision model names must be non-empty and "
                "unique");
        }
        if (!std::all_of(
                named.name.begin(), named.name.end(),
                [](unsigned char character) {
                    return std::isalnum(character) ||
                           character == '_' ||
                           character == '-';
                })) {
            throw std::runtime_error(
                "named collision model names may contain only "
                "letters, digits, '_', and '-'");
        }
        validate_collision(
            named.config,
            "collision model '" + named.name + "'", true);
        if (named.config.enabled &&
            !target_names.insert(
                named.config.species).second) {
            throw std::runtime_error(
                "multiple enabled collision models target species '" +
                named.config.species + "'");
        }
    }
    validate_runtime_policy(cfg.runtime);
    validate_initialization_acceptance(
        cfg.initialization_acceptance,
        "1D initialization acceptance config");
    std::unordered_set<std::string> species_names;
    for (const auto& s : cfg.species) {
        if (s.name.empty()) throw std::runtime_error("species name must not be empty");
        if (!species_names.insert(s.name).second) {
            throw std::runtime_error(
                "1D species names must be unique: " + s.name);
        }
        validate_positive(s.mass, "species '" + s.name + "' mass");
        validate_positive(s.weight, "species '" + s.name + "' weight");
        validate_positive(s.density, "species '" + s.name + "' density");
        if (!std::isfinite(s.charge)) throw std::runtime_error("species '" + s.name + "' charge must be finite");
        if (s.particles == 0) throw std::runtime_error("species '" + s.name + "' particles must be positive");
        if (s.timestep_multiplier == 0) {
            throw std::runtime_error(
                "species '" + s.name +
                "' timestep_multiplier must be positive");
        }
        if (s.particles > cfg.max_particles_per_species) {
            throw std::runtime_error(
                "species '" + s.name +
                "' initial particles exceed "
                "max_particles_per_species");
        }
        if (!std::isfinite(s.drift_velocity) ||
            !std::isfinite(s.drift_velocity_y) ||
            !std::isfinite(s.drift_velocity_z)) {
            throw std::runtime_error(
                "species '" + s.name +
                "' drift velocities must be finite");
        }
        if (cfg.velocity_dimensions == 1 &&
            (s.drift_velocity_y != 0.0 ||
             s.drift_velocity_z != 0.0)) {
            throw std::runtime_error(
                "transverse drift velocities require "
                "velocity_dimensions = 3");
        }
        validate_non_negative(s.thermal_velocity, "species '" + s.name + "' thermal_velocity");
        validate_particle_initialization(
            s.initialization, cfg.velocity_dimensions,
            s.thermal_velocity,
            "species '" + s.name + "'");
        validate_density_profile(
            s.initialization, 1, s.particles,
            "species '" + s.name + "'");
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
    (void)cfg.units.permittivity();
    if (!cfg.restart_path.empty() &&
        !cfg.initial_state_path.empty()) {
        throw std::runtime_error(
            "2D restart_path and initial_state_path are mutually exclusive");
    }
    if (cfg.initial_state_signature &&
        cfg.initial_state_path.empty()) {
        throw std::runtime_error(
            "2D initial_state_signature requires initial_state_path");
    }
    if (cfg.nx < 3) throw std::runtime_error("2D nx must be at least 3");
    if (cfg.ny < 3) throw std::runtime_error("2D ny must be at least 3");
    validate_positive(cfg.length_x, "length_x");
    validate_positive(cfg.length_y, "length_y");
    validate_positive(
        cfg.out_of_plane_depth, "out_of_plane_depth");
    validate_positive(cfg.dt, "dt");
    if (cfg.output_interval == 0) throw std::runtime_error("output_interval must be positive");
    validate_positive(cfg.steady_tolerance, "steady_tolerance");
    if (cfg.steady_window == 0) throw std::runtime_error("steady_window must be positive");
    if (cfg.mode == RunMode::SteadyState && cfg.max_steps == 0) {
        throw std::runtime_error("max_steps must be positive for steady-state mode");
    }
    if (cfg.particle_output_stride == 0) throw std::runtime_error("particle_output_stride must be positive");
    if (cfg.max_particles_per_species == 0) {
        throw std::runtime_error(
            "2D max_particles_per_species must be positive");
    }
    if (!std::isfinite(cfg.magnetic_field_x) ||
        !std::isfinite(cfg.magnetic_field_y) ||
        !std::isfinite(cfg.magnetic_field_z)) {
        throw std::runtime_error(
            "2D magnetic_field components must be finite");
    }
    if (cfg.magnetic_field_profile) {
        if (cfg.magnetic_field_x != 0.0 ||
            cfg.magnetic_field_y != 0.0 ||
            cfg.magnetic_field_z != 0.0) {
            throw std::runtime_error(
                "2D uniform magnetic_field components and magnetic field profile are mutually exclusive");
        }
        cfg.magnetic_field_profile->validate_domain(
            {0.0, 0.0, 0.0},
            {cfg.length_x, cfg.length_y, 0.0},
            "2D config");
    }
    validate_runtime_policy(cfg.runtime);
    validate_initialization_acceptance(
        cfg.initialization_acceptance,
        "2D initialization acceptance config");
    validate_boundary_side(cfg.boundary_config.left, "left");
    validate_boundary_side(cfg.boundary_config.right, "right");
    validate_boundary_side(cfg.boundary_config.bottom, "bottom");
    validate_boundary_side(cfg.boundary_config.top, "top");
    const Boundary boundary_x =
        cfg.boundary_x.value_or(cfg.boundary);
    const Boundary boundary_y =
        cfg.boundary_y.value_or(cfg.boundary);
    if (boundary_x == Boundary::Periodic &&
        (cfg.boundary_config.left.potential != 0.0 ||
         cfg.boundary_config.right.potential != 0.0)) {
        throw std::runtime_error(
            "2D phi_left/phi_right must be zero when boundary_x is periodic");
    }
    if (boundary_y == Boundary::Periodic &&
        (cfg.boundary_config.bottom.potential != 0.0 ||
         cfg.boundary_config.top.potential != 0.0)) {
        throw std::runtime_error(
            "2D phi_bottom/phi_top must be zero when boundary_y is periodic");
    }
    if (cfg.potential_reference) {
        const auto& reference = *cfg.potential_reference;
        if (reference.axis == CoordinateAxis::Z ||
            !std::isfinite(reference.coordinate) ||
            !std::isfinite(reference.target)) {
            throw std::runtime_error(
                "2D potential reference requires a finite x/y coordinate and target");
        }
        const double length =
            reference.axis == CoordinateAxis::X
                ? cfg.length_x
                : cfg.length_y;
        if (reference.coordinate < 0.0 ||
            reference.coordinate > length) {
            throw std::runtime_error(
                "2D potential reference coordinate lies outside the domain");
        }
    }
    for (const auto& s : cfg.species) {
        if (s.name.empty()) throw std::runtime_error("2D species name must not be empty");
        validate_positive(s.mass, "2D species '" + s.name + "' mass");
        validate_positive(s.weight, "2D species '" + s.name + "' weight");
        if (!std::isfinite(s.charge)) throw std::runtime_error("2D species '" + s.name + "' charge must be finite");
        if (s.particles == 0) throw std::runtime_error("2D species '" + s.name + "' particles must be positive");
        if (!std::isfinite(s.drift_velocity_x) ||
            !std::isfinite(s.drift_velocity_y) ||
            !std::isfinite(s.drift_velocity_z)) {
            throw std::runtime_error("2D species '" + s.name + "' drift velocities must be finite");
        }
        validate_non_negative(s.thermal_velocity, "2D species '" + s.name + "' thermal_velocity");
        validate_particle_initialization(
            s.initialization, 3, s.thermal_velocity,
            "2D species '" + s.name + "'");
        validate_density_profile(
            s.initialization, 2, s.particles,
            "2D species '" + s.name + "'");
        if (s.init_x_min < 0.0) throw std::runtime_error("2D species '" + s.name + "' init_x_min must be non-negative");
        if (s.init_y_min < 0.0) throw std::runtime_error("2D species '" + s.name + "' init_y_min must be non-negative");
        const double xmax = s.init_x_max < 0.0 ? cfg.length_x : s.init_x_max;
        const double ymax = s.init_y_max < 0.0 ? cfg.length_y : s.init_y_max;
        if (xmax > cfg.length_x) throw std::runtime_error("2D species '" + s.name + "' init_x_max exceeds domain length_x");
        if (ymax > cfg.length_y) throw std::runtime_error("2D species '" + s.name + "' init_y_max exceeds domain length_y");
        if (!(s.init_x_min < xmax)) throw std::runtime_error("2D species '" + s.name + "' x initialization interval must have positive width");
        if (!(s.init_y_min < ymax)) throw std::runtime_error("2D species '" + s.name + "' y initialization interval must have positive width");
    }
    if (cfg.current_regulated_source) {
        const auto& source = *cfg.current_regulated_source;
        if (source.species.empty() ||
            std::none_of(
                cfg.species.begin(), cfg.species.end(),
                [&](const Species2DConfig& species) {
                    return species.name == source.species;
                })) {
            throw std::runtime_error(
                "2D current-regulated source must reference a configured species");
        }
        const double normal_length =
            source.emission_boundary == BoundarySide2DName::Left ||
                    source.emission_boundary == BoundarySide2DName::Right
                ? cfg.length_x
                : cfg.length_y;
        if (!std::isfinite(source.emission_inset) ||
            source.emission_inset < 0.0 ||
            !(source.emission_inset < normal_length) ||
            !std::isfinite(source.drift.x) ||
            !std::isfinite(source.drift.y) ||
            !std::isfinite(source.drift.z) ||
            !std::isfinite(source.thermal_velocity) ||
            source.thermal_velocity < 0.0) {
            throw std::runtime_error(
                "2D current-regulated source has invalid position or velocity controls");
        }
    }
    for (const auto& source : cfg.sources) {
        const std::string label =
            "2D source '" + source.name + "'";
        if (source.name.empty() || source.first_species.empty() ||
            source.second_species.empty()) {
            throw std::runtime_error(
                label + " must name both source species");
        }
        if (!std::all_of(
                source.name.begin(), source.name.end(),
                [](unsigned char character) {
                    return std::isalnum(character) ||
                           character == '_' || character == '-';
                })) {
            throw std::runtime_error(
                "2D source names may contain only letters, digits, '_' and '-'");
        }
        if (source.first_species == source.second_species) {
            throw std::runtime_error(
                label + " source species must be distinct");
        }
        const bool fixed_rate = source.pairs_per_step != 0;
        const bool physical_rate =
            source.represented_pair_rate.has_value();
        const bool volumetric_rate =
            source.peak_volumetric_pair_rate.has_value();
        if (static_cast<unsigned>(fixed_rate) +
                static_cast<unsigned>(physical_rate) +
                static_cast<unsigned>(volumetric_rate) !=
            1U) {
            throw std::runtime_error(
                label + " requires exactly one of positive pairs_per_step, represented_pair_rate, or peak_volumetric_pair_rate");
        }
        if (source.represented_pair_rate) {
            validate_positive(
                *source.represented_pair_rate,
                label + " represented_pair_rate");
        }
        if (source.peak_volumetric_pair_rate) {
            validate_positive(
                *source.peak_volumetric_pair_rate,
                label + " peak_volumetric_pair_rate");
        }
        if (source.end_step != 0 &&
            source.end_step <= source.start_step) {
            throw std::runtime_error(
                label + " end_step must exceed start_step");
        }
        const double xmax =
            source.x_max < 0.0 ? cfg.length_x : source.x_max;
        const double ymax =
            source.y_max < 0.0 ? cfg.length_y : source.y_max;
        if (source.x_min < 0.0 || source.y_min < 0.0 ||
            xmax > cfg.length_x || ymax > cfg.length_y ||
            !(source.x_min < xmax) || !(source.y_min < ymax)) {
            throw std::runtime_error(
                label + " region must have positive area inside the domain");
        }
        const auto finite_vec = [](Vec3 value) {
            return std::isfinite(value.x) &&
                   std::isfinite(value.y) &&
                   std::isfinite(value.z);
        };
        if (!finite_vec(source.first_drift) ||
            !finite_vec(source.second_drift)) {
            throw std::runtime_error(
                label + " drift velocities must be finite");
        }
        validate_non_negative(
            source.first_thermal_velocity,
            label + " first_thermal_velocity");
        validate_non_negative(
            source.second_thermal_velocity,
            label + " second_thermal_velocity");
        validate_density_profile(
            source.spatial_profile, 2, 1,
            label + " spatial profile");
    }
}

void validate_config_3d(const Simulation3DConfig& cfg) {
    (void)cfg.units.permittivity();
    if (!cfg.restart_path.empty() &&
        !cfg.initial_state_path.empty()) {
        throw std::runtime_error(
            "3D restart_path and initial_state_path are mutually exclusive");
    }
    if (cfg.initial_state_signature &&
        cfg.initial_state_path.empty()) {
        throw std::runtime_error(
            "3D initial_state_signature requires initial_state_path");
    }
    if (cfg.nx < 3) throw std::runtime_error("3D nx must be at least 3");
    if (cfg.ny < 3) throw std::runtime_error("3D ny must be at least 3");
    if (cfg.nz < 3) throw std::runtime_error("3D nz must be at least 3");
    validate_positive(cfg.length_x, "length_x");
    validate_positive(cfg.length_y, "length_y");
    validate_positive(cfg.length_z, "length_z");
    validate_positive(cfg.dt, "dt");
    if (cfg.output_interval == 0) throw std::runtime_error("output_interval must be positive");
    validate_positive(cfg.steady_tolerance, "steady_tolerance");
    if (cfg.steady_window == 0) throw std::runtime_error("steady_window must be positive");
    if (cfg.mode == RunMode::SteadyState && cfg.max_steps == 0) {
        throw std::runtime_error("max_steps must be positive for steady-state mode");
    }
    if (cfg.particle_output_stride == 0) throw std::runtime_error("particle_output_stride must be positive");
    if (!std::isfinite(cfg.magnetic_field.x) || !std::isfinite(cfg.magnetic_field.y) || !std::isfinite(cfg.magnetic_field.z)) {
        throw std::runtime_error("magnetic_field components must be finite");
    }
    if (cfg.magnetic_field_profile) {
        if (cfg.magnetic_field.x != 0.0 ||
            cfg.magnetic_field.y != 0.0 ||
            cfg.magnetic_field.z != 0.0) {
            throw std::runtime_error(
                "3D uniform magnetic_field components and magnetic field profile are mutually exclusive");
        }
        cfg.magnetic_field_profile->validate_domain(
            {0.0, 0.0, 0.0},
            {cfg.length_x, cfg.length_y, cfg.length_z},
            "3D config");
    }
    validate_runtime_policy(cfg.runtime);
    validate_initialization_acceptance(
        cfg.initialization_acceptance,
        "3D initialization acceptance config");
    for (const auto& s : cfg.species) {
        if (s.name.empty()) throw std::runtime_error("3D species name must not be empty");
        validate_positive(s.mass, "3D species '" + s.name + "' mass");
        validate_positive(s.weight, "3D species '" + s.name + "' weight");
        if (!std::isfinite(s.charge)) throw std::runtime_error("3D species '" + s.name + "' charge must be finite");
        if (s.particles == 0) throw std::runtime_error("3D species '" + s.name + "' particles must be positive");
        if (!std::isfinite(s.drift_velocity_x) || !std::isfinite(s.drift_velocity_y) || !std::isfinite(s.drift_velocity_z)) {
            throw std::runtime_error("3D species '" + s.name + "' drift velocities must be finite");
        }
        validate_non_negative(s.thermal_velocity, "3D species '" + s.name + "' thermal_velocity");
        validate_particle_initialization(
            s.initialization, 3, s.thermal_velocity,
            "3D species '" + s.name + "'");
        validate_density_profile(
            s.initialization, 3, s.particles,
            "3D species '" + s.name + "'");
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
    struct NamedBlock {
        std::string name{};
        KeyValue values{};
    };
    struct CollisionModelBlock {
        std::string name{};
        KeyValue values{};
        std::vector<KeyValue> channel_blocks{};
    };
    KeyValue global;
    KeyValue collisions;
    std::vector<KeyValue> collision_channel_blocks;
    std::vector<CollisionModelBlock> collision_model_blocks;
    std::vector<KeyValue> species_blocks;
    std::vector<NamedBlock> source_blocks;
};

ParsedBlocks parse_config_blocks(const std::string& path,
                                 const std::unordered_set<std::string>& global_keys,
                                 const std::unordered_set<std::string>& species_keys,
                                 const std::unordered_set<std::string>* collision_keys,
                                 const std::unordered_set<std::string>*
                                     collision_channel_keys,
                                 const std::unordered_set<std::string>*
                                     source_keys,
                                 const std::string& loader_name) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open config: " + path);

    ParsedBlocks blocks;
    std::string section = "global", line;
    KeyValue* current = &blocks.global;
    const std::unordered_set<std::string>* allowed = &global_keys;
    std::size_t line_number = 0;
    const auto collision_model_block =
        [&](const std::string& name) -> ParsedBlocks::CollisionModelBlock& {
        const auto found = std::find_if(
            blocks.collision_model_blocks.begin(),
            blocks.collision_model_blocks.end(),
            [&](const auto& block) {
                return block.name == name;
            });
        if (found != blocks.collision_model_blocks.end()) {
            return *found;
        }
        blocks.collision_model_blocks.push_back({name, {}, {}});
        return blocks.collision_model_blocks.back();
    };
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
            } else if (collision_keys &&
                       starts_with(section, "collisions.")) {
                const std::string suffix =
                    trim(section.substr(
                        std::string("collisions.").size()));
                const std::string channel_marker = ".channel.";
                const auto marker = suffix.find(channel_marker);
                if (suffix.empty()) {
                    throw config_error(
                        line_number,
                        "empty named collision model suffix");
                }
                if (marker == std::string::npos) {
                    if (suffix.find('.') != std::string::npos) {
                        throw config_error(
                            line_number,
                            "named collision model section must be "
                            "[collisions.<model>] or "
                            "[collisions.<model>.channel.<channel>]");
                    }
                    auto& model = collision_model_block(suffix);
                    current = &model.values;
                    allowed = collision_keys;
                } else {
                    const std::string model_name =
                        trim(suffix.substr(0, marker));
                    const std::string channel_name =
                        trim(suffix.substr(
                            marker + channel_marker.size()));
                    if (model_name.empty() || channel_name.empty() ||
                        model_name.find('.') != std::string::npos ||
                        channel_name.find('.') != std::string::npos) {
                        throw config_error(
                            line_number,
                            "invalid named collision channel section");
                    }
                    auto& model =
                        collision_model_block(model_name);
                    model.channel_blocks.emplace_back();
                    current = &model.channel_blocks.back();
                    allowed = collision_channel_keys;
                    assign_key(
                        *current, "name", channel_name,
                        line_number);
                }
            } else if (collision_channel_keys &&
                       (starts_with(section, "collision.") ||
                        starts_with(section, "collision_channel."))) {
                const std::string prefix =
                    starts_with(section, "collision.")
                        ? "collision."
                        : "collision_channel.";
                const std::string channel_name =
                    trim(section.substr(prefix.size()));
                if (channel_name.empty()) {
                    throw config_error(
                        line_number,
                        "empty collision channel section suffix");
                }
                blocks.collision_channel_blocks.emplace_back();
                current = &blocks.collision_channel_blocks.back();
                allowed = collision_channel_keys;
                assign_key(
                    *current, "name", channel_name, line_number);
            } else if (section == "species" || starts_with(section, "species.")) {
                blocks.species_blocks.emplace_back();
                current = &blocks.species_blocks.back();
                allowed = &species_keys;
                if (starts_with(section, "species.")) {
                    const std::string species_name = trim(section.substr(std::string("species.").size()));
                    if (species_name.empty()) throw config_error(line_number, "empty species section suffix");
                    assign_key(*current, "name", species_name, line_number);
                }
            } else if (source_keys && starts_with(section, "source.")) {
                const std::string source_name =
                    trim(section.substr(std::string("source.").size()));
                if (source_name.empty() ||
                    source_name.find('.') != std::string::npos) {
                    throw config_error(
                        line_number, "invalid source section suffix");
                }
                const auto duplicate = std::find_if(
                    blocks.source_blocks.begin(),
                    blocks.source_blocks.end(),
                    [&](const auto& block) {
                        return block.name == source_name;
                    });
                if (duplicate != blocks.source_blocks.end()) {
                    throw config_error(
                        line_number,
                        "duplicate source section '" + source_name + "'");
                }
                blocks.source_blocks.push_back({source_name, {}});
                current = &blocks.source_blocks.back().values;
                allowed = source_keys;
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

void validate_spatial_average_1d(const Config& cfg) {
    const auto& average = cfg.spatial_average;
    const auto& eedf = cfg.phase_eedf;
    const auto& flux = cfg.phase_surface_flux;
    if (!eedf.enabled) {
        if (eedf.history_enabled || !eedf.species.empty() ||
            eedf.energy_bins != 0 ||
            eedf.energy_max != 0.0 || eedf.tail_threshold != 0.0 ||
            !eedf.regions.empty()) {
            throw std::runtime_error(
                "disabled phase_eedf cannot configure history, species, "
                "bins, energy maximum, tail threshold, or regions");
        }
    } else {
        if (!average.enabled || average.phase_bins == 0) {
            throw std::runtime_error(
                "phase_eedf requires spatial averaging with phase bins");
        }
        if (eedf.species.empty() || eedf.energy_bins == 0 ||
            eedf.energy_bins > 1000000 ||
            !std::isfinite(eedf.energy_max) || eedf.energy_max <= 0.0 ||
            eedf.regions.empty()) {
            throw std::runtime_error(
                "phase_eedf requires a species, 1..1000000 bins, a "
                "positive energy maximum, and at least one region");
        }
        if (!std::isfinite(eedf.tail_threshold) ||
            eedf.tail_threshold < 0.0 ||
            eedf.tail_threshold >= eedf.energy_max) {
            throw std::runtime_error(
                "phase_eedf tail threshold must be finite, non-negative, "
                "and below its energy maximum");
        }
        if (std::none_of(cfg.species.begin(), cfg.species.end(),
                         [&](const auto& species) {
                             return species.name == eedf.species;
                         })) {
            throw std::runtime_error(
                "phase_eedf species does not exist: " + eedf.species);
        }
        std::unordered_set<std::string> names;
        for (const auto& region : eedf.regions) {
            if (region.name.empty() || !names.insert(region.name).second ||
                !std::all_of(region.name.begin(), region.name.end(),
                    [](unsigned char c) {
                        return std::isalnum(c) || c == '_' || c == '-';
                    }) ||
                !std::isfinite(region.x_min) ||
                !std::isfinite(region.x_max) ||
                region.x_min < 0.0 || region.x_max > cfg.length ||
                region.x_max <= region.x_min) {
                throw std::runtime_error(
                    "phase_eedf regions require unique safe names and "
                    "0 <= x_min < x_max <= length");
            }
        }
    }
    if (!flux.enabled) {
        if (flux.reset_on_restart || !flux.species.empty() ||
            !flux.positions.empty() || flux.energy_bins != 0 ||
            flux.energy_max != 0.0) {
            throw std::runtime_error(
                "disabled phase_surface_flux cannot configure restart "
                "reset, species, positions, bins, or energy maximum");
        }
    } else {
        if (!average.enabled || average.phase_bins == 0 ||
            average.rf_frequency <= 0.0) {
            throw std::runtime_error(
                "phase_surface_flux requires phase-resolved spatial averaging");
        }
        if (cfg.boundary == Boundary::Periodic) {
            throw std::runtime_error(
                "phase_surface_flux does not yet support periodic boundaries");
        }
        if (flux.species.empty() || flux.positions.empty() ||
            flux.energy_bins == 0 || flux.energy_bins > 1000000 ||
            !std::isfinite(flux.energy_max) || flux.energy_max <= 0.0 ||
            std::none_of(cfg.species.begin(), cfg.species.end(),
                         [&](const auto& species) {
                             return species.name == flux.species;
                         })) {
            throw std::runtime_error(
                "phase_surface_flux requires an existing species, at least "
                "one position, 1..1000000 bins, and a positive energy maximum");
        }
        double previous = 0.0;
        for (std::size_t index = 0; index < flux.positions.size(); ++index) {
            const double position = flux.positions[index];
            if (!std::isfinite(position) || position <= 0.0 ||
                position >= cfg.length ||
                (index != 0 && position <= previous)) {
                throw std::runtime_error(
                    "phase_surface_flux_positions must be finite, strictly "
                    "increasing, and internal to the domain");
            }
            previous = position;
        }
    }
    if (!std::isfinite(average.rf_frequency) ||
        average.rf_frequency < 0.0) {
        throw std::runtime_error(
            "spatial_average_rf_frequency must be non-negative "
            "and finite");
    }
    if (!average.enabled) {
        if (average.reset_on_restart ||
            average.rf_cycles != 0 ||
            average.rf_frequency != 0.0 ||
            average.phase_bins != 0 ||
            average.sampling_order !=
                SpatialAverageSamplingOrder1D::PostCollision) {
            throw std::runtime_error(
                "disabled spatial_average cannot configure restart "
                "reset, sampling order, or an RF contract");
        }
        return;
    }
    if (average.interval == 0) {
        throw std::runtime_error(
            "spatial_average_interval must be positive");
    }
    if (average.start_step == 0 ||
        average.end_step < average.start_step) {
        throw std::runtime_error(
            "spatial_average requires 1 <= start_step <= end_step");
    }
    if (average.rf_frequency == 0.0) {
        if (average.rf_cycles != 0 || average.phase_bins != 0) {
            throw std::runtime_error(
                "spatial_average_rf_cycles requires positive "
                "spatial_average_rf_frequency");
        }
        return;
    }
    if (average.rf_cycles == 0) {
        throw std::runtime_error(
            "positive spatial_average_rf_frequency requires positive "
            "spatial_average_rf_cycles");
    }
    const double steps_per_cycle_value =
        1.0 / (average.rf_frequency * cfg.dt);
    if (!std::isfinite(steps_per_cycle_value) ||
        steps_per_cycle_value >
            static_cast<double>(
                std::numeric_limits<long long>::max())) {
        throw std::runtime_error(
            "spatial-average RF steps per cycle are invalid");
    }
    const auto steps_per_cycle = static_cast<std::size_t>(
        std::llround(steps_per_cycle_value));
    if (steps_per_cycle == 0 ||
        std::abs(steps_per_cycle_value -
                 static_cast<double>(steps_per_cycle)) >
            1e-10 * std::max(1.0, steps_per_cycle_value)) {
        throw std::runtime_error(
            "spatial-average RF period must contain an integer "
            "number of timesteps");
    }
    if (average.rf_cycles >
        std::numeric_limits<std::size_t>::max() /
            steps_per_cycle) {
        throw std::runtime_error(
            "spatial-average RF window step count overflows");
    }
    const std::size_t window_steps =
        average.end_step - average.start_step + 1;
    if (window_steps !=
        average.rf_cycles * steps_per_cycle) {
        throw std::runtime_error(
            "spatial-average step window does not span the "
            "configured whole RF cycles");
    }
    if (average.end_step % steps_per_cycle != 0) {
        throw std::runtime_error(
            "spatial-average RF window must end at the time-zero "
            "drive phase");
    }
    if (steps_per_cycle % average.interval != 0) {
        throw std::runtime_error(
            "spatial_average_interval must divide RF steps per cycle");
    }
    const std::size_t samples_per_cycle =
        steps_per_cycle / average.interval;
    if (average.phase_bins > samples_per_cycle ||
        (average.phase_bins != 0 &&
         samples_per_cycle % average.phase_bins != 0)) {
        throw std::runtime_error(
            "spatial_average_phase_bins must divide the sampled RF cycle");
    }
    for (const auto& drive :
         {cfg.phi_left_drive, cfg.phi_right_drive}) {
        if (drive.amplitude == 0.0) continue;
        if (std::abs(
                drive.frequency - average.rf_frequency) >
            1e-12 * std::max(
                drive.frequency, average.rf_frequency)) {
            throw std::runtime_error(
                "spatial-average RF frequency does not match "
                "the electrode drive");
        }
    }
}

Config load_config(const std::string& path) {
    static const std::unordered_set<std::string> global_keys{
        "nx", "length", "velocity_dimensions", "dt", "steps",
        "output_interval", "output_dir", "seed",
        "spatial_average", "spatial_average_reset_on_restart",
        "spatial_average_interval",
        "spatial_average_start_step", "spatial_average_end_step",
        "spatial_average_rf_frequency", "spatial_average_rf_cycles",
        "spatial_average_phase_bins", "spatial_average_sampling_order",
        "phase_eedf", "phase_eedf_history", "phase_eedf_species",
        "phase_eedf_energy_bins",
        "phase_eedf_energy_max", "phase_eedf_tail_threshold",
        "phase_eedf_regions",
        "phase_surface_flux", "phase_surface_flux_reset_on_restart",
        "phase_surface_flux_species", "phase_surface_flux_positions",
        "phase_surface_flux_energy_bins", "phase_surface_flux_energy_max",
        "wall_impact_spectrum", "wall_impact_reset_on_restart",
        "wall_impact_energy_bins",
        "wall_impact_energy_max",
        "max_particles_per_species",
        "collision_velocity_sampling",
        "phi_left", "phi_right", "steady_tolerance", "steady_window", "max_steps",
        "phi_left_amplitude", "phi_left_frequency", "phi_left_phase",
        "phi_right_amplitude", "phi_right_frequency", "phi_right_phase",
        "boundary", "mode", "dimension", "config_version", "checkpoint_output", "checkpoint_interval",
        "checkpoint_path", "restart_path", "initial_state_path",
        "initial_state_signature",
        "runtime_backend", "runtime_threads",
        "units", "relative_permittivity",
        "initialization_max_relative_charge_imbalance",
        "initialization_max_relative_current_imbalance",
        "initialization_max_relative_pair_imbalance",
        "initialization_charge_pairs"
    };
    static const std::unordered_set<std::string> collision_keys{
        "enabled", "model", "frequency", "neutral_temperature_velocity",
        "neutral_density", "species", "max_frequency",
        "max_candidates_per_particle", "neutral_mass",
        "neutral_temperature", "gas_data_file"
    };
    static const std::unordered_set<std::string> collision_channel_keys{
        "name", "type", "cross_section_file", "threshold_energy",
        "energy_scale", "cross_section_scale",
        "secondary_species", "ion_species", "angular_model",
        "mean_cosine_file", "mean_cosine_energy_scale",
        "energy_frame", "ionization_kinematics",
        "ionization_ejected_energy_scale",
        "cross_section_interpolation", "inelastic_transform"
    };
    static const std::unordered_set<std::string> species_keys{
        "name", "charge", "mass", "weight", "particles", "density",
        "drift_velocity", "drift_velocity_y", "drift_velocity_z",
        "thermal_velocity", "thermal_velocity_x",
        "thermal_velocity_y", "thermal_velocity_z",
        "initialization_version",
        "loading", "density_profile", "profile_center_x",
        "profile_scale_x", "profile_amplitude", "profile_phase",
        "profile_mode_x", "max_profile_sampling_attempts",
        "init_x_min", "init_x_max", "timestep_multiplier"
    };

    auto blocks = parse_config_blocks(
        path, global_keys, species_keys, &collision_keys,
        &collision_channel_keys, nullptr, "1D");
    const auto& global = blocks.global;
    const auto& collision = blocks.collisions;

    (void)parse_config_version(global, "1D");
    const auto dimension = as<std::size_t>(global, "dimension", 1);
    if (dimension != 1) throw std::runtime_error("1D config loader requires dimension = 1 or no dimension key");

    Config cfg;
    cfg.units = parse_units(global, cfg.units);
    cfg.nx = as<std::size_t>(global, "nx", cfg.nx);
    cfg.velocity_dimensions = as<std::size_t>(
        global, "velocity_dimensions", cfg.velocity_dimensions);
    cfg.length = as<double>(global, "length", cfg.length);
    cfg.dt = as<double>(global, "dt", cfg.dt);
    cfg.steps = as<std::size_t>(global, "steps", cfg.steps);
    cfg.output_interval = as<std::size_t>(global, "output_interval", cfg.output_interval);
    cfg.spatial_average.enabled = parse_bool(
        global, "spatial_average", cfg.spatial_average.enabled);
    cfg.spatial_average.reset_on_restart = parse_bool(
        global, "spatial_average_reset_on_restart",
        cfg.spatial_average.reset_on_restart);
    cfg.spatial_average.interval = as<std::size_t>(
        global, "spatial_average_interval",
        cfg.spatial_average.interval);
    cfg.spatial_average.start_step = as<std::size_t>(
        global, "spatial_average_start_step",
        cfg.spatial_average.start_step);
    cfg.spatial_average.end_step = as<std::size_t>(
        global, "spatial_average_end_step",
        cfg.spatial_average.end_step);
    cfg.spatial_average.rf_frequency = as<double>(
        global, "spatial_average_rf_frequency",
        cfg.spatial_average.rf_frequency);
    cfg.spatial_average.rf_cycles = as<std::size_t>(
        global, "spatial_average_rf_cycles",
        cfg.spatial_average.rf_cycles);
    cfg.spatial_average.phase_bins = as<std::size_t>(
        global, "spatial_average_phase_bins",
        cfg.spatial_average.phase_bins);
    cfg.spatial_average.sampling_order =
        parse_spatial_average_sampling_order(
            global, cfg.spatial_average.sampling_order);
    cfg.phase_eedf.enabled = parse_bool(
        global, "phase_eedf", cfg.phase_eedf.enabled);
    cfg.phase_eedf.history_enabled = parse_bool(
        global, "phase_eedf_history",
        cfg.phase_eedf.history_enabled);
    cfg.phase_eedf.species = as<std::string>(
        global, "phase_eedf_species", cfg.phase_eedf.species);
    cfg.phase_eedf.energy_bins = as<std::size_t>(
        global, "phase_eedf_energy_bins", cfg.phase_eedf.energy_bins);
    cfg.phase_eedf.energy_max = as<double>(
        global, "phase_eedf_energy_max", cfg.phase_eedf.energy_max);
    cfg.phase_eedf.tail_threshold = as<double>(
        global, "phase_eedf_tail_threshold",
        cfg.phase_eedf.tail_threshold);
    cfg.phase_eedf.regions = parse_phase_eedf_regions(global);
    cfg.phase_surface_flux.enabled = parse_bool(
        global, "phase_surface_flux", cfg.phase_surface_flux.enabled);
    cfg.phase_surface_flux.reset_on_restart = parse_bool(
        global, "phase_surface_flux_reset_on_restart",
        cfg.phase_surface_flux.reset_on_restart);
    cfg.phase_surface_flux.species = as<std::string>(
        global, "phase_surface_flux_species",
        cfg.phase_surface_flux.species);
    cfg.phase_surface_flux.positions =
        parse_phase_surface_flux_positions(global);
    cfg.phase_surface_flux.energy_bins = as<std::size_t>(
        global, "phase_surface_flux_energy_bins",
        cfg.phase_surface_flux.energy_bins);
    cfg.phase_surface_flux.energy_max = as<double>(
        global, "phase_surface_flux_energy_max",
        cfg.phase_surface_flux.energy_max);
    cfg.wall_impact_spectrum.enabled = parse_bool(
        global, "wall_impact_spectrum",
        cfg.wall_impact_spectrum.enabled);
    cfg.wall_impact_spectrum.reset_on_restart = parse_bool(
        global, "wall_impact_reset_on_restart",
        cfg.wall_impact_spectrum.reset_on_restart);
    cfg.wall_impact_spectrum.energy_bins = as<std::size_t>(
        global, "wall_impact_energy_bins",
        cfg.wall_impact_spectrum.energy_bins);
    cfg.wall_impact_spectrum.energy_max = as<double>(
        global, "wall_impact_energy_max",
        cfg.wall_impact_spectrum.energy_max);
    cfg.output_dir = as<std::string>(global, "output_dir", cfg.output_dir);
    cfg.seed = as<unsigned>(global, "seed", cfg.seed);
    cfg.phi_left = as<double>(global, "phi_left", cfg.phi_left);
    cfg.phi_right = as<double>(global, "phi_right", cfg.phi_right);
    cfg.phi_left_drive.amplitude = as<double>(
        global, "phi_left_amplitude",
        cfg.phi_left_drive.amplitude);
    cfg.phi_left_drive.frequency = as<double>(
        global, "phi_left_frequency",
        cfg.phi_left_drive.frequency);
    cfg.phi_left_drive.phase = as<double>(
        global, "phi_left_phase", cfg.phi_left_drive.phase);
    cfg.phi_right_drive.amplitude = as<double>(
        global, "phi_right_amplitude",
        cfg.phi_right_drive.amplitude);
    cfg.phi_right_drive.frequency = as<double>(
        global, "phi_right_frequency",
        cfg.phi_right_drive.frequency);
    cfg.phi_right_drive.phase = as<double>(
        global, "phi_right_phase", cfg.phi_right_drive.phase);
    cfg.steady_tolerance = as<double>(global, "steady_tolerance", cfg.steady_tolerance);
    cfg.steady_window = as<std::size_t>(global, "steady_window", cfg.steady_window);
    cfg.max_steps = as<std::size_t>(global, "max_steps", cfg.max_steps);
    cfg.max_particles_per_species = as<std::size_t>(
        global, "max_particles_per_species",
        cfg.max_particles_per_species);
    cfg.collision_velocity_sampling =
        parse_collision_velocity_sampling(
            global, cfg.collision_velocity_sampling);
    cfg.boundary = parse_boundary(global, cfg.boundary);
    cfg.mode = parse_mode(global, cfg.mode);
    cfg.collisions.enabled = parse_bool(collision, "enabled", cfg.collisions.enabled);
    cfg.collisions.model =
        parse_collision_model(collision, cfg.collisions.model);
    cfg.collisions.frequency = as<double>(collision, "frequency", cfg.collisions.frequency);
    cfg.collisions.neutral_temperature_velocity = as<double>(collision, "neutral_temperature_velocity", cfg.collisions.neutral_temperature_velocity);
    cfg.collisions.neutral_density =
        as<double>(
            collision, "neutral_density",
            cfg.collisions.neutral_density);
    cfg.collisions.species =
        as<std::string>(
            collision, "species", cfg.collisions.species);
    cfg.collisions.max_frequency =
        as<double>(
            collision, "max_frequency",
            cfg.collisions.max_frequency);
    cfg.collisions.max_candidates_per_particle =
        as<std::size_t>(
            collision, "max_candidates_per_particle",
            cfg.collisions.max_candidates_per_particle);
    cfg.collisions.neutral_mass =
        as<double>(
            collision, "neutral_mass",
            cfg.collisions.neutral_mass);
    cfg.collisions.neutral_temperature =
        as<double>(
            collision, "neutral_temperature",
            cfg.collisions.neutral_temperature);
    if (collision.count("gas_data_file")) {
        throw std::runtime_error(
            "1D gas_data_file requires a named "
            "[collisions.<model>] section");
    }
    cfg.collisions.gas_data_units = cfg.units.system;
    cfg.checkpoint_output = parse_bool(global, "checkpoint_output", cfg.checkpoint_output);
    cfg.checkpoint_interval = as<std::size_t>(global, "checkpoint_interval", cfg.checkpoint_interval);
    cfg.checkpoint_path = as<std::string>(global, "checkpoint_path", cfg.checkpoint_path);
    cfg.restart_path = as<std::string>(global, "restart_path", cfg.restart_path);
    if (global.count("initial_state_path")) {
        cfg.initial_state_path = resolved_input_path(
            path, as<std::string>(
                      global, "initial_state_path", ""));
    }
    cfg.initial_state_signature = parse_optional_uint64(
        global, "initial_state_signature");
    cfg.runtime = parse_runtime_policy(global, cfg.runtime);
    cfg.initialization_acceptance =
        parse_initialization_acceptance(global);

    cfg.species.clear();
    for (const auto& block : blocks.species_blocks) {
        SpeciesConfig s;
        s.name = as<std::string>(block, "name", s.name);
        s.charge = as<double>(block, "charge", s.charge);
        s.mass = as<double>(block, "mass", s.mass);
        s.particles = as<std::size_t>(block, "particles", s.particles);
        s.density = as<double>(block, "density", s.density);
        s.drift_velocity = as<double>(block, "drift_velocity", s.drift_velocity);
        s.drift_velocity_y = as<double>(
            block, "drift_velocity_y", s.drift_velocity_y);
        s.drift_velocity_z = as<double>(
            block, "drift_velocity_z", s.drift_velocity_z);
        s.thermal_velocity = as<double>(block, "thermal_velocity", s.thermal_velocity);
        s.initialization.version = as<std::size_t>(
            block, "initialization_version", s.initialization.version);
        if (block.count("loading")) {
            s.initialization.loading = particle_loading_from_string(
                as<std::string>(block, "loading", "random"));
        }
        if (block.count("thermal_velocity_x")) {
            s.initialization.thermal_velocity_x =
                as<double>(block, "thermal_velocity_x", 0.0);
        }
        if (block.count("thermal_velocity_y")) {
            s.initialization.thermal_velocity_y =
                as<double>(block, "thermal_velocity_y", 0.0);
        }
        if (block.count("thermal_velocity_z")) {
            s.initialization.thermal_velocity_z =
                as<double>(block, "thermal_velocity_z", 0.0);
        }
        parse_density_profile(block, s.initialization, 1);
        s.init_x_min = as<double>(block, "init_x_min", s.init_x_min);
        s.init_x_max = as<double>(block, "init_x_max", s.init_x_max);
        s.timestep_multiplier = as<std::size_t>(
            block, "timestep_multiplier", s.timestep_multiplier);
        require_species_scale_source(block, s.name, "");
        if (block.count("weight")) {
            s.weight = as<double>(block, "weight", s.weight);
        } else if (block.count("density")) {
            const double xmax = s.init_x_max < 0.0 ? cfg.length : s.init_x_max;
            s.weight = s.density * (xmax - s.init_x_min) / static_cast<double>(s.particles);
        }
        cfg.species.push_back(s);
    }
    if (cfg.species.empty()) cfg.species.push_back(SpeciesConfig{});
    cfg.collisions.channels.clear();
    const auto config_directory =
        std::filesystem::absolute(std::filesystem::path(path))
            .parent_path();
    const auto parse_channel =
        [&](const KeyValue& block) {
        CollisionChannelConfig channel;
        channel.name =
            as<std::string>(block, "name", channel.name);
        channel.process =
            parse_collision_process(block, channel.process);
        const auto cross_section_file =
            as<std::string>(block, "cross_section_file", "");
        if (!cross_section_file.empty()) {
            const std::filesystem::path configured(cross_section_file);
            channel.cross_section_file =
                (configured.is_absolute()
                     ? configured
                     : config_directory / configured)
                    .lexically_normal();
        }
        channel.threshold_energy =
            as<double>(
                block, "threshold_energy",
                channel.threshold_energy);
        channel.energy_scale =
            as<double>(
                block, "energy_scale", channel.energy_scale);
        channel.cross_section_scale =
            as<double>(
                block, "cross_section_scale",
                channel.cross_section_scale);
        channel.cross_section_interpolation =
            parse_cross_section_interpolation(
                block, channel.cross_section_interpolation);
        channel.inelastic_transform =
            parse_inelastic_transform(
                block, channel.inelastic_transform);
        channel.angular_scattering =
            parse_angular_scattering(
                block, channel.angular_scattering);
        channel.energy_frame =
            parse_collision_energy_frame(
                block, channel.energy_frame);
        channel.ionization_kinematics =
            parse_ionization_kinematics(
                block, channel.ionization_kinematics);
        channel.ionization_ejected_energy_scale =
            as<double>(
                block, "ionization_ejected_energy_scale",
                channel.ionization_ejected_energy_scale);
        const auto mean_cosine_file =
            as<std::string>(block, "mean_cosine_file", "");
        if (!mean_cosine_file.empty()) {
            const std::filesystem::path configured(mean_cosine_file);
            channel.mean_cosine_file =
                (configured.is_absolute()
                     ? configured
                     : config_directory / configured)
                    .lexically_normal();
        }
        channel.mean_cosine_energy_scale =
            as<double>(
                block, "mean_cosine_energy_scale",
                channel.mean_cosine_energy_scale);
        channel.secondary_species =
            as<std::string>(
                block, "secondary_species",
                channel.secondary_species);
        channel.ion_species =
            as<std::string>(
                block, "ion_species",
                channel.ion_species);
        return channel;
    };
    for (const auto& block : blocks.collision_channel_blocks) {
        cfg.collisions.channels.push_back(
            parse_channel(block));
    }
    if (!blocks.collision_model_blocks.empty() &&
        (!blocks.collisions.empty() ||
         !blocks.collision_channel_blocks.empty())) {
        throw std::runtime_error(
            "named [collisions.<model>] sections cannot be mixed with "
            "legacy [collisions] or [collision.<channel>] sections");
    }
    cfg.collision_models.clear();
    for (const auto& block : blocks.collision_model_blocks) {
        NamedCollisionConfig named;
        named.name = block.name;
        named.config.enabled =
            parse_bool(block.values, "enabled", true);
        named.config.model =
            parse_collision_model(
                block.values,
                CollisionModelKind::NullCollision);
        named.config.frequency =
            as<double>(
                block.values, "frequency",
                named.config.frequency);
        named.config.neutral_temperature_velocity =
            as<double>(
                block.values, "neutral_temperature_velocity",
                named.config.neutral_temperature_velocity);
        named.config.neutral_density =
            as<double>(
                block.values, "neutral_density",
                named.config.neutral_density);
        named.config.species =
            as<std::string>(
                block.values, "species",
                named.config.species);
        named.config.max_frequency =
            as<double>(
                block.values, "max_frequency",
                named.config.max_frequency);
        named.config.max_candidates_per_particle =
            as<std::size_t>(
                block.values, "max_candidates_per_particle",
                named.config.max_candidates_per_particle);
        named.config.neutral_mass =
            as<double>(
                block.values, "neutral_mass",
                named.config.neutral_mass);
        named.config.neutral_temperature =
            as<double>(
                block.values, "neutral_temperature",
                named.config.neutral_temperature);
        named.config.gas_data_units = cfg.units.system;
        const auto gas_data_file =
            as<std::string>(
                block.values, "gas_data_file", "");
        if (!gas_data_file.empty()) {
            if (block.values.count("neutral_mass")) {
                throw std::runtime_error(
                    "named collision model '" + named.name +
                    "' gas_data_file cannot be combined with "
                    "neutral_mass");
            }
            const std::filesystem::path configured(gas_data_file);
            named.config.gas_data_file =
                (configured.is_absolute()
                     ? configured
                     : config_directory / configured)
                    .lexically_normal();
            const auto dataset =
                load_gas_dataset(named.config.gas_data_file);
            if (dataset.unit_system != cfg.units.system) {
                throw std::runtime_error(
                    "named collision model '" + named.name +
                    "' gas dataset units '" +
                    to_string(dataset.unit_system) +
                    "' do not match simulation units '" +
                    to_string(cfg.units.system) + "'");
            }
            named.config.gas_data_version =
                dataset.format_version;
            named.config.gas_data_units =
                dataset.unit_system;
            named.config.gas_name = dataset.gas_name;
            named.config.neutral_mass =
                dataset.neutral_mass;
            named.config.data_provenance =
                dataset.data_provenance;
            named.config.dataset_id = dataset.dataset_id;
            named.config.dataset_version =
                dataset.dataset_version;
            named.config.citation = dataset.citation;
            named.config.retrieved = dataset.retrieved;
            named.config.license = dataset.license;
            named.config.channels = dataset.channels;
            for (const auto& override : block.channel_blocks) {
                const std::string channel_name =
                    as<std::string>(override, "name", "");
                const auto found = std::find_if(
                    named.config.channels.begin(),
                    named.config.channels.end(),
                    [&](const auto& channel) {
                        return channel.name == channel_name;
                    });
                if (found == named.config.channels.end()) {
                    throw std::runtime_error(
                        "named collision model '" + named.name +
                        "' gas dataset has no channel '" +
                        channel_name + "'");
                }
                for (const auto& [key, unused] : override) {
                    (void)unused;
                    if (key != "name" &&
                        key != "secondary_species" &&
                        key != "ion_species") {
                        throw std::runtime_error(
                            "named collision model '" + named.name +
                            "' gas dataset channel '" + channel_name +
                            "' permits only reactive product mappings");
                    }
                }
                if (found->process !=
                    CollisionProcessKind::Ionization) {
                    throw std::runtime_error(
                        "named collision model '" + named.name +
                        "' gas dataset channel '" + channel_name +
                        "' is not ionization and cannot map products");
                }
                found->secondary_species =
                    as<std::string>(
                        override, "secondary_species",
                        found->secondary_species);
                found->ion_species =
                    as<std::string>(
                        override, "ion_species",
                        found->ion_species);
            }
        } else {
            for (const auto& channel : block.channel_blocks) {
                named.config.channels.push_back(
                    parse_channel(channel));
            }
        }
        cfg.collision_models.push_back(std::move(named));
    }
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
        "dimension", "config_version", "nx", "ny", "length_x", "length_y",
        "out_of_plane_depth", "dt", "steps",
        "mode", "steady_tolerance", "steady_window", "max_steps",
        "output_interval", "output_dir", "seed", "boundary",
        "boundary_x", "boundary_y", "vtk_output", "vtk_format",
        "particle_output", "particle_output_interval", "particle_output_stride", "particle_sample_count",
        "resolved_diagnostics", "resolved_diagnostic_interval",
        "resolved_diagnostic_start_step", "resolved_profile_axis",
        "resolved_mode_axis", "resolved_max_mode",
        "max_particles_per_species",
        "checkpoint_output", "checkpoint_interval", "checkpoint_path",
        "restart_path", "initial_state_path",
        "initial_state_signature",
        "runtime_backend", "runtime_threads", "magnetic_field_x",
        "magnetic_field_y", "magnetic_field_z",
        "magnetic_field_profile_file", "magnetic_field_profile_axis",
        "particle_boundary", "particle_boundary_left", "particle_boundary_right",
        "particle_boundary_bottom", "particle_boundary_top",
        "phi_left", "phi_right", "phi_bottom", "phi_top",
        "boundary_left_tag", "boundary_right_tag", "boundary_bottom_tag",
        "boundary_top_tag", "units", "relative_permittivity",
        "current_source_species", "current_source_monitor_boundary",
        "current_source_control_mode",
        "current_source_emission_boundary", "current_source_emission_inset",
        "current_source_drift_velocity_x",
        "current_source_drift_velocity_y",
        "current_source_drift_velocity_z",
        "current_source_thermal_velocity",
        "current_source_temperature_ev",
        "potential_reference_axis", "potential_reference_coordinate",
        "potential_reference_target", "potential_reference_correction",
        "initialization_max_relative_charge_imbalance",
        "initialization_max_relative_current_imbalance",
        "initialization_max_relative_pair_imbalance",
        "initialization_charge_pairs"
    };
    static const std::unordered_set<std::string> species_keys{
        "name", "charge", "mass", "weight", "density", "particles", "drift_velocity_x",
        "drift_velocity_y", "drift_velocity_z", "thermal_velocity",
        "temperature_ev",
        "thermal_velocity_x", "thermal_velocity_y", "thermal_velocity_z",
        "initialization_version", "loading",
        "density_profile", "profile_center_x", "profile_center_y",
        "profile_scale_x", "profile_scale_y", "profile_amplitude",
        "profile_phase", "profile_mode_x", "profile_mode_y",
        "max_profile_sampling_attempts",
        "init_x_min", "init_x_max",
        "init_y_min", "init_y_max"
    };
    static const std::unordered_set<std::string> source_keys{
        "first_species", "second_species", "pairs_per_step",
        "represented_pair_rate", "peak_volumetric_pair_rate",
        "start_step", "end_step", "x_min", "x_max", "y_min", "y_max",
        "first_drift_velocity_x", "first_drift_velocity_y",
        "first_drift_velocity_z", "second_drift_velocity_x",
        "second_drift_velocity_y", "second_drift_velocity_z",
        "first_thermal_velocity", "second_thermal_velocity",
        "first_temperature_ev", "second_temperature_ev",
        "density_profile", "profile_center_x", "profile_center_y",
        "profile_scale_x", "profile_scale_y", "profile_amplitude",
        "profile_phase", "profile_mode_x", "profile_mode_y",
        "max_profile_sampling_attempts"
    };

    auto blocks = parse_config_blocks(
        path, global_keys, species_keys, nullptr, nullptr,
        &source_keys, "2D");
    const auto& global = blocks.global;

    (void)parse_config_version(global, "2D");
    if (!global.count("dimension")) throw std::runtime_error("2D config loader requires dimension = 2");
    const auto dimension = as<std::size_t>(global, "dimension", 2);
    if (dimension != 2) throw std::runtime_error("2D config loader requires dimension = 2");

    Simulation2DConfig cfg;
    cfg.units = parse_units(global, cfg.units);
    cfg.nx = as<std::size_t>(global, "nx", cfg.nx);
    cfg.ny = as<std::size_t>(global, "ny", cfg.ny);
    cfg.length_x = as<double>(global, "length_x", cfg.length_x);
    cfg.length_y = as<double>(global, "length_y", cfg.length_y);
    cfg.out_of_plane_depth = as<double>(
        global, "out_of_plane_depth",
        cfg.out_of_plane_depth);
    cfg.dt = as<double>(global, "dt", cfg.dt);
    cfg.steps = as<std::size_t>(global, "steps", cfg.steps);
    cfg.mode = parse_mode(global, cfg.mode);
    cfg.steady_tolerance = as<double>(global, "steady_tolerance", cfg.steady_tolerance);
    cfg.steady_window = as<std::size_t>(global, "steady_window", cfg.steady_window);
    cfg.max_steps = as<std::size_t>(global, "max_steps", cfg.max_steps);
    cfg.output_interval = as<std::size_t>(global, "output_interval", cfg.output_interval);
    cfg.output_dir = as<std::string>(global, "output_dir", cfg.output_dir.string());
    cfg.seed = as<unsigned>(global, "seed", cfg.seed);
    cfg.boundary = parse_boundary(global, cfg.boundary);
    if (global.count("boundary_x")) {
        cfg.boundary_x =
            parse_boundary_key(global, "boundary_x", cfg.boundary);
    }
    if (global.count("boundary_y")) {
        cfg.boundary_y =
            parse_boundary_key(global, "boundary_y", cfg.boundary);
    }
    cfg.vtk_output = parse_bool(global, "vtk_output", cfg.vtk_output);
    cfg.vtk_format = parse_vtk_output_format(global, cfg.vtk_format);
    cfg.particle_output = parse_bool(global, "particle_output", cfg.particle_output);
    cfg.particle_output_interval = as<std::size_t>(global, "particle_output_interval", cfg.particle_output_interval);
    cfg.particle_output_stride = as<std::size_t>(global, "particle_output_stride", cfg.particle_output_stride);
    cfg.particle_sample_count = as<std::size_t>(global, "particle_sample_count", cfg.particle_sample_count);
    cfg.resolved_diagnostics.enabled = parse_bool(
        global, "resolved_diagnostics",
        cfg.resolved_diagnostics.enabled);
    cfg.resolved_diagnostics.interval = as<std::size_t>(
        global, "resolved_diagnostic_interval",
        cfg.resolved_diagnostics.interval);
    cfg.resolved_diagnostics.start_step = as<std::size_t>(
        global, "resolved_diagnostic_start_step",
        cfg.resolved_diagnostics.start_step);
    if (global.count("resolved_profile_axis")) {
        cfg.resolved_diagnostics.profile_axis =
            parse_coordinate_axis(as<std::string>(
                global, "resolved_profile_axis", "x"));
    }
    if (global.count("resolved_mode_axis")) {
        cfg.resolved_diagnostics.mode_axis =
            parse_coordinate_axis(as<std::string>(
                global, "resolved_mode_axis", "y"));
    }
    cfg.resolved_diagnostics.max_mode = as<std::size_t>(
        global, "resolved_max_mode",
        cfg.resolved_diagnostics.max_mode);
    cfg.checkpoint_output = parse_bool(global, "checkpoint_output", cfg.checkpoint_output);
    cfg.checkpoint_interval = as<std::size_t>(global, "checkpoint_interval", cfg.checkpoint_interval);
    cfg.checkpoint_path = as<std::string>(global, "checkpoint_path", cfg.checkpoint_path.string());
    cfg.restart_path = as<std::string>(global, "restart_path", cfg.restart_path.string());
    if (global.count("initial_state_path")) {
        cfg.initial_state_path = resolved_input_path(
            path, as<std::string>(
                      global, "initial_state_path", ""));
    }
    cfg.initial_state_signature = parse_optional_uint64(
        global, "initial_state_signature");
    cfg.max_particles_per_species = as<std::size_t>(
        global, "max_particles_per_species",
        cfg.max_particles_per_species);
    cfg.runtime = parse_runtime_policy(global, cfg.runtime);
    cfg.initialization_acceptance =
        parse_initialization_acceptance(global);
    cfg.magnetic_field_x =
        as<double>(
            global, "magnetic_field_x", cfg.magnetic_field_x);
    cfg.magnetic_field_y =
        as<double>(
            global, "magnetic_field_y", cfg.magnetic_field_y);
    cfg.magnetic_field_z = as<double>(global, "magnetic_field_z", cfg.magnetic_field_z);
    cfg.magnetic_field_profile =
        parse_magnetic_field_profile(global, path);
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
    std::optional<double> current_source_temperature_ev;
    const bool has_current_source = std::any_of(
        global.begin(), global.end(), [](const auto& entry) {
            return entry.first.starts_with("current_source_");
        });
    if (has_current_source) {
        if (!global.count("current_source_species")) {
            throw std::runtime_error(
                "2D current_source_* controls require current_source_species");
        }
        CurrentRegulatedSource2DConfig source;
        source.species = as<std::string>(
            global, "current_source_species", "");
        source.control_mode = parse_current_source_control_mode(
            global, "current_source_control_mode",
            source.control_mode);
        source.monitor_boundary = parse_boundary_side_2d(
            global, "current_source_monitor_boundary",
            source.monitor_boundary);
        source.emission_boundary = parse_boundary_side_2d(
            global, "current_source_emission_boundary",
            source.emission_boundary);
        source.emission_inset = as<double>(
            global, "current_source_emission_inset",
            source.emission_inset);
        source.drift.x = as<double>(
            global, "current_source_drift_velocity_x", 0.0);
        source.drift.y = as<double>(
            global, "current_source_drift_velocity_y", 0.0);
        source.drift.z = as<double>(
            global, "current_source_drift_velocity_z", 0.0);
        if (global.count("current_source_temperature_ev")) {
            if (global.count("current_source_thermal_velocity")) {
                throw std::runtime_error(
                    "2D current source temperature_ev and thermal_velocity are mutually exclusive");
            }
            if (cfg.units.system != UnitSystem::SI) {
                throw std::runtime_error(
                    "2D current_source_temperature_ev requires units = si");
            }
            current_source_temperature_ev = as<double>(
                global, "current_source_temperature_ev", 0.0);
        } else {
            source.thermal_velocity = as<double>(
                global, "current_source_thermal_velocity", 0.0);
        }
        cfg.current_regulated_source = std::move(source);
    }
    const bool has_potential_reference = std::any_of(
        global.begin(), global.end(), [](const auto& entry) {
            return entry.first.starts_with("potential_reference_");
        });
    if (has_potential_reference) {
        if (!global.count("potential_reference_axis") ||
            !global.count("potential_reference_coordinate")) {
            throw std::runtime_error(
                "2D potential reference requires potential_reference_axis and potential_reference_coordinate");
        }
        PotentialReference2DConfig reference;
        reference.axis = parse_coordinate_axis(
            as<std::string>(
                global, "potential_reference_axis", ""));
        reference.coordinate = as<double>(
            global, "potential_reference_coordinate", 0.0);
        reference.target = as<double>(
            global, "potential_reference_target", 0.0);
        reference.correction = parse_potential_reference_correction(
            global, "potential_reference_correction",
            reference.correction);
        cfg.potential_reference = reference;
    }

    cfg.species.clear();
    for (const auto& block : blocks.species_blocks) {
        Species2DConfig s;
        s.name = as<std::string>(block, "name", s.name);
        s.charge = as<double>(block, "charge", s.charge);
        s.mass = as<double>(block, "mass", s.mass);
        s.particles = as<std::size_t>(block, "particles", s.particles);
        s.drift_velocity_x = as<double>(block, "drift_velocity_x", s.drift_velocity_x);
        s.drift_velocity_y = as<double>(block, "drift_velocity_y", s.drift_velocity_y);
        s.drift_velocity_z =
            as<double>(
                block, "drift_velocity_z",
                s.drift_velocity_z);
        if (block.count("temperature_ev")) {
            if (block.count("thermal_velocity") ||
                block.count("thermal_velocity_x") ||
                block.count("thermal_velocity_y") ||
                block.count("thermal_velocity_z")) {
                throw std::runtime_error(
                    "2D species '" + s.name +
                    "' temperature_ev and thermal_velocity controls are mutually exclusive");
            }
            if (cfg.units.system != UnitSystem::SI) {
                throw std::runtime_error(
                    "2D species '" + s.name +
                    "' temperature_ev requires units = si");
            }
            s.thermal_velocity =
                maxwellian_thermal_velocity_from_ev(
                    as<double>(block, "temperature_ev", 0.0),
                    s.mass);
        } else {
            s.thermal_velocity = as<double>(
                block, "thermal_velocity", s.thermal_velocity);
        }
        s.initialization.version = as<std::size_t>(
            block, "initialization_version", s.initialization.version);
        if (block.count("loading")) {
            s.initialization.loading = particle_loading_from_string(
                as<std::string>(block, "loading", "random"));
        }
        if (block.count("thermal_velocity_x")) {
            s.initialization.thermal_velocity_x =
                as<double>(block, "thermal_velocity_x", 0.0);
        }
        if (block.count("thermal_velocity_y")) {
            s.initialization.thermal_velocity_y =
                as<double>(block, "thermal_velocity_y", 0.0);
        }
        if (block.count("thermal_velocity_z")) {
            s.initialization.thermal_velocity_z =
                as<double>(block, "thermal_velocity_z", 0.0);
        }
        parse_density_profile(block, s.initialization, 2);
        s.init_x_min = as<double>(block, "init_x_min", s.init_x_min);
        s.init_x_max = as<double>(block, "init_x_max", s.init_x_max);
        s.init_y_min = as<double>(block, "init_y_min", s.init_y_min);
        s.init_y_max = as<double>(block, "init_y_max", s.init_y_max);
        require_species_scale_source(block, s.name, "2D");
        if (block.count("density")) {
            const double density = as<double>(block, "density", 1.0);
            validate_positive(density, "2D species '" + s.name + "' density");
            if (!block.count("weight")) {
                const double xmax = s.init_x_max < 0.0 ? cfg.length_x : s.init_x_max;
                const double ymax = s.init_y_max < 0.0 ? cfg.length_y : s.init_y_max;
                s.weight = density * (xmax - s.init_x_min) *
                    (ymax - s.init_y_min) *
                    cfg.out_of_plane_depth /
                    static_cast<double>(s.particles);
            }
        }
        if (block.count("weight")) s.weight = as<double>(block, "weight", s.weight);
        cfg.species.push_back(s);
    }
    if (cfg.species.empty()) cfg.species.push_back(Species2DConfig{});
    if (current_source_temperature_ev) {
        const auto& source = *cfg.current_regulated_source;
        const auto species = std::find_if(
            cfg.species.begin(), cfg.species.end(),
            [&](const Species2DConfig& candidate) {
                return candidate.name == source.species;
            });
        if (species == cfg.species.end()) {
            throw std::runtime_error(
                "2D current source temperature references unknown species '" +
                source.species + "'");
        }
        cfg.current_regulated_source->thermal_velocity =
            maxwellian_thermal_velocity_from_ev(
                *current_source_temperature_ev, species->mass);
    }
    for (const auto& block : blocks.source_blocks) {
        VolumetricPairSource2DConfig source;
        source.name = block.name;
        source.first_species = as<std::string>(
            block.values, "first_species", "");
        source.second_species = as<std::string>(
            block.values, "second_species", "");
        source.pairs_per_step = as<std::size_t>(
            block.values, "pairs_per_step", 0);
        if (block.values.count("represented_pair_rate")) {
            source.represented_pair_rate = as<double>(
                block.values, "represented_pair_rate", 0.0);
        }
        if (block.values.count("peak_volumetric_pair_rate")) {
            source.peak_volumetric_pair_rate = as<double>(
                block.values, "peak_volumetric_pair_rate", 0.0);
        }
        source.start_step = as<std::size_t>(
            block.values, "start_step", source.start_step);
        source.end_step = as<std::size_t>(
            block.values, "end_step", source.end_step);
        source.x_min = as<double>(
            block.values, "x_min", source.x_min);
        source.x_max = as<double>(
            block.values, "x_max", source.x_max);
        source.y_min = as<double>(
            block.values, "y_min", source.y_min);
        source.y_max = as<double>(
            block.values, "y_max", source.y_max);
        source.first_drift.x = as<double>(
            block.values, "first_drift_velocity_x", 0.0);
        source.first_drift.y = as<double>(
            block.values, "first_drift_velocity_y", 0.0);
        source.first_drift.z = as<double>(
            block.values, "first_drift_velocity_z", 0.0);
        source.second_drift.x = as<double>(
            block.values, "second_drift_velocity_x", 0.0);
        source.second_drift.y = as<double>(
            block.values, "second_drift_velocity_y", 0.0);
        source.second_drift.z = as<double>(
            block.values, "second_drift_velocity_z", 0.0);
        const auto species_mass =
            [&](const std::string& name) {
                const auto species = std::find_if(
                    cfg.species.begin(), cfg.species.end(),
                    [&](const Species2DConfig& candidate) {
                        return candidate.name == name;
                    });
                if (species == cfg.species.end()) {
                    throw std::runtime_error(
                        "2D source '" + source.name +
                        "' temperature references unknown species '" +
                        name + "'");
                }
                return species->mass;
            };
        if (block.values.count("first_temperature_ev")) {
            if (block.values.count("first_thermal_velocity")) {
                throw std::runtime_error(
                    "2D source '" + source.name +
                    "' first_temperature_ev and first_thermal_velocity are mutually exclusive");
            }
            if (cfg.units.system != UnitSystem::SI) {
                throw std::runtime_error(
                    "2D source '" + source.name +
                    "' first_temperature_ev requires units = si");
            }
            source.first_thermal_velocity =
                maxwellian_thermal_velocity_from_ev(
                    as<double>(
                        block.values,
                        "first_temperature_ev", 0.0),
                    species_mass(source.first_species));
        } else {
            source.first_thermal_velocity = as<double>(
                block.values, "first_thermal_velocity", 0.0);
        }
        if (block.values.count("second_temperature_ev")) {
            if (block.values.count("second_thermal_velocity")) {
                throw std::runtime_error(
                    "2D source '" + source.name +
                    "' second_temperature_ev and second_thermal_velocity are mutually exclusive");
            }
            if (cfg.units.system != UnitSystem::SI) {
                throw std::runtime_error(
                    "2D source '" + source.name +
                    "' second_temperature_ev requires units = si");
            }
            source.second_thermal_velocity =
                maxwellian_thermal_velocity_from_ev(
                    as<double>(
                        block.values,
                        "second_temperature_ev", 0.0),
                    species_mass(source.second_species));
        } else {
            source.second_thermal_velocity = as<double>(
                block.values, "second_thermal_velocity", 0.0);
        }
        parse_density_profile(
            block.values, source.spatial_profile, 2);
        cfg.sources.push_back(std::move(source));
    }
    if (cfg.resolved_diagnostics.enabled &&
        cfg.resolved_diagnostics.interval == 0) {
        cfg.resolved_diagnostics.interval = cfg.output_interval;
    }
    if (cfg.checkpoint_output && cfg.checkpoint_interval == 0) cfg.checkpoint_interval = cfg.output_interval;
    validate_config_2d(cfg);
    return cfg;
}

Simulation3DConfig load_config_3d(const std::string& path) {
    static const std::unordered_set<std::string> global_keys{
        "dimension", "config_version", "nx", "ny", "nz", "length_x", "length_y", "length_z", "dt", "steps",
        "mode", "steady_tolerance", "steady_window", "max_steps",
        "output_interval", "output_dir", "seed", "boundary", "vtk_output", "vtk_format",
        "particle_output", "particle_output_interval", "particle_output_stride", "particle_sample_count",
        "checkpoint_output", "checkpoint_interval", "checkpoint_path",
        "restart_path", "initial_state_path",
        "initial_state_signature",
        "runtime_backend", "runtime_threads", "magnetic_field_x", "magnetic_field_y", "magnetic_field_z",
        "magnetic_field_profile_file", "magnetic_field_profile_axis",
        "particle_boundary", "particle_boundary_left", "particle_boundary_right",
        "particle_boundary_bottom", "particle_boundary_top",
        "particle_boundary_back", "particle_boundary_front", "units",
        "relative_permittivity",
        "initialization_max_relative_charge_imbalance",
        "initialization_max_relative_current_imbalance",
        "initialization_max_relative_pair_imbalance",
        "initialization_charge_pairs"
    };
    static const std::unordered_set<std::string> species_keys{
        "name", "charge", "mass", "weight", "density", "particles", "drift_velocity_x",
        "drift_velocity_y", "drift_velocity_z", "thermal_velocity",
        "thermal_velocity_x", "thermal_velocity_y", "thermal_velocity_z",
        "initialization_version", "loading",
        "density_profile", "profile_center_x", "profile_center_y",
        "profile_center_z", "profile_scale_x", "profile_scale_y",
        "profile_scale_z", "profile_amplitude", "profile_phase",
        "profile_mode_x", "profile_mode_y", "profile_mode_z",
        "max_profile_sampling_attempts", "init_x_min", "init_x_max",
        "init_y_min", "init_y_max", "init_z_min", "init_z_max"
    };

    auto blocks = parse_config_blocks(
        path, global_keys, species_keys, nullptr, nullptr, nullptr, "3D");
    const auto& global = blocks.global;

    (void)parse_config_version(global, "3D");
    if (!global.count("dimension")) throw std::runtime_error("3D config loader requires dimension = 3");
    const auto dimension = as<std::size_t>(global, "dimension", 3);
    if (dimension != 3) throw std::runtime_error("3D config loader requires dimension = 3");

    Simulation3DConfig cfg;
    cfg.units = parse_units(global, cfg.units);
    cfg.nx = as<std::size_t>(global, "nx", cfg.nx);
    cfg.ny = as<std::size_t>(global, "ny", cfg.ny);
    cfg.nz = as<std::size_t>(global, "nz", cfg.nz);
    cfg.length_x = as<double>(global, "length_x", cfg.length_x);
    cfg.length_y = as<double>(global, "length_y", cfg.length_y);
    cfg.length_z = as<double>(global, "length_z", cfg.length_z);
    cfg.dt = as<double>(global, "dt", cfg.dt);
    cfg.steps = as<std::size_t>(global, "steps", cfg.steps);
    cfg.mode = parse_mode(global, cfg.mode);
    cfg.steady_tolerance = as<double>(global, "steady_tolerance", cfg.steady_tolerance);
    cfg.steady_window = as<std::size_t>(global, "steady_window", cfg.steady_window);
    cfg.max_steps = as<std::size_t>(global, "max_steps", cfg.max_steps);
    cfg.output_interval = as<std::size_t>(global, "output_interval", cfg.output_interval);
    cfg.output_dir = as<std::string>(global, "output_dir", cfg.output_dir.string());
    cfg.seed = as<unsigned>(global, "seed", cfg.seed);
    cfg.boundary = parse_boundary(global, cfg.boundary);
    cfg.vtk_output = parse_bool(global, "vtk_output", cfg.vtk_output);
    cfg.vtk_format = parse_vtk_output_format(global, cfg.vtk_format);
    cfg.particle_output = parse_bool(global, "particle_output", cfg.particle_output);
    cfg.particle_output_interval = as<std::size_t>(global, "particle_output_interval", cfg.particle_output_interval);
    cfg.particle_output_stride = as<std::size_t>(global, "particle_output_stride", cfg.particle_output_stride);
    cfg.particle_sample_count = as<std::size_t>(global, "particle_sample_count", cfg.particle_sample_count);
    cfg.checkpoint_output = parse_bool(global, "checkpoint_output", cfg.checkpoint_output);
    cfg.checkpoint_interval = as<std::size_t>(global, "checkpoint_interval", cfg.checkpoint_interval);
    cfg.checkpoint_path = as<std::string>(global, "checkpoint_path", cfg.checkpoint_path.string());
    cfg.restart_path = as<std::string>(global, "restart_path", cfg.restart_path.string());
    if (global.count("initial_state_path")) {
        cfg.initial_state_path = resolved_input_path(
            path, as<std::string>(
                      global, "initial_state_path", ""));
    }
    cfg.initial_state_signature = parse_optional_uint64(
        global, "initial_state_signature");
    cfg.runtime = parse_runtime_policy(global, cfg.runtime);
    cfg.initialization_acceptance =
        parse_initialization_acceptance(global);
    cfg.magnetic_field.x = as<double>(global, "magnetic_field_x", cfg.magnetic_field.x);
    cfg.magnetic_field.y = as<double>(global, "magnetic_field_y", cfg.magnetic_field.y);
    cfg.magnetic_field.z = as<double>(global, "magnetic_field_z", cfg.magnetic_field.z);
    cfg.magnetic_field_profile =
        parse_magnetic_field_profile(global, path);
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
        s.initialization.version = as<std::size_t>(
            block, "initialization_version", s.initialization.version);
        if (block.count("loading")) {
            s.initialization.loading = particle_loading_from_string(
                as<std::string>(block, "loading", "random"));
        }
        if (block.count("thermal_velocity_x")) {
            s.initialization.thermal_velocity_x =
                as<double>(block, "thermal_velocity_x", 0.0);
        }
        if (block.count("thermal_velocity_y")) {
            s.initialization.thermal_velocity_y =
                as<double>(block, "thermal_velocity_y", 0.0);
        }
        if (block.count("thermal_velocity_z")) {
            s.initialization.thermal_velocity_z =
                as<double>(block, "thermal_velocity_z", 0.0);
        }
        parse_density_profile(block, s.initialization, 3);
        s.init_x_min = as<double>(block, "init_x_min", s.init_x_min);
        s.init_x_max = as<double>(block, "init_x_max", s.init_x_max);
        s.init_y_min = as<double>(block, "init_y_min", s.init_y_min);
        s.init_y_max = as<double>(block, "init_y_max", s.init_y_max);
        s.init_z_min = as<double>(block, "init_z_min", s.init_z_min);
        s.init_z_max = as<double>(block, "init_z_max", s.init_z_max);
        require_species_scale_source(block, s.name, "3D");
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
