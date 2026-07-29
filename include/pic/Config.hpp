#pragma once
#include "pic/Initialization.hpp"
#include "pic/Runtime.hpp"
#include "pic/Types.hpp"
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace pic {
struct Simulation2DConfig;
struct Simulation3DConfig;

struct SpeciesConfig {
    std::string name{"electrons"};
    double charge{-1.0};
    double mass{1.0};
    double weight{1.0};
    std::size_t particles{1000};
    double density{1.0};
    double drift_velocity{0.0};
    double thermal_velocity{0.1};
    double init_x_min{0.0};
    double init_x_max{-1.0}; // negative means full domain
    ParticleInitializationConfig initialization{};
    double drift_velocity_y{0.0};
    double drift_velocity_z{0.0};
};

enum class CollisionModelKind { BGK, NullCollision };
enum class AngularScatteringKind {
    Isotropic,
    HenyeyGreenstein
};
enum class CollisionProcessKind {
    Elastic,
    Excitation,
    Ionization,
    Attachment,
    ChargeExchange
};

inline std::string to_string(CollisionModelKind model) {
    switch (model) {
        case CollisionModelKind::BGK: return "bgk";
        case CollisionModelKind::NullCollision: return "null_collision";
    }
    return "unknown";
}

inline std::string to_string(CollisionProcessKind process) {
    switch (process) {
        case CollisionProcessKind::Elastic: return "elastic";
        case CollisionProcessKind::Excitation: return "excitation";
        case CollisionProcessKind::Ionization: return "ionization";
        case CollisionProcessKind::Attachment: return "attachment";
        case CollisionProcessKind::ChargeExchange:
            return "charge_exchange";
    }
    return "unknown";
}

inline std::string to_string(AngularScatteringKind model) {
    switch (model) {
        case AngularScatteringKind::Isotropic:
            return "isotropic";
        case AngularScatteringKind::HenyeyGreenstein:
            return "henyey_greenstein";
    }
    return "unknown";
}

struct CollisionChannelConfig {
    std::string name{};
    CollisionProcessKind process{CollisionProcessKind::Elastic};
    std::filesystem::path cross_section_file{};
    double threshold_energy{0.0};
    double energy_scale{1.0};
    double cross_section_scale{1.0};
    std::string secondary_species{};
    std::string ion_species{};
    std::string attachment_species{};
    AngularScatteringKind angular_scattering{
        AngularScatteringKind::Isotropic};
    std::filesystem::path mean_cosine_file{};
    double mean_cosine_energy_scale{1.0};
};

struct CollisionConfig {
    bool enabled{false};
    CollisionModelKind model{CollisionModelKind::BGK};
    double frequency{0.0};
    double neutral_temperature_velocity{0.0};
    double neutral_density{0.0};
    std::string species{};
    double max_frequency{0.0};
    std::size_t max_candidates_per_particle{64};
    std::vector<CollisionChannelConfig> channels{};
    std::string gas_name{};
    double neutral_mass{0.0};
    double neutral_temperature{0.0};
    std::string data_provenance{};
    std::filesystem::path gas_data_file{};
    std::size_t gas_data_version{0};
    UnitSystem gas_data_units{UnitSystem::Normalized};
    std::string dataset_id{};
    std::string dataset_version{};
    std::string citation{};
    std::string retrieved{};
    std::string license{};
};

struct NamedCollisionConfig {
    std::string name{};
    CollisionConfig config{};
};

struct SinusoidalVoltageConfig {
    double amplitude{0.0};
    double frequency{0.0};
    double phase{0.0};
};

struct Config {
    UnitSystemConfig units{};
    std::size_t velocity_dimensions{1};
    std::size_t nx{128};
    double length{1.0};
    double dt{0.02};
    std::size_t steps{100};
    std::size_t output_interval{10};
    Boundary boundary{Boundary::Periodic};
    RunMode mode{RunMode::Transient};
    double phi_left{0.0};
    double phi_right{0.0};
    SinusoidalVoltageConfig phi_left_drive{};
    SinusoidalVoltageConfig phi_right_drive{};
    unsigned seed{12345};
    std::string output_dir{"output"};
    double steady_tolerance{1e-6};
    std::size_t steady_window{25};
    std::size_t max_steps{10000};
    std::size_t max_particles_per_species{10000000};
    CollisionConfig collisions{};
    std::vector<NamedCollisionConfig> collision_models{};
    bool checkpoint_output{false};
    std::size_t checkpoint_interval{0}; // zero inherits output_interval
    std::string checkpoint_path{};      // empty writes output_dir/checkpoint_<step>.apc
    std::string restart_path{};
    std::filesystem::path initial_state_path{};
    std::optional<std::uint64_t> initial_state_signature{};
    RuntimePolicy runtime{};
    InitializationAcceptanceConfig initialization_acceptance{};
    std::vector<SpeciesConfig> species{};
};

Config load_config(const std::string& path);
Simulation2DConfig load_config_2d(const std::string& path);
Simulation3DConfig load_config_3d(const std::string& path);
unsigned detect_config_dimension(const std::string& path);
std::string to_string(Boundary b);
std::string to_string(RunMode m);
}
