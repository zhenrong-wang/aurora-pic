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
    std::size_t timestep_multiplier{1};
};

enum class CollisionModelKind { BGK, NullCollision };
enum class AngularScatteringKind {
    Isotropic,
    Backward,
    HenyeyGreenstein
};
enum class CollisionEnergyFrame {
    Projectile,
    CenterOfMass
};
enum class CrossSectionInterpolationKind {
    Linear,
    LowerBin
};
enum class IonizationKinematicsKind {
    EqualEnergyIsotropic,
    OpalBeatyPeterson
};
enum class InelasticTransformKind {
    HeavyTarget,
    FiniteMassCenterOfMass
};
enum class CollisionProcessKind {
    Elastic,
    Excitation,
    Ionization,
    Attachment,
    ChargeExchange
};

enum class CollisionVelocitySampling1D {
    TimeCentered,
    LeapfrogHalfStep
};

inline std::string to_string(CollisionVelocitySampling1D sampling) {
    switch (sampling) {
        case CollisionVelocitySampling1D::TimeCentered:
            return "time_centered";
        case CollisionVelocitySampling1D::LeapfrogHalfStep:
            return "leapfrog_half_step";
    }
    return "unknown";
}

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
        case AngularScatteringKind::Backward:
            return "backward";
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
    CollisionEnergyFrame energy_frame{
        CollisionEnergyFrame::Projectile};
    IonizationKinematicsKind ionization_kinematics{
        IonizationKinematicsKind::EqualEnergyIsotropic};
    double ionization_ejected_energy_scale{0.0};
    CrossSectionInterpolationKind cross_section_interpolation{
        CrossSectionInterpolationKind::Linear};
    InelasticTransformKind inelastic_transform{
        InelasticTransformKind::HeavyTarget};
};

inline std::string to_string(CollisionEnergyFrame frame) {
    switch (frame) {
        case CollisionEnergyFrame::Projectile:
            return "projectile";
        case CollisionEnergyFrame::CenterOfMass:
            return "center_of_mass";
    }
    return "unknown";
}

inline std::string to_string(IonizationKinematicsKind model) {
    switch (model) {
        case IonizationKinematicsKind::EqualEnergyIsotropic:
            return "equal_energy_isotropic";
        case IonizationKinematicsKind::OpalBeatyPeterson:
            return "opal_beaty_peterson";
    }
    return "unknown";
}

inline std::string to_string(CrossSectionInterpolationKind interpolation) {
    switch (interpolation) {
        case CrossSectionInterpolationKind::Linear:
            return "linear";
        case CrossSectionInterpolationKind::LowerBin:
            return "lower_bin";
    }
    return "unknown";
}

inline std::string to_string(InelasticTransformKind transform) {
    switch (transform) {
        case InelasticTransformKind::HeavyTarget:
            return "heavy_target";
        case InelasticTransformKind::FiniteMassCenterOfMass:
            return "finite_mass_center_of_mass";
    }
    return "unknown";
}

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

enum class SpatialAverageSamplingOrder1D {
    PostCollision,
    PreCollision
};

inline std::string to_string(SpatialAverageSamplingOrder1D order) {
    switch (order) {
        case SpatialAverageSamplingOrder1D::PostCollision:
            return "post_collision";
        case SpatialAverageSamplingOrder1D::PreCollision:
            return "pre_collision";
    }
    return "unknown";
}

struct SpatialAverage1DConfig {
    bool enabled{false};
    bool reset_on_restart{false};
    std::size_t interval{1};
    std::size_t start_step{1};
    std::size_t end_step{0};
    double rf_frequency{0.0};
    std::size_t rf_cycles{0};
    std::size_t phase_bins{0};
    SpatialAverageSamplingOrder1D sampling_order{
        SpatialAverageSamplingOrder1D::PostCollision};
};

struct PhaseEedfRegion1DConfig {
    std::string name{};
    double x_min{0.0};
    double x_max{0.0};
};

struct PhaseEedf1DConfig {
    bool enabled{false};
    bool history_enabled{false};
    std::string species{};
    std::size_t energy_bins{0};
    double energy_max{0.0};
    double tail_threshold{0.0};
    std::vector<PhaseEedfRegion1DConfig> regions{};
};

struct PhaseSurfaceFlux1DConfig {
    bool enabled{false};
    bool reset_on_restart{false};
    std::string species{};
    std::vector<double> positions{};
    std::size_t energy_bins{0};
    double energy_max{0.0};
};

struct WallImpactSpectrum1DConfig {
    bool enabled{false};
    bool reset_on_restart{false};
    std::size_t energy_bins{0};
    double energy_max{0.0};
};

struct Config {
    UnitSystemConfig units{};
    std::size_t velocity_dimensions{1};
    std::size_t nx{128};
    double length{1.0};
    double dt{0.02};
    std::size_t steps{100};
    std::size_t output_interval{10};
    SpatialAverage1DConfig spatial_average{};
    PhaseEedf1DConfig phase_eedf{};
    PhaseSurfaceFlux1DConfig phase_surface_flux{};
    WallImpactSpectrum1DConfig wall_impact_spectrum{};
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
    CollisionVelocitySampling1D collision_velocity_sampling{
        CollisionVelocitySampling1D::TimeCentered};
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

void validate_spatial_average_1d(const Config& cfg);
Config load_config(const std::string& path);
Simulation2DConfig load_config_2d(const std::string& path);
Simulation3DConfig load_config_3d(const std::string& path);
unsigned detect_config_dimension(const std::string& path);
std::string to_string(Boundary b);
std::string to_string(RunMode m);
}
