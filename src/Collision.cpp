#include "pic/Collision.hpp"

#include <algorithm>
#include <bit>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace pic {
namespace {

constexpr double boltzmann_constant_si = 1.380649e-23;

void hash_bytes(std::uint64_t& hash, const void* data, std::size_t size) {
    constexpr std::uint64_t prime = 1099511628211ULL;
    const auto* bytes = static_cast<const unsigned char*>(data);
    for (std::size_t i = 0; i < size; ++i) {
        hash ^= bytes[i];
        hash *= prime;
    }
}

void hash_string(std::uint64_t& hash, const std::string& value) {
    hash_bytes(hash, value.data(), value.size());
    const unsigned char separator = 0xff;
    hash_bytes(hash, &separator, 1);
}

void hash_uint64(std::uint64_t& hash, std::uint64_t value) {
    unsigned char bytes[8]{};
    for (unsigned byte = 0; byte < 8; ++byte) {
        bytes[byte] =
            static_cast<unsigned char>(value >> (byte * 8));
    }
    hash_bytes(hash, bytes, sizeof(bytes));
}

void hash_double(std::uint64_t& hash, double value) {
    hash_uint64(hash, std::bit_cast<std::uint64_t>(value));
}

double open_unit_interval(std::mt19937_64& rng) {
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    double value = unit(rng);
    while (!(value > 0.0)) value = unit(rng);
    return value;
}

Vec3 isotropic_velocity(double speed, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    const double cosine = 2.0 * unit(rng) - 1.0;
    const double sine =
        std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
    const double azimuth =
        2.0 * std::acos(-1.0) * unit(rng);
    return {
        speed * sine * std::cos(azimuth),
        speed * sine * std::sin(azimuth),
        speed * cosine};
}

Vec3 angular_velocity(
    const Vec3& incoming,
    double speed,
    double mean_cosine,
    std::mt19937_64& rng) {
    const double incoming_speed = std::sqrt(
        incoming.x * incoming.x +
        incoming.y * incoming.y +
        incoming.z * incoming.z);
    if (!(incoming_speed > 0.0)) {
        return isotropic_velocity(speed, rng);
    }
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    const double sample = unit(rng);
    double cosine = 2.0 * sample - 1.0;
    if (mean_cosine != 0.0) {
        const double numerator =
            1.0 - mean_cosine * mean_cosine;
        const double denominator =
            1.0 - mean_cosine + 2.0 * mean_cosine * sample;
        const double ratio = numerator / denominator;
        cosine =
            (1.0 + mean_cosine * mean_cosine - ratio * ratio) /
            (2.0 * mean_cosine);
        cosine = std::clamp(cosine, -1.0, 1.0);
    }
    const double sine =
        std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
    const double azimuth =
        2.0 * std::acos(-1.0) * unit(rng);
    const Vec3 direction{
        incoming.x / incoming_speed,
        incoming.y / incoming_speed,
        incoming.z / incoming_speed};
    const Vec3 helper =
        std::abs(direction.z) < 0.9
            ? Vec3{0.0, 0.0, 1.0}
            : Vec3{0.0, 1.0, 0.0};
    Vec3 tangent{
        helper.y * direction.z - helper.z * direction.y,
        helper.z * direction.x - helper.x * direction.z,
        helper.x * direction.y - helper.y * direction.x};
    const double tangent_norm = std::sqrt(
        tangent.x * tangent.x +
        tangent.y * tangent.y +
        tangent.z * tangent.z);
    tangent.x /= tangent_norm;
    tangent.y /= tangent_norm;
    tangent.z /= tangent_norm;
    const Vec3 bitangent{
        direction.y * tangent.z - direction.z * tangent.y,
        direction.z * tangent.x - direction.x * tangent.z,
        direction.x * tangent.y - direction.y * tangent.x};
    const double tangent_factor = sine * std::cos(azimuth);
    const double bitangent_factor = sine * std::sin(azimuth);
    return {
        speed * (
            cosine * direction.x +
            tangent_factor * tangent.x +
            bitangent_factor * bitangent.x),
        speed * (
            cosine * direction.y +
            tangent_factor * tangent.y +
            bitangent_factor * bitangent.y),
        speed * (
            cosine * direction.z +
            tangent_factor * tangent.z +
            bitangent_factor * bitangent.z)};
}

double speed(const Vec3& value) {
    return std::sqrt(
        value.x * value.x +
        value.y * value.y +
        value.z * value.z);
}

Vec3 subtract(const Vec3& first, const Vec3& second) {
    return {
        first.x - second.x,
        first.y - second.y,
        first.z - second.z};
}

Vec3 add(const Vec3& first, const Vec3& second) {
    return {
        first.x + second.x,
        first.y + second.y,
        first.z + second.z};
}

struct IonizationVelocityPair {
    Vec3 primary{};
    Vec3 secondary{};
};

IonizationVelocityPair opal_ionization_velocities(
    const Vec3& incoming, double available_energy, double particle_mass,
    double ejected_energy_scale, std::mt19937_64& rng) {
    if (!(available_energy > 0.0)) return {};
    const double incoming_speed = speed(incoming);
    if (!(incoming_speed > 0.0)) return {};
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    const double secondary_energy = ejected_energy_scale * std::tan(
        unit(rng) * std::atan(
            available_energy / (2.0 * ejected_energy_scale)));
    const double primary_energy =
        std::max(0.0, available_energy - secondary_energy);
    const double primary_speed =
        std::sqrt(2.0 * primary_energy / particle_mass);
    const double secondary_speed =
        std::sqrt(2.0 * secondary_energy / particle_mass);
    const double primary_cosine =
        std::sqrt(primary_energy / available_energy);
    const double secondary_cosine =
        std::sqrt(secondary_energy / available_energy);
    const double primary_sine = secondary_cosine;
    const double secondary_sine = primary_cosine;
    const double azimuth = 2.0 * std::acos(-1.0) * unit(rng);
    const Vec3 direction{
        incoming.x / incoming_speed,
        incoming.y / incoming_speed,
        incoming.z / incoming_speed};
    const Vec3 helper =
        std::abs(direction.z) < 0.9
            ? Vec3{0.0, 0.0, 1.0}
            : Vec3{0.0, 1.0, 0.0};
    Vec3 tangent{
        helper.y * direction.z - helper.z * direction.y,
        helper.z * direction.x - helper.x * direction.z,
        helper.x * direction.y - helper.y * direction.x};
    const double tangent_norm = speed(tangent);
    tangent.x /= tangent_norm;
    tangent.y /= tangent_norm;
    tangent.z /= tangent_norm;
    const Vec3 bitangent{
        direction.y * tangent.z - direction.z * tangent.y,
        direction.z * tangent.x - direction.x * tangent.z,
        direction.x * tangent.y - direction.y * tangent.x};
    const Vec3 transverse{
        std::cos(azimuth) * tangent.x +
            std::sin(azimuth) * bitangent.x,
        std::cos(azimuth) * tangent.y +
            std::sin(azimuth) * bitangent.y,
        std::cos(azimuth) * tangent.z +
            std::sin(azimuth) * bitangent.z};
    return {
        {
            primary_speed *
                (primary_cosine * direction.x +
                 primary_sine * transverse.x),
            primary_speed *
                (primary_cosine * direction.y +
                 primary_sine * transverse.y),
            primary_speed *
                (primary_cosine * direction.z +
                 primary_sine * transverse.z),
        },
        {
            secondary_speed *
                (secondary_cosine * direction.x -
                 secondary_sine * transverse.x),
            secondary_speed *
                (secondary_cosine * direction.y -
                 secondary_sine * transverse.y),
            secondary_speed *
                (secondary_cosine * direction.z -
                 secondary_sine * transverse.z),
        }};
}

} // namespace

CrossSectionTable::CrossSectionTable(
    const std::filesystem::path& path,
    double energy_scale,
    double cross_section_scale,
    CrossSectionInterpolationKind interpolation)
    : interpolation_(interpolation) {
    if (!std::isfinite(energy_scale) || !(energy_scale > 0.0)) {
        throw std::invalid_argument(
            "collision cross-section energy_scale must be positive and finite");
    }
    if (!std::isfinite(cross_section_scale) ||
        !(cross_section_scale > 0.0)) {
        throw std::invalid_argument(
            "collision cross_section_scale must be positive and finite");
    }
    if (interpolation != CrossSectionInterpolationKind::Linear &&
        interpolation != CrossSectionInterpolationKind::LowerBin) {
        throw std::invalid_argument(
            "unsupported cross-section interpolation mode");
    }
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error(
            "cannot open collision cross-section table: " + path.string());
    }
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        const auto comment = line.find_first_of("#;");
        if (comment != std::string::npos) line.resize(comment);
        std::istringstream row(line);
        double energy = 0.0;
        double cross_section = 0.0;
        if (!(row >> energy)) continue;
        if (!(row >> cross_section)) {
            throw std::runtime_error(
                "collision table " + path.string() + ":" +
                std::to_string(line_number) +
                " must contain energy and cross section");
        }
        std::string trailing;
        if (row >> trailing) {
            throw std::runtime_error(
                "collision table " + path.string() + ":" +
                std::to_string(line_number) +
                " has trailing columns");
        }
        energy *= energy_scale;
        cross_section *= cross_section_scale;
        if (!std::isfinite(energy) || energy < 0.0 ||
            !std::isfinite(cross_section) || cross_section < 0.0) {
            throw std::runtime_error(
                "collision table values must be finite and non-negative");
        }
        if (!energies_.empty() && !(energy > energies_.back())) {
            throw std::runtime_error(
                "collision table energies must be strictly increasing");
        }
        energies_.push_back(energy);
        cross_sections_.push_back(cross_section);
    }
    if (energies_.size() < 2) {
        throw std::runtime_error(
            "collision cross-section table requires at least two rows: " +
            path.string());
    }
    const double candidate_spacing =
        (energies_.back() - energies_.front()) /
        static_cast<double>(energies_.size() - 1);
    bool uniformly_spaced =
        std::isfinite(candidate_spacing) && candidate_spacing > 0.0;
    for (std::size_t index = 1;
         uniformly_spaced && index + 1 < energies_.size(); ++index) {
        const double expected = energies_.front() +
            static_cast<double>(index) * candidate_spacing;
        const double tolerance =
            64.0 * std::numeric_limits<double>::epsilon() *
            std::max({std::abs(expected), candidate_spacing,
                      std::numeric_limits<double>::min()});
        uniformly_spaced =
            std::abs(energies_[index] - expected) <= tolerance;
    }
    if (uniformly_spaced) {
        uniform_energy_spacing_ = candidate_spacing;
    }
    maximum_tree_leaf_count_ = 1;
    while (maximum_tree_leaf_count_ < cross_sections_.size()) {
        maximum_tree_leaf_count_ *= 2;
    }
    maximum_tree_.assign(2 * maximum_tree_leaf_count_, 0.0);
    std::copy(
        cross_sections_.begin(), cross_sections_.end(),
        maximum_tree_.begin() +
            static_cast<std::ptrdiff_t>(maximum_tree_leaf_count_));
    for (std::size_t node = maximum_tree_leaf_count_ - 1;
         node > 0; --node) {
        maximum_tree_[node] = std::max(
            maximum_tree_[2 * node],
            maximum_tree_[2 * node + 1]);
    }
}

std::size_t CrossSectionTable::lower_bound_index(double energy) const {
    if (uniform_energy_spacing_ == 0.0) {
        return static_cast<std::size_t>(
            std::lower_bound(energies_.begin(), energies_.end(), energy) -
            energies_.begin());
    }
    if (energy <= energies_.front()) return 0;
    if (energy > energies_.back()) return energies_.size();
    const double coordinate =
        (energy - energies_.front()) / uniform_energy_spacing_;
    std::size_t index = static_cast<std::size_t>(std::floor(coordinate));
    while (index > 0 && energies_[index - 1] >= energy) --index;
    while (index < energies_.size() && energies_[index] < energy) ++index;
    return index;
}

std::size_t CrossSectionTable::upper_bound_index(double energy) const {
    if (uniform_energy_spacing_ == 0.0) {
        return static_cast<std::size_t>(
            std::upper_bound(energies_.begin(), energies_.end(), energy) -
            energies_.begin());
    }
    if (energy < energies_.front()) return 0;
    if (energy >= energies_.back()) return energies_.size();
    const double coordinate =
        (energy - energies_.front()) / uniform_energy_spacing_;
    std::size_t index =
        static_cast<std::size_t>(std::floor(coordinate)) + 1;
    while (index > 0 && energies_[index - 1] > energy) --index;
    while (index < energies_.size() && energies_[index] <= energy) ++index;
    return index;
}

double CrossSectionTable::evaluate(double energy) const {
    if (!std::isfinite(energy) || energy < 0.0) {
        throw std::invalid_argument(
            "collision lookup energy must be finite and non-negative");
    }
    if (energy <= energies_.front()) return cross_sections_.front();
    if (energy >= energies_.back()) return cross_sections_.back();
    const std::size_t high = upper_bound_index(energy);
    const std::size_t low = high - 1;
    if (interpolation_ == CrossSectionInterpolationKind::LowerBin) {
        return cross_sections_[low];
    }
    const double fraction =
        (energy - energies_[low]) / (energies_[high] - energies_[low]);
    return cross_sections_[low] +
           fraction * (cross_sections_[high] - cross_sections_[low]);
}

double CrossSectionTable::maximum_between(
    double minimum_energy, double maximum_energy) const {
    if (!std::isfinite(minimum_energy) || minimum_energy < 0.0 ||
        !std::isfinite(maximum_energy) ||
        maximum_energy < minimum_energy) {
        throw std::invalid_argument(
            "collision maximum lookup requires finite non-negative "
            "ordered energies");
    }
    double result = std::max(
        evaluate(minimum_energy), evaluate(maximum_energy));
    std::size_t first = lower_bound_index(minimum_energy);
    std::size_t last = upper_bound_index(maximum_energy);
    first += maximum_tree_leaf_count_;
    last += maximum_tree_leaf_count_;
    while (first < last) {
        if ((first & 1U) != 0U) {
            result = std::max(result, maximum_tree_[first]);
            ++first;
        }
        if ((last & 1U) != 0U) {
            --last;
            result = std::max(result, maximum_tree_[last]);
        }
        first /= 2;
        last /= 2;
    }
    return result;
}

MeanCosineTable::MeanCosineTable(
    const std::filesystem::path& path,
    double energy_scale) {
    if (!std::isfinite(energy_scale) || !(energy_scale > 0.0)) {
        throw std::invalid_argument(
            "mean-cosine energy_scale must be positive and finite");
    }
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error(
            "cannot open mean-cosine table: " + path.string());
    }
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        const auto comment = line.find_first_of("#;");
        if (comment != std::string::npos) line.resize(comment);
        std::istringstream row(line);
        double energy = 0.0;
        double mean_cosine = 0.0;
        if (!(row >> energy)) continue;
        if (!(row >> mean_cosine)) {
            throw std::runtime_error(
                "mean-cosine table " + path.string() + ":" +
                std::to_string(line_number) +
                " must contain energy and mean cosine");
        }
        std::string trailing;
        if (row >> trailing) {
            throw std::runtime_error(
                "mean-cosine table " + path.string() + ":" +
                std::to_string(line_number) +
                " has trailing columns");
        }
        energy *= energy_scale;
        if (!std::isfinite(energy) || energy < 0.0 ||
            !std::isfinite(mean_cosine) ||
            !(std::abs(mean_cosine) < 1.0)) {
            throw std::runtime_error(
                "mean-cosine table requires finite non-negative energies "
                "and mean cosines strictly between -1 and 1");
        }
        if (!energies_.empty() && !(energy > energies_.back())) {
            throw std::runtime_error(
                "mean-cosine table energies must be strictly increasing");
        }
        energies_.push_back(energy);
        mean_cosines_.push_back(mean_cosine);
    }
    if (energies_.size() < 2) {
        throw std::runtime_error(
            "mean-cosine table requires at least two rows: " +
            path.string());
    }
}

double MeanCosineTable::evaluate(double energy) const {
    if (!std::isfinite(energy) || energy < 0.0) {
        throw std::invalid_argument(
            "mean-cosine lookup energy must be finite and non-negative");
    }
    if (energy <= energies_.front()) return mean_cosines_.front();
    if (energy >= energies_.back()) return mean_cosines_.back();
    const auto upper =
        std::upper_bound(energies_.begin(), energies_.end(), energy);
    const std::size_t high =
        static_cast<std::size_t>(upper - energies_.begin());
    const std::size_t low = high - 1;
    const double fraction =
        (energy - energies_[low]) /
        (energies_[high] - energies_[low]);
    return mean_cosines_[low] +
           fraction * (mean_cosines_[high] - mean_cosines_[low]);
}

NullCollisionModel::NullCollisionModel(
    CollisionConfig config,
    double particle_mass)
    : config_(std::move(config)),
      particle_mass_(particle_mass) {
    if (!std::isfinite(particle_mass_) || !(particle_mass_ > 0.0)) {
        throw std::invalid_argument(
            "MCC particle mass must be positive and finite");
    }
    if (!std::isfinite(config_.neutral_density) ||
        !(config_.neutral_density > 0.0)) {
        throw std::invalid_argument(
            "MCC neutral_density must be positive and finite");
    }
    if (!std::isfinite(config_.max_frequency) ||
        !(config_.max_frequency > 0.0)) {
        throw std::invalid_argument(
            "MCC max_frequency must be positive and finite");
    }
    if (config_.max_candidates_per_particle == 0) {
        throw std::invalid_argument(
            "MCC max_candidates_per_particle must be positive");
    }
    if (config_.channels.empty()) {
        throw std::invalid_argument(
            "MCC requires at least one collision channel");
    }
    if (!std::isfinite(config_.neutral_mass) ||
        config_.neutral_mass < 0.0 ||
        !std::isfinite(config_.neutral_temperature) ||
        config_.neutral_temperature < 0.0) {
        throw std::invalid_argument(
            "MCC neutral mass and temperature must be finite and non-negative");
    }
    if (config_.gas_data_units == UnitSystem::SI &&
        config_.neutral_temperature > 0.0) {
        if (!(config_.neutral_mass > 0.0)) {
            throw std::invalid_argument(
                "thermal-neutral MCC requires positive neutral_mass");
        }
        neutral_velocity_stddev_ = std::sqrt(
            boltzmann_constant_si * config_.neutral_temperature /
            config_.neutral_mass);
        if (!std::isfinite(neutral_velocity_stddev_) ||
            !(neutral_velocity_stddev_ > 0.0)) {
            throw std::invalid_argument(
                "thermal-neutral velocity scale is not finite and positive");
        }
    }

    constexpr std::uint64_t fnv_offset = 1469598103934665603ULL;
    signature_ = fnv_offset;
    hash_string(signature_, to_string(config_.model));
    hash_string(signature_, config_.species);
    hash_double(signature_, config_.neutral_density);
    hash_double(signature_, config_.max_frequency);
    hash_double(signature_, particle_mass_);
    // Preserve the historical 1D MCC signature when the optional gas
    // metadata is absent. Imported 2D MCC requires these fields and therefore
    // fingerprints them without invalidating existing 1D v3 checkpoints.
    if (!config_.gas_name.empty() ||
        config_.neutral_mass != 0.0 ||
        config_.neutral_temperature != 0.0 ||
        !config_.data_provenance.empty()) {
        hash_string(signature_, config_.gas_name);
        hash_double(signature_, config_.neutral_mass);
        hash_double(signature_, config_.neutral_temperature);
        hash_string(signature_, config_.data_provenance);
        if (neutral_velocity_stddev_ > 0.0) {
            hash_string(
                signature_,
                "bounded_maxwellian_neutral_kinematics_v1");
            hash_double(signature_, neutral_velocity_stddev_);
            hash_double(
                signature_, neutral_speed_limit_sigma_);
        } else if (config_.neutral_mass > 0.0) {
            hash_string(
                signature_,
                "stationary_finite_mass_neutral_kinematics_v2");
        }
    }
    if (!config_.dataset_id.empty() ||
        !config_.dataset_version.empty() ||
        !config_.citation.empty() ||
        !config_.retrieved.empty() ||
        !config_.license.empty()) {
        hash_string(signature_, config_.dataset_id);
        hash_string(signature_, config_.dataset_version);
        hash_string(signature_, config_.citation);
        hash_string(signature_, config_.retrieved);
        hash_string(signature_, config_.license);
    }
    if (!config_.gas_data_file.empty()) {
        hash_uint64(
            signature_,
            static_cast<std::uint64_t>(
                config_.gas_data_version));
        hash_string(
            signature_, to_string(config_.gas_data_units));
    }
    const std::uint64_t candidate_limit =
        config_.max_candidates_per_particle;
    hash_uint64(signature_, candidate_limit);
    for (const auto& channel_config : config_.channels) {
        if (channel_config.name.empty()) {
            throw std::invalid_argument(
                "MCC collision channel name must not be empty");
        }
        if (!std::all_of(
                channel_config.name.begin(),
                channel_config.name.end(),
                [](unsigned char character) {
                    return std::isalnum(character) ||
                           character == '_' || character == '-' ||
                           character == '.';
                })) {
            throw std::invalid_argument(
                "MCC collision channel names may contain only letters, "
                "digits, '_', '-', and '.'");
        }
        if (!std::isfinite(channel_config.threshold_energy) ||
            channel_config.threshold_energy < 0.0) {
            throw std::invalid_argument(
                "MCC threshold_energy must be finite and non-negative");
        }
        if (channel_config.inelastic_transform !=
                InelasticTransformKind::HeavyTarget &&
            channel_config.inelastic_transform !=
                InelasticTransformKind::FiniteMassCenterOfMass) {
            throw std::invalid_argument(
                "unsupported MCC inelastic transform");
        }
        if ((channel_config.process == CollisionProcessKind::Elastic ||
             channel_config.process == CollisionProcessKind::Attachment ||
             channel_config.process ==
                 CollisionProcessKind::ChargeExchange) &&
            channel_config.threshold_energy != 0.0) {
            throw std::invalid_argument(
                "elastic, attachment, and charge-exchange MCC channel "
                "threshold_energy must be zero");
        }
        if ((channel_config.process == CollisionProcessKind::Excitation ||
             channel_config.process == CollisionProcessKind::Ionization) &&
            !(channel_config.threshold_energy > 0.0)) {
            throw std::invalid_argument(
                "inelastic MCC channel threshold_energy must be positive");
        }
        if (channel_config.angular_scattering !=
                AngularScatteringKind::Isotropic &&
            channel_config.process != CollisionProcessKind::Elastic) {
            throw std::invalid_argument(
                "anisotropic MCC scattering is supported only for "
                "elastic channels");
        }
        if (channel_config.angular_scattering ==
            AngularScatteringKind::HenyeyGreenstein) {
            if (channel_config.mean_cosine_file.empty()) {
                throw std::invalid_argument(
                    "Henyey-Greenstein scattering requires "
                    "mean_cosine_file");
            }
            if (!std::isfinite(
                    channel_config.mean_cosine_energy_scale) ||
                !(channel_config.mean_cosine_energy_scale > 0.0)) {
                throw std::invalid_argument(
                    "mean_cosine_energy_scale must be positive and finite");
            }
        } else if (!channel_config.mean_cosine_file.empty() ||
                   channel_config.mean_cosine_energy_scale != 1.0) {
            throw std::invalid_argument(
                "mean-cosine data requires angular_model = "
                "henyey_greenstein");
        }
        if (channel_config.energy_frame ==
                CollisionEnergyFrame::CenterOfMass &&
            !(config_.neutral_mass > 0.0)) {
            throw std::invalid_argument(
                "center-of-mass collision energy requires positive "
                "neutral_mass");
        }
        if (channel_config.energy_frame ==
                CollisionEnergyFrame::CenterOfMass &&
            channel_config.process != CollisionProcessKind::Elastic &&
            channel_config.process !=
                CollisionProcessKind::ChargeExchange) {
            throw std::invalid_argument(
                "center-of-mass collision energy is supported only for "
                "elastic and charge-exchange channels");
        }
        if (channel_config.process == CollisionProcessKind::Ionization &&
            (channel_config.secondary_species.empty() ||
             channel_config.ion_species.empty())) {
            throw std::invalid_argument(
                "ionization MCC channel requires secondary and ion species");
        }
        if (channel_config.process == CollisionProcessKind::Ionization) {
            if (channel_config.ionization_kinematics ==
                    IonizationKinematicsKind::OpalBeatyPeterson) {
                if (!std::isfinite(
                        channel_config.ionization_ejected_energy_scale) ||
                    !(channel_config.ionization_ejected_energy_scale > 0.0)) {
                    throw std::invalid_argument(
                        "Opal-Beaty-Peterson ionization requires positive "
                        "finite ionization_ejected_energy_scale");
                }
            } else if (
                channel_config.ionization_ejected_energy_scale != 0.0) {
                throw std::invalid_argument(
                    "ionization_ejected_energy_scale requires "
                    "Opal-Beaty-Peterson ionization kinematics");
            }
        } else if (
            channel_config.ionization_kinematics !=
                IonizationKinematicsKind::EqualEnergyIsotropic ||
            channel_config.ionization_ejected_energy_scale != 0.0) {
            throw std::invalid_argument(
                "ionization kinematics is valid only for ionization channels");
        }
        const bool inelastic =
            channel_config.process == CollisionProcessKind::Excitation ||
            channel_config.process == CollisionProcessKind::Ionization;
        if (channel_config.inelastic_transform !=
                InelasticTransformKind::HeavyTarget &&
            (!inelastic || !(config_.neutral_mass > 0.0))) {
            throw std::invalid_argument(
                "finite-mass inelastic transform requires an "
                "excitation or ionization channel and positive "
                "neutral_mass");
        }
        if (channel_config.process == CollisionProcessKind::Attachment &&
            channel_config.attachment_species.empty()) {
            throw std::invalid_argument(
                "attachment MCC channel requires an attachment species");
        }
        if (channel_config.process ==
            CollisionProcessKind::ChargeExchange) {
            const double mass_tolerance =
                64.0 * std::numeric_limits<double>::epsilon() *
                std::max(particle_mass_, config_.neutral_mass);
            if (!(config_.neutral_mass > 0.0) ||
                std::abs(particle_mass_ - config_.neutral_mass) >
                    mass_tolerance) {
                throw std::invalid_argument(
                    "resonant charge exchange requires projectile mass "
                    "equal to neutral_mass");
            }
        }
        if (std::find(channel_names_.begin(), channel_names_.end(),
                      channel_config.name) != channel_names_.end()) {
            throw std::invalid_argument(
                "duplicate MCC collision channel name: " +
                channel_config.name);
        }
        channels_.emplace_back(channel_config);
        if (channels_.back().mean_cosine.has_value() &&
            (channels_.back().mean_cosine->energies().front() >
                 channels_.back().table.energies().front() ||
             channels_.back().mean_cosine->energies().back() <
                 channels_.back().table.energies().back())) {
            throw std::invalid_argument(
                "mean-cosine energy range must cover the "
                "cross-section table for channel " +
                channel_config.name);
        }
        channel_names_.push_back(channel_config.name);
        hash_string(signature_, channel_config.name);
        hash_string(signature_, to_string(channel_config.process));
        hash_string(
            signature_,
            to_string(channel_config.angular_scattering));
        if (channel_config.cross_section_interpolation !=
            CrossSectionInterpolationKind::Linear) {
            hash_string(
                signature_,
                to_string(channel_config.cross_section_interpolation));
        }
        if (channel_config.energy_frame !=
            CollisionEnergyFrame::Projectile) {
            hash_string(
                signature_, to_string(channel_config.energy_frame));
        }
        if (channel_config.inelastic_transform !=
            InelasticTransformKind::HeavyTarget) {
            hash_string(
                signature_,
                to_string(channel_config.inelastic_transform));
        }
        hash_double(signature_, channel_config.threshold_energy);
        if (channel_config.process ==
            CollisionProcessKind::Ionization) {
            hash_string(signature_, channel_config.secondary_species);
            hash_string(signature_, channel_config.ion_species);
            if (channel_config.ionization_kinematics !=
                IonizationKinematicsKind::EqualEnergyIsotropic) {
                hash_string(
                    signature_,
                    to_string(channel_config.ionization_kinematics));
                hash_double(
                    signature_,
                    channel_config.ionization_ejected_energy_scale);
            }
        } else if (channel_config.process ==
                   CollisionProcessKind::Attachment) {
            hash_string(signature_, channel_config.attachment_species);
        }
        for (std::size_t i = 0;
             i < channels_.back().table.energies().size(); ++i) {
            hash_double(signature_, channels_.back().table.energies()[i]);
            hash_double(
                signature_, channels_.back().table.cross_sections()[i]);
        }
        if (channels_.back().mean_cosine.has_value()) {
            for (std::size_t i = 0;
                 i < channels_.back().mean_cosine->energies().size();
                 ++i) {
                hash_double(
                    signature_,
                    channels_.back().mean_cosine->energies()[i]);
                hash_double(
                    signature_,
                    channels_.back().mean_cosine->mean_cosines()[i]);
            }
        }
    }
}

double NullCollisionModel::collision_energy(
    const Channel& channel, double relative_speed) const {
    double energy_mass = particle_mass_;
    if (channel.config.energy_frame ==
        CollisionEnergyFrame::CenterOfMass) {
        energy_mass =
            particle_mass_ * config_.neutral_mass /
            (particle_mass_ + config_.neutral_mass);
    }
    const double energy =
        0.5 * energy_mass * relative_speed * relative_speed;
    if (!std::isfinite(energy)) {
        throw std::overflow_error("MCC collision energy overflow");
    }
    return energy;
}

double NullCollisionModel::rates_for_speed(
    double speed, std::vector<double>& rates) const {
    if (!std::isfinite(speed) || speed < 0.0) {
        throw std::invalid_argument(
            "MCC particle speed must be finite and non-negative");
    }
    rates.resize(channels_.size());
    double total = 0.0;
    for (std::size_t index = 0; index < channels_.size(); ++index) {
        const auto& channel = channels_[index];
        const double energy = collision_energy(channel, speed);
        double rate = 0.0;
        if (energy >= channel.config.threshold_energy) {
            rate = config_.neutral_density *
                   channel.table.evaluate(energy) * speed;
        }
        if (!std::isfinite(rate) || rate < 0.0) {
            throw std::overflow_error("MCC collision rate overflow");
        }
        rates[index] = rate;
        total += rate;
    }
    const double tolerance =
        64.0 * std::numeric_limits<double>::epsilon() *
        std::max(config_.max_frequency, total);
    if (total > config_.max_frequency + tolerance) {
        throw std::runtime_error(
            "MCC total collision frequency exceeds configured max_frequency");
    }
    return total;
}

void NullCollisionModel::validate_frequency_bound(
    double projectile_speed,
    std::vector<double>& rate_scratch) const {
    if (neutral_velocity_stddev_ == 0.0) {
        (void)rates_for_speed(projectile_speed, rate_scratch);
        return;
    }
    const double thermal_speed_limit =
        neutral_speed_limit_sigma_ * neutral_velocity_stddev_;
    const double maximum_relative_speed =
        projectile_speed + thermal_speed_limit;
    if (!std::isfinite(maximum_relative_speed)) {
        throw std::overflow_error(
            "thermal-neutral relative speed overflow");
    }
    double frequency_bound = 0.0;
    for (const auto& channel : channels_) {
        const double maximum_energy =
            collision_energy(channel, maximum_relative_speed);
        if (maximum_energy < channel.config.threshold_energy) {
            continue;
        }
        const double minimum_relative_speed =
            std::max(
                0.0,
                projectile_speed - thermal_speed_limit);
        const double minimum_energy =
            std::max(
                channel.config.threshold_energy,
                collision_energy(channel, minimum_relative_speed));
        const double maximum_cross_section =
            channel.table.maximum_between(
                minimum_energy, maximum_energy);
        frequency_bound +=
            config_.neutral_density *
            maximum_cross_section *
            maximum_relative_speed;
    }
    if (!std::isfinite(frequency_bound)) {
        throw std::overflow_error(
            "thermal-neutral collision-frequency bound overflow");
    }
    const double tolerance =
        64.0 * std::numeric_limits<double>::epsilon() *
        std::max(config_.max_frequency, frequency_bound);
    if (frequency_bound > config_.max_frequency + tolerance) {
        throw std::runtime_error(
            "MCC thermal-neutral collision-frequency bound exceeds "
            "configured max_frequency");
    }
}

double NullCollisionModel::sample_neutral_velocity(
    std::mt19937_64& rng) const {
    if (neutral_velocity_stddev_ == 0.0) return 0.0;
    std::normal_distribution<double> normal(
        0.0, neutral_velocity_stddev_);
    const double limit =
        neutral_speed_limit_sigma_ * neutral_velocity_stddev_;
    for (std::size_t attempt = 0; attempt < 1024; ++attempt) {
        const double value = normal(rng);
        if (std::abs(value) <= limit) return value;
    }
    throw std::runtime_error(
        "thermal-neutral scalar sampler exceeded its rejection limit");
}

Vec3 NullCollisionModel::sample_neutral_velocity_3v(
    std::mt19937_64& rng) const {
    if (neutral_velocity_stddev_ == 0.0) return {};
    std::normal_distribution<double> normal(
        0.0, neutral_velocity_stddev_);
    const double limit =
        neutral_speed_limit_sigma_ * neutral_velocity_stddev_;
    for (std::size_t attempt = 0; attempt < 1024; ++attempt) {
        const Vec3 value{normal(rng), normal(rng), normal(rng)};
        if (speed(value) <= limit) return value;
    }
    throw std::runtime_error(
        "thermal-neutral 3V sampler exceeded its rejection limit");
}

void NullCollisionModel::apply_channel(
    std::size_t channel_index,
    double& velocity,
    double neutral_velocity,
    std::mt19937_64& rng) const {
    const auto& channel = channels_.at(channel_index);
    if (channel.config.process == CollisionProcessKind::Attachment) {
        velocity = 0.0;
        return;
    }
    if (channel.config.process ==
        CollisionProcessKind::ChargeExchange) {
        velocity = neutral_velocity;
        return;
    }
    const double relative_velocity = velocity - neutral_velocity;
    double energy =
        0.5 * particle_mass_ *
        relative_velocity * relative_velocity;
    if (channel.config.process == CollisionProcessKind::Excitation ||
        channel.config.process == CollisionProcessKind::Ionization) {
        energy = std::max(
            0.0, energy - channel.config.threshold_energy);
    }
    if (channel.config.process == CollisionProcessKind::Ionization) {
        energy *= 0.5;
    }
    const double speed = std::sqrt(2.0 * energy / particle_mass_);
    std::uniform_int_distribution<int> direction(0, 1);
    const double scattered =
        direction(rng) == 0 ? -speed : speed;
    if (channel.config.inelastic_transform ==
        InelasticTransformKind::FiniteMassCenterOfMass) {
        const double total_mass =
            particle_mass_ + config_.neutral_mass;
        velocity =
            (particle_mass_ * velocity +
             config_.neutral_mass * neutral_velocity) /
                total_mass +
            (config_.neutral_mass / total_mass) * scattered;
    } else if (channel.config.process == CollisionProcessKind::Elastic &&
        config_.neutral_mass > 0.0) {
        const double total_mass =
            particle_mass_ + config_.neutral_mass;
        velocity =
            (particle_mass_ * velocity +
             config_.neutral_mass * neutral_velocity) /
                total_mass +
            (config_.neutral_mass / total_mass) * scattered;
    } else {
        velocity = neutral_velocity + scattered;
    }
}

void NullCollisionModel::apply_channel(
    std::size_t channel_index,
    Vec3& velocity,
    const Vec3& neutral_velocity,
    std::mt19937_64& rng,
    std::optional<Vec3>& ionization_secondary_relative) const {
    const auto& channel = channels_.at(channel_index);
    if (channel.config.process == CollisionProcessKind::Attachment) {
        velocity = {};
        return;
    }
    if (channel.config.process ==
        CollisionProcessKind::ChargeExchange) {
        velocity = neutral_velocity;
        return;
    }
    const Vec3 initial_relative =
        subtract(velocity, neutral_velocity);
    const double initial_speed = speed(initial_relative);
    double energy =
        0.5 * particle_mass_ * initial_speed * initial_speed;
    if (channel.config.process == CollisionProcessKind::Excitation ||
        channel.config.process == CollisionProcessKind::Ionization) {
        energy = std::max(
            0.0, energy - channel.config.threshold_energy);
    }
    if (channel.config.process == CollisionProcessKind::Ionization &&
        channel.config.ionization_kinematics ==
            IonizationKinematicsKind::OpalBeatyPeterson) {
        const auto pair = opal_ionization_velocities(
            initial_relative, energy, particle_mass_,
            channel.config.ionization_ejected_energy_scale, rng);
        if (channel.config.inelastic_transform ==
            InelasticTransformKind::FiniteMassCenterOfMass) {
            const double total_mass =
                particle_mass_ + config_.neutral_mass;
            const double initial_factor =
                particle_mass_ / total_mass;
            const double relative_factor =
                config_.neutral_mass / total_mass;
            const Vec3 drift_relative{
                initial_factor * initial_relative.x,
                initial_factor * initial_relative.y,
                initial_factor * initial_relative.z};
            const Vec3 primary_relative{
                drift_relative.x + relative_factor * pair.primary.x,
                drift_relative.y + relative_factor * pair.primary.y,
                drift_relative.z + relative_factor * pair.primary.z};
            ionization_secondary_relative = Vec3{
                drift_relative.x + relative_factor * pair.secondary.x,
                drift_relative.y + relative_factor * pair.secondary.y,
                drift_relative.z + relative_factor * pair.secondary.z};
            velocity = add(neutral_velocity, primary_relative);
        } else {
            velocity = add(neutral_velocity, pair.primary);
            ionization_secondary_relative = pair.secondary;
        }
        return;
    }
    if (channel.config.process == CollisionProcessKind::Ionization) {
        energy *= 0.5;
    }
    const double speed = std::sqrt(2.0 * energy / particle_mass_);
    const double mean_cosine =
        channel.mean_cosine.has_value()
            ? channel.mean_cosine->evaluate(
                  collision_energy(channel, initial_speed))
            : 0.0;
    const Vec3 scattered_relative =
        channel.config.angular_scattering ==
                AngularScatteringKind::Backward
            ? Vec3{
                  -initial_relative.x,
                  -initial_relative.y,
                  -initial_relative.z}
            : channel.mean_cosine.has_value()
            ? angular_velocity(
                  initial_relative, speed, mean_cosine, rng)
            : isotropic_velocity(speed, rng);
    if (channel.config.inelastic_transform ==
        InelasticTransformKind::FiniteMassCenterOfMass) {
        const double total_mass =
            particle_mass_ + config_.neutral_mass;
        const double initial_factor =
            particle_mass_ / total_mass;
        const double relative_factor =
            config_.neutral_mass / total_mass;
        const Vec3 drift_relative{
            initial_factor * initial_relative.x,
            initial_factor * initial_relative.y,
            initial_factor * initial_relative.z};
        const Vec3 primary_relative{
            drift_relative.x + relative_factor * scattered_relative.x,
            drift_relative.y + relative_factor * scattered_relative.y,
            drift_relative.z + relative_factor * scattered_relative.z};
        velocity = add(neutral_velocity, primary_relative);
        if (channel.config.process ==
            CollisionProcessKind::Ionization) {
            const Vec3 secondary_scattered =
                isotropic_velocity(speed, rng);
            ionization_secondary_relative = Vec3{
                drift_relative.x +
                    relative_factor * secondary_scattered.x,
                drift_relative.y +
                    relative_factor * secondary_scattered.y,
                drift_relative.z +
                    relative_factor * secondary_scattered.z};
        }
    } else if (channel.config.process == CollisionProcessKind::Elastic &&
        config_.neutral_mass > 0.0) {
        const double total_mass =
            particle_mass_ + config_.neutral_mass;
        const double relative_factor =
            config_.neutral_mass / total_mass;
        velocity = {
            (particle_mass_ * velocity.x +
             config_.neutral_mass * neutral_velocity.x) /
                    total_mass +
                relative_factor * scattered_relative.x,
            (particle_mass_ * velocity.y +
             config_.neutral_mass * neutral_velocity.y) /
                    total_mass +
                relative_factor * scattered_relative.y,
            (particle_mass_ * velocity.z +
             config_.neutral_mass * neutral_velocity.z) /
                    total_mass +
                relative_factor * scattered_relative.z};
    } else {
        velocity = add(neutral_velocity, scattered_relative);
    }
}

CollisionStepStatistics NullCollisionModel::collide(
    double& velocity,
    double timestep,
    std::mt19937_64& rng) const {
    CollisionWorkspace workspace;
    collide_reusing_storage(
        velocity, timestep, rng, workspace);
    return std::move(workspace.statistics);
}

CollisionStepStatistics& NullCollisionModel::collide_reusing_storage(
    double& velocity,
    double timestep,
    std::mt19937_64& rng,
    CollisionWorkspace& workspace) const {
    auto& statistics = workspace.statistics;
    auto& channel_rates = workspace.channel_rates;
    if (std::any_of(
            channels_.begin(), channels_.end(),
            [](const Channel& channel) {
                return channel.config.process ==
                           CollisionProcessKind::Ionization ||
                       channel.config.angular_scattering !=
                           AngularScatteringKind::Isotropic;
            })) {
        throw std::logic_error(
            "ionization and anisotropic MCC require the 3V "
            "collision interface");
    }
    if (!std::isfinite(timestep) || !(timestep > 0.0)) {
        throw std::invalid_argument(
            "MCC timestep must be positive and finite");
    }
    statistics.candidates = 0;
    statistics.null_collisions = 0;
    statistics.channel_collisions.assign(channels_.size(), 0);
    statistics.channel_projectile_energy_change.assign(
        channels_.size(), 0.0);
    statistics.secondaries.clear();
    statistics.primary_removal_channel.reset();
    statistics.primary_removal_product_velocity.reset();
    double elapsed = 0.0;
    while (true) {
        validate_frequency_bound(
            std::abs(velocity), channel_rates);
        const double waiting_time =
            -std::log(open_unit_interval(rng)) / config_.max_frequency;
        if (!std::isfinite(waiting_time) ||
            waiting_time >= timestep - elapsed) {
            break;
        }
        elapsed += waiting_time;
        ++statistics.candidates;
        if (statistics.candidates >
            config_.max_candidates_per_particle) {
            throw std::runtime_error(
                "MCC candidate limit exceeded; reduce dt or increase "
                "max_candidates_per_particle");
        }
        const double neutral_velocity =
            sample_neutral_velocity(rng);
        if (neutral_velocity_stddev_ != 0.0) {
            (void)rates_for_speed(
                std::abs(velocity - neutral_velocity),
                channel_rates);
        }
        const double selection =
            std::generate_canonical<double, 64>(rng) *
            config_.max_frequency;
        double cumulative = 0.0;
        bool accepted = false;
        for (std::size_t channel = 0;
             channel < channel_rates.size(); ++channel) {
            cumulative += channel_rates[channel];
            if (selection < cumulative) {
                const double energy_before =
                    0.5 * particle_mass_ * velocity * velocity;
                apply_channel(
                    channel, velocity, neutral_velocity, rng);
                statistics.channel_projectile_energy_change[channel] +=
                    0.5 * particle_mass_ * velocity * velocity - energy_before;
                ++statistics.channel_collisions[channel];
                if (channels_[channel].config.process ==
                    CollisionProcessKind::Attachment) {
                    statistics.primary_removal_channel = channel;
                    statistics.primary_removal_product_velocity =
                        Vec3{neutral_velocity, 0.0, 0.0};
                }
                accepted = true;
                break;
            }
        }
        if (!accepted) ++statistics.null_collisions;
        if (statistics.primary_removal_channel) break;
    }
    return statistics;
}

CollisionStepStatistics NullCollisionModel::collide(
    Vec3& velocity,
    double timestep,
    std::mt19937_64& rng) const {
    CollisionWorkspace workspace;
    collide_reusing_storage(
        velocity, timestep, rng, workspace);
    return std::move(workspace.statistics);
}

CollisionStepStatistics& NullCollisionModel::collide_reusing_storage(
    Vec3& velocity,
    double timestep,
    std::mt19937_64& rng,
    CollisionWorkspace& workspace) const {
    auto& statistics = workspace.statistics;
    auto& channel_rates = workspace.channel_rates;
    if (!std::isfinite(velocity.x) ||
        !std::isfinite(velocity.y) ||
        !std::isfinite(velocity.z)) {
        throw std::invalid_argument(
            "MCC particle velocity must be finite");
    }
    if (!std::isfinite(timestep) || !(timestep > 0.0)) {
        throw std::invalid_argument(
            "MCC timestep must be positive and finite");
    }
    statistics.candidates = 0;
    statistics.null_collisions = 0;
    statistics.channel_collisions.assign(channels_.size(), 0);
    statistics.channel_projectile_energy_change.assign(
        channels_.size(), 0.0);
    statistics.secondaries.clear();
    statistics.primary_removal_channel.reset();
    statistics.primary_removal_product_velocity.reset();
    const auto projectile_speed = [&]() {
        return speed(velocity);
    };
    double elapsed = 0.0;
    while (true) {
        validate_frequency_bound(
            projectile_speed(), channel_rates);
        const double waiting_time =
            -std::log(open_unit_interval(rng)) / config_.max_frequency;
        if (!std::isfinite(waiting_time) ||
            waiting_time >= timestep - elapsed) {
            break;
        }
        elapsed += waiting_time;
        ++statistics.candidates;
        if (statistics.candidates >
            config_.max_candidates_per_particle) {
            throw std::runtime_error(
                "MCC candidate limit exceeded; reduce dt or increase "
                "max_candidates_per_particle");
        }
        const Vec3 neutral_velocity =
            sample_neutral_velocity_3v(rng);
        if (neutral_velocity_stddev_ != 0.0) {
            (void)rates_for_speed(
                speed(subtract(velocity, neutral_velocity)),
                channel_rates);
        }
        const double selection =
            std::generate_canonical<double, 64>(rng) *
            config_.max_frequency;
        double cumulative = 0.0;
        bool accepted = false;
        for (std::size_t channel = 0;
             channel < channel_rates.size(); ++channel) {
            cumulative += channel_rates[channel];
            if (selection < cumulative) {
                const double energy_before = 0.5 * particle_mass_ *
                    (velocity.x * velocity.x + velocity.y * velocity.y +
                     velocity.z * velocity.z);
                std::optional<Vec3> ionization_secondary_relative;
                apply_channel(
                    channel, velocity, neutral_velocity, rng,
                    ionization_secondary_relative);
                statistics.channel_projectile_energy_change[channel] +=
                    0.5 * particle_mass_ *
                        (velocity.x * velocity.x + velocity.y * velocity.y +
                         velocity.z * velocity.z) - energy_before;
                if (channels_[channel].config.process ==
                    CollisionProcessKind::Ionization) {
                    const double secondary_speed = speed(
                        subtract(velocity, neutral_velocity));
                    const Vec3 secondary_relative =
                        ionization_secondary_relative
                            ? *ionization_secondary_relative
                            : isotropic_velocity(secondary_speed, rng);
                    statistics.secondaries.push_back({
                        channel,
                        add(neutral_velocity, secondary_relative),
                        neutral_velocity});
                }
                ++statistics.channel_collisions[channel];
                if (channels_[channel].config.process ==
                    CollisionProcessKind::Attachment) {
                    statistics.primary_removal_channel = channel;
                    statistics.primary_removal_product_velocity =
                        neutral_velocity;
                }
                accepted = true;
                break;
            }
        }
        if (!accepted) ++statistics.null_collisions;
        if (statistics.primary_removal_channel) break;
    }
    return statistics;
}

} // namespace pic
