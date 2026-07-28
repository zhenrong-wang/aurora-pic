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

} // namespace

CrossSectionTable::CrossSectionTable(
    const std::filesystem::path& path,
    double energy_scale,
    double cross_section_scale) {
    if (!std::isfinite(energy_scale) || !(energy_scale > 0.0)) {
        throw std::invalid_argument(
            "collision cross-section energy_scale must be positive and finite");
    }
    if (!std::isfinite(cross_section_scale) ||
        !(cross_section_scale > 0.0)) {
        throw std::invalid_argument(
            "collision cross_section_scale must be positive and finite");
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
}

double CrossSectionTable::evaluate(double energy) const {
    if (!std::isfinite(energy) || energy < 0.0) {
        throw std::invalid_argument(
            "collision lookup energy must be finite and non-negative");
    }
    if (energy <= energies_.front()) return cross_sections_.front();
    if (energy >= energies_.back()) return cross_sections_.back();
    const auto upper =
        std::upper_bound(energies_.begin(), energies_.end(), energy);
    const std::size_t high =
        static_cast<std::size_t>(upper - energies_.begin());
    const std::size_t low = high - 1;
    const double fraction =
        (energy - energies_[low]) / (energies_[high] - energies_[low]);
    return cross_sections_[low] +
           fraction * (cross_sections_[high] - cross_sections_[low]);
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
        if (config_.neutral_mass > 0.0) {
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
        if ((channel_config.process == CollisionProcessKind::Elastic ||
             channel_config.process ==
                 CollisionProcessKind::ChargeExchange) &&
            channel_config.threshold_energy != 0.0) {
            throw std::invalid_argument(
                "elastic and charge-exchange MCC channel threshold_energy "
                "must be zero");
        }
        if ((channel_config.process == CollisionProcessKind::Excitation ||
             channel_config.process == CollisionProcessKind::Ionization) &&
            !(channel_config.threshold_energy > 0.0)) {
            throw std::invalid_argument(
                "inelastic MCC channel threshold_energy must be positive");
        }
        if (channel_config.process == CollisionProcessKind::Ionization &&
            (channel_config.secondary_species.empty() ||
             channel_config.ion_species.empty())) {
            throw std::invalid_argument(
                "ionization MCC channel requires secondary and ion species");
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
        channel_names_.push_back(channel_config.name);
        hash_string(signature_, channel_config.name);
        hash_string(signature_, to_string(channel_config.process));
        hash_double(signature_, channel_config.threshold_energy);
        if (channel_config.process ==
            CollisionProcessKind::Ionization) {
            hash_string(signature_, channel_config.secondary_species);
            hash_string(signature_, channel_config.ion_species);
        }
        for (std::size_t i = 0;
             i < channels_.back().table.energies().size(); ++i) {
            hash_double(signature_, channels_.back().table.energies()[i]);
            hash_double(
                signature_, channels_.back().table.cross_sections()[i]);
        }
    }
}

std::vector<double> NullCollisionModel::rates(double velocity) const {
    if (!std::isfinite(velocity)) {
        throw std::invalid_argument("MCC particle velocity must be finite");
    }
    return rates_for_speed(std::abs(velocity));
}

std::vector<double> NullCollisionModel::rates_for_speed(double speed) const {
    if (!std::isfinite(speed) || speed < 0.0) {
        throw std::invalid_argument(
            "MCC particle speed must be finite and non-negative");
    }
    const double energy = 0.5 * particle_mass_ * speed * speed;
    if (!std::isfinite(energy)) {
        throw std::overflow_error("MCC particle energy overflow");
    }
    std::vector<double> result;
    result.reserve(channels_.size());
    double total = 0.0;
    for (const auto& channel : channels_) {
        double rate = 0.0;
        if (energy >= channel.config.threshold_energy) {
            rate = config_.neutral_density *
                   channel.table.evaluate(energy) * speed;
        }
        if (!std::isfinite(rate) || rate < 0.0) {
            throw std::overflow_error("MCC collision rate overflow");
        }
        result.push_back(rate);
        total += rate;
    }
    const double tolerance =
        64.0 * std::numeric_limits<double>::epsilon() *
        std::max(config_.max_frequency, total);
    if (total > config_.max_frequency + tolerance) {
        throw std::runtime_error(
            "MCC total collision frequency exceeds configured max_frequency");
    }
    return result;
}

void NullCollisionModel::apply_channel(
    std::size_t channel_index,
    double& velocity,
    std::mt19937_64& rng) const {
    const auto& channel = channels_.at(channel_index);
    if (channel.config.process ==
        CollisionProcessKind::ChargeExchange) {
        velocity = 0.0;
        return;
    }
    double energy = 0.5 * particle_mass_ * velocity * velocity;
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
    if (channel.config.process == CollisionProcessKind::Elastic &&
        config_.neutral_mass > 0.0) {
        const double total_mass =
            particle_mass_ + config_.neutral_mass;
        velocity =
            (particle_mass_ / total_mass) * velocity +
            (config_.neutral_mass / total_mass) * scattered;
    } else {
        velocity = scattered;
    }
}

void NullCollisionModel::apply_channel(
    std::size_t channel_index,
    Vec3& velocity,
    std::mt19937_64& rng) const {
    const auto& channel = channels_.at(channel_index);
    if (channel.config.process ==
        CollisionProcessKind::ChargeExchange) {
        velocity = {};
        return;
    }
    const double initial_speed = std::sqrt(
        velocity.x * velocity.x +
        velocity.y * velocity.y +
        velocity.z * velocity.z);
    double energy =
        0.5 * particle_mass_ * initial_speed * initial_speed;
    if (channel.config.process == CollisionProcessKind::Excitation ||
        channel.config.process == CollisionProcessKind::Ionization) {
        energy = std::max(
            0.0, energy - channel.config.threshold_energy);
    }
    if (channel.config.process == CollisionProcessKind::Ionization) {
        energy *= 0.5;
    }
    const double speed = std::sqrt(2.0 * energy / particle_mass_);
    const Vec3 scattered = isotropic_velocity(speed, rng);
    if (channel.config.process == CollisionProcessKind::Elastic &&
        config_.neutral_mass > 0.0) {
        const double total_mass =
            particle_mass_ + config_.neutral_mass;
        const double center_factor = particle_mass_ / total_mass;
        const double relative_factor =
            config_.neutral_mass / total_mass;
        velocity = {
            center_factor * velocity.x +
                relative_factor * scattered.x,
            center_factor * velocity.y +
                relative_factor * scattered.y,
            center_factor * velocity.z +
                relative_factor * scattered.z};
    } else {
        velocity = scattered;
    }
}

CollisionStepStatistics NullCollisionModel::collide(
    double& velocity,
    double timestep,
    std::mt19937_64& rng) const {
    if (std::any_of(
            channels_.begin(), channels_.end(),
            [](const Channel& channel) {
                return channel.config.process ==
                       CollisionProcessKind::Ionization;
            })) {
        throw std::logic_error(
            "ionization MCC requires the 3V collision interface");
    }
    if (!std::isfinite(timestep) || !(timestep > 0.0)) {
        throw std::invalid_argument(
            "MCC timestep must be positive and finite");
    }
    CollisionStepStatistics statistics;
    statistics.channel_collisions.assign(channels_.size(), 0);
    (void)rates(velocity);
    double elapsed = 0.0;
    while (true) {
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
        const auto channel_rates = rates(velocity);
        const double selection =
            std::generate_canonical<double, 64>(rng) *
            config_.max_frequency;
        double cumulative = 0.0;
        bool accepted = false;
        for (std::size_t channel = 0;
             channel < channel_rates.size(); ++channel) {
            cumulative += channel_rates[channel];
            if (selection < cumulative) {
                apply_channel(channel, velocity, rng);
                ++statistics.channel_collisions[channel];
                accepted = true;
                break;
            }
        }
        if (!accepted) ++statistics.null_collisions;
    }
    return statistics;
}

CollisionStepStatistics NullCollisionModel::collide(
    Vec3& velocity,
    double timestep,
    std::mt19937_64& rng) const {
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
    CollisionStepStatistics statistics;
    statistics.channel_collisions.assign(channels_.size(), 0);
    const auto speed = [&]() {
        return std::sqrt(
            velocity.x * velocity.x +
            velocity.y * velocity.y +
            velocity.z * velocity.z);
    };
    (void)rates_for_speed(speed());
    double elapsed = 0.0;
    while (true) {
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
        const auto channel_rates = rates_for_speed(speed());
        const double selection =
            std::generate_canonical<double, 64>(rng) *
            config_.max_frequency;
        double cumulative = 0.0;
        bool accepted = false;
        for (std::size_t channel = 0;
             channel < channel_rates.size(); ++channel) {
            cumulative += channel_rates[channel];
            if (selection < cumulative) {
                apply_channel(channel, velocity, rng);
                if (channels_[channel].config.process ==
                    CollisionProcessKind::Ionization) {
                    const double secondary_speed = speed();
                    statistics.secondaries.push_back({
                        channel,
                        isotropic_velocity(secondary_speed, rng)});
                }
                ++statistics.channel_collisions[channel];
                accepted = true;
                break;
            }
        }
        if (!accepted) ++statistics.null_collisions;
    }
    return statistics;
}

} // namespace pic
