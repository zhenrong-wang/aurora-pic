#pragma once

#include "pic/Config.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <random>
#include <string>
#include <vector>

namespace pic {

class CrossSectionTable {
public:
    CrossSectionTable(
        const std::filesystem::path& path,
        double energy_scale = 1.0,
        double cross_section_scale = 1.0,
        CrossSectionInterpolationKind interpolation =
            CrossSectionInterpolationKind::Linear);

    double evaluate(double energy) const;
    double maximum_between(
        double minimum_energy, double maximum_energy) const;
    const std::vector<double>& energies() const { return energies_; }
    const std::vector<double>& cross_sections() const {
        return cross_sections_;
    }

private:
    std::vector<double> energies_{};
    std::vector<double> cross_sections_{};
    CrossSectionInterpolationKind interpolation_{
        CrossSectionInterpolationKind::Linear};
    std::size_t maximum_tree_leaf_count_{0};
    std::vector<double> maximum_tree_{};
};

class MeanCosineTable {
public:
    MeanCosineTable(
        const std::filesystem::path& path,
        double energy_scale = 1.0);

    double evaluate(double energy) const;
    const std::vector<double>& energies() const { return energies_; }
    const std::vector<double>& mean_cosines() const {
        return mean_cosines_;
    }

private:
    std::vector<double> energies_{};
    std::vector<double> mean_cosines_{};
};

struct CollisionStepStatistics {
    std::uint64_t candidates{0};
    std::uint64_t null_collisions{0};
    std::vector<std::uint64_t> channel_collisions{};
    std::vector<double> channel_projectile_energy_change{};
    struct Secondary {
        std::size_t channel{0};
        Vec3 velocity{};
        Vec3 ion_velocity{};
    };
    std::vector<Secondary> secondaries{};
    std::optional<std::size_t> primary_removal_channel{};
    std::optional<Vec3> primary_removal_product_velocity{};
};

struct CollisionWorkspace {
    CollisionStepStatistics statistics{};
    std::vector<double> channel_rates{};
};

struct CollisionDiagnostics {
    std::uint64_t candidates{0};
    std::uint64_t null_collisions{0};
    std::vector<std::string> channel_names{};
    std::vector<std::uint64_t> channel_collisions{};
    std::vector<double> channel_energy_change{};
};

class NullCollisionModel {
public:
    NullCollisionModel(CollisionConfig config, double particle_mass);

    CollisionStepStatistics collide(
        double& velocity,
        double timestep,
        std::mt19937_64& rng) const;
    CollisionStepStatistics collide(
        Vec3& velocity,
        double timestep,
        std::mt19937_64& rng) const;
    CollisionStepStatistics& collide_reusing_storage(
        double& velocity,
        double timestep,
        std::mt19937_64& rng,
        CollisionWorkspace& workspace) const;
    CollisionStepStatistics& collide_reusing_storage(
        Vec3& velocity,
        double timestep,
        std::mt19937_64& rng,
        CollisionWorkspace& workspace) const;

    const CollisionConfig& config() const { return config_; }
    const std::vector<std::string>& channel_names() const {
        return channel_names_;
    }
    std::uint64_t signature() const { return signature_; }
    double neutral_velocity_stddev() const {
        return neutral_velocity_stddev_;
    }
    double neutral_speed_limit_sigma() const {
        return neutral_velocity_stddev_ > 0.0
                   ? neutral_speed_limit_sigma_
                   : 0.0;
    }

private:
    struct Channel {
        CollisionChannelConfig config{};
        CrossSectionTable table;

        explicit Channel(const CollisionChannelConfig& channel)
            : config(channel),
              table(channel.cross_section_file,
                    channel.energy_scale,
                    channel.cross_section_scale,
                    channel.cross_section_interpolation) {
            if (channel.angular_scattering ==
                AngularScatteringKind::HenyeyGreenstein) {
                mean_cosine.emplace(
                    channel.mean_cosine_file,
                    channel.mean_cosine_energy_scale);
            }
        }

        std::optional<MeanCosineTable> mean_cosine{};
    };

    double rates_for_speed(
        double speed, std::vector<double>& rates) const;
    double collision_energy(
        const Channel& channel, double relative_speed) const;
    void validate_frequency_bound(
        double projectile_speed,
        std::vector<double>& rate_scratch) const;
    double sample_neutral_velocity(std::mt19937_64& rng) const;
    Vec3 sample_neutral_velocity_3v(std::mt19937_64& rng) const;
    void apply_channel(
        std::size_t channel,
        double& velocity,
        double neutral_velocity,
        std::mt19937_64& rng) const;
    void apply_channel(
        std::size_t channel,
        Vec3& velocity,
        const Vec3& neutral_velocity,
        std::mt19937_64& rng,
        std::optional<Vec3>& ionization_secondary_relative) const;

    CollisionConfig config_{};
    double particle_mass_{0.0};
    std::vector<Channel> channels_{};
    std::vector<std::string> channel_names_{};
    std::uint64_t signature_{0};
    double neutral_velocity_stddev_{0.0};
    static constexpr double neutral_speed_limit_sigma_{8.0};
};

} // namespace pic
