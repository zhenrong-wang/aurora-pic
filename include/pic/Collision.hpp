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
        double cross_section_scale = 1.0);

    double evaluate(double energy) const;
    const std::vector<double>& energies() const { return energies_; }
    const std::vector<double>& cross_sections() const {
        return cross_sections_;
    }

private:
    std::vector<double> energies_{};
    std::vector<double> cross_sections_{};
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
    struct Secondary {
        std::size_t channel{0};
        Vec3 velocity{};
    };
    std::vector<Secondary> secondaries{};
};

struct CollisionDiagnostics {
    std::uint64_t candidates{0};
    std::uint64_t null_collisions{0};
    std::vector<std::string> channel_names{};
    std::vector<std::uint64_t> channel_collisions{};
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

    const CollisionConfig& config() const { return config_; }
    const std::vector<std::string>& channel_names() const {
        return channel_names_;
    }
    std::uint64_t signature() const { return signature_; }

private:
    struct Channel {
        CollisionChannelConfig config{};
        CrossSectionTable table;

        explicit Channel(const CollisionChannelConfig& channel)
            : config(channel),
              table(channel.cross_section_file,
                    channel.energy_scale,
                    channel.cross_section_scale) {
            if (channel.angular_scattering ==
                AngularScatteringKind::HenyeyGreenstein) {
                mean_cosine.emplace(
                    channel.mean_cosine_file,
                    channel.mean_cosine_energy_scale);
            }
        }

        std::optional<MeanCosineTable> mean_cosine{};
    };

    std::vector<double> rates(double velocity) const;
    std::vector<double> rates_for_speed(double speed) const;
    void apply_channel(
        std::size_t channel,
        double& velocity,
        std::mt19937_64& rng) const;
    void apply_channel(
        std::size_t channel,
        Vec3& velocity,
        std::mt19937_64& rng) const;

    CollisionConfig config_{};
    double particle_mass_{0.0};
    std::vector<Channel> channels_{};
    std::vector<std::string> channel_names_{};
    std::uint64_t signature_{0};
};

} // namespace pic
