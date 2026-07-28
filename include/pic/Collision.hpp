#pragma once

#include "pic/Config.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
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

struct CollisionStepStatistics {
    std::uint64_t candidates{0};
    std::uint64_t null_collisions{0};
    std::vector<std::uint64_t> channel_collisions{};
};

class NullCollisionModel {
public:
    NullCollisionModel(CollisionConfig config, double particle_mass);

    CollisionStepStatistics collide(
        double& velocity,
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
                    channel.cross_section_scale) {}
    };

    std::vector<double> rates(double velocity) const;
    void apply_channel(
        std::size_t channel,
        double& velocity,
        std::mt19937_64& rng) const;

    CollisionConfig config_{};
    double particle_mass_{0.0};
    std::vector<Channel> channels_{};
    std::vector<std::string> channel_names_{};
    std::uint64_t signature_{0};
};

} // namespace pic
