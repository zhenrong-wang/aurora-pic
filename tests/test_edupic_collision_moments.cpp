#include "pic/Collision.hpp"

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

namespace {
constexpr std::size_t target_collisions = 12000;
constexpr std::size_t maximum_attempts = 1500000;

double norm(const pic::Vec3& value) {
    return std::sqrt(
        value.x * value.x + value.y * value.y + value.z * value.z);
}

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void require_near(
    double actual, double expected, double tolerance,
    const std::string& message) {
    if (std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(
            message + ": actual=" + std::to_string(actual) +
            ", expected=" + std::to_string(expected));
    }
}

pic::CollisionConfig base_config(
    const std::filesystem::path& table, double neutral_mass,
    double maximum_frequency) {
    pic::CollisionConfig config;
    config.enabled = true;
    config.model = pic::CollisionModelKind::NullCollision;
    config.neutral_density = 1.0;
    config.neutral_mass = neutral_mass;
    config.species = "frozen_test";
    config.max_frequency = maximum_frequency;
    config.max_candidates_per_particle = 8;
    pic::CollisionChannelConfig channel;
    channel.name = "channel";
    channel.cross_section_file = table;
    config.channels = {channel};
    return config;
}

struct AngularMoments {
    std::size_t accepted{0};
    double mean_cosine{0.0};
    double mean_squared_cosine{0.0};
};

AngularMoments sample_isotropic(
    const pic::NullCollisionModel& model, double timestep,
    const pic::Vec3& incoming, double drift_factor,
    double relative_factor, double post_relative_speed,
    std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    double cosine_sum = 0.0;
    double cosine_squared_sum = 0.0;
    std::size_t accepted = 0;
    for (std::size_t attempt = 0;
         attempt < maximum_attempts && accepted < target_collisions;
         ++attempt) {
        pic::Vec3 projectile = incoming;
        const auto statistics = model.collide(projectile, timestep, rng);
        if (statistics.channel_collisions[0] != 1) continue;
        const pic::Vec3 relative{
            (projectile.x - drift_factor * incoming.x) / relative_factor,
            (projectile.y - drift_factor * incoming.y) / relative_factor,
            (projectile.z - drift_factor * incoming.z) / relative_factor};
        require_near(
            norm(relative), post_relative_speed, 2e-12,
            "post-collision relative speed differs from eduPIC law");
        const double cosine =
            relative.x / post_relative_speed;
        cosine_sum += cosine;
        cosine_squared_sum += cosine * cosine;
        ++accepted;
    }
    require(
        accepted == target_collisions,
        "insufficient accepted collisions for moment audit");
    return {
        accepted,
        cosine_sum / static_cast<double>(accepted),
        cosine_squared_sum / static_cast<double>(accepted)};
}

struct IonizationMoments {
    std::size_t accepted{0};
    double mean_ejected_energy{0.0};
    double expected_mean_ejected_energy{0.0};
    double maximum_pair_energy_error{0.0};
    double maximum_pair_momentum_error{0.0};
};

IonizationMoments sample_ionization(
    const pic::NullCollisionModel& model, double timestep,
    const pic::Vec3& incoming, double drift_factor,
    double relative_factor, double threshold, double scale,
    std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    const double initial_energy = 0.5 * norm(incoming) * norm(incoming);
    const double available_energy = initial_energy - threshold;
    double ejected_energy_sum = 0.0;
    double maximum_energy_error = 0.0;
    double maximum_momentum_error = 0.0;
    std::size_t accepted = 0;
    for (std::size_t attempt = 0;
         attempt < maximum_attempts && accepted < target_collisions;
         ++attempt) {
        pic::Vec3 primary = incoming;
        const auto statistics = model.collide(primary, timestep, rng);
        if (statistics.channel_collisions[0] != 1 ||
            statistics.secondaries.size() != 1) {
            continue;
        }
        const auto recover = [&](const pic::Vec3& transformed) {
            return pic::Vec3{
                (transformed.x - drift_factor * incoming.x) /
                    relative_factor,
                (transformed.y - drift_factor * incoming.y) /
                    relative_factor,
                (transformed.z - drift_factor * incoming.z) /
                    relative_factor};
        };
        const pic::Vec3 primary_relative = recover(primary);
        const pic::Vec3 secondary_relative = recover(
            statistics.secondaries.front().velocity);
        const double primary_energy =
            0.5 * norm(primary_relative) * norm(primary_relative);
        const double secondary_energy =
            0.5 * norm(secondary_relative) * norm(secondary_relative);
        ejected_energy_sum += secondary_energy;
        maximum_energy_error = std::max(
            maximum_energy_error,
            std::abs(primary_energy + secondary_energy - available_energy));
        const pic::Vec3 momentum_error{
            primary_relative.x + secondary_relative.x -
                std::sqrt(2.0 * available_energy),
            primary_relative.y + secondary_relative.y,
            primary_relative.z + secondary_relative.z};
        maximum_momentum_error = std::max(
            maximum_momentum_error, norm(momentum_error));
        require(
            secondary_energy <= 0.5 * available_energy + 2e-12,
            "Opal ejected energy exceeds the eduPIC support");
        ++accepted;
    }
    require(
        accepted == target_collisions,
        "insufficient ionization collisions for moment audit");
    const double alpha = std::atan(available_energy / (2.0 * scale));
    const double expected_mean =
        scale * (-std::log(std::cos(alpha))) / alpha;
    const double sampled_mean =
        ejected_energy_sum / static_cast<double>(accepted);
    require_near(
        sampled_mean, expected_mean, 0.02 * expected_mean,
        "Opal mean ejected energy differs from the eduPIC distribution");
    require(
        maximum_energy_error < 2e-12,
        "Opal electron-pair energy does not close");
    require(
        maximum_momentum_error < 2e-12,
        "Opal paired-angle momentum invariant does not close");
    return {
        accepted, sampled_mean, expected_mean, maximum_energy_error,
        maximum_momentum_error};
}
} // namespace

int main() {
    const std::filesystem::path table =
        "test_edupic_collision_moments_constant.dat";
    try {
        {
            std::ofstream output(table);
            output << "0 1\n100 1\n";
        }
        constexpr pic::Vec3 incoming{4.0, 0.0, 0.0};
        constexpr double projectile_mass = 1.0;
        constexpr double neutral_mass = 3.0;
        constexpr double drift_factor =
            projectile_mass / (projectile_mass + neutral_mass);
        constexpr double relative_factor =
            neutral_mass / (projectile_mass + neutral_mass);

        auto elastic_config = base_config(table, neutral_mass, 4.0);
        elastic_config.channels.front().process =
            pic::CollisionProcessKind::Elastic;
        const pic::NullCollisionModel elastic(
            elastic_config, projectile_mass);
        const auto elastic_moments = sample_isotropic(
            elastic, 0.003, incoming, drift_factor, relative_factor,
            norm(incoming), 51949);

        auto excitation_config = base_config(table, neutral_mass, 4.0);
        excitation_config.channels.front().process =
            pic::CollisionProcessKind::Excitation;
        excitation_config.channels.front().threshold_energy = 0.5;
        excitation_config.channels.front().inelastic_transform =
            pic::InelasticTransformKind::FiniteMassCenterOfMass;
        const pic::NullCollisionModel excitation(
            excitation_config, projectile_mass);
        const auto excitation_moments = sample_isotropic(
            excitation, 0.003, incoming, drift_factor, relative_factor,
            std::sqrt(15.0), 63059);

        for (const auto& [label, moments] : {
                 std::pair{"elastic", elastic_moments},
                 std::pair{"excitation", excitation_moments}}) {
            require_near(
                moments.mean_cosine, 0.0, 0.02,
                std::string(label) + " mean scattering cosine differs");
            require_near(
                moments.mean_squared_cosine, 1.0 / 3.0, 0.02,
                std::string(label) +
                    " second scattering-cosine moment differs");
        }

        auto ionization_config = base_config(table, neutral_mass, 4.0);
        ionization_config.channels.front().process =
            pic::CollisionProcessKind::Ionization;
        ionization_config.channels.front().threshold_energy = 0.5;
        ionization_config.channels.front().secondary_species = "electrons";
        ionization_config.channels.front().ion_species = "ions";
        ionization_config.channels.front().ionization_kinematics =
            pic::IonizationKinematicsKind::OpalBeatyPeterson;
        ionization_config.channels.front().ionization_ejected_energy_scale =
            1.0;
        ionization_config.channels.front().inelastic_transform =
            pic::InelasticTransformKind::FiniteMassCenterOfMass;
        const pic::NullCollisionModel ionization(
            ionization_config, projectile_mass);
        const auto ionization_moments = sample_ionization(
            ionization, 0.003, incoming, drift_factor, relative_factor,
            0.5, 1.0, 2718281);

        constexpr pic::Vec3 ion_incoming{2.0, 0.0, 0.0};
        auto ion_isotropic_config = base_config(table, 1.0, 2.0);
        ion_isotropic_config.channels.front().process =
            pic::CollisionProcessKind::Elastic;
        ion_isotropic_config.channels.front().energy_frame =
            pic::CollisionEnergyFrame::CenterOfMass;
        const pic::NullCollisionModel ion_isotropic(
            ion_isotropic_config, 1.0);
        const auto ion_moments = sample_isotropic(
            ion_isotropic, 0.006, ion_incoming, 0.5, 0.5,
            norm(ion_incoming), 173205);
        require_near(
            ion_moments.mean_cosine, 0.0, 0.02,
            "ion isotropic mean scattering cosine differs");
        require_near(
            ion_moments.mean_squared_cosine, 1.0 / 3.0, 0.02,
            "ion isotropic second scattering-cosine moment differs");

        auto backward_config = ion_isotropic_config;
        backward_config.channels.front().angular_scattering =
            pic::AngularScatteringKind::Backward;
        const pic::NullCollisionModel backward(backward_config, 1.0);
        std::mt19937_64 backward_rng(141421);
        bool observed_backward = false;
        for (std::size_t attempt = 0;
             attempt < maximum_attempts && !observed_backward; ++attempt) {
            pic::Vec3 projectile = ion_incoming;
            const auto statistics = backward.collide(
                projectile, 0.1, backward_rng);
            if (statistics.channel_collisions[0] != 1) continue;
            require_near(
                norm(projectile), 0.0, 2e-14,
                "equal-mass backward ion collision did not exchange velocity");
            observed_backward = true;
        }
        require(observed_backward, "backward ion collision was not sampled");

        std::filesystem::remove(table);
        std::cout << std::setprecision(17)
                  << "{\n"
                  << "  \"samples_per_stochastic_channel\": "
                  << target_collisions << ",\n"
                  << "  \"elastic_mean_cosine\": "
                  << elastic_moments.mean_cosine << ",\n"
                  << "  \"elastic_mean_squared_cosine\": "
                  << elastic_moments.mean_squared_cosine << ",\n"
                  << "  \"excitation_mean_cosine\": "
                  << excitation_moments.mean_cosine << ",\n"
                  << "  \"excitation_mean_squared_cosine\": "
                  << excitation_moments.mean_squared_cosine << ",\n"
                  << "  \"opal_mean_ejected_energy\": "
                  << ionization_moments.mean_ejected_energy << ",\n"
                  << "  \"opal_expected_mean_ejected_energy\": "
                  << ionization_moments.expected_mean_ejected_energy
                  << ",\n"
                  << "  \"opal_mean_relative_error\": "
                  << std::abs(
                         ionization_moments.mean_ejected_energy /
                             ionization_moments.expected_mean_ejected_energy -
                         1.0)
                  << ",\n"
                  << "  \"opal_maximum_pair_energy_error\": "
                  << ionization_moments.maximum_pair_energy_error << ",\n"
                  << "  \"opal_maximum_pair_momentum_error\": "
                  << ionization_moments.maximum_pair_momentum_error << ",\n"
                  << "  \"ion_isotropic_mean_cosine\": "
                  << ion_moments.mean_cosine << ",\n"
                  << "  \"ion_isotropic_mean_squared_cosine\": "
                  << ion_moments.mean_squared_cosine << "\n"
                  << "}\n";
        return 0;
    } catch (const std::exception& error) {
        std::filesystem::remove(table);
        std::cerr << error.what() << '\n';
        return 1;
    }
}
