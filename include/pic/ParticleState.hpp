#pragma once

#include "pic/Types.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace pic {

struct ExternalParticleRecord {
    Vec3 position{};
    Vec3 velocity{};
};

struct ExternalParticleState {
    std::size_t version{1};
    std::size_t spatial_dimension{0};
    std::size_t velocity_dimensions{0};
    UnitSystem unit_system{UnitSystem::Normalized};
    std::map<std::string, std::vector<ExternalParticleRecord>>
        species;
    std::size_t particle_count{0};
    std::uint64_t signature{0};
};

struct ExternalParticleStateMetadata {
    std::size_t version{1};
    std::size_t spatial_dimension{0};
    std::size_t velocity_dimensions{0};
    UnitSystem unit_system{UnitSystem::Normalized};
    std::size_t particle_count{0};
    std::uint64_t signature{0};
};

struct ExternalSpeciesExpectation {
    std::string name;
    std::size_t particle_count{0};
};

using ExternalParticleRecordConsumer = std::function<void(
    std::size_t species_index,
    std::size_t record_index,
    const ExternalParticleRecord& record)>;

ExternalParticleState load_external_particle_state(
    const std::filesystem::path& path,
    std::size_t max_particles);

void validate_external_particle_state(
    const ExternalParticleState& state,
    std::size_t spatial_dimension,
    std::size_t velocity_dimensions,
    UnitSystem unit_system,
    const std::vector<ExternalSpeciesExpectation>& expected_species,
    const std::string& context);

inline void validate_external_particle_state(
    const ExternalParticleState& state,
    std::size_t spatial_dimension,
    UnitSystem unit_system,
    const std::vector<ExternalSpeciesExpectation>& expected_species,
    const std::string& context) {
    validate_external_particle_state(
        state, spatial_dimension,
        spatial_dimension == 1 ? 1 : 3,
        unit_system, expected_species, context);
}

ExternalParticleState load_validated_external_particle_state(
    const std::filesystem::path& path,
    std::size_t spatial_dimension,
    std::size_t velocity_dimensions,
    UnitSystem unit_system,
    const std::vector<ExternalSpeciesExpectation>& expected_species,
    const std::string& context,
    std::optional<std::uint64_t> expected_signature = {});

inline ExternalParticleState load_validated_external_particle_state(
    const std::filesystem::path& path,
    std::size_t spatial_dimension,
    UnitSystem unit_system,
    const std::vector<ExternalSpeciesExpectation>& expected_species,
    const std::string& context,
    std::optional<std::uint64_t> expected_signature = {}) {
    return load_validated_external_particle_state(
        path, spatial_dimension,
        spatial_dimension == 1 ? 1 : 3,
        unit_system, expected_species, context,
        expected_signature);
}

ExternalParticleStateMetadata
load_validated_external_particle_state_bounded(
    const std::filesystem::path& path,
    std::size_t spatial_dimension,
    std::size_t velocity_dimensions,
    UnitSystem unit_system,
    const std::vector<ExternalSpeciesExpectation>& expected_species,
    const std::string& context,
    const ExternalParticleRecordConsumer& consumer,
    std::optional<std::uint64_t> expected_signature = {});

inline ExternalParticleStateMetadata
load_validated_external_particle_state_bounded(
    const std::filesystem::path& path,
    std::size_t spatial_dimension,
    UnitSystem unit_system,
    const std::vector<ExternalSpeciesExpectation>& expected_species,
    const std::string& context,
    const ExternalParticleRecordConsumer& consumer,
    std::optional<std::uint64_t> expected_signature = {}) {
    return load_validated_external_particle_state_bounded(
        path, spatial_dimension,
        spatial_dimension == 1 ? 1 : 3,
        unit_system, expected_species, context, consumer,
        expected_signature);
}

std::uint64_t external_particle_state_signature(
    const ExternalParticleState& state);

ExternalParticleStateMetadata external_particle_state_metadata(
    const ExternalParticleState& state);

void write_external_particle_state(
    const std::filesystem::path& path,
    const ExternalParticleState& state,
    bool overwrite = false);

void write_external_particle_state_metadata(
    const std::filesystem::path& path,
    const std::filesystem::path& source_path,
    const ExternalParticleStateMetadata& metadata,
    std::optional<std::uint64_t> expected_signature = {});

} // namespace pic
