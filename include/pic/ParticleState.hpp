#pragma once

#include "pic/Types.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
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
    UnitSystem unit_system{UnitSystem::Normalized};
    std::map<std::string, std::vector<ExternalParticleRecord>>
        species;
    std::size_t particle_count{0};
    std::uint64_t signature{0};
};

struct ExternalParticleStateMetadata {
    std::size_t version{1};
    std::size_t spatial_dimension{0};
    UnitSystem unit_system{UnitSystem::Normalized};
    std::size_t particle_count{0};
    std::uint64_t signature{0};
};

struct ExternalSpeciesExpectation {
    std::string name;
    std::size_t particle_count{0};
};

ExternalParticleState load_external_particle_state(
    const std::filesystem::path& path,
    std::size_t max_particles);

void validate_external_particle_state(
    const ExternalParticleState& state,
    std::size_t spatial_dimension,
    UnitSystem unit_system,
    const std::vector<ExternalSpeciesExpectation>& expected_species,
    const std::string& context);

ExternalParticleState load_validated_external_particle_state(
    const std::filesystem::path& path,
    std::size_t spatial_dimension,
    UnitSystem unit_system,
    const std::vector<ExternalSpeciesExpectation>& expected_species,
    const std::string& context,
    std::optional<std::uint64_t> expected_signature = {});

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
