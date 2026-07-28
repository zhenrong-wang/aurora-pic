#include "pic/ParticleState.hpp"

#include <cmath>
#include <fstream>
#include <limits>
#include <set>
#include <stdexcept>

namespace pic {
namespace {

constexpr const char* particle_state_magic =
    "AuroraPIC-particle-state-v1";

void require_key(
    std::istream& input,
    const std::string& expected,
    const std::filesystem::path& path) {
    std::string key;
    input >> key;
    if (!input || key != expected) {
        throw std::runtime_error(
            "external particle state '" + path.string() +
            "' expected key '" + expected + "'");
    }
}

UnitSystem parse_unit_system(
    const std::string& value,
    const std::filesystem::path& path) {
    if (value == "normalized") return UnitSystem::Normalized;
    if (value == "si") return UnitSystem::SI;
    throw std::runtime_error(
        "external particle state '" + path.string() +
        "' has invalid units '" + value + "'");
}

bool finite(const ExternalParticleRecord& record) {
    return std::isfinite(record.position.x) &&
           std::isfinite(record.position.y) &&
           std::isfinite(record.position.z) &&
           std::isfinite(record.velocity.x) &&
           std::isfinite(record.velocity.y) &&
           std::isfinite(record.velocity.z);
}

} // namespace

ExternalParticleState load_external_particle_state(
    const std::filesystem::path& path,
    std::size_t max_particles) {
    if (max_particles == 0) {
        throw std::invalid_argument(
            "external particle-state limit must be positive");
    }
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error(
            "cannot open external particle state: " +
            path.string());
    }
    std::string magic;
    input >> magic;
    if (magic != particle_state_magic) {
        throw std::runtime_error(
            "invalid external particle-state magic in: " +
            path.string());
    }

    ExternalParticleState state;
    require_key(input, "dimension", path);
    input >> state.spatial_dimension;
    if (!input || state.spatial_dimension == 0 ||
        state.spatial_dimension > 3) {
        throw std::runtime_error(
            "external particle state has invalid dimension");
    }

    require_key(input, "units", path);
    std::string units;
    input >> units;
    state.unit_system = parse_unit_system(units, path);

    require_key(input, "weighting", path);
    std::string weighting;
    input >> weighting;
    if (weighting != "species_constant") {
        throw std::runtime_error(
            "external particle state supports only weighting species_constant");
    }

    require_key(input, "velocity_staggering", path);
    std::string velocity_staggering;
    input >> velocity_staggering;
    if (velocity_staggering != "time_centered") {
        throw std::runtime_error(
            "external particle state supports only time_centered velocities");
    }

    require_key(input, "particle_count", path);
    input >> state.particle_count;
    if (!input || state.particle_count == 0 ||
        state.particle_count > max_particles) {
        throw std::runtime_error(
            "external particle-state count is zero, invalid, or exceeds the configured limit");
    }

    require_key(input, "records", path);
    for (std::size_t index = 0;
         index < state.particle_count; ++index) {
        require_key(input, "particle", path);
        std::string species;
        ExternalParticleRecord record;
        input >> species >>
            record.position.x >> record.position.y >>
            record.position.z >> record.velocity.x >>
            record.velocity.y >> record.velocity.z;
        if (!input || species.empty() || !finite(record)) {
            throw std::runtime_error(
                "invalid external particle record " +
                std::to_string(index) + " in: " +
                path.string());
        }
        if (state.spatial_dimension < 3 &&
            record.position.z != 0.0) {
            throw std::runtime_error(
                "external particle state has nonzero inactive z position");
        }
        if (state.spatial_dimension == 1 &&
            (record.position.y != 0.0 ||
             record.velocity.y != 0.0 ||
             record.velocity.z != 0.0)) {
            throw std::runtime_error(
                "1D1V external particle state has nonzero inactive components");
        }
        state.species[species].push_back(record);
    }
    require_key(input, "end", path);
    std::string trailing;
    if (input >> trailing) {
        throw std::runtime_error(
            "external particle state contains trailing data: " +
            path.string());
    }
    return state;
}

void validate_external_particle_state(
    const ExternalParticleState& state,
    std::size_t spatial_dimension,
    UnitSystem unit_system,
    const std::vector<ExternalSpeciesExpectation>& expected_species,
    const std::string& context) {
    if (state.version != 1) {
        throw std::invalid_argument(
            context + " external particle-state version is unsupported");
    }
    if (state.spatial_dimension != spatial_dimension) {
        throw std::invalid_argument(
            context + " external particle-state dimension does not match");
    }
    if (state.unit_system != unit_system) {
        throw std::invalid_argument(
            context + " external particle-state unit system does not match");
    }
    std::set<std::string> expected_names;
    std::size_t expected_total = 0;
    for (const auto& expectation : expected_species) {
        if (expectation.name.empty() ||
            !expected_names.insert(expectation.name).second) {
            throw std::invalid_argument(
                context +
                " requires unique non-empty configured species names");
        }
        if (expectation.particle_count >
            std::numeric_limits<std::size_t>::max() -
                expected_total) {
            throw std::overflow_error(
                context + " configured particle count overflow");
        }
        expected_total += expectation.particle_count;
        const auto records =
            state.species.find(expectation.name);
        if (records == state.species.end()) {
            throw std::invalid_argument(
                context + " external particle state is missing species '" +
                expectation.name + "'");
        }
        if (records->second.size() !=
            expectation.particle_count) {
            throw std::invalid_argument(
                context + " external particle count for species '" +
                expectation.name + "' does not match config");
        }
    }
    for (const auto& [name, records] : state.species) {
        (void)records;
        if (!expected_names.contains(name)) {
            throw std::invalid_argument(
                context + " external particle state contains unknown species '" +
                name + "'");
        }
    }
    if (state.particle_count != expected_total) {
        throw std::invalid_argument(
            context + " external total particle count does not match config");
    }
}

ExternalParticleState load_validated_external_particle_state(
    const std::filesystem::path& path,
    std::size_t spatial_dimension,
    UnitSystem unit_system,
    const std::vector<ExternalSpeciesExpectation>& expected_species,
    const std::string& context) {
    std::size_t expected_total = 0;
    for (const auto& species : expected_species) {
        if (species.particle_count >
            std::numeric_limits<std::size_t>::max() -
                expected_total) {
            throw std::overflow_error(
                context + " configured particle count overflow");
        }
        expected_total += species.particle_count;
    }
    if (expected_total == 0) {
        throw std::invalid_argument(
            context + " external state requires configured particles");
    }
    auto state =
        load_external_particle_state(path, expected_total);
    validate_external_particle_state(
        state, spatial_dimension, unit_system,
        expected_species, context);
    return state;
}

} // namespace pic
