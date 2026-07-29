#include "pic/ParticleState.hpp"

#include <bit>
#include <cmath>
#include <fstream>
#include <iomanip>
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

void validate_structure(
    const ExternalParticleState& state,
    const std::string& context) {
    if (state.version != 1 ||
        state.spatial_dimension == 0 ||
        state.spatial_dimension > 3 ||
        state.particle_count == 0) {
        throw std::invalid_argument(
            context + " has invalid particle-state metadata");
    }
    std::size_t count = 0;
    for (const auto& [species, records] : state.species) {
        if (species.empty() || records.empty()) {
            throw std::invalid_argument(
                context + " has an empty species name or population");
        }
        if (records.size() >
            std::numeric_limits<std::size_t>::max() - count) {
            throw std::overflow_error(
                context + " particle count overflow");
        }
        count += records.size();
        for (const auto& record : records) {
            if (!finite(record) ||
                (state.spatial_dimension < 3 &&
                 record.position.z != 0.0) ||
                (state.spatial_dimension == 1 &&
                 (record.position.y != 0.0 ||
                  record.velocity.y != 0.0 ||
                  record.velocity.z != 0.0))) {
                throw std::invalid_argument(
                    context + " has invalid or active out-of-dimension components");
            }
        }
    }
    if (count != state.particle_count) {
        throw std::invalid_argument(
            context + " particle count does not match its records");
    }
}

void hash_uint64(std::uint64_t& hash, std::uint64_t value) {
    constexpr std::uint64_t prime = 1099511628211ULL;
    for (unsigned byte = 0; byte < 8; ++byte) {
        hash ^= static_cast<unsigned char>(
            value >> (byte * 8));
        hash *= prime;
    }
}

void hash_string(std::uint64_t& hash, const std::string& value) {
    hash_uint64(hash, value.size());
    for (const unsigned char character : value) {
        hash ^= character;
        hash *= 1099511628211ULL;
    }
}

void hash_double(std::uint64_t& hash, double value) {
    hash_uint64(
        hash, std::bit_cast<std::uint64_t>(value));
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
    state.signature =
        external_particle_state_signature(state);
    return state;
}

void validate_external_particle_state(
    const ExternalParticleState& state,
    std::size_t spatial_dimension,
    UnitSystem unit_system,
    const std::vector<ExternalSpeciesExpectation>& expected_species,
    const std::string& context) {
    validate_structure(state, context);
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
    const std::string& context,
    std::optional<std::uint64_t> expected_signature) {
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
    if (expected_signature &&
        state.signature != *expected_signature) {
        throw std::invalid_argument(
            context +
            " external particle-state signature mismatch: expected " +
            std::to_string(*expected_signature) + ", got " +
            std::to_string(state.signature));
    }
    return state;
}

std::uint64_t external_particle_state_signature(
    const ExternalParticleState& state) {
    validate_structure(
        state, "external particle state");
    std::uint64_t hash = 14695981039346656037ULL;
    hash_string(hash, particle_state_magic);
    hash_uint64(hash, state.version);
    hash_uint64(hash, state.spatial_dimension);
    hash_string(hash, to_string(state.unit_system));
    hash_uint64(hash, state.particle_count);
    hash_uint64(hash, state.species.size());
    for (const auto& [species, records] : state.species) {
        hash_string(hash, species);
        hash_uint64(hash, records.size());
        for (const auto& record : records) {
            hash_double(hash, record.position.x);
            hash_double(hash, record.position.y);
            hash_double(hash, record.position.z);
            hash_double(hash, record.velocity.x);
            hash_double(hash, record.velocity.y);
            hash_double(hash, record.velocity.z);
        }
    }
    return hash;
}

ExternalParticleStateMetadata external_particle_state_metadata(
    const ExternalParticleState& state) {
    const auto realized_signature =
        external_particle_state_signature(state);
    if (state.signature != realized_signature) {
        throw std::invalid_argument(
            "external particle-state metadata received stale signature data");
    }
    return {
        state.version,
        state.spatial_dimension,
        state.unit_system,
        state.particle_count,
        realized_signature};
}

void write_external_particle_state(
    const std::filesystem::path& path,
    const ExternalParticleState& state,
    bool overwrite) {
    validate_structure(
        state, "external particle state writer");
    if (!overwrite && std::filesystem::exists(path)) {
        throw std::runtime_error(
            "external particle-state output already exists: " +
            path.string());
    }
    const auto parent = path.parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }
    std::ofstream output(path, std::ios::trunc);
    if (!output) {
        throw std::runtime_error(
            "cannot write external particle state: " +
            path.string());
    }
    output << particle_state_magic << '\n'
           << "dimension " << state.spatial_dimension << '\n'
           << "units " << to_string(state.unit_system) << '\n'
           << "weighting species_constant\n"
           << "velocity_staggering time_centered\n"
           << "particle_count " << state.particle_count << '\n'
           << "records\n"
           << std::setprecision(17);
    for (const auto& [species, records] : state.species) {
        for (const auto& record : records) {
            output << "particle " << species << ' '
                   << record.position.x << ' '
                   << record.position.y << ' '
                   << record.position.z << ' '
                   << record.velocity.x << ' '
                   << record.velocity.y << ' '
                   << record.velocity.z << '\n';
        }
    }
    output << "end\n";
    if (!output) {
        throw std::runtime_error(
            "failed while writing external particle state: " +
            path.string());
    }
}

void write_external_particle_state_metadata(
    const std::filesystem::path& path,
    const std::filesystem::path& source_path,
    const ExternalParticleStateMetadata& metadata,
    std::optional<std::uint64_t> expected_signature) {
    if (metadata.version != 1 ||
        metadata.spatial_dimension == 0 ||
        metadata.spatial_dimension > 3 ||
        metadata.particle_count == 0) {
        throw std::invalid_argument(
            "external particle-state metadata is invalid");
    }
    const auto parent = path.parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error(
            "cannot write external particle-state metadata: " +
            path.string());
    }
    output << "format aurorapic_particle_state\n"
           << "version " << metadata.version << '\n'
           << "source_path "
           << std::quoted(source_path.string()) << '\n'
           << "dimension " << metadata.spatial_dimension << '\n'
           << "units " << to_string(metadata.unit_system) << '\n'
           << "particle_count " << metadata.particle_count << '\n'
           << "signature " << metadata.signature << '\n'
           << "expected_signature ";
    if (expected_signature) output << *expected_signature;
    else output << "none";
    output << '\n';
    if (!output) {
        throw std::runtime_error(
            "failed while writing external particle-state metadata: " +
            path.string());
    }
}

} // namespace pic
