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

constexpr const char* particle_state_magic_v1 =
    "AuroraPIC-particle-state-v1";
constexpr const char* particle_state_magic_v2 =
    "AuroraPIC-particle-state-v2";

std::size_t inferred_velocity_dimensions(
    std::size_t version, std::size_t spatial_dimension,
    std::size_t velocity_dimensions) {
    if (version == 1) {
        return spatial_dimension == 1 ? 1 : 3;
    }
    return velocity_dimensions;
}

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
    const auto velocity_dimensions = inferred_velocity_dimensions(
        state.version, state.spatial_dimension,
        state.velocity_dimensions);
    if ((state.version != 1 && state.version != 2) ||
        state.spatial_dimension == 0 ||
        state.spatial_dimension > 3 ||
        (velocity_dimensions != 1 && velocity_dimensions != 3) ||
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
                 record.position.y != 0.0) ||
                (velocity_dimensions == 1 &&
                 (record.velocity.y != 0.0 ||
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

struct ParticleStateScan {
    ExternalParticleStateMetadata metadata;
    std::map<std::string, std::size_t> species_counts;
    std::uint64_t file_order_fingerprint{0};
};

template <typename Consumer>
ParticleStateScan scan_external_particle_state(
    const std::filesystem::path& path,
    std::size_t max_particles,
    Consumer&& consumer) {
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
    ParticleStateScan scan;
    std::string magic;
    input >> magic;
    if (magic == particle_state_magic_v1) {
        scan.metadata.version = 1;
    } else if (magic == particle_state_magic_v2) {
        scan.metadata.version = 2;
    } else {
        throw std::runtime_error(
            "invalid external particle-state magic in: " +
            path.string());
    }
    require_key(input, "dimension", path);
    input >> scan.metadata.spatial_dimension;
    if (!input || scan.metadata.spatial_dimension == 0 ||
        scan.metadata.spatial_dimension > 3) {
        throw std::runtime_error(
            "external particle state has invalid dimension");
    }
    if (scan.metadata.version == 2) {
        require_key(input, "velocity_dimensions", path);
        input >> scan.metadata.velocity_dimensions;
        if (!input ||
            (scan.metadata.velocity_dimensions != 1 &&
             scan.metadata.velocity_dimensions != 3)) {
            throw std::runtime_error(
                "external particle state has invalid velocity dimensions");
        }
    } else {
        scan.metadata.velocity_dimensions =
            inferred_velocity_dimensions(
                1, scan.metadata.spatial_dimension, 0);
    }

    require_key(input, "units", path);
    std::string units;
    input >> units;
    scan.metadata.unit_system =
        parse_unit_system(units, path);

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
    input >> scan.metadata.particle_count;
    if (!input || scan.metadata.particle_count == 0 ||
        scan.metadata.particle_count > max_particles) {
        throw std::runtime_error(
            "external particle-state count is zero, invalid, or exceeds the configured limit");
    }

    std::uint64_t fingerprint =
        14695981039346656037ULL;
    hash_string(fingerprint, magic);
    hash_uint64(fingerprint, scan.metadata.version);
    hash_uint64(
        fingerprint, scan.metadata.spatial_dimension);
    if (scan.metadata.version == 2) {
        hash_uint64(
            fingerprint, scan.metadata.velocity_dimensions);
    }
    hash_string(
        fingerprint, to_string(scan.metadata.unit_system));
    hash_uint64(
        fingerprint, scan.metadata.particle_count);

    require_key(input, "records", path);
    for (std::size_t index = 0;
         index < scan.metadata.particle_count; ++index) {
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
        if ((scan.metadata.spatial_dimension < 3 &&
             record.position.z != 0.0) ||
            (scan.metadata.spatial_dimension == 1 &&
             record.position.y != 0.0)) {
            throw std::runtime_error(
                "external particle state has nonzero inactive position");
        }
        if (scan.metadata.velocity_dimensions == 1 &&
            (record.velocity.y != 0.0 ||
             record.velocity.z != 0.0)) {
            throw std::runtime_error(
                "external particle state has nonzero inactive velocity");
        }
        auto& species_count =
            scan.species_counts[species];
        const auto species_index = species_count++;
        hash_string(fingerprint, species);
        hash_double(fingerprint, record.position.x);
        hash_double(fingerprint, record.position.y);
        hash_double(fingerprint, record.position.z);
        hash_double(fingerprint, record.velocity.x);
        hash_double(fingerprint, record.velocity.y);
        hash_double(fingerprint, record.velocity.z);
        consumer(species, species_index, record);
    }
    require_key(input, "end", path);
    std::string trailing;
    if (input >> trailing) {
        throw std::runtime_error(
            "external particle state contains trailing data: " +
            path.string());
    }
    scan.file_order_fingerprint = fingerprint;
    return scan;
}

} // namespace

ExternalParticleState load_external_particle_state(
    const std::filesystem::path& path,
    std::size_t max_particles) {
    ExternalParticleState state;
    const auto scan = scan_external_particle_state(
        path, max_particles,
        [&](const std::string& species,
            std::size_t,
            const ExternalParticleRecord& record) {
            state.species[species].push_back(record);
        });
    state.version = scan.metadata.version;
    state.spatial_dimension =
        scan.metadata.spatial_dimension;
    state.velocity_dimensions =
        scan.metadata.velocity_dimensions;
    state.unit_system = scan.metadata.unit_system;
    state.particle_count =
        scan.metadata.particle_count;
    state.signature =
        external_particle_state_signature(state);
    return state;
}

void validate_external_particle_state(
    const ExternalParticleState& state,
    std::size_t spatial_dimension,
    std::size_t velocity_dimensions,
    UnitSystem unit_system,
    const std::vector<ExternalSpeciesExpectation>& expected_species,
    const std::string& context) {
    validate_structure(state, context);
    if (state.version != 1 && state.version != 2) {
        throw std::invalid_argument(
            context + " external particle-state version is unsupported");
    }
    if (state.spatial_dimension != spatial_dimension) {
        throw std::invalid_argument(
            context + " external particle-state dimension does not match");
    }
    if (inferred_velocity_dimensions(
            state.version, state.spatial_dimension,
            state.velocity_dimensions) != velocity_dimensions) {
        throw std::invalid_argument(
            context +
            " external particle-state velocity dimensions do not match");
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
    std::size_t velocity_dimensions,
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
        state, spatial_dimension, velocity_dimensions, unit_system,
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

ExternalParticleStateMetadata
load_validated_external_particle_state_bounded(
    const std::filesystem::path& path,
    std::size_t spatial_dimension,
    std::size_t velocity_dimensions,
    UnitSystem unit_system,
    const std::vector<ExternalSpeciesExpectation>& expected_species,
    const std::string& context,
    const ExternalParticleRecordConsumer& consumer,
    std::optional<std::uint64_t> expected_signature) {
    if (!consumer) {
        throw std::invalid_argument(
            context +
            " external particle-state consumer is empty");
    }

    struct ExpectedSpecies {
        std::size_t configured_index{0};
        std::size_t particle_count{0};
    };
    std::map<std::string, ExpectedSpecies> canonical_species;
    std::size_t expected_total = 0;
    for (std::size_t index = 0;
         index < expected_species.size(); ++index) {
        const auto& species = expected_species[index];
        if (species.name.empty() ||
            !canonical_species.emplace(
                species.name,
                ExpectedSpecies{
                    index, species.particle_count}).second) {
            throw std::invalid_argument(
                context +
                " requires unique non-empty configured species names");
        }
        if (species.particle_count >
            std::numeric_limits<std::size_t>::max() -
                expected_total) {
            throw std::overflow_error(
                context +
                " configured particle count overflow");
        }
        expected_total += species.particle_count;
    }
    if (expected_total == 0) {
        throw std::invalid_argument(
            context +
            " external state requires configured particles");
    }

    const auto baseline = scan_external_particle_state(
        path, expected_total,
        [](const std::string&, std::size_t,
           const ExternalParticleRecord&) {});
    if ((baseline.metadata.version != 1 &&
         baseline.metadata.version != 2) ||
        baseline.metadata.spatial_dimension != spatial_dimension) {
        throw std::invalid_argument(
            context +
            " external particle-state dimension or version does not match");
    }
    if (baseline.metadata.velocity_dimensions != velocity_dimensions) {
        throw std::invalid_argument(
            context +
            " external particle-state velocity dimensions do not match");
    }
    if (baseline.metadata.unit_system != unit_system) {
        throw std::invalid_argument(
            context +
            " external particle-state unit system does not match");
    }
    if (baseline.metadata.particle_count != expected_total) {
        throw std::invalid_argument(
            context +
            " external total particle count does not match config");
    }
    for (const auto& [name, count] : baseline.species_counts) {
        const auto expected = canonical_species.find(name);
        if (expected == canonical_species.end()) {
            throw std::invalid_argument(
                context +
                " external particle state contains unknown species '" +
                name + "'");
        }
        if (count != expected->second.particle_count) {
            throw std::invalid_argument(
                context + " external particle count for species '" +
                name + "' does not match config");
        }
    }
    for (const auto& [name, expectation] : canonical_species) {
        (void)expectation;
        if (!baseline.species_counts.contains(name)) {
            throw std::invalid_argument(
                context + " external particle state is missing species '" +
                name + "'");
        }
    }

    std::uint64_t signature = 14695981039346656037ULL;
    hash_string(
        signature,
        baseline.metadata.version == 1
            ? particle_state_magic_v1
            : particle_state_magic_v2);
    hash_uint64(signature, baseline.metadata.version);
    hash_uint64(signature, spatial_dimension);
    if (baseline.metadata.version == 2) {
        hash_uint64(signature, velocity_dimensions);
    }
    hash_string(signature, to_string(unit_system));
    hash_uint64(signature, expected_total);
    hash_uint64(signature, canonical_species.size());

    for (const auto& [species_name, expectation] :
         canonical_species) {
        hash_string(signature, species_name);
        hash_uint64(
            signature, expectation.particle_count);
        const auto scan = scan_external_particle_state(
            path, expected_total,
            [&](const std::string& record_species,
                std::size_t,
                const ExternalParticleRecord& record) {
                if (record_species != species_name) return;
                hash_double(signature, record.position.x);
                hash_double(signature, record.position.y);
                hash_double(signature, record.position.z);
                hash_double(signature, record.velocity.x);
                hash_double(signature, record.velocity.y);
                hash_double(signature, record.velocity.z);
            });
        if (scan.file_order_fingerprint !=
                baseline.file_order_fingerprint ||
            scan.species_counts !=
                baseline.species_counts) {
            throw std::runtime_error(
                context +
                " external particle state changed while it was being read");
        }
    }

    if (expected_signature &&
        signature != *expected_signature) {
        throw std::invalid_argument(
            context +
            " external particle-state signature mismatch: expected " +
            std::to_string(*expected_signature) + ", got " +
            std::to_string(signature));
    }

    const auto final_scan = scan_external_particle_state(
        path, expected_total,
        [&](const std::string& species_name,
            std::size_t record_index,
            const ExternalParticleRecord& record) {
            const auto expected =
                canonical_species.find(species_name);
            if (expected == canonical_species.end()) {
                throw std::runtime_error(
                    context +
                    " external particle state changed while it was being read");
            }
            consumer(
                expected->second.configured_index,
                record_index, record);
        });
    if (final_scan.file_order_fingerprint !=
            baseline.file_order_fingerprint ||
        final_scan.species_counts !=
            baseline.species_counts) {
        throw std::runtime_error(
            context +
            " external particle state changed while it was being read");
    }

    auto metadata = baseline.metadata;
    metadata.signature = signature;
    return metadata;
}

std::uint64_t external_particle_state_signature(
    const ExternalParticleState& state) {
    validate_structure(
        state, "external particle state");
    std::uint64_t hash = 14695981039346656037ULL;
    hash_string(
        hash,
        state.version == 1
            ? particle_state_magic_v1
            : particle_state_magic_v2);
    hash_uint64(hash, state.version);
    hash_uint64(hash, state.spatial_dimension);
    if (state.version == 2) {
        hash_uint64(hash, state.velocity_dimensions);
    }
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
        inferred_velocity_dimensions(
            state.version, state.spatial_dimension,
            state.velocity_dimensions),
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
    output << (state.version == 1
                   ? particle_state_magic_v1
                   : particle_state_magic_v2) << '\n'
           << "dimension " << state.spatial_dimension << '\n';
    if (state.version == 2) {
        output << "velocity_dimensions "
               << state.velocity_dimensions << '\n';
    }
    output
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
    if ((metadata.version != 1 && metadata.version != 2) ||
        metadata.spatial_dimension == 0 ||
        metadata.spatial_dimension > 3 ||
        (metadata.velocity_dimensions != 1 &&
         metadata.velocity_dimensions != 3) ||
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
           << "velocity_dimensions "
           << metadata.velocity_dimensions << '\n'
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
