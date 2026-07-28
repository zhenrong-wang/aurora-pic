#include "pic/Swarm.hpp"

#include "pic/Collision.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace pic {
namespace {

constexpr double electron_mass_kg = 9.1093837139e-31;
constexpr double elementary_charge_c = 1.602176634e-19;
constexpr double ev_to_j = elementary_charge_c;
constexpr double townsend_v_m2 = 1.0e-21;

std::string trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

template <typename T>
T parse_number(const std::string& value, const std::string& context) {
    if constexpr (std::is_unsigned_v<T>) {
        if (!value.empty() && value.front() == '-') {
            throw std::runtime_error(
                context + " must not be negative");
        }
    }
    std::istringstream input(value);
    T result{};
    if (!(input >> result)) {
        throw std::runtime_error(context + " must be numeric");
    }
    std::string trailing;
    if (input >> trailing) {
        throw std::runtime_error(context + " has trailing content");
    }
    return result;
}

template <typename T>
T required_number(
    const std::map<std::string, std::string>& values,
    const std::string& key,
    const std::string& context) {
    const auto found = values.find(key);
    if (found == values.end() || found->second.empty()) {
        throw std::runtime_error(context + " requires '" + key + "'");
    }
    return parse_number<T>(
        found->second, context + " key '" + key + "'");
}

std::vector<double> parse_fields(
    const std::string& value,
    const std::string& context) {
    if (value.empty() || value.front() == ',' || value.back() == ',') {
        throw std::runtime_error(
            context + " reduced_fields_td contains an empty value");
    }
    std::vector<double> result;
    std::istringstream input(value);
    std::string field;
    while (std::getline(input, field, ',')) {
        field = trim(field);
        if (field.empty()) {
            throw std::runtime_error(
                context + " reduced_fields_td contains an empty value");
        }
        result.push_back(parse_number<double>(
            field, context + " reduced_fields_td value"));
    }
    if (result.empty()) {
        throw std::runtime_error(
            context + " requires at least one reduced field");
    }
    return result;
}

double kinetic_energy_ev(const Vec3& velocity) {
    const double speed_squared =
        velocity.x * velocity.x +
        velocity.y * velocity.y +
        velocity.z * velocity.z;
    return 0.5 * electron_mass_kg * speed_squared / ev_to_j;
}

std::pair<double, double> block_mean_and_error(
    const std::vector<double>& samples,
    std::size_t blocks) {
    if (samples.empty() || blocks == 0 || blocks > samples.size()) {
        throw std::invalid_argument(
            "swarm block statistics require at least one sample per block");
    }
    std::vector<double> means;
    means.reserve(blocks);
    for (std::size_t block = 0; block < blocks; ++block) {
        const std::size_t begin = block * samples.size() / blocks;
        const std::size_t end = (block + 1) * samples.size() / blocks;
        const double total = std::accumulate(
            samples.begin() + static_cast<std::ptrdiff_t>(begin),
            samples.begin() + static_cast<std::ptrdiff_t>(end), 0.0);
        means.push_back(total / static_cast<double>(end - begin));
    }
    const double mean = std::accumulate(
        means.begin(), means.end(), 0.0) /
        static_cast<double>(means.size());
    if (means.size() == 1) return {mean, 0.0};
    double squared_deviation = 0.0;
    for (const double value : means) {
        const double delta = value - mean;
        squared_deviation += delta * delta;
    }
    const double standard_error = std::sqrt(
        squared_deviation /
        (static_cast<double>(means.size()) *
         static_cast<double>(means.size() - 1)));
    return {mean, standard_error};
}

std::string csv_cell(const std::string& value) {
    if (value.find_first_of(",\"\r\n") == std::string::npos) {
        return value;
    }
    std::string result{"\""};
    for (const char character : value) {
        if (character == '"') result.push_back('"');
        result.push_back(character);
    }
    result.push_back('"');
    return result;
}

void validate_config(
    const SwarmBenchmarkConfig& config,
    const GasDataset& dataset) {
    if (dataset.format_version != 2 ||
        dataset.unit_system != UnitSystem::SI) {
        throw std::runtime_error(
            "swarm benchmarks require a version-2 SI gas dataset");
    }
    const auto positive_finite = [](double value) {
        return std::isfinite(value) && value > 0.0;
    };
    if (!positive_finite(config.neutral_density)) {
        throw std::runtime_error(
            "swarm neutral_density must be positive and finite");
    }
    if (!positive_finite(config.max_frequency)) {
        throw std::runtime_error(
            "swarm max_frequency must be positive and finite");
    }
    if (!positive_finite(config.timestep)) {
        throw std::runtime_error(
            "swarm timestep must be positive and finite");
    }
    if (config.max_frequency * config.timestep > 0.1) {
        throw std::runtime_error(
            "swarm max_frequency * timestep must not exceed 0.1");
    }
    if (config.steps == 0 ||
        config.burn_in_steps >= config.steps) {
        throw std::runtime_error(
            "swarm steps must exceed burn_in_steps");
    }
    if (config.particles == 0) {
        throw std::runtime_error(
            "swarm particles must be positive");
    }
    constexpr std::size_t particle_limit = 10000000;
    if (config.particles > particle_limit) {
        throw std::runtime_error(
            "swarm particles exceeds the safety limit of 10000000");
    }
    if (config.work_item_limit == 0) {
        throw std::runtime_error(
            "swarm work_item_limit must be positive");
    }
    const auto exceeds_product_limit = [](
        std::uint64_t first,
        std::uint64_t second,
        std::uint64_t limit) {
        return first != 0 && second > limit / first;
    };
    const auto fields = static_cast<std::uint64_t>(
        config.reduced_fields_td.size());
    const auto particles =
        static_cast<std::uint64_t>(config.particles);
    const auto steps = static_cast<std::uint64_t>(config.steps);
    if (exceeds_product_limit(fields, particles,
                              config.work_item_limit) ||
        exceeds_product_limit(
            fields * particles, steps,
            config.work_item_limit)) {
        throw std::runtime_error(
            "swarm scan exceeds work_item_limit; raise it explicitly "
            "only on an appropriate compute host");
    }
    const std::size_t sampling_steps =
        config.steps - config.burn_in_steps;
    if (config.uncertainty_blocks == 0 ||
        config.uncertainty_blocks > sampling_steps) {
        throw std::runtime_error(
            "swarm uncertainty_blocks must be between 1 and the "
            "number of sampling steps");
    }
    if (sampling_steps % config.uncertainty_blocks != 0) {
        throw std::runtime_error(
            "swarm sampling steps must be divisible by "
            "uncertainty_blocks");
    }
    if (!std::isfinite(config.initial_mean_energy_ev) ||
        config.initial_mean_energy_ev < 0.0) {
        throw std::runtime_error(
            "swarm initial_mean_energy_ev must be finite and non-negative");
    }
    if (!positive_finite(config.max_energy_ev)) {
        throw std::runtime_error(
            "swarm max_energy_ev must be positive and finite");
    }
    if (config.reduced_fields_td.empty()) {
        throw std::runtime_error(
            "swarm requires at least one reduced field");
    }
    for (const double field : config.reduced_fields_td) {
        if (!std::isfinite(field) || field < 0.0) {
            throw std::runtime_error(
                "swarm reduced fields must be finite and non-negative");
        }
        const double electric_field =
            field * townsend_v_m2 * config.neutral_density;
        if (!std::isfinite(electric_field)) {
            throw std::runtime_error(
                "swarm reduced field and neutral density overflow "
                "the electric field");
        }
    }
    for (std::size_t first = 0;
         first < config.reduced_fields_td.size(); ++first) {
        for (std::size_t second = first + 1;
             second < config.reduced_fields_td.size(); ++second) {
            if (config.reduced_fields_td[first] ==
                config.reduced_fields_td[second]) {
                throw std::runtime_error(
                    "swarm reduced fields must not contain duplicates");
            }
        }
    }
    if (dataset.channels.empty()) {
        throw std::runtime_error(
            "swarm gas dataset has no collision channels");
    }
    bool has_elastic = false;
    for (const auto& channel : dataset.channels) {
        has_elastic =
            has_elastic ||
            channel.process == CollisionProcessKind::Elastic;
        if (channel.process == CollisionProcessKind::ChargeExchange) {
            throw std::runtime_error(
                "electron swarm benchmark does not support "
                "charge-exchange channels");
        }
        const CrossSectionTable table(
            channel.cross_section_file,
            channel.energy_scale,
            channel.cross_section_scale);
        const double table_max_ev =
            table.energies().back() / ev_to_j;
        if (config.max_energy_ev > table_max_ev) {
            throw std::runtime_error(
                "swarm max_energy_ev exceeds channel '" +
                channel.name + "' table coverage");
        }
        if (channel.angular_scattering ==
            AngularScatteringKind::HenyeyGreenstein) {
            const MeanCosineTable angular_table(
                channel.mean_cosine_file,
                channel.mean_cosine_energy_scale);
            const double angular_max_ev =
                angular_table.energies().back() / ev_to_j;
            if (config.max_energy_ev > angular_max_ev) {
                throw std::runtime_error(
                    "swarm max_energy_ev exceeds channel '" +
                    channel.name + "' angular table coverage");
            }
        }
    }
    if (!has_elastic) {
        throw std::runtime_error(
            "electron swarm benchmark requires an elastic channel");
    }
}

CollisionConfig collision_config(
    const SwarmBenchmarkConfig& config,
    const GasDataset& dataset) {
    CollisionConfig result;
    result.enabled = true;
    result.model = CollisionModelKind::NullCollision;
    result.neutral_density = config.neutral_density;
    result.species = "electrons";
    result.max_frequency = config.max_frequency;
    result.gas_name = dataset.gas_name;
    result.neutral_mass = dataset.neutral_mass;
    result.data_provenance = dataset.data_provenance;
    result.gas_data_file = config.gas_data_file;
    result.gas_data_version = dataset.format_version;
    result.gas_data_units = dataset.unit_system;
    result.dataset_id = dataset.dataset_id;
    result.dataset_version = dataset.dataset_version;
    result.citation = dataset.citation;
    result.retrieved = dataset.retrieved;
    result.license = dataset.license;
    result.channels = dataset.channels;
    for (auto& channel : result.channels) {
        if (channel.process == CollisionProcessKind::Ionization) {
            channel.secondary_species = "fixed_population_electron";
            channel.ion_species = "diagnostic_ion";
        }
    }
    return result;
}

SwarmBenchmarkResult run_field(
    const SwarmBenchmarkConfig& config,
    const GasDataset& dataset,
    double reduced_field_td,
    std::uint64_t seed) {
    const CollisionConfig collisions =
        collision_config(config, dataset);
    const NullCollisionModel model(collisions, electron_mass_kg);
    std::mt19937_64 rng(seed);
    const double component_stddev = std::sqrt(
        2.0 * config.initial_mean_energy_ev * ev_to_j /
        (3.0 * electron_mass_kg));
    std::normal_distribution<double> normal(0.0, component_stddev);
    std::vector<Vec3> velocities(config.particles);
    for (auto& velocity : velocities) {
        velocity = {normal(rng), normal(rng), normal(rng)};
    }
    std::vector<Vec3> positions(config.particles);
    std::vector<Vec3> sampling_start_positions(config.particles);

    const double electric_field =
        reduced_field_td * townsend_v_m2 * config.neutral_density;
    const double acceleration =
        -elementary_charge_c * electric_field / electron_mass_kg;
    const double half_kick = 0.5 * acceleration * config.timestep;
    const std::size_t sampling_steps =
        config.steps - config.burn_in_steps;
    std::vector<double> velocity_samples;
    std::vector<double> energy_samples;
    velocity_samples.reserve(sampling_steps);
    energy_samples.reserve(sampling_steps);

    SwarmBenchmarkResult result;
    result.collision_model_signature = model.signature();
    result.reduced_field_td = reduced_field_td;
    result.electric_field_v_m = electric_field;
    result.channels.reserve(model.channel_names().size());
    for (const auto& name : model.channel_names()) {
        result.channels.push_back({name});
    }

    for (std::size_t step = 0; step < config.steps; ++step) {
        if (step == config.burn_in_steps) {
            sampling_start_positions = positions;
        }
        for (std::size_t particle = 0;
             particle < velocities.size(); ++particle) {
            auto& velocity = velocities[particle];
            velocity.x += half_kick;
            const double pre_collision_energy_ev =
                kinetic_energy_ev(velocity);
            result.maximum_observed_energy_ev = std::max(
                result.maximum_observed_energy_ev,
                pre_collision_energy_ev);
            if (pre_collision_energy_ev > config.max_energy_ev) {
                throw std::runtime_error(
                    "swarm particle energy exceeded max_energy_ev before "
                    "collision lookup at E/N=" +
                    std::to_string(reduced_field_td) + " Td");
            }
            positions[particle].x +=
                velocity.x * config.timestep;
            positions[particle].y +=
                velocity.y * config.timestep;
            positions[particle].z +=
                velocity.z * config.timestep;
            const auto statistics =
                model.collide(velocity, config.timestep, rng);
            velocity.x += half_kick;

            const double energy_ev = kinetic_energy_ev(velocity);
            result.maximum_observed_energy_ev = std::max(
                result.maximum_observed_energy_ev, energy_ev);
            if (energy_ev > config.max_energy_ev) {
                throw std::runtime_error(
                    "swarm particle energy exceeded max_energy_ev at "
                    "E/N=" + std::to_string(reduced_field_td) + " Td");
            }
            if (step < config.burn_in_steps) continue;
            result.collision_candidates += statistics.candidates;
            result.null_collisions += statistics.null_collisions;
            for (std::size_t channel = 0;
                 channel < result.channels.size(); ++channel) {
                result.channels[channel].collisions +=
                    statistics.channel_collisions[channel];
            }
        }
        if (step < config.burn_in_steps) continue;
        double velocity_sum = 0.0;
        double energy_sum = 0.0;
        for (const auto& velocity : velocities) {
            velocity_sum += velocity.x;
            energy_sum += kinetic_energy_ev(velocity);
        }
        velocity_samples.push_back(
            velocity_sum / static_cast<double>(velocities.size()));
        energy_samples.push_back(
            energy_sum / static_cast<double>(velocities.size()));
    }

    const auto velocity_statistics = block_mean_and_error(
        velocity_samples, config.uncertainty_blocks);
    const auto energy_statistics = block_mean_and_error(
        energy_samples, config.uncertainty_blocks);
    result.mean_velocity_x_m_s = velocity_statistics.first;
    result.mean_velocity_x_standard_error_m_s =
        velocity_statistics.second;
    result.electron_drift_velocity_m_s =
        -result.mean_velocity_x_m_s;
    result.reduced_mobility_1_v_m_s =
        electric_field == 0.0
            ? 0.0
            : result.electron_drift_velocity_m_s *
                  config.neutral_density / electric_field;
    result.mean_energy_ev = energy_statistics.first;
    result.mean_energy_standard_error_ev = energy_statistics.second;

    double mean_dx = 0.0;
    double mean_dy = 0.0;
    double mean_dz = 0.0;
    for (std::size_t particle = 0;
         particle < positions.size(); ++particle) {
        mean_dx += positions[particle].x -
                   sampling_start_positions[particle].x;
        mean_dy += positions[particle].y -
                   sampling_start_positions[particle].y;
        mean_dz += positions[particle].z -
                   sampling_start_positions[particle].z;
    }
    const double particle_count =
        static_cast<double>(positions.size());
    mean_dx /= particle_count;
    mean_dy /= particle_count;
    mean_dz /= particle_count;
    double variance_x = 0.0;
    double variance_yz = 0.0;
    for (std::size_t particle = 0;
         particle < positions.size(); ++particle) {
        const double dx =
            positions[particle].x -
            sampling_start_positions[particle].x - mean_dx;
        const double dy =
            positions[particle].y -
            sampling_start_positions[particle].y - mean_dy;
        const double dz =
            positions[particle].z -
            sampling_start_positions[particle].z - mean_dz;
        variance_x += dx * dx;
        variance_yz += dy * dy + dz * dz;
    }
    variance_x /= particle_count;
    variance_yz /= particle_count;
    const double sampling_time =
        static_cast<double>(sampling_steps) * config.timestep;
    result.longitudinal_diffusion_m2_s =
        variance_x / (2.0 * sampling_time);
    result.transverse_diffusion_m2_s =
        variance_yz / (4.0 * sampling_time);

    const double exposure =
        particle_count * sampling_time;
    for (auto& channel : result.channels) {
        channel.rate_per_electron_s =
            static_cast<double>(channel.collisions) / exposure;
        channel.poisson_standard_error_s =
            std::sqrt(static_cast<double>(channel.collisions)) /
            exposure;
    }
    return result;
}

} // namespace

SwarmBenchmarkConfig load_swarm_benchmark_config(
    const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error(
            "cannot open swarm config: " + path.string());
    }
    static const std::set<std::string> allowed{
        "swarm_config_version",
        "gas_data_file",
        "neutral_density",
        "reduced_fields_td",
        "max_frequency",
        "timestep",
        "steps",
        "burn_in_steps",
        "particles",
        "uncertainty_blocks",
        "work_item_limit",
        "initial_mean_energy_ev",
        "max_energy_ev",
        "seed",
        "output_file",
    };
    std::map<std::string, std::string> values;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        const auto comment = line.find_first_of("#;");
        if (comment != std::string::npos) line.resize(comment);
        line = trim(line);
        if (line.empty()) continue;
        const auto equals = line.find('=');
        if (equals == std::string::npos) {
            throw std::runtime_error(
                path.string() + ":" + std::to_string(line_number) +
                ": expected key = value");
        }
        const std::string key = trim(line.substr(0, equals));
        const std::string value = trim(line.substr(equals + 1));
        if (!allowed.contains(key)) {
            throw std::runtime_error(
                path.string() + ":" + std::to_string(line_number) +
                ": unknown key '" + key + "'");
        }
        if (!values.emplace(key, value).second) {
            throw std::runtime_error(
                path.string() + ":" + std::to_string(line_number) +
                ": duplicate key '" + key + "'");
        }
    }

    const std::string context =
        "swarm config '" + path.string() + "'";
    const std::size_t version =
        required_number<std::size_t>(
            values, "swarm_config_version", context);
    if (version != 1) {
        throw std::runtime_error(
            context + " requires swarm_config_version = 1");
    }
    const auto required_string = [&](const std::string& key) {
        const auto found = values.find(key);
        if (found == values.end() || found->second.empty()) {
            throw std::runtime_error(
                context + " requires '" + key + "'");
        }
        return found->second;
    };

    SwarmBenchmarkConfig result;
    const auto base = path.parent_path();
    result.gas_data_file =
        (base / required_string("gas_data_file")).lexically_normal();
    result.neutral_density =
        required_number<double>(
            values, "neutral_density", context);
    result.reduced_fields_td = parse_fields(
        required_string("reduced_fields_td"), context);
    result.max_frequency =
        required_number<double>(
            values, "max_frequency", context);
    result.timestep =
        required_number<double>(values, "timestep", context);
    result.steps =
        required_number<std::size_t>(values, "steps", context);
    result.burn_in_steps =
        required_number<std::size_t>(
            values, "burn_in_steps", context);
    result.particles =
        required_number<std::size_t>(
            values, "particles", context);
    if (values.contains("uncertainty_blocks")) {
        result.uncertainty_blocks =
            parse_number<std::size_t>(
                values.at("uncertainty_blocks"),
                context + " key 'uncertainty_blocks'");
    }
    if (values.contains("work_item_limit")) {
        result.work_item_limit =
            parse_number<std::uint64_t>(
                values.at("work_item_limit"),
                context + " key 'work_item_limit'");
    }
    if (values.contains("initial_mean_energy_ev")) {
        result.initial_mean_energy_ev =
            parse_number<double>(
                values.at("initial_mean_energy_ev"),
                context + " key 'initial_mean_energy_ev'");
    }
    result.max_energy_ev =
        required_number<double>(
            values, "max_energy_ev", context);
    if (values.contains("seed")) {
        result.seed = parse_number<std::uint64_t>(
            values.at("seed"), context + " key 'seed'");
    }
    if (values.contains("output_file")) {
        result.output_file =
            (base / required_string("output_file")).lexically_normal();
    } else {
        result.output_file = base / result.output_file;
    }
    const auto absolute_config =
        std::filesystem::absolute(path).lexically_normal();
    const auto absolute_gas =
        std::filesystem::absolute(
            result.gas_data_file).lexically_normal();
    const auto absolute_output =
        std::filesystem::absolute(
            result.output_file).lexically_normal();
    if (absolute_output == absolute_config ||
        absolute_output == absolute_gas) {
        throw std::runtime_error(
            context + " output_file must not overwrite an input file");
    }
    return result;
}

std::vector<SwarmBenchmarkResult> run_swarm_benchmark(
    const SwarmBenchmarkConfig& config) {
    const GasDataset dataset =
        load_gas_dataset(config.gas_data_file);
    validate_config(config, dataset);
    std::vector<SwarmBenchmarkResult> results;
    results.reserve(config.reduced_fields_td.size());
    for (std::size_t field = 0;
         field < config.reduced_fields_td.size(); ++field) {
        results.push_back(run_field(
            config, dataset, config.reduced_fields_td[field],
            config.seed + static_cast<std::uint64_t>(field)));
    }
    return results;
}

void write_swarm_benchmark_csv(
    const std::filesystem::path& path,
    const SwarmBenchmarkConfig& config,
    const GasDataset& dataset,
    const std::vector<SwarmBenchmarkResult>& results) {
    if (results.empty()) {
        throw std::invalid_argument(
            "cannot write an empty swarm result set");
    }
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error(
            "cannot write swarm output: " + path.string());
    }
    output << std::setprecision(17);
    output
        << "dataset_id,dataset_version,gas,retrieved,provenance,"
        << "citation,license,gas_data_file,population_model,"
        << "neutral_density_m3,timestep_s,steps,burn_in_steps,"
        << "particles,work_item_limit,seed,reduced_field_td,"
        << "collision_model_signature,"
        << "electric_field_v_m,mean_velocity_x_m_s,"
        << "mean_velocity_x_standard_error_m_s,"
        << "electron_drift_velocity_m_s,"
        << "reduced_mobility_1_v_m_s,mean_energy_ev,"
        << "mean_energy_standard_error_ev,"
        << "longitudinal_diffusion_m2_s,"
        << "transverse_diffusion_m2_s,maximum_observed_energy_ev,"
        << "collision_candidates,null_collisions";
    for (const auto& channel : results.front().channels) {
        output << ',' << channel.name << "_collisions"
               << ',' << channel.name << "_rate_per_electron_s"
               << ',' << channel.name << "_rate_standard_error_s";
    }
    output << '\n';
    for (std::size_t row = 0; row < results.size(); ++row) {
        const auto& result = results[row];
        if (result.channels.size() != results.front().channels.size()) {
            throw std::invalid_argument(
                "swarm result channel counts are inconsistent");
        }
        output
            << csv_cell(dataset.dataset_id) << ','
            << csv_cell(dataset.dataset_version) << ','
            << csv_cell(dataset.gas_name) << ','
            << csv_cell(dataset.retrieved) << ','
            << csv_cell(dataset.data_provenance) << ','
            << csv_cell(dataset.citation) << ','
            << csv_cell(dataset.license) << ','
            << csv_cell(config.gas_data_file.string()) << ','
            << "fixed_population_no_avalanche,"
            << config.neutral_density << ','
            << config.timestep << ','
            << config.steps << ','
            << config.burn_in_steps << ','
            << config.particles << ','
            << config.work_item_limit << ','
            << config.seed + static_cast<std::uint64_t>(row) << ','
            << result.reduced_field_td << ','
            << result.collision_model_signature << ','
            << result.electric_field_v_m << ','
            << result.mean_velocity_x_m_s << ','
            << result.mean_velocity_x_standard_error_m_s << ','
            << result.electron_drift_velocity_m_s << ','
            << result.reduced_mobility_1_v_m_s << ','
            << result.mean_energy_ev << ','
            << result.mean_energy_standard_error_ev << ','
            << result.longitudinal_diffusion_m2_s << ','
            << result.transverse_diffusion_m2_s << ','
            << result.maximum_observed_energy_ev << ','
            << result.collision_candidates << ','
            << result.null_collisions;
        for (const auto& channel : result.channels) {
            output << ',' << channel.collisions
                   << ',' << channel.rate_per_electron_s
                   << ',' << channel.poisson_standard_error_s;
        }
        output << '\n';
    }
    if (!output) {
        throw std::runtime_error(
            "failed while writing swarm output: " + path.string());
    }
}

} // namespace pic
