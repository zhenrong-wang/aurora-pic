#include "pic/GasDataset.hpp"
#include "pic/Collision.hpp"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace pic {
namespace {

using Values = std::map<std::string, std::string>;

struct ChannelBlock {
    std::string name;
    Values values;
};

std::string trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

std::string lower(std::string value) {
    std::transform(
        value.begin(), value.end(), value.begin(),
        [](unsigned char character) {
            return static_cast<char>(std::tolower(character));
        });
    return value;
}

[[noreturn]] void dataset_error(
    const std::filesystem::path& path,
    std::size_t line,
    const std::string& message) {
    throw std::runtime_error(
        "gas dataset '" + path.string() + "' line " +
        std::to_string(line) + ": " + message);
}

void assign(
    Values& values,
    std::string key,
    std::string value,
    const std::filesystem::path& path,
    std::size_t line) {
    key = lower(trim(std::move(key)));
    value = trim(std::move(value));
    if (key.empty() || value.empty()) {
        dataset_error(path, line, "empty key or value");
    }
    if (!values.emplace(key, value).second) {
        dataset_error(path, line, "duplicate key '" + key + "'");
    }
}

std::string required(
    const Values& values,
    const std::string& key,
    const std::string& context) {
    const auto found = values.find(key);
    if (found == values.end() || found->second.empty()) {
        throw std::runtime_error(
            context + " requires key '" + key + "'");
    }
    return found->second;
}

template <typename Number>
Number number(
    const Values& values,
    const std::string& key,
    Number fallback,
    const std::string& context) {
    const auto found = values.find(key);
    if (found == values.end()) return fallback;
    std::size_t used = 0;
    try {
        if constexpr (std::is_same_v<Number, double>) {
            const double value = std::stod(found->second, &used);
            if (used != found->second.size() || !std::isfinite(value)) {
                throw std::invalid_argument("invalid");
            }
            return value;
        } else {
            const auto value = std::stoull(found->second, &used);
            if (used != found->second.size()) {
                throw std::invalid_argument("invalid");
            }
            return static_cast<Number>(value);
        }
    } catch (const std::exception&) {
        throw std::runtime_error(
            context + " has invalid numeric key '" + key + "'");
    }
}

std::filesystem::path resolved_path(
    const std::filesystem::path& manifest,
    const std::string& value) {
    std::filesystem::path result(value);
    if (result.is_relative()) result = manifest.parent_path() / result;
    return std::filesystem::absolute(result).lexically_normal();
}

bool valid_iso_date(const std::string& value) {
    if (value.size() != 10 || value[4] != '-' || value[7] != '-') {
        return false;
    }
    for (std::size_t index = 0; index < value.size(); ++index) {
        if (index == 4 || index == 7) continue;
        if (!std::isdigit(static_cast<unsigned char>(value[index]))) {
            return false;
        }
    }
    const int year = std::stoi(value.substr(0, 4));
    const unsigned month =
        static_cast<unsigned>(std::stoul(value.substr(5, 2)));
    const unsigned day =
        static_cast<unsigned>(std::stoul(value.substr(8, 2)));
    return std::chrono::year_month_day{
               std::chrono::year{year},
               std::chrono::month{month},
               std::chrono::day{day}}
        .ok();
}

bool valid_identifier(const std::string& value) {
    return !value.empty() &&
           std::all_of(
               value.begin(), value.end(),
               [](unsigned char character) {
                   return std::isalnum(character) ||
                          character == '_' || character == '-' ||
                          character == '.' || character == ':';
               });
}

} // namespace

GasDataset load_gas_dataset(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error(
            "cannot open gas dataset: " + path.string());
    }

    static const std::set<std::string> global_keys{
        "gas_data_version", "units", "gas", "neutral_mass", "dataset_id",
        "dataset_version", "data_provenance", "citation", "retrieved",
        "license",
    };
    static const std::set<std::string> channel_keys{
        "type", "cross_section_file", "threshold_energy",
        "energy_scale", "cross_section_scale",
        "angular_model", "mean_cosine_file",
        "mean_cosine_energy_scale",
    };

    Values global;
    std::vector<ChannelBlock> channels;
    Values* current = &global;
    const std::set<std::string>* allowed = &global_keys;
    std::set<std::string> channel_names;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        const auto comment = line.find_first_of("#;");
        if (comment != std::string::npos) line.resize(comment);
        line = trim(line);
        if (line.empty()) continue;
        if (line.front() == '[' && line.back() == ']') {
            const std::string section =
                trim(line.substr(1, line.size() - 2));
            const std::string lowered = lower(section);
            constexpr const char* prefix = "collision.";
            if (lowered.rfind(prefix, 0) != 0) {
                dataset_error(
                    path, line_number,
                    "unknown section '" + section + "'");
            }
            const std::string name =
                trim(section.substr(
                    std::char_traits<char>::length(prefix)));
            if (name.empty() || !channel_names.insert(name).second) {
                dataset_error(
                    path, line_number,
                    "empty or duplicate collision channel");
            }
            channels.push_back({name, {}});
            current = &channels.back().values;
            allowed = &channel_keys;
            continue;
        }
        const auto equals = line.find('=');
        if (equals == std::string::npos) {
            dataset_error(path, line_number, "expected key = value");
        }
        const std::string key = lower(trim(line.substr(0, equals)));
        if (!allowed->contains(key)) {
            dataset_error(
                path, line_number, "unknown key '" + key + "'");
        }
        assign(
            *current, line.substr(0, equals), line.substr(equals + 1),
            path, line_number);
    }

    const std::string context = "gas dataset '" + path.string() + "'";
    const auto format_version =
        number<std::size_t>(
            global, "gas_data_version", 0, context);
    if (format_version != 1 && format_version != 2) {
        throw std::runtime_error(
            context + " supports gas_data_version 1 or 2");
    }

    GasDataset result;
    result.format_version = format_version;
    if (format_version == 2) {
        const std::string units =
            lower(required(global, "units", context));
        if (units == "si") {
            result.unit_system = UnitSystem::SI;
        } else if (units == "normalized") {
            result.unit_system = UnitSystem::Normalized;
        } else {
            throw std::runtime_error(
                context + " units must be si or normalized");
        }
    } else if (global.contains("units")) {
        throw std::runtime_error(
            context + " units requires gas_data_version = 2");
    }
    result.gas_name = required(global, "gas", context);
    (void)required(global, "neutral_mass", context);
    result.neutral_mass =
        number<double>(global, "neutral_mass", 0.0, context);
    if (!(result.neutral_mass > 0.0)) {
        throw std::runtime_error(
            context + " requires positive neutral_mass");
    }
    result.dataset_id = required(global, "dataset_id", context);
    if (!valid_identifier(result.dataset_id)) {
        throw std::runtime_error(
            context + " dataset_id contains unsupported characters");
    }
    result.dataset_version =
        required(global, "dataset_version", context);
    result.data_provenance =
        required(global, "data_provenance", context);
    result.citation = required(global, "citation", context);
    result.retrieved = required(global, "retrieved", context);
    result.license = required(global, "license", context);
    if (!valid_iso_date(result.retrieved)) {
        throw std::runtime_error(
            context + " retrieved must use a valid YYYY-MM-DD date");
    }
    if (channels.empty()) {
        throw std::runtime_error(
            context + " requires collision channel sections");
    }

    result.channels.reserve(channels.size());
    for (const auto& channel : channels) {
        const std::string channel_context =
            context + " channel '" + channel.name + "'";
        if (!valid_identifier(channel.name)) {
            throw std::runtime_error(
                channel_context + " name contains unsupported characters");
        }
        CollisionChannelConfig value;
        value.name = channel.name;
        const std::string type = lower(
            required(channel.values, "type", channel_context));
        if (type == "elastic") {
            value.process = CollisionProcessKind::Elastic;
        } else if (type == "excitation") {
            value.process = CollisionProcessKind::Excitation;
        } else if (type == "ionization") {
            value.process = CollisionProcessKind::Ionization;
        } else if (type == "attachment") {
            value.process = CollisionProcessKind::Attachment;
        } else if (type == "charge_exchange" ||
                   type == "charge-exchange") {
            value.process =
                CollisionProcessKind::ChargeExchange;
        } else {
            throw std::runtime_error(
                channel_context +
                " type must be elastic, excitation, ionization, "
                "attachment, or charge_exchange");
        }
        value.cross_section_file = resolved_path(
            path, required(
                      channel.values, "cross_section_file",
                      channel_context));
        value.threshold_energy = number<double>(
            channel.values, "threshold_energy",
            value.threshold_energy, channel_context);
        value.energy_scale = number<double>(
            channel.values, "energy_scale",
            value.energy_scale, channel_context);
        value.cross_section_scale = number<double>(
            channel.values, "cross_section_scale",
            value.cross_section_scale, channel_context);
        const bool has_angular_data =
            channel.values.contains("angular_model") ||
            channel.values.contains("mean_cosine_file") ||
            channel.values.contains("mean_cosine_energy_scale");
        if (has_angular_data && format_version != 2) {
            throw std::runtime_error(
                channel_context +
                " angular scattering requires gas_data_version = 2");
        }
        const std::string angular_model = lower(
            channel.values.contains("angular_model")
                ? channel.values.at("angular_model")
                : "isotropic");
        if (angular_model == "isotropic") {
            value.angular_scattering =
                AngularScatteringKind::Isotropic;
        } else if (
            angular_model == "henyey_greenstein" ||
            angular_model == "henyey-greenstein") {
            value.angular_scattering =
                AngularScatteringKind::HenyeyGreenstein;
        } else {
            throw std::runtime_error(
                channel_context +
                " angular_model must be isotropic or "
                "henyey_greenstein");
        }
        if (channel.values.contains("mean_cosine_file")) {
            value.mean_cosine_file = resolved_path(
                path, channel.values.at("mean_cosine_file"));
        }
        value.mean_cosine_energy_scale = number<double>(
            channel.values, "mean_cosine_energy_scale",
            value.mean_cosine_energy_scale, channel_context);
        if ((value.process == CollisionProcessKind::Elastic ||
             value.process == CollisionProcessKind::Attachment ||
             value.process ==
                 CollisionProcessKind::ChargeExchange) &&
            value.threshold_energy != 0.0) {
            throw std::runtime_error(
                channel_context +
                " elastic, attachment, and charge-exchange "
                "threshold_energy must be zero");
        }
        if ((value.process == CollisionProcessKind::Excitation ||
             value.process == CollisionProcessKind::Ionization) &&
            !(value.threshold_energy > 0.0)) {
            throw std::runtime_error(
                channel_context +
                " inelastic threshold_energy must be positive");
        }
        if (value.angular_scattering !=
                AngularScatteringKind::Isotropic &&
            value.process != CollisionProcessKind::Elastic) {
            throw std::runtime_error(
                channel_context +
                " anisotropic scattering is valid only for elastic channels");
        }
        if (value.angular_scattering ==
            AngularScatteringKind::HenyeyGreenstein) {
            if (value.mean_cosine_file.empty()) {
                throw std::runtime_error(
                    channel_context +
                    " Henyey-Greenstein scattering requires "
                    "mean_cosine_file");
            }
        } else if (!value.mean_cosine_file.empty() ||
                   value.mean_cosine_energy_scale != 1.0) {
            throw std::runtime_error(
                channel_context +
                " mean-cosine data requires angular_model = "
                "henyey_greenstein");
        }
        const CrossSectionTable cross_section(
            value.cross_section_file,
            value.energy_scale,
            value.cross_section_scale);
        if (value.angular_scattering ==
            AngularScatteringKind::HenyeyGreenstein) {
            const MeanCosineTable mean_cosine(
                value.mean_cosine_file,
                value.mean_cosine_energy_scale);
            if (mean_cosine.energies().front() >
                    cross_section.energies().front() ||
                mean_cosine.energies().back() <
                    cross_section.energies().back()) {
                throw std::runtime_error(
                    channel_context +
                    " mean-cosine energy range must cover the "
                    "cross-section table");
            }
        }
        result.channels.push_back(std::move(value));
    }
    return result;
}

} // namespace pic
