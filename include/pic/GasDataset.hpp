#pragma once

#include "pic/Config.hpp"

#include <filesystem>
#include <string>
#include <vector>

namespace pic {

struct GasDataset {
    std::string gas_name{};
    double neutral_mass{0.0};
    std::string dataset_id{};
    std::string dataset_version{};
    std::string data_provenance{};
    std::string citation{};
    std::string retrieved{};
    std::string license{};
    std::vector<CollisionChannelConfig> channels{};
};

GasDataset load_gas_dataset(const std::filesystem::path& path);

} // namespace pic
