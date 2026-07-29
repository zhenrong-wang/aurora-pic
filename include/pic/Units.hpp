#pragma once

#include "pic/Types.hpp"

#include <cstddef>
#include <filesystem>
#include <optional>

namespace pic {

void write_unit_metadata(
    const std::filesystem::path& output_dir,
    const UnitSystemConfig& units,
    std::size_t spatial_dimension,
    std::optional<double> out_of_plane_depth = {});

} // namespace pic
