#pragma once

#include "pic/Types.hpp"

#include <cstddef>
#include <filesystem>

namespace pic {

void write_unit_metadata(
    const std::filesystem::path& output_dir,
    const UnitSystemConfig& units,
    std::size_t spatial_dimension);

} // namespace pic
