#pragma once
#include "pic/Mesh2D.hpp"
#include <filesystem>
#include <string>

namespace pic {
void write_legacy_vtk(const Mesh2D& mesh, const std::filesystem::path& path, const std::string& title = "AuroraPIC 2D fields");
}
