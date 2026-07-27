#pragma once
#include "pic/Mesh2D.hpp"
#include "pic/Mesh3D.hpp"
#include "pic/UnstructuredMesh2D.hpp"
#include <filesystem>
#include <string>

namespace pic {
void write_legacy_vtk(const Mesh2D& mesh, const std::filesystem::path& path, const std::string& title = "AuroraPIC 2D fields");
void write_legacy_vtk(const Mesh3D& mesh, const std::filesystem::path& path, const std::string& title = "AuroraPIC 3D fields");
void write_vtk_xml(const Mesh2D& mesh, const std::filesystem::path& path);
void write_vtk_xml(const Mesh3D& mesh, const std::filesystem::path& path);
void write_vtk_xml(const UnstructuredMesh2D& mesh, const std::filesystem::path& path);
}
