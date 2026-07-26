#include "pic/VTKWriter.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace pic {
void write_legacy_vtk(const Mesh2D& mesh, const std::filesystem::path& path, const std::string& title) {
    if (path.has_parent_path()) std::filesystem::create_directories(path.parent_path());

    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot open VTK output: " + path.string());

    out << "# vtk DataFile Version 3.0\n";
    out << title << "\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_GRID\n";
    out << "DIMENSIONS " << mesh.nx() << ' ' << mesh.ny() << " 1\n";
    out << "POINTS " << mesh.size() << " double\n";
    out << std::setprecision(17);
    for (std::size_t j = 0; j < mesh.ny(); ++j) {
        for (std::size_t i = 0; i < mesh.nx(); ++i) {
            out << mesh.node_x(i) << ' ' << mesh.node_y(j) << " 0\n";
        }
    }

    out << "POINT_DATA " << mesh.size() << "\n";
    out << "SCALARS rho double 1\nLOOKUP_TABLE default\n";
    for (double value : mesh.rho()) out << value << '\n';

    out << "SCALARS phi double 1\nLOOKUP_TABLE default\n";
    for (double value : mesh.phi()) out << value << '\n';

    out << "VECTORS electric double\n";
    for (std::size_t idx = 0; idx < mesh.size(); ++idx) {
        out << mesh.electric_x()[idx] << ' ' << mesh.electric_y()[idx] << " 0\n";
    }
}

void write_legacy_vtk(const Mesh3D& mesh, const std::filesystem::path& path, const std::string& title) {
    if (path.has_parent_path()) std::filesystem::create_directories(path.parent_path());

    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot open VTK output: " + path.string());

    out << "# vtk DataFile Version 3.0\n";
    out << title << "\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_GRID\n";
    out << "DIMENSIONS " << mesh.nx() << ' ' << mesh.ny() << ' ' << mesh.nz() << "\n";
    out << "POINTS " << mesh.size() << " double\n";
    out << std::setprecision(17);
    for (std::size_t k = 0; k < mesh.nz(); ++k) {
        for (std::size_t j = 0; j < mesh.ny(); ++j) {
            for (std::size_t i = 0; i < mesh.nx(); ++i) {
                out << mesh.node_x(i) << ' ' << mesh.node_y(j) << ' ' << mesh.node_z(k) << "\n";
            }
        }
    }

    out << "POINT_DATA " << mesh.size() << "\n";
    out << "SCALARS rho double 1\nLOOKUP_TABLE default\n";
    for (double value : mesh.rho()) out << value << '\n';

    out << "SCALARS phi double 1\nLOOKUP_TABLE default\n";
    for (double value : mesh.phi()) out << value << '\n';

    out << "VECTORS electric double\n";
    for (std::size_t idx = 0; idx < mesh.size(); ++idx) {
        out << mesh.electric_x()[idx] << ' ' << mesh.electric_y()[idx] << ' ' << mesh.electric_z()[idx] << "\n";
    }
}
}
