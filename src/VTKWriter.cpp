#include "pic/VTKWriter.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

namespace pic {
namespace {
void ensure_parent_directory(const std::filesystem::path& path) {
    if (path.has_parent_path()) std::filesystem::create_directories(path.parent_path());
}

void require_open(const std::ofstream& out, const std::filesystem::path& path, const char* format) {
    if (!out) throw std::runtime_error(std::string("cannot open ") + format + " VTK output: " + path.string());
}

void write_scalar_data_array(std::ofstream& out, const std::string& name, const std::vector<double>& values) {
    out << "        <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
    out << "          ";
    for (double value : values) out << value << ' ';
    out << "\n";
    out << "        </DataArray>\n";
}

void write_vtk_xml_header(std::ofstream& out,
                          std::size_t nx,
                          std::size_t ny,
                          std::size_t nz) {
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <StructuredGrid WholeExtent=\"0 " << (nx - 1)
        << " 0 " << (ny - 1)
        << " 0 " << (nz - 1) << "\">\n";
    out << "    <Piece Extent=\"0 " << (nx - 1)
        << " 0 " << (ny - 1)
        << " 0 " << (nz - 1) << "\">\n";
}

void write_vtk_xml_footer(std::ofstream& out) {
    out << "    </Piece>\n";
    out << "  </StructuredGrid>\n";
    out << "</VTKFile>\n";
}
} // namespace

void write_legacy_vtk(const Mesh2D& mesh, const std::filesystem::path& path, const std::string& title) {
    ensure_parent_directory(path);

    std::ofstream out(path);
    require_open(out, path, "legacy");

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
    ensure_parent_directory(path);

    std::ofstream out(path);
    require_open(out, path, "legacy");

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

void write_vtk_xml(const Mesh2D& mesh, const std::filesystem::path& path) {
    ensure_parent_directory(path);

    std::ofstream out(path);
    require_open(out, path, "XML");
    out << std::setprecision(17);

    write_vtk_xml_header(out, mesh.nx(), mesh.ny(), 1);
    out << "      <PointData Scalars=\"rho\" Vectors=\"electric\">\n";
    write_scalar_data_array(out, "rho", mesh.rho());
    write_scalar_data_array(out, "phi", mesh.phi());
    out << "        <DataArray type=\"Float64\" Name=\"electric\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    out << "          ";
    for (std::size_t idx = 0; idx < mesh.size(); ++idx) {
        out << mesh.electric_x()[idx] << ' ' << mesh.electric_y()[idx] << " 0 ";
    }
    out << "\n";
    out << "        </DataArray>\n";
    out << "      </PointData>\n";
    out << "      <CellData/>\n";
    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    out << "          ";
    for (std::size_t j = 0; j < mesh.ny(); ++j) {
        for (std::size_t i = 0; i < mesh.nx(); ++i) {
            out << mesh.node_x(i) << ' ' << mesh.node_y(j) << " 0 ";
        }
    }
    out << "\n";
    out << "        </DataArray>\n";
    out << "      </Points>\n";
    write_vtk_xml_footer(out);
}

void write_vtk_xml(const Mesh3D& mesh, const std::filesystem::path& path) {
    ensure_parent_directory(path);

    std::ofstream out(path);
    require_open(out, path, "XML");
    out << std::setprecision(17);

    write_vtk_xml_header(out, mesh.nx(), mesh.ny(), mesh.nz());
    out << "      <PointData Scalars=\"rho\" Vectors=\"electric\">\n";
    write_scalar_data_array(out, "rho", mesh.rho());
    write_scalar_data_array(out, "phi", mesh.phi());
    out << "        <DataArray type=\"Float64\" Name=\"electric\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    out << "          ";
    for (std::size_t idx = 0; idx < mesh.size(); ++idx) {
        out << mesh.electric_x()[idx] << ' ' << mesh.electric_y()[idx] << ' ' << mesh.electric_z()[idx] << ' ';
    }
    out << "\n";
    out << "        </DataArray>\n";
    out << "      </PointData>\n";
    out << "      <CellData/>\n";
    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    out << "          ";
    for (std::size_t k = 0; k < mesh.nz(); ++k) {
        for (std::size_t j = 0; j < mesh.ny(); ++j) {
            for (std::size_t i = 0; i < mesh.nx(); ++i) {
                out << mesh.node_x(i) << ' ' << mesh.node_y(j) << ' ' << mesh.node_z(k) << ' ';
            }
        }
    }
    out << "\n";
    out << "        </DataArray>\n";
    out << "      </Points>\n";
    write_vtk_xml_footer(out);
}
} // namespace pic
