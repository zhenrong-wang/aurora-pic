#include "pic/ImportedMesh2D.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace pic {
namespace {
constexpr double Z_TOLERANCE = 1e-12;

std::string trim(std::string text) {
    const auto begin = text.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos) return "";
    const auto end = text.find_last_not_of(" \t\r\n");
    return text.substr(begin, end - begin + 1);
}

std::string physical_fallback(int dimension, int tag) {
    if (tag == 0) return "unlabeled";
    return std::string{dimension == 1 ? "boundary" : "region"} + "_physical_" + std::to_string(tag);
}

std::string read_required_line(std::istream& input, const std::string& context) {
    std::string line;
    if (!std::getline(input, line)) {
        throw std::runtime_error("unexpected end of Gmsh file while reading " + context);
    }
    return trim(line);
}

void require_section_end(std::istream& input, const std::string& section_name) {
    const std::string expected = "$End" + section_name.substr(1);
    const std::string actual = read_required_line(input, expected);
    if (actual != expected) {
        throw std::runtime_error("expected " + expected + ", found " + actual);
    }
}

std::string parse_physical_name_line(const std::string& line, int& dimension, int& tag) {
    std::istringstream stream(line);
    if (!(stream >> dimension >> tag)) {
        throw std::runtime_error("invalid $PhysicalNames row: " + line);
    }
    const auto quote_begin = line.find('"');
    const auto quote_end = line.find_last_of('"');
    if (quote_begin == std::string::npos || quote_end == quote_begin) {
        throw std::runtime_error("physical name must be quoted: " + line);
    }
    return line.substr(quote_begin + 1, quote_end - quote_begin - 1);
}

bool is_supported_cell_type(int type) {
    return type == 2 || type == 3; // Gmsh v2 triangle and quadrilateral.
}

std::size_t required_node_count(int type) {
    switch (type) {
        case 1: return 2; // line
        case 2: return 3; // triangle
        case 3: return 4; // quadrilateral
        default: return 0;
    }
}

ImportedCellShape2D cell_shape_for(int type) {
    if (type == 2) return ImportedCellShape2D::Triangle;
    if (type == 3) return ImportedCellShape2D::Quadrilateral;
    throw std::runtime_error("unsupported 2D cell type");
}
} // namespace

void ImportedMesh2D::add_node(ImportedMeshNode2D node) {
    if (node.id == 0) throw std::invalid_argument("imported mesh node id must be positive");
    const auto duplicate = std::find_if(nodes_.begin(), nodes_.end(), [&](const auto& existing) {
        return existing.id == node.id;
    });
    if (duplicate != nodes_.end()) {
        throw std::invalid_argument("duplicate imported mesh node id: " + std::to_string(node.id));
    }
    nodes_.push_back(node);
}

void ImportedMesh2D::add_cell(ImportedCell2D cell) {
    if (cell.id == 0) throw std::invalid_argument("imported mesh cell id must be positive");
    const std::size_t expected = cell.shape == ImportedCellShape2D::Triangle ? 3 : 4;
    if (cell.node_ids.size() != expected) {
        throw std::invalid_argument("imported mesh cell has wrong node count");
    }
    if (cell.label.empty()) cell.label = label_for_physical_tag(2, cell.physical_tag);
    cells_.push_back(std::move(cell));
}

void ImportedMesh2D::add_boundary_face(ImportedBoundaryFace2D face) {
    if (face.id == 0) throw std::invalid_argument("imported boundary face id must be positive");
    if (face.node_ids[0] == face.node_ids[1]) {
        throw std::invalid_argument("imported boundary face has duplicate node ids");
    }
    if (face.label.empty()) face.label = label_for_physical_tag(1, face.physical_tag);
    boundary_faces_.push_back(std::move(face));
}

void ImportedMesh2D::set_physical_name(int dimension, int tag, std::string name) {
    if (dimension != 1 && dimension != 2) {
        throw std::invalid_argument("physical name dimension must be 1 or 2 for imported 2D meshes");
    }
    if (tag <= 0) throw std::invalid_argument("physical tag must be positive");
    if (name.empty()) throw std::invalid_argument("physical name must not be empty");
    physical_names_[{dimension, tag}] = std::move(name);
}

const ImportedMeshNode2D& ImportedMesh2D::node_by_id(std::size_t id) const {
    const auto it = std::find_if(nodes_.begin(), nodes_.end(), [&](const auto& node) { return node.id == id; });
    if (it == nodes_.end()) throw std::out_of_range("imported mesh node id not found: " + std::to_string(id));
    return *it;
}

std::string ImportedMesh2D::label_for_physical_tag(int dimension, int tag) const {
    const auto it = physical_names_.find({dimension, tag});
    return it == physical_names_.end() ? physical_fallback(dimension, tag) : it->second;
}

std::vector<std::string> ImportedMesh2D::boundary_labels() const {
    std::set<std::string> labels;
    for (const auto& face : boundary_faces_) labels.insert(face.label);
    return {labels.begin(), labels.end()};
}

Vec2 ImportedMesh2D::min_corner() const {
    if (nodes_.empty()) throw std::runtime_error("imported mesh has no nodes");
    Vec2 result{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    for (const auto& node : nodes_) {
        result.x = std::min(result.x, node.position.x);
        result.y = std::min(result.y, node.position.y);
    }
    return result;
}

Vec2 ImportedMesh2D::max_corner() const {
    if (nodes_.empty()) throw std::runtime_error("imported mesh has no nodes");
    Vec2 result{-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    for (const auto& node : nodes_) {
        result.x = std::max(result.x, node.position.x);
        result.y = std::max(result.y, node.position.y);
    }
    return result;
}

void ImportedMesh2D::validate() const {
    if (nodes_.empty()) throw std::runtime_error("imported mesh must contain nodes");
    if (cells_.empty()) throw std::runtime_error("imported 2D mesh must contain triangle or quadrilateral cells");
    if (boundary_faces_.empty()) throw std::runtime_error("imported 2D mesh must contain tagged line boundary faces");

    std::set<std::size_t> node_ids;
    for (const auto& node : nodes_) {
        node_ids.insert(node.id);
    }
    const auto require_node = [&](std::size_t id, const std::string& owner) {
        if (!node_ids.contains(id)) throw std::runtime_error(owner + " references missing node id " + std::to_string(id));
    };
    for (const auto& cell : cells_) {
        for (const auto id : cell.node_ids) require_node(id, "imported cell " + std::to_string(cell.id));
        if (cell.label.empty()) throw std::runtime_error("imported cell has empty label");
    }
    for (const auto& face : boundary_faces_) {
        require_node(face.node_ids[0], "imported boundary face " + std::to_string(face.id));
        require_node(face.node_ids[1], "imported boundary face " + std::to_string(face.id));
        if (face.label.empty()) throw std::runtime_error("imported boundary face has empty label");
    }
}

ImportedMesh2D load_gmsh2_ascii_mesh2d(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open Gmsh mesh: " + path.string());

    ImportedMesh2D mesh;
    bool saw_mesh_format = false;
    bool saw_nodes = false;
    bool saw_elements = false;

    std::string section;
    while (std::getline(input, section)) {
        section = trim(section);
        if (section.empty()) continue;

        if (section == "$MeshFormat") {
            const std::string row = read_required_line(input, section);
            std::istringstream stream(row);
            double version = 0.0;
            int file_type = -1;
            int data_size = 0;
            if (!(stream >> version >> file_type >> data_size)) {
                throw std::runtime_error("invalid $MeshFormat row: " + row);
            }
            if (version < 2.0 || version >= 3.0) {
                throw std::runtime_error("only Gmsh v2 ASCII meshes are supported");
            }
            if (file_type != 0) {
                throw std::runtime_error("binary Gmsh meshes are not supported");
            }
            saw_mesh_format = true;
            require_section_end(input, section);
        } else if (section == "$PhysicalNames") {
            const std::size_t count = static_cast<std::size_t>(std::stoull(read_required_line(input, section)));
            for (std::size_t i = 0; i < count; ++i) {
                int dimension = 0;
                int tag = 0;
                const std::string name = parse_physical_name_line(read_required_line(input, section), dimension, tag);
                if (dimension == 1 || dimension == 2) mesh.set_physical_name(dimension, tag, name);
            }
            require_section_end(input, section);
        } else if (section == "$Nodes") {
            const std::size_t count = static_cast<std::size_t>(std::stoull(read_required_line(input, section)));
            for (std::size_t i = 0; i < count; ++i) {
                const std::string row = read_required_line(input, section);
                std::istringstream stream(row);
                std::size_t id = 0;
                double x = 0.0;
                double y = 0.0;
                double z = 0.0;
                if (!(stream >> id >> x >> y >> z)) throw std::runtime_error("invalid $Nodes row: " + row);
                if (std::abs(z) > Z_TOLERANCE) {
                    throw std::runtime_error("2D Gmsh importer requires z=0 nodes; node " + std::to_string(id) + " has z=" + std::to_string(z));
                }
                mesh.add_node(ImportedMeshNode2D{id, Vec2{x, y}});
            }
            saw_nodes = true;
            require_section_end(input, section);
        } else if (section == "$Elements") {
            const std::size_t count = static_cast<std::size_t>(std::stoull(read_required_line(input, section)));
            for (std::size_t i = 0; i < count; ++i) {
                const std::string row = read_required_line(input, section);
                std::istringstream stream(row);
                std::size_t id = 0;
                int type = 0;
                int tag_count = 0;
                if (!(stream >> id >> type >> tag_count)) throw std::runtime_error("invalid $Elements row: " + row);
                if (tag_count < 0) throw std::runtime_error("invalid negative Gmsh element tag count: " + row);
                std::vector<int> tags(static_cast<std::size_t>(tag_count));
                for (auto& tag : tags) {
                    if (!(stream >> tag)) throw std::runtime_error("invalid $Elements tag list: " + row);
                }
                const int physical_tag = tags.empty() ? 0 : tags.front();
                const std::size_t node_count = required_node_count(type);
                std::vector<std::size_t> node_ids(node_count);
                for (auto& node_id : node_ids) {
                    if (!(stream >> node_id)) throw std::runtime_error("invalid $Elements node list: " + row);
                }
                if (type == 1) {
                    mesh.add_boundary_face(ImportedBoundaryFace2D{id, {node_ids[0], node_ids[1]},
                                                                 physical_tag, mesh.label_for_physical_tag(1, physical_tag)});
                } else if (is_supported_cell_type(type)) {
                    mesh.add_cell(ImportedCell2D{id, cell_shape_for(type), std::move(node_ids),
                                                 physical_tag, mesh.label_for_physical_tag(2, physical_tag)});
                }
            }
            saw_elements = true;
            require_section_end(input, section);
        } else if (section.starts_with("$End")) {
            throw std::runtime_error("unexpected Gmsh section terminator: " + section);
        } else if (section.starts_with("$")) {
            const std::string end_section = "$End" + section.substr(1);
            std::string row;
            bool closed = false;
            while (std::getline(input, row)) {
                if (trim(row) == end_section) {
                    closed = true;
                    break;
                }
            }
            if (!closed) throw std::runtime_error("unterminated unsupported Gmsh section: " + section);
        }
    }

    if (!saw_mesh_format) throw std::runtime_error("Gmsh file is missing $MeshFormat section");
    if (!saw_nodes) throw std::runtime_error("Gmsh file is missing $Nodes section");
    if (!saw_elements) throw std::runtime_error("Gmsh file is missing $Elements section");
    mesh.validate();
    return mesh;
}

} // namespace pic
