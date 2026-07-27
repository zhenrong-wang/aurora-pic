#include "pic/ImportedMesh2D.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace pic {
namespace {
constexpr double Z_TOLERANCE = 1e-12;
using Edge = std::pair<std::size_t, std::size_t>;

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

void require_no_trailing_tokens(std::istringstream& stream, const std::string& context) {
    std::string extra;
    if (stream >> extra) throw std::runtime_error(context + " has trailing data");
}

std::size_t parse_bounded_count(const std::string& row, std::size_t limit, const std::string& context) {
    if (row.empty() || !std::all_of(row.begin(), row.end(),
                                    [](unsigned char character) { return std::isdigit(character) != 0; })) {
        throw std::runtime_error("invalid " + context + " count: " + row);
    }
    std::size_t consumed = 0;
    unsigned long long value = 0;
    try {
        value = std::stoull(row, &consumed);
    } catch (const std::exception&) {
        throw std::runtime_error("invalid " + context + " count: " + row);
    }
    if (consumed != row.size() || value > std::numeric_limits<std::size_t>::max()) {
        throw std::runtime_error("invalid " + context + " count: " + row);
    }
    const auto result = static_cast<std::size_t>(value);
    if (result > limit) {
        throw std::runtime_error(context + " count exceeds configured import limit");
    }
    return result;
}

std::size_t parse_positive_id(const std::string& token, const std::string& context) {
    const std::size_t id = parse_bounded_count(token, std::numeric_limits<std::size_t>::max(), context);
    if (id == 0) throw std::runtime_error(context + " must be positive");
    return id;
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
    std::string quoted_name;
    std::getline(stream, quoted_name);
    quoted_name = trim(quoted_name);
    if (quoted_name.size() < 2 || quoted_name.front() != '"' || quoted_name.back() != '"' ||
        quoted_name.find('"', 1) != quoted_name.size() - 1) {
        throw std::runtime_error("physical name must be quoted: " + line);
    }
    const std::string name = quoted_name.substr(1, quoted_name.size() - 2);
    if (name.empty()) throw std::runtime_error("physical name must not be empty");
    return name;
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

Edge canonical_edge(std::size_t first, std::size_t second) {
    return std::minmax(first, second);
}

double cross(Vec2 first, Vec2 second, Vec2 third) {
    return (second.x - first.x) * (third.y - first.y) -
           (second.y - first.y) * (third.x - first.x);
}

double signed_polygon_area(const ImportedMesh2D& mesh, const ImportedCell2D& cell) {
    double twice_area = 0.0;
    for (std::size_t i = 0; i < cell.node_ids.size(); ++i) {
        const auto& current = mesh.node_by_id(cell.node_ids[i]).position;
        const auto& next = mesh.node_by_id(cell.node_ids[(i + 1) % cell.node_ids.size()]).position;
        twice_area += current.x * next.y - next.x * current.y;
    }
    return 0.5 * twice_area;
}
} // namespace

void ImportedMesh2D::add_node(ImportedMeshNode2D node) {
    if (node.id == 0) throw std::invalid_argument("imported mesh node id must be positive");
    if (!std::isfinite(node.position.x) || !std::isfinite(node.position.y)) {
        throw std::invalid_argument("imported mesh node coordinates must be finite");
    }
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
    if (std::any_of(cells_.begin(), cells_.end(), [&](const auto& existing) { return existing.id == cell.id; }) ||
        std::any_of(boundary_faces_.begin(), boundary_faces_.end(),
                    [&](const auto& existing) { return existing.id == cell.id; })) {
        throw std::invalid_argument("duplicate imported element id: " + std::to_string(cell.id));
    }
    const std::size_t expected = cell.shape == ImportedCellShape2D::Triangle ? 3 : 4;
    if (cell.node_ids.size() != expected) {
        throw std::invalid_argument("imported mesh cell has wrong node count");
    }
    const std::set<std::size_t> unique_nodes(cell.node_ids.begin(), cell.node_ids.end());
    if (unique_nodes.size() != cell.node_ids.size()) {
        throw std::invalid_argument("imported mesh cell has duplicate node ids");
    }
    if (cell.label.empty()) cell.label = label_for_physical_tag(2, cell.physical_tag);
    cells_.push_back(std::move(cell));
}

void ImportedMesh2D::add_boundary_face(ImportedBoundaryFace2D face) {
    if (face.id == 0) throw std::invalid_argument("imported boundary face id must be positive");
    if (std::any_of(cells_.begin(), cells_.end(), [&](const auto& existing) { return existing.id == face.id; }) ||
        std::any_of(boundary_faces_.begin(), boundary_faces_.end(),
                    [&](const auto& existing) { return existing.id == face.id; })) {
        throw std::invalid_argument("duplicate imported element id: " + std::to_string(face.id));
    }
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
    const bool inserted = physical_names_.emplace(std::make_pair(dimension, tag), std::move(name)).second;
    if (!inserted) {
        throw std::invalid_argument("duplicate imported physical tag " + std::to_string(tag) +
                                    " for dimension " + std::to_string(dimension));
    }
}

const ImportedMeshNode2D& ImportedMesh2D::node_by_id(std::size_t id) const {
    const auto it = std::find_if(nodes_.begin(), nodes_.end(), [&](const auto& node) { return node.id == id; });
    if (it == nodes_.end()) throw std::out_of_range("imported mesh node id not found: " + std::to_string(id));
    return *it;
}

const ImportedCell2D& ImportedMesh2D::cell_by_id(std::size_t id) const {
    const auto it = std::find_if(cells_.begin(), cells_.end(), [&](const auto& cell) { return cell.id == id; });
    if (it == cells_.end()) throw std::out_of_range("imported mesh cell id not found: " + std::to_string(id));
    return *it;
}

const ImportedBoundaryFace2D& ImportedMesh2D::boundary_face_by_id(std::size_t id) const {
    const auto it = std::find_if(boundary_faces_.begin(), boundary_faces_.end(),
                                 [&](const auto& face) { return face.id == id; });
    if (it == boundary_faces_.end()) {
        throw std::out_of_range("imported boundary face id not found: " + std::to_string(id));
    }
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

double ImportedMesh2D::cell_area(std::size_t id) const {
    return std::abs(signed_polygon_area(*this, cell_by_id(id)));
}

Vec2 ImportedMesh2D::cell_centroid(std::size_t id) const {
    const auto& cell = cell_by_id(id);
    double twice_area = 0.0;
    Vec2 weighted{};
    for (std::size_t i = 0; i < cell.node_ids.size(); ++i) {
        const auto& current = node_by_id(cell.node_ids[i]).position;
        const auto& next = node_by_id(cell.node_ids[(i + 1) % cell.node_ids.size()]).position;
        const double factor = current.x * next.y - next.x * current.y;
        twice_area += factor;
        weighted.x += (current.x + next.x) * factor;
        weighted.y += (current.y + next.y) * factor;
    }
    if (twice_area == 0.0 || !std::isfinite(twice_area)) {
        throw std::runtime_error("cannot compute centroid of degenerate imported cell " + std::to_string(id));
    }
    return {weighted.x / (3.0 * twice_area), weighted.y / (3.0 * twice_area)};
}

double ImportedMesh2D::boundary_face_length(std::size_t id) const {
    const auto& face = boundary_face_by_id(id);
    const auto& first = node_by_id(face.node_ids[0]).position;
    const auto& second = node_by_id(face.node_ids[1]).position;
    return std::hypot(second.x - first.x, second.y - first.y);
}

double ImportedMesh2D::total_area() const {
    double result = 0.0;
    for (const auto& cell : cells_) result += cell_area(cell.id);
    return result;
}

void ImportedMesh2D::validate() const {
    if (nodes_.empty()) throw std::runtime_error("imported mesh must contain nodes");
    if (cells_.empty()) throw std::runtime_error("imported 2D mesh must contain triangle or quadrilateral cells");
    if (boundary_faces_.empty()) throw std::runtime_error("imported 2D mesh must contain tagged line boundary faces");

    const Vec2 minimum = min_corner();
    const Vec2 maximum = max_corner();
    const double characteristic_length = std::max(maximum.x - minimum.x, maximum.y - minimum.y);
    if (!(characteristic_length > 0.0) || !std::isfinite(characteristic_length)) {
        throw std::runtime_error("imported mesh bounding box must have finite nonzero extent");
    }
    const double length_tolerance =
        128.0 * std::numeric_limits<double>::epsilon() * characteristic_length;
    const double area_tolerance =
        128.0 * std::numeric_limits<double>::epsilon() * characteristic_length * characteristic_length;

    std::set<std::size_t> node_ids;
    for (const auto& node : nodes_) {
        if (node.id == 0 || !node_ids.insert(node.id).second) {
            throw std::runtime_error("imported mesh contains an invalid or duplicate node id");
        }
        if (!std::isfinite(node.position.x) || !std::isfinite(node.position.y)) {
            throw std::runtime_error("imported mesh contains a non-finite node coordinate");
        }
    }
    const auto require_node = [&](std::size_t id, const std::string& owner) {
        if (!node_ids.contains(id)) throw std::runtime_error(owner + " references missing node id " + std::to_string(id));
    };
    std::set<std::size_t> element_ids;
    std::set<std::size_t> used_node_ids;
    std::map<Edge, std::size_t> cell_edge_counts;
    int orientation = 0;
    for (const auto& cell : cells_) {
        if (cell.id == 0 || !element_ids.insert(cell.id).second) {
            throw std::runtime_error("imported mesh contains an invalid or duplicate element id");
        }
        const std::size_t expected = cell.shape == ImportedCellShape2D::Triangle ? 3 : 4;
        if (cell.node_ids.size() != expected ||
            std::set<std::size_t>(cell.node_ids.begin(), cell.node_ids.end()).size() != expected) {
            throw std::runtime_error("imported cell " + std::to_string(cell.id) +
                                     " has invalid node connectivity");
        }
        for (const auto id : cell.node_ids) require_node(id, "imported cell " + std::to_string(cell.id));
        used_node_ids.insert(cell.node_ids.begin(), cell.node_ids.end());
        if (cell.label.empty()) throw std::runtime_error("imported cell has empty label");

        const double signed_area = signed_polygon_area(*this, cell);
        if (!std::isfinite(signed_area) || std::abs(signed_area) <= area_tolerance) {
            throw std::runtime_error("imported cell " + std::to_string(cell.id) + " is degenerate");
        }
        const int cell_orientation = signed_area > 0.0 ? 1 : -1;
        if (orientation == 0) orientation = cell_orientation;
        if (orientation != cell_orientation) {
            throw std::runtime_error("imported cells must use a consistent node orientation");
        }

        if (cell.shape == ImportedCellShape2D::Quadrilateral) {
            int corner_orientation = 0;
            for (std::size_t i = 0; i < cell.node_ids.size(); ++i) {
                const auto& first = node_by_id(cell.node_ids[i]).position;
                const auto& second = node_by_id(cell.node_ids[(i + 1) % 4]).position;
                const auto& third = node_by_id(cell.node_ids[(i + 2) % 4]).position;
                const double turn = cross(first, second, third);
                if (!std::isfinite(turn) || std::abs(turn) <= area_tolerance) {
                    throw std::runtime_error("imported quadrilateral cell " + std::to_string(cell.id) +
                                             " has a degenerate corner");
                }
                const int turn_orientation = turn > 0.0 ? 1 : -1;
                if (corner_orientation == 0) corner_orientation = turn_orientation;
                if (corner_orientation != turn_orientation) {
                    throw std::runtime_error("imported quadrilateral cell " + std::to_string(cell.id) +
                                             " must be convex and non-self-intersecting");
                }
            }
        }

        for (std::size_t i = 0; i < cell.node_ids.size(); ++i) {
            const Edge edge = canonical_edge(cell.node_ids[i], cell.node_ids[(i + 1) % cell.node_ids.size()]);
            auto& count = cell_edge_counts[edge];
            ++count;
            if (count > 2) throw std::runtime_error("imported mesh contains a non-manifold cell edge");
        }
    }

    std::set<Edge> actual_boundary_edges;
    for (const auto& face : boundary_faces_) {
        if (face.id == 0 || !element_ids.insert(face.id).second) {
            throw std::runtime_error("imported mesh contains an invalid or duplicate element id");
        }
        require_node(face.node_ids[0], "imported boundary face " + std::to_string(face.id));
        require_node(face.node_ids[1], "imported boundary face " + std::to_string(face.id));
        if (face.node_ids[0] == face.node_ids[1]) {
            throw std::runtime_error("imported boundary face has duplicate node ids");
        }
        if (face.label.empty()) throw std::runtime_error("imported boundary face has empty label");
        if (boundary_face_length(face.id) <= length_tolerance) {
            throw std::runtime_error("imported boundary face " + std::to_string(face.id) + " is degenerate");
        }
        if (!actual_boundary_edges.insert(canonical_edge(face.node_ids[0], face.node_ids[1])).second) {
            throw std::runtime_error("imported mesh contains duplicate boundary faces");
        }
    }

    if (used_node_ids.size() != node_ids.size()) {
        throw std::runtime_error("imported mesh contains nodes not referenced by any cell");
    }

    std::set<Edge> expected_boundary_edges;
    for (const auto& [edge, count] : cell_edge_counts) {
        if (count == 1) expected_boundary_edges.insert(edge);
    }
    if (actual_boundary_edges != expected_boundary_edges) {
        throw std::runtime_error("tagged boundary faces do not exactly close the imported cell domain");
    }
}

ImportedMesh2D load_gmsh2_ascii_mesh2d(const std::filesystem::path& path, Gmsh2ImportLimits limits) {
    if (limits.max_physical_names == 0 || limits.max_nodes == 0 || limits.max_elements == 0 ||
        limits.max_tags_per_element == 0) {
        throw std::invalid_argument("Gmsh import limits must be positive");
    }
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open Gmsh mesh: " + path.string());

    ImportedMesh2D mesh;
    bool saw_mesh_format = false;
    bool saw_physical_names = false;
    bool saw_nodes = false;
    bool saw_elements = false;

    std::string section;
    while (std::getline(input, section)) {
        section = trim(section);
        if (section.empty()) continue;

        if (section == "$MeshFormat") {
            if (saw_mesh_format) throw std::runtime_error("duplicate $MeshFormat section");
            const std::string row = read_required_line(input, section);
            std::istringstream stream(row);
            double version = 0.0;
            int file_type = -1;
            int data_size = 0;
            if (!(stream >> version >> file_type >> data_size)) {
                throw std::runtime_error("invalid $MeshFormat row: " + row);
            }
            require_no_trailing_tokens(stream, "$MeshFormat row");
            if (!std::isfinite(version) || version < 2.0 || version >= 3.0) {
                throw std::runtime_error("only Gmsh v2 ASCII meshes are supported");
            }
            if (file_type != 0) {
                throw std::runtime_error("binary Gmsh meshes are not supported");
            }
            if (data_size != 8) throw std::runtime_error("Gmsh mesh must use 8-byte floating-point values");
            saw_mesh_format = true;
            require_section_end(input, section);
        } else if (section == "$PhysicalNames") {
            if (saw_physical_names) throw std::runtime_error("duplicate $PhysicalNames section");
            saw_physical_names = true;
            const std::size_t count = parse_bounded_count(read_required_line(input, section),
                                                          limits.max_physical_names, "$PhysicalNames");
            for (std::size_t i = 0; i < count; ++i) {
                int dimension = 0;
                int tag = 0;
                const std::string name = parse_physical_name_line(read_required_line(input, section), dimension, tag);
                if (dimension == 1 || dimension == 2) mesh.set_physical_name(dimension, tag, name);
            }
            require_section_end(input, section);
        } else if (section == "$Nodes") {
            if (saw_nodes) throw std::runtime_error("duplicate $Nodes section");
            const std::size_t count = parse_bounded_count(read_required_line(input, section),
                                                          limits.max_nodes, "$Nodes");
            for (std::size_t i = 0; i < count; ++i) {
                const std::string row = read_required_line(input, section);
                std::istringstream stream(row);
                std::string id_token;
                double x = 0.0;
                double y = 0.0;
                double z = 0.0;
                if (!(stream >> id_token >> x >> y >> z)) throw std::runtime_error("invalid $Nodes row: " + row);
                require_no_trailing_tokens(stream, "$Nodes row");
                const std::size_t id = parse_positive_id(id_token, "Gmsh node id");
                if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
                    throw std::runtime_error("Gmsh node coordinates must be finite");
                }
                if (std::abs(z) > Z_TOLERANCE) {
                    throw std::runtime_error("2D Gmsh importer requires z=0 nodes; node " + std::to_string(id) + " has z=" + std::to_string(z));
                }
                mesh.add_node(ImportedMeshNode2D{id, Vec2{x, y}});
            }
            saw_nodes = true;
            require_section_end(input, section);
        } else if (section == "$Elements") {
            if (saw_elements) throw std::runtime_error("duplicate $Elements section");
            const std::size_t count = parse_bounded_count(read_required_line(input, section),
                                                          limits.max_elements, "$Elements");
            std::set<std::size_t> element_ids;
            for (std::size_t i = 0; i < count; ++i) {
                const std::string row = read_required_line(input, section);
                std::istringstream stream(row);
                std::string id_token;
                int type = 0;
                int tag_count = 0;
                if (!(stream >> id_token >> type >> tag_count)) {
                    throw std::runtime_error("invalid $Elements row: " + row);
                }
                const std::size_t id = parse_positive_id(id_token, "Gmsh element id");
                if (!element_ids.insert(id).second) {
                    throw std::runtime_error("Gmsh element ids must be positive and unique");
                }
                if (tag_count < 0) throw std::runtime_error("invalid negative Gmsh element tag count: " + row);
                if (static_cast<std::size_t>(tag_count) > limits.max_tags_per_element) {
                    throw std::runtime_error("Gmsh element tag count exceeds configured import limit");
                }
                std::vector<int> tags(static_cast<std::size_t>(tag_count));
                for (auto& tag : tags) {
                    if (!(stream >> tag)) throw std::runtime_error("invalid $Elements tag list: " + row);
                }
                const int physical_tag = tags.empty() ? 0 : tags.front();
                const std::size_t node_count = required_node_count(type);
                std::vector<std::size_t> node_ids(node_count);
                for (auto& node_id : node_ids) {
                    std::string node_id_token;
                    if (!(stream >> node_id_token)) {
                        throw std::runtime_error("invalid $Elements node list: " + row);
                    }
                    node_id = parse_positive_id(node_id_token, "Gmsh element node id");
                }
                if (node_count != 0) require_no_trailing_tokens(stream, "$Elements row");
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
