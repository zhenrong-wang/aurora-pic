#pragma once

#include "pic/Types.hpp"
#include <array>
#include <cstddef>
#include <filesystem>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace pic {

struct ImportedMeshNode2D {
    std::size_t id{0};
    Vec2 position{};
};

enum class ImportedCellShape2D { Triangle, Quadrilateral };

struct ImportedCell2D {
    std::size_t id{0};
    ImportedCellShape2D shape{ImportedCellShape2D::Triangle};
    std::vector<std::size_t> node_ids;
    int physical_tag{0};
    std::string label;
};

struct ImportedBoundaryFace2D {
    std::size_t id{0};
    std::array<std::size_t, 2> node_ids{};
    int physical_tag{0};
    std::string label;
};

struct Gmsh2ImportLimits {
    std::size_t max_physical_names{100000};
    std::size_t max_nodes{10000000};
    std::size_t max_elements{20000000};
    std::size_t max_tags_per_element{64};
};

struct ImportedPointLocation2D {
    std::size_t cell_id{0};
    std::vector<std::size_t> node_ids;
    std::vector<double> shape_weights;
};

struct ImportedMeshQuality2D {
    double minimum_cell_area{0.0};
    double maximum_cell_area{0.0};
    double minimum_edge_length{0.0};
    double maximum_edge_length{0.0};
    double minimum_corner_angle_degrees{0.0};
    double maximum_cell_edge_ratio{0.0};
};

class ImportedMesh2D {
public:
    void add_node(ImportedMeshNode2D node);
    void add_cell(ImportedCell2D cell);
    void add_boundary_face(ImportedBoundaryFace2D face);
    void set_physical_name(int dimension, int tag, std::string name);

    const std::vector<ImportedMeshNode2D>& nodes() const { return nodes_; }
    const std::vector<ImportedCell2D>& cells() const { return cells_; }
    const std::vector<ImportedBoundaryFace2D>& boundary_faces() const { return boundary_faces_; }
    const std::map<std::pair<int, int>, std::string>& physical_names() const { return physical_names_; }

    const ImportedMeshNode2D& node_by_id(std::size_t id) const;
    const ImportedCell2D& cell_by_id(std::size_t id) const;
    const ImportedBoundaryFace2D& boundary_face_by_id(std::size_t id) const;
    std::string label_for_physical_tag(int dimension, int tag) const;
    std::vector<std::string> boundary_labels() const;
    Vec2 min_corner() const;
    Vec2 max_corner() const;
    double cell_area(std::size_t id) const;
    Vec2 cell_centroid(std::size_t id) const;
    double boundary_face_length(std::size_t id) const;
    double total_area() const;
    ImportedMeshQuality2D quality() const;
    std::optional<ImportedPointLocation2D> cell_coordinates(std::size_t cell_id, Vec2 point,
                                                            double relative_tolerance = 1e-12) const;
    std::optional<ImportedPointLocation2D> locate_point(Vec2 point,
                                                        double relative_tolerance = 1e-12) const;
    void validate() const;

private:
    std::vector<ImportedMeshNode2D> nodes_;
    std::vector<ImportedCell2D> cells_;
    std::vector<ImportedBoundaryFace2D> boundary_faces_;
    std::map<std::pair<int, int>, std::string> physical_names_;
    std::map<std::size_t, std::size_t> node_indices_;
    std::map<std::size_t, std::size_t> cell_indices_;
    std::map<std::size_t, std::size_t> boundary_face_indices_;
    Vec2 minimum_corner_{};
    Vec2 maximum_corner_{};
    bool has_bounds_{false};
};

ImportedMesh2D load_gmsh2_ascii_mesh2d(const std::filesystem::path& path,
                                       Gmsh2ImportLimits limits = {});

} // namespace pic
