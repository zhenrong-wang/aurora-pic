#pragma once

#include "pic/ImportedMesh2D.hpp"
#include "pic/Runtime.hpp"

#include <cstddef>
#include <map>
#include <optional>
#include <vector>

namespace pic {

struct UnstructuredDepositSummary2D {
    std::size_t deposited_particles{0};
    std::size_t outside_particles{0};
    std::size_t location_cache_hits{0};
    std::size_t location_searches{0};
    double deposited_charge{0.0};
};

struct UnstructuredParticleLocation2D {
    ImportedPointLocation2D location;
    bool valid{false};
};

class UnstructuredMesh2D {
public:
    explicit UnstructuredMesh2D(ImportedMesh2D topology);

    const ImportedMesh2D& topology() const { return topology_; }
    std::size_t size() const { return topology_.nodes().size(); }
    std::size_t node_index(std::size_t node_id) const;
    double node_control_area(std::size_t node_id) const;
    std::optional<ImportedPointLocation2D> locate_point(Vec2 point,
                                                        double relative_tolerance = 1e-12) const;

    std::vector<double>& rho() { return rho_; }
    std::vector<double>& phi() { return phi_; }
    std::vector<Vec2>& electric() { return electric_; }
    const std::vector<double>& rho() const { return rho_; }
    const std::vector<double>& phi() const { return phi_; }
    const std::vector<Vec2>& electric() const { return electric_; }
    const std::vector<double>& node_control_areas() const { return node_control_areas_; }

    void clear_charge();

private:
    ImportedMesh2D topology_;
    std::map<std::size_t, std::size_t> node_indices_;
    std::vector<double> node_control_areas_;
    std::vector<double> rho_;
    std::vector<double> phi_;
    std::vector<Vec2> electric_;
    Vec2 minimum_corner_{};
    Vec2 maximum_corner_{};
    std::size_t spatial_bins_x_{0};
    std::size_t spatial_bins_y_{0};
    std::vector<std::vector<std::size_t>> spatial_cell_ids_;
    std::vector<std::size_t> spatial_global_cell_ids_;
};

UnstructuredDepositSummary2D deposit_charge_shape(UnstructuredMesh2D& mesh,
                                                   const std::vector<Particle2D>& particles,
                                                   double charge, double weight);
UnstructuredDepositSummary2D deposit_charge_shape(
    UnstructuredMesh2D& mesh,
    const std::vector<Particle2D>& particles,
    double charge, double weight,
    const RuntimePolicy& runtime);
UnstructuredDepositSummary2D deposit_charge_shape(
    UnstructuredMesh2D& mesh,
    const std::vector<Particle2D>& particles,
    double charge, double weight,
    const RuntimePolicy& runtime,
    std::vector<UnstructuredParticleLocation2D>& locations);
std::optional<Vec2> interpolate_electric(const UnstructuredMesh2D& mesh, Vec2 position);
std::optional<Vec2> interpolate_electric(
    const UnstructuredMesh2D& mesh, Vec2 position,
    UnstructuredParticleLocation2D& location,
    bool* cache_hit = nullptr);

} // namespace pic
