#pragma once

#include "pic/Mesh2D.hpp"
#include "pic/PrescribedField.hpp"
#include "pic/Species2D.hpp"
#include "pic/Types.hpp"

#include <cstddef>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace pic {

struct ResolvedDiagnostics2DConfig {
    bool enabled{false};
    std::size_t interval{0}; // zero inherits output_interval
    std::size_t start_step{0};
    CoordinateAxis profile_axis{CoordinateAxis::X};
    CoordinateAxis mode_axis{CoordinateAxis::Y};
    std::size_t max_mode{8};
};

struct ResolvedFieldProfilePoint2D {
    double coordinate{0.0};
    double potential{0.0};
    double electric_x{0.0};
    double electric_y{0.0};
    double charge_density{0.0};
};

struct ResolvedSpeciesProfilePoint2D {
    std::string species{};
    double coordinate{0.0};
    double macro_particle_equivalent{0.0};
    double represented_number{0.0};
    double number_density{0.0};
    double mean_velocity_x{0.0};
    double mean_velocity_y{0.0};
    double mean_velocity_z{0.0};
    double thermal_speed_x{0.0};
    double thermal_speed_y{0.0};
    double thermal_speed_z{0.0};
    double temperature_ev{0.0};
    double current_density_x{0.0};
    double current_density_y{0.0};
    double current_density_z{0.0};
    double density_velocity_x{0.0};
    double density_velocity_y{0.0};
    double density_velocity_z{0.0};
    double density_velocity2_x{0.0};
    double density_velocity2_y{0.0};
    double density_velocity2_z{0.0};
};

struct ResolvedModeCoefficient2D {
    std::size_t mode{0};
    double wavenumber{0.0};
    std::string quantity{};
    std::string species{};
    double real{0.0};
    double imaginary{0.0};
    double amplitude{0.0};
};

struct ResolvedDiagnosticSnapshot2D {
    std::size_t step{0};
    double time{0.0};
    std::vector<ResolvedFieldProfilePoint2D> fields{};
    std::vector<ResolvedSpeciesProfilePoint2D> species{};
    std::vector<ResolvedModeCoefficient2D> modes{};
};

ResolvedDiagnosticSnapshot2D compute_resolved_diagnostics_2d(
    std::size_t step,
    double time,
    const Mesh2D& mesh,
    const std::vector<Species2D>& species,
    const ResolvedDiagnostics2DConfig& config,
    UnitSystem unit_system,
    double out_of_plane_depth = 1.0);

class ResolvedDiagnostics2D {
public:
    ResolvedDiagnostics2D(
        std::filesystem::path output_dir,
        ResolvedDiagnostics2DConfig config,
        const Mesh2D& mesh,
        const std::vector<Species2D>& species,
        UnitSystem unit_system,
        double out_of_plane_depth = 1.0);
    ~ResolvedDiagnostics2D();

    ResolvedDiagnostics2D(const ResolvedDiagnostics2D&) = delete;
    ResolvedDiagnostics2D& operator=(const ResolvedDiagnostics2D&) = delete;
    ResolvedDiagnostics2D(ResolvedDiagnostics2D&&) noexcept;
    ResolvedDiagnostics2D& operator=(ResolvedDiagnostics2D&&) noexcept;

    ResolvedDiagnosticSnapshot2D sample(
        std::size_t step,
        double time,
        const Mesh2D& mesh,
        const std::vector<Species2D>& species);
    void write_time_averages();
    std::size_t sample_count() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace pic
