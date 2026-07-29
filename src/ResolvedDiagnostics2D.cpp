#include "pic/ResolvedDiagnostics2D.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <utility>

namespace pic {
namespace {

constexpr double TWO_PI = 6.283185307179586476925286766559;

double axis_length(const Mesh2D& mesh, CoordinateAxis axis) {
    if (axis == CoordinateAxis::X) return mesh.length_x();
    if (axis == CoordinateAxis::Y) return mesh.length_y();
    throw std::invalid_argument(
        "resolved 2D diagnostics axis must be x or y");
}

std::size_t axis_nodes(const Mesh2D& mesh, CoordinateAxis axis) {
    if (axis == CoordinateAxis::X) return mesh.nx();
    if (axis == CoordinateAxis::Y) return mesh.ny();
    throw std::invalid_argument(
        "resolved 2D diagnostics axis must be x or y");
}

double axis_spacing(const Mesh2D& mesh, CoordinateAxis axis) {
    return axis == CoordinateAxis::X ? mesh.dx() : mesh.dy();
}

Boundary axis_boundary(const Mesh2D& mesh, CoordinateAxis axis) {
    return axis == CoordinateAxis::X
        ? mesh.boundary_x()
        : mesh.boundary_y();
}

double coordinate(
    const Particle2D& particle, CoordinateAxis axis) {
    return axis == CoordinateAxis::X
        ? particle.position.x
        : particle.position.y;
}

double node_coordinate(
    const Mesh2D& mesh, CoordinateAxis axis, std::size_t node) {
    return axis == CoordinateAxis::X
        ? mesh.node_x(node)
        : mesh.node_y(node);
}

double node_width(
    const Mesh2D& mesh, CoordinateAxis axis, std::size_t node) {
    const double spacing = axis_spacing(mesh, axis);
    if (axis_boundary(mesh, axis) == Boundary::Periodic) {
        return spacing;
    }
    return node == 0 || node + 1 == axis_nodes(mesh, axis)
        ? 0.5 * spacing
        : spacing;
}

std::pair<std::size_t, std::size_t> node_indices(
    CoordinateAxis profile_axis,
    std::size_t profile_node,
    std::size_t transverse_node) {
    return profile_axis == CoordinateAxis::X
        ? std::pair{profile_node, transverse_node}
        : std::pair{transverse_node, profile_node};
}

struct DepositionWeights1D {
    std::size_t lower{0};
    std::size_t upper{0};
    double lower_weight{1.0};
    double upper_weight{0.0};
};

DepositionWeights1D profile_weights(
    const Mesh2D& mesh,
    CoordinateAxis axis,
    double position) {
    const double length = axis_length(mesh, axis);
    const double spacing = axis_spacing(mesh, axis);
    const std::size_t nodes = axis_nodes(mesh, axis);
    const bool periodic =
        axis_boundary(mesh, axis) == Boundary::Periodic;
    double resolved = position;
    if (periodic) {
        resolved = std::fmod(
            std::fmod(position, length) + length, length);
    } else {
        resolved = std::clamp(position, 0.0, length);
    }
    const double grid_coordinate = resolved / spacing;
    std::size_t lower =
        static_cast<std::size_t>(std::floor(grid_coordinate));
    double fraction =
        grid_coordinate - static_cast<double>(lower);
    if (periodic) {
        return {
            lower % nodes,
            (lower + 1) % nodes,
            1.0 - fraction,
            fraction};
    }
    lower = std::min(lower, nodes - 2);
    fraction = std::clamp(
        grid_coordinate - static_cast<double>(lower),
        0.0, 1.0);
    return {lower, lower + 1, 1.0 - fraction, fraction};
}

struct RawSpeciesPoint {
    double macro_particles{0.0};
    double represented_number{0.0};
    double velocity_x{0.0};
    double velocity_y{0.0};
    double velocity_z{0.0};
    double velocity2_x{0.0};
    double velocity2_y{0.0};
    double velocity2_z{0.0};
};

void add_particle(
    RawSpeciesPoint& point,
    const Particle2D& particle,
    double macro_weight,
    double shape_weight) {
    point.macro_particles += shape_weight;
    const double represented = macro_weight * shape_weight;
    point.represented_number += represented;
    point.velocity_x += represented * particle.velocity.x;
    point.velocity_y += represented * particle.velocity.y;
    point.velocity_z += represented * particle.velocity_z;
    point.velocity2_x +=
        represented * particle.velocity.x * particle.velocity.x;
    point.velocity2_y +=
        represented * particle.velocity.y * particle.velocity.y;
    point.velocity2_z +=
        represented * particle.velocity_z * particle.velocity_z;
}

double nonnegative_variance(
    double first_moment,
    double second_moment,
    double density) {
    if (!(density > 0.0)) return 0.0;
    const double mean = first_moment / density;
    return std::max(0.0, second_moment / density - mean * mean);
}

ResolvedSpeciesProfilePoint2D resolve_species_point(
    const Species2D& species,
    double coordinate_value,
    double volume,
    const RawSpeciesPoint& raw,
    UnitSystem unit_system) {
    ResolvedSpeciesProfilePoint2D result;
    result.species = species.name();
    result.coordinate = coordinate_value;
    result.macro_particle_equivalent = raw.macro_particles;
    result.represented_number = raw.represented_number;
    if (unit_system != UnitSystem::SI) {
        result.temperature_ev =
            std::numeric_limits<double>::quiet_NaN();
    }
    if (!(raw.represented_number > 0.0)) return result;
    result.number_density = raw.represented_number / volume;
    result.density_velocity_x = raw.velocity_x / volume;
    result.density_velocity_y = raw.velocity_y / volume;
    result.density_velocity_z = raw.velocity_z / volume;
    result.density_velocity2_x = raw.velocity2_x / volume;
    result.density_velocity2_y = raw.velocity2_y / volume;
    result.density_velocity2_z = raw.velocity2_z / volume;
    result.mean_velocity_x =
        result.density_velocity_x / result.number_density;
    result.mean_velocity_y =
        result.density_velocity_y / result.number_density;
    result.mean_velocity_z =
        result.density_velocity_z / result.number_density;
    const double variance_x = nonnegative_variance(
        result.density_velocity_x,
        result.density_velocity2_x,
        result.number_density);
    const double variance_y = nonnegative_variance(
        result.density_velocity_y,
        result.density_velocity2_y,
        result.number_density);
    const double variance_z = nonnegative_variance(
        result.density_velocity_z,
        result.density_velocity2_z,
        result.number_density);
    result.thermal_speed_x = std::sqrt(variance_x);
    result.thermal_speed_y = std::sqrt(variance_y);
    result.thermal_speed_z = std::sqrt(variance_z);
    if (unit_system == UnitSystem::SI) {
        result.temperature_ev =
            species.mass() *
            (variance_x + variance_y + variance_z) /
            (3.0 * ELEMENTARY_CHARGE_SI);
    } else {
        result.temperature_ev =
            std::numeric_limits<double>::quiet_NaN();
    }
    result.current_density_x =
        species.charge() * result.density_velocity_x;
    result.current_density_y =
        species.charge() * result.density_velocity_y;
    result.current_density_z =
        species.charge() * result.density_velocity_z;
    return result;
}

std::complex<double> phase_factor(
    std::size_t mode, double value, double length) {
    const double phase =
        -TWO_PI * static_cast<double>(mode) * value / length;
    return {std::cos(phase), std::sin(phase)};
}

double mode_amplitude(
    std::complex<double> coefficient,
    std::size_t mode,
    std::size_t node_count) {
    const double factor =
        mode == 0 || 2 * mode == node_count ? 1.0 : 2.0;
    return factor * std::abs(coefficient);
}

void append_mode(
    ResolvedDiagnosticSnapshot2D& snapshot,
    std::size_t mode,
    double wavenumber,
    std::size_t node_count,
    std::string quantity,
    std::string species,
    std::complex<double> coefficient) {
    snapshot.modes.push_back({
        mode,
        wavenumber,
        std::move(quantity),
        std::move(species),
        coefficient.real(),
        coefficient.imag(),
        mode_amplitude(coefficient, mode, node_count)});
}

void require_finite_snapshot(
    const ResolvedDiagnosticSnapshot2D& snapshot) {
    const auto finite = [](double value) {
        return std::isfinite(value);
    };
    for (const auto& field : snapshot.fields) {
        if (!finite(field.coordinate) ||
            !finite(field.potential) ||
            !finite(field.electric_x) ||
            !finite(field.electric_y) ||
            !finite(field.charge_density)) {
            throw std::runtime_error(
                "resolved 2D field profile is not finite");
        }
    }
    for (const auto& point : snapshot.species) {
        if (!finite(point.coordinate) ||
            !finite(point.macro_particle_equivalent) ||
            !finite(point.represented_number) ||
            !finite(point.number_density) ||
            !finite(point.mean_velocity_x) ||
            !finite(point.mean_velocity_y) ||
            !finite(point.mean_velocity_z) ||
            !finite(point.thermal_speed_x) ||
            !finite(point.thermal_speed_y) ||
            !finite(point.thermal_speed_z) ||
            (std::isnan(point.temperature_ev)
                 ? false
                 : !finite(point.temperature_ev)) ||
            !finite(point.current_density_x) ||
            !finite(point.current_density_y) ||
            !finite(point.current_density_z)) {
            throw std::runtime_error(
                "resolved 2D species profile is not finite");
        }
    }
    for (const auto& mode : snapshot.modes) {
        if (!finite(mode.wavenumber) ||
            !finite(mode.real) ||
            !finite(mode.imaginary) ||
            !finite(mode.amplitude)) {
            throw std::runtime_error(
                "resolved 2D mode coefficient is not finite");
        }
    }
}

} // namespace

ResolvedDiagnosticSnapshot2D compute_resolved_diagnostics_2d(
    std::size_t step,
    double time,
    const Mesh2D& mesh,
    const std::vector<Species2D>& species,
    const ResolvedDiagnostics2DConfig& config,
    UnitSystem unit_system,
    double out_of_plane_depth) {
    if (!config.enabled) {
        throw std::invalid_argument(
            "cannot compute disabled resolved 2D diagnostics");
    }
    if (!std::isfinite(time) || time < 0.0) {
        throw std::invalid_argument(
            "resolved 2D diagnostic time must be non-negative and finite");
    }
    if (!std::isfinite(out_of_plane_depth) ||
        !(out_of_plane_depth > 0.0)) {
        throw std::invalid_argument(
            "resolved 2D diagnostic depth must be positive and finite");
    }
    if (config.profile_axis == CoordinateAxis::Z ||
        config.mode_axis == CoordinateAxis::Z ||
        config.profile_axis == config.mode_axis) {
        throw std::invalid_argument(
            "resolved 2D diagnostics require distinct x/y profile and mode axes");
    }
    if (axis_boundary(mesh, config.mode_axis) !=
        Boundary::Periodic) {
        throw std::invalid_argument(
            "resolved 2D mode axis must be periodic");
    }
    const std::size_t mode_nodes =
        axis_nodes(mesh, config.mode_axis);
    if (config.max_mode > mode_nodes / 2) {
        throw std::invalid_argument(
            "resolved 2D max_mode exceeds the mode-axis Nyquist limit");
    }

    ResolvedDiagnosticSnapshot2D snapshot;
    snapshot.step = step;
    snapshot.time = time;

    const std::size_t profile_nodes =
        axis_nodes(mesh, config.profile_axis);
    const std::size_t transverse_nodes =
        axis_nodes(mesh, config.mode_axis);
    const double transverse_length =
        axis_length(mesh, config.mode_axis);
    snapshot.fields.resize(profile_nodes);
    for (std::size_t profile_node = 0;
         profile_node < profile_nodes; ++profile_node) {
        auto& output = snapshot.fields[profile_node];
        output.coordinate =
            node_coordinate(
                mesh, config.profile_axis, profile_node);
        for (std::size_t transverse_node = 0;
             transverse_node < transverse_nodes;
             ++transverse_node) {
            const auto [i, j] = node_indices(
                config.profile_axis,
                profile_node, transverse_node);
            const std::size_t index = mesh.index(i, j);
            const double weight =
                node_width(
                    mesh, config.mode_axis,
                    transverse_node) /
                transverse_length;
            output.potential += mesh.phi()[index] * weight;
            output.electric_x +=
                mesh.electric_x()[index] * weight;
            output.electric_y +=
                mesh.electric_y()[index] * weight;
            output.charge_density +=
                mesh.rho()[index] * weight;
        }
    }

    snapshot.species.reserve(species.size() * profile_nodes);
    for (const auto& current_species : species) {
        std::vector<RawSpeciesPoint> raw(profile_nodes);
        for (const auto& particle : current_species.particles()) {
            if (!particle.alive) continue;
            const auto weights = profile_weights(
                mesh, config.profile_axis,
                coordinate(particle, config.profile_axis));
            add_particle(
                raw[weights.lower], particle,
                current_species.weight(),
                weights.lower_weight);
            add_particle(
                raw[weights.upper], particle,
                current_species.weight(),
                weights.upper_weight);
        }
        for (std::size_t profile_node = 0;
             profile_node < profile_nodes; ++profile_node) {
            const double volume =
                node_width(
                    mesh, config.profile_axis,
                    profile_node) *
                transverse_length * out_of_plane_depth;
            snapshot.species.push_back(
                resolve_species_point(
                    current_species,
                    node_coordinate(
                        mesh, config.profile_axis,
                        profile_node),
                    volume, raw[profile_node], unit_system));
        }
    }

    const double area =
        mesh.length_x() * mesh.length_y();
    const double mode_length =
        axis_length(mesh, config.mode_axis);
    struct SpeciesModeGrid {
        std::vector<double> number;
        std::vector<double> current_x;
        std::vector<double> current_y;
        std::vector<double> current_z;
    };
    std::vector<SpeciesModeGrid> species_mode_grids;
    species_mode_grids.reserve(species.size());
    for (const auto& current_species : species) {
        SpeciesModeGrid grid{
            std::vector<double>(mode_nodes, 0.0),
            std::vector<double>(mode_nodes, 0.0),
            std::vector<double>(mode_nodes, 0.0),
            std::vector<double>(mode_nodes, 0.0)};
        for (const auto& particle :
             current_species.particles()) {
            if (!particle.alive) continue;
            const auto weights = profile_weights(
                mesh, config.mode_axis,
                coordinate(particle, config.mode_axis));
            const auto deposit =
                [&](std::size_t node, double shape) {
                    const double represented =
                        current_species.weight() * shape;
                    const double charge_weight =
                        current_species.charge() * represented;
                    grid.number[node] += represented;
                    grid.current_x[node] +=
                        charge_weight * particle.velocity.x;
                    grid.current_y[node] +=
                        charge_weight * particle.velocity.y;
                    grid.current_z[node] +=
                        charge_weight * particle.velocity_z;
                };
            deposit(weights.lower, weights.lower_weight);
            deposit(weights.upper, weights.upper_weight);
        }
        species_mode_grids.push_back(std::move(grid));
    }
    const double volume =
        area * out_of_plane_depth;
    for (std::size_t mode = 0;
         mode <= config.max_mode; ++mode) {
        const double wavenumber =
            TWO_PI * static_cast<double>(mode) / mode_length;
        std::complex<double> rho{};
        std::complex<double> electric_x{};
        std::complex<double> electric_y{};
        for (std::size_t j = 0; j < mesh.ny(); ++j) {
            for (std::size_t i = 0; i < mesh.nx(); ++i) {
                const double mode_coordinate =
                    config.mode_axis == CoordinateAxis::X
                        ? mesh.node_x(i)
                        : mesh.node_y(j);
                const auto phase =
                    phase_factor(
                        mode, mode_coordinate, mode_length);
                const double weight =
                    mesh.node_area(i, j) / area;
                const std::size_t index = mesh.index(i, j);
                rho += mesh.rho()[index] * weight * phase;
                electric_x +=
                    mesh.electric_x()[index] * weight * phase;
                electric_y +=
                    mesh.electric_y()[index] * weight * phase;
            }
        }
        append_mode(
            snapshot, mode, wavenumber, mode_nodes,
            "charge_density", "", rho);
        append_mode(
            snapshot, mode, wavenumber, mode_nodes,
            "electric_x", "", electric_x);
        append_mode(
            snapshot, mode, wavenumber, mode_nodes,
            "electric_y", "", electric_y);

        for (std::size_t species_id = 0;
             species_id < species.size(); ++species_id) {
            const auto& current_species = species[species_id];
            const auto& grid = species_mode_grids[species_id];
            std::complex<double> number_density{};
            std::complex<double> current_x{};
            std::complex<double> current_y{};
            std::complex<double> current_z{};
            for (std::size_t node = 0;
                 node < mode_nodes; ++node) {
                const auto phase = phase_factor(
                    mode,
                    node_coordinate(
                        mesh, config.mode_axis, node),
                    mode_length);
                number_density +=
                    grid.number[node] / volume * phase;
                current_x +=
                    grid.current_x[node] / volume * phase;
                current_y +=
                    grid.current_y[node] / volume * phase;
                current_z +=
                    grid.current_z[node] / volume * phase;
            }
            append_mode(
                snapshot, mode, wavenumber, mode_nodes,
                "number_density", current_species.name(),
                number_density);
            append_mode(
                snapshot, mode, wavenumber, mode_nodes,
                "current_x", current_species.name(), current_x);
            append_mode(
                snapshot, mode, wavenumber, mode_nodes,
                "current_y", current_species.name(), current_y);
            append_mode(
                snapshot, mode, wavenumber, mode_nodes,
                "current_z", current_species.name(), current_z);
        }
    }
    require_finite_snapshot(snapshot);
    return snapshot;
}

struct ResolvedDiagnostics2D::Impl {
    std::filesystem::path output_dir;
    ResolvedDiagnostics2DConfig config;
    UnitSystem unit_system{UnitSystem::Normalized};
    double out_of_plane_depth{1.0};
    std::vector<double> masses;
    std::vector<double> charges;
    std::ofstream field_output;
    std::ofstream species_output;
    std::ofstream mode_output;
    std::size_t samples{0};
    double start_time{0.0};
    double last_time{0.0};
    ResolvedDiagnosticSnapshot2D previous;
    ResolvedDiagnosticSnapshot2D integral;

    Impl(
        std::filesystem::path directory,
        ResolvedDiagnostics2DConfig requested,
        const Mesh2D& mesh,
        const std::vector<Species2D>& species,
        UnitSystem units,
        double depth)
        : output_dir(std::move(directory)),
          config(requested),
          unit_system(units),
          out_of_plane_depth(depth) {
        (void)compute_resolved_diagnostics_2d(
            0, 0.0, mesh, species, config,
            unit_system, out_of_plane_depth);
        masses.reserve(species.size());
        charges.reserve(species.size());
        for (const auto& current_species : species) {
            masses.push_back(current_species.mass());
            charges.push_back(current_species.charge());
        }
        std::filesystem::create_directories(output_dir);
        field_output.open(
            output_dir / "resolved_field_profiles.csv");
        species_output.open(
            output_dir / "resolved_species_profiles.csv");
        mode_output.open(
            output_dir / "resolved_modes.csv");
        if (!field_output || !species_output || !mode_output) {
            throw std::runtime_error(
                "cannot open resolved 2D diagnostic output");
        }
        field_output
            << "step,time,profile_axis,coordinate,potential,"
               "electric_x,electric_y,charge_density\n";
        species_output
            << "step,time,profile_axis,coordinate,species,"
               "macro_particle_equivalent,represented_number,"
               "number_density,mean_velocity_x,mean_velocity_y,"
               "mean_velocity_z,thermal_speed_x,thermal_speed_y,"
               "thermal_speed_z,temperature_ev,current_density_x,"
               "current_density_y,current_density_z\n";
        mode_output
            << "step,time,mode_axis,mode,wavenumber,quantity,"
               "species,real,imaginary,amplitude\n";
    }

    void write_snapshot(
        const ResolvedDiagnosticSnapshot2D& snapshot) {
        for (const auto& point : snapshot.fields) {
            field_output
                << snapshot.step << ',' << std::setprecision(17)
                << snapshot.time << ','
                << to_string(config.profile_axis) << ','
                << point.coordinate << ',' << point.potential
                << ',' << point.electric_x << ','
                << point.electric_y << ','
                << point.charge_density << '\n';
        }
        for (const auto& point : snapshot.species) {
            species_output
                << snapshot.step << ',' << std::setprecision(17)
                << snapshot.time << ','
                << to_string(config.profile_axis) << ','
                << point.coordinate << ',' << point.species
                << ',' << point.macro_particle_equivalent
                << ',' << point.represented_number
                << ',' << point.number_density
                << ',' << point.mean_velocity_x
                << ',' << point.mean_velocity_y
                << ',' << point.mean_velocity_z
                << ',' << point.thermal_speed_x
                << ',' << point.thermal_speed_y
                << ',' << point.thermal_speed_z
                << ',' << point.temperature_ev
                << ',' << point.current_density_x
                << ',' << point.current_density_y
                << ',' << point.current_density_z << '\n';
        }
        for (const auto& mode : snapshot.modes) {
            mode_output
                << snapshot.step << ',' << std::setprecision(17)
                << snapshot.time << ','
                << to_string(config.mode_axis) << ','
                << mode.mode << ',' << mode.wavenumber
                << ',' << mode.quantity << ',' << mode.species
                << ',' << mode.real << ',' << mode.imaginary
                << ',' << mode.amplitude << '\n';
        }
        field_output.flush();
        species_output.flush();
        mode_output.flush();
    }

    static void initialize_integral(
        ResolvedDiagnosticSnapshot2D& target,
        const ResolvedDiagnosticSnapshot2D& source) {
        target = source;
        for (auto& field : target.fields) {
            field.potential = 0.0;
            field.electric_x = 0.0;
            field.electric_y = 0.0;
            field.charge_density = 0.0;
        }
        for (auto& point : target.species) {
            point.macro_particle_equivalent = 0.0;
            point.represented_number = 0.0;
            point.number_density = 0.0;
            point.mean_velocity_x = 0.0;
            point.mean_velocity_y = 0.0;
            point.mean_velocity_z = 0.0;
            point.thermal_speed_x = 0.0;
            point.thermal_speed_y = 0.0;
            point.thermal_speed_z = 0.0;
            point.temperature_ev = 0.0;
            point.current_density_x = 0.0;
            point.current_density_y = 0.0;
            point.current_density_z = 0.0;
            point.density_velocity_x = 0.0;
            point.density_velocity_y = 0.0;
            point.density_velocity_z = 0.0;
            point.density_velocity2_x = 0.0;
            point.density_velocity2_y = 0.0;
            point.density_velocity2_z = 0.0;
        }
        target.modes.clear();
    }

    void accumulate(
        const ResolvedDiagnosticSnapshot2D& snapshot) {
        if (samples == 0) {
            samples = 1;
            start_time = snapshot.time;
            last_time = snapshot.time;
            previous = snapshot;
            initialize_integral(integral, snapshot);
            return;
        }
        if (!(snapshot.time > last_time)) {
            throw std::runtime_error(
                "resolved 2D diagnostic sample times must increase");
        }
        if (snapshot.fields.size() != previous.fields.size() ||
            snapshot.species.size() != previous.species.size()) {
            throw std::runtime_error(
                "resolved 2D diagnostic shape changed during a run");
        }
        const double half_dt =
            0.5 * (snapshot.time - last_time);
        for (std::size_t index = 0;
             index < snapshot.fields.size(); ++index) {
            auto& total = integral.fields[index];
            const auto& before = previous.fields[index];
            const auto& after = snapshot.fields[index];
            total.potential +=
                half_dt * (before.potential + after.potential);
            total.electric_x +=
                half_dt *
                (before.electric_x + after.electric_x);
            total.electric_y +=
                half_dt *
                (before.electric_y + after.electric_y);
            total.charge_density +=
                half_dt *
                (before.charge_density + after.charge_density);
        }
        for (std::size_t index = 0;
             index < snapshot.species.size(); ++index) {
            auto& total = integral.species[index];
            const auto& before = previous.species[index];
            const auto& after = snapshot.species[index];
            const auto integrate =
                [&](double& destination,
                    double first, double second) {
                    destination += half_dt * (first + second);
                };
            integrate(
                total.macro_particle_equivalent,
                before.macro_particle_equivalent,
                after.macro_particle_equivalent);
            integrate(
                total.represented_number,
                before.represented_number,
                after.represented_number);
            integrate(
                total.number_density,
                before.number_density,
                after.number_density);
            integrate(
                total.density_velocity_x,
                before.density_velocity_x,
                after.density_velocity_x);
            integrate(
                total.density_velocity_y,
                before.density_velocity_y,
                after.density_velocity_y);
            integrate(
                total.density_velocity_z,
                before.density_velocity_z,
                after.density_velocity_z);
            integrate(
                total.density_velocity2_x,
                before.density_velocity2_x,
                after.density_velocity2_x);
            integrate(
                total.density_velocity2_y,
                before.density_velocity2_y,
                after.density_velocity2_y);
            integrate(
                total.density_velocity2_z,
                before.density_velocity2_z,
                after.density_velocity2_z);
        }
        ++samples;
        last_time = snapshot.time;
        previous = snapshot;
    }
};

ResolvedDiagnostics2D::ResolvedDiagnostics2D(
    std::filesystem::path output_dir,
    ResolvedDiagnostics2DConfig config,
    const Mesh2D& mesh,
    const std::vector<Species2D>& species,
    UnitSystem unit_system,
    double out_of_plane_depth)
    : impl_(std::make_unique<Impl>(
          std::move(output_dir), config, mesh, species,
          unit_system, out_of_plane_depth)) {}

ResolvedDiagnostics2D::~ResolvedDiagnostics2D() = default;
ResolvedDiagnostics2D::ResolvedDiagnostics2D(
    ResolvedDiagnostics2D&&) noexcept = default;
ResolvedDiagnostics2D& ResolvedDiagnostics2D::operator=(
    ResolvedDiagnostics2D&&) noexcept = default;

ResolvedDiagnosticSnapshot2D ResolvedDiagnostics2D::sample(
    std::size_t step,
    double time,
    const Mesh2D& mesh,
    const std::vector<Species2D>& species) {
    auto snapshot = compute_resolved_diagnostics_2d(
        step, time, mesh, species, impl_->config,
        impl_->unit_system, impl_->out_of_plane_depth);
    impl_->write_snapshot(snapshot);
    impl_->accumulate(snapshot);
    return snapshot;
}

void ResolvedDiagnostics2D::write_time_averages() {
    if (impl_->samples == 0) return;
    const double duration = impl_->last_time - impl_->start_time;
    ResolvedDiagnosticSnapshot2D average =
        duration > 0.0 ? impl_->integral : impl_->previous;
    if (duration > 0.0) {
        for (auto& field : average.fields) {
            field.potential /= duration;
            field.electric_x /= duration;
            field.electric_y /= duration;
            field.charge_density /= duration;
        }
        for (auto& point : average.species) {
            point.macro_particle_equivalent /= duration;
            point.represented_number /= duration;
            point.number_density /= duration;
            point.density_velocity_x /= duration;
            point.density_velocity_y /= duration;
            point.density_velocity_z /= duration;
            point.density_velocity2_x /= duration;
            point.density_velocity2_y /= duration;
            point.density_velocity2_z /= duration;
        }
    }
    std::ofstream field_average(
        impl_->output_dir /
        "resolved_field_time_average.csv");
    std::ofstream species_average(
        impl_->output_dir /
        "resolved_species_time_average.csv");
    if (!field_average || !species_average) {
        throw std::runtime_error(
            "cannot open resolved 2D time-average output");
    }
    field_average
        << "start_time,end_time,duration,samples,profile_axis,"
           "coordinate,potential,electric_x,electric_y,"
           "charge_density\n";
    species_average
        << "start_time,end_time,duration,samples,profile_axis,"
           "coordinate,species,macro_particle_equivalent,"
           "represented_number,number_density,mean_velocity_x,"
           "mean_velocity_y,mean_velocity_z,thermal_speed_x,"
           "thermal_speed_y,thermal_speed_z,temperature_ev,"
           "current_density_x,current_density_y,"
           "current_density_z\n";
    for (const auto& field : average.fields) {
        field_average
            << std::setprecision(17)
            << impl_->start_time << ',' << impl_->last_time
            << ',' << duration << ',' << impl_->samples
            << ',' << to_string(impl_->config.profile_axis)
            << ',' << field.coordinate << ',' << field.potential
            << ',' << field.electric_x << ','
            << field.electric_y << ','
            << field.charge_density << '\n';
    }
    const std::size_t points_per_species =
        average.fields.size();
    for (std::size_t index = 0;
         index < average.species.size(); ++index) {
        auto point = average.species[index];
        const std::size_t species_id =
            index / points_per_species;
        const double density = point.number_density;
        if (density > 0.0) {
            point.mean_velocity_x =
                point.density_velocity_x / density;
            point.mean_velocity_y =
                point.density_velocity_y / density;
            point.mean_velocity_z =
                point.density_velocity_z / density;
            const double variance_x = nonnegative_variance(
                point.density_velocity_x,
                point.density_velocity2_x, density);
            const double variance_y = nonnegative_variance(
                point.density_velocity_y,
                point.density_velocity2_y, density);
            const double variance_z = nonnegative_variance(
                point.density_velocity_z,
                point.density_velocity2_z, density);
            point.thermal_speed_x = std::sqrt(variance_x);
            point.thermal_speed_y = std::sqrt(variance_y);
            point.thermal_speed_z = std::sqrt(variance_z);
            point.temperature_ev =
                impl_->unit_system == UnitSystem::SI
                    ? impl_->masses.at(species_id) *
                          (variance_x + variance_y + variance_z) /
                          (3.0 * ELEMENTARY_CHARGE_SI)
                    : std::numeric_limits<double>::quiet_NaN();
            point.current_density_x =
                impl_->charges.at(species_id) *
                point.density_velocity_x;
            point.current_density_y =
                impl_->charges.at(species_id) *
                point.density_velocity_y;
            point.current_density_z =
                impl_->charges.at(species_id) *
                point.density_velocity_z;
        }
        species_average
            << std::setprecision(17)
            << impl_->start_time << ',' << impl_->last_time
            << ',' << duration << ',' << impl_->samples
            << ',' << to_string(impl_->config.profile_axis)
            << ',' << point.coordinate << ',' << point.species
            << ',' << point.macro_particle_equivalent
            << ',' << point.represented_number
            << ',' << point.number_density
            << ',' << point.mean_velocity_x
            << ',' << point.mean_velocity_y
            << ',' << point.mean_velocity_z
            << ',' << point.thermal_speed_x
            << ',' << point.thermal_speed_y
            << ',' << point.thermal_speed_z
            << ',' << point.temperature_ev
            << ',' << point.current_density_x
            << ',' << point.current_density_y
            << ',' << point.current_density_z << '\n';
    }
}

std::size_t ResolvedDiagnostics2D::sample_count() const {
    return impl_->samples;
}

} // namespace pic
