#include "pic/Units.hpp"

#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdexcept>

namespace pic {

double maxwellian_thermal_velocity_from_ev(
    double temperature_ev, double mass_kg) {
    if (!std::isfinite(temperature_ev) ||
        temperature_ev < 0.0) {
        throw std::invalid_argument(
            "temperature_ev must be non-negative and finite");
    }
    if (!std::isfinite(mass_kg) || !(mass_kg > 0.0)) {
        throw std::invalid_argument(
            "SI particle mass must be positive and finite");
    }
    const double velocity = std::sqrt(
        temperature_ev * ELEMENTARY_CHARGE_SI / mass_kg);
    if (!std::isfinite(velocity)) {
        throw std::overflow_error(
            "temperature_ev to thermal-velocity conversion overflow");
    }
    return velocity;
}

void write_unit_metadata(
    const std::filesystem::path& output_dir,
    const UnitSystemConfig& units,
    std::size_t spatial_dimension,
    std::optional<double> out_of_plane_depth) {
    if (spatial_dimension < 1 || spatial_dimension > 3) {
        throw std::invalid_argument(
            "unit metadata spatial dimension must be 1, 2, or 3");
    }
    if (out_of_plane_depth &&
        (!std::isfinite(*out_of_plane_depth) ||
         !(*out_of_plane_depth > 0.0))) {
        throw std::invalid_argument(
            "unit metadata extrusion depth must be positive and finite");
    }
    std::filesystem::create_directories(output_dir);
    std::ofstream output(output_dir / "units.txt");
    if (!output) {
        throw std::runtime_error("cannot open unit metadata output");
    }
    output << std::setprecision(17)
           << "unit_system " << to_string(units.system) << '\n'
           << "relative_permittivity "
           << units.relative_permittivity << '\n'
           << "permittivity " << units.permittivity() << '\n'
           << "spatial_dimension " << spatial_dimension << '\n';
    if (spatial_dimension == 2 && out_of_plane_depth) {
        output << "out_of_plane_depth "
               << *out_of_plane_depth << '\n';
    }
    if (units.system == UnitSystem::SI) {
        output << "length m\n"
               << "time s\n"
               << "velocity m/s\n"
               << "potential V\n"
               << "electric_field V/m\n"
               << "magnetic_field T\n"
               << "charge C\n"
               << "mass kg\n";
        if (spatial_dimension == 1) {
            output << "macro_weight particles/m^2\n"
                   << "energy J/m^2\n";
        } else if (spatial_dimension == 2) {
            if (out_of_plane_depth) {
                output << "macro_weight particles\n"
                       << "energy J\n";
            } else {
                output << "macro_weight particles/m\n"
                       << "energy J/m\n";
            }
        } else {
            output << "macro_weight particles\n"
                   << "energy J\n";
        }
    } else {
        output << "quantity_basis dimensionless_consistent_normalization\n"
               << "macro_weight represented_normalized_particles\n"
               << "energy normalized_energy\n";
    }
    if (!output) {
        throw std::runtime_error("failed while writing unit metadata");
    }
}

} // namespace pic
