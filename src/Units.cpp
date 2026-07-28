#include "pic/Units.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace pic {

void write_unit_metadata(
    const std::filesystem::path& output_dir,
    const UnitSystemConfig& units,
    std::size_t spatial_dimension) {
    if (spatial_dimension < 1 || spatial_dimension > 3) {
        throw std::invalid_argument(
            "unit metadata spatial dimension must be 1, 2, or 3");
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
            output << "macro_weight particles/m\n"
                   << "energy J/m\n";
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
