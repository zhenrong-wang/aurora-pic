#pragma once

#include "pic/Types.hpp"
#include <filesystem>
#include <string>
#include <vector>

namespace pic {

enum class CoordinateAxis { X, Y, Z };

CoordinateAxis parse_coordinate_axis(const std::string& value);
std::string to_string(CoordinateAxis axis);

class TabulatedVectorField1D {
public:
    TabulatedVectorField1D(
        CoordinateAxis axis,
        std::vector<double> coordinates,
        std::vector<Vec3> values);

    CoordinateAxis axis() const { return axis_; }
    const std::vector<double>& coordinates() const { return coordinates_; }
    const std::vector<Vec3>& values() const { return values_; }

    Vec3 evaluate(Vec3 position) const;
    void validate_domain(Vec3 minimum, Vec3 maximum,
                         const std::string& context) const;

private:
    CoordinateAxis axis_{CoordinateAxis::X};
    std::vector<double> coordinates_;
    std::vector<Vec3> values_;
};

TabulatedVectorField1D load_tabulated_vector_field_1d(
    const std::filesystem::path& path,
    CoordinateAxis axis);

} // namespace pic
