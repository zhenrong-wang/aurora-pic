#include "pic/PrescribedField.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace pic {
namespace {

double coordinate(Vec3 position, CoordinateAxis axis) {
    switch (axis) {
        case CoordinateAxis::X: return position.x;
        case CoordinateAxis::Y: return position.y;
        case CoordinateAxis::Z: return position.z;
    }
    throw std::logic_error("unknown coordinate axis");
}

std::string lower(std::string value) {
    std::transform(
        value.begin(), value.end(), value.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

Vec3 interpolate(Vec3 left, Vec3 right, double fraction) {
    return {
        left.x + fraction * (right.x - left.x),
        left.y + fraction * (right.y - left.y),
        left.z + fraction * (right.z - left.z)};
}

bool finite(Vec3 value) {
    return std::isfinite(value.x) &&
           std::isfinite(value.y) &&
           std::isfinite(value.z);
}

} // namespace

CoordinateAxis parse_coordinate_axis(const std::string& value) {
    const auto normalized = lower(value);
    if (normalized == "x") return CoordinateAxis::X;
    if (normalized == "y") return CoordinateAxis::Y;
    if (normalized == "z") return CoordinateAxis::Z;
    throw std::invalid_argument(
        "coordinate axis must be x, y, or z");
}

std::string to_string(CoordinateAxis axis) {
    switch (axis) {
        case CoordinateAxis::X: return "x";
        case CoordinateAxis::Y: return "y";
        case CoordinateAxis::Z: return "z";
    }
    return "unknown";
}

TabulatedVectorField1D::TabulatedVectorField1D(
    CoordinateAxis axis,
    std::vector<double> coordinates,
    std::vector<Vec3> values)
    : axis_(axis),
      coordinates_(std::move(coordinates)),
      values_(std::move(values)) {
    if (coordinates_.size() != values_.size()) {
        throw std::invalid_argument(
            "tabulated vector field coordinate/value sizes differ");
    }
    if (coordinates_.size() < 2) {
        throw std::invalid_argument(
            "tabulated vector field requires at least two rows");
    }
    for (std::size_t i = 0; i < coordinates_.size(); ++i) {
        if (!std::isfinite(coordinates_[i]) || !finite(values_[i])) {
            throw std::invalid_argument(
                "tabulated vector field values must be finite");
        }
        if (i != 0 && !(coordinates_[i] > coordinates_[i - 1])) {
            throw std::invalid_argument(
                "tabulated vector field coordinates must be strictly increasing");
        }
    }
}

Vec3 TabulatedVectorField1D::evaluate(Vec3 position) const {
    const double query = coordinate(position, axis_);
    if (!std::isfinite(query)) {
        throw std::invalid_argument(
            "tabulated vector field query coordinate must be finite");
    }
    if (query < coordinates_.front() || query > coordinates_.back()) {
        throw std::out_of_range(
            "tabulated vector field query on axis " + to_string(axis_) +
            " lies outside [" + std::to_string(coordinates_.front()) +
            ", " + std::to_string(coordinates_.back()) + "]");
    }
    if (query == coordinates_.front()) return values_.front();
    if (query == coordinates_.back()) return values_.back();

    const auto upper =
        std::upper_bound(coordinates_.begin(), coordinates_.end(), query);
    const std::size_t right =
        static_cast<std::size_t>(upper - coordinates_.begin());
    const std::size_t left = right - 1;
    const double fraction =
        (query - coordinates_[left]) /
        (coordinates_[right] - coordinates_[left]);
    return interpolate(values_[left], values_[right], fraction);
}

void TabulatedVectorField1D::validate_domain(
    Vec3 minimum, Vec3 maximum,
    const std::string& context) const {
    const double domain_minimum = coordinate(minimum, axis_);
    const double domain_maximum = coordinate(maximum, axis_);
    if (!std::isfinite(domain_minimum) ||
        !std::isfinite(domain_maximum) ||
        domain_minimum > domain_maximum) {
        throw std::invalid_argument(
            context + " has an invalid domain on axis " + to_string(axis_));
    }
    if (coordinates_.front() > domain_minimum ||
        coordinates_.back() < domain_maximum) {
        throw std::invalid_argument(
            context + " tabulated field on axis " + to_string(axis_) +
            " does not cover the simulation domain [" +
            std::to_string(domain_minimum) + ", " +
            std::to_string(domain_maximum) + "]");
    }
}

TabulatedVectorField1D load_tabulated_vector_field_1d(
    const std::filesystem::path& path,
    CoordinateAxis axis) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error(
            "cannot open tabulated vector field: " + path.string());
    }

    std::vector<double> coordinates;
    std::vector<Vec3> values;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        const auto comment = line.find_first_of("#;");
        if (comment != std::string::npos) line.erase(comment);
        std::istringstream row(line);
        row >> std::ws;
        if (row.eof()) continue;

        double location = 0.0;
        Vec3 value{};
        if (!(row >> location >> value.x >> value.y >> value.z)) {
            throw std::runtime_error(
                "tabulated vector field " + path.string() +
                " line " + std::to_string(line_number) +
                " must contain coordinate field_x field_y field_z");
        }
        row >> std::ws;
        if (!row.eof()) {
            throw std::runtime_error(
                "tabulated vector field " + path.string() +
                " line " + std::to_string(line_number) +
                " has trailing columns");
        }
        coordinates.push_back(location);
        values.push_back(value);
    }
    try {
        return TabulatedVectorField1D(
            axis, std::move(coordinates), std::move(values));
    } catch (const std::invalid_argument& error) {
        throw std::runtime_error(
            "invalid tabulated vector field " + path.string() +
            ": " + error.what());
    }
}

} // namespace pic
