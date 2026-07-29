#pragma once
#include <cstddef>
#include <cmath>
#include <stdexcept>
#include <string>
namespace pic {
constexpr double NORMALIZED_PERMITTIVITY = 1.0;
constexpr double VACUUM_PERMITTIVITY_SI = 8.8541878128e-12;
constexpr double EPS0 = NORMALIZED_PERMITTIVITY; // compatibility alias

enum class UnitSystem { Normalized, SI };

struct UnitSystemConfig {
    UnitSystem system{UnitSystem::Normalized};
    double relative_permittivity{1.0};

    double permittivity() const {
        if (!std::isfinite(relative_permittivity) ||
            !(relative_permittivity > 0.0)) {
            throw std::invalid_argument(
                "relative_permittivity must be positive and finite");
        }
        const double base = system == UnitSystem::SI
                                ? VACUUM_PERMITTIVITY_SI
                                : NORMALIZED_PERMITTIVITY;
        return base * relative_permittivity;
    }
};

enum class Boundary { Periodic, Dirichlet };
enum class RunMode { Transient, SteadyState };
enum class ParticleBoundary { Auto, Absorbing, Reflecting, Periodic };
enum class VTKOutputFormat { Legacy, Xml, Both };

inline std::string to_string(UnitSystem system) {
    switch (system) {
        case UnitSystem::Normalized: return "normalized";
        case UnitSystem::SI: return "si";
    }
    return "unknown";
}

inline std::string to_string(Boundary boundary) {
    switch (boundary) {
        case Boundary::Periodic: return "periodic";
        case Boundary::Dirichlet: return "dirichlet";
    }
    return "unknown";
}

inline std::string to_string(RunMode mode) {
    switch (mode) {
        case RunMode::Transient: return "transient";
        case RunMode::SteadyState: return "steady_state";
    }
    return "unknown";
}

inline std::string to_string(ParticleBoundary boundary) {
    switch (boundary) {
        case ParticleBoundary::Auto: return "auto";
        case ParticleBoundary::Absorbing: return "absorbing";
        case ParticleBoundary::Reflecting: return "reflecting";
        case ParticleBoundary::Periodic: return "periodic";
    }
    return "unknown";
}

inline std::string to_string(VTKOutputFormat format) {
    switch (format) {
        case VTKOutputFormat::Legacy: return "legacy";
        case VTKOutputFormat::Xml: return "xml";
        case VTKOutputFormat::Both: return "both";
    }
    return "unknown";
}

struct Vec2 {
    double x{0.0};
    double y{0.0};
};

struct Vec3 {
    double x{0.0};
    double y{0.0};
    double z{0.0};
};

struct Particle {
    double x{0.0};
    double v{0.0};          // time-centered velocity used for diagnostics/output
    bool alive{true};
    double v_half{0.0};     // leapfrog velocity at the adjacent half time step
    double velocity_y{0.0};
    double velocity_z{0.0};
};

struct Particle2D {
    Vec2 position{};
    Vec2 velocity{};        // time-centered velocity used for diagnostics/output
    bool alive{true};
    Vec2 velocity_half{};   // leapfrog velocity at the adjacent half time step
    double velocity_z{0.0};
    double velocity_half_z{0.0};
};

struct Particle3D {
    Vec3 position{};
    Vec3 velocity{};        // time-centered velocity used for diagnostics/output
    bool alive{true};
    Vec3 velocity_half{};   // leapfrog velocity at the adjacent half time step
};
}
