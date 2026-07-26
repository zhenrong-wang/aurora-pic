#pragma once
#include <cstddef>
#include <string>
namespace pic {
constexpr double EPS0 = 1.0; // normalized permittivity; SI scaling can be layered via config later.

enum class Boundary { Periodic, Dirichlet };
enum class RunMode { Transient, SteadyState };
enum class ParticleBoundary { Auto, Absorbing, Reflecting, Periodic };

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

struct Vec2 {
    double x{0.0};
    double y{0.0};
};

struct Particle {
    double x{0.0};
    double v{0.0};
    bool alive{true};
};

struct Particle2D {
    Vec2 position{};
    Vec2 velocity{};
    bool alive{true};
};
}
