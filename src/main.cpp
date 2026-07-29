#include "pic/Config.hpp"
#include "pic/Simulation.hpp"
#include "pic/Simulation2D.hpp"
#include "pic/Simulation3D.hpp"
#include "pic/UnstructuredSimulation2D.hpp"
#include <algorithm>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace {
constexpr std::string_view large_run_acknowledgement =
    "I_UNDERSTAND_THIS_IS_A_LARGE_RUN";
constexpr std::size_t local_particle_update_limit = 100'000'000;

template <typename Species>
bool exceeds_local_run_limit(
    const std::vector<Species>& species,
    std::size_t steps) {
    std::size_t particles = 0;
    for (const auto& item : species) {
        if (item.particles >
            std::numeric_limits<std::size_t>::max() - particles) {
            return true;
        }
        particles += item.particles;
    }
    return steps != 0 &&
        particles > local_particle_update_limit / steps;
}

void require_large_run_acknowledgement(
    bool exceeds_limit,
    bool acknowledged) {
    if (exceeds_limit && !acknowledged) {
        throw std::runtime_error(
            "estimated initial particle updates exceed the conservative "
            "100,000,000-update CLI limit; inspect with --validate-only, "
            "or deliberately launch with --allow-large-run "
            "I_UNDERSTAND_THIS_IS_A_LARGE_RUN");
    }
}
} // namespace

int main(int argc, char** argv) {
    bool validate_only = false;
    bool large_run_acknowledged = false;
    const char* config_path = nullptr;
    for (int index = 1; index < argc; ++index) {
        const std::string_view argument(argv[index]);
        if (argument == "--validate-only") {
            validate_only = true;
        } else if (argument == "--allow-large-run" && index + 1 < argc) {
            large_run_acknowledged =
                std::string_view(argv[++index]) ==
                large_run_acknowledgement;
        } else if (!config_path && !argument.starts_with("--")) {
            config_path = argv[index];
        } else {
            config_path = nullptr;
            break;
        }
    }
    if (!config_path) {
        std::cerr
            << "usage: aurorapic_cli [--validate-only] "
               "[--allow-large-run I_UNDERSTAND_THIS_IS_A_LARGE_RUN] "
               "<config.cfg>\n";
        return 2;
    }
    try {
        const unsigned dimension = pic::detect_config_dimension(config_path);
        if (dimension == 2) {
            if (pic::config_uses_unstructured_mesh_2d(config_path)) {
                auto cfg = pic::load_unstructured_config_2d(config_path);
                std::cout << "AuroraPIC imported 2D: mesh=" << cfg.mesh_path.string()
                          << " dt=" << cfg.dt
                          << " mode=" << pic::to_string(cfg.mode)
                          << " units=" << pic::to_string(cfg.units.system)
                          << " permittivity=" << cfg.units.permittivity()
                          << " vtk_output=" << (cfg.vtk_output ? "yes" : "no")
                          << "\n";
                if (validate_only) {
                    std::cout << "configuration valid; simulation not launched\n";
                    return 0;
                }
                require_large_run_acknowledgement(
                    exceeds_local_run_limit(
                        cfg.species,
                        cfg.mode == pic::RunMode::Transient
                            ? cfg.steps : cfg.max_steps),
                    large_run_acknowledged);
                pic::UnstructuredSimulation2D sim(std::move(cfg));
                const auto quality = sim.mesh().topology().quality();
                std::cout << "mesh nodes=" << sim.mesh().topology().nodes().size()
                          << " cells=" << sim.mesh().topology().cells().size()
                          << " boundary_faces="
                          << sim.mesh().topology().boundary_faces().size()
                          << " min_angle_degrees="
                          << quality.minimum_corner_angle_degrees
                          << " max_edge_ratio="
                          << quality.maximum_cell_edge_ratio << "\n";
                const auto summary = sim.run();
                std::cout << "completed steps=" << summary.steps_completed
                          << " time=" << summary.final_time
                          << " live_particles=" << summary.final_sample.live_particles
                          << " steady="
                          << (summary.steady_state_reached ? "yes" : "no")
                          << " poisson_residual="
                          << summary.final_sample.poisson.final_residual
                          << " total_energy=" << summary.final_sample.total_energy
                          << "\n";
                return summary.steady_state_reached ||
                               summary.steps_completed > 0
                           ? 0
                           : 1;
            }
            auto cfg = pic::load_config_2d(config_path);
            const pic::Boundary boundary_x =
                cfg.boundary_x.value_or(cfg.boundary);
            const pic::Boundary boundary_y =
                cfg.boundary_y.value_or(cfg.boundary);
            const char* field_solver =
                boundary_x == pic::Boundary::Periodic &&
                        boundary_y == pic::Boundary::Periodic
                    ? "spectral"
                    : boundary_x != boundary_y
                          ? "mixed_spectral_tridiagonal"
                          : "sor";
            std::cout << "AuroraPIC 2D: nx=" << cfg.nx << " ny=" << cfg.ny
                      << " length_x=" << cfg.length_x << " length_y=" << cfg.length_y
                      << " out_of_plane_depth=" << cfg.out_of_plane_depth
                      << " dt=" << cfg.dt << " mode=" << pic::to_string(cfg.mode)
                      << " units=" << pic::to_string(cfg.units.system)
                      << " permittivity=" << cfg.units.permittivity()
                      << " boundary_x="
                      << pic::to_string(boundary_x)
                      << " boundary_y="
                      << pic::to_string(boundary_y)
                      << " field_solver=" << field_solver
                      << " current_source="
                      << (cfg.current_regulated_source ? "yes" : "no")
                      << " potential_reference="
                      << (cfg.potential_reference ? "yes" : "no")
                      << " vtk_output=" << (cfg.vtk_output ? "yes" : "no") << "\n";
            if (validate_only) {
                std::cout << "configuration valid; simulation not launched\n";
                return 0;
            }
            require_large_run_acknowledgement(
                exceeds_local_run_limit(
                    cfg.species,
                    cfg.mode == pic::RunMode::Transient
                        ? cfg.steps : cfg.max_steps),
                large_run_acknowledged);
            pic::Simulation2D sim(std::move(cfg));
            auto summary = sim.run();
            std::cout << "completed steps=" << summary.steps_completed << " time=" << summary.final_time
                      << " live_particles=" << summary.final_sample.live_particles
                      << " steady=" << (summary.steady_state_reached ? "yes" : "no")
                      << " total_energy=" << summary.final_sample.total_energy << "\n";
            return summary.steps_completed > 0 ? 0 : 1;
        }
        if (dimension == 3) {
            auto cfg = pic::load_config_3d(config_path);
            std::cout << "AuroraPIC 3D: nx=" << cfg.nx << " ny=" << cfg.ny << " nz=" << cfg.nz
                      << " length_x=" << cfg.length_x << " length_y=" << cfg.length_y
                      << " length_z=" << cfg.length_z << " dt=" << cfg.dt
                      << " mode=" << pic::to_string(cfg.mode)
                      << " units=" << pic::to_string(cfg.units.system)
                      << " permittivity=" << cfg.units.permittivity()
                      << " boundary=" << pic::to_string(cfg.boundary)
                      << " vtk_output=" << (cfg.vtk_output ? "yes" : "no") << "\n";
            if (validate_only) {
                std::cout << "configuration valid; simulation not launched\n";
                return 0;
            }
            require_large_run_acknowledgement(
                exceeds_local_run_limit(
                    cfg.species,
                    cfg.mode == pic::RunMode::Transient
                        ? cfg.steps : cfg.max_steps),
                large_run_acknowledged);
            pic::Simulation3D sim(std::move(cfg));
            auto summary = sim.run();
            std::cout << "completed steps=" << summary.steps_completed << " time=" << summary.final_time
                      << " live_particles=" << summary.final_sample.live_particles
                      << " steady=" << (summary.steady_state_reached ? "yes" : "no")
                      << " total_energy=" << summary.final_sample.total_energy << "\n";
            return summary.steps_completed > 0 ? 0 : 1;
        }

        auto cfg = pic::load_config(config_path);
        std::cout << "AuroraPIC 1D"
                  << cfg.velocity_dimensions << "V: nx=" << cfg.nx
                  << " length=" << cfg.length << " dt=" << cfg.dt
                  << " mode=" << pic::to_string(cfg.mode)
                  << " units=" << pic::to_string(cfg.units.system)
                  << " permittivity=" << cfg.units.permittivity()
                  << " boundary=" << pic::to_string(cfg.boundary)
                  << " collisions="
                  << (!cfg.collision_models.empty()
                          ? "named_mcc(" +
                                std::to_string(
                                    std::count_if(
                                        cfg.collision_models.begin(),
                                        cfg.collision_models.end(),
                                        [](const auto& model) {
                                            return model.config.enabled;
                                        })) +
                                ")"
                          : (cfg.collisions.enabled
                                 ? pic::to_string(
                                       cfg.collisions.model)
                                 : "off"))
                  << "\n";
        if (validate_only) {
            std::cout << "configuration valid; simulation not launched\n";
            return 0;
        }
        require_large_run_acknowledgement(
            exceeds_local_run_limit(
                cfg.species,
                cfg.mode == pic::RunMode::Transient
                    ? cfg.steps : cfg.max_steps),
            large_run_acknowledged);
        pic::Simulation sim(std::move(cfg));
        auto summary = sim.run();
        std::cout << "completed steps=" << summary.steps_completed << " time=" << summary.final_time
                  << " steady=" << (summary.steady_state_reached ? "yes" : "no")
                  << " total_energy=" << summary.final_sample.total_energy << "\n";
        return summary.steady_state_reached || summary.steps_completed > 0 ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
