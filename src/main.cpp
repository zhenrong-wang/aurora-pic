#include "pic/Config.hpp"
#include "pic/Simulation.hpp"
#include "pic/Simulation2D.hpp"
#include "pic/Simulation3D.hpp"
#include "pic/UnstructuredSimulation2D.hpp"
#include <algorithm>
#include <exception>
#include <iostream>
#include <string_view>
#include <utility>

int main(int argc, char** argv) {
    const bool validate_only =
        argc == 3 && std::string_view(argv[1]) == "--validate-only";
    if ((!validate_only && argc != 2) || (validate_only && argc != 3)) {
        std::cerr
            << "usage: aurorapic_cli [--validate-only] <config.cfg>\n";
        return 2;
    }
    const char* config_path = argv[validate_only ? 2 : 1];
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
