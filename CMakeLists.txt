cmake_minimum_required(VERSION 3.20)
project(AuroraPIC VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
option(AURORA_BUILD_TESTS "Build AuroraPIC tests" ON)
option(AURORA_ENABLE_OPENMP "Enable OpenMP runtime backend when available" ON)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
set(AURORAPIC_USES_OPENMP OFF)

add_library(aurorapic
  src/Collision.cpp
  src/Config.cpp
  src/Diagnostics.cpp
  src/ResolvedDiagnostics2D.cpp
  src/FieldSolver.cpp
  src/GasDataset.cpp
  src/Grid.cpp
  src/ImportedMesh2D.cpp
  src/Initialization.cpp
  src/Mesh2D.cpp
  src/Mesh3D.cpp
  src/ParticleState.cpp
  src/PrescribedField.cpp
  src/Pusher.cpp
  src/Runtime.cpp
  src/Simulation.cpp
  src/Simulation2D.cpp
  src/Simulation3D.cpp
  src/Species.cpp
  src/Species2D.cpp
  src/Species3D.cpp
  src/Swarm.cpp
  src/UnstructuredConfig2D.cpp
  src/UnstructuredFieldSolver2D.cpp
  src/UnstructuredMesh2D.cpp
  src/UnstructuredSimulation2D.cpp
  src/Units.cpp
  src/VTKWriter.cpp
)
target_include_directories(aurorapic PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_features(aurorapic PUBLIC cxx_std_20)
target_compile_options(aurorapic PRIVATE -Wall -Wextra -Wpedantic -O3)

if(AURORA_ENABLE_OPENMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(aurorapic PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(aurorapic PUBLIC AURORA_HAVE_OPENMP=1)
    set(AURORAPIC_USES_OPENMP ON)
  endif()
endif()

add_executable(aurorapic_cli src/main.cpp)
target_link_libraries(aurorapic_cli PRIVATE aurorapic)
target_compile_options(aurorapic_cli PRIVATE -Wall -Wextra -Wpedantic -O3)

add_executable(aurorapic_swarm src/swarm_main.cpp)
target_link_libraries(aurorapic_swarm PRIVATE aurorapic)
target_compile_options(aurorapic_swarm PRIVATE -Wall -Wextra -Wpedantic -O3)

install(TARGETS aurorapic aurorapic_cli aurorapic_swarm
  EXPORT AuroraPICTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY examples/ DESTINATION ${CMAKE_INSTALL_DATADIR}/aurorapic/examples
  FILES_MATCHING PATTERN "*.cfg" PATTERN "*.msh" PATTERN "*.geo"
  PATTERN "*.dat" PATTERN "*.gas" PATTERN "*.swarm"
  PATTERN "*.aps" PATTERN "*.case" PATTERN "*.sources")
install(PROGRAMS scripts/compare_swarm.py scripts/compare_hall.py
        scripts/compare_turner.py
        scripts/analyze_turner_balance.py
        scripts/analyze_turner_density_blocks.py
        scripts/digitize_charoy_figure6.py
  scripts/aggregate_hall_ensemble.py scripts/analyze_hall_pilot.py
  scripts/preflight_hall.py scripts/qualify_hall_runtime.py
  scripts/prepare_hall_campaign.py scripts/lock_hall_source.py
  scripts/verify_turner_source.py scripts/normalize_turner_source.py
  scripts/prepare_turner_case.py
  scripts/prepare_turner_ensemble.py
  scripts/attach_turner_ensemble_result.py
  scripts/analyze_turner_ensemble.py
  scripts/audit_turner_credibility.py
  scripts/compare_edupic_phase_space.py
  scripts/analyze_edupic_heating_trend.py
  scripts/analyze_edupic_post_trend_horizon.py
  scripts/prepare_turner_sensitivity.py
  scripts/qualify_turner_runtime.py
  scripts/run_turner_startup.py
  scripts/extend_turner_horizon.py
  scripts/analyze_turner_stationarity.py
  scripts/analyze_turner_spatial_structure.py
  scripts/analyze_1d_energy_budget.py
  scripts/analyze_1d_spatial_collision.py
  scripts/run_aurorapic_edupic_pilot.py
  scripts/extend_aurorapic_edupic_horizon.py
  scripts/export_checkpoint_particle_state.py
  scripts/run_aurorapic_initialization_ab.py
  scripts/run_aurorapic_initialization_followup.py
  scripts/run_aurorapic_density_bracket.py
  scripts/prepare_quasineutral_particle_state.py
  scripts/augment_neutral_particle_pairs.py
  scripts/analyze_phase_eedf.py
  scripts/analyze_edupic_eedf_discrepancy.py
  scripts/analyze_pic_resolution.py
  scripts/analyze_aurorapic_timestep_refinement.py
  scripts/analyze_aurorapic_mesh_refinement.py
  scripts/analyze_aurorapic_particle_refinement.py
  scripts/analyze_aurorapic_region_matched_collisions.py
  scripts/analyze_aurorapic_heating_ledger.py
  scripts/analyze_aurorapic_matched_heating.py
  scripts/analyze_aurorapic_heating_localization.py
  scripts/analyze_aurorapic_heating_seed_ensemble.py
  scripts/analyze_aurorapic_heating_factor_attribution.py
  scripts/analyze_aurorapic_long_window_heating.py
  scripts/audit_edupic_normalization.py
  scripts/analyze_aurorapic_microstate_ensemble.py
  scripts/analyze_edupic_transport_balance.py
  scripts/run_aurorapic_balance_stationarity_block.py
  scripts/compare_aurorapic_measurement_windows.py
  scripts/run_aurorapic_timestep_refinement.py
  scripts/split_particle_state.py
  scripts/randomize_particle_state.py
  scripts/export_phase_eedf.py
  scripts/compare_phase_eedf.py
  scripts/run_edupic_stage.py
  scripts/run_edupic_measurement_stage.py
  scripts/advance_edupic_measurement.py
  scripts/analyze_edupic_measurement_blocks.py
  scripts/evaluate_edupic_measurement_stationarity.py
  scripts/analyze_edupic_convergence.py
  scripts/advance_edupic_equilibration.py
  scripts/normalize_hall_reference.py scripts/prepare_hall_ensemble.py
  scripts/prepare_hall_convergence.py scripts/analyze_hall_convergence.py
  scripts/analyze_hall_flux_stationarity.py
  scripts/compare_hall_flux_blocks.py
  scripts/prepare_hall_horizon_stage.py
  scripts/analyze_hall_horizon_stage.py
  scripts/run_swarm_campaign.py
  scripts/generate_real_case_mesh.py scripts/import_lxcat.py
  scripts/validate_kinetic_benchmarks.py scripts/validate_hall_case.py
  DESTINATION ${CMAKE_INSTALL_DATADIR}/aurorapic/scripts)
install(FILES README.md docs/ccp-validation.md docs/collisions.md docs/external-particle-state.md docs/gas-data-workflow.md docs/hall-thruster-validation.md docs/kinetic-validation.md docs/methodology.md docs/multidimensional-roadmap.md docs/performance-envelope.md docs/phase-eedf-interchange.md docs/real-case-validation.md docs/swarm-validation.md docs/units.md DESTINATION ${CMAKE_INSTALL_DOCDIR})
install(EXPORT AuroraPICTargets
  FILE AuroraPICTargets.cmake
  NAMESPACE AuroraPIC::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AuroraPIC
)
configure_package_config_file(
  cmake/AuroraPICConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/AuroraPICConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AuroraPIC
)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/AuroraPICConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/AuroraPICConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/AuroraPICConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AuroraPIC
)

set(CPACK_PACKAGE_NAME "AuroraPIC")
set(CPACK_PACKAGE_VENDOR "AuroraPIC")
set(CPACK_PACKAGE_CONTACT "AuroraPIC maintainers")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Electrostatic PIC research prototype with 1D/2D/3D examples")
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set(CPACK_GENERATOR "TGZ")
include(CPack)

if(AURORA_BUILD_TESTS)
  enable_testing()
  add_executable(aurorapic_tests tests/test_core.cpp)
  target_link_libraries(aurorapic_tests PRIVATE aurorapic)
  target_compile_options(aurorapic_tests PRIVATE -Wall -Wextra -Wpedantic -O2)
  target_compile_definitions(aurorapic_tests PRIVATE AURORA_TEST_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_test(NAME aurorapic_core COMMAND aurorapic_tests)

  add_executable(aurorapic_initialization_tests tests/test_initialization.cpp)
  target_link_libraries(aurorapic_initialization_tests PRIVATE aurorapic)
  target_compile_options(aurorapic_initialization_tests PRIVATE -Wall -Wextra -Wpedantic -O2)
  target_compile_definitions(aurorapic_initialization_tests PRIVATE AURORA_TEST_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_test(NAME aurorapic_initialization COMMAND aurorapic_initialization_tests)
endif()
