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
  src/Config.cpp
  src/Diagnostics.cpp
  src/FieldSolver.cpp
  src/Grid.cpp
  src/ImportedMesh2D.cpp
  src/Mesh2D.cpp
  src/Mesh3D.cpp
  src/Pusher.cpp
  src/Runtime.cpp
  src/Simulation.cpp
  src/Simulation2D.cpp
  src/Simulation3D.cpp
  src/Species.cpp
  src/Species2D.cpp
  src/Species3D.cpp
  src/UnstructuredFieldSolver2D.cpp
  src/UnstructuredMesh2D.cpp
  src/UnstructuredSimulation2D.cpp
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

install(TARGETS aurorapic aurorapic_cli
  EXPORT AuroraPICTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY examples/ DESTINATION ${CMAKE_INSTALL_DATADIR}/aurorapic/examples FILES_MATCHING PATTERN "*.cfg")
install(FILES README.md docs/methodology.md docs/multidimensional-roadmap.md docs/performance-envelope.md DESTINATION ${CMAKE_INSTALL_DOCDIR})
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
endif()
