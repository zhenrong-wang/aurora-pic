cmake_minimum_required(VERSION 3.20)
project(AuroraPIC VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(AURORA_BUILD_TESTS "Build AuroraPIC tests" ON)
add_library(aurorapic
  src/Config.cpp
  src/Diagnostics.cpp
  src/FieldSolver.cpp
  src/Grid.cpp
  src/Mesh2D.cpp
  src/Mesh3D.cpp
  src/Pusher.cpp
  src/Simulation.cpp
  src/Simulation2D.cpp
  src/Species.cpp
  src/Species2D.cpp
  src/Species3D.cpp
  src/VTKWriter.cpp
)
target_include_directories(aurorapic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(aurorapic PRIVATE -Wall -Wextra -Wpedantic -O3)

add_executable(aurorapic_cli src/main.cpp)
target_link_libraries(aurorapic_cli PRIVATE aurorapic)
target_compile_options(aurorapic_cli PRIVATE -Wall -Wextra -Wpedantic -O3)

if(AURORA_BUILD_TESTS)
  enable_testing()
  add_executable(aurorapic_tests tests/test_core.cpp)
  target_link_libraries(aurorapic_tests PRIVATE aurorapic)
  target_compile_options(aurorapic_tests PRIVATE -Wall -Wextra -Wpedantic -O2)
  add_test(NAME aurorapic_core COMMAND aurorapic_tests)
endif()
