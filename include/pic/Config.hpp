#pragma once
#include "pic/Types.hpp"
#include <cstddef>
#include <string>
#include <vector>

namespace pic {
struct Simulation2DConfig;
struct Simulation3DConfig;

struct SpeciesConfig {
    std::string name{"electrons"};
    double charge{-1.0};
    double mass{1.0};
    double weight{1.0};
    std::size_t particles{1000};
    double density{1.0};
    double drift_velocity{0.0};
    double thermal_velocity{0.1};
    double init_x_min{0.0};
    double init_x_max{-1.0}; // negative means full domain
};

struct CollisionConfig {
    bool enabled{false};
    double frequency{0.0};
    double neutral_temperature_velocity{0.0};
};

struct Config {
    std::size_t nx{128};
    double length{1.0};
    double dt{0.02};
    std::size_t steps{100};
    std::size_t output_interval{10};
    Boundary boundary{Boundary::Periodic};
    RunMode mode{RunMode::Transient};
    double phi_left{0.0};
    double phi_right{0.0};
    unsigned seed{12345};
    std::string output_dir{"output"};
    double steady_tolerance{1e-6};
    std::size_t steady_window{25};
    std::size_t max_steps{10000};
    CollisionConfig collisions{};
    bool checkpoint_output{false};
    std::size_t checkpoint_interval{0}; // zero inherits output_interval
    std::string checkpoint_path{};      // empty writes output_dir/checkpoint_<step>.apc
    std::string restart_path{};
    std::vector<SpeciesConfig> species{};
};

Config load_config(const std::string& path);
Simulation2DConfig load_config_2d(const std::string& path);
Simulation3DConfig load_config_3d(const std::string& path);
unsigned detect_config_dimension(const std::string& path);
std::string to_string(Boundary b);
std::string to_string(RunMode m);
}
