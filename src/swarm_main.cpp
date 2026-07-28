#include "pic/GasDataset.hpp"
#include "pic/Swarm.hpp"

#include <exception>
#include <iostream>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: aurorapic_swarm <swarm.cfg>\n";
        return 2;
    }
    try {
        const auto config =
            pic::load_swarm_benchmark_config(argv[1]);
        const auto dataset =
            pic::load_gas_dataset(config.gas_data_file);
        const auto results = pic::run_swarm_benchmark(config);
        pic::write_swarm_benchmark_csv(
            config.output_file, config, dataset, results);
        for (const auto& result : results) {
            std::cout
                << "E/N=" << result.reduced_field_td
                << " Td electron_drift_velocity="
                << result.electron_drift_velocity_m_s
                << " m/s mean_energy=" << result.mean_energy_ev
                << " eV max_energy="
                << result.maximum_observed_energy_ev << " eV";
            if (config.population_model ==
                pic::SwarmPopulationModel::BranchingResampled) {
                std::cout
                    << " growth_rate="
                    << result.temporal_growth_rate_s << " 1/s"
                    << " electron_weight="
                    << result.final_total_electron_weight;
                if (result.townsend_available) {
                    std::cout
                        << " growth_over_flux_drift_townsend="
                        << result
                               .growth_over_flux_drift_townsend_1_m
                        << " 1/m";
                }
            }
            std::cout << '\n';
        }
        std::cout << "wrote " << config.output_file.string() << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "swarm error: " << error.what() << '\n';
        return 1;
    }
}
