#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace pic {

template <typename Sample>
bool adjacent_energy_windows_converged(const std::vector<Sample>& history,
                                       std::size_t window,
                                       double relative_tolerance) {
    if (window == 0 || history.size() / 2 < window ||
        !std::isfinite(relative_tolerance) || relative_tolerance <= 0.0) {
        return false;
    }

    long double current_window_sum = 0.0L;
    long double previous_window_sum = 0.0L;
    for (std::size_t i = history.size() - window; i < history.size(); ++i) {
        if (!std::isfinite(history[i].total_energy)) return false;
        current_window_sum += static_cast<long double>(history[i].total_energy);
    }
    for (std::size_t i = history.size() - 2 * window; i < history.size() - window; ++i) {
        if (!std::isfinite(history[i].total_energy)) return false;
        previous_window_sum += static_cast<long double>(history[i].total_energy);
    }

    const long double scale = static_cast<long double>(window);
    const long double current_window_mean = current_window_sum / scale;
    const long double previous_window_mean = previous_window_sum / scale;
    const long double relative_change =
        std::abs(current_window_mean - previous_window_mean) /
        std::max(1.0e-30L, std::abs(previous_window_mean));
    return relative_change < static_cast<long double>(relative_tolerance);
}

// Statistical readiness criteria for block-averaged observables from a
// periodically driven simulation. Blocks must contain complete drive cycles;
// constructing those blocks is deliberately separate from evaluating them.
struct BlockConvergenceCriteria {
    std::size_t minimum_blocks{16};
    double minimum_effective_blocks{8.0};
    double maximum_absolute_projected_fractional_drift{0.01};
    double maximum_absolute_split_half_fractional_change{0.01};
    double maximum_relative_standard_error{
        std::numeric_limits<double>::infinity()};
};

enum class BlockConvergenceClassification {
    Invalid,
    HorizonIncomplete,
    EffectiveSampleIncomplete,
    NonStationary,
    Converged
};

struct BlockConvergenceResult {
    BlockConvergenceClassification classification{
        BlockConvergenceClassification::Invalid};
    std::size_t blocks{0};
    double mean{std::numeric_limits<double>::quiet_NaN()};
    double sample_standard_deviation{
        std::numeric_limits<double>::quiet_NaN()};
    std::optional<double> lag_one_correlation{};
    std::optional<double> effective_blocks{};
    std::optional<double> projected_fractional_drift{};
    std::optional<double> split_half_fractional_change{};
    std::optional<double> relative_standard_error{};
    bool minimum_blocks_passed{false};
    bool minimum_effective_blocks_passed{false};
    bool drift_passed{false};
    bool split_half_passed{false};
    bool relative_standard_error_passed{false};

    bool converged() const {
        return classification == BlockConvergenceClassification::Converged;
    }
};

BlockConvergenceResult evaluate_block_convergence(
    const std::vector<double>& block_means,
    const BlockConvergenceCriteria& criteria = {});

std::string to_string(BlockConvergenceClassification classification);

struct PeriodicBlockState {
    std::size_t steps_per_cycle{0};
    std::size_t cycles_per_block{0};
    std::size_t origin_step{0};
    std::size_t last_step{0};
    std::size_t samples_in_current_block{0};
    std::vector<std::string> observable_names{};
    std::vector<long double> current_sums{};
    std::vector<std::vector<double>> completed_block_means{};
};

// Collects one vector of observables after every simulation step. The origin
// must be an RF phase boundary; strict step contiguity prevents incomplete or
// phase-shifted cycles from entering the statistical decision.
class PeriodicBlockConvergence {
public:
    PeriodicBlockConvergence(
        std::size_t steps_per_cycle,
        std::size_t cycles_per_block,
        std::size_t origin_step,
        std::vector<std::string> observable_names,
        BlockConvergenceCriteria criteria = {});

    PeriodicBlockConvergence(
        PeriodicBlockState state,
        BlockConvergenceCriteria criteria = {});

    bool observe(std::size_t step, const std::vector<double>& values);
    std::vector<BlockConvergenceResult> evaluate() const;
    bool converged() const;
    std::size_t completed_blocks() const;
    std::size_t samples_per_block() const;
    const PeriodicBlockState& state() const { return state_; }

private:
    void validate_state() const;
    PeriodicBlockState state_{};
    BlockConvergenceCriteria criteria_{};
};

} // namespace pic
