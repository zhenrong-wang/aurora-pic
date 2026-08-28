#include "pic/Convergence.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void require_near(double actual, double expected, double tolerance,
                  const std::string& message) {
    if (std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(
            message + ": actual=" + std::to_string(actual) +
            " expected=" + std::to_string(expected));
    }
}

void test_constant_series_converges() {
    const std::vector<double> blocks(16, 4.0);
    const auto result = pic::evaluate_block_convergence(blocks);
    require(result.converged(), "constant series should converge");
    require(result.effective_blocks.has_value(),
            "constant series should have an effective count");
    require_near(*result.effective_blocks, 16.0, 1e-12,
                 "constant effective count");
    require_near(*result.projected_fractional_drift, 0.0, 1e-12,
                 "constant drift");
    require_near(*result.relative_standard_error, 0.0, 1e-12,
                 "constant uncertainty");
}

void test_short_series_is_horizon_incomplete() {
    const std::vector<double> blocks(8, 2.0);
    const auto result = pic::evaluate_block_convergence(blocks);
    require(result.classification ==
                pic::BlockConvergenceClassification::HorizonIncomplete,
            "short series should not claim convergence");
}

void test_slow_correlated_mode_is_rejected() {
    // Deterministic AR-like relaxation with small alternating perturbations.
    // Ordinary adjacent changes are small, but serial correlation reduces the
    // independent information below the eight-block gate.
    std::vector<double> blocks;
    double value = 1.0;
    for (std::size_t i = 0; i < 16; ++i) {
        value = 0.92 * value + 0.08 * 2.0 +
            (i % 2 == 0 ? 0.001 : -0.001);
        blocks.push_back(value);
    }
    pic::BlockConvergenceCriteria criteria;
    criteria.maximum_absolute_projected_fractional_drift = 1.0;
    criteria.maximum_absolute_split_half_fractional_change = 1.0;
    const auto result = pic::evaluate_block_convergence(blocks, criteria);
    require(result.classification ==
                pic::BlockConvergenceClassification::EffectiveSampleIncomplete,
            "correlated mode should fail effective-sample gate");
    require(result.lag_one_correlation &&
                *result.lag_one_correlation > 0.5,
            "correlated fixture should retain positive lag-one correlation");
    require(result.effective_blocks && *result.effective_blocks < 8.0,
            "correlated fixture should have fewer than eight effective blocks");
}

void test_drift_is_rejected_after_effective_sample_gate() {
    std::vector<double> blocks;
    for (std::size_t i = 0; i < 16; ++i) {
        // Alternation suppresses lag-one correlation while the linear term
        // creates an unambiguous secular drift.
        blocks.push_back(10.0 + 0.03 * static_cast<double>(i) +
                         (i % 2 == 0 ? 0.3 : -0.3));
    }
    pic::BlockConvergenceCriteria criteria;
    criteria.minimum_effective_blocks = 4.0;
    criteria.maximum_absolute_projected_fractional_drift = 0.01;
    criteria.maximum_absolute_split_half_fractional_change = 1.0;
    const auto result = pic::evaluate_block_convergence(blocks, criteria);
    require(result.minimum_effective_blocks_passed,
            "drift fixture should clear effective-sample gate");
    require(!result.drift_passed,
            "drift fixture should fail projected-drift gate");
    require(result.classification ==
                pic::BlockConvergenceClassification::NonStationary,
            "drifting series should be nonstationary");
}

void test_nonfinite_and_zero_mean_are_invalid() {
    auto values = std::vector<double>(16, 1.0);
    values[4] = std::numeric_limits<double>::quiet_NaN();
    require(pic::evaluate_block_convergence(values).classification ==
                pic::BlockConvergenceClassification::Invalid,
            "nonfinite blocks must be invalid");
    require(pic::evaluate_block_convergence(
                std::vector<double>(16, 0.0)).classification ==
                pic::BlockConvergenceClassification::Invalid,
            "zero-mean relative convergence must be invalid");
}

void test_periodic_blocks_are_complete_and_restart_safe() {
    pic::BlockConvergenceCriteria criteria;
    criteria.minimum_blocks = 4;
    criteria.minimum_effective_blocks = 4.0;
    pic::PeriodicBlockConvergence continuous(
        4, 2, 100, {"electron_population", "total_energy"}, criteria);
    pic::PeriodicBlockConvergence split(
        4, 2, 100, {"electron_population", "total_energy"}, criteria);

    for (std::size_t step = 101; step <= 117; ++step) {
        const std::vector<double> values{5.0, 9.0};
        const bool continuous_closed = continuous.observe(step, values);
        const bool split_closed = split.observe(step, values);
        require(continuous_closed == (step == 108 || step == 116),
                "block must close only after two complete RF cycles");
        require(split_closed == continuous_closed,
                "split prefix must match continuous closure");
    }
    require(continuous.completed_blocks() == 2,
            "continuous prefix should contain two blocks");
    require(continuous.state().samples_in_current_block == 1,
            "continuous prefix should retain partial block");

    pic::PeriodicBlockConvergence restarted(split.state(), criteria);
    for (std::size_t step = 118; step <= 132; ++step) {
        const std::vector<double> values{5.0, 9.0};
        continuous.observe(step, values);
        restarted.observe(step, values);
    }
    require(continuous.completed_blocks() == 4,
            "four complete blocks expected");
    require(restarted.state().last_step == continuous.state().last_step,
            "restart should preserve last step");
    require(restarted.state().samples_in_current_block ==
                continuous.state().samples_in_current_block,
            "restart should preserve partial-block sample count");
    require(restarted.state().current_sums == continuous.state().current_sums,
            "restart should preserve partial-block sums");
    require(restarted.state().completed_block_means ==
                continuous.state().completed_block_means,
            "restart should reproduce block means exactly");
    require(restarted.converged(),
            "constant complete-cycle observables should converge");
}

void test_periodic_blocks_reject_gaps_and_bad_state() {
    pic::PeriodicBlockConvergence blocks(4, 1, 20, {"density"});
    bool gap_rejected = false;
    try {
        blocks.observe(22, {1.0});
    } catch (const std::invalid_argument&) {
        gap_rejected = true;
    }
    require(gap_rejected, "step gaps must be rejected");

    auto state = blocks.state();
    state.samples_in_current_block = 4;
    bool state_rejected = false;
    try {
        pic::PeriodicBlockConvergence invalid(state);
    } catch (const std::invalid_argument&) {
        state_rejected = true;
    }
    require(state_rejected, "invalid restored state must be rejected");
}

void test_turner_slow_mode_golden_series() {
    // Checksum-pinned seed 162631810 continuation from the corrected Turner
    // Case 1 campaign (16 consecutive 32-cycle density blocks).
    const std::vector<double> blocks{
        4873376247863.77, 4910235080871.582, 4864251791687.014,
        4866666924133.303, 4942945692138.675, 4955437344055.173,
        4984151939392.09, 4948202026977.542, 4966138553466.799,
        4899038262023.926, 4875369195251.466, 4874410220947.267,
        4881444525756.835, 4901569961547.853, 4888318977661.131,
        4916364329528.81};
    const auto result = pic::evaluate_block_convergence(blocks);
    require_near(*result.lag_one_correlation, 0.562893888984654,
                 1e-14, "Turner lag-one correlation parity");
    require_near(*result.effective_blocks, 4.474838519452556,
                 1e-13, "Turner effective-block parity");
    require_near(*result.projected_fractional_drift,
                 -0.00096194577209566, 1e-15,
                 "Turner projected-drift parity");
    require_near(*result.split_half_fractional_change,
                 -0.0036312360400212714, 1e-15,
                 "Turner split-half parity");
    require(result.classification ==
                pic::BlockConvergenceClassification::EffectiveSampleIncomplete,
            "Turner golden series must reproduce slow-mode classification");
}

} // namespace

int main() {
    try {
        test_constant_series_converges();
        test_short_series_is_horizon_incomplete();
        test_slow_correlated_mode_is_rejected();
        test_drift_is_rejected_after_effective_sample_gate();
        test_nonfinite_and_zero_mean_are_invalid();
        test_periodic_blocks_are_complete_and_restart_safe();
        test_periodic_blocks_reject_gaps_and_bad_state();
        test_turner_slow_mode_golden_series();
        std::cout << "convergence tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "convergence test failed: " << error.what() << '\n';
        return 1;
    }
}
