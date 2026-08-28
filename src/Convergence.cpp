#include "pic/Convergence.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace pic {
namespace {

bool valid_criteria(const BlockConvergenceCriteria& value) {
    return value.minimum_blocks >= 2 &&
        std::isfinite(value.minimum_effective_blocks) &&
        value.minimum_effective_blocks > 0.0 &&
        std::isfinite(value.maximum_absolute_projected_fractional_drift) &&
        value.maximum_absolute_projected_fractional_drift >= 0.0 &&
        std::isfinite(
            value.maximum_absolute_split_half_fractional_change) &&
        value.maximum_absolute_split_half_fractional_change >= 0.0 &&
        !std::isnan(value.maximum_relative_standard_error) &&
        value.maximum_relative_standard_error >= 0.0;
}

double arithmetic_mean(const std::vector<double>& values) {
    const long double sum = std::accumulate(
        values.begin(), values.end(), 0.0L,
        [](long double total, double value) {
            return total + static_cast<long double>(value);
        });
    return static_cast<double>(sum / values.size());
}

double linear_slope(const std::vector<double>& values) {
    const long double x_mean =
        static_cast<long double>(values.size() - 1) / 2.0L;
    const long double y_mean = arithmetic_mean(values);
    long double numerator = 0.0L;
    long double denominator = 0.0L;
    for (std::size_t i = 0; i < values.size(); ++i) {
        const long double dx = static_cast<long double>(i) - x_mean;
        numerator += dx * (static_cast<long double>(values[i]) - y_mean);
        denominator += dx * dx;
    }
    return static_cast<double>(numerator / denominator);
}

std::optional<double> lag_one_correlation(
    const std::vector<double>& values, double mean) {
    if (values.size() < 3) return std::nullopt;
    long double numerator = 0.0L;
    long double variance = 0.0L;
    for (double value : values) {
        const long double centered =
            static_cast<long double>(value) - mean;
        variance += centered * centered;
    }
    for (std::size_t i = 0; i + 1 < values.size(); ++i) {
        const long double left =
            static_cast<long double>(values[i]) - mean;
        const long double right =
            static_cast<long double>(values[i + 1]) - mean;
        numerator += left * right;
    }
    if (variance == 0.0L) {
        // An exactly constant series has no serial uncertainty. Treating its
        // correlation as zero gives it the full nominal effective count.
        return 0.0;
    }
    return static_cast<double>(numerator / variance);
}

} // namespace

BlockConvergenceResult evaluate_block_convergence(
    const std::vector<double>& block_means,
    const BlockConvergenceCriteria& criteria) {
    BlockConvergenceResult result;
    result.blocks = block_means.size();
    if (!valid_criteria(criteria) || block_means.size() < 2 ||
        !std::all_of(block_means.begin(), block_means.end(),
                     [](double value) { return std::isfinite(value); })) {
        return result;
    }

    result.mean = arithmetic_mean(block_means);
    if (!std::isfinite(result.mean) || result.mean == 0.0) return result;

    long double squared_deviation = 0.0L;
    for (double value : block_means) {
        const long double delta =
            static_cast<long double>(value) - result.mean;
        squared_deviation += delta * delta;
    }
    result.sample_standard_deviation = static_cast<double>(std::sqrt(
        squared_deviation /
        static_cast<long double>(block_means.size() - 1)));
    result.lag_one_correlation =
        lag_one_correlation(block_means, result.mean);
    if (result.lag_one_correlation) {
        const double bounded = std::clamp(
            *result.lag_one_correlation, -0.99, 0.99);
        result.effective_blocks = std::clamp(
            static_cast<double>(block_means.size()) *
                (1.0 - bounded) / (1.0 + bounded),
            1.0, static_cast<double>(block_means.size()));
    }

    result.projected_fractional_drift =
        linear_slope(block_means) *
        static_cast<double>(block_means.size() - 1) / result.mean;
    const std::size_t half = block_means.size() / 2;
    if (half > 0) {
        const std::vector<double> first(
            block_means.begin(), block_means.begin() + half);
        const std::vector<double> second(
            block_means.end() - half, block_means.end());
        result.split_half_fractional_change =
            (arithmetic_mean(second) - arithmetic_mean(first)) / result.mean;
    }
    if (result.effective_blocks) {
        result.relative_standard_error =
            result.sample_standard_deviation /
            std::sqrt(*result.effective_blocks) / std::abs(result.mean);
    }

    result.minimum_blocks_passed =
        block_means.size() >= criteria.minimum_blocks;
    result.minimum_effective_blocks_passed =
        result.effective_blocks &&
        *result.effective_blocks >= criteria.minimum_effective_blocks;
    result.drift_passed = result.projected_fractional_drift &&
        std::abs(*result.projected_fractional_drift) <=
            criteria.maximum_absolute_projected_fractional_drift;
    result.split_half_passed = result.split_half_fractional_change &&
        std::abs(*result.split_half_fractional_change) <=
            criteria.maximum_absolute_split_half_fractional_change;
    result.relative_standard_error_passed =
        result.relative_standard_error &&
        *result.relative_standard_error <=
            criteria.maximum_relative_standard_error;

    if (!result.minimum_blocks_passed) {
        result.classification =
            BlockConvergenceClassification::HorizonIncomplete;
    } else if (!result.minimum_effective_blocks_passed) {
        result.classification =
            BlockConvergenceClassification::EffectiveSampleIncomplete;
    } else if (!result.drift_passed || !result.split_half_passed ||
               !result.relative_standard_error_passed) {
        result.classification =
            BlockConvergenceClassification::NonStationary;
    } else {
        result.classification = BlockConvergenceClassification::Converged;
    }
    return result;
}

std::string to_string(BlockConvergenceClassification classification) {
    switch (classification) {
        case BlockConvergenceClassification::Invalid: return "invalid";
        case BlockConvergenceClassification::HorizonIncomplete:
            return "horizon_incomplete";
        case BlockConvergenceClassification::EffectiveSampleIncomplete:
            return "effective_sample_incomplete";
        case BlockConvergenceClassification::NonStationary:
            return "nonstationary";
        case BlockConvergenceClassification::Converged: return "converged";
    }
    throw std::logic_error("unknown block convergence classification");
}

PeriodicBlockConvergence::PeriodicBlockConvergence(
    std::size_t steps_per_cycle,
    std::size_t cycles_per_block,
    std::size_t origin_step,
    std::vector<std::string> observable_names,
    BlockConvergenceCriteria criteria)
    : criteria_(criteria) {
    state_.steps_per_cycle = steps_per_cycle;
    state_.cycles_per_block = cycles_per_block;
    state_.origin_step = origin_step;
    state_.last_step = origin_step;
    state_.observable_names = std::move(observable_names);
    state_.current_sums.assign(state_.observable_names.size(), 0.0L);
    state_.completed_block_means.resize(state_.observable_names.size());
    validate_state();
}

PeriodicBlockConvergence::PeriodicBlockConvergence(
    PeriodicBlockState state,
    BlockConvergenceCriteria criteria)
    : state_(std::move(state)), criteria_(criteria) {
    validate_state();
}

void PeriodicBlockConvergence::validate_state() const {
    const std::size_t block_samples = samples_per_block();
    const std::size_t elapsed = state_.last_step >= state_.origin_step
        ? state_.last_step - state_.origin_step : 0;
    if (!valid_criteria(criteria_) || state_.steps_per_cycle == 0 ||
        state_.cycles_per_block == 0 ||
        state_.observable_names.empty() ||
        state_.current_sums.size() != state_.observable_names.size() ||
        state_.completed_block_means.size() !=
            state_.observable_names.size() ||
        state_.last_step < state_.origin_step ||
        state_.samples_in_current_block >= block_samples ||
        elapsed / block_samples != completed_blocks() ||
        elapsed % block_samples != state_.samples_in_current_block) {
        throw std::invalid_argument(
            "invalid periodic block convergence state");
    }
    const std::size_t blocks = completed_blocks();
    for (std::size_t i = 0; i < state_.observable_names.size(); ++i) {
        if (state_.observable_names[i].empty() ||
            state_.completed_block_means[i].size() != blocks ||
            !std::isfinite(state_.current_sums[i]) ||
            !std::all_of(state_.completed_block_means[i].begin(),
                         state_.completed_block_means[i].end(),
                         [](double value) { return std::isfinite(value); })) {
            throw std::invalid_argument(
                "invalid periodic block convergence observable state");
        }
        if (std::find(state_.observable_names.begin(),
                      state_.observable_names.begin() + i,
                      state_.observable_names[i]) !=
            state_.observable_names.begin() + i) {
            throw std::invalid_argument(
                "periodic convergence observable names must be unique");
        }
    }
}

bool PeriodicBlockConvergence::observe(
    std::size_t step, const std::vector<double>& values) {
    if (state_.last_step == std::numeric_limits<std::size_t>::max() ||
        step != state_.last_step + 1) {
        throw std::invalid_argument(
            "periodic convergence samples must be step-contiguous");
    }
    if (values.size() != state_.observable_names.size() ||
        !std::all_of(values.begin(), values.end(),
                     [](double value) { return std::isfinite(value); })) {
        throw std::invalid_argument(
            "periodic convergence sample is invalid");
    }
    for (std::size_t i = 0; i < values.size(); ++i) {
        state_.current_sums[i] += static_cast<long double>(values[i]);
    }
    state_.last_step = step;
    ++state_.samples_in_current_block;
    if (state_.samples_in_current_block != samples_per_block()) return false;

    const long double count =
        static_cast<long double>(state_.samples_in_current_block);
    for (std::size_t i = 0; i < values.size(); ++i) {
        state_.completed_block_means[i].push_back(
            static_cast<double>(state_.current_sums[i] / count));
        state_.current_sums[i] = 0.0L;
    }
    state_.samples_in_current_block = 0;
    return true;
}

std::vector<BlockConvergenceResult>
PeriodicBlockConvergence::evaluate() const {
    std::vector<BlockConvergenceResult> result;
    result.reserve(state_.observable_names.size());
    for (const auto& blocks : state_.completed_block_means) {
        result.push_back(evaluate_block_convergence(blocks, criteria_));
    }
    return result;
}

bool PeriodicBlockConvergence::converged() const {
    const auto results = evaluate();
    return !results.empty() &&
        std::all_of(results.begin(), results.end(),
                    [](const BlockConvergenceResult& value) {
                        return value.converged();
                    });
}

std::size_t PeriodicBlockConvergence::completed_blocks() const {
    return state_.completed_block_means.empty()
        ? 0 : state_.completed_block_means.front().size();
}

std::size_t PeriodicBlockConvergence::samples_per_block() const {
    if (state_.steps_per_cycle != 0 &&
        state_.cycles_per_block >
            std::numeric_limits<std::size_t>::max() /
                state_.steps_per_cycle) {
        throw std::invalid_argument(
            "periodic convergence block length overflows");
    }
    return state_.steps_per_cycle * state_.cycles_per_block;
}

} // namespace pic
