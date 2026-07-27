#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
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

} // namespace pic
