#pragma once

#include <algorithm>
#include <cstddef>
#include <limits>
#include <string>
namespace pic {

enum class RuntimeBackend { Serial, OpenMP, MPI, GPU };

struct RuntimePolicy {
    RuntimeBackend backend{RuntimeBackend::Serial};
    std::size_t threads{1};
};

struct RuntimeInfo {
    RuntimeBackend backend{RuntimeBackend::Serial};
    std::size_t requested_threads{1};
    std::size_t active_threads{1};
    bool openmp_compiled{false};
    bool mpi_compiled{false};
    bool gpu_compiled{false};
};

std::string to_string(RuntimeBackend backend);
RuntimeInfo runtime_info(const RuntimePolicy& policy);
void validate_runtime_policy(const RuntimePolicy& policy);

template <typename Body>
void runtime_parallel_for(std::size_t begin, std::size_t end, const RuntimePolicy& policy, Body&& body) {
    validate_runtime_policy(policy);
    if (begin >= end) return;
#ifdef AURORA_HAVE_OPENMP
    if (policy.backend == RuntimeBackend::OpenMP && policy.threads > 1) {
        const int requested_threads = static_cast<int>(policy.threads);
        #pragma omp parallel for schedule(static) num_threads(requested_threads)
        for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(begin); i < static_cast<std::ptrdiff_t>(end); ++i) {
            body(static_cast<std::size_t>(i));
        }
        return;
    }
#else
    (void)policy;
#endif
    for (std::size_t i = begin; i < end; ++i) body(i);
}

template <typename Body>
void runtime_static_chunks(std::size_t begin, std::size_t end, const RuntimePolicy& policy, Body&& body) {
    validate_runtime_policy(policy);
    if (begin >= end) return;
    const std::size_t active_threads = runtime_info(policy).active_threads;
    const std::size_t total = end - begin;
    const std::size_t chunk_count = active_threads == 0 ? 1 : active_threads;
    const std::size_t chunk_size = (total + chunk_count - 1) / chunk_count;
    runtime_parallel_for(std::size_t{0}, chunk_count, policy, [&](std::size_t chunk) {
        const std::size_t chunk_begin = begin + chunk * chunk_size;
        const std::size_t chunk_end = std::min(end, chunk_begin + chunk_size);
        if (chunk_begin < chunk_end) body(chunk, chunk_begin, chunk_end);
    });
}

} // namespace pic
