#include "pic/Runtime.hpp"

#ifdef AURORA_HAVE_OPENMP
#include <omp.h>
#endif

#include <limits>
#include <stdexcept>

namespace pic {

std::string to_string(RuntimeBackend backend) {
    switch (backend) {
        case RuntimeBackend::Serial: return "serial";
        case RuntimeBackend::OpenMP: return "openmp";
        case RuntimeBackend::MPI: return "mpi";
        case RuntimeBackend::GPU: return "gpu";
    }
    return "unknown";
}

RuntimeInfo runtime_info(const RuntimePolicy& policy) {
    RuntimeInfo info;
    info.backend = policy.backend;
    info.requested_threads = policy.threads == 0 ? 1 : policy.threads;
#ifdef AURORA_HAVE_OPENMP
    info.openmp_compiled = true;
#endif
    info.mpi_compiled = false;
    info.gpu_compiled = false;

    switch (policy.backend) {
        case RuntimeBackend::Serial:
            info.active_threads = 1;
            break;
        case RuntimeBackend::OpenMP:
#ifdef AURORA_HAVE_OPENMP
            info.active_threads = info.requested_threads;
#else
            info.active_threads = 1;
#endif
            break;
        case RuntimeBackend::MPI:
        case RuntimeBackend::GPU:
            info.active_threads = 1;
            break;
    }
    return info;
}

void validate_runtime_policy(const RuntimePolicy& policy) {
    if (policy.threads == 0) {
        throw std::invalid_argument("runtime_threads must be positive");
    }
    if (policy.threads > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument("runtime_threads exceeds supported thread-count range");
    }
    switch (policy.backend) {
        case RuntimeBackend::Serial:
            if (policy.threads != 1) {
                throw std::invalid_argument("serial runtime_backend requires runtime_threads = 1");
            }
            return;
        case RuntimeBackend::OpenMP:
#ifndef AURORA_HAVE_OPENMP
            if (policy.threads > 1) {
                throw std::invalid_argument("OpenMP runtime_backend requested but AuroraPIC was built without OpenMP support");
            }
#endif
            return;
        case RuntimeBackend::MPI:
            throw std::invalid_argument("MPI runtime_backend is reserved for future distributed runs and is not implemented yet");
        case RuntimeBackend::GPU:
            throw std::invalid_argument("GPU runtime_backend is reserved for future accelerator runs and is not implemented yet");
    }
}

} // namespace pic
