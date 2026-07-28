#!/usr/bin/env sh
set -eu

# Keep local verification responsive by default. Override deliberately on a
# dedicated build host, for example AURORA_BUILD_JOBS=8 ./scripts/verify.sh.
BUILD_JOBS=${AURORA_BUILD_JOBS:-1}
TEST_JOBS=${AURORA_TEST_JOBS:-1}
OPENMP_THREADS=${AURORA_OPENMP_THREADS:-1}
case "$BUILD_JOBS:$TEST_JOBS:$OPENMP_THREADS" in
  *[!0-9:]*|0:*|*:0:*|*:0)
    echo "AURORA_BUILD_JOBS, AURORA_TEST_JOBS, and AURORA_OPENMP_THREADS must be positive integers" >&2
    exit 2
    ;;
esac
export CMAKE_BUILD_PARALLEL_LEVEL=$BUILD_JOBS
export CTEST_PARALLEL_LEVEL=$TEST_JOBS
export OMP_NUM_THREADS=$OPENMP_THREADS

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel "$BUILD_JOBS"
python3 scripts/validate_milestones.py
python3 scripts/validate_release_artifacts.py
ctest --test-dir build --parallel "$TEST_JOBS" --output-on-failure
python3 scripts/validate_pushers.py
python3 scripts/verify_examples.py build/aurorapic_cli
python3 scripts/benchmark_unstructured.py build/aurorapic_cli --repeats 1
python3 scripts/verify_install_package.py build --jobs "$BUILD_JOBS"
