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
python3 scripts/validate_hall_case.py
python3 scripts/test_import_lxcat.py
python3 scripts/test_swarm_cli.py build/aurorapic_swarm
python3 scripts/test_compare_swarm.py
python3 scripts/test_swarm_campaign.py
python3 scripts/test_hall_comparison.py
python3 scripts/test_hall_ensemble.py
python3 scripts/test_hall_convergence.py build/aurorapic_cli
python3 scripts/test_hall_flux_stationarity.py
python3 scripts/test_hall_flux_block_comparison.py
python3 scripts/test_hall_horizon_stage.py build/aurorapic_cli
python3 scripts/test_lock_hall_source.py
python3 scripts/test_normalize_hall_reference.py
python3 scripts/test_hall_pilot.py build/aurorapic_cli
python3 scripts/test_hall_runtime_qualification.py build/aurorapic_cli
python3 scripts/test_turner_balance.py
python3 scripts/test_compare_turner.py
python3 scripts/test_turner_density_blocks.py
python3 scripts/test_prepare_turner_ensemble.py
python3 scripts/test_attach_turner_ensemble_result.py
python3 scripts/test_analyze_turner_ensemble.py
python3 scripts/test_prepare_turner_sensitivity.py
python3 scripts/test_spatial_average_1d.py build/aurorapic_cli
python3 scripts/test_analyze_turner_spatial_structure.py
python3 scripts/test_phase_eedf_interchange.py
python3 scripts/test_edupic_stage.py
PYTHONPATH=scripts python3 scripts/test_analyze_edupic_convergence.py
PYTHONPATH=scripts python3 scripts/test_advance_edupic_equilibration.py
ctest --test-dir build --parallel "$TEST_JOBS" --output-on-failure
python3 scripts/validate_pushers.py
python3 scripts/test_kinetic_benchmarks.py
python3 scripts/validate_kinetic_benchmarks.py build/aurorapic_cli
python3 scripts/verify_examples.py build/aurorapic_cli
python3 scripts/benchmark_unstructured.py build/aurorapic_cli --repeats 1
python3 scripts/verify_install_package.py build --jobs "$BUILD_JOBS"
