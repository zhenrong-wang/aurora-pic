#!/usr/bin/env sh
set -eu
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
python3 scripts/validate_milestones.py
python3 scripts/validate_release_artifacts.py
ctest --test-dir build --output-on-failure
python3 scripts/validate_pushers.py
python3 scripts/verify_examples.py build/aurorapic_cli
