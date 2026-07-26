#!/usr/bin/env sh
set -eu
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
./build/aurorapic_cli examples/two_stream.cfg
./build/aurorapic_cli examples/sheath_steady.cfg
./build/aurorapic_cli examples/plasma_2d.cfg
./build/aurorapic_cli examples/electrode_2d.cfg
