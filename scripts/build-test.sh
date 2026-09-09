#!/usr/bin/env bash
set -e
mkdir -p build-test
cd build-test
cmake .. -DCMAKE_BUILD_TYPE=Test -DCMAKE_CXX_COMPILER=g++-15
cmake --build . --parallel "$(nproc)"
ctest --output-on-failure
