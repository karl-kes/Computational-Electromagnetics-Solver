#pragma once

#include "utilities/macros.cuh"

#if defined(__CUDACC__)
  #include <cstdio>
  #include <cuda_runtime.h>

inline constexpr int cuda_test_skip_return_code{77};

inline bool cuda_test_device_available() {
  int device_count{};
  const cudaError_t error{cudaGetDeviceCount(&device_count)};
  if (error != cudaSuccess || device_count <= 0) {
    std::fprintf(stderr, "CUDA test skipped: %s\n", cudaGetErrorString(error));
    return false;
  }

  CUDA_CHECK(cudaSetDevice(0));
  return true;
}

  #define HEAT_SOLVER_SKIP_CUDA_TEST_IF_UNAVAILABLE() \
    do { \
      if (!cuda_test_device_available()) { return cuda_test_skip_return_code; } \
    } while (0)
#else
  #define HEAT_SOLVER_SKIP_CUDA_TEST_IF_UNAVAILABLE() static_cast<void>(0)
#endif
