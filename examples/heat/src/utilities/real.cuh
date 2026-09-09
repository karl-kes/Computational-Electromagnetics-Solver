#pragma once

#include "macros.cuh"

#include <cmath>
#include <type_traits>

#if defined(HEAT_SOLVER_PRECISION_DOUBLE)
using real_t = double;
#else
using real_t = float;
#endif

CUDA_CALLABLE
inline real_t real_t_exp(real_t value) {
#if defined(__CUDACC__)
  #if defined(HEAT_SOLVER_PRECISION_DOUBLE)
    return exp(value);
  #else
    return expf(value);
  #endif
#else
  return std::exp(value);
#endif
}

CUDA_CALLABLE
inline real_t real_t_cos(real_t value) {
#if defined(__CUDACC__)
  #if defined(HEAT_SOLVER_PRECISION_DOUBLE)
    return cos(value);
  #else
    return cosf(value);
  #endif
#else
  return std::cos(value);
#endif
}
