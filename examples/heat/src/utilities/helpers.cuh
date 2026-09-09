#pragma once

#include "real.cuh"

#include <cmath>
#include <cstddef>
#include <numbers>

// Per-axis factor of finite-domain Neumann cosine eigenmode at cell-centered index on an n-cell axis (n-2 interior).
CUDA_CALLABLE
inline real_t cosine_mode(std::size_t index, std::size_t n) {
  const real_t arg{
    std::numbers::pi_v<real_t> *
      (static_cast<real_t>(index) - static_cast<real_t>(0.5)) /
      static_cast<real_t>(n - 2)
  };
  return real_t_cos(arg);
}

inline real_t neumann_decay_rate(
  real_t alpha,
  real_t dx, real_t dy, real_t dz,
  std::size_t nx, std::size_t ny, std::size_t nz
) {
  constexpr real_t pi_sq{
    std::numbers::pi_v<real_t> * std::numbers::pi_v<real_t>
  };
  const real_t lx{static_cast<real_t>(nx - 2) * dx};
  const real_t ly{static_cast<real_t>(ny - 2) * dy};
  const real_t lz{static_cast<real_t>(nz - 2) * dz};
  
  return alpha * pi_sq * (
    static_cast<real_t>(1.0) / (lx * lx) +
    static_cast<real_t>(1.0) / (ly * ly) +
    static_cast<real_t>(1.0) / (lz * lz)
  );
}
