#pragma once

#include "../utilities/real.cuh"

#include <cstddef>
#include <cstdint>

enum class InitCondition { Gaussian, NeumannCosine };

struct Config {
  std::size_t total_steps{1000};
  std::size_t output_interval{0};
  InitCondition ic{InitCondition::Gaussian};

  std::size_t nx{64};
  std::size_t ny{64};
  std::size_t nz{64};

  real_t dx{static_cast<real_t>(1.0)};
  real_t dy{static_cast<real_t>(1.0)};
  real_t dz{static_cast<real_t>(1.0)};

  real_t alpha{static_cast<real_t>(1.0)};
  real_t dt{stable_dt(alpha, dx, dy, dz)};

  static constexpr real_t stable_dt(real_t alpha, real_t dx, real_t dy, real_t dz) {
    return static_cast<real_t>(0.85) / (
      static_cast<real_t>(2.0) * alpha *
      (static_cast<real_t>(1.0)/(dx*dx) +
       static_cast<real_t>(1.0)/(dy*dy) +
       static_cast<real_t>(1.0)/(dz*dz))
    );
  }

  static Config parse(int argc, char** argv);
};