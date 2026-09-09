#include "config/config.hpp"
#include "simulation/simulation.hpp"
#include "cuda_test_guard.cuh"
#include "utilities/helpers.cuh"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

std::vector<real_t> field_snapshot(const Grid& grid) {
  std::vector<real_t> field(grid.total_size());
  grid.copy_to_host(field.data());
  return field;
}

bool nearly_equal(real_t actual, real_t expected) {
  const double diff{std::fabs(static_cast<double>(actual) - static_cast<double>(expected))};
  const double scale{std::max(1.0, std::fabs(static_cast<double>(expected)))};
  return diff <= 1e-6 * scale;
}

} // namespace

int main() {
  HEAT_SOLVER_SKIP_CUDA_TEST_IF_UNAVAILABLE();

  Config cfg{};
  cfg.nx = 17;
  cfg.ny = 19;
  cfg.nz = 21;
  cfg.total_steps = 0;
  cfg.output_interval = 0;

  {
    cfg.ic = InitCondition::Gaussian;
    Simulation sim{cfg};
    const Grid& grid{sim.grid()};
    const std::vector<real_t> field{field_snapshot(grid)};

    const real_t center_x{static_cast<real_t>(0.5) * static_cast<real_t>(cfg.nx - 1)};
    const real_t center_y{static_cast<real_t>(0.5) * static_cast<real_t>(cfg.ny - 1)};
    const real_t center_z{static_cast<real_t>(0.5) * static_cast<real_t>(cfg.nz - 1)};
    const real_t sigma{static_cast<real_t>(0.1) * static_cast<real_t>(cfg.nx)};
    const real_t inv_two_sigma_sq{static_cast<real_t>(1.0) / (static_cast<real_t>(2.0) * sigma * sigma)};

    for (std::size_t k{}; k < cfg.nz; ++k) {
      for (std::size_t j{}; j < cfg.ny; ++j) {
        for (std::size_t i{}; i < cfg.nx; ++i) {
          const real_t rx{static_cast<real_t>(i) - center_x};
          const real_t ry{static_cast<real_t>(j) - center_y};
          const real_t rz{static_cast<real_t>(k) - center_z};
          const real_t expected{real_t_exp(-(rx * rx + ry * ry + rz * rz) * inv_two_sigma_sq)};
          const real_t actual{field[grid.idx(i, j, k)]};

          if (!nearly_equal(actual, expected)) {
            std::fprintf(stderr, "gaussian initialization mismatch at (%zu,%zu,%zu)\n", i, j, k);
            return 1;
          }
        }
      }
    }
  }

  {
    cfg.ic = InitCondition::NeumannCosine;
    Simulation sim{cfg};
    const Grid& grid{sim.grid()};
    const std::vector<real_t> field{field_snapshot(grid)};

    for (std::size_t k{}; k < cfg.nz; ++k) {
      for (std::size_t j{}; j < cfg.ny; ++j) {
        for (std::size_t i{}; i < cfg.nx; ++i) {
          const real_t expected{cosine_mode(i, cfg.nx) * cosine_mode(j, cfg.ny) * cosine_mode(k, cfg.nz)};
          const real_t actual{field[grid.idx(i, j, k)]};

          if (!nearly_equal(actual, expected)) {
            std::fprintf(stderr, "cosine initialization mismatch at (%zu,%zu,%zu)\n", i, j, k);
            return 1;
          }
        }
      }
    }
  }

  return 0;
}
