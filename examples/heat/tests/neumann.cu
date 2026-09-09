#include "config/config.hpp"
#include "simulation/simulation.hpp"
#include "cuda_test_guard.cuh"
#include "utilities/helpers.cuh"

#include <cmath>
#include <cstdio>
#include <vector>

namespace {

std::vector<real_t> field_snapshot(const Grid& grid) {
  std::vector<real_t> field(grid.p_nx() * grid.p_ny() * grid.p_nz());
  grid.copy_to_host(field.data());
  return field;
}

} // namespace

int main() {
  HEAT_SOLVER_SKIP_CUDA_TEST_IF_UNAVAILABLE();

  Config cfg{};
  cfg.nx = cfg.ny = cfg.nz = 64;
  cfg.dx = cfg.dy = cfg.dz = static_cast<real_t>(1.0);
  cfg.alpha = static_cast<real_t>(1.0);
  cfg.ic = InitCondition::NeumannCosine;
  cfg.output_interval = 0;
  cfg.total_steps = 200;
  cfg.dt = Config::stable_dt(cfg.alpha, cfg.dx, cfg.dy, cfg.dz);

  Simulation sim{cfg};
  sim.run();
  const Grid& grid{sim.grid()};
  const std::vector<real_t> field{field_snapshot(grid)};

  // Exact finite-box solution: u = prod cos(pi r/L) * exp(-lambda t), du/dn=0 on every face.
  const real_t t{static_cast<real_t>(cfg.total_steps) * cfg.dt};
  const real_t lambda{neumann_decay_rate(cfg.alpha, cfg.dx, cfg.dy, cfg.dz, cfg.nx, cfg.ny, cfg.nz)};
  const real_t decay{std::exp(-lambda * t)};

  // Check 1: global L2 relative error against the analytic solution over the interior.
  double sq_error_sum{};
  double sq_expected_sum{};
  for (std::size_t k{1}; k < cfg.nz - 1; ++k) {

    const real_t cz{cosine_mode(k, cfg.nz)};
    for (std::size_t j{1}; j < cfg.ny - 1; ++j) {

      const real_t cy{cosine_mode(j, cfg.ny)};
      for (std::size_t i{1}; i < cfg.nx - 1; ++i) {
        const real_t cx{cosine_mode(i, cfg.nx)};
        const real_t expected{cx * cy * cz * decay};
        const real_t actual{field[grid.idx(i, j, k)]};
        const double diff{static_cast<double>(actual) - static_cast<double>(expected)};

        sq_error_sum += diff * diff;
        sq_expected_sum += static_cast<double>(expected) * static_cast<double>(expected);
      }
    }
  }

  const double l2_rel_error{std::sqrt(sq_error_sum / sq_expected_sum)};

  if (l2_rel_error > 5e-4) {
    std::fprintf(stderr,
      "neumann: L2 relative error too large (%.3e)\n", l2_rel_error);
    return 1;
  }

  // Check 2: du/dn = 0 still holds exactly at the faces under the cosine condition.
  for (std::size_t k{1}; k < cfg.nz - 1; ++k) {
    for (std::size_t j{1}; j < cfg.ny - 1; ++j) {
      if (field[grid.idx(0, j, k)] != field[grid.idx(1, j, k)] ||
          field[grid.idx(cfg.nx - 1, j, k)] != field[grid.idx(cfg.nx - 2, j, k)]) {
            
        std::fprintf(stderr, "neumann: x-boundary not insulated\n");
        return 1;
      }
    }
  }

  return 0;
}
