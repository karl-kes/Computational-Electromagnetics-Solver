#include "config/config.hpp"
#include "simulation/simulation.hpp"
#include "cuda_test_guard.cuh"

#include <cmath>
#include <cstdio>
#include <vector>

namespace {

bool nearly_equal(real_t actual, real_t expected, real_t rel_tol) {
  const real_t diff{std::fabs(actual - expected)};
  const real_t scale{std::fabs(expected) > static_cast<real_t>(1e-6) ? std::fabs(expected) : static_cast<real_t>(1.0)};
  return diff <= rel_tol * scale;
}

std::vector<real_t> field_snapshot(const Grid& grid) {
  std::vector<real_t> field(grid.p_nx() * grid.p_ny() * grid.p_nz());
  grid.copy_to_host(field.data());
  return field;
}

} // namespace

int main() {
  HEAT_SOLVER_SKIP_CUDA_TEST_IF_UNAVAILABLE();

  Config cfg{};
  cfg.nx = cfg.ny = cfg.nz = 32;
  cfg.dx = cfg.dy = cfg.dz = static_cast<real_t>(1.0);
  cfg.alpha = static_cast<real_t>(1.0);
  cfg.output_interval = 0;
  cfg.total_steps = 0;
  cfg.dt = Config::stable_dt(cfg.alpha, cfg.dx, cfg.dy, cfg.dz);

  const std::size_t total_steps{40};

  Simulation sim0{cfg};
  const std::vector<real_t> field0{field_snapshot(sim0.grid())};
  const Grid& grid0{sim0.grid()};

  cfg.total_steps = total_steps;
  Simulation sim{cfg};
  sim.run();
  const std::vector<real_t> field{field_snapshot(sim.grid())};
  const Grid& grid{sim.grid()};

  // Check 1: center value vs. u(r,t) = amplitude * (sigma0/sigma_t)^3 * exp(-r^2/(2*sigma_t^2))
  const real_t center_x{static_cast<real_t>(0.5) * static_cast<real_t>(cfg.nx - 1)};
  const real_t center_y{static_cast<real_t>(0.5) * static_cast<real_t>(cfg.ny - 1)};
  const real_t center_z{static_cast<real_t>(0.5) * static_cast<real_t>(cfg.nz - 1)};

  const std::size_t cx{cfg.nx / 2};
  const std::size_t cy{cfg.ny / 2};
  const std::size_t cz{cfg.nz / 2};

  const real_t rx{static_cast<real_t>(cx) - center_x};
  const real_t ry{static_cast<real_t>(cy) - center_y};
  const real_t rz{static_cast<real_t>(cz) - center_z};
  const real_t r_sq{rx * rx + ry * ry + rz * rz};

  const real_t sigma0{static_cast<real_t>(0.1) * static_cast<real_t>(cfg.nx)};
  const real_t t{static_cast<real_t>(total_steps) * cfg.dt};
  const real_t sigma_t_sq{sigma0 * sigma0 + static_cast<real_t>(2.0) * cfg.alpha * t};
  const real_t amplitude_t{std::pow(sigma0 * sigma0 / sigma_t_sq, static_cast<real_t>(1.5))};
  const real_t expected{amplitude_t * std::exp(-r_sq / (static_cast<real_t>(2.0) * sigma_t_sq))};

  const real_t actual{field[grid.idx(cx, cy, cz)]};

  if (!nearly_equal(actual, expected, static_cast<real_t>(0.05))) {
    std::fprintf(stderr,
      "physical_correctness: center value mismatch (expected %.6f, got %.6f)\n",
      expected, actual);
    return 1;
  }

  // Check 2: conservation of energy over interior cells (ghost layer excluded)
  real_t sum0{};
  real_t sum{};
  for (std::size_t k{1}; k < cfg.nz - 1; ++k) {
    for (std::size_t j{1}; j < cfg.ny - 1; ++j) {
      for (std::size_t i{1}; i < cfg.nx - 1; ++i) {
        sum0 += field0[grid0.idx(i, j, k)];
        sum += field[grid.idx(i, j, k)];
      }
    }
  }

  if (!nearly_equal(sum, sum0, static_cast<real_t>(0.02))) {
    std::fprintf(stderr,
      "physical_correctness: sum not conserved (t=0 sum %.6f, final sum %.6f)\n",
      sum0, sum);
    return 1;
  }

  // Check 3: du/dn = 0 -- ghost cells must exactly mirror the adjacent interior cell.
  for (std::size_t k{1}; k < cfg.nz - 1; ++k) {
    for (std::size_t j{1}; j < cfg.ny - 1; ++j) {
      const real_t left_ghost{field[grid.idx(0, j, k)]};
      const real_t left_inner{field[grid.idx(1, j, k)]};
      const real_t right_ghost{field[grid.idx(cfg.nx - 1, j, k)]};
      const real_t right_inner{field[grid.idx(cfg.nx - 2, j, k)]};

      if (left_ghost != left_inner || right_ghost != right_inner) {
        std::fprintf(stderr,
          "physical_correctness: x-boundary not insulated at (j=%zu, k=%zu)\n", j, k);
        return 1;
      }
    }
  }

  for (std::size_t k{1}; k < cfg.nz - 1; ++k) {
    for (std::size_t i{}; i < cfg.nx; ++i) {
      const real_t front_ghost{field[grid.idx(i, 0, k)]};
      const real_t front_inner{field[grid.idx(i, 1, k)]};
      const real_t back_ghost{field[grid.idx(i, cfg.ny - 1, k)]};
      const real_t back_inner{field[grid.idx(i, cfg.ny - 2, k)]};

      if (front_ghost != front_inner || back_ghost != back_inner) {
        std::fprintf(stderr,
          "physical_correctness: y-boundary not insulated at (i=%zu, k=%zu)\n", i, k);
        return 1;
      }
    }
  }

  for (std::size_t j{}; j < cfg.ny; ++j) {
    for (std::size_t i{}; i < cfg.nx; ++i) {
      const real_t bottom_ghost{field[grid.idx(i, j, 0)]};
      const real_t bottom_inner{field[grid.idx(i, j, 1)]};
      const real_t top_ghost{field[grid.idx(i, j, cfg.nz - 1)]};
      const real_t top_inner{field[grid.idx(i, j, cfg.nz - 2)]};

      if (bottom_ghost != bottom_inner || top_ghost != top_inner) {
        std::fprintf(stderr,
          "physical_correctness: z-boundary not insulated at (i=%zu, j=%zu)\n", i, j);
        return 1;
      }
    }
  }

  return 0;
}
