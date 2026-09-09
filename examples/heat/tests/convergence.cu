#include "config/config.hpp"
#include "simulation/simulation.hpp"
#include "cuda_test_guard.cuh"
#include "utilities/helpers.cuh"

#include <cmath>
#include <cstdio>
#include <vector>

int main(int argc, char** argv) {
  HEAT_SOLVER_SKIP_CUDA_TEST_IF_UNAVAILABLE();

  const Config cfg{Config::parse(argc, argv)};

  Simulation sim{cfg};
  sim.run();

  const Grid& grid{sim.grid()};
  std::vector<real_t> field(grid.total_size());
  grid.copy_to_host(field.data());

  if (cfg.ic == InitCondition::NeumannCosine) {
    const real_t t{static_cast<real_t>(cfg.total_steps) * cfg.dt};
    const real_t lambda{neumann_decay_rate(cfg.alpha, cfg.dx, cfg.dy, cfg.dz, cfg.nx, cfg.ny, cfg.nz)};
    const real_t decay{std::exp(-lambda * t)};

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
    std::printf(
      "%zu,%.6f,%zu,%.8f,%.6f,%.10f,%.10f\n",
      cfg.nx, cfg.dx, cfg.total_steps, cfg.dt, t, l2_rel_error, 0.0
    );
    return 0;
  }

  const real_t center_x{static_cast<real_t>(0.5) * static_cast<real_t>(cfg.nx - 1)};
  const real_t center_y{static_cast<real_t>(0.5) * static_cast<real_t>(cfg.ny - 1)};
  const real_t center_z{static_cast<real_t>(0.5) * static_cast<real_t>(cfg.nz - 1)};

  const real_t sigma0{static_cast<real_t>(0.1) * static_cast<real_t>(cfg.nx) * cfg.dx};
  const real_t t{static_cast<real_t>(cfg.total_steps) * cfg.dt};
  const real_t sigma_t_sq{sigma0 * sigma0 + static_cast<real_t>(2.0) * cfg.alpha * t};
  const real_t amplitude_t{std::pow(sigma0 * sigma0 / sigma_t_sq, static_cast<real_t>(1.5))};

  double sq_error_sum{};
  double sq_expected_sum{};

  for (std::size_t k{}; k < cfg.nz; ++k) {
    for (std::size_t j{}; j < cfg.ny; ++j) {
      for (std::size_t i{}; i < cfg.nx; ++i) {
        const real_t rx{(static_cast<real_t>(i) - center_x) * cfg.dx};
        const real_t ry{(static_cast<real_t>(j) - center_y) * cfg.dy};
        const real_t rz{(static_cast<real_t>(k) - center_z) * cfg.dz};
        const real_t r_sq{rx * rx + ry * ry + rz * rz};
        const real_t expected{amplitude_t * std::exp(-r_sq / (static_cast<real_t>(2.0) * sigma_t_sq))};
        const real_t actual{field[grid.idx(i, j, k)]};

        const double diff{static_cast<double>(actual) - static_cast<double>(expected)};
        sq_error_sum += diff * diff;
        sq_expected_sum += static_cast<double>(expected) * static_cast<double>(expected);
      }
    }
  }

  const double l2_rel_error{std::sqrt(sq_error_sum / sq_expected_sum)};

  const std::size_t cx{cfg.nx / 2};
  const std::size_t cy{cfg.ny / 2};
  const std::size_t cz{cfg.nz / 2};
  const real_t rx_c{(static_cast<real_t>(cx) - center_x) * cfg.dx};
  const real_t ry_c{(static_cast<real_t>(cy) - center_y) * cfg.dy};
  const real_t rz_c{(static_cast<real_t>(cz) - center_z) * cfg.dz};
  const real_t r_sq_c{rx_c * rx_c + ry_c * ry_c + rz_c * rz_c};
  const real_t expected_center{amplitude_t * std::exp(-r_sq_c / (static_cast<real_t>(2.0) * sigma_t_sq))};
  const real_t actual_center{field[grid.idx(cx, cy, cz)]};
  const double center_rel_error{
    std::fabs(static_cast<double>(actual_center) - static_cast<double>(expected_center))
    / static_cast<double>(expected_center)
  };

  std::printf(
    "%zu,%.6f,%zu,%.8f,%.6f,%.10f,%.10f\n",
    cfg.nx, cfg.dx, cfg.total_steps, cfg.dt, t, l2_rel_error, center_rel_error
  );

  return 0;
}
