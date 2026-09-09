#include "config/config.hpp"
#include "cuda_test_guard.cuh"

#include <vector>

namespace {

Config parse(std::vector<const char*> args) {
  std::vector<char*> argv;
  argv.push_back(const_cast<char*>("heat_solver"));
  for (const char* arg : args) { argv.push_back(const_cast<char*>(arg)); }

  return Config::parse(static_cast<int>(argv.size()), argv.data());
}

} // namespace

int main() {
  HEAT_SOLVER_SKIP_CUDA_TEST_IF_UNAVAILABLE();

  // no arguments leaves the compile-time defaults untouched
  {
    const Config cfg{parse({})};
    if (cfg.nx != 64 || cfg.ny != 64 || cfg.nz != 64) { return 1; }
    if (cfg.total_steps != 1000 || cfg.output_interval != 0) { return 1; }
  }

  // explicit flags override the defaults
  {
    const Config cfg{parse({"--nx", "16", "--steps", "5"})};
    if (cfg.nx != 16 || cfg.total_steps != 5) { return 1; }
    if (cfg.ny != 64 || cfg.nz != 64) { return 1; }
  }

  // dt is always re-derived from alpha/dx/dy/dz, never left stale
  {
    const Config cfg{parse({"--alpha", "2.0"})};
    const real_t expected_dt{Config::stable_dt(cfg.alpha, cfg.dx, cfg.dy, cfg.dz)};
    if (cfg.dt != expected_dt) { return 1; }
  }

  return 0;
}
