#pragma once

#include "../grid/grid.cuh"
#include "../integrator/integrator.cuh"

#include <utility>
#include <memory>

class Simulation {
private:
  std::size_t total_steps_;
  std::size_t output_interval_;
  InitCondition ic_;

  Grid grid_a_;
  Grid grid_b_;
  Grid* current_grid_{&grid_a_};
  std::unique_ptr<Integrator> integrator_;

public:
  Simulation(const Config& config);

  void run();

  [[nodiscard]] const Grid& grid() const { return *current_grid_; }

private:
  void initialize();
};