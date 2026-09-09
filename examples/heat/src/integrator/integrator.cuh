#pragma once

#include "../grid/grid.cuh"
#include "../config/config.hpp"
#include "../utilities/macros.cuh"

class Integrator {
private:
  real_t dt_;
  real_t alpha_;

public:
  explicit Integrator(const Config& config);
  virtual ~Integrator() = default;

  virtual void integrate(const Grid& old_grid, Grid& new_grid) = 0;

  real_t dt() const { return dt_; }
  real_t alpha() const { return alpha_; }
};

class ExplicitEuler : public Integrator {
public:
  ExplicitEuler(const Config& config);
  void integrate(const Grid& old_grid, Grid& new_grid) override;
};