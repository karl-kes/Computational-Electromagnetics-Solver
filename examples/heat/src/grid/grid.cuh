#pragma once

#include "../utilities/aligned_soa.cuh"
#include "../utilities/macros.cuh"
#include "../config/config.hpp"

#include <cstdint>

class Grid {
private:
  std::size_t nx_, ny_, nz_;
  std::size_t p_nx_, p_ny_, p_nz_;
  real_t dx_, dy_, dz_;
  real_t inv_dx_sq_, inv_dy_sq_, inv_dz_sq_;

  AlignedSoA<real_t> data_;

  enum : std::size_t {
    U, NUM_SUB_ARR
  };

public:
  Grid(const Config& config);
  ~Grid();

  [[nodiscard]] real_t* field() { return data_[U]; }
  [[nodiscard]] const real_t* field() const { return data_[U]; }

  void copy_to_host(real_t* dst) const;

  [[nodiscard]] std::size_t nx() const { return nx_; }
  [[nodiscard]] std::size_t ny() const { return ny_; }
  [[nodiscard]] std::size_t nz() const { return nz_; }

  [[nodiscard]] std::size_t p_nx() const { return p_nx_; }
  [[nodiscard]] std::size_t p_ny() const { return p_ny_; }
  [[nodiscard]] std::size_t p_nz() const { return p_nz_; }
  
  [[nodiscard]] std::size_t total_size() const { return p_nx_*p_ny_*p_nz_; }

  [[nodiscard]] real_t dx() const { return dx_; }
  [[nodiscard]] real_t dy() const { return dy_; }
  [[nodiscard]] real_t dz() const { return dz_; }

  [[nodiscard]] real_t inv_dx_sq() const { return inv_dx_sq_; }
  [[nodiscard]] real_t inv_dy_sq() const { return inv_dy_sq_; }
  [[nodiscard]] real_t inv_dz_sq() const { return inv_dz_sq_; }

  [[nodiscard]]
  std::size_t idx(std::size_t x, std::size_t y, std::size_t z) const {
    return x + p_nx_ * (y + p_ny_ * z);
  }
};