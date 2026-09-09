#include "simulation.hpp"
#include "../utilities/helpers.cuh"
#include "../io/vtk_writer.hpp"

#include <cmath>

Simulation::Simulation(const Config& config)
: total_steps_{config.total_steps}
, output_interval_{config.output_interval}
, ic_{config.ic}
, grid_a_{config}
, grid_b_{config}
, integrator_{std::make_unique<ExplicitEuler>(config)}
{ initialize(); }

namespace {

#if defined(__CUDACC__)
__global__
void gpuInitializeGaussianKernel(
  std::size_t nx, std::size_t ny, std::size_t nz,
  std::size_t p_nx, std::size_t p_ny,
  real_t* RESTRICT u
) {
  const std::size_t i{blockIdx.x * blockDim.x + threadIdx.x};
  const std::size_t j{blockIdx.y * blockDim.y + threadIdx.y};
  const std::size_t k{blockIdx.z * blockDim.z + threadIdx.z};

  if (i >= nx || j >= ny || k >= nz) { return; }

  const real_t center_x{static_cast<real_t>(0.5) * static_cast<real_t>(nx - 1)};
  const real_t center_y{static_cast<real_t>(0.5) * static_cast<real_t>(ny - 1)};
  const real_t center_z{static_cast<real_t>(0.5) * static_cast<real_t>(nz - 1)};

  const real_t amplitude{static_cast<real_t>(1.0)};
  const real_t sigma{static_cast<real_t>(0.1) * static_cast<real_t>(nx)};
  const real_t two_sigma_sq{static_cast<real_t>(2.0) * sigma * sigma};
  const real_t inv_two_sig_sq{static_cast<real_t>(1.0) / two_sigma_sq};

  const real_t rx{static_cast<real_t>(i) - center_x};
  const real_t ry{static_cast<real_t>(j) - center_y};
  const real_t rz{static_cast<real_t>(k) - center_z};
  const real_t r_sq{rx*rx + ry*ry + rz*rz};

  const std::size_t point{i + p_nx * (j + p_ny * k)};
  u[point] = amplitude * real_t_exp(-r_sq * inv_two_sig_sq);
}

__global__
void gpuInitializeCosineKernel(
  std::size_t nx, std::size_t ny, std::size_t nz,
  std::size_t p_nx, std::size_t p_ny,
  real_t* RESTRICT u
) {
  const std::size_t i{blockIdx.x * blockDim.x + threadIdx.x};
  const std::size_t j{blockIdx.y * blockDim.y + threadIdx.y};
  const std::size_t k{blockIdx.z * blockDim.z + threadIdx.z};

  if (i >= nx || j >= ny || k >= nz) { return; }

  const std::size_t point{i + p_nx * (j + p_ny * k)};
  u[point] = cosine_mode(i, nx) * cosine_mode(j, ny) * cosine_mode(k, nz);
}
#else
  void cpuInitializeGaussianKernel(
  std::size_t nx, std::size_t ny, std::size_t nz,
  std::size_t p_nx, std::size_t p_ny,
  real_t* RESTRICT u
) {
  ASSUME_ALIGNED(u, SIMD_BYTES);

  const real_t center_x{static_cast<real_t>(0.5) * static_cast<real_t>(nx - 1)};
  const real_t center_y{static_cast<real_t>(0.5) * static_cast<real_t>(ny - 1)};
  const real_t center_z{static_cast<real_t>(0.5) * static_cast<real_t>(nz - 1)};

  const real_t amplitude{static_cast<real_t>(1.0)};
  const real_t sigma{static_cast<real_t>(0.1) * static_cast<real_t>(nx)};
  const real_t two_sigma_sq{static_cast<real_t>(2.0) * sigma * sigma};
  const real_t inv_two_sig_sq{static_cast<real_t>(1.0) / two_sigma_sq};

  #pragma omp parallel for collapse(2)
  for (std::ptrdiff_t k = 0; k < static_cast<std::ptrdiff_t>(nz); ++k) {
    for (std::ptrdiff_t j = 0; j < static_cast<std::ptrdiff_t>(ny); ++j) {

      #pragma omp simd
      for (std::size_t i = 0; i < nx; ++i) {
        const real_t rx{static_cast<real_t>(i) - center_x};
        const real_t ry{static_cast<real_t>(j) - center_y};
        const real_t rz{static_cast<real_t>(k) - center_z};
        const real_t r_sq{rx*rx + ry*ry + rz*rz};

        const std::size_t point{i + p_nx * (j + p_ny * k)};
        u[point] = amplitude * real_t_exp(-r_sq * inv_two_sig_sq);
      }
    }
  }
}

void cpuInitializeCosineKernel(
  std::size_t nx, std::size_t ny, std::size_t nz,
  std::size_t p_nx, std::size_t p_ny,
  real_t* RESTRICT u
) {
  ASSUME_ALIGNED(u, SIMD_BYTES);

  #pragma omp parallel for collapse(2)
  for (std::ptrdiff_t k = 0; k < static_cast<std::ptrdiff_t>(nz); ++k) {
    for (std::ptrdiff_t j = 0; j < static_cast<std::ptrdiff_t>(ny); ++j) {

      const real_t cy{cosine_mode(static_cast<std::size_t>(j), ny)};
      const real_t cz{cosine_mode(static_cast<std::size_t>(k), nz)};

      #pragma omp simd
      for (std::size_t i = 0; i < nx; ++i) {
        const std::size_t point{i + p_nx * (j + p_ny * k)};
        u[point] = cosine_mode(i, nx) * cy * cz;
      }
    }
  }
}
#endif

} // namespace

void Simulation::initialize() {
#if defined(__CUDACC__)
  dim3 threads(8, 8, 8);

  if (ic_ == InitCondition::NeumannCosine) {
    dim3 blocks(
      (static_cast<unsigned>(grid_a_.nx()) + threads.x - 1) / threads.x,
      (static_cast<unsigned>(grid_a_.ny()) + threads.y - 1) / threads.y,
      (static_cast<unsigned>(grid_a_.nz()) + threads.z - 1) / threads.z
    );

    gpuInitializeCosineKernel<<<blocks, threads>>>(
      grid_a_.nx(), grid_a_.ny(), grid_a_.nz(),
      grid_a_.p_nx(), grid_a_.p_ny(),
      grid_a_.field()
    );
    CUDA_CHECK(cudaGetLastError());

    return;
  }

  dim3 blocks(
    (static_cast<unsigned>(grid_a_.nx()) + threads.x - 1) / threads.x,
    (static_cast<unsigned>(grid_a_.ny()) + threads.y - 1) / threads.y,
    (static_cast<unsigned>(grid_a_.nz()) + threads.z - 1) / threads.z
  );

  gpuInitializeGaussianKernel<<<blocks, threads>>>(
    grid_a_.nx(), grid_a_.ny(), grid_a_.nz(),
    grid_a_.p_nx(), grid_a_.p_ny(),
    grid_a_.field()
  );
  CUDA_CHECK(cudaGetLastError());
#else
  if (ic_ == InitCondition::NeumannCosine) {
    cpuInitializeCosineKernel(
      grid_a_.nx(), grid_a_.ny(), grid_a_.nz(),
      grid_a_.p_nx(), grid_a_.p_ny(),
      grid_a_.field()
    );
    return;
  }

  cpuInitializeGaussianKernel(
    grid_a_.nx(), grid_a_.ny(), grid_a_.nz(),
    grid_a_.p_nx(), grid_a_.p_ny(),
    grid_a_.field()
  );
#endif
}

void Simulation::run() {
  auto curr_grid{current_grid_};
  auto next_grid{(curr_grid == &grid_a_) ? &grid_b_ : &grid_a_};

  const bool enable_vtk{output_interval_ > 0};
  for (std::size_t step{}; step < total_steps_; ++step) {
    const bool output{
      enable_vtk           &&
      output_interval_ > 0 &&
      step % output_interval_ == 0
    };

    if (output) { vtk::write(*curr_grid, step); }

    integrator_->integrate(*curr_grid, *next_grid);
    std::swap(curr_grid, next_grid);
  }
  if (enable_vtk) { vtk::write(*curr_grid, total_steps_); }

  #if defined(__CUDACC__)
    CUDA_CHECK(cudaDeviceSynchronize());
  #endif

  current_grid_ = curr_grid;
}