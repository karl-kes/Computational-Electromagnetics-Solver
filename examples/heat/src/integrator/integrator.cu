#include "integrator.cuh"

#include <omp.h>

Integrator::Integrator(const Config& config)
: dt_{config.dt}
, alpha_{config.alpha}
{ }

ExplicitEuler::ExplicitEuler(const Config& config)
: Integrator(config)
{ }

namespace {

#if defined(__CUDACC__)
__global__
void gpuBoundaryXFacesKernel(
  std::size_t nx, std::size_t ny, std::size_t nz,
  std::size_t p_nx, std::size_t p_ny,
  real_t* RESTRICT u_new
) {
  const std::size_t j{blockIdx.x * blockDim.x + threadIdx.x + 1};
  const std::size_t k{blockIdx.y * blockDim.y + threadIdx.y + 1};

  if (j >= ny-1 || k >= nz-1) { return; }

  const std::size_t left_ghost{(0) + p_nx * ((j) + p_ny * (k))};
  const std::size_t left_inner{(1) + p_nx * ((j) + p_ny * (k))};
  const std::size_t right_ghost{(nx-1) + p_nx * ((j) + p_ny * (k))};
  const std::size_t right_inner{(nx-2) + p_nx * ((j) + p_ny * (k))};

  u_new[left_ghost] = u_new[left_inner];
  u_new[right_ghost] = u_new[right_inner];
}

__global__
void gpuBoundaryYFacesKernel(
  std::size_t nx, std::size_t ny, std::size_t nz,
  std::size_t p_nx, std::size_t p_ny,
  real_t* RESTRICT u_new
) {
  const std::size_t i{blockIdx.x * blockDim.x + threadIdx.x};
  const std::size_t k{blockIdx.y * blockDim.y + threadIdx.y + 1};

  if (i >= nx || k >= nz-1) { return; }

  const std::size_t front_ghost{(i) + p_nx * ((0) + p_ny * (k))};
  const std::size_t front_inner{(i) + p_nx * ((1) + p_ny * (k))};
  const std::size_t back_ghost{(i) + p_nx * ((ny-1) + p_ny * (k))};
  const std::size_t back_inner{(i) + p_nx * ((ny-2) + p_ny * (k))};

  u_new[front_ghost] = u_new[front_inner];
  u_new[back_ghost] = u_new[back_inner];
}

__global__
void gpuBoundaryZFacesKernel(
  std::size_t nx, std::size_t ny, std::size_t nz,
  std::size_t p_nx, std::size_t p_ny,
  real_t* RESTRICT u_new
) {
  const std::size_t i{blockIdx.x * blockDim.x + threadIdx.x};
  const std::size_t j{blockIdx.y * blockDim.y + threadIdx.y};

  if (i >= nx || j >= ny) { return; }

  const std::size_t bottom_ghost{(i) + p_nx * ((j) + p_ny * (0))};
  const std::size_t bottom_inner{(i) + p_nx * ((j) + p_ny * (1))};
  const std::size_t top_ghost{(i) + p_nx * ((j) + p_ny * (nz-1))};
  const std::size_t top_inner{(i) + p_nx * ((j) + p_ny * (nz-2))};

  u_new[bottom_ghost] = u_new[bottom_inner];
  u_new[top_ghost] = u_new[top_inner];
}

__global__
void gpuIntegrateGridKernel(
  std::size_t nx, std::size_t ny, std::size_t nz,
  std::size_t p_nx, std::size_t p_ny,
  real_t inv_dx_sq, real_t inv_dy_sq, real_t inv_dz_sq,
  real_t alpha_dt,
  const real_t* RESTRICT u_old,
  real_t* RESTRICT u_new
) {
  const std::size_t i{blockIdx.x * blockDim.x + threadIdx.x + 1};
  const std::size_t j{blockIdx.y * blockDim.y + threadIdx.y + 1};
  const std::size_t k{blockIdx.z * blockDim.z + threadIdx.z + 1};

  if (i >= nx-1 || j >= ny-1 || k >= nz-1) { return; }

  const std::size_t x_low{(i-1) + p_nx * (j + p_ny * k)};
  const std::size_t x_high{(i+1) + p_nx * (j + p_ny * k)};

  const std::size_t y_low{i + p_nx * ((j-1) + p_ny * k)};
  const std::size_t y_high{i + p_nx * ((j+1) + p_ny * k)};

  const std::size_t z_low{i + p_nx * (j + p_ny * (k-1))};
  const std::size_t z_high{i + p_nx * (j + p_ny * (k+1))};

  const std::size_t point{i + p_nx * (j + p_ny * k)};

  const real_t laplacian{
    (u_old[x_low] - static_cast<real_t>(2.0) * u_old[point] + u_old[x_high]) * inv_dx_sq +
    (u_old[y_low] - static_cast<real_t>(2.0) * u_old[point] + u_old[y_high]) * inv_dy_sq +
    (u_old[z_low] - static_cast<real_t>(2.0) * u_old[point] + u_old[z_high]) * inv_dz_sq
  };

  u_new[point] = u_old[point] + alpha_dt * laplacian;
}
#else
void cpuBoundaryConditionKernel(
  std::size_t nx, std::size_t ny, std::size_t nz,
  std::size_t p_nx, std::size_t p_ny,
  real_t* RESTRICT u_new
) {
  ASSUME_ALIGNED(u_new, SIMD_BYTES);

  // X-faces:
  for (std::size_t k = 1; k < nz-1; ++k) {

    #pragma omp simd
    for (std::size_t j = 1; j < ny-1; ++j) {
      const std::size_t left_ghost{(0) + p_nx * ((j) + p_ny * (k))};
      const std::size_t left_inner{(1) + p_nx * ((j) + p_ny * (k))};
      const std::size_t right_ghost{(nx-1) + p_nx * ((j) + p_ny * (k))};
      const std::size_t right_inner{(nx-2) + p_nx * ((j) + p_ny * (k))};

      u_new[left_ghost] = u_new[left_inner];
      u_new[right_ghost] = u_new[right_inner];
    }
  }

  // Y-faces:
  for (std::size_t k = 1; k < nz-1; ++k) {

    #pragma omp simd
    for (std::size_t i = 0; i < nx; ++i) {
      const std::size_t front_ghost{(i) + p_nx * ((0) + p_ny * (k))};
      const std::size_t front_inner{(i) + p_nx * ((1) + p_ny * (k))};
      const std::size_t back_ghost{(i) + p_nx * ((ny-1) + p_ny * (k))};
      const std::size_t back_inner{(i) + p_nx * ((ny-2) + p_ny * (k))};

      u_new[front_ghost] = u_new[front_inner];
      u_new[back_ghost] = u_new[back_inner];
    }
  }

  // Z-faces:
  for (std::size_t j = 0; j < ny; ++j) {

    #pragma omp simd
    for (std::size_t i = 0; i < nx; ++i) {
      const std::size_t bottom_ghost{(i) + p_nx * ((j) + p_ny * (0))};
      const std::size_t bottom_inner{(i) + p_nx * ((j) + p_ny * (1))};
      const std::size_t top_ghost{(i) + p_nx * ((j) + p_ny * (nz-1))};
      const std::size_t top_inner{(i) + p_nx * ((j) + p_ny * (nz-2))};

      u_new[bottom_ghost] = u_new[bottom_inner];
      u_new[top_ghost] = u_new[top_inner];
    }
  }
}

void cpuIntegrateGridKernel(
  std::size_t nx, std::size_t ny, std::size_t nz,
  std::size_t p_nx, std::size_t p_ny,
  real_t inv_dx_sq, real_t inv_dy_sq, real_t inv_dz_sq,
  real_t alpha_dt,
  const real_t* RESTRICT u_old,
  real_t* RESTRICT u_new
) {
  ASSUME_ALIGNED(u_new, SIMD_BYTES);
  ASSUME_ALIGNED(u_old, SIMD_BYTES);

  #pragma omp parallel for collapse(2)
  for (std::ptrdiff_t k = 1; k < static_cast<std::ptrdiff_t>(nz-1); ++k) {
    for (std::ptrdiff_t j = 1; j < static_cast<std::ptrdiff_t>(ny-1); ++j) {

      #pragma omp simd
      for (std::size_t i = 1; i < nx-1; ++i) {
        const std::size_t x_low{(i-1) + p_nx * (j + p_ny * k)};
        const std::size_t x_high{(i+1) + p_nx * (j + p_ny * k)};

        const std::size_t y_low{i + p_nx * ((j-1) + p_ny * k)};
        const std::size_t y_high{i + p_nx * ((j+1) + p_ny * k)};

        const std::size_t z_low{i + p_nx * (j + p_ny * (k-1))};
        const std::size_t z_high{i + p_nx * (j + p_ny * (k+1))};

        const std::size_t point{i + p_nx * (j + p_ny * k)};

        const real_t laplacian{
          (u_old[x_low] - static_cast<real_t>(2.0) * u_old[point] + u_old[x_high]) * inv_dx_sq +
          (u_old[y_low] - static_cast<real_t>(2.0) * u_old[point] + u_old[y_high]) * inv_dy_sq +
          (u_old[z_low] - static_cast<real_t>(2.0) * u_old[point] + u_old[z_high]) * inv_dz_sq
        };

        u_new[point] = u_old[point] + alpha_dt * laplacian;
      }
    }
  }
}
#endif

} // namespace

void ExplicitEuler::integrate(const Grid& old_grid, Grid& new_grid) {
#if defined(__CUDACC__)
  dim3 threads(8, 8, 8);
  dim3 blocks(
    (static_cast<unsigned>(old_grid.nx()-2) + threads.x - 1) / threads.x,
    (static_cast<unsigned>(old_grid.ny()-2) + threads.y - 1) / threads.y,
    (static_cast<unsigned>(old_grid.nz()-2) + threads.z - 1) / threads.z
  );

  gpuIntegrateGridKernel<<<blocks, threads>>>(
    old_grid.nx(), old_grid.ny(), old_grid.nz(),
    old_grid.p_nx(), old_grid.p_ny(),
    old_grid.inv_dx_sq(), old_grid.inv_dy_sq(), old_grid.inv_dz_sq(),
    alpha() * dt(),
    old_grid.field(),
    new_grid.field()
  );
  CUDA_CHECK(cudaGetLastError());

  dim3 face_threads(16, 16);

  dim3 x_blocks(
    (static_cast<unsigned>(old_grid.ny()-2) + face_threads.x - 1) / face_threads.x,
    (static_cast<unsigned>(old_grid.nz()-2) + face_threads.y - 1) / face_threads.y
  );
  gpuBoundaryXFacesKernel<<<x_blocks, face_threads>>>(
    old_grid.nx(), old_grid.ny(), old_grid.nz(),
    old_grid.p_nx(), old_grid.p_ny(),
    new_grid.field()
  );
  CUDA_CHECK(cudaGetLastError());

  dim3 y_blocks(
    (static_cast<unsigned>(old_grid.nx()) + face_threads.x - 1) / face_threads.x,
    (static_cast<unsigned>(old_grid.nz()-2) + face_threads.y - 1) / face_threads.y
  );
  gpuBoundaryYFacesKernel<<<y_blocks, face_threads>>>(
    old_grid.nx(), old_grid.ny(), old_grid.nz(),
    old_grid.p_nx(), old_grid.p_ny(),
    new_grid.field()
  );
  CUDA_CHECK(cudaGetLastError());

  dim3 z_blocks(
    (static_cast<unsigned>(old_grid.nx()) + face_threads.x - 1) / face_threads.x,
    (static_cast<unsigned>(old_grid.ny()) + face_threads.y - 1) / face_threads.y
  );
  gpuBoundaryZFacesKernel<<<z_blocks, face_threads>>>(
    old_grid.nx(), old_grid.ny(), old_grid.nz(),
    old_grid.p_nx(), old_grid.p_ny(),
    new_grid.field()
  );
  CUDA_CHECK(cudaGetLastError());
#else
  cpuIntegrateGridKernel(
    old_grid.nx(), old_grid.ny(), old_grid.nz(),
    old_grid.p_nx(), old_grid.p_ny(),
    old_grid.inv_dx_sq(), old_grid.inv_dy_sq(), old_grid.inv_dz_sq(),
    alpha() * dt(),
    old_grid.field(),
    new_grid.field()
  );

  cpuBoundaryConditionKernel(
    old_grid.nx(), old_grid.ny(), old_grid.nz(),
    old_grid.p_nx(), old_grid.p_ny(),
    new_grid.field()
  );
#endif
}