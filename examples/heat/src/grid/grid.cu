#include "grid.cuh"

Grid::Grid(const Config& config)
: nx_{config.nx}, ny_{config.ny}, nz_{config.nz}
, p_nx_{AlignedSoA<real_t>::round_up(config.nx)}
, p_ny_{AlignedSoA<real_t>::round_up(config.ny)}
, p_nz_{AlignedSoA<real_t>::round_up(config.nz)}
, dx_{config.dx}, dy_{config.dy}, dz_{config.dz}
, inv_dx_sq_{static_cast<real_t>(1.0)/(dx_*dx_)}
, inv_dy_sq_{static_cast<real_t>(1.0)/(dy_*dy_)}
, inv_dz_sq_{static_cast<real_t>(1.0)/(dz_*dz_)}
, data_{total_size(), NUM_SUB_ARR}
{ }

Grid::~Grid() = default;

void Grid::copy_to_host(real_t* dst) const {
  const std::size_t total{total_size()};

#if defined(__CUDACC__)
  CUDA_CHECK(cudaMemcpy(dst, field(), total * sizeof(real_t), cudaMemcpyDeviceToHost));
#else
  std::copy_n(field(), total, dst);
#endif
}