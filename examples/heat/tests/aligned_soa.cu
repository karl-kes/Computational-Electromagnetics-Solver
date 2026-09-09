#include "cuda_test_guard.cuh"
#include "utilities/aligned_soa.cuh"

#include <algorithm>
#include <vector>

namespace {

std::vector<real_t> copy_array(const AlignedSoA<real_t>& soa, std::size_t array_index) {
  std::vector<real_t> values(soa.stride());
#if defined(__CUDACC__)
  CUDA_CHECK(cudaMemcpy(values.data(), soa[array_index], values.size() * sizeof(real_t), cudaMemcpyDeviceToHost));
#else
  std::copy_n(soa[array_index], values.size(), values.data());
#endif
  return values;
}

} // namespace

int main() {
  HEAT_SOLVER_SKIP_CUDA_TEST_IF_UNAVAILABLE();

  constexpr std::size_t elems_per_align{SIMD_BYTES / sizeof(real_t)};

  // round_up pads up to the next multiple of the SIMD alignment
  if (AlignedSoA<real_t>::round_up(0) != 0) { return 1; }
  if (AlignedSoA<real_t>::round_up(elems_per_align) != elems_per_align) { return 1; }
  if (AlignedSoA<real_t>::round_up(1) != elems_per_align) { return 1; }
  if (AlignedSoA<real_t>::round_up(elems_per_align + 1) != 2 * elems_per_align) { return 1; }

  // allocated storage is correctly strided and zero-initialized
  AlignedSoA<real_t> soa{5, 2};
  if (soa.stride() != AlignedSoA<real_t>::round_up(5)) { return 1; }

  for (std::size_t arr{}; arr < 2; ++arr) {
    const std::vector<real_t> values{copy_array(soa, arr)};
    for (const real_t value : values) {
      if (value != real_t{}) { return 1; }
    }
  }

  return 0;
}
