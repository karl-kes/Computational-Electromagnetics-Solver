#pragma once

#include "../grid/grid.cuh"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <vector>

namespace vtk {

namespace detail {
  inline uint32_t bswap32(uint32_t v) {
    return ((v & 0xFF000000u) >> 24) |
           ((v & 0x00FF0000u) >>  8) |
           ((v & 0x0000FF00u) <<  8) |
           ((v & 0x000000FFu) << 24);
  }

  inline uint64_t bswap64(uint64_t v) {
    return ((v & 0xFF00000000000000ull) >> 56) |
           ((v & 0x00FF000000000000ull) >> 40) |
           ((v & 0x0000FF0000000000ull) >> 24) |
           ((v & 0x000000FF00000000ull) >>  8) |
           ((v & 0x00000000FF000000ull) <<  8) |
           ((v & 0x0000000000FF0000ull) << 24) |
           ((v & 0x000000000000FF00ull) << 40) |
           ((v & 0x00000000000000FFull) << 56);
  }

  inline const std::filesystem::path& output_dir() {
    static const std::filesystem::path dir{"out"};

    static const bool created{
      (std::filesystem::create_directories(dir), true)
    };
    static_cast<void>(created);

    return dir;
  }
}

inline void write(const Grid& grid, std::size_t step) {
  char name[32];
  
  std::snprintf(name, sizeof(name), "step_%04zu.vtk", step);
  const auto filename{detail::output_dir() / name};

  std::ofstream out(filename, std::ios::binary);
  if (!out) { throw std::runtime_error("failed to open vtk output file"); }

  const std::size_t nx{grid.nx()};
  const std::size_t ny{grid.ny()};
  const std::size_t nz{grid.nz()};

  out << "# vtk DataFile Version 3.0\n"
      << "Heat solver step " << step << "\n"
      << "BINARY\n"
      << "DATASET STRUCTURED_POINTS\n"
      << "DIMENSIONS " << nx << ' ' << ny << ' ' << nz << '\n'
      << "ORIGIN 0 0 0\n"
      << "SPACING " << grid.dx() << ' ' << grid.dy() << ' ' << grid.dz() << '\n'
      << "POINT_DATA " << nx * ny * nz << '\n'
      << "SCALARS temperature " << (sizeof(real_t) == sizeof(float) ? "float" : "double") << " 1\n"
      << "LOOKUP_TABLE default\n";

  const std::size_t padded_total{grid.total_size()};
  std::vector<real_t> host_field(padded_total);
  grid.copy_to_host(host_field.data());
  const real_t* u{host_field.data()};

  for (std::size_t k{}; k < nz; ++k) {
    for (std::size_t j{}; j < ny; ++j) {
      if constexpr (sizeof(real_t) == sizeof(float)) {
        std::vector<uint32_t> row(nx);
        for (std::size_t i{}; i < nx; ++i) {
          const float val{static_cast<float>(u[grid.idx(i, j, k)])};
          uint32_t bits;
          std::memcpy(&bits, &val, sizeof(bits));
          row[i] = detail::bswap32(bits);
        }

        out.write(
          reinterpret_cast<const char*>(row.data()),
          static_cast<std::streamsize>(nx * sizeof(uint32_t))
        );
      } else {
        std::vector<uint64_t> row(nx);
        for (std::size_t i{}; i < nx; ++i) {
          const double val{static_cast<double>(u[grid.idx(i, j, k)])};
          uint64_t bits;
          std::memcpy(&bits, &val, sizeof(bits));
          row[i] = detail::bswap64(bits);
        }

        out.write(
          reinterpret_cast<const char*>(row.data()),
          static_cast<std::streamsize>(nx * sizeof(uint64_t))
        );
      }
    }
  }
}

} // namespace vtk