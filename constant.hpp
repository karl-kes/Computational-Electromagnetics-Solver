#pragma once

#include <cstdint>

namespace constant {
    static constexpr std::size_t Nx{ 5 }, Ny{ 5 }, Nz{ 5 };
    static constexpr double cfl_factor{ 0.02 };

    static constexpr double inject{ 10.0 };
    static constexpr int elapsed_time{ 1000 };

    static constexpr char B_field{ 'B' };
    static constexpr char E_field{ 'E' };
}