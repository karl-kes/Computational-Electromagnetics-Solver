#pragma once

#include <cstdint>

namespace constant {
    static constexpr std::size_t Nx{ 8 }, Ny{ 8 }, Nz{ 8 };
    static constexpr double cfl_factor{ 0.02 };
    static constexpr int elapsed_time{ 1000 };

    static constexpr double inject{ 10.0 };
    static constexpr double amp_one{ 10.0 };
    static constexpr double amp_two{ 10.0 };
    static constexpr double freq_one{ 1.0 };
    static constexpr double freq_two{ 1.0 };

    static constexpr char B_field{ 'B' };
    static constexpr char E_field{ 'E' };
}