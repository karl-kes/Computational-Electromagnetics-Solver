#pragma once

#include <cstdint>

namespace constant {
    static constexpr std::size_t Nx{ 12 }, Ny{ 12 }, Nz{ 12 };
    static constexpr double dx{ 5.0 }, dy{ 5.0 }, dz{ 5.0 };
    static constexpr double eps{ 1.0 }, mu{ 1.0 };
    static constexpr double cfl_factor{ 0.025 };
    static constexpr std::size_t total_time{ 1000 };

    static constexpr double inject{ 10.0 };
    static constexpr double amp_one{ 10.0 };
    static constexpr double amp_two{ 10.0 };
    static constexpr double freq_one{ 1.0 };
    static constexpr double freq_two{ 1.0 };

    static constexpr char B_field{ 'b' };
    static constexpr char E_field{ 'e' };
}