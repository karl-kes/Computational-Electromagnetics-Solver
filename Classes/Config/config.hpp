#pragma once

#include <memory>
#include <fstream>
#include <iostream>
#include <cctype>
#include <string>
#include <omp.h>
#include <vector>
#include <filesystem>
#include <cstdint>
#include <cmath>
#include <stdexcept>

class Simulation_Config {
public:
    // Grid Dimensions:
    std::size_t Nx{ 100 };
    std::size_t Ny{ 100 };
    std::size_t Nz{ 100 };
    std::size_t size{ ( Nx + 1 ) * ( Ny + 1 ) * ( Nz + 1 ) };

    // Spatial Step Size:
    double dx{ 5.0 };
    double dy{ 5.0 };
    double dz{ 5.0 };

    // Material Properties:
    double mu{ 1.0 };
    double eps{ 1.0 };
    double c{ 1.0 / ( std::sqrt ( mu * eps ) ) };

    // Time Stepping:
    double cfl_factor{ 0.5 };
    double total_time{ 1000.0 };
    double dt{ cfl_factor / ( c * std::sqrt( 1.0 / (dx*dx) + 1.0 / (dy*dy) + 1.0 / (dz*dz) ) ) };

    // Validation:
    bool run_validation{ true };

    std::size_t output_interval() const { return static_cast<std::size_t>( total_time / 50) ; }
};

enum class Field {
    ELECTRIC,
    MAGNETIC
};

enum class Component {
    X,
    Y,
    Z
};

namespace config {
    static constexpr double PI{ 3.14159265358979323846264338327950288419716939937510 };
}