#include "grid.hpp"

Grid::Grid( Simulation_Config const &config )
: Nx_{ config.Nx + 1 }
, Ny_{ config.Ny + 1 }
, Nz_{ config.Nz + 1 }
, dx_{ config.dx }
, dy_{ config.dy }
, dz_{ config.dz }
, eps_{ config.eps }
, mu_{ config.mu }
, c_{ config.c }
, c_sq_{ config.c * config.c }
, dt_{ config.dt } {
    auto allocate = [size = config.size]() {
        return std::make_unique<double[]>( size );
    };

    Ex_ = allocate();
    Ey_ = allocate();
    Ez_ = allocate();

    Bx_ = allocate();
    By_ = allocate();
    Bz_ = allocate();

    Jx_ = allocate();
    Jy_ = allocate();
    Jz_ = allocate();
}