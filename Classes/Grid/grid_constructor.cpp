#include "grid.hpp"

Grid::Grid( Simulation_Config const &config ):
Nx_{ config.Nx + 1 }, Ny_{ config.Ny + 1 }, Nz_{ config.Nz + 1 },
dx_{ config.dx }, dy_{ config.dy }, dz_{ config.dz },
eps_{ config.eps }, mu_{ config.mu },
c_{ config.c }, c_sq_{ config.c * config.c }, dt_{ config.dt } {
    Ex_ = std::make_unique<double[]>( config.size );
    Ey_ = std::make_unique<double[]>( config.size );
    Ez_ = std::make_unique<double[]>( config.size );

    Bx_ = std::make_unique<double[]>( config.size );
    By_ = std::make_unique<double[]>( config.size );
    Bz_ = std::make_unique<double[]>( config.size );

    Jx_ = std::make_unique<double[]>( config.size );
    Jy_ = std::make_unique<double[]>( config.size );
    Jz_ = std::make_unique<double[]>( config.size );
}