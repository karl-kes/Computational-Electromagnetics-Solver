#include "grid.hpp"

std::size_t Grid::idx( std::size_t const x, std::size_t const y, std::size_t const z ) const {
    return x + Nx() * ( y + Ny() * z );
}

double Grid::curl_x(
    double const Y_0, double const Y_1,
    double const Z_0, double const Z_1 ) const {

    double dz_dy{ ( Z_1 - Z_0 ) / dy() };
    double dy_dz{ ( Y_1 - Y_0 ) / dz() };

    return dy_dz- dz_dy;
}

double Grid::curl_y(
    double const X_0, double const X_1,
    double const Z_0, double const Z_1 ) const {

    double dz_dx{ ( Z_1 - Z_0 ) / dx() };
    double dx_dz{ ( X_1 - X_0 ) / dz() };

    return dz_dx - dx_dz;
}

double Grid::curl_z(
    double const Y_0, double const Y_1,
    double const X_0, double const X_1 ) const {

    double dy_dx{ ( Y_1 - Y_0 ) / dx() };
    double dx_dy{ ( X_1 - X_0 ) / dy() };

    return dx_dy - dy_dx;
}

double Grid::total_energy() const {
    double energy{};
    double dV{ dx() * dy() * dz() };

    #pragma omp parallel for collapse( 2 ) reduction( +:energy )
    for ( std::size_t z = 0; z < Nz(); ++z ) {
        for ( std::size_t y = 0; y < Ny(); ++y ) {
            #pragma omp simd
            for ( std::size_t x = 0; x < Nx(); ++x ) {
                    double E_sq{ Ex(x,y,z)*Ex(x,y,z) + Ey(x,y,z)*Ey(x,y,z) + Ez(x,y,z)*Ez(x,y,z) };
                    double B_sq{ Bx(x,y,z)*Bx(x,y,z) + By(x,y,z)*By(x,y,z) + Bz(x,y,z)*Bz(x,y,z) };

                    energy += 0.5 * ( eps() * E_sq + B_sq / mu() );
            }
        }
    }
    return energy * dV;
}