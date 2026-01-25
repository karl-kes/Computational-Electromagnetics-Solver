#include "grid.hpp"

void Grid::add_source( std::unique_ptr<Source> source ) {
    sources_.push_back( std::move( source ) );
}

void Grid::apply_sources( double time_step ) {
    for ( auto const &source : sources_ ) {
        source->apply( *this, time_step );
    }
}

void Grid::update_B() {
    // ∂B/∂t = -curl( E )
    #pragma omp parallel for collapse( 2 )
    // Start at 0; stagger for Yee-Cell grid.
    for ( std::size_t z = 0; z < Nz() - 1; ++z ) {
        for ( std::size_t y = 0; y < Ny() - 1; ++y ) {
            #pragma omp simd
            for ( std::size_t x = 0; x < Nx() - 1; ++x ) {
                // Take curl of components and apply B -= ∂B:
                // ∂B_x = ∂t * curl_x( E )
                Bx(x,y,z) -= dt() * curl_x( Ey(x,y,z), Ey(x,y,z+1),
                                            Ez(x,y,z), Ez(x,y+1,z) );

                // ∂B_y = ∂t * curl_y( E )
                By(x,y,z) -= dt() * curl_y( Ex(x,y,z), Ex(x,y,z+1),
                                            Ez(x,y,z), Ez(x+1,y,z) );

                // ∂B_z = ∂t * curl_z( E )
                Bz(x,y,z) -= dt() * curl_z( Ey(x,y,z), Ey(x+1,y,z),
                                            Ex(x,y,z), Ex(x,y+1,z) );
            }
        }
    }
}

void Grid::update_E() {
    // ∂E/∂t = c*c * curl(B)
    #pragma omp parallel for collapse( 2 )
    // Start at 1; stagger for Yee-Cell grid.
    for ( std::size_t z = 1; z < Nz(); ++z ) {
        for ( std::size_t y = 1; y < Ny(); ++y ) {
            #pragma omp simd
            for ( std::size_t x = 1; x < Nx(); ++x ) {
                // Curl of components and apply E += ∂E:
                // ∂E_x = ∂t * c*c * (∂E_z/∂y - ∂E_y/∂z)
                Ex(x,y,z) += dt() * (
                                      c_sq() *
                                      curl_x( By(x,y,z-1), By(x,y,z),
                                              Bz(x,y-1,z), Bz(x,y,z) ) -
                                      Jx(x,y,z) / eps()
                                    );

                // ∂E_y = ∂t * c*c * (∂Ex/∂z - ∂Ez/∂x)
                Ey(x,y,z) += dt() * (
                                      c_sq() *
                                      curl_y( Bx(x,y,z-1), Bx(x,y,z),
                                              Bz(x-1,y,z), Bz(x,y,z) ) -
                                      Jy(x,y,z) / eps()
                                    );

                // ∂E_z = ∂t * c*c * (∂Ex/∂y - ∂Ey/∂x)
                Ez(x,y,z) += dt() * (
                                      c_sq() *
                                      curl_z( By(x-1,y,z), By(x,y,z),
                                              Bx(x,y-1,z), Bx(x,y,z) ) -
                                      Jz(x,y,z) / eps()
                                    );
            }
        }
    }
}

void Grid::step() {
    update_B();
    update_E();
}

double Grid::field( char field, char component, 
                    std::size_t x, std::size_t y, std::size_t z ) const {
    if ( field == 'e' ) {
        switch ( component ) {
            case 'x': return Ex(x,y,z);
            case 'y': return Ey(x,y,z);
            case 'z': return Ez(x,y,z);
        }
    } else if ( field == 'b' ) {
        switch ( component ) {
            case 'x': return Bx(x,y,z);
            case 'y': return By(x,y,z);
            case 'z': return Bz(x,y,z);
        }
    }
    throw std::invalid_argument{ "Invalid field or component specifier" };
}

double &Grid::field( char field, char component,
                     std::size_t x, std::size_t y, std::size_t z ) {
    if ( field == 'e' ) {
        switch ( component ) {
            case 'x': return Ex(x,y,z);
            case 'y': return Ey(x,y,z);
            case 'z': return Ez(x,y,z);
        }
    } else if ( field == 'b' ) {
        switch ( component ) {
            case 'x': return Bx(x,y,z);
            case 'y': return By(x,y,z);
            case 'z': return Bz(x,y,z);
        }
    }
    throw std::invalid_argument{ "Invalid field or component specifier" };
}

double Grid::field_magnitude(char field, std::size_t x, std::size_t y, std::size_t z) const {
    double Fx{ this->field( field, 'x', x, y, z ) };
    double Fy{ this->field( field, 'y', x, y, z ) };
    double Fz{ this->field( field, 'z', x, y, z ) };

    return std::sqrt( Fx*Fx + Fy*Fy + Fz*Fz );
}