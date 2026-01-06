#include "Classes/Grid/grid.hpp"
#include "constant.hpp"

/* 
    To compile and run.

    For No Parallel (NOT RECOMMENDED):
    g++ -std=c++17 main.cpp Classes/Grid/grid.cpp -o main.exe
    ./main.exe
    ./render.py

    For Parallel (RECOMMENDED):
    g++ -std=c++17 main.cpp Classes/Grid/grid.cpp -o main.exe -fopenmp
    ./main.exe
    ./render.py
*/

int main() {

    // Construct the grid and create its directories
    Grid grid{ constant::Nx+1, constant::Ny+1, constant::Nz+1 };
    grid.create_directories();

    // Inject the original source
    grid.inject_source( constant::Nx/2, constant::Ny/2, constant::Nz/2, constant::inject );

    // Calculate initial values
    double initial_energy{ grid.total_energy() };
    double initial_div_B{ grid.div_B() };

    // Initialize max as initial
    double max_energy_drift{ initial_energy };
    double max_div_B{ initial_div_B };

    auto start{ std::chrono::high_resolution_clock::now() };
    for ( int t{}; t <= constant::elapsed_time; ++t ) {
        // Step simulation forward
        grid.step();

        // Max between current and previous max
        max_energy_drift = std::max( grid.total_energy(), max_energy_drift );
        max_div_B = std::max( grid.div_B(), max_div_B );

        // Output simulation data
        if ( t % 10 == 0 ) {
            std::string t_sec{ std::to_string( t ) };
            grid.vector_volume( "output/E/E" + t_sec + ".bin", constant::E_field );
            grid.vector_volume( "output/B/B" + t_sec + ".bin", constant::B_field );
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    // Calculate max energy drift from initial system energy
    max_energy_drift = ( 100.0 * ( max_energy_drift - initial_energy ) / initial_energy );

    // Output final metrics
    std::cout << std::endl;
    std::cout << "Duration of Simulation: " << duration.count() << " ms\n" << std::endl;
    std::cout << "Physical Time Simulated: " << constant::elapsed_time * grid.dt() << " s" << std::endl;
    std::cout << "Max Energy Drift: " << max_energy_drift << "%" << std::endl;
    std::cout << "Max Magnetic Divergence: " << max_div_B << std::endl << std::endl;

    return 0;
}