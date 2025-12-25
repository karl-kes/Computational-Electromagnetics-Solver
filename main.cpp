#include <chrono>
#include <filesystem>

#include "Classes/grid.hpp"

int main() {
    /* 
        To compile and run.

        For No Parallel (NOT RECOMMENDED):
        g++ -std=c++17 main.cpp Classes/grid.cpp -o main.exe
        ./main.exe
        ./render.py

        For Parallel (RECOMMENDED):
        g++ -std=c++17 main.cpp Classes/grid.cpp -o main.exe -fopenmp
        ./main.exe
        ./render.py
    */

    // Clear previous and create new output folder.
    std::filesystem::remove_all("output");
    std::filesystem::create_directories("output");

    static constexpr int elapsed_time{ 5000 };
    static constexpr std::size_t Nx{ 20 }, Ny{ 20 }, Nz{ 20 };
    static constexpr double inject{ 10.0 };

    Grid grid{ Nx, Ny, Nz };

    grid.inject_source( 3*Nx/4, 3*Ny/4,3*Nz/4, inject );
    grid.inject_source( Nx/4, Ny/4, Nz/4, inject );

    double initial_energy{ grid.total_energy() };
    double max_energy_drift{ initial_energy };

    auto start = std::chrono::high_resolution_clock::now();
    for ( int t{}; t <= elapsed_time; ++t ) {
        grid.step();

        // Max between current and previous max.
        max_energy_drift = std::max( grid.total_energy(), max_energy_drift );

        if ( t % 10 == 0 ) {
            std::string t_sec{ std::to_string( t ) };
            grid.vector_volume( "output/E" + t_sec + ".csv", 'E' );
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    double drift{ 100 * ( max_energy_drift - initial_energy ) / initial_energy };

    std::cout << "Physical Time Simulated: " << elapsed_time * grid.dt() << " s" << std::endl;
    std::cout << "Duration of Simulation: " << duration.count() << " ms" << std::endl;
    std::cout << "Max Energy Drift: " << drift << "%" << std::endl;

    return 0;
}