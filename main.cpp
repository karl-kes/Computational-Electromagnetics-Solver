#include "Classes/Config/config.hpp"
#include "Classes/Grid/grid.hpp"
#include "Classes/Source/source.hpp"
#include "Classes/Write_Output/output.hpp"
#include "Simulation.hpp"
#include "Validation.hpp"

/* 
    To compile and run.

    For No Parallel (NOT RECOMMENDED):
    g++ -std=c++17 main.cpp Classes/Grid/*.cpp Classes/Source/*.cpp Classes/Write_Output/*.cpp -o main.exe
    ./main.exe
    ./render.py

    For Parallel (RECOMMENDED):
    g++ -std=c++17 -O3 -march=native -fopenmp main.cpp Classes/Grid/*.cpp Classes/Source/*.cpp Classes/Write_Output/*.cpp -o main.exe
    ./main.exe
    ./render.py
*/

int main() {
    Simulation_Config config{};
    Simulation sim{ config };

    // Validation:
    if ( config.run_validation ) {
        Plane_Wave_Test test{ config };
        Validation_Result result{ test.run() };
        test.print_report( result );
    }

    // Simulation:
    sim.initialize();
    sim.run();

    return 0;
}