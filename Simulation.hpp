#pragma once

#include "Classes/Config/config.hpp"
#include "Classes/Grid/grid.hpp"
#include "Classes/Source/source.hpp"
#include "Classes/Write_Output/output.hpp"

class Output;
class Simulation_Config;

class Simulation {
public:
    Simulation_Config config;
    Grid grid;
    Output output;

    explicit Simulation( Simulation_Config const &new_config )
    : config{ new_config }
    , grid{ config }
    , output{ "output" }
    { }

    void initialize() {
        output.initialize();

        // Add sources:
        grid.add_source( std::make_unique<Straight_Wire_X>(
            10.0,                           // amplitude
            1.0,                            // frequency
            config.Ny / 2,                  // y position
            config.Nz / 2,                  // z position
            config.Nx / 4,                  // x start
            3 * config.Nx / 4               // x end
        ) );

        // grid.add_source( std::make_unique<Point_Source>(
        //     10.0,
        //     config.Nx / 2,
        //     config.Ny / 2,
        //     config.Nz / 2
        // ) );

        // grid.add_source( std::make_unique<Gaussian_Pulse>(
        //     10.0,                   // Amplitude
        //     30 * 3 * grid.dt(),     // Center Time
        //     30 * grid.dt(),         // Width
        //     config.Nx / 2,          // x
        //     config.Ny / 2,          // y
        //     config.Nz / 2           // z
        // ) );
    }

    void run() {
        std::cout << "<-----Maxwell Simulation----->" << std::endl;

        // Run @ t = 0:
        grid.apply_sources();
        grid.step( config, output, 0 );

        // Run simulation and start timer:
        std::size_t const output_interval{ config.output_interval() };
        auto const start_time{ std::chrono::high_resolution_clock::now() };

        // Simulation Loop:
        for ( std::size_t curr_time{1}; curr_time <= config.total_time; ++curr_time ) {
            grid.apply_sources( curr_time );
            grid.step( config, output, curr_time );

            if ( ( curr_time % config.output_interval() ) == 0 ) {
                output.write_field( grid, curr_time );
                grid.print_progress( curr_time, config.total_time );
            }
        }

        // End Timer:
        auto const end_time{ std::chrono::high_resolution_clock::now() };
        auto const duration{ std::chrono::duration_cast<std::chrono::milliseconds>( end_time - start_time ) };

        // Report Results:
        std::cout << "\n\nResults:\n";
        std::cout << "--------\n";
        std::cout << "Duration: " << duration.count() << " ms\n";
        std::cout << "Physical time: " << config.total_time * grid.dt() << " s\n";
    }
};