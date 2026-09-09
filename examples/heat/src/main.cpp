#include "simulation/simulation.hpp"
#include "utilities/timer.hpp"

#include <iostream>

int main(int argc, char** argv) {
  const Config cfg{Config::parse(argc, argv)};
  Simulation sim{cfg};

  Timer timer{};
  sim.run();
  std::cout << timer.elapsed_ms() << "ms\n";
  std::cout.flush();

  return 0;
}