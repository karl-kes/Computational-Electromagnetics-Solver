#include "config.hpp"

#include <charconv>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string_view>

namespace {

void print_usage() {
  std::cout <<
    "Usage: kestrel-heat [options]\n"
    "  --nx <n>                   grid size in x (default 64)\n"
    "  --ny <n>                   grid size in y (default 64)\n"
    "  --nz <n>                   grid size in z (default 64)\n"
    "  --steps <n>                total integration steps (default 1000)\n"
    "  --output-interval <n>      VTK output interval, 0 = disabled (default 0)\n"
    "  --ic <gaussian|cosine>     initial condition (default gaussian)\n"
    "  --alpha <f>                thermal diffusivity (default 1.0)\n"
    "  --dx <f> --dy <f> --dz <f> grid spacing (default 1.0)\n"
    "  --help                     show this message\n";
}

template <typename T>
T parse_number(std::string_view flag, std::string_view raw) {
  T value{};
  const auto result{std::from_chars(raw.data(), raw.data() + raw.size(), value)};
  if (result.ec != std::errc{} || result.ptr != raw.data() + raw.size()) {
    std::cerr << "invalid value for " << flag << ": " << raw << '\n';
    std::exit(1);
  }
  return value;
}

void require(bool condition, const char* message) {
  if (!condition) {
    std::cerr << message << '\n';
    std::exit(1);
  }
}

} // namespace

Config Config::parse(int argc, char** argv) {
  Config cfg{};

  for (int i{1}; i < argc; ++i) {
    const std::string_view flag{argv[i]};

    if (flag == "--help" || flag == "-h") {
      print_usage();
      std::exit(0);
    }

    if (i + 1 >= argc) {
      std::cerr << "missing value for " << flag << '\n';
      std::exit(1);
    }

    const std::string_view value{argv[++i]};

    if (flag == "--nx") { cfg.nx = parse_number<std::size_t>(flag, value); }
    else if (flag == "--ny") { cfg.ny = parse_number<std::size_t>(flag, value); }
    else if (flag == "--nz") { cfg.nz = parse_number<std::size_t>(flag, value); }
    else if (flag == "--steps") { cfg.total_steps = parse_number<std::size_t>(flag, value); }
    else if (flag == "--output-interval") { cfg.output_interval = parse_number<std::size_t>(flag, value); }
    else if (flag == "--ic") {
      if (value == "gaussian") { cfg.ic = InitCondition::Gaussian; }
      else if (value == "cosine") { cfg.ic = InitCondition::NeumannCosine; }
      else {
        std::cerr << "invalid value for --ic: " << value << '\n';
        std::exit(1);
      }
    }
    else if (flag == "--alpha") { cfg.alpha = parse_number<real_t>(flag, value); }
    else if (flag == "--dx") { cfg.dx = parse_number<real_t>(flag, value); }
    else if (flag == "--dy") { cfg.dy = parse_number<real_t>(flag, value); }
    else if (flag == "--dz") { cfg.dz = parse_number<real_t>(flag, value); }
    else {
      std::cerr << "unknown option: " << flag << '\n';
      print_usage();
      std::exit(1);
    }
  }

  require(cfg.nx >= 3, "--nx must be at least 3");
  require(cfg.ny >= 3, "--ny must be at least 3");
  require(cfg.nz >= 3, "--nz must be at least 3");
  require(std::isfinite(cfg.alpha) && cfg.alpha > static_cast<real_t>(0), "--alpha must be positive");
  require(std::isfinite(cfg.dx) && cfg.dx > static_cast<real_t>(0), "--dx must be positive");
  require(std::isfinite(cfg.dy) && cfg.dy > static_cast<real_t>(0), "--dy must be positive");
  require(std::isfinite(cfg.dz) && cfg.dz > static_cast<real_t>(0), "--dz must be positive");

  cfg.dt = stable_dt(cfg.alpha, cfg.dx, cfg.dy, cfg.dz);

  return cfg;
}
