#pragma once

#include <chrono>

class Timer {
  using Clock = std::chrono::high_resolution_clock;
  using TimePoint = Clock::time_point;

  TimePoint start_;

public:
  Timer() : start_{Clock::now()} {}

  void reset() { start_ = Clock::now(); }

  double elapsed_ms() const {
    return std::chrono::duration<double, std::milli>(Clock::now() - start_).count();
  }

  double elapsed_s() const {
    return std::chrono::duration<double>(Clock::now() - start_).count();
  }
};
