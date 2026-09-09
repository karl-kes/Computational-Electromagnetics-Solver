#pragma once

#include <chrono>

namespace kestrel {

inline auto timer(auto&& callable) -> std::uint64_t {
  const auto begin{std::chrono::steady_clock::now()};
  callable();
  const auto end{std::chrono::steady_clock::now()};

  const auto duration{
  std::chrono::duration_cast<std::chrono::nanoseconds>(
      end - begin
    )
  };

  return static_cast<std::uint64_t>(duration.count());
}

} // namespace kestrel
