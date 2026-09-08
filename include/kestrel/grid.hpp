#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>

#include <kestrel/precision.hpp>

namespace kestrel {

template <std::size_t Dimensions>
class grid {
  static_assert(Dimensions > 0, "A grid must have at least one dimension");

public:
  static constexpr auto dimensions{Dimensions};

  using extent = std::array<std::size_t, Dimensions>;
  using spacing = std::array<fp_t, Dimensions>;

private:
  extent cells_;
  spacing cell_spacing_;

public:
  constexpr explicit grid(extent cells, spacing cell_spacing) noexcept
    : cells_{cells}
    , cell_spacing_{cell_spacing}
  { }

  [[nodiscard]]
  constexpr auto total_cells() const noexcept -> std::size_t {
    return std::ranges::fold_left(cells_, 1uz, std::multiplies<>());
  }

  [[nodiscard]]
  constexpr auto cells() const noexcept -> extent {
    return cells_;
  }

  [[nodiscard]]
  constexpr auto cell_spacing() const noexcept -> spacing {
    return cell_spacing_;
  }
};

} // namespace kestrel
