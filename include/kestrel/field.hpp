#pragma once

#include <cstddef>

#include <kestrel/concepts.hpp>
#include <xpu/soa.hpp>

namespace kestrel {

template <arithmetic T, std::size_t Components>
struct field_view {
  static_assert(Components > 0, "A field must have at least one component.");

  xpu::soa_view<T, Components> data;
};

template <arithmetic T, std::size_t Components>
class field {
  static_assert(Components > 0);

private:
  xpu::soa<T, Components> storage_;

public:
  explicit field(std::size_t count)
    : storage_{count}
  { }

  using const_view_t = field_view<const T, Components>;
  using view_t = field_view<T, Components>;

  [[nodiscard]]
  auto view() const noexcept -> const_view_t {
    return {storage_.view()};
  }

  [[nodiscard]]
  auto view() noexcept -> view_t {
    return {storage_.view()};
  }
};

template <arithmetic T>
using scalar_field = field<T, 1uz>;

template <arithmetic T>
using vector_field = field<T, 3uz>;

} // namespace kestrel
