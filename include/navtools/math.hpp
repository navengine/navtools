#ifndef NAVENGINE_INCLUDE_NAVTOOLS_MATH_HPP
#define NAVENGINE_INCLUDE_NAVTOOLS_MATH_HPP

#include <cmath>

namespace navtools {

template<typename RealType>
constexpr RealType circular_fmod(RealType& x, const RealType y)
{
  assert(y > static_cast<RealType>(0));
  RealType mult = std::floor(x / y);
  x -= mult * y;
  return mult;
}

} // namespace navtools
#endif
