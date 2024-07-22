#ifndef NAVENGINE_INCLUDE_NAVTOOLS_CONSTANTS_HPP
#define NAVENGINE_INCLUDE_NAVTOOLS_CONSTANTS_HPP

#include <numbers>
#include <complex>

// FP is floating point
#define NAV_DEFINE_FP_CONSTANT(name,value) \
  template<typename Scalar = double> \
  inline constexpr Scalar name = static_cast<Scalar>(value)

namespace navtools {

template<typename T>
inline constexpr T TwoPi = static_cast<T>(2) * std::numbers::pi_v<T>;

template<typename T>
inline constexpr std::complex<T> ComplexI = {static_cast<T>(0),static_cast<T>(1)};

NAV_DEFINE_FP_CONSTANT(LIGHT_SPEED, 299792458.0);

} // namespace navtools
#endif
