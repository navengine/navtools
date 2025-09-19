/**
 * *constants.hpp*
 *
 * =======  ========================================================================================
 * @file    include/navtools/constants.hpp
 * @brief   Useful navigation constants.
 * @author  Blake Baker, Daniel Sturdivant
 * @ref     Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems
 *            - (2013) Paul D. Groves
 * @date    March 2025
 * =======  ========================================================================================
 */

#ifndef NAVTOOLS_CONSTANTS_HPP
#define NAVTOOLS_CONSTANTS_HPP

#include <cmath>
#include <complex>
#include <numbers>
#include <cassert>

#include <Eigen/Dense>

#include "navtools/types.hpp"

namespace navtools {

//* ===== Physical Constants =================================================================== *//

// DEFINE_FP_CONSTANT(PI, 3.141592653589793238462643383279502884);  //! pi
template <typename T = double>
inline constexpr std::complex<T> COMPLEX_I{
    static_cast<T>(0), static_cast<T>(1)};                         //! imaginary value
DEFINE_FP_CONSTANT(PI, std::numbers::pi_v<T>);                     //! pi
DEFINE_FP_CONSTANT(HALF_PI, 0.5 * PI<T>);                          //! pi/2
DEFINE_FP_CONSTANT(TWO_PI, 2.0 * PI<T>);                           //! 2*pi
DEFINE_FP_CONSTANT(PI_SQU, PI<T>* PI<T>);                          //! pi^2
DEFINE_FP_CONSTANT(SQRT_PI, 1.0 / std::numbers::inv_sqrtpi_v<T>);  //! sqrt(pi)
DEFINE_FP_CONSTANT(RAD2DEG, 180.0 / PI<T>);                        //! radians to degrees
DEFINE_FP_CONSTANT(DEG2RAD, PI<T> / 180.0);                        //! degrees to radians
template <typename T = double>
inline static const Eigen::Vector<T, 3> LLA_RAD2DEG{RAD2DEG<T>, RAD2DEG<T>, 1.0};
template <typename T = double>
inline static const Eigen::Vector<T, 3> LLA_DEG2RAD{DEG2RAD<T>, DEG2RAD<T>, 1.0};

template <int Numerator, int Denominator = 1, typename T = double>
constexpr T PiFraction() {
  return static_cast<T>(Numerator) * PI<T> / static_cast<T>(Denominator);
}

template <int Numerator, int Denominator = 1, typename T = double>
constexpr T Fraction() {
  return static_cast<T>(Numerator) / static_cast<T>(Denominator);
}

// CAUTION: max val allowable for Pow depends on the bit size of intmax_t. It should be at least 64 on most systems.
template<int Pow, typename Float = double>
static constexpr Float PowerOfTwo()
{
  static_assert((Pow < 64) && (Pow > -64), "PowerOfTwo: Pow of magnitude greater than 63 cannot be used. Use a different method.");
  if constexpr (Pow == 0) return Float(1);
  else if constexpr (Pow < 0) {
    return Float(1) / Float(uint64_t(1) << -Pow);
  }
  else {
    return Float(uint64_t(1) << Pow);
  }
}

DEFINE_FP_CONSTANT(LIGHT_SPEED, 299792458.0);  //! speed of light [m/s]
DEFINE_FP_CONSTANT(BOLTZMANN, 1.38e-23);       //! Boltsman constant [J/K]
DEFINE_FP_CONSTANT(GAUSS_TO_TESLA, 1e-4);      //! Gauss to Tesla
DEFINE_FP_CONSTANT(METERS_PER_FOOT, 0.3048);   //! Feet to meters
DEFINE_FP_CONSTANT(MINUTES_PER_DAY, 1440.0);   //! minutes per day

//* ===== Earth Info =========================================================================== *//

DEFINE_FP_CONSTANT(GRAVITY, 9.80665);     //! Earth's gravity constant [m/s^2]
DEFINE_FP_CONSTANT(RE, 6378.1363e3);      //! Earth equatorial radius [m]
DEFINE_FP_CONSTANT(J2, 1.0826269e-03);    //! Earth second zonal harmonic coefficient
DEFINE_FP_CONSTANT(J3, -2.5323000e-06);   //! Earth third zonal harmonic coefficient
DEFINE_FP_CONSTANT(J4, -1.6204000e-06);   //! Earth fourth zonal harmonic coefficient
DEFINE_FP_CONSTANT(F, -4.442807633e-10);  //! Relativistic coefficient

//* ===== WGS84 Models ========================================================================= *//

// DEFINE_FP_CONSTANT(WGS84_MU, 3.986004415e14);     //! (mu) Gravitational constant [m^3/s^2]
DEFINE_FP_CONSTANT(WGS84_MU, 3.986005e14);
DEFINE_FP_CONSTANT(WGS84_R0, 6378137.0);          //! WGS84 Equatorial radius (semi-major axis) [m]
DEFINE_FP_CONSTANT(WGS84_RP, 6356752.314245);     //! WGS84 Polar radius (semi-minor axis) [m]
DEFINE_FP_CONSTANT(WGS84_E, 0.0818191908429654);  //! WGS84 eccentricity
DEFINE_FP_CONSTANT(WGS84_E2, 0.00669437999019758);  //! WGS84 eccentricity squared
DEFINE_FP_CONSTANT(WGS84_F, 0.00335281066477569);   //! WGS84 inverse flattening
DEFINE_FP_CONSTANT(WGS84_OMEGA, 7.2921151467e-5);   //! WGS84 Earth rotational constant
template <typename T = double>
inline static const Eigen::Vector<T, 3> WGS84_OMEGA_VEC{0.0, 0.0, WGS84_OMEGA<T>};
template <typename T = double>
inline static const Eigen::Matrix<T, 3, 3> WGS84_OMEGA_SKEW{
    {0.0, 0.0, 0.0},
    {0.0, 0.0, +WGS84_OMEGA<T>},
    {0.0, -WGS84_OMEGA<T>, 0.0},
};

}  // namespace navtools

#endif
