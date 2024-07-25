/**
|======================================== constants.hpp ===========================================|
|                                                                                                  |
|   @file     include/navtools/constants.hpp                                                       |
|   @brief    Useful navigation constants.                                                         |
|   @ref      Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems          |
|               - (2013) Paul D. Groves                                                            |
|                                                                                                  |
|==================================================================================================|
*/

#ifndef NAVTOOLS_CONSTANTS_HPP
#define NAVTOOLS_CONSTANTS_HPP

#include <cmath>
#include <complex>
#include <numbers>

#include "navtools/types.hpp"

namespace navtools {

//* ===== Physical Constants =================================================================== *//

template <typename T = double>
inline constexpr std::complex<T> J{static_cast<T>(0), static_cast<T>(1)};  //! imaginary value
template <typename T = double>                                             //
inline constexpr T PI = std::numbers::pi_v<T>;                             //! pi
template <typename T = double>                                             //
inline constexpr T HALF_PI = std::numbers::pi_v<T> / static_cast<T>(2);    //! pi/2
template <typename T = double>                                             //
inline constexpr T TWO_PI = static_cast<T>(2) * std::numbers::pi_v<T>;     //! 2*pi
template <typename T = double>                                             //
inline constexpr T PI_SQU = std::pow(std::numbers::pi_v<T>, 2);            //! pi^2
template <typename T = double>                                             //
inline constexpr T SQRT_PI = std::sqrt(std::numbers::pi_v<T>);             //! sqrt(pi)
template <typename T = double>                                             //
inline constexpr T RAD2DEG = 180.0 / std::numbers::pi_v<T>;                //! radians to degrees
template <typename T = double>                                             //
inline constexpr T DEG2RAD = std::numbers::pi_v<T> / 180.0;                //! degrees to radians
template <typename T = double>
inline static const Vec3<T> LLA_RAD2DEG{RAD2DEG<T>, RAD2DEG<T>, 1.0};
template <typename T = double>
inline static const Vec3<T> LLA_DEG2RAD{DEG2RAD<T>, DEG2RAD<T>, 1.0};

template <int Numerator, int Denominator, typename T = double>
constexpr T PiFraction() {
    return static_cast<T>(Numerator) * PI<T> / static_cast<T>(Denominator);
}

template <int Numerator, int Denominator, typename T = double>
constexpr T Fraction() {
    return static_cast<T>(Numerator) / static_cast<T>(Denominator);
}

DEFINE_FP_CONSTANT(LIGHT_SPEED, 299792458.0);  //! speed of light [m/s]
DEFINE_FP_CONSTANT(BOLTZMANN, 1.38e-23);       //! Boltsman constant [J/K]
DEFINE_FP_CONSTANT(GAUSS_TO_TESLA, 1e-4);      //! Gauss to Tesla
DEFINE_FP_CONSTANT(METERS_PER_FOOT, 0.3048);   //! Feet to meters
DEFINE_FP_CONSTANT(MINUTES_PER_DAY, 1440.0);   //! minutes per day

//* ===== Earth Info =========================================================================== *//

DEFINE_FP_CONSTANT(GRAVITY, 9.80665);     //! Earth's gravity constant [m/s^2]
DEFINE_FP_CONSTANT(RE, 6378.1363e3);      //! Earth equatorial radius [m]
DEFINE_FP_CONSTANT(GM, 3.986004415e14);   //! (mu) Gravitational constant [m^3/s^2]
DEFINE_FP_CONSTANT(J2, 1.0826269e-03);    //! Earth second zonal harmonic coefficient
DEFINE_FP_CONSTANT(J3, -2.5323000e-06);   //! Earth third zonal harmonic coefficient
DEFINE_FP_CONSTANT(J4, -1.6204000e-06);   //! Earth fourth zonal harmonic coefficient
DEFINE_FP_CONSTANT(F, -4.442807633e-10);  //! Relativistic coefficient

//* ===== WGS84 Models ========================================================================= *//

DEFINE_FP_CONSTANT(WGS84_R0, 6378137.0);          //! WGS84 Equatorial radius (semi-major axis) [m]
DEFINE_FP_CONSTANT(WGS84_RP, 6356752.314245);     //! WGS84 Polar radius (semi-minor axis) [m]
DEFINE_FP_CONSTANT(WGS84_E, 0.0818191908429654);  //! WGS84 eccentricity
DEFINE_FP_CONSTANT(WGS84_E2, 0.00669437999019758);  //! WGS84 eccentricity squared
DEFINE_FP_CONSTANT(WGS84_F, 0.00335281066477569);   //! WGS84 inverse flattening
DEFINE_FP_CONSTANT(WGS84_OMEGA, 7.2921151467e-5);   //! WGS84 Earth rotational constant
template <typename T = double>
inline static const Vec3<T> WGS84_OMEGA_VEC{0.0, 0.0, WGS84_OMEGA<T>};
template <typename T = double>
inline static const Mat3x3<T> WGS84_OMEGA_SKEW{
        {0.0, -WGS84_OMEGA<T>, 0.0},
        {WGS84_OMEGA<T>, 0.0, 0.0},
        {0.0, 0.0, 0.0},
};

//* ===== SGP Models =========================================================================== *//

DEFINE_FP_CONSTANT(SGP_AE, 1.0);                      // distance per earth radii
DEFINE_FP_CONSTANT(SGP_XKMPER, 6378.135);             // km per earth radii
DEFINE_FP_CONSTANT(SGP_S, 1.01222928);                // s
DEFINE_FP_CONSTANT(SGP_QOMS2T, 1.88027916e-9);        // (q0 - s)^4
DEFINE_FP_CONSTANT(SGP_XKE, 7.43669161331734132e-2);  // sqrt(G*M)
template <typename T = double>                        //
inline constexpr T SGP_CK2 = 1.0 / 2.0 * J2<T> * std::pow(SGP_AE<T>, 2);   // 1/2 * J2 * a_E^2
template <typename T = double>                                             //
inline constexpr T SGP_CK4 = -3.0 / 8.0 * J4<T> * std::pow(SGP_AE<T>, 4);  // -3/8 * J4 * a_E^4
template <typename T = double>                                             //
inline constexpr T SGP_A3OVK2 = -J3<T> / SGP_CK2<T> * SGP_AE<T> * std::pow(SGP_AE<T>, 2);  //

}  // namespace navtools

#endif
