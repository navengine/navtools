#ifndef NAVTOOLS_CORE_CONSTANTS_HPP
#define NAVTOOLS_CORE_CONSTANTS_HPP

#include "navtools/core/macros.hpp"
#include <Eigen/Dense>
#include <numbers>

namespace nt {

namespace num = std::numbers;

// *=== Physical constants ===*
NEW_FP_CONST(SQRT2, num::sqrt2_v<T>);               // sqrt(2)
NEW_FP_CONST(SQRT3, num::sqrt3_v<T>);               // sqrt(2)
NEW_FP_CONST(SQRT_HALF, num::sqrt2_v<T> / 2.0);     // sqrt(0.5)
NEW_FP_CONST(PI, num::pi_v<T>);                     // pi
NEW_FP_CONST(HALF_PI, 0.5 * num::pi_v<T>);          // pi/2
NEW_FP_CONST(TWO_PI, 2.0 * num::pi_v<T>);           // 2*pi
NEW_FP_CONST(PI_SQU, num::pi_v<T>* num::pi_v<T>);   // pi^2
NEW_FP_CONST(SQRT_PI, 1.0 / num::inv_sqrtpi_v<T>);  // sqrt(pi)
NEW_FP_CONST(R2D, 180.0 * num::inv_pi_v<T>);        // radians to degrees
NEW_FP_CONST(D2R, num::pi_v<T> / 180.0);            // degrees to radians
NEW_FP_CONST(BOLTZMANN, 1.38e-23);                  // Boltzman constant (J/K)
NEW_FP_CONST(GAUSS2TESLA, 1e-4);                    // Gauss to Tesla
NEW_FP_CONST(FT2M, 381.0 / 1250.0);                 // Feet to meters
NEW_FP_CONST(M2FT, 1250.0 / 381.0);                 // meters to feet
NEW_FP_CONST(MIN2DAY, 1440.0);                      // minutes to day
NEW_FP_CONST(MIN2HR, 60.0);                         // minutes to hour
NEW_FP_CONST(SEC2DAY, 86400.0);                     // seconds to day
NEW_FP_CONST(SEC2HR, 3600.0);                       // seconds to hour
NEW_FP_CONST(SEC2MIN, 60.0);                        // seconds to minutes
NEW_FP_CONST(LIGHT_SPEED, 299792458.0);             // speed of light (m/s)
NEW_FP_CONST(GRAVITY, 9.80665);                     // gravity (m/s^2)

template <int Pow, typename Float = double>
static constexpr Float PowerOfTwo() {
  static_assert(
      (Pow < 8 * sizeof(Float)) && (Pow > -8 * sizeof(Float)),
      "PowerOfTwo: Pow of magnitude greater than the number of bits cannot be used. Use a "
      "different method.");
  if constexpr (Pow == 0)
    return Float(1);
  else if constexpr (Pow < 0) {
    return Float(1) / Float(uint64_t(1) << -Pow);
  } else {
    return Float(uint64_t(1) << Pow);
  }
}

// *=== WGS84 constants ===*
NEW_FP_CONST(WGS84_A, 6378137.0);             // WGS84 Equatorial radius (semi-major axis) [m]
NEW_FP_CONST(WGS84_B, 6356752.3142);          // WGS84 Polar radius (semi-minor axis) [m]
NEW_FP_CONST(WGS84_E, 0.081819190842622);     // WGS84 eccentricity
NEW_FP_CONST(WGS84_E2, 6.69437999014e-3);     // WGS84 eccentricity squared
NEW_FP_CONST(WGS84_FINV, 298.257223563);      // WGS84 inverse geometrical flattening
NEW_FP_CONST(WGS84_F, 0.082094437949696);     // WGS84 geometrical flattening
NEW_FP_CONST(WGS84_M, 0.00344978650684);      // WGS84 physical flattening (gravity ratio)
NEW_FP_CONST(WGS84_OMEGA, 7.292115e-5);       // WGS84 Earth rotational constant
NEW_FP_CONST(WGS84_GM, 3.986004418e14);       // WGS84 gravitational constant [m^3/s^2]
NEW_FP_CONST(WGS84_EM, 5.9733328e24);         // WGS84 mass of the Earth [kg]
NEW_FP_CONST(WGS84_GA, 9.7803253359);         // WGS84 gravity at equator [m/s^2]
NEW_FP_CONST(WGS84_GB, 9.8321849378);         // WGS84 gravity at pole [m/s^2]
NEW_FP_CONST(WGS84_SOMGLIANA, 0.001931853);   // WGS84 Somigliana constant
NEW_FP_CONST(WGS84_GRAVITY, 9.7803253359);    // Earth's gravity constant [m/s^2]
NEW_FP_CONST(WGS84_RE, 6378.1363e3);          // Earth equatorial radius [m]
NEW_FP_CONST(WGS84_J2, 1.0826269e-03);        // Earth second zonal harmonic coefficient
NEW_FP_CONST(WGS84_J3, -2.5323000e-06);       // Earth third zonal harmonic coefficient
NEW_FP_CONST(WGS84_J4, -1.6204000e-06);       // Earth fourth zonal harmonic coefficient
NEW_FP_CONST(WGS84_REL_F, -4.442807633e-10);  // Relativistic coefficient

NEW_EIGEN_CONST(OMEGA_ECEF, 3, 1, 0.0, 0.0, WGS84_OMEGA<T>);
NEW_EIGEN_CONST(
    OMEGA_ECEF_SKEW, 3, 3, 0.0, -WGS84_OMEGA<T>, 0.0, WGS84_OMEGA<T>, 0.0, 0.0, 0.0, 0.0, 0.0);

}  // namespace nt

#endif