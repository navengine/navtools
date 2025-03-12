/**
 * *models.hpp*
 *
 * =======  ========================================================================================
 * @file    include/navtools/models.hpp
 * @brief   Simple Earth models commonly used in navigation equations.
 * @author  Daniel Sturdivant, Blake Baker
 * @ref     Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems
 *            - (2013) Paul D. Groves
 * @date    March 2025
 * =======  ========================================================================================
 */
// TODO: maybe split into earth_models

#ifndef NAVTOOLS_MODELS_HPP
#define NAVTOOLS_MODELS_HPP

#include <Eigen/Dense>

#include "navtools/constants.hpp"

namespace navtools {

//* ===== Radii of Curvature =================================================================== *//

//! === TRANSVERSERADIUS ===
/// @brief      Calculates the transverse radius relative to user latitude
/// @param phi  Latitude [rad]
/// @param Re   Transverse radius [m]
/// @returns    Earth's transverse radius at Latitude
template <typename T = double>
void TransverseRadius(T &Re, const T &phi) {
  T sin_phi = std::sin(phi);
  T t = 1.0 - WGS84_E2<T> * sin_phi * sin_phi;
  Re = WGS84_R0<T> / std::sqrt(t);
}
template <typename T = double>
T TransverseRadius(const T &phi) {
  T Re;
  TransverseRadius<T>(Re, phi);
  return Re;
}

//! === MERIDIANRADIUS ===
/// @brief      Calculates the meridian radius relative to user latitude
/// @param phi  Latitude [rad]
/// @param Rn   Meridian radius [m]
/// @returns    Earth's meridian radius at Latitude
template <typename T = double>
void MeridianRadius(T &Rn, const T &phi) {
  T sin_phi = std::sin(phi);
  T t = 1.0 - WGS84_E2<T> * sin_phi * sin_phi;
  Rn = WGS84_R0<T> * (1.0 - WGS84_E2<T>) / std::pow(t, 1.5);
}
template <typename T = double>
T MeridianRadius(const T &phi) {
  T Rn;
  MeridianRadius<T>(Rn, phi);
  return Rn;
}

//! === GEOCENTRICRADIUS ===
/// @brief      Calculates the geocentric radius relative to user latitude
/// @param phi      Latitude [rad]
/// @param R_es_e   Geocentric radius [m]
/// @returns    Earth's geocentric radius at Latitude
template <typename T = double>
void GeocentricRadius(T &R_es_e, const T &phi) {
  T sin_phi2 = std::sin(phi);
  sin_phi2 *= sin_phi2;
  T cos_phi = std::cos(phi);
  T t = 1.0 - WGS84_E2<T> * sin_phi2;
  T Re = WGS84_R0<T> / std::sqrt(t);
  T o_e2 = 1.0 - WGS84_E2<T>;
  R_es_e = Re * std::sqrt(cos_phi * cos_phi + o_e2 * o_e2 * sin_phi2);
}
template <typename T = double>
T GeocentricRadius(const T &phi) {
  T R_es_e;
  GeocentricRadius<T>(R_es_e, phi);
  return R_es_e;
}

//! === TRANSANDMERRADII ===
/// @brief      Calculates the {Transverse, Meridian} radii relative to user latitude
/// @param phi      Latitude [rad]
/// @param radii    vector containing {Transverse, Meridian} radii [m]
/// @returns    Earth's {Transverse, Meridian} radii at Latitude
template <typename T = double>
void TransAndMerRadii(Eigen::Ref<Eigen::Vector<T, 2>> radii, const T &phi) {
  T sin_phi = std::sin(phi);
  T t = 1.0 - WGS84_E2<T> * sin_phi * sin_phi;

  radii(0) = WGS84_R0<T> / std::sqrt(t);
  radii(1) = WGS84_R0<T> * (1.0 - WGS84_E2<T>) / std::pow(t, 1.5);
}
template <typename T = double>
void TransAndMerRadii(T &Re, T &Rn, const T &phi) {
  T sin_phi = std::sin(phi);
  T t = 1.0 - WGS84_E2<T> * sin_phi * sin_phi;

  Re = WGS84_R0<T> / std::sqrt(t);
  Rn = WGS84_R0<T> * (1.0 - WGS84_E2<T>) / std::pow(t, 1.5);
}
template <typename T = double>
Eigen::Vector<T, 2> TransAndMerRadii(const T &phi) {
  Eigen::Vector<T, 2> radii;
  TransAndMerRadii<T>(radii, phi);
  return radii;
}

//! === RADIIOFCURVATURE ===
/// @brief      Calculates the {Transverse, Meridian, Geocentric} radii relative to user latitude
/// @param phi      Latitude [rad]
/// @param radii    vector containing {Transverse, Meridian, Geocentric} radii [m]
/// @returns    Earth's {Transverse, Meridian, Geocentric} radii at Latitude
template <typename T = double>
void RadiiOfCurvature(Eigen::Ref<Eigen::Vector<T, 3>> radii, const T &phi) {
  T sin_phi2 = std::sin(phi);
  sin_phi2 *= sin_phi2;
  T cos_phi = std::cos(phi);
  T t = 1.0 - WGS84_E2<T> * sin_phi2;
  T o_e2 = 1.0 - WGS84_E2<T>;

  radii(0) = WGS84_R0<T> / std::sqrt(t);
  radii(1) = WGS84_R0<T> * o_e2 / std::pow(t, 1.5);
  radii(2) = radii(0) * std::sqrt(cos_phi * cos_phi + o_e2 * o_e2 * sin_phi2);
}
template <typename T = double>
void RadiiOfCurvature(T &Re, T &Rn, T &R_es_e, const T &phi) {
  T sin_phi2 = std::sin(phi);
  sin_phi2 *= sin_phi2;
  T cos_phi = std::cos(phi);
  T t = 1.0 - WGS84_E2<T> * sin_phi2;
  T o_e2 = 1.0 - WGS84_E2<T>;

  Re = WGS84_R0<T> / std::sqrt(t);
  Rn = WGS84_R0<T> * o_e2 / std::pow(t, 1.5);
  R_es_e = Re * std::sqrt(cos_phi * cos_phi + o_e2 * o_e2 * sin_phi2);
}
template <typename T = double>
Eigen::Vector<T, 3> RadiiOfCurvature(const T &phi) {
  Eigen::Vector<T, 3> radii;
  RadiiOfCurvature<T>(radii, phi);
  return radii;
}

//* ===== Coriolis ============================================================================= *//

//! === EARTHRATE ===
/// @brief      Rotation rate of the earth relative to the 'NAV' frame
/// @param phi      Latitude [rad]
/// @param frame    string representing the NAV-frame to rotate into
/// @param w_ie_n   size 3 vector of earth's rotation in the 'NAV' frame
/// @param W_ie_n   3x3 skew-symmetric matrix of earth's rotation in the 'NAV' frame
/// @returns    earth's rotation in the 'NAV' frame
template <typename T = double>
void EarthRate(Eigen::Ref<Eigen::Vector<T, 3>> w_ie_n, const T &phi, const bool isNed = true) {
  if (isNed) {
    w_ie_n(0) = WGS84_OMEGA<T> * std::cos(phi);
    w_ie_n(1) = 0.0;
    w_ie_n(2) = WGS84_OMEGA<T> * std::sin(phi);
  } else {
    w_ie_n(0) = 0.0;
    w_ie_n(1) = WGS84_OMEGA<T> * std::cos(phi);
    w_ie_n(2) = -WGS84_OMEGA<T> * std::sin(phi);
  }
}
template <typename T = double>
Eigen::Vector<T, 3> EarthRate(const T &phi, const bool isNed = true) {
  Eigen::Vector<T, 3> w_ie_n;
  EarthRate<T>(w_ie_n, phi, isNed);
  return w_ie_n;
}
template <typename T = double>
void EarthRateSkew(Eigen::Matrix<T, 3, 3> &W_ie_n, const T &phi, const bool isNed = true) {
  W_ie_n = Skew<T>(EarthRate<T>(phi, isNed));
}
template <typename T = double>
Eigen::Matrix<T, 3, 3> EarthRateSkew(const T &phi, const bool isNed = true) {
  Eigen::Matrix<T, 3, 3> W_ie_n;
  EarthRateSkew<T>(W_ie_n, phi, isNed);
  return W_ie_n;
}

//! === TRANSPORTRATE ===
/// @brief      Transport rate of the 'ECEF' frame relative to the 'NAV' frame
/// @param lla      Latitude, Longitude, Height [rad, rad, m]
/// @param v_nb_e   size 3 velocity vector in the 'NAV' coordinate system
/// @param frame    string representing the NAV-frame to rotate into
/// @param w_en_n   size 3 vector of earth's rotation in the 'NAV' frame
/// @param W_en_n   3x3 skew-symmetric matrix of earth's rotation in the 'NAV' frame
/// @returns    Transport rate in the 'NAV' frame
template <typename T = double>
void TransportRate(
    Eigen::Ref<Eigen::Vector<T, 3>> w_en_n,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_nb_e,
    const bool isNed = true) {
  Eigen::Vector<T, 2> radii = TransAndMerRadii<T>(lla(0));
  if (isNed) {
    T ve_Reh = v_nb_e(1) / (radii(0) + lla(2));    // ve / (Re + h)
    w_en_n(0) = ve_Reh;                            // ve / (Re + h)
    w_en_n(1) = -v_nb_e(0) / (radii(1) + lla(2));  // -vn / (Rn + h)
    w_en_n(2) = -ve_Reh * std::tan(lla(0));        // -ve * tan(phi) / (Re + h)
  } else {
    T ve_Reh = v_nb_e(0) / (radii(0) + lla(2));    // ve / (Re + h)
    w_en_n(0) = -v_nb_e(1) / (radii(1) + lla(2));  // -vn / (Rn + h)
    w_en_n(1) = ve_Reh;                            // ve / (Re + h)
    w_en_n(2) = ve_Reh * std::tan(lla(0));         // ve * tan(phi) / (Re + h)
  }
}
template <typename T = double>
Eigen::Vector<T, 3> TransportRate(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_nb_e,
    const bool isNed = true) {
  Eigen::Vector<T, 3> w_en_n;
  TransportRate<T>(w_en_n, lla, v_nb_e, isNed);
  return w_en_n;
}
template <typename T = double>
void TransportRateSkew(
    Eigen::Matrix<T, 3, 3> &W_en_n,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_nb_e,
    const bool isNed = true) {
  W_en_n = Skew<T>(TransportRate<T>(lla, v_nb_e, isNed));
}
template <typename T = double>
Eigen::Matrix<T, 3, 3> TransportRateSkew(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_nb_e,
    const bool isNed = true) {
  Eigen::Matrix<T, 3, 3> W_en_n;
  TransportRateSkew<T>(W_en_n, lla, v_nb_e, isNed);
  return W_en_n;
}

//! === CORIOLISRATE ===
/// @brief      Coriolis effect perceived in the "NAV" frame
/// @param lla      Latitude, Longitude, Height [rad, rad, m]
/// @param v_nb_e   size 3 velocity vector in the 'NAV' coordinate system
/// @param frame    string representing the NAV-frame to rotate into
/// @param coriolis size 3 coriolis effect
/// @returns    Coriolis effect
template <typename T = double>
void CoriolisRate(
    Eigen::Ref<Eigen::Vector<T, 3>> coriolis,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_nb_e,
    const bool isNed = true) {
  Eigen::Vector<T, 3> w_ie_n = EarthRate<T>(lla(0), isNed);
  Eigen::Vector<T, 3> w_en_n = TransportRate<T>(lla, v_nb_e, isNed);
  coriolis = Skew<T>(w_en_n + 2.0 * w_ie_n) * v_nb_e;
}
template <typename T = double>
Eigen::Vector<T, 3> CoriolisRate(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_nb_e,
    const bool isNed = true) {
  Eigen::Vector<T, 3> coriolis;
  CoriolisRate<T>(coriolis, lla, v_nb_e, isNed);
  return coriolis;
}

//* ===== Gravity ============================================================================== *//

//! === SOMIGLIANA ===
/// @brief      Calculates the somilgiana model reference gravity
/// @param phi  Latitude [rad]
/// @param g0   Somigliana gravity
/// @returns Somgiliana gravity
template <typename T = double>
void Somigliana(T &g0, const T &phi) {
  T sin_phi2 = std::sin(phi);
  sin_phi2 *= sin_phi2;
  g0 = 9.7803253359 * ((1.0 + 0.001931853 * sin_phi2) / std::sqrt(1.0 - WGS84_E2<T> * sin_phi2));
}
template <typename T = double>
T Somigliana(const T &phi) {
  T g0;
  Somigliana<T>(g0, phi);
  return g0;
}

//! === LOCALGRAVITY ===
/// @brief      Calculates gravity in the Local/NAV (ENU or NED) frame
/// @param lla      Latitude, Longitude, Height [rad, rad, m]
/// @param frame    string representing the NAV-frame to rotate into
/// @param g        size 3 Local/NAV frame gravity vector
/// @returns    Local/NAV frame gravity
template <typename T = double>
void LocalGravity(
    Eigen::Ref<Eigen::Vector<T, 3>> g,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla,
    const bool isNed = true) {
  T sin_phi2 = std::sin(lla(0));
  sin_phi2 *= sin_phi2;
  T g0 = 9.7803253359 * ((1.0 + 0.001931853 * sin_phi2) / std::sqrt(1.0 - WGS84_E2<T> * sin_phi2));
  T R02 = WGS84_R0<T> * WGS84_R0<T>;
  T OMEGA2 = WGS84_OMEGA<T> * WGS84_OMEGA<T>;
  T h2 = lla(2) * lla(2);
  if (isNed) {
    // clang-format off
        g(0) = -8.08e-9 * lla(2) * std::sin(2.0 * lla(0));
        g(1) = 0.0;
        g(2) = g0 * (1.0 - 
                    (2.0 / WGS84_R0<T>) *
                    (1.0 + WGS84_F<T> * (1.0 - 2.0 * sin_phi2) + (OMEGA2 * R02 * WGS84_RP<T> / WGS84_MU<T>)) *
                    lla(2) + (3.0 * h2 / R02));
    // clang-format on
  } else {
    // clang-format off
        g(0) = 0.0;
        g(1) = -8.08e-9 * lla(2) * std::sin(2.0 * lla(0));
        g(2) = -g0 * (1.0 - 
                     (2.0 / WGS84_R0<T>) *
                     (1.0 + WGS84_F<T> * (1.0 - 2.0 * sin_phi2) + (OMEGA2 * R02 * WGS84_RP<T> / WGS84_MU<T>)) *
                     lla(2) + (3.0 * h2 / R02));
    // clang-format on
  }
}
template <typename T = double>
Eigen::Vector<T, 3> LocalGravity(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla, const bool isNed = true) {
  Eigen::Vector<T, 3> g;
  LocalGravity<T>(g, lla, isNed);
  return g;
}

//! === ECEFGRAVITY ===
/// @brief      Calculates gravity in the Earth-Centered-Earth-Fixed frame
/// @param xyz      ECEF position [m]
/// @param g        size 3 ECEF frame gravity vector
/// @param gamma    size 3 ECEF frame gravitational acceleration
/// @returns    ECEF frame gravity
template <typename T = double>
void EcefGravity(
    Eigen::Ref<Eigen::Vector<T, 3>> g,
    Eigen::Ref<Eigen::Vector<T, 3>> gamma,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &xyz) {
  T mag_r = xyz.norm();
  if (mag_r == 0.0) {
    g << 0.0, 0.0, 0.0;
    gamma << 0.0, 0.0, 0.0;
  } else {
    T zeta = 5.0 * std::pow(xyz(2) / mag_r, 2.0);
    T omega2 = WGS84_OMEGA<T> * WGS84_OMEGA<T>;
    Eigen::Vector<T, 3> v{1.0 - zeta, 1.0 - zeta, 3.0 - zeta};

    gamma =
        -WGS84_MU<T> / std::pow(mag_r, 3.0) *
        (xyz.array() + 1.5 * J2<T> * std::pow(WGS84_R0<T> / mag_r, 2.0) * v.array() * xyz.array());
    v << xyz(0) * omega2, xyz(1) * omega2, 0.0;
    g = gamma + v;
  }
}
template <typename T = double>
void EcefGravity(
    Eigen::Ref<Eigen::Vector<T, 3>> g, const Eigen::Ref<const Eigen::Vector<T, 3>> &xyz) {
  Eigen::Vector<T, 3> gamma;
  EcefGravity<T>(g, gamma, xyz);
}
template <typename T = double>
Eigen::Vector<T, 3> EcefGravity(const Eigen::Ref<const Eigen::Vector<T, 3>> &xyz) {
  Eigen::Vector<T, 3> g;
  EcefGravity<T>(g, xyz);
  return g;
}

//* ===== 3d Ranging =========================================================================== *//

//! === CALCRANGE ===
/// @brief      Computes the range from user to satellite
/// @param sv_xyz   Satellite position [m]
/// @param user_xyz User position [m]
/// @param r        Calculated range [m]
/// @returns    Range
template <typename T = double>
void CalcRange(
    T &r,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &sv_xyz,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &user_xyz) {
  r = (user_xyz - sv_xyz).norm();
}
template <typename T = double>
T CalcRange(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &sv_xyz,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &user_xyz) {
  T r;
  CalcRange<T>(r, sv_xyz, user_xyz);
  return r;
}

//! === CALCUNITVEC ===
/// @brief      Computes the unit vector from user to satellite
/// @param sv_xyz   Satellite position [m]
/// @param user_xyz User position [m]
/// @param u        Calculated unit vector [m]
/// @returns    Unit vector
template <typename T = double>
void CalcUnitVec(
    Eigen::Ref<Eigen::Vector<T, 3>> u,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &sv_xyz,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &user_xyz) {
  Eigen::Vector<T, 3> dr = user_xyz - sv_xyz;
  u = dr / dr.norm();
}
template <typename T = double>
Eigen::Vector<T, 3> CalcUnitVec(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &sv_xyz,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &user_xyz) {
  Eigen::Vector<T, 3> u;
  CalcUnitVec<T>(u, sv_xyz, user_xyz);
  return u;
}

//! === CALCRANGEANDUNITVEC ===
/// @brief      Computes the range and unit vector from user to satellite
/// @param sv_xyz   Satellite position [m]
/// @param user_xyz User position [m]
/// @param r        Calculated range [m]
/// @param u        Calculated unit vector [m]
/// @returns    Range and unit vector
template <typename T = double>
void CalcRangeAndUnitVec(
    T &r,
    Eigen::Ref<Eigen::Vector<T, 3>> u,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &sv_xyz,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &user_xyz) {
  Eigen::Vector<T, 3> dr = user_xyz - sv_xyz;
  r = dr.norm();
  u = dr / r;
}

//! === CALCRANGERATE ===
/// @brief      Computes the range-rate between the user and satellite
/// @param u        Unit vector [m]
/// @param sv_vel   Satellite velocity [m/s]
/// @param user_vel User velocity [m/s]
/// @param rr       Calculated range-rate [m/s]
template <typename T = double>
void CalcRangeRate(
    T &rr,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &u,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &sv_vel,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &user_vel) {
  Eigen::Vector<T, 3> dv = user_vel - sv_vel;
  rr = u(0) * dv(0) + u(1) * dv(1) + u(2) * dv(2);
}
template <typename T = double>
T CalcRangeRate(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &u,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &sv_vel,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &user_vel) {
  T rr;
  CalcRangeRate<T>(rr, u, sv_vel, user_vel);
  return rr;
}

}  // namespace navtools

#endif
