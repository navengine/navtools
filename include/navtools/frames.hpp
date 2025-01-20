/**
|========================================== frames.hpp ============================================|
|                                                                                                  |
|   @file     include/navtools/frames.hpp                                                          |
|   @brief    Common coordinate frame transformations.                                             |
|   @ref      Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems          |
|               - (2013) Paul D. Groves                                                            |
|   @date     July 2024                                                                            |
|                                                                                                  |
|==================================================================================================|
*/

#ifndef NAVTOOLS_FRAMES_HPP
#define NAVTOOLS_FRAMES_HPP

#include <Eigen/Dense>
#include <cmath>
#include <navtools/constants.hpp>
#include <navtools/types.hpp>

namespace navtools {

//* ===== Direction Cosine Matrices ============================================================ *//

//! --- ECI2ECEFDCM ---
/// @brief      Earth-Centered-Inertial to Earth-Centered-Earth-Fixed direction cosine matrix
/// @param dt   time elapsed between frames [s]
/// @returns    3x3 ECI->ECEF direction cosine matrix
template <typename T = double>
void eci2ecefDcm(Eigen::Ref<Eigen::Matrix<T, 3, 3>> C, const T &dt) {
  T omega_dt = WGS84_OMEGA<T> * dt;
  T sin_omega_dt = std::sin(omega_dt);
  T cos_omega_dt = std::cos(omega_dt);
  // clang-format off
  C <<  cos_omega_dt, sin_omega_dt, 0.0, 
       -sin_omega_dt, cos_omega_dt, 0.0, 
                 0.0,          0.0, 1.0;
  // clang-format on
}
template <typename T = double>
Eigen::Matrix<T, 3, 3> eci2ecefDcm(const T &dt) {
  Eigen::Matrix<T, 3, 3> C;
  eci2ecefDcm<T>(C, dt);
  return C;
}

//! --- ECI2NEDDCM ---
/// @brief      Earth-Centered-Inertial to North-East-Down direction cosine matrix
/// @param lla  3x1 Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @param dt   time elapsed between frames [s]
/// @returns    3x3 ECI->NED direction cosine matrix
template <typename T = double>
void eci2nedDcm(
    Eigen::Ref<Eigen::Matrix<T, 3, 3>> C,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla,
    const T &dt) {
  T omega_dt = WGS84_OMEGA<T> * dt;
  T sin_phi = std::sin(lla(0));
  T cos_phi = std::cos(lla(0));
  T sin_lam_omega_dt = std::sin(lla(1) + omega_dt);
  T cos_lam_omega_dt = std::cos(lla(1) + omega_dt);
  // clang-format off
  C << -sin_phi * cos_lam_omega_dt, -sin_phi * sin_lam_omega_dt,  cos_phi,
                 -sin_lam_omega_dt,            cos_lam_omega_dt,      0.0,
       -cos_phi * cos_lam_omega_dt, -cos_phi * sin_lam_omega_dt, -sin_phi;
  // clang-format on
}
template <typename T = double>
Eigen::Matrix<T, 3, 3> eci2nedDcm(const Eigen::Ref<const Eigen::Vector<T, 3>> &lla, const T &dt) {
  Eigen::Matrix<T, 3, 3> C;
  eci2nedDcm<T>(C, lla, dt);
  return C;
}

//! --- ECI2ENUDCM ---
/// @brief      Earth-Centered-Inertial to East-North-Up direction cosine matrix
/// @param lla  3x1 Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @param dt   time elapsed between frames [s]
/// @returns    3x3 ECI->ENU direction cosine matrix
template <typename T = double>
void eci2enuDcm(
    Eigen::Ref<Eigen::Matrix<T, 3, 3>> C,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla,
    const T &dt) {
  T omega_dt = WGS84_OMEGA<T> * dt;
  T sin_phi = std::sin(lla(0));
  T cos_phi = std::cos(lla(0));
  T sin_lam_omega_dt = std::sin(lla(1) + omega_dt);
  T cos_lam_omega_dt = std::cos(lla(1) + omega_dt);
  // clang-format off
  C <<           -sin_lam_omega_dt,            cos_lam_omega_dt,     0.0, 
       -sin_phi * cos_lam_omega_dt, -sin_phi * sin_lam_omega_dt, cos_phi, 
        cos_phi * cos_lam_omega_dt,  cos_phi * sin_lam_omega_dt, sin_phi;
          // clang-format off
}
template <typename T = double>
Eigen::Matrix<T,3,3> eci2enuDcm(const Eigen::Vector<T,3> &lla, const T &dt)
{
  Eigen::Matrix<T,3,3> C;
  eci2enuDcm<T>(C, lla, dt);
  return C;
}

//! --- ECEF2ECIDCM ---
/// @brief      Earth-Centered-Earth-Fixed to Earth-Centered-Inertial direction cosine matrix
/// @param dt   time elapsed between frames [s]
/// @returns    3x3 ECEF->ECI direction cosine matrix
template <typename T = double>
void ecef2eciDcm(Eigen::Matrix<T,3,3> &C, const T &dt)
{
  T omega_dt = WGS84_OMEGA<T> * dt;
  T sin_omega_dt = std::sin(omega_dt);
  T cos_omega_dt = std::cos(omega_dt);
  // clang-format off
  C << cos_omega_dt, -sin_omega_dt, 0.0, 
       sin_omega_dt,  cos_omega_dt, 0.0, 
                0.0,          0.0,  1.0;
  // clang-format on
}
template <typename T = double>
Eigen::Matrix<T, 3, 3> ecef2eciDcm(const T &dt) {
  Eigen::Matrix<T, 3, 3> C;
  ecef2eciDcm<T>(C, dt);
  return C;
}

//! --- ECEF2NEDDCM ---
/// @brief      Earth-Centered-Earth-Fixed to North-East-Down direction cosine matrix
/// @param lla  3x1 Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @returns    3x3 ECEF->NED direction cosine matrix
template <typename T = double>
void ecef2nedDcm(
    Eigen::Ref<Eigen::Matrix<T, 3, 3>> C, const Eigen::Ref<const Eigen::Vector<T, 3>> &lla) {
  T sin_phi = std::sin(lla(0));
  T cos_phi = std::cos(lla(0));
  T sin_lam = std::sin(lla(1));
  T cos_lam = std::cos(lla(1));
  // clang-format off
  C << -sin_phi * cos_lam, -sin_phi * sin_lam,  cos_phi,
                 -sin_lam,            cos_lam,      0.0,
       -cos_phi * cos_lam, -cos_phi * sin_lam, -sin_phi;
  // clang-format on
}
template <typename T = double>
Eigen::Matrix<T, 3, 3> ecef2nedDcm(const Eigen::Ref<const Eigen::Vector<T, 3>> &lla) {
  Eigen::Matrix<T, 3, 3> C;
  ecef2nedDcm<T>(C, lla);
  return C;
}

//! --- ECEF2ENUDCM ---
/// @brief      Earth-Centered-Earth-Fixed to East-North-Up direction cosine matrix
/// @param lla  3x1 Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @returns    3x3 ECEF->ENU direction cosine matrix
template <typename T = double>
void ecef2enuDcm(
    Eigen::Ref<Eigen::Matrix<T, 3, 3>> C, const Eigen::Ref<const Eigen::Vector<T, 3>> &lla) {
  T sin_phi = std::sin(lla(0));
  T cos_phi = std::cos(lla(0));
  T sin_lam = std::sin(lla(1));
  T cos_lam = std::cos(lla(1));
  // clang-format off
  C <<           -sin_lam,            cos_lam,     0.0,
       -sin_phi * cos_lam, -sin_phi * sin_lam, cos_phi,
        cos_phi * cos_lam,  cos_phi * sin_lam, sin_phi;
  // clang-format on
}
template <typename T = double>
Eigen::Matrix<T, 3, 3> ecef2enuDcm(const Eigen::Ref<const Eigen::Vector<T, 3>> &lla) {
  Eigen::Matrix<T, 3, 3> C;
  ecef2enuDcm<T>(C, lla);
  return C;
}

//! --- NED2ECIDCM ---
/// @brief      North-East-Down to Earth-Centered-Inertial direction cosine matrix
/// @param lla  3x1 Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @param dt   time elapsed between frames [s]
/// @returns    3x3 ECEF->ECI direction cosine matrix
template <typename T = double>
void ned2eciDcm(
    Eigen::Ref<Eigen::Matrix<T, 3, 3>> C,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla,
    const T &dt) {
  T omega_dt = WGS84_OMEGA<T> * dt;
  T sin_phi = std::sin(lla(0));
  T cos_phi = std::cos(lla(0));
  T sin_lam_omega_dt = std::sin(lla(1) + omega_dt);
  T cos_lam_omega_dt = std::cos(lla(1) + omega_dt);
  // clang-format off
  C << -sin_phi * cos_lam_omega_dt, -sin_lam_omega_dt, -cos_phi * cos_lam_omega_dt,
       -sin_phi * sin_lam_omega_dt,  cos_lam_omega_dt, -cos_phi * sin_lam_omega_dt,
                           cos_phi,               0.0,                    -sin_phi;
  // clang-format on
}
template <typename T = double>
Eigen::Matrix<T, 3, 3> ned2eciDcm(const Eigen::Ref<const Eigen::Vector<T, 3>> &lla, const T &dt) {
  Eigen::Matrix<T, 3, 3> C;
  ned2eciDcm<T>(C, lla, dt);
  return C;
}

//! --- NED2ECEFDCM ---
/// @brief      North-East-Down to Earth-Centered-Earth-Fixed direction cosine matrix
/// @param lla  3x1 Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @returns    3x3 ECEF->NED direction cosine matrix
template <typename T = double>
void ned2ecefDcm(
    Eigen::Ref<Eigen::Matrix<T, 3, 3>> C, const Eigen::Ref<const Eigen::Vector<T, 3>> &lla) {
  T sin_phi = std::sin(lla(0));
  T cos_phi = std::cos(lla(0));
  T sin_lam = std::sin(lla(1));
  T cos_lam = std::cos(lla(1));
  // clang-format off
  C << -sin_phi * cos_lam, -sin_lam, -cos_phi * cos_lam,
       -sin_phi * sin_lam,  cos_lam, -cos_phi * sin_lam,
                  cos_phi,      0.0,           -sin_phi;
  // clang-format on
}
template <typename T = double>
Eigen::Matrix<T, 3, 3> ned2ecefDcm(const Eigen::Ref<const Eigen::Vector<T, 3>> &lla) {
  Eigen::Matrix<T, 3, 3> C;
  ned2ecefDcm<T>(C, lla);
  return C;
}

//! --- NED2ENUDCM ---
/// @brief      North-East-Down to East-North-Up direction cosine matrix
/// @returns    3x3 ECEF->ENU direction cosine matrix
template <typename T = double>
void ned2enuDcm(Eigen::Ref<Eigen::Matrix<T, 3, 3>> C) {
  C << 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0;
}
template <typename T = double>
Eigen::Matrix<T, 3, 3> ned2enuDcm() {
  Eigen::Matrix<T, 3, 3> C;
  ned2enuDcm<T>(C);
  return C;
}

//! --- ENU2ECIDCM ---
/// @brief      East-North-Up to Earth-Centered-Inertial direction cosine matrix
/// @param lla  3x1 Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @param dt   time elapsed between frames [s]
/// @returns    3x3 ECEF->ECI direction cosine matrix
template <typename T = double>
void enu2eciDcm(
    Eigen::Ref<Eigen::Matrix<T, 3, 3>> C,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla,
    const T &dt) {
  T omega_dt = WGS84_OMEGA<T> * dt;
  T sin_phi = std::sin(lla(0));
  T cos_phi = std::cos(lla(0));
  T sin_lam_omega_dt = std::sin(lla(1) + omega_dt);
  T cos_lam_omega_dt = std::cos(lla(1) + omega_dt);
  // clang-format off
  C << -sin_lam_omega_dt, -sin_phi * cos_lam_omega_dt, cos_phi * cos_lam_omega_dt,
        cos_lam_omega_dt, -sin_phi * sin_lam_omega_dt, cos_phi * sin_lam_omega_dt,
                     0.0,                     cos_phi,                    sin_phi;
  // clang-format on
}
template <typename T = double>
Eigen::Matrix<T, 3, 3> enu2eciDcm(const Eigen::Ref<const Eigen::Vector<T, 3>> &lla, const T &dt) {
  Eigen::Matrix<T, 3, 3> C;
  enu2eciDcm<T>(C, lla, dt);
  return C;
}

//! --- ENU2ECEFDCM ---
/// @brief      East-North-Up to Earth-Centered-Earth-Fixed direction cosine matrix
/// @param lla  3x1 Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @returns    3x3 ECEF->NED direction cosine matrix
template <typename T = double>
void enu2ecefDcm(
    Eigen::Ref<Eigen::Matrix<T, 3, 3>> C, const Eigen::Ref<const Eigen::Vector<T, 3>> &lla) {
  T sin_phi = std::sin(lla(0));
  T cos_phi = std::cos(lla(0));
  T sin_lam = std::sin(lla(1));
  T cos_lam = std::cos(lla(1));
  // clang-format off
  C << -sin_lam, -cos_lam * sin_phi, cos_lam * cos_phi,
        cos_lam, -sin_lam * sin_phi, sin_lam * cos_phi,
            0.0,            cos_phi,           sin_phi;
  // clang-format on
}
template <typename T = double>
Eigen::Matrix<T, 3, 3> enu2ecefDcm(const Eigen::Ref<const Eigen::Vector<T, 3>> &lla) {
  Eigen::Matrix<T, 3, 3> C;
  enu2ecefDcm<T>(C, lla);
  return C;
}

//! --- ENU2NEDDCM ---
/// @brief      East-North-Up to North-East-Down direction cosine matrix
/// @returns    3x3 ECEF->ENU direction cosine matrix
template <typename T = double>
void enu2nedDcm(Eigen::Ref<Eigen::Matrix<T, 3, 3>> C) {
  C << 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0;
}
template <typename T = double>
Eigen::Matrix<T, 3, 3> enu2nedDcm() {
  Eigen::Matrix<T, 3, 3> C;
  enu2nedDcm<T>(C);
  return C;
}

//* ===== Position Transformations ============================================================= *//

//! --- LLA2ECI ---
/// @brief      Latitude-Longitude-Height to Earth-Centered-Inertial position coordinates
/// @param eci  3x1 ECI position [m]
/// @param lla  3x1 Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @param dt   time elapsed between frames [s]
/// @returns    ECI position
template <typename T = double>
void lla2eci(
    Eigen::Ref<Eigen::Vector<T, 3>> eci,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla,
    const T &dt) {
  Eigen::Vector<T, 3> xyz = lla2ecef<T>(lla);
  Eigen::Matrix<T, 3, 3> C_e_i = ecef2eciDcm<T>(dt);
  eci = C_e_i * xyz;
}
template <typename T = double>
Eigen::Vector<T, 3> lla2eci(const Eigen::Ref<const Eigen::Vector<T, 3>> &lla, const T &dt) {
  Eigen::Vector<T, 3> eci;
  lla2eci<T>(eci, lla, dt);
  return eci;
}

//! --- LLA2ECEF ---
/// @brief      Latitude-Longitude-Height to Earth-Centered-Earth-Fixed position coordinates
/// @param xyz  3x1 ECEF position [m]
/// @param lla  3x1 Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @returns    ECEF position
template <typename T = double>
void lla2ecef(
    Eigen::Ref<Eigen::Vector<T, 3>> xyz, const Eigen::Ref<const Eigen::Vector<T, 3>> &lla) {
  T sin_phi = std::sin(lla(0));
  T cos_phi = std::cos(lla(0));
  T sin_lam = std::sin(lla(1));
  T cos_lam = std::cos(lla(1));
  T h = lla(2);

  T Re = WGS84_R0<T> / std::sqrt(1.0 - WGS84_E2<T> * sin_phi * sin_phi);
  xyz(0) = (Re + h) * cos_phi * cos_lam;
  xyz(1) = (Re + h) * cos_phi * sin_lam;
  xyz(2) = (Re * (1.0 - WGS84_E2<T>)+h) * sin_phi;
}
template <typename T = double>
Eigen::Vector<T, 3> lla2ecef(const Eigen::Ref<const Eigen::Vector<T, 3>> &lla) {
  Eigen::Vector<T, 3> xyz;
  lla2ecef<T>(xyz, lla);
  return xyz;
}

//! --- LLA2NED ---
/// @brief      Latitude-Longitude-Height to North-East-Down position coordinates
/// @param ned  3x1 NED position [m]
/// @param lla  3x1 Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @param lla0 3x1 Reference Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @returns    NED position
template <typename T = double>
void lla2ned(
    Eigen::Ref<Eigen::Vector<T, 3>> ned,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Matrix<T, 3, 3> C_e_n = ecef2nedDcm<T>(lla0);
  Eigen::Vector<T, 3> xyz0 = lla2ecef<T>(lla0);
  Eigen::Vector<T, 3> xyz = lla2ecef<T>(lla);
  ned = C_e_n * (xyz - xyz0);
}
template <typename T = double>
Eigen::Vector<T, 3> lla2ned(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Vector<T, 3> ned;
  lla2ned<T>(ned, lla, lla0);
  return ned;
}

//! --- LLA2ENU ---
/// @brief      Latitude-Longitude-Height to East-North-Up position coordinates
/// @param enu  3x1 ENU position [m]
/// @param lla  3x1 Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @param lla0 3x1 Reference Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @returns    ENU position
template <typename T = double>
void lla2enu(
    Eigen::Ref<Eigen::Vector<T, 3>> enu,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Matrix<T, 3, 3> C_e_n = ecef2enuDcm<T>(lla0);
  Eigen::Vector<T, 3> xyz0 = lla2ecef<T>(lla0);
  Eigen::Vector<T, 3> xyz = lla2ecef<T>(lla);
  enu = C_e_n * (xyz - xyz0);
}
template <typename T = double>
Eigen::Vector<T, 3> lla2enu(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Vector<T, 3> enu;
  lla2enu<T>(enu, lla, lla0);
  return enu;
}

//! --- LLA2AER ---
/// @brief      Latitude-Longitude-Height to Azimuth-Elevation-Range position coordinates
/// @param aer  3x1 AER position [m]
/// @param llaR 3x1 Reference Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @param llaT 3x1 Target Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @returns    AER from reference to target
template <typename T = double>
void lla2aer(
    Eigen::Ref<Eigen::Vector<T, 3>> aer,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &llaR,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &llaT) {
  Eigen::Vector<T, 3> enu = lla2enu<T>(llaT, llaR);
  aer(2) = enu.norm();
  aer(1) = std::asin(enu(2) / aer(2));
  aer(0) = std::atan2(enu(0), enu(1));
}
template <typename T = double>
Eigen::Vector<T, 3> lla2aer(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &llaR,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &llaT) {
  Eigen::Vector<T, 3> aer;
  lla2aer<T>(aer, llaR, llaT);
  return aer;
}

//! --- ECI2ECEF ---
/// @brief      Earth-Centered-Inertial to Earth-Centered-Earth-Fixed position coordinates
/// @param eci  3x1 ECI position [m]
/// @param xyz  3x1 ECEF position [m]
/// @param dt   time elapsed between frames [s]
/// @returns    ECEF position
template <typename T = double>
void eci2ecef(
    Eigen::Ref<Eigen::Vector<T, 3>> xyz,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &eci,
    const T &dt) {
  Eigen::Matrix<T, 3, 3> C_i_e = eci2ecefDcm<T>(dt);
  xyz = C_i_e * eci;
}
template <typename T = double>
Eigen::Vector<T, 3> eci2ecef(const Eigen::Ref<const Eigen::Vector<T, 3>> &eci, const T &dt) {
  Eigen::Vector<T, 3> xyz;
  eci2ecef<T>(xyz, eci, dt);
  return xyz;
}

//! --- ECI2LLA ---
/// @brief      Earth-Centered-Inertial to Latitude-Longitude-Height position coordinates
/// @param eci  3x1 ECI position [m]
/// @param lla  3x1 LLA position [rad, rad, m]
/// @param dt   time elapsed between frames [s]
/// @returns    lla position
template <typename T = double>
void eci2lla(
    Eigen::Ref<Eigen::Vector<T, 3>> lla,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &eci,
    const T &dt) {
  Eigen::Vector<T, 3> xyz = eci2ecef<T>(eci, dt);
  lla = ecef2lla<T>(xyz);
}
template <typename T = double>
Eigen::Vector<T, 3> eci2lla(const Eigen::Ref<const Eigen::Vector<T, 3>> &eci, const T &dt) {
  Eigen::Vector<T, 3> lla;
  eci2lla<T>(lla, eci, dt);
  return lla;
}

//! --- ECI2NED ---
/// @brief      Earth-Centered-Inertial to North-East-Down position coordinates
/// @param eci  3x1 ECI position [m]
/// @param lla0 3x1 Reference Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @param ned  3x1 NED position [m]
/// @param dt   time elapsed between frames [s]
/// @returns    lla position
template <typename T = double>
void eci2ned(
    Eigen::Ref<Eigen::Vector<T, 3>> ned,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &eci,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Vector<T, 3> xyz = eci2ecef<T>(eci, dt);
  ned = ecef2ned<T>(xyz, lla0);
}
template <typename T = double>
Eigen::Vector<T, 3> eci2ned(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &eci,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Vector<T, 3> ned;
  eci2ned<T>(ned, eci, lla0, dt);
  return ned;
}

//! --- ECI2ENU ---
/// @brief      Earth-Centered-Inertial to East-North-Up position coordinates
/// @param eci  3x1 ECI position [m]
/// @param lla0 3x1 Reference Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @param enu  3x1 ENU position [m]
/// @param dt   time elapsed between frames [s]
/// @returns    lla position
template <typename T = double>
void eci2enu(
    Eigen::Ref<Eigen::Vector<T, 3>> enu,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &eci,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Vector<T, 3> xyz = eci2ecef<T>(eci, dt);
  enu = ecef2enu<T>(xyz, lla0);
}
template <typename T = double>
Eigen::Vector<T, 3> eci2enu(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &eci,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Vector<T, 3> enu;
  eci2enu<T>(enu, eci, lla0, dt);
  return enu;
}

//! --- ECI2AER ---
/// @brief      Earth-Centered-Inertial to Azimuth-Elevation-Range position coordinates
/// @param aer  3x1 AER position [m]
/// @param eciR 3x1 Reference ECI position [m]
/// @param eciT 3x1 Target ECI position [m]
/// @param dt   time elapsed between frames [s]
/// @returns    AER from reference to target
template <typename T = double>
void eci2aer(
    Eigen::Ref<Eigen::Vector<T, 3>> aer,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &eciR,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &eciT,
    const T &dt) {
  aer = ecef2aer<T>(eci2ecef<T>(eciT, dt), eci2ecef<T>(eciR, dt));
}
template <typename T = double>
Eigen::Vector<T, 3> eci2aer(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &eciR,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &eciT,
    const T &dt) {
  Eigen::Vector<T, 3> aer;
  eci2aer<T>(aer, eciR, eciT, dt);
  return aer;
}

//! --- ECEF2ECI ---
/// @brief      Earth-Centered-Earth-Fixed to Earth-Centered-Inertial position coordinates
/// @param xyz  3x1 ECEF position [m]
/// @param eci  3x1 ECI position [m]
/// @param dt   time elapsed between frames [s]
/// @returns    ECEF position
template <typename T = double>
void ecef2eci(
    Eigen::Ref<Eigen::Vector<T, 3>> eci,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &xyz,
    const T &dt) {
  Eigen::Matrix<T, 3, 3> C_e_i = ecef2eciDcm<T>(dt);
  eci = C_e_i * xyz;
}
template <typename T = double>
Eigen::Vector<T, 3> ecef2eci(const Eigen::Ref<const Eigen::Vector<T, 3>> &xyz, const T &dt) {
  Eigen::Vector<T, 3> eci;
  ecef2eci<T>(eci, xyz, dt);
  return eci;
}

//! --- ECEF2LLA ---
/// @brief      Earth-Centered-Earth-Fixed to Latitude-Longitude-Height position coordinates
/// @param xyz  3x1 ECEF position [m]
/// @param lla  3x1 LLA position [rad, rad, m]
/// @returns    lla position
template <typename T = double>
void ecef2lla(
    Eigen::Ref<Eigen::Vector<T, 3>> lla, const Eigen::Ref<const Eigen::Vector<T, 3>> &xyz) {
  const T &x = xyz(0);
  const T &y = xyz(1);
  const T &z = xyz(2);

  T sign_z = std::copysign(1.0, z);
  T sqrt_1_e2 = std::sqrt(1.0 - WGS84_E2<T>);

  T beta = std::sqrt(x * x + y * y);  // (Groves C.18)
  T a = sqrt_1_e2 * std::abs(z);
  T b = WGS84_E2<T> * WGS84_R0<T>;
  T E = (a - b) / beta;             // (Groves C.29)
  T F = (a + b) / beta;             // (Groves C.30)
  T P = 4.0 / 3.0 * (E * F + 1.0);  // (Groves C.31)
  T Q = 2.0 * (E * E - F * F);      // (Groves C.32)
  T D = P * P * P + Q * Q;          // (Groves C.33)
  T sqrt_D = std::sqrt(D);
  T V = std::pow(sqrt_D - Q, 1.0 / 3.0) - std::pow(sqrt_D + Q, 1.0 / 3.0);  // (Groves C.34)
  T G = 0.5 * (std::sqrt(E * E + V) + E);                                   // (Groves C.35)
  T t = std::sqrt(G * G + ((F - V * G) / (2.0 * G - E))) - G;               // (Groves C.36)
  lla(0) = sign_z * std::atan((1.0 - t * t) / (2.0 * t * sqrt_1_e2));       // (Groves C.37)
  lla(1) = std::atan2(y, x);
  lla(2) = (beta - WGS84_R0<T> * t) * std::cos(lla(0)) +
           (z - sign_z * WGS84_R0<T> * sqrt_1_e2) * std::sin(lla(0));  // (Groves C.38)
}
template <typename T = double>
Eigen::Vector<T, 3> ecef2lla(const Eigen::Ref<const Eigen::Vector<T, 3>> &xyz) {
  Eigen::Vector<T, 3> lla;
  ecef2lla<T>(lla, xyz);
  return lla;
}

//! --- ECEF2NED ---
/// @brief      Earth-Centered-Earth-Fixed to North-East-Down position coordinates
/// @param xyz  3x1 ECEF position [m]
/// @param ned  3x1 NED position [m]
/// @param lla0 3x1 Reference LLA position [rad, rad, m]
/// @returns    NED position
template <typename T = double>
void ecef2ned(
    Eigen::Ref<Eigen::Vector<T, 3>> ned,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &xyz,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Matrix<T, 3, 3> C_e_n = ecef2nedDcm<T>(lla0);
  Eigen::Vector<T, 3> xyz0 = lla2ecef<T>(lla0);
  ned = C_e_n * (xyz - xyz0);
}
template <typename T = double>
Eigen::Vector<T, 3> ecef2ned(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &xyz,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Vector<T, 3> ned;
  ecef2ned<T>(ned, xyz, lla0);
  return ned;
}

//! --- ECEF2ENU ---
/// @brief      Earth-Centered-Earth-Fixed to East-North-Up position coordinates
/// @param xyz  3x1 ECEF position [m]
/// @param enu  3x1 ENU position [m]
/// @param lla0 3x1 Reference LLA position [rad, rad, m]
/// @returns    ENU position
template <typename T = double>
void ecef2enu(
    Eigen::Ref<Eigen::Vector<T, 3>> enu,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &xyz,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Matrix<T, 3, 3> C_e_n = ecef2enuDcm<T>(lla0);
  Eigen::Vector<T, 3> xyz0 = lla2ecef<T>(lla0);
  enu = C_e_n * (xyz - xyz0);
}
template <typename T = double>
Eigen::Vector<T, 3> ecef2enu(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &xyz,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Vector<T, 3> enu;
  ecef2enu<T>(enu, xyz, lla0);
  return enu;
}

//! --- ECEF2AER ---
/// @brief      Earth-Centered-Earth-Fixed to Azimuth-Elevation-Range position coordinates
/// @param aer  3x1 AER position [rad, rad, m]
/// @param xyzR 3x1 Reference ECEF position [m]
/// @param xyzT 3x1 Target ECEF position [m]
/// @returns    AER position
template <typename T = double>
void ecef2aer(
    Eigen::Ref<Eigen::Vector<T, 3>> aer,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &xyzR,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &xyzT) {
  Eigen::Vector<T, 3> lla0 = ecef2lla<T>(xyzR);
  Eigen::Matrix<T, 3, 3> C_e_n = ecef2enuDcm<T>(lla0);
  Eigen::Vector<T, 3> enu = C_e_n * (xyzT - xyzR);

  aer(2) = enu.norm();
  aer(1) = std::asin(enu(2) / aer(2));
  aer(0) = std::atan2(enu(0), enu(1));
}
template <typename T = double>
Eigen::Vector<T, 3> ecef2aer(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &xyzR,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &xyzT) {
  Eigen::Vector<T, 3> aer;
  ecef2aer<T>(aer, xyzR, xyzT);
  return aer;
}

//! --- NED2ECI ---
/// @brief      North-East-Down to Earth-Centered-Inertial position coordinates
/// @param ned  3x1 NED position [m]
/// @param eci  3x1 ECI position [m]
/// @param lla0 3x1 Reference LLA position [rad, rad, m]
/// @param dt   time elapsed between frames [s]
/// @returns    ECI position
template <typename T = double>
void ned2eci(
    Eigen::Ref<Eigen::Vector<T, 3>> eci,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &ned,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Vector<T, 3> xyz = ned2ecef<T>(ned, lla0);
  Eigen::Matrix<T, 3, 3> C_e_i = ecef2eciDcm<T>(dt);
  eci = C_e_i * xyz;
}
template <typename T = double>
Eigen::Vector<T, 3> ned2eci(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &ned,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Vector<T, 3> eci;
  ned2eci<T>(eci, ned, lla0, dt);
  return eci;
}

//! --- NED2ECEF ---
/// @brief      North-East-Down to Earth-Centered-Earth-Fixed position coordinates
/// @param ned  3x1 NED position [m]
/// @param xyz  3x1 ECEF position [m]
/// @param lla0 3x1 Reference LLA position [rad, rad, m]
/// @returns    ECEF position
template <typename T = double>
void ned2ecef(
    Eigen::Ref<Eigen::Vector<T, 3>> xyz,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &ned,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Matrix<T, 3, 3> C_n_e = ned2ecefDcm<T>(lla0);
  xyz = lla2ecef<T>(lla0) + C_n_e * ned;
}
template <typename T = double>
Eigen::Vector<T, 3> ned2ecef(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &ned,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Vector<T, 3> xyz;
  ned2ecef<T>(xyz, ned, lla0);
  return xyz;
}

//! --- NED2LLA ---
/// @brief      North-East-Down to Latitude-Longitude-Height position coordinates
/// @param ned  3x1 NED position [m]
/// @param lla  3x1 LLA position [m]
/// @param lla0 3x1 Reference LLA position [rad, rad, m]
/// @returns    LLA position
template <typename T = double>
void ned2lla(
    Eigen::Ref<Eigen::Vector<T, 3>> lla,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &ned,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Vector<T, 3> xyz = ned2ecef<T>(ned, lla0);
  lla = ecef2lla<T>(xyz);
}
template <typename T = double>
Eigen::Vector<T, 3> ned2lla(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &ned,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Vector<T, 3> lla;
  ned2ecef<T>(lla, ned, lla0);
  return lla;
}

//! --- NED2AER ---
/// @brief      North-East-Down to Azimuth-Elevation-Range position coordinates
/// @param aer  3x1 AER position [m]
/// @param nedR 3x1 Reference NED position [m]
/// @param nedT 3x1 Target NED position [m]
/// @returns    AER position
template <typename T = double>
void ned2aer(
    Eigen::Ref<Eigen::Vector<T, 3>> aer,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &nedR,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &nedT) {
  Eigen::Vector<T, 3> d_ned = nedT - nedR;
  aer(2) = d_ned.norm();
  aer(1) = std::asin(-d_ned(2) / aer(2));
  aer(0) = std::atan2(d_ned(1), d_ned(0));
}
template <typename T = double>
Eigen::Vector<T, 3> ned2aer(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &nedR,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &nedT) {
  Eigen::Vector<T, 3> aer;
  ned2aer<T>(aer, nedR, nedT);
  return aer;
}

//! --- ENU2ECI ---
/// @brief      East-North-Up to Earth-Centered-Inertial position coordinates
/// @param enu  3x1 ENU position [m]
/// @param eci  3x1 ECI position [m]
/// @param lla0 3x1 Reference LLA position [rad, rad, m]
/// @param dt   time elapsed between frames [s]
/// @returns    ECI position
template <typename T = double>
void enu2eci(
    Eigen::Ref<Eigen::Vector<T, 3>> eci,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &enu,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Vector<T, 3> xyz = enu2ecef<T>(enu, lla0);
  Eigen::Matrix<T, 3, 3> C_e_i = ecef2eciDcm<T>(dt);
  eci = C_e_i * xyz;
}
template <typename T = double>
Eigen::Vector<T, 3> enu2eci(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &enu,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Vector<T, 3> eci;
  enu2eci<T>(eci, enu, lla0, dt);
  return eci;
}

//! --- ENU2ECEF ---
/// @brief      East-North-Up to Earth-Centered-Earth-Fixed position coordinates
/// @param enu  3x1 ENU position [m]
/// @param xyz  3x1 ECEF position [m]
/// @param lla0 3x1 Reference LLA position [rad, rad, m]
/// @returns    ECEF position
template <typename T = double>
void enu2ecef(
    Eigen::Ref<Eigen::Vector<T, 3>> xyz,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &enu,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Matrix<T, 3, 3> C_n_e = enu2ecefDcm<T>(lla0);
  xyz = lla2ecef<T>(lla0) + C_n_e * enu;
}
template <typename T = double>
Eigen::Vector<T, 3> enu2ecef(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &enu,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Vector<T, 3> xyz;
  enu2ecef<T>(xyz, enu, lla0);
  return xyz;
}

//! --- ENU2LLA ---
/// @brief      East-North-Up to Latitude-Longitude-Height position coordinates
/// @param enu  3x1 ENU position [m]
/// @param lla  3x1 LLA position [m]
/// @param lla0 3x1 Reference LLA position [rad, rad, m]
/// @returns    LLA position
template <typename T = double>
void enu2lla(
    Eigen::Ref<Eigen::Vector<T, 3>> lla,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &enu,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Vector<T, 3> xyz = enu2ecef<T>(enu, lla0);
  lla = ecef2lla<T>(xyz);
}
template <typename T = double>
Eigen::Vector<T, 3> enu2lla(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &enu,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Vector<T, 3> lla;
  enu2lla<T>(lla, enu, lla0);
  return lla;
}

//! --- ENU2AER ---
/// @brief      East-North-Up to Azimuth-Elevation-Range position coordinates
/// @param aer  3x1 AER position [m]
/// @param enuR 3x1 Reference NED position [m]
/// @param enuT 3x1 Target NED position [m]
/// @returns    AER position
template <typename T = double>
void enu2aer(
    Eigen::Ref<Eigen::Vector<T, 3>> aer,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &enuR,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &enuT) {
  Eigen::Vector<T, 3> d_enu = enuT - enuR;
  aer(2) = d_enu.norm();
  aer(1) = std::asin(d_enu(2) / aer(2));
  aer(0) = std::atan2(d_enu(0), d_enu(1));
}
template <typename T = double>
Eigen::Vector<T, 3> enu2aer(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &enuR,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &enuT) {
  Eigen::Vector<T, 3> aer;
  enu2aer<T>(aer, enuR, enuT);
  return aer;
}

//* ===== Velocity Transformations ============================================================= *//

//! --- ECI2ECEFV ---
/// @brief      Converts Earth-Centered-Inertial to Earth-Centered-Earth-Fixed velocity
/// @param r_ib_i   3x1 ECI position [m]
/// @param v_ib_i   3x1 ECI velocity [m/s]
/// @param xyz      3x1 ECEF velocity [m/s]
/// @param dt       time elapsed between frames [s]
/// @returns    ECEF velocity
template <typename T = double>
void eci2ecefv(
    Eigen::Ref<Eigen::Vector<T, 3>> xyz,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_ib_i,
    const T &dt) {
  Eigen::Matrix<T, 3, 3> C_i_e = eci2ecefDcm<T>(dt);
  xyz = C_i_e * (v_ib_i - WGS84_OMEGA_SKEW<T> * r_ib_i);
}
template <typename T = double>
Eigen::Vector<T, 3> eci2ecefv(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_ib_i,
    const T &dt) {
  Eigen::Vector<T, 3> xyz;
  eci2ecefv<T>(xyz, r_ib_i, v_ib_i, dt);
  return xyz;
}

//! --- ECI2NEDV ---
/// @brief      Converts Earth-Centered-Inertial to North-East-Down velocity
/// @param r_ib_i   3x1 ECI position [m]
/// @param v_ib_i   3x1 ECI velocity [m/s]
/// @param lla0     3x1 Reference LLA position [rad, rad, m]
/// @param ned      3x1 NED velocity [m/s]
/// @param dt       time elapsed between frames [s]
/// @returns    NED velocity
template <typename T = double>
void eci2nedv(
    Eigen::Ref<Eigen::Vector<T, 3>> ned,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Matrix<T, 3, 3> C_i_n = eci2nedDcm<T>(lla0, dt);
  ned = C_i_n * (v_ib_i - WGS84_OMEGA_SKEW<T> * r_ib_i);
}
template <typename T = double>
Eigen::Vector<T, 3> eci2nedv(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Vector<T, 3> ned;
  eci2nedv<T>(ned, r_ib_i, v_ib_i, lla0, dt);
  return ned;
}

//! --- ECI2ENUV ---
/// @brief      Converts Earth-Centered-Inertial to East-North-Up velocity
/// @param r_ib_i   3x1 ECI position [m]
/// @param v_ib_i   3x1 ECI velocity [m/s]
/// @param lla0     3x1 Reference LLA position [rad, rad, m]
/// @param enu      3x1 ENU velocity [m/s]
/// @param dt       time elapsed between frames [s]
/// @returns    ENU velocity
template <typename T = double>
void eci2enuv(
    Eigen::Ref<Eigen::Vector<T, 3>> enu,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Matrix<T, 3, 3> C_i_n = eci2enuDcm<T>(lla0, dt);
  enu = C_i_n * (v_ib_i - WGS84_OMEGA_SKEW<T> * r_ib_i);
}
template <typename T = double>
Eigen::Vector<T, 3> eci2enuv(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Vector<T, 3> enu;
  eci2enuv<T>(enu, r_ib_i, v_ib_i, lla0, dt);
  return enu;
}

//! --- ECEF2ECIV ---
/// @brief      Converts Earth-Centered-Earth-Fixed to Earth-Centered-Inertial velocity
/// @param r_eb_e   3x1 ECEF position [m]
/// @param v_eb_e   3x1 ECEF velocity [m/s]
/// @param eci      3x1 ECI velocity [m/s]
/// @param dt       time elapsed between frames [s]
/// @returns    ECI velocity
template <typename T = double>
void ecef2eciv(
    Eigen::Ref<Eigen::Vector<T, 3>> eci,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_eb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_eb_e,
    const T &dt) {
  Eigen::Vector<T, 3> C_e_i = ecef2eciDcm<T>(dt);
  eci = C_e_i * (v_eb_e - WGS84_OMEGA_SKEW<T> * r_eb_e);
}
template <typename T = double>
Eigen::Vector<T, 3> ecef2eciv(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_eb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_eb_e,
    const T &dt) {
  Eigen::Vector<T, 3> eci;
  ecef2eciv<T>(eci, r_eb_e, v_eb_e, dt);
  return eci;
}

//! --- ECEF2NEDV ---
/// @brief      Converts Earth-Centered-Earth-Fixed to North-East-Down velocity
/// @param v_eb_e   3x1 ECEF velocity [m/s]
/// @param lla0     3x1 Reference LLA position [rad, rad, m]
/// @param ned      3x1 NED velocity [m/s]
/// @returns    NED velocity
template <typename T = double>
void ecef2nedv(
    Eigen::Ref<Eigen::Vector<T, 3>> ned,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_eb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Matrix<T, 3, 3> C_e_n = ecef2nedDcm<T>(lla0);
  ned = C_e_n * v_eb_e;
}
template <typename T = double>
Eigen::Vector<T, 3> ecef2nedv(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_eb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Vector<T, 3> ned;
  ecef2nedv<T>(ned, v_eb_e, lla0);
  return ned;
}

//! --- ECEF2ENUV ---
/// @brief      Converts Earth-Centered-Earth-Fixed to East-North0Up velocity
/// @param v_eb_e   3x1 ECEF velocity [m/s]
/// @param lla0     3x1 Reference LLA position [rad, rad, m]
/// @param enu      3x1 NED velocity [m/s]
/// @returns    ENU velocity
template <typename T = double>
void ecef2enuv(
    Eigen::Ref<Eigen::Vector<T, 3>> enu,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_eb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Matrix<T, 3, 3> C_e_n = ecef2enuDcm<T>(lla0);
  enu = C_e_n * v_eb_e;
}
template <typename T = double>
Eigen::Vector<T, 3> ecef2enuv(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_eb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Vector<T, 3> enu;
  ecef2enuv<T>(enu, v_eb_e, lla0);
  return enu;
}

//! --- NED2ECIV ---
/// @brief      Converts North-East-Down to Earth-Centered-Inertial velocity
/// @param r_nb_e   3x1 NED position [m]
/// @param v_nb_e   3x1 NED velocity [m/s]
/// @param lla0     3x1 Reference LLA position [rad, rad, m]
/// @param eci      3x1 NED velocity [m/s]
/// @param dt       time elapsed between frames [s]
/// @returns    ECI velocity
template <typename T = double>
void ned2eciv(
    Eigen::Ref<Eigen::Vector<T, 3>> eci,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Matrix<T, 3, 3> C_n_i = ned2eciDcm<T>(lla0, dt);
  Eigen::Matrix<T, 3, 3> C_e_i = ecef2eciDcm<T>(dt);
  Eigen::Vector<T, 3> xyz = ned2ecef<T>(r_nb_e, lla0);
  eci = C_n_i * v_nb_e + C_e_i * WGS84_OMEGA_SKEW<T> * xyz;
}
template <typename T = double>
Eigen::Vector<T, 3> ned2eciv(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Vector<T, 3> eci;
  ned2eciv<T>(eci, r_nb_e, v_nb_e, lla0, dt);
  return eci;
}

//! --- NED2ECEFV ---
/// @brief      Converts North-East-Down to Earth-Centered-Earth-Fixed velocity
/// @param v_nb_e   3x1 NED velocity [m/s]
/// @param lla0     3x1 Reference LLA position [rad, rad, m]
/// @param xyz      3x1 ECEF velocity [m/s]
/// @returns    ENU velocity
template <typename T = double>
void ned2ecefv(
    Eigen::Ref<Eigen::Vector<T, 3>> xyz,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Matrix<T, 3, 3> C_n_e = ned2ecefDcm<T>(lla0);
  xyz = C_n_e * v_nb_e;
}
template <typename T = double>
Eigen::Vector<T, 3> ned2ecefv(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Vector<T, 3> xyz;
  ned2ecefv<T>(xyz, v_nb_e, lla0);
  return xyz;
}

//! --- ENU2ECIV ---
/// @brief      Converts North-East-Down to Earth-Centered-Inertial velocity
/// @param r_nb_e   3x1 ENU position [m]
/// @param v_nb_e   3x1 ENU velocity [m/s]
/// @param lla0     3x1 Reference LLA position [rad, rad, m]
/// @param eci      3x1 ECI velocity [m/s]
/// @param dt       time elapsed between frames [s]
/// @returns    ECI velocity
template <typename T = double>
void enu2eciv(
    Eigen::Ref<Eigen::Vector<T, 3>> eci,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Matrix<T, 3, 3> C_n_i = enu2eciDcm<T>(lla0, dt);
  Eigen::Matrix<T, 3, 3> C_e_i = ecef2eciDcm<T>(dt);
  Eigen::Vector<T, 3> xyz = ned2ecef<T>(r_nb_e, lla0);
  eci = C_n_i * v_nb_e + C_e_i * WGS84_OMEGA_SKEW<T> * xyz;
}
template <typename T = double>
Eigen::Vector<T, 3> enu2eciv(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Vector<T, 3> eci;
  enu2eciv<T>(eci, r_nb_e, v_nb_e, lla0, dt);
  return eci;
}

//! --- ENU2ECEFV ---
/// @brief      Converts East-North-Up to Earth-Centered-Earth-Fixed velocity
/// @param v_nb_e   3x1 ENU velocity [m/s]
/// @param lla0     3x1 Reference LLA position [rad, rad, m]
/// @param xyz      3x1 ECEF velocity [m/s]
/// @returns    ENU velocity
template <typename T = double>
void enu2ecefv(
    Eigen::Ref<Eigen::Vector<T, 3>> xyz,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Matrix<T, 3, 3> C_n_e = enu2ecefDcm<T>(lla0);
  xyz = C_n_e * v_nb_e;
}
template <typename T = double>
Eigen::Vector<T, 3> enu2ecefv(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Vector<T, 3> xyz;
  enu2ecefv<T>(xyz, v_nb_e, lla0);
  return xyz;
}

//* ===== Angular Velocity Transformations ===================================================== *//

//! --- ECI2ECEFW ---
/// @brief      Converts Earth-Centered-Inertial to Earth-Centered-Earth-Fixed angular velocity
/// @param w_ib_i   3x1 ECI angular velocity [rad/s]
/// @param dt       time elapsed between frames [s]
/// @param xyz      3x1 ECEF angular velocity [rad/s]
/// @returns    ECEF angular velocity
template <typename T = double>
void eci2ecefw(
    Eigen::Ref<Eigen::Vector<T, 3>> xyz,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &w_ib_i,
    const T &dt) {
  Eigen::Matrix<T, 3, 3> C_i_e = eci2ecefDcm<T>(dt);
  xyz = C_i_e * (w_ib_i + WGS84_OMEGA_VEC<T>);
}
template <typename T = double>
Eigen::Vector<T, 3> eci2ecefw(const Eigen::Ref<const Eigen::Vector<T, 3>> &w_ib_i, const T &dt) {
  Eigen::Vector<T, 3> xyz;
  eci2ecefw<T>(xyz, w_ib_i, dt);
  return xyz;
}

//! --- ECI2NEDW ---
/// @brief      Converts Earth-Centered-Inertial to North-East-Down angular velocity
/// @param w_ib_i   3x1 ECI angular velocity [rad/s]
/// @param r_ib_i   3x1 ECI position [m]
/// @param v_ib_i   3x1 ECI velocity [m/s]
/// @param lla0     Reference LLA [rad, rad, m]
/// @param dt       time elapsed between frames [s]
/// @param ned      3x1 NED angular velocity [rad/s]
/// @returns    NED angular velocity
template <typename T = double>
void eci2nedw(
    Eigen::Ref<Eigen::Vector<T, 3>> ned,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &w_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Vector<T, 3> C_i_n = eci2nedDcm<T>(lla0, dt);

  Eigen::Vector<T, 3> v_nb_e = eci2nedv(r_ib_i, v_ib_i, lla0, dt);
  T vn = v_nb_e(0);
  T ve = v_nb_e(1);
  T phi = lla0(0);
  T h = lla0(2);
  T sin_phi = std::sin(phi);

  T trans = 1.0 - WGS84_E2<T> * sin_phi * sin_phi;
  T Re = WGS84_R0<T> / std::sqrt(trans);
  T Rn = WGS84_R0<T> * (1.0 - WGS84_E2<T>) / std::pow(trans, 1.5);
  Eigen::Vector<T, 3> w_en_n{ve / (Re + h), -vn / (Rn + h), -ve * std::tan(phi) / (Re + h)};

  ned = -w_en_n + C_i_n * (w_ib_i - WGS84_OMEGA_VEC<T>);
}
template <typename T = double>
Eigen::Vector<T, 3> eci2nedw(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &w_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Vector<T, 3> ned;
  eci2nedw(ned, w_ib_i, r_ib_i, v_ib_i, lla0, dt);
  return ned;
}

//! --- ECI2ENUW ---
/// @brief      Converts Earth-Centered-Inertial to East-North-Up angular velocity
/// @param w_ib_i   3x1 ECI angular velocity [rad/s]
/// @param r_ib_i   3x1 ECI position [m]
/// @param v_ib_i   3x1 ECI velocity [m/s]
/// @param lla0     Reference LLA [rad, rad, m]
/// @param dt       time elapsed between frames [s]
/// @param enu      3x1 ECU angular velocity [rad/s]
/// @returns    ECU angular velocity
template <typename T = double>
void eci2enuw(
    Eigen::Ref<Eigen::Vector<T, 3>> enu,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &w_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Vector<T, 3> C_i_n = eci2enuDcm<T>(lla0, dt);

  Eigen::Vector<T, 3> v_nb_e = eci2enuv(r_ib_i, v_ib_i, lla0, dt);
  T vn = v_nb_e(0);
  T ve = v_nb_e(1);
  T phi = lla0(0);
  T h = lla0(2);
  T sin_phi = std::sin(phi);

  T trans = 1.0 - WGS84_E2<T> * sin_phi * sin_phi;
  T Re = WGS84_R0<T> / std::sqrt(trans);
  T Rn = WGS84_R0<T> * (1.0 - WGS84_E2<T>) / std::pow(trans, 1.5);
  Eigen::Vector<T, 3> w_en_n{-vn / (Rn + h), ve / (Re + h), ve * std::tan(phi) / (Re + h)};

  enu = -w_en_n + C_i_n * (w_ib_i - WGS84_OMEGA_VEC<T>);
}
template <typename T = double>
Eigen::Vector<T, 3> eci2enuw(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &w_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Vector<T, 3> enu;
  eci2enuw(enu, w_ib_i, r_ib_i, v_ib_i, lla0, dt);
  return enu;
}

//! --- ECEF2ECIW ---
/// @brief      Converts Earth-Centered-Earth-Fixed to Earth-Centered-Inertial angular velocity
/// @param w_eb_e   3x1 ECEF angular velocity [rad/s]
/// @param dt       time elapsed between frames [s]
/// @param eci      3x1 ECI angular velocity [rad/s]
/// @returns    ECI angular velocity
template <typename T = double>
void ecef2eciw(
    Eigen::Ref<Eigen::Vector<T, 3>> eci,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &w_eb_e,
    const T &dt) {
  Eigen::Matrix<T, 3, 3> C_e_i = ecef2eciDcm<T>(dt);
  eci = C_e_i * (w_eb_e + WGS84_OMEGA_VEC<T>);
}
template <typename T = double>
Eigen::Vector<T, 3> ecef2eciw(const Eigen::Ref<const Eigen::Vector<T, 3>> &w_eb_e, const T &dt) {
  Eigen::Vector<T, 3> eci;
  ecef2eciw(eci, w_eb_e, dt);
  return eci;
}

//! --- ECEF2NEDW ---
/// @brief      Converts Earth-Centered-Earth-Fixed to North-East-Down angular velocity
/// @param w_eb_e   3x1 ECEF angular velocity [rad/s]
/// @param lla0     Reference LLA [rad, rad, m]
/// @param ned      3x1 NED angular velocity [rad/s]
/// @returns    NED angular velocity
template <typename T = double>
void ecef2nedw(
    Eigen::Ref<Eigen::Vector<T, 3>> ned,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &w_eb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Matrix<T, 3, 3> C_e_n = ecef2nedDcm<T>(lla0);
  ned = C_e_n * w_eb_e;
}
template <typename T = double>
Eigen::Vector<T, 3> ecef2nedw(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &w_eb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Vector<T, 3> ned;
  ecef2nedw(ned, w_eb_e, lla0);
  return ned;
}

//! --- ECEF2ENUW ---
/// @brief      Converts Earth-Centered-Earth-Fixed to East-North-Up angular velocity
/// @param w_eb_e   3x1 ECEF angular velocity [rad/s]
/// @param lla0     Reference LLA [rad, rad, m]
/// @param ned      3x1 ECU angular velocity [rad/s]
/// @returns    ENU angular velocity
template <typename T = double>
void ecef2enuw(
    Eigen::Ref<Eigen::Vector<T, 3>> enu,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &w_eb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Matrix<T, 3, 3> C_e_n = ecef2enuDcm<T>(lla0);
  enu = C_e_n * w_eb_e;
}
template <typename T = double>
Eigen::Vector<T, 3> ecef2enuw(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &w_eb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Vector<T, 3> enu;
  ecef2nedw(enu, w_eb_e, lla0);
  return enu;
}

//! --- NED2ECIW ---
/// @brief      Converts North-East-Down to Earth-Centered-Inertial angular velocity
/// @param w_nb_e   3x1 NED angular velocity [rad/s]
/// @param v_nb_e   3x1 NED velocity [m/s]
/// @param lla0     Reference LLA [rad, rad, m]
/// @param dt       time elapsed between frames [s]
/// @param eci      3x1 ECI angular velocity [rad/s]
/// @returns    ECI angular velocity
template <typename T = double>
void ned2eciw(
    Eigen::Ref<Eigen::Vector<T, 3>> eci,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &w_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Vector<T, 3> C_n_i = ned2eciDcm<T>(lla0, dt);

  T vn = v_nb_e(0);
  T ve = v_nb_e(1);
  T phi = lla0(0);
  T h = lla0(2);
  T sin_phi = std::sin(phi);

  T trans = 1.0 - WGS84_E2<T> * sin_phi * sin_phi;
  T Re = WGS84_R0<T> / std::sqrt(trans);
  T Rn = WGS84_R0<T> * (1.0 - WGS84_E2<T>) / std::pow(trans, 1.5);
  Eigen::Vector<T, 3> w_en_n{ve / (Re + h), -vn / (Rn + h), -ve * std::tan(phi) / (Re + h)};

  eci = C_n_i * (w_nb_e + w_en_n) + WGS84_OMEGA_VEC<T>;
}
template <typename T = double>
Eigen::Vector<T, 3> ned2eciw(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &w_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Vector<T, 3> eci;
  ned2eciw(eci, w_nb_e, v_nb_e, lla0, dt);
  return eci;
}

//! --- NED2ECEFW ---
/// @brief      Converts North-East-Down to Earth-Centered-Earth-Fixed angular velocity
/// @param w_nb_e   3x1 NED angular velocity [rad/s]
/// @param lla0     Reference LLA [rad, rad, m]
/// @param xyz      3x1 ECEF angular velocity [rad/s]
/// @returns    NED angular velocity
template <typename T = double>
void ned2ecefw(
    Eigen::Ref<Eigen::Vector<T, 3>> xyz,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &w_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Matrix<T, 3, 3> C_n_e = ned2ecefDcm<T>(lla0);
  xyz = C_n_e * w_nb_e;
}
template <typename T = double>
Eigen::Vector<T, 3> ned2ecefw(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &w_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Vector<T, 3> xyz;
  ned2ecefw(xyz, w_nb_e, lla0);
  return xyz;
}

//! --- ENU2ECIW ---
/// @brief      Converts East-North-Up to Earth-Centered-Inertial angular velocity
/// @param w_nb_e   3x1 ENU angular velocity [rad/s]
/// @param v_nb_e   3x1 ENU velocity [m/s]
/// @param lla0     Reference LLA [rad, rad, m]
/// @param dt       time elapsed between frames [s]
/// @param eci      3x1 ECI angular velocity [rad/s]
/// @returns    ECI angular velocity
template <typename T = double>
void enu2eciw(
    Eigen::Ref<Eigen::Vector<T, 3>> eci,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &w_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Vector<T, 3> C_n_i = ned2eciDcm<T>(lla0, dt);

  T vn = v_nb_e(0);
  T ve = v_nb_e(1);
  T phi = lla0(0);
  T h = lla0(2);
  T sin_phi = std::sin(phi);

  T trans = 1.0 - WGS84_E2<T> * sin_phi * sin_phi;
  T Re = WGS84_R0<T> / std::sqrt(trans);
  T Rn = WGS84_R0<T> * (1.0 - WGS84_E2<T>) / std::pow(trans, 1.5);
  Eigen::Vector<T, 3> w_en_n{-vn / (Rn + h), ve / (Re + h), ve * std::tan(phi) / (Re + h)};

  eci = C_n_i * (w_nb_e + w_en_n) + WGS84_OMEGA_VEC<T>;
}
template <typename T = double>
Eigen::Vector<T, 3> enu2eciw(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &w_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Vector<T, 3> eci;
  enu2eciw(eci, w_nb_e, v_nb_e, lla0, dt);
  return eci;
}

//! --- ENU2ECEFW ---
/// @brief      Converts North-East-Down to Earth-Centered-Earth-Fixed angular velocity
/// @param w_nb_e   3x1 NED angular velocity [rad/s]
/// @param lla0     Reference LLA [rad, rad, m]
/// @param xyz      3x1 ECEF angular velocity [rad/s]
/// @returns    NED angular velocity
template <typename T = double>
void enu2ecefw(
    Eigen::Ref<Eigen::Vector<T, 3>> xyz,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &w_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Matrix<T, 3, 3> C_n_e = enu2ecefDcm<T>(lla0);
  xyz = C_n_e * w_nb_e;
}
template <typename T = double>
Eigen::Vector<T, 3> enu2ecefw(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &w_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Vector<T, 3> xyz;
  enu2ecefw(xyz, w_nb_e, lla0);
  return xyz;
}

//* ===== Acceleration Transformations ========================================================= *//

//! --- ECI2ECEFA ---
/// @brief      Converts Earth-Centered-Inertial to Earth-Centered-Earth-Fixed acceleration
/// @param a_ib_i   3x1 ECI acceleration [rad/s]
/// @param r_ib_i   3x1 ECI position [m]
/// @param v_ib_i   3x1 ECI velocity [m/s]
/// @param dt       time elapsed between frames [s]
/// @param xyz      3x1 ECEF acceleration [rad/s]
/// @returns    ECEF acceleration
template <typename T = double>
void eci2ecefa(
    Eigen::Ref<Eigen::Vector<T, 3>> xyz,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &a_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_ib_i,
    const T &dt) {
  Eigen::Matrix<T, 3, 3> C_i_e = eci2ecefDcm<T>(dt);
  xyz = C_i_e * (a_ib_i - 2.0 * WGS84_OMEGA_SKEW<T> * v_ib_i +
                 WGS84_OMEGA_SKEW<T> * WGS84_OMEGA_SKEW<T> * r_ib_i);
}
template <typename T = double>
Eigen::Vector<T, 3> eci2ecefa(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &a_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_ib_i,
    const T &dt) {
  Eigen::Vector<T, 3> xyz;
  eci2ecefa(xyz, a_ib_i, r_ib_i, v_ib_i, dt);
  return xyz;
}

//! --- ECI2NEDA ---
/// @brief      Converts Earth-Centered-Inertial to North-East-Down acceleration
/// @param a_ib_i   3x1 ECI acceleration [rad/s]
/// @param r_ib_i   3x1 ECI position [m]
/// @param v_ib_i   3x1 ECI velocity [m/s]
/// @param lla0     Reference LLA [rad, rad, m]
/// @param dt       time elapsed between frames [s]
/// @param ned      3x1 NED acceleration [rad/s]
/// @returns    NED angular velocity
template <typename T = double>
void eci2neda(
    Eigen::Ref<Eigen::Vector<T, 3>> ned,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &a_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Matrix<T, 3, 3> C_i_n = eci2nedDcm<T>(lla0, dt);
  ned = C_i_n * (a_ib_i + 2.0 * WGS84_OMEGA_SKEW<T> * v_ib_i +
                 WGS84_OMEGA_SKEW<T> * WGS84_OMEGA_SKEW<T> * r_ib_i);
}
template <typename T = double>
Eigen::Vector<T, 3> eci2neda(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &a_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Vector<T, 3> ned;
  eci2enua(ned, a_ib_i, r_ib_i, v_ib_i, lla0, dt);
  return ned;
}

//! --- ECI2ENUA ---
/// @brief      Converts Earth-Centered-Inertial to East-North-Up acceleration
/// @param a_ib_i   3x1 ECI acceleration [rad/s]
/// @param r_ib_i   3x1 ECI position [m]
/// @param v_ib_i   3x1 ECI velocity [m/s]
/// @param lla0     Reference LLA [rad, rad, m]
/// @param dt       time elapsed between frames [s]
/// @param enu      3x1 ECU acceleration [rad/s]
/// @returns    ECU acceleration
template <typename T = double>
void eci2enua(
    Eigen::Ref<Eigen::Vector<T, 3>> enu,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &a_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Matrix<T, 3, 3> C_i_n = eci2enuDcm<T>(lla0, dt);
  enu = C_i_n * (a_ib_i + 2.0 * WGS84_OMEGA_SKEW<T> * v_ib_i +
                 WGS84_OMEGA_SKEW<T> * WGS84_OMEGA_SKEW<T> * r_ib_i);
}
template <typename T = double>
Eigen::Vector<T, 3> eci2enua(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &a_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_ib_i,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Vector<T, 3> enu;
  eci2enua<T>(enu, a_ib_i, r_ib_i, v_ib_i, lla0, dt);
  return enu;
}

//! --- ECEF2ECIA ---
/// @brief      Converts Earth-Centered-Earth-Fixed to Earth-Centered-Inertial acceleration
/// @param a_eb_e   3x1 ECEF acceleration [rad/s]
/// @param r_eb_e   3x1 ECEF position [m]
/// @param v_eb_e   3x1 ECEF velocity [m/s]
/// @param dt       time elapsed between frames [s]
/// @param eci      3x1 ECI acceleration [rad/s]
/// @returns    ECI acceleration
template <typename T = double>
void ecef2ecia(
    Eigen::Ref<Eigen::Vector<T, 3>> eci,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &a_eb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_eb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_eb_e,
    const T &dt) {
  Eigen::Matrix<T, 3, 3> C_i_e = eci2ecefDcm<T>(dt);
  eci = C_i_e * (a_eb_e - 2.0 * WGS84_OMEGA_SKEW<T> * v_eb_e +
                 WGS84_OMEGA_SKEW<T> * WGS84_OMEGA_SKEW<T> * r_eb_e);
}
template <typename T = double>
Eigen::Vector<T, 3> ecef2ecia(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &a_eb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_eb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_eb_e,
    const T &dt) {
  Eigen::Vector<T, 3> eci;
  ecef2ecia<T>(eci, a_eb_e, r_eb_e, v_eb_e, dt);
  return eci;
}

//! --- ECEF2NEDA ---
/// @brief      Converts Earth-Centered-Earth-Fixed to North-East-Down acceleration
/// @param a_eb_e   3x1 ECEF acceleration [rad/s]
/// @param lla0     Reference LLA [rad, rad, m]
/// @param ned      3x1 NED acceleration [rad/s]
/// @returns    ENU acceleration
template <typename T = double>
void ecef2neda(
    Eigen::Ref<Eigen::Vector<T, 3>> ned,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &a_eb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Matrix<T, 3, 3> C_e_n = ecef2nedDcm<T>(lla0);
  ned = C_e_n * a_eb_e;
}
template <typename T = double>
Eigen::Vector<T, 3> ecef2neda(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &a_eb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Vector<T, 3> ned;
  ecef2neda<T>(ned, a_eb_e, lla0);
  return ned;
}

//! --- ECEF2ENUA ---
/// @brief      Converts Earth-Centered-Earth-Fixed to East-North-Up acceleration
/// @param a_eb_e   3x1 ECEF acceleration [rad/s]
/// @param lla0     Reference LLA [rad, rad, m]
/// @param enu      3x1 ENU acceleration [rad/s]
/// @returns    ENU acceleration
template <typename T = double>
void ecef2enua(
    Eigen::Ref<Eigen::Vector<T, 3>> enu,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &a_eb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Matrix<T, 3, 3> C_e_n = ecef2enuDcm<T>(lla0);
  enu = C_e_n * a_eb_e;
}
template <typename T = double>
Eigen::Vector<T, 3> ecef2enua(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &a_eb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Vector<T, 3> enu;
  ecef2enua<T>(enu, a_eb_e, lla0);
  return enu;
}

//! --- NED2ECIA ---
/// @brief      Converts North-East-Down to Earth-Centered-Inertial acceleration
/// @param a_nb_e   3x1 NED acceleration [rad/s]
/// @param r_nb_e   3x1 NED position [m]
/// @param v_nb_e   3x1 NED velocity [m/s]
/// @param lla0     Reference LLA [rad, rad, m]
/// @param dt       time elapsed between frames [s]
/// @param eci      3x1 ECI acceleration [rad/s]
/// @returns    ECI acceleration
template <typename T = double>
void ned2ecia(
    Eigen::Ref<Eigen::Vector<T, 3>> eci,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &a_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Vector<T, 3> r_eb_e = ned2ecef<T>(r_nb_e, lla0);
  Eigen::Matrix<T, 3, 3> C_n_e = ned2ecefDcm<T>(lla0);
  T omega_sin_phi = WGS84_OMEGA<T> * std::sin(lla0(0));
  T omega_cos_phi = WGS84_OMEGA<T> * std::cos(lla0(0));
  Eigen::Matrix<T, 3, 3> omega_ie_n{
      {0.0, omega_cos_phi, 0.0}, {-omega_cos_phi, 0.0, -omega_sin_phi}, {0.0, omega_sin_phi, 0.0}};
  Eigen::Matrix<T, 3, 3> C_n_i = ned2eciDcm<T>(lla0, dt);
  Eigen::Matrix<T, 3, 3> C_e_i = ecef2eciDcm<T>(dt);

  Eigen::Vector<T, 3> a_eb_e = C_n_e * a_nb_e;
  Eigen::Vector<T, 3> v_eb_e = C_n_e * v_nb_e;
  eci = C_n_i * (a_eb_e + 2.0 * omega_ie_n * v_eb_e +
                 C_e_i * (WGS84_OMEGA_SKEW<T> * WGS84_OMEGA_SKEW<T>, *r_eb_e));
}
template <typename T = double>
Eigen::Vector<T, 3> ned2ecia(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &a_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Vector<T, 3> eci;
  ned2ecia<T>(eci, a_nb_e, r_nb_e, v_nb_e, lla0, dt);
  return eci;
}

//! --- NED2ECEFA ---
/// @brief      Converts North-East-Down to Earth-Centered-Earth-Fixed acceleration
/// @param a_nb_e   3x1 NED acceleration [rad/s]
/// @param lla0     Reference LLA [rad, rad, m]
/// @param xyz      3x1 ECEF acceleration [rad/s]
/// @returns    ECEF acceleration
template <typename T = double>
void ned2ecefa(
    Eigen::Ref<Eigen::Vector<T, 3>> xyz,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &a_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Matrix<T, 3, 3> C_n_e = ned2ecefDcm<T>(lla0);
  xyz = C_n_e * a_nb_e;
}
template <typename T = double>
Eigen::Vector<T, 3> ned2ecefa(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &a_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Vector<T, 3> xyz;
  ned2ecefa<T>(xyz, a_nb_e, lla0);
  return xyz;
}

//! --- ENU2ECIA ---
/// @brief      Converts East-North-Up to Earth-Centered-Inertial acceleration
/// @param a_nb_e   3x1 ENU acceleration [rad/s]
/// @param r_nb_e   3x1 ENU position [m]
/// @param v_nb_e   3x1 ENU velocity [m/s]
/// @param lla0     Reference LLA [rad, rad, m]
/// @param dt       time elapsed between frames [s]
/// @param eci      3x1 ECI acceleration [rad/s]
/// @returns    ECI acceleration
template <typename T = double>
void enu2ecia(
    Eigen::Ref<Eigen::Vector<T, 3>> eci,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &a_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Vector<T, 3> r_eb_e = enu2ecef<T>(r_nb_e, lla0);
  Eigen::Matrix<T, 3, 3> C_n_e = enu2ecefDcm<T>(lla0);
  T omega_sin_phi = WGS84_OMEGA<T> * std::sin(lla0(0));
  T omega_cos_phi = WGS84_OMEGA<T> * std::cos(lla0(0));
  Eigen::Matrix<T, 3, 3> omega_ie_n{
      {0.0, 0.0, omega_cos_phi}, {0.0, 0.0, omega_sin_phi}, {-omega_cos_phi, -omega_sin_phi, 0.0}};
  Eigen::Matrix<T, 3, 3> C_n_i = enu2eciDcm<T>(lla0, dt);
  Eigen::Matrix<T, 3, 3> C_e_i = ecef2eciDcm<T>(dt);

  Eigen::Vector<T, 3> a_eb_e = C_n_e * a_nb_e;
  Eigen::Vector<T, 3> v_eb_e = C_n_e * v_nb_e;
  eci = C_n_i * (a_eb_e + 2.0 * omega_ie_n * v_eb_e +
                 C_e_i * (WGS84_OMEGA_SKEW<T> * WGS84_OMEGA_SKEW<T>, r_eb_e));
}
template <typename T = double>
Eigen::Vector<T, 3> enu2ecia(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &a_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &r_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &v_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0,
    const T &dt) {
  Eigen::Vector<T, 3> eci;
  enu2ecia<T>(eci, a_nb_e, r_nb_e, v_nb_e, lla0, dt);
  return eci;
}

//! --- ENU2ECEFA ---
/// @brief      Converts East-North-Up to Earth-Centered-Earth-Fixed acceleration
/// @param a_nb_e   3x1 ENU acceleration [rad/s]
/// @param lla0     Reference LLA [rad, rad, m]
/// @param xyz      3x1 ECEF acceleration [rad/s]
/// @returns    ECEF acceleration
template <typename T = double>
void enu2ecefa(
    Eigen::Ref<Eigen::Vector<T, 3>> xyz,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &a_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Matrix<T, 3, 3> C_n_e = enu2ecefDcm<T>(lla0);
  xyz = C_n_e * a_nb_e;
}
template <typename T = double>
Eigen::Vector<T, 3> enu2ecefa(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &a_nb_e,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &lla0) {
  Eigen::Vector<T, 3> xyz;
  enu2ecefa<T>(xyz, a_nb_e, lla0);
  return xyz;
}

//! --- VECEXP ---
/// @brief

//! --- VECLOG ---

}  // namespace navtools

#endif
