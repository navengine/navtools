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

#include "navtools/constants.hpp"
#include "navtools/types.hpp"

namespace navtools {

//* ===== Direction Cosine Matrices ============================================================ *//

//! --- ECI2ECEFDCM ---
/// @brief      Earth-Centered-Inertial to Earth-Centered-Earth-Fixed direction cosine matrix
/// @param dt   time elapsed between frames [s]
/// @returns    3x3 ECI->ECEF direction cosine matrix
template <typename Float = double>
void eci2ecefDcm(Mat3x3<Float> &C, const Float &dt) {
    Float omega_dt = WGS84_OMEGA<Float> * dt;
    Float sin_omega_dt = std::sin(omega_dt);
    Float cos_omega_dt = std::cos(omega_dt);
    // clang-format off
    C <<  cos_omega_dt, sin_omega_dt, 0.0, 
         -sin_omega_dt, cos_omega_dt, 0.0, 
                   0.0,          0.0, 1.0;
    // clang-format on
}
template <typename Float = double>
Mat3x3<Float> eci2ecefDcm(const Float &dt) {
    Mat3x3<Float> C;
    eci2ecefDcm(C, dt);
    return C;
}

//! --- ECI2NEDDCM ---
/// @brief      Earth-Centered-Inertial to North-East-Down direction cosine matrix
/// @param lla  3x1 Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @param dt   time elapsed between frames [s]
/// @returns    3x3 ECI->NED direction cosine matrix
template <typename Float = double>
void eci2nedDcm(Mat3x3<Float> &C, const Vec3<Float> &lla, const Float &dt) {
    Float omega_dt = WGS84_OMEGA<Float> * dt;
    Float sin_phi = std::sin(lla(0));
    Float cos_phi = std::cos(lla(0));
    Float sin_lam_omega_dt = std::sin(lla(1) + omega_dt);
    Float cos_lam_omega_dt = std::cos(lla(1) + omega_dt);
    // clang-format off
    C << -sin_phi * cos_lam_omega_dt, -sin_phi * sin_lam_omega_dt,  cos_phi,
                   -sin_lam_omega_dt,            cos_lam_omega_dt,      0.0,
         -cos_phi * cos_lam_omega_dt, -cos_phi * sin_lam_omega_dt, -sin_phi;
    // clang-format on
}
template <typename Float = double>
Mat3x3<Float> eci2nedDcm(const Vec3<Float> &lla, const Float &dt) {
    Mat3x3<Float> C;
    eci2nedDcm(C, lla, dt);
    return C;
}

//! --- ECI2ENUDCM ---
/// @brief      Earth-Centered-Inertial to East-North-Up direction cosine matrix
/// @param lla  3x1 Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @param dt   time elapsed between frames [s]
/// @returns    3x3 ECI->ENU direction cosine matrix
template <typename Float = double>
void eci2enuDcm(Mat3x3<Float> &C, const Vec3<Float> &lla, const Float &dt) {
    Float omega_dt = WGS84_OMEGA<Float> * dt;
    Float sin_phi = std::sin(lla(0));
    Float cos_phi = std::cos(lla(0));
    Float sin_lam_omega_dt = std::sin(lla(1) + omega_dt);
    Float cos_lam_omega_dt = std::cos(lla(1) + omega_dt);
    // clang-format off
    C <<           -sin_lam_omega_dt,            cos_lam_omega_dt,     0.0, 
         -sin_phi * cos_lam_omega_dt, -sin_phi * sin_lam_omega_dt, cos_phi, 
          cos_phi * cos_lam_omega_dt,  cos_phi * sin_lam_omega_dt, sin_phi;
            // clang-format off
}
template <typename Float = double>
Mat3x3<Float> eci2enuDcm(const Vec3<Float> &lla, const Float &dt) {
    Mat3x3<Float> C;
    eci2enuDcm(C, lla, dt);
    return C;
}

//! --- ECEF2ECIDCM ---
/// @brief      Earth-Centered-Earth-Fixed to Earth-Centered-Inertial direction cosine matrix
/// @param dt   time elapsed between frames [s]
/// @returns    3x3 ECEF->ECI direction cosine matrix
template <typename Float = double>
void ecef2eciDcm(Mat3x3<Float> &C, const Float &dt) {
    Float omega_dt = WGS84_OMEGA<Float> * dt;
    Float sin_omega_dt = std::sin(omega_dt);
    Float cos_omega_dt = std::cos(omega_dt);
    // clang-format off
    C << cos_omega_dt, -sin_omega_dt, 0.0, 
         sin_omega_dt,  cos_omega_dt, 0.0, 
                  0.0,          0.0,  1.0;
    // clang-format on
}
template <typename Float = double>
Mat3x3<Float> ecef2eciDcm(const Float &dt) {
    Mat3x3<Float> C;
    ecef2eciDcm(C, dt);
    return C;
}

//! --- ECEF2NEDDCM ---
/// @brief      Earth-Centered-Earth-Fixed to North-East-Down direction cosine matrix
/// @param lla  3x1 Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @returns    3x3 ECEF->NED direction cosine matrix
template <typename Float = double>
void ecef2nedDcm(Mat3x3<Float> &C, const Vec3<Float> &lla) {
    Float sin_phi = std::sin(lla(0));
    Float cos_phi = std::cos(lla(0));
    Float sin_lam = std::sin(lla(1));
    Float cos_lam = std::cos(lla(1));
    // clang-format off
    C << -sin_phi * cos_lam, -sin_phi * sin_lam,  cos_phi,
                   -sin_lam,            cos_lam,      0.0,
         -cos_phi * cos_lam, -cos_phi * sin_lam, -sin_phi;
    // clang-format on
}
template <typename Float = double>
Mat3x3<Float> ecef2nedDcm(const Vec3<Float> &lla) {
    Mat3x3<Float> C;
    ecef2nedDcm(C, lla);
    return C;
}

//! --- ECEF2ENUDCM ---
/// @brief      Earth-Centered-Earth-Fixed to East-North-Up direction cosine matrix
/// @param lla  3x1 Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @returns    3x3 ECEF->ENU direction cosine matrix
template <typename Float = double>
void ecef2enuDcm(Mat3x3<Float> &C, const Vec3<Float> &lla) {
    Float sin_phi = std::sin(lla(0));
    Float cos_phi = std::cos(lla(0));
    Float sin_lam = std::sin(lla(1));
    Float cos_lam = std::cos(lla(1));
    // clang-format off
    C <<           -sin_lam,            cos_lam,     0.0,
         -sin_phi * cos_lam, -sin_phi * sin_lam, cos_phi,
          cos_phi * cos_lam,  cos_phi * sin_lam, sin_phi;
    // clang-format on
}
template <typename Float = double>
Mat3x3<Float> ecef2enuDcm(const Vec3<Float> &lla) {
    Mat3x3<Float> C;
    ecef2enuDcm(C, lla);
    return C;
}

//! --- NED2ECIDCM ---
/// @brief      North-East-Down to Earth-Centered-Inertial direction cosine matrix
/// @param lla  3x1 Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @param dt   time elapsed between frames [s]
/// @returns    3x3 ECEF->ECI direction cosine matrix
template <typename Float = double>
void ned2eciDcm(Mat3x3<Float> &C, const Vec3<Float> &lla, const Float &dt) {
    Float omega_dt = WGS84_OMEGA<Float> * dt;
    Float sin_phi = std::sin(lla(0));
    Float cos_phi = std::cos(lla(0));
    Float sin_lam_omega_dt = std::sin(lla(1) + omega_dt);
    Float cos_lam_omega_dt = std::cos(lla(1) + omega_dt);
    // clang-format off
    C << -sin_phi * cos_lam_omega_dt, -sin_lam_omega_dt, -cos_phi * cos_lam_omega_dt,
         -sin_phi * sin_lam_omega_dt,  cos_lam_omega_dt, -cos_phi * sin_lam_omega_dt,
                             cos_phi,               0.0,                    -sin_phi;
    // clang-format on
}
template <typename Float = double>
Mat3x3<Float> ned2eciDcm(const Vec3<Float> &lla, const Float &dt) {
    Mat3x3<Float> C;
    ned2eciDcm(C, lla, dt);
    return C;
}

//! --- NED2ECEFDCM ---
/// @brief      North-East-Down to Earth-Centered-Earth-Fixed direction cosine matrix
/// @param lla  3x1 Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @returns    3x3 ECEF->NED direction cosine matrix
template <typename Float = double>
void ned2ecefDcm(Mat3x3<Float> &C, const Vec3<Float> &lla) {
    Float sin_phi = std::sin(lla(0));
    Float cos_phi = std::cos(lla(0));
    Float sin_lam = std::sin(lla(1));
    Float cos_lam = std::cos(lla(1));
    // clang-format off
    C << -sin_phi * cos_lam, -sin_lam, -cos_phi * cos_lam,
         -sin_phi * sin_lam,  cos_lam, -cos_phi * sin_lam,
                    cos_phi,      0.0,           -sin_phi;
    // clang-format on
}
template <typename Float = double>
Mat3x3<Float> ned2ecefDcm(const Vec3<Float> &lla) {
    Mat3x3<Float> C;
    ned2ecefDcm(C, lla);
    return C;
}

//! --- NED2ENUDCM ---
/// @brief      North-East-Down to East-North-Up direction cosine matrix
/// @returns    3x3 ECEF->ENU direction cosine matrix
template <typename Float = double>
void ned2enuDcm(Mat3x3<Float> &C) {
    C << 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0;
}
template <typename Float = double>
Mat3x3<Float> ned2enuDcm() {
    Mat3x3<Float> C;
    ned2enuDcm(C);
    return C;
}

//! --- ENU2ECIDCM ---
/// @brief      East-North-Up to Earth-Centered-Inertial direction cosine matrix
/// @param lla  3x1 Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @param dt   time elapsed between frames [s]
/// @returns    3x3 ECEF->ECI direction cosine matrix
template <typename Float = double>
void enu2eciDcm(Mat3x3<Float> &C, const Vec3<Float> &lla, const Float &dt) {
    Float omega_dt = WGS84_OMEGA<Float> * dt;
    Float sin_phi = std::sin(lla(0));
    Float cos_phi = std::cos(lla(0));
    Float sin_lam_omega_dt = std::sin(lla(1) + omega_dt);
    Float cos_lam_omega_dt = std::cos(lla(1) + omega_dt);
    // clang-format off
    C << -sin_lam_omega_dt, -sin_phi * cos_lam_omega_dt, cos_phi * cos_lam_omega_dt,
          cos_lam_omega_dt, -sin_phi * sin_lam_omega_dt, cos_phi * sin_lam_omega_dt,
                       0.0,                     cos_phi,                    sin_phi;
    // clang-format on
}
template <typename Float = double>
Mat3x3<Float> enu2eciDcm(const Vec3<Float> &lla, const Float &dt) {
    Mat3x3<Float> C;
    enu2eciDcm(C, lla, dt);
    return C;
}

//! --- ENU2ECEFDCM ---
/// @brief      East-North-Up to Earth-Centered-Earth-Fixed direction cosine matrix
/// @param lla  3x1 Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @returns    3x3 ECEF->NED direction cosine matrix
template <typename Float = double>
void enu2ecefDcm(Mat3x3<Float> &C, const Vec3<Float> &lla) {
    Float sin_phi = std::sin(lla(0));
    Float cos_phi = std::cos(lla(0));
    Float sin_lam = std::sin(lla(1));
    Float cos_lam = std::cos(lla(1));
    // clang-format off
    C << -sin_lam, -cos_lam * sin_phi, cos_lam * cos_phi,
          cos_lam, -sin_lam * sin_phi, sin_lam * cos_phi,
              0.0,            cos_phi,           sin_phi;
    // clang-format on
}
template <typename Float = double>
Mat3x3<Float> enu2ecefDcm(const Vec3<Float> &lla) {
    Mat3x3<Float> C;
    enu2ecefDcm(C, lla);
    return C;
}

//! --- ENU2NEDDCM ---
/// @brief      East-North-Up to North-East-Down direction cosine matrix
/// @returns    3x3 ECEF->ENU direction cosine matrix
template <typename Float = double>
void enu2nedDcm(Mat3x3<Float> &C) {
    C << 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0;
}
template <typename Float = double>
Mat3x3<Float> enu2nedDcm() {
    Mat3x3<Float> C;
    enu2nedDcm(C);
    return C;
}

//* ===== Position Transformations ============================================================= *//

//! --- LLA2ECI ---
/// @brief      Latitude-Longitude-Height to Earth-Centered-Inertial position coordinates
/// @param eci  3x1 ECI position [m]
/// @param lla  3x1 Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @param dt   time elapsed between frames [s]
/// @returns    ECI position
template <typename Float = double>
void lla2eci(Vec3<Float> &eci, const Vec3<Float> &lla, const Float &dt) {
    Vec3<Float> xyz = lla2ecef(lla);
    Mat3x3<Float> C_e_i = ecef2eciDcm(dt);
    eci = C_e_i * xyz;
}
template <typename Float = double>
Vec3<Float> lla2eci(const Vec3<Float> &lla, const Float &dt) {
    Vec3<Float> eci;
    lla2eci(eci, lla, dt);
    return eci;
}

//! --- LLA2ECEF ---
/// @brief      Latitude-Longitude-Height to Earth-Centered-Earth-Fixed position coordinates
/// @param xyz  3x1 ECEF position [m]
/// @param lla  3x1 Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @returns    ECEF position
template <typename Float = double>
void lla2ecef(Vec3<Float> &xyz, const Vec3<Float> &lla) {
    Float sin_phi = std::sin(lla(0));
    Float cos_phi = std::cos(lla(0));
    Float sin_lam = std::sin(lla(1));
    Float cos_lam = std::cos(lla(1));
    Float h = lla(2);

    Float Re = WGS84_R0<Float> / std::sqrt(1.0 - WGS84_E2<Float> * sin_phi * sin_phi);
    xyz(0) = (Re + h) * cos_phi * cos_lam;
    xyz(1) = (Re + h) * cos_phi * sin_lam;
    xyz(2) = (Re * (1.0 - WGS84_E2<Float>)+h) * sin_phi;
}
template <typename Float = double>
Vec3<Float> lla2ecef(const Vec3<Float> &lla) {
    Vec3<Float> xyz;
    lla2ecef(xyz, lla);
    return xyz;
}

//! --- LLA2NED ---
/// @brief      Latitude-Longitude-Height to North-East-Down position coordinates
/// @param ned  3x1 NED position [m]
/// @param lla  3x1 Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @param lla0 3x1 Reference Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @returns    NED position
template <typename Float = double>
void lla2ned(Vec3<Float> &ned, const Vec3<Float> &lla, const Vec3<Float> &lla0) {
    Mat3x3<Float> C_e_n = ecef2nedDcm(lla0);
    Vec3<Float> xyz0 = lla2ecef(lla0);
    Vec3<Float> xyz = lla2ecef(lla);
    ned = C_e_n * (xyz - xyz0);
}
template <typename Float = double>
Vec3<Float> lla2ned(const Vec3<Float> &lla, const Vec3<Float> &lla0) {
    Vec3<Float> ned;
    lla2ned(ned, lla, lla0);
    return ned;
}

//! --- LLA2ENU ---
/// @brief      Latitude-Longitude-Height to East-North-Up position coordinates
/// @param enu  3x1 ENU position [m]
/// @param lla  3x1 Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @param lla0 3x1 Reference Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @returns    ENU position
template <typename Float = double>
void lla2enu(Vec3<Float> &enu, const Vec3<Float> &lla, const Vec3<Float> &lla0) {
    Mat3x3<Float> C_e_n = ecef2enuDcm(lla0);
    Vec3<Float> xyz0 = lla2ecef(lla0);
    Vec3<Float> xyz = lla2ecef(lla);
    enu = C_e_n * (xyz - xyz0);
}
template <typename Float = double>
Vec3<Float> lla2enu(const Vec3<Float> &lla, const Vec3<Float> &lla0) {
    Vec3<Float> enu;
    lla2ned(enu, lla, lla0);
    return enu;
}

//! --- LLA2AER ---
/// @brief      Latitude-Longitude-Height to Azimuth-Elevation-Range position coordinates
/// @param aer  3x1 AER position [m]
/// @param llaR 3x1 Reference Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @param llaT 3x1 Target Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @returns    AER from reference to target
template <typename Float = double>
void lla2aer(Vec3<Float> &aer, const Vec3<Float> &llaR, const Vec3<Float> &llaT) {
    Vec3<Float> enu = lla2enu(llaT, llaR);
    aer(2) = enu.norm();
    aer(1) = std::asin(enu(2) / aer(2));
    aer(0) = std::atan2(enu(0), enu(1));
}
template <typename Float = double>
Vec3<Float> lla2aer(const Vec3<Float> &llaR, const Vec3<Float> &llaT) {
    Vec3<Float> aer;
    lla2aer(aer, llaR, llaT);
    return aer;
}

//! --- ECI2ECEF ---
/// @brief      Earth-Centered-Inertial to Earth-Centered-Earth-Fixed position coordinates
/// @param eci  3x1 ECI position [m]
/// @param xyz  3x1 ECEF position [m]
/// @param dt   time elapsed between frames [s]
/// @returns    ECEF position
template <typename Float = double>
void eci2ecef(Vec3<Float> &xyz, const Vec3<Float> &eci, const Float &dt) {
    Mat3x3<Float> C_i_e = eci2ecefDcm(dt);
    xyz = C_i_e * eci;
}
template <typename Float = double>
Vec3<Float> eci2ecef(const Vec3<Float> &eci, const Float &dt) {
    Vec3<Float> xyz;
    eci2ecef(xyz, eci, dt);
    return xyz;
}

//! --- ECI2LLA ---
/// @brief      Earth-Centered-Inertial to Latitude-Longitude-Height position coordinates
/// @param eci  3x1 ECI position [m]
/// @param lla  3x1 LLA position [rad, rad, m]
/// @param dt   time elapsed between frames [s]
/// @returns    lla position
template <typename Float = double>
void eci2lla(Vec3<Float> &lla, const Vec3<Float> &eci, const Float &dt) {
    Vec3<Float> xyz = eci2ecef(eci, dt);
    lla = ecef2lla(xyz);
}
template <typename Float = double>
Vec3<Float> eci2lla(const Vec3<Float> &eci, const Float &dt) {
    Vec3<Float> lla;
    eci2lla(lla, eci, dt);
    return lla;
}

//! --- ECI2NED ---
/// @brief      Earth-Centered-Inertial to North-East-Down position coordinates
/// @param eci  3x1 ECI position [m]
/// @param lla0 3x1 Reference Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @param ned  3x1 NED position [m]
/// @param dt   time elapsed between frames [s]
/// @returns    lla position
template <typename Float = double>
void eci2ned(Vec3<Float> &ned, const Vec3<Float> &eci, const Vec3<Float> &lla0, const Float &dt) {
    Vec3<Float> xyz = eci2ecef(eci, dt);
    ned = ecef2ned(xyz, lla0);
}
template <typename Float = double>
Vec3<Float> eci2ned(const Vec3<Float> &eci, const Vec3<Float> &lla0, const Float &dt) {
    Vec3<Float> ned;
    eci2ned(ned, eci, lla0, dt);
    return ned;
}

//! --- ECI2ENU ---
/// @brief      Earth-Centered-Inertial to East-North-Up position coordinates
/// @param eci  3x1 ECI position [m]
/// @param lla0 3x1 Reference Geodetic Latitude, Longitude, Height [rad, rad, m]
/// @param enu  3x1 ENU position [m]
/// @param dt   time elapsed between frames [s]
/// @returns    lla position
template <typename Float = double>
void eci2enu(Vec3<Float> &enu, const Vec3<Float> &eci, const Vec3<Float> &lla0, const Float &dt) {
    Vec3<Float> xyz = eci2ecef(eci, dt);
    enu = ecef2enu(xyz, lla0);
}
template <typename Float = double>
Vec3<Float> eci2enu(const Vec3<Float> &eci, const Vec3<Float> &lla0, const Float &dt) {
    Vec3<Float> enu;
    eci2enu(enu, eci, lla0, dt);
    return enu;
}

//! --- ECI2AER ---
/// @brief      Earth-Centered-Inertial to Azimuth-Elevation-Range position coordinates
/// @param aer  3x1 AER position [m]
/// @param eciR 3x1 Reference ECI position [m]
/// @param eciT 3x1 Target ECI position [m]
/// @param dt   time elapsed between frames [s]
/// @returns    AER from reference to target
template <typename Float = double>
void eci2aer(Vec3<Float> &aer, const Vec3<Float> &eciR, const Vec3<Float> &eciT, const Float &dt) {
    aer = ecef2aer(eci2ecef(eciT, dt), eci2ecef(eciR, dt));
}
template <typename Float = double>
Vec3<Float> eci2aer(const Vec3<Float> &eciR, const Vec3<Float> &eciT, const Float &dt) {
    Vec3<Float> aer;
    eci2aer(aer, eciR, eciT, dt);
    return aer;
}

//! --- ECEF2ECI ---
/// @brief      Earth-Centered-Earth-Fixed to Earth-Centered-Inertial position coordinates
/// @param xyz  3x1 ECEF position [m]
/// @param eci  3x1 ECI position [m]
/// @param dt   time elapsed between frames [s]
/// @returns    ECEF position
template <typename Float = double>
void ecef2eci(Vec3<Float> &eci, const Vec3<Float> &xyz, const Float &dt) {
    Mat3x3<Float> C_e_i = ecef2eciDcm(dt);
    eci = C_e_i * xyz;
}
template <typename Float = double>
Vec3<Float> ecef2eci(const Vec3<Float> &xyz, const Float &dt) {
    Vec3<Float> eci;
    ecef2eci(eci, xyz, dt);
    return eci;
}

//! --- ECEF2LLA ---
/// @brief      Earth-Centered-Earth-Fixed to Latitude-Longitude-Height position coordinates
/// @param xyz  3x1 ECEF position [m]
/// @param lla  3x1 LLA position [rad, rad, m]
/// @returns    lla position
template <typename Float = double>
void ecef2lla(Vec3<Float> &lla, const Vec3<Float> &xyz) {
    Float x = xyz(0);
    Float y = xyz(1);
    Float z = xyz(2);

    Float sign_z = std::copysign(1.0, z);
    Float sqrt_1_e2 = std::sqrt(1.0 - WGS84_E2<Float>);

    Float beta = std::sqrt(x * x + y * y);  // (Groves C.18)
    Float a = sqrt_1_e2 * std::abs(z);
    Float b = WGS84_E2<Float> * WGS84_R0<Float>;
    Float E = (a - b) / beta;             // (Groves C.29)
    Float F = (a + b) / beta;             // (Groves C.30)
    Float P = 4.0 / 3.0 * (E * F + 1.0);  // (Groves C.31)
    Float Q = 2.0 * (E * E - F * F);      // (Groves C.32)
    Float D = P * P * P + Q * Q;          // (Groves C.33)
    Float sqrt_D = std::sqrt(D);
    Float V = std::pow(sqrt_D - Q, 1.0 / 3.0) - std::pow(sqrt_D + Q, 1.0 / 3.0);  // (Groves C.34)
    Float G = 0.5 * (std::sqrt(E * E + V) + E);                                   // (Groves C.35)
    Float T = std::sqrt(G * G + ((F - V * G) / (2.0 * G - E))) - G;               // (Groves C.36)
    lla(0) = sign_z * std::atan((1.0 - T * T) / (2.0 * T * sqrt_1_e2));           // (Groves C.37)
    lla(1) = std::atan2(y, x);
    lla(2) = (beta - WGS84_R0<Float> * T) * std::cos(lla(0)) +
             (z - sign_z * WGS84_R0<Float> * sqrt_1_e2) * std::sin(lla(0));  // (Groves C.38)
}
template <typename Float = double>
Vec3<Float> ecef2lla(const Vec3<Float> &xyz) {
    Vec3<Float> lla;
    ecef2lla(lla, xyz);
    return lla;
}

//! --- ECEF2NED ---
/// @brief      Earth-Centered-Earth-Fixed to North-East-Down position coordinates
/// @param xyz  3x1 ECEF position [m]
/// @param ned  3x1 NED position [m]
/// @param lla0 3x1 Reference LLA position [rad, rad, m]
/// @returns    NED position
template <typename Float = double>
void ecef2ned(Vec3<Float> &ned, const Vec3<Float> &xyz, const Vec3<Float> &lla0) {
    Mat3x3<Float> C_e_n = ecef2nedDcm(lla0);
    Vec3<Float> xyz0 = lla2ecef(lla0);
    ned = C_e_n * (xyz - xyz0);
}
template <typename Float = double>
Vec3<Float> ecef2ned(const Vec3<Float> &xyz, const Vec3<Float> &lla0) {
    Vec3<Float> ned;
    ecef2ned(ned, xyz, lla0);
    return ned;
}

//! --- ECEF2ENU ---
/// @brief      Earth-Centered-Earth-Fixed to East-North-Up position coordinates
/// @param xyz  3x1 ECEF position [m]
/// @param enu  3x1 ENU position [m]
/// @param lla0 3x1 Reference LLA position [rad, rad, m]
/// @returns    ENU position
template <typename Float = double>
void ecef2enu(Vec3<Float> &enu, const Vec3<Float> &xyz, const Vec3<Float> &lla0) {
    Mat3x3<Float> C_e_n = ecef2enuDcm(lla0);
    Vec3<Float> xyz0 = lla2ecef(lla0);
    enu = C_e_n * (xyz - xyz0);
}
template <typename Float = double>
Vec3<Float> ecef2enu(const Vec3<Float> &xyz, const Vec3<Float> &lla0) {
    Vec3<Float> enu;
    ecef2enu(enu, xyz, lla0);
    return enu;
}

//! --- ECEF2AER ---
/// @brief      Earth-Centered-Earth-Fixed to Azimuth-Elevation-Range position coordinates
/// @param aer  3x1 AER position [m]
/// @param xyzR 3x1 Reference ECEF position [m]
/// @param xyzT 3x1 Target ECEF position [m]
/// @returns    AER position
template <typename Float = double>
void ecef2aer(Vec3<Float> &aer, const Vec3<Float> &xyzR, const Vec3<Float> &xyzT) {
    Vec3<Float> lla0 = ecef2lla(xyzR);
    Mat3x3<Float> C_e_n = ecef2enuDcm(lla0);
    Vec3<Float> enu = C_e_n * (xyzT - xyzR);

    aer(2) = enu.norm();
    aer(1) = std::asin(enu(2) / aer(2));
    aer(0) = std::atan2(enu(0), enu(1));
}
template <typename Float = double>
Vec3<Float> ecef2aer(const Vec3<Float> &xyzR, const Vec3<Float> &xyzT) {
    Vec3<Float> aer;
    ecef2aer(aer, xyzR, xyzT);
    return aer;
}

//! --- NED2ECI ---
/// @brief      North-East-Down to Earth-Centered-Inertial position coordinates
/// @param ned  3x1 NED position [m]
/// @param eci  3x1 ECI position [m]
/// @param lla0 3x1 Reference LLA position [rad, rad, m]
/// @param dt   time elapsed between frames [s]
/// @returns    ECI position
template <typename Float = double>
void ned2eci(Vec3<Float> &eci, const Vec3<Float> &ned, const Vec3<Float> &lla0, const Float &dt) {
    Vec3<Float> xyz = ned2ecef(ned, lla0);
    Mat3x3<Float> C_e_i = ecef2eciDcm(dt);
    eci = C_e_i * xyz;
}
template <typename Float = double>
Vec3<Float> ned2eci(const Vec3<Float> &ned, const Vec3<Float> &lla0, const Float &dt) {
    Vec3<Float> eci;
    ned2eci(eci, ned, lla0, dt);
    return eci;
}

//! --- NED2ECEF ---
/// @brief      North-East-Down to Earth-Centered-Earth-Fixed position coordinates
/// @param ned  3x1 NED position [m]
/// @param xyz  3x1 ECEF position [m]
/// @param lla0 3x1 Reference LLA position [rad, rad, m]
/// @returns    ECEF position
template <typename Float = double>
void ned2ecef(Vec3<Float> &xyz, const Vec3<Float> &ned, const Vec3<Float> &lla0) {
    Mat3x3<Float> C_n_e = ned2ecefDcm(lla0);
    xyz = lla2ecef(lla0) + C_n_e * ned;
}
template <typename Float = double>
Vec3<Float> ned2ecef(const Vec3<Float> &ned, const Vec3<Float> &lla0) {
    Vec3<Float> xyz;
    ned2ecef(xyz, ned, lla0);
    return xyz;
}

//! --- NED2LLA ---
/// @brief      North-East-Down to Latitude-Longitude-Height position coordinates
/// @param ned  3x1 NED position [m]
/// @param lla  3x1 LLA position [m]
/// @param lla0 3x1 Reference LLA position [rad, rad, m]
/// @returns    LLA position
template <typename Float = double>
void ned2lla(Vec3<Float> &lla, const Vec3<Float> &ned, const Vec3<Float> &lla0) {
    Vec3<Float> xyz = ned2ecef(ned, lla0);
    lla = ecef2lla(xyz);
}
template <typename Float = double>
Vec3<Float> ned2lla(const Vec3<Float> &ned, const Vec3<Float> &lla0) {
    Vec3<Float> lla;
    ned2ecef(lla, ned, lla0);
    return lla;
}

//! --- NED2AER ---
/// @brief      North-East-Down to Azimuth-Elevation-Range position coordinates
/// @param aer  3x1 AER position [m]
/// @param nedR 3x1 Reference NED position [m]
/// @param nedT 3x1 Target NED position [m]
/// @returns    AER position
template <typename Float = double>
void ned2aer(Vec3<Float> &aer, const Vec3<Float> &nedR, const Vec3<Float> &nedT) {
    Vec3<Float> d_ned = nedT - nedR;
    aer(2) = d_ned.norm();
    aer(1) = std::asin(-d_ned(2) / aer(2));
    aer(0) = std::atan2(d_ned(1), d_ned(0));
}
template <typename Float = double>
Vec3<Float> ned2aer(const Vec3<Float> &nedR, const Vec3<Float> &nedT) {
    Vec3<Float> aer;
    ned2aer(aer, nedR, nedT);
    return aer;
}

//! --- ENU2ECI ---
/// @brief      East-North-Up to Earth-Centered-Inertial position coordinates
/// @param enu  3x1 ENU position [m]
/// @param eci  3x1 ECI position [m]
/// @param lla0 3x1 Reference LLA position [rad, rad, m]
/// @param dt   time elapsed between frames [s]
/// @returns    ECI position
template <typename Float = double>
void enu2eci(Vec3<Float> &eci, const Vec3<Float> &enu, const Vec3<Float> &lla0, const Float &dt) {
    Vec3<Float> xyz = enu2ecef(enu, lla0);
    Mat3x3<Float> C_e_i = ecef2eciDcm(dt);
    eci = C_e_i * xyz;
}
template <typename Float = double>
Vec3<Float> enu2eci(const Vec3<Float> &enu, const Vec3<Float> &lla0, const Float &dt) {
    Vec3<Float> eci;
    enu2eci(eci, enu, lla0, dt);
    return eci;
}

//! --- ENU2ECEF ---
/// @brief      East-North-Up to Earth-Centered-Earth-Fixed position coordinates
/// @param enu  3x1 ENU position [m]
/// @param xyz  3x1 ECEF position [m]
/// @param lla0 3x1 Reference LLA position [rad, rad, m]
/// @returns    ECEF position
template <typename Float = double>
void enu2ecef(Vec3<Float> &xyz, const Vec3<Float> &enu, const Vec3<Float> &lla0) {
    Mat3x3<Float> C_n_e = enu2ecefDcm(lla0);
    xyz = lla2ecef(lla0) + C_n_e * enu;
}
template <typename Float = double>
Vec3<Float> enu2ecef(const Vec3<Float> &enu, const Vec3<Float> &lla0) {
    Vec3<Float> xyz;
    enu2ecef(xyz, enu, lla0);
    return xyz;
}

//! --- ENU2LLA ---
/// @brief      East-North-Up to Latitude-Longitude-Height position coordinates
/// @param enu  3x1 ENU position [m]
/// @param lla  3x1 LLA position [m]
/// @param lla0 3x1 Reference LLA position [rad, rad, m]
/// @returns    LLA position
template <typename Float = double>
void enu2lla(Vec3<Float> &lla, const Vec3<Float> &enu, const Vec3<Float> &lla0) {
    Vec3<Float> xyz = enu2ecef(enu, lla0);
    lla = ecef2lla(xyz);
}
template <typename Float = double>
Vec3<Float> enu2lla(const Vec3<Float> &enu, const Vec3<Float> &lla0) {
    Vec3<Float> lla;
    enu2ecef(lla, enu, lla0);
    return lla;
}

//! --- ENU2AER ---
/// @brief      East-North-Up to Azimuth-Elevation-Range position coordinates
/// @param aer  3x1 AER position [m]
/// @param enuR 3x1 Reference NED position [m]
/// @param enuT 3x1 Target NED position [m]
/// @returns    AER position
template <typename Float = double>
void enu2aer(Vec3<Float> &aer, const Vec3<Float> &enuR, const Vec3<Float> &enuT) {
    Vec3<Float> d_enu = enuT - enuR;
    aer(2) = d_enu.norm();
    aer(1) = std::asin(d_enu(2) / aer(2));
    aer(0) = std::atan2(d_enu(0), d_enu(1));
}
template <typename Float = double>
Vec3<Float> enu2aer(const Vec3<Float> &enuR, const Vec3<Float> &enuT) {
    Vec3<Float> aer;
    ned2aer(aer, enuR, enuT);
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
template <typename Float = double>
void eci2ecefv(
        Vec3<Float> &xyz, const Vec3<Float> &r_ib_i, const Vec3<Float> &v_ib_i, const Float &dt) {
    Mat3x3<Float> C_i_e = eci2ecefDcm(dt);
    xyz = C_i_e * (v_ib_i - WGS84_OMEGA_SKEW<Float> * r_ib_i);
}
template <typename Float = double>
Vec3<Float> eci2ecefv(const Vec3<Float> &r_ib_i, const Vec3<Float> &v_ib_i, const Float &dt) {
    Vec3<Float> xyz;
    eci2ecefv(xyz, r_ib_i, v_ib_i, dt);
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
template <typename Float = double>
void eci2nedv(
        Vec3<Float> &ned,
        const Vec3<Float> &r_ib_i,
        const Vec3<Float> &v_ib_i,
        const Vec3<Float> &lla0,
        const Float &dt) {
    Mat3x3<Float> C_i_n = eci2nedDcm(lla0, dt);
    ned = C_i_n * (v_ib_i - WGS84_OMEGA_SKEW<Float> * r_ib_i);
}
template <typename Float = double>
Vec3<Float> eci2nedv(
        const Vec3<Float> &r_ib_i,
        const Vec3<Float> &v_ib_i,
        const Vec3<Float> &lla0,
        const Float &dt) {
    Vec3<Float> ned;
    eci2nedv(ned, r_ib_i, v_ib_i, lla0, dt);
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
template <typename Float = double>
void eci2enuv(
        Vec3<Float> &enu,
        const Vec3<Float> &r_ib_i,
        const Vec3<Float> &v_ib_i,
        const Vec3<Float> &lla0,
        const Float &dt) {
    Mat3x3<Float> C_i_n = eci2enuDcm(lla0, dt);
    enu = C_i_n * (v_ib_i - WGS84_OMEGA_SKEW<Float> * r_ib_i);
}
template <typename Float = double>
Vec3<Float> eci2enuv(
        const Vec3<Float> &r_ib_i,
        const Vec3<Float> &v_ib_i,
        const Vec3<Float> &lla0,
        const Float &dt) {
    Vec3<Float> enu;
    eci2enuv(enu, r_ib_i, v_ib_i, lla0, dt);
    return enu;
}

//! --- ECEF2ECIV ---
/// @brief      Converts Earth-Centered-Earth-Fixed to Earth-Centered-Inertial velocity
/// @param r_eb_e   3x1 ECEF position [m]
/// @param v_eb_e   3x1 ECEF velocity [m/s]
/// @param eci      3x1 ECI velocity [m/s]
/// @param dt       time elapsed between frames [s]
/// @returns    ECI velocity
template <typename Float = double>
void ecef2eciv(
        Vec3<Float> &eci, const Vec3<Float> &r_eb_e, const Vec3<Float> &v_eb_e, const Float &dt) {
    Vec3<Float> C_e_i = ecef2eciDcm(dt);
    eci = C_e_i * (v_eb_e - WGS84_OMEGA_SKEW<Float> * r_eb_e);
}
template <typename Float = double>
Vec3<Float> ecef2eciv(const Vec3<Float> &r_eb_e, const Vec3<Float> &v_eb_e, const Float &dt) {
    Vec3<Float> eci;
    ecef2eciv(eci, r_eb_e, v_eb_e, dt);
    return eci;
}

//! --- ECEF2NEDV ---
/// @brief      Converts Earth-Centered-Earth-Fixed to North-East-Down velocity
/// @param v_eb_e   3x1 ECEF velocity [m/s]
/// @param lla0     3x1 Reference LLA position [rad, rad, m]
/// @param ned      3x1 NED velocity [m/s]
/// @returns    NED velocity
template <typename Float = double>
void ecef2nedv(Vec3<Float> &ned, const Vec3<Float> &v_eb_e, const Vec3<Float> &lla0) {
    Mat3x3<Float> C_e_n = ecef2nedDcm(lla0);
    ned = C_e_n * v_eb_e;
}
template <typename Float = double>
Vec3<Float> ecef2nedv(const Vec3<Float> &v_eb_e, const Vec3<Float> &lla0) {
    Vec3<Float> ned;
    ecef2nedv(ned, v_eb_e, lla0);
    return ned;
}

//! --- ECEF2ENUV ---
/// @brief      Converts Earth-Centered-Earth-Fixed to East-North0Up velocity
/// @param v_eb_e   3x1 ECEF velocity [m/s]
/// @param lla0     3x1 Reference LLA position [rad, rad, m]
/// @param enu      3x1 NED velocity [m/s]
/// @returns    ENU velocity
template <typename Float = double>
void ecef2enuv(Vec3<Float> &enu, const Vec3<Float> &v_eb_e, const Vec3<Float> &lla0) {
    Mat3x3<Float> C_e_n = ecef2enuDcm(lla0);
    enu = C_e_n * v_eb_e;
}
template <typename Float = double>
Vec3<Float> ecef2enuv(const Vec3<Float> &v_eb_e, const Vec3<Float> &lla0) {
    Vec3<Float> enu;
    ecef2enuv(enu, v_eb_e, lla0);
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
template <typename Float = double>
void ned2eciv(
        Vec3<Float> &eci,
        const Vec3<Float> &r_nb_e,
        const Vec3<Float> &v_nb_e,
        const Vec3<Float> &lla0,
        const Float &dt) {
    Mat3x3<Float> C_n_i = ned2eciDcm(lla0, dt);
    Mat3x3<Float> C_e_i = ecef2eciDcm(dt);
    Vec3<Float> xyz = ned2ecef(r_nb_e, lla0);
    eci = C_n_i * v_nb_e + C_e_i * WGS84_OMEGA_SKEW<Float> * xyz;
}
template <typename Float = double>
Vec3<Float> ned2eciv(
        const Vec3<Float> &r_nb_e,
        const Vec3<Float> &v_nb_e,
        const Vec3<Float> &lla0,
        const Float &dt) {
    Vec3<Float> eci;
    ned2eciv(eci, r_nb_e, v_nb_e, lla0, dt);
    return eci;
}

//! --- NED2ECEFV ---
/// @brief      Converts North-East-Down to Earth-Centered-Earth-Fixed velocity
/// @param v_nb_e   3x1 NED velocity [m/s]
/// @param lla0     3x1 Reference LLA position [rad, rad, m]
/// @param xyz      3x1 ECEF velocity [m/s]
/// @returns    ENU velocity
template <typename Float = double>
void ned2ecefv(Vec3<Float> &xyz, const Vec3<Float> &v_nb_e, const Vec3<Float> &lla0) {
    Mat3x3<Float> C_n_e = ned2ecefDcm(lla0);
    xyz = C_n_e * v_nb_e;
}
template <typename Float = double>
Vec3<Float> ned2ecefv(const Vec3<Float> &v_nb_e, const Vec3<Float> &lla0) {
    Vec3<Float> xyz;
    ned2ecefv(xyz, v_nb_e, lla0);
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
template <typename Float = double>
void enu2eciv(
        Vec3<Float> &eci,
        const Vec3<Float> &r_nb_e,
        const Vec3<Float> &v_nb_e,
        const Vec3<Float> &lla0,
        const Float &dt) {
    Mat3x3<Float> C_n_i = enu2eciDcm(lla0, dt);
    Mat3x3<Float> C_e_i = ecef2eciDcm(dt);
    Vec3<Float> xyz = ned2ecef(r_nb_e, lla0);
    eci = C_n_i * v_nb_e + C_e_i * WGS84_OMEGA_SKEW<Float> * xyz;
}
template <typename Float = double>
Vec3<Float> enu2eciv(
        const Vec3<Float> &r_nb_e,
        const Vec3<Float> &v_nb_e,
        const Vec3<Float> &lla0,
        const Float &dt) {
    Vec3<Float> eci;
    enu2eciv(eci, r_nb_e, v_nb_e, lla0, dt);
    return eci;
}

//! --- ENU2ECEFV ---
/// @brief      Converts East-North-Up to Earth-Centered-Earth-Fixed velocity
/// @param v_nb_e   3x1 ENU velocity [m/s]
/// @param lla0     3x1 Reference LLA position [rad, rad, m]
/// @param xyz      3x1 ECEF velocity [m/s]
/// @returns    ENU velocity
template <typename Float = double>
void enu2ecefv(Vec3<Float> &xyz, const Vec3<Float> &v_nb_e, const Vec3<Float> &lla0) {
    Mat3x3<Float> C_n_e = enu2ecefDcm(lla0);
    xyz = C_n_e * v_nb_e;
}
template <typename Float = double>
Vec3<Float> enu2ecefv(const Vec3<Float> &v_nb_e, const Vec3<Float> &lla0) {
    Vec3<Float> xyz;
    enu2ecefv(xyz, v_nb_e, lla0);
    return xyz;
}

//* ===== Angular Velocity Transformations ===================================================== *//

//! --- ECI2ECEFW ---
/// @brief      Converts Earth-Centered-Inertial to Earth-Centered-Earth-Fixed angular velocity
/// @param w_ib_i   3x1 ECI angular velocity [rad/s]
/// @param dt       time elapsed between frames [s]
/// @param xyz      3x1 ECEF angular velocity [rad/s]
/// @returns    ECEF angular velocity
template <typename Float = double>
void eci2ecefw(Vec3<Float> &xyz, const Vec3<Float> &w_ib_i, const Float &dt) {
    Mat3x3<Float> C_i_e = eci2ecefDcm(dt);
    xyz = C_i_e * (w_ib_i + WGS84_OMEGA_VEC<Float>);
}
template <typename Float = double>
Vec3<Float> eci2ecefw(const Vec3<Float> &w_ib_i, const Float &dt) {
    Vec3<Float> xyz;
    eci2ecefw(xyz, w_ib_i, dt);
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
template <typename Float = double>
void eci2nedw(
        Vec3<Float> &ned,
        const Vec3<Float> &w_ib_i,
        const Vec3<Float> &r_ib_i,
        const Vec3<Float> &v_ib_i,
        const Vec3<Float> &lla0,
        const Float &dt) {
    Vec3<Float> C_i_n = eci2nedDcm(lla0, dt);

    Vec3<Float> v_nb_e = eci2nedv(r_ib_i, v_ib_i, lla0, dt);
    Float vn = v_nb_e(0);
    Float ve = v_nb_e(1);
    Float phi = lla0(0);
    Float h = lla0(2);
    Float sin_phi = std::sin(phi);

    Float trans = 1.0 - WGS84_E2<Float> * sin_phi * sin_phi;
    Float Re = WGS84_R0<Float> / std::sqrt(trans);
    Float Rn = WGS84_R0<Float> * (1.0 - WGS84_E2<Float>) / std::pow(trans, 1.5);
    Vec3<Float> w_en_n{ve / (Re + h), -vn / (Rn + h), -ve * std::tan(phi) / (Re + h)};

    ned = -w_en_n + C_i_n * (w_ib_i - WGS84_OMEGA_VEC<Float>);
}
template <typename Float = double>
Vec3<Float> eci2nedw(
        const Vec3<Float> &w_ib_i,
        const Vec3<Float> &r_ib_i,
        const Vec3<Float> &v_ib_i,
        const Vec3<Float> &lla0,
        const Float &dt) {
    Vec3<Float> ned;
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
template <typename Float = double>
void eci2enuw(
        Vec3<Float> &enu,
        const Vec3<Float> &w_ib_i,
        const Vec3<Float> &r_ib_i,
        const Vec3<Float> &v_ib_i,
        const Vec3<Float> &lla0,
        const Float &dt) {
    Vec3<Float> C_i_n = eci2enuDcm(lla0, dt);

    Vec3<Float> v_nb_e = eci2enuv(r_ib_i, v_ib_i, lla0, dt);
    Float vn = v_nb_e(0);
    Float ve = v_nb_e(1);
    Float phi = lla0(0);
    Float h = lla0(2);
    Float sin_phi = std::sin(phi);

    Float trans = 1.0 - WGS84_E2<Float> * sin_phi * sin_phi;
    Float Re = WGS84_R0<Float> / std::sqrt(trans);
    Float Rn = WGS84_R0<Float> * (1.0 - WGS84_E2<Float>) / std::pow(trans, 1.5);
    Vec3<Float> w_en_n{-vn / (Rn + h), ve / (Re + h), ve * std::tan(phi) / (Re + h)};

    enu = -w_en_n + C_i_n * (w_ib_i - WGS84_OMEGA_VEC<Float>);
}
template <typename Float = double>
Vec3<Float> eci2enuw(
        const Vec3<Float> &w_ib_i,
        const Vec3<Float> &r_ib_i,
        const Vec3<Float> &v_ib_i,
        const Vec3<Float> &lla0,
        const Float &dt) {
    Vec3<Float> enu;
    eci2enuw(enu, w_ib_i, r_ib_i, v_ib_i, lla0, dt);
    return enu;
}

//! --- ECEF2ECIW ---
/// @brief      Converts Earth-Centered-Earth-Fixed to Earth-Centered-Inertial angular velocity
/// @param w_eb_e   3x1 ECEF angular velocity [rad/s]
/// @param dt       time elapsed between frames [s]
/// @param eci      3x1 ECI angular velocity [rad/s]
/// @returns    ECI angular velocity
template <typename Float = double>
void ecef2eciw(Vec3<Float> &eci, const Vec3<Float> &w_eb_e, const Float &dt) {
    Mat3x3<Float> C_e_i = ecef2eciDcm(dt);
    eci = C_e_i * (w_eb_e + WGS84_OMEGA_VEC<Float>);
}
template <typename Float = double>
Vec3<Float> ecef2eciw(const Vec3<Float> &w_eb_e, const Float &dt) {
    Vec3<Float> eci;
    ecef2eciw(eci, w_eb_e, dt);
    return eci;
}

//! --- ECEF2NEDW ---
/// @brief      Converts Earth-Centered-Earth-Fixed to North-East-Down angular velocity
/// @param w_eb_e   3x1 ECEF angular velocity [rad/s]
/// @param lla0     Reference LLA [rad, rad, m]
/// @param ned      3x1 NED angular velocity [rad/s]
/// @returns    NED angular velocity
template <typename Float = double>
void ecef2nedw(Vec3<Float> &ned, const Vec3<Float> &w_eb_e, const Vec3<Float> &lla0) {
    Mat3x3<Float> C_e_n = ecef2nedDcm(lla0);
    ned = C_e_n * w_eb_e;
}
template <typename Float = double>
Vec3<Float> ecef2nedw(const Vec3<Float> &w_eb_e, const Vec3<Float> &lla0) {
    Vec3<Float> ned;
    ecef2nedw(ned, w_eb_e, lla0);
    return ned;
}

//! --- ECEF2ENUW ---
/// @brief      Converts Earth-Centered-Earth-Fixed to East-North-Up angular velocity
/// @param w_eb_e   3x1 ECEF angular velocity [rad/s]
/// @param lla0     Reference LLA [rad, rad, m]
/// @param ned      3x1 ECU angular velocity [rad/s]
/// @returns    ENU angular velocity
template <typename Float = double>
void ecef2enuw(Vec3<Float> &enu, const Vec3<Float> &w_eb_e, const Vec3<Float> &lla0) {
    Mat3x3<Float> C_e_n = ecef2enuDcm(lla0);
    enu = C_e_n * w_eb_e;
}
template <typename Float = double>
Vec3<Float> ecef2enuw(const Vec3<Float> &w_eb_e, const Vec3<Float> &lla0) {
    Vec3<Float> enu;
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
template <typename Float = double>
void ned2eciw(
        Vec3<Float> &eci,
        const Vec3<Float> &w_nb_e,
        const Vec3<Float> &v_nb_e,
        const Vec3<Float> &lla0,
        const Float &dt) {
    Vec3<Float> C_n_i = ned2eciDcm(lla0, dt);

    Float vn = v_nb_e(0);
    Float ve = v_nb_e(1);
    Float phi = lla0(0);
    Float h = lla0(2);
    Float sin_phi = std::sin(phi);

    Float trans = 1.0 - WGS84_E2<Float> * sin_phi * sin_phi;
    Float Re = WGS84_R0<Float> / std::sqrt(trans);
    Float Rn = WGS84_R0<Float> * (1.0 - WGS84_E2<Float>) / std::pow(trans, 1.5);
    Vec3<Float> w_en_n{ve / (Re + h), -vn / (Rn + h), -ve * std::tan(phi) / (Re + h)};

    eci = C_n_i * (w_nb_e + w_en_n) + WGS84_OMEGA_VEC<Float>;
}
template <typename Float = double>
Vec3<Float> ned2eciw(
        const Vec3<Float> &w_nb_e,
        const Vec3<Float> &v_nb_e,
        const Vec3<Float> &lla0,
        const Float &dt) {
    Vec3<Float> eci;
    ned2eciw(eci, w_nb_e, v_nb_e, lla0, dt);
    return eci;
}

//! --- NED2ECEFW ---
/// @brief      Converts North-East-Down to Earth-Centered-Earth-Fixed angular velocity
/// @param w_nb_e   3x1 NED angular velocity [rad/s]
/// @param lla0     Reference LLA [rad, rad, m]
/// @param xyz      3x1 ECEF angular velocity [rad/s]
/// @returns    NED angular velocity
template <typename Float = double>
void ned2ecefw(Vec3<Float> &xyz, const Vec3<Float> &w_nb_e, const Vec3<Float> &lla0) {
    Mat3x3<Float> C_n_e = ned2ecefDcm(lla0);
    xyz = C_n_e * w_nb_e;
}
template <typename Float = double>
Vec3<Float> ned2ecefw(const Vec3<Float> &w_nb_e, const Vec3<Float> &lla0) {
    Vec3<Float> xyz;
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
template <typename Float = double>
void enu2eciw(
        Vec3<Float> &eci,
        const Vec3<Float> &w_nb_e,
        const Vec3<Float> &v_nb_e,
        const Vec3<Float> &lla0,
        const Float &dt) {
    Vec3<Float> C_n_i = ned2eciDcm(lla0, dt);

    Float vn = v_nb_e(0);
    Float ve = v_nb_e(1);
    Float phi = lla0(0);
    Float h = lla0(2);
    Float sin_phi = std::sin(phi);

    Float trans = 1.0 - WGS84_E2<Float> * sin_phi * sin_phi;
    Float Re = WGS84_R0<Float> / std::sqrt(trans);
    Float Rn = WGS84_R0<Float> * (1.0 - WGS84_E2<Float>) / std::pow(trans, 1.5);
    Vec3<Float> w_en_n{-vn / (Rn + h), ve / (Re + h), ve * std::tan(phi) / (Re + h)};

    eci = C_n_i * (w_nb_e + w_en_n) + WGS84_OMEGA_VEC<Float>;
}
template <typename Float = double>
Vec3<Float> enu2eciw(
        const Vec3<Float> &w_nb_e,
        const Vec3<Float> &v_nb_e,
        const Vec3<Float> &lla0,
        const Float &dt) {
    Vec3<Float> eci;
    enu2eciw(eci, w_nb_e, v_nb_e, lla0, dt);
    return eci;
}

//! --- ENU2ECEFW ---
/// @brief      Converts North-East-Down to Earth-Centered-Earth-Fixed angular velocity
/// @param w_nb_e   3x1 NED angular velocity [rad/s]
/// @param lla0     Reference LLA [rad, rad, m]
/// @param xyz      3x1 ECEF angular velocity [rad/s]
/// @returns    NED angular velocity
template <typename Float = double>
void enu2ecefw(Vec3<Float> &xyz, const Vec3<Float> &w_nb_e, const Vec3<Float> &lla0) {
    Mat3x3<Float> C_n_e = enu2ecefDcm(lla0);
    xyz = C_n_e * w_nb_e;
}
template <typename Float = double>
Vec3<Float> enu2ecefw(const Vec3<Float> &w_nb_e, const Vec3<Float> &lla0) {
    Vec3<Float> xyz;
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
template <typename Float = double>
void eci2ecefa(
        Vec3<Float> &xyz,
        const Vec3<Float> &a_ib_i,
        const Vec3<Float> &r_ib_i,
        const Vec3<Float> &v_ib_i,
        const Float &dt) {
    Mat3x3<Float> C_i_e = eci2ecefDcm(dt);
    xyz = C_i_e * (a_ib_i - 2.0 * WGS84_OMEGA_SKEW<Float> * v_ib_i +
                   WGS84_OMEGA_SKEW<Float> * WGS84_OMEGA_SKEW<Float> * r_ib_i);
}
template <typename Float = double>
Vec3<Float> eci2ecefa(
        const Vec3<Float> &a_ib_i,
        const Vec3<Float> &r_ib_i,
        const Vec3<Float> &v_ib_i,
        const Float &dt) {
    Vec3<Float> xyz;
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
template <typename Float = double>
void eci2neda(
        Vec3<Float> &ned,
        const Vec3<Float> &a_ib_i,
        const Vec3<Float> &r_ib_i,
        const Vec3<Float> &v_ib_i,
        const Vec3<Float> &lla0,
        const Float &dt) {
    Mat3x3<Float> C_i_n = eci2nedDcm(lla0, dt);
    ned = C_i_n * (a_ib_i + 2.0 * WGS84_OMEGA_SKEW<Float> * v_ib_i +
                   WGS84_OMEGA_SKEW<Float> * WGS84_OMEGA_SKEW<Float> * r_ib_i);
}
template <typename Float = double>
Vec3<Float> eci2neda(
        const Vec3<Float> &a_ib_i,
        const Vec3<Float> &r_ib_i,
        const Vec3<Float> &v_ib_i,
        const Vec3<Float> &lla0,
        const Float &dt) {
    Vec3<Float> ned;
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
template <typename Float = double>
void eci2enua(
        Vec3<Float> &enu,
        const Vec3<Float> &a_ib_i,
        const Vec3<Float> &r_ib_i,
        const Vec3<Float> &v_ib_i,
        const Vec3<Float> &lla0,
        const Float &dt) {
    Mat3x3<Float> C_i_n = eci2enuDcm(lla0, dt);
    enu = C_i_n * (a_ib_i + 2.0 * WGS84_OMEGA_SKEW<Float> * v_ib_i +
                   WGS84_OMEGA_SKEW<Float> * WGS84_OMEGA_SKEW<Float> * r_ib_i);
}
template <typename Float = double>
Vec3<Float> eci2enua(
        const Vec3<Float> &a_ib_i,
        const Vec3<Float> &r_ib_i,
        const Vec3<Float> &v_ib_i,
        const Vec3<Float> &lla0,
        const Float &dt) {
    Vec3<Float> enu;
    eci2enua(enu, a_ib_i, r_ib_i, v_ib_i, lla0, dt);
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
template <typename Float = double>
void ecef2ecia(
        Vec3<Float> &eci,
        const Vec3<Float> &a_eb_e,
        const Vec3<Float> &r_eb_e,
        const Vec3<Float> &v_eb_e,
        const Float &dt) {
    Mat3x3<Float> C_i_e = eci2ecefDcm(dt);
    eci = C_i_e * (a_eb_e - 2.0 * WGS84_OMEGA_SKEW<Float> * v_eb_e +
                   WGS84_OMEGA_SKEW<Float> * WGS84_OMEGA_SKEW<Float> * r_eb_e);
}
template <typename Float = double>
Vec3<Float> ecef2ecia(
        const Vec3<Float> &a_eb_e,
        const Vec3<Float> &r_eb_e,
        const Vec3<Float> &v_eb_e,
        const Float &dt) {
    Vec3<Float> eci;
    ecef2ecia(eci, a_eb_e, r_eb_e, v_eb_e, dt);
    return eci;
}

//! --- ECEF2NEDA ---
/// @brief      Converts Earth-Centered-Earth-Fixed to North-East-Down acceleration
/// @param a_eb_e   3x1 ECEF acceleration [rad/s]
/// @param lla0     Reference LLA [rad, rad, m]
/// @param ned      3x1 NED acceleration [rad/s]
/// @returns    ENU acceleration
template <typename Float = double>
void ecef2neda(Vec3<Float> &ned, const Vec3<Float> &a_eb_e, const Vec3<Float> &lla0) {
    Mat3x3<Float> C_e_n = ecef2nedDcm(lla0);
    ned = C_e_n * a_eb_e;
}
template <typename Float = double>
Vec3<Float> ecef2neda(const Vec3<Float> &a_eb_e, const Vec3<Float> &lla0) {
    Vec3<Float> ned;
    ecef2neda(ned, a_eb_e, lla0);
    return ned;
}

//! --- ECEF2ENUA ---
/// @brief      Converts Earth-Centered-Earth-Fixed to East-North-Up acceleration
/// @param a_eb_e   3x1 ECEF acceleration [rad/s]
/// @param lla0     Reference LLA [rad, rad, m]
/// @param enu      3x1 ENU acceleration [rad/s]
/// @returns    ENU acceleration
template <typename Float = double>
void ecef2enua(Vec3<Float> &enu, const Vec3<Float> &a_eb_e, const Vec3<Float> &lla0) {
    Mat3x3<Float> C_e_n = ecef2enuDcm(lla0);
    enu = C_e_n * a_eb_e;
}
template <typename Float = double>
Vec3<Float> ecef2enua(const Vec3<Float> &a_eb_e, const Vec3<Float> &lla0) {
    Vec3<Float> enu;
    ecef2enua(enu, a_eb_e, lla0);
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
template <typename Float = double>
void ned2ecia(
        Vec3<Float> &eci,
        const Vec3<Float> &a_nb_e,
        const Vec3<Float> &r_nb_e,
        const Vec3<Float> &v_nb_e,
        const Vec3<Float> &lla0,
        const Float &dt) {
    Vec3<Float> r_eb_e = ned2ecef(r_nb_e, lla0);
    Mat3x3<Float> C_n_e = ned2ecefDcm(lla0);
    Float omega_sin_phi = WGS84_OMEGA<Float> * std::sin(lla0(0));
    Float omega_cos_phi = WGS84_OMEGA<Float> * std::cos(lla0(0));
    Mat3x3<Float> omega_ie_n{
            {0.0, omega_cos_phi, 0.0},
            {-omega_cos_phi, 0.0, -omega_sin_phi},
            {0.0, omega_sin_phi, 0.0}};
    Mat3x3<Float> C_n_i = ned2eciDcm(lla0, dt);
    Mat3x3<Float> C_e_i = ecef2eciDcm(dt);

    auto a_eb_e = C_n_e * a_nb_e;
    auto v_eb_e = C_n_e * v_nb_e;
    eci = C_n_i * (a_eb_e + 2.0 * omega_ie_n * v_eb_e +
                   C_e_i * (WGS84_OMEGA_SKEW<Float> * WGS84_OMEGA_SKEW<Float>, *r_eb_e));
}
template <typename Float = double>
Vec3<Float> ned2ecia(
        const Vec3<Float> &a_nb_e,
        const Vec3<Float> &r_nb_e,
        const Vec3<Float> &v_nb_e,
        const Vec3<Float> &lla0,
        const Float &dt) {
    Vec3<Float> eci;
    ned2ecia(eci, a_nb_e, r_nb_e, v_nb_e, lla0, dt);
    return eci;
}

//! --- NED2ECEFA ---
/// @brief      Converts North-East-Down to Earth-Centered-Earth-Fixed acceleration
/// @param a_nb_e   3x1 NED acceleration [rad/s]
/// @param lla0     Reference LLA [rad, rad, m]
/// @param xyz      3x1 ECEF acceleration [rad/s]
/// @returns    ECEF acceleration
template <typename Float = double>
void ned2ecefa(Vec3<Float> &xyz, const Vec3<Float> &a_nb_e, const Vec3<Float> &lla0) {
    Mat3x3<Float> C_n_e = ned2ecefDcm(lla0);
    xyz = C_n_e * a_nb_e;
}
template <typename Float = double>
Vec3<Float> ned2ecefa(const Vec3<Float> &a_nb_e, const Vec3<Float> &lla0) {
    Vec3<Float> xyz;
    ned2ecefa(xyz, a_nb_e, lla0);
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
template <typename Float = double>
void enu2ecia(
        Vec3<Float> &eci,
        const Vec3<Float> &a_nb_e,
        const Vec3<Float> &r_nb_e,
        const Vec3<Float> &v_nb_e,
        const Vec3<Float> &lla0,
        const Float &dt) {
    Vec3<Float> r_eb_e = enu2ecef(r_nb_e, lla0);
    Mat3x3<Float> C_n_e = enu2ecefDcm(lla0);
    Float omega_sin_phi = WGS84_OMEGA<Float> * std::sin(lla0(0));
    Float omega_cos_phi = WGS84_OMEGA<Float> * std::cos(lla0(0));
    Mat3x3<Float> omega_ie_n{
            {0.0, 0.0, omega_cos_phi},
            {0.0, 0.0, omega_sin_phi},
            {-omega_cos_phi, -omega_sin_phi, 0.0}};
    Mat3x3<Float> C_n_i = enu2eciDcm(lla0, dt);
    Mat3x3<Float> C_e_i = ecef2eciDcm(dt);

    Vec3<Float> a_eb_e = C_n_e * a_nb_e;
    Vec3<Float> v_eb_e = C_n_e * v_nb_e;
    eci = C_n_i * (a_eb_e + 2.0 * omega_ie_n * v_eb_e +
                   C_e_i * (WGS84_OMEGA_SKEW<Float> * WGS84_OMEGA_SKEW<Float>, r_eb_e));
}
template <typename Float = double>
Vec3<Float> enu2ecia(
        const Vec3<Float> &a_nb_e,
        const Vec3<Float> &r_nb_e,
        const Vec3<Float> &v_nb_e,
        const Vec3<Float> &lla0,
        const Float &dt) {
    Vec3<Float> eci;
    enu2ecia(eci, a_nb_e, r_nb_e, v_nb_e, lla0, dt);
    return eci;
}

//! --- ENU2ECEFA ---
/// @brief      Converts East-North-Up to Earth-Centered-Earth-Fixed acceleration
/// @param a_nb_e   3x1 ENU acceleration [rad/s]
/// @param lla0     Reference LLA [rad, rad, m]
/// @param xyz      3x1 ECEF acceleration [rad/s]
/// @returns    ECEF acceleration
template <typename Float = double>
void enu2ecefa(Vec3<Float> &xyz, const Vec3<Float> &a_nb_e, const Vec3<Float> &lla0) {
    Mat3x3<Float> C_n_e = enu2ecefDcm(lla0);
    xyz = C_n_e * a_nb_e;
}
template <typename Float = double>
Vec3<Float> enu2ecefa(const Vec3<Float> &a_nb_e, const Vec3<Float> &lla0) {
    Vec3<Float> xyz;
    enu2ecefa(xyz, a_nb_e, lla0);
    return xyz;
}

//! --- VECEXP ---
/// @brief

//! --- VECLOG ---

}  // namespace navtools

#endif
