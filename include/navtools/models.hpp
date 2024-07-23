/**
|========================================== models.hpp ============================================|
|                                                                                                  |
|   Property of Daniel Sturdivant. Unauthorized copying of this file via any medium is would be    |
|   super sad and unfortunate for me. Proprietary and confidential.                                |
|                                                                                                  |
|--------------------------------------------------------------------------------------------------|
|                                                                                                  |
|   @file     include/navtools/models.hpp                                                          |
|   @brief    Simple models commonly used in navigation equations.                                 |
|   @ref      Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems          |
|               - (2013) Paul D. Groves                                                            |
|   @author   Daniel Sturdivant <sturdivant20@gmail.com>                                           |
|   @date     July 2024                                                                            |
|                                                                                                  |
|==================================================================================================|
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
template <typename Float = double>
void transverseRadius(Float &Re, const Float &phi) {
    Float sin_phi = std::sin(phi);
    Float t = 1.0 - WGS84_E2<Float> * sin_phi * sin_phi;
    Re = WGS84_R0<Float> / std::sqrt(t);
}
template <typename Float = double>
Float transverseRadius(const Float &phi) {
    Float Re;
    transverseRadius(Re, phi);
    return Re;
}

//! === MERIDIANRADIUS ===
/// @brief      Calculates the meridian radius relative to user latitude
/// @param phi  Latitude [rad]
/// @param Rn   Meridian radius [m]
/// @returns    Earth's meridian radius at Latitude
template <typename Float = double>
void meridianRadius(Float &Rn, const Float &phi) {
    Float sin_phi = std::sin(phi);
    Float t = 1.0 - WGS84_E2<Float> * sin_phi * sin_phi;
    Rn = WGS84_R0<Float> * (1.0 - WGS84_E2<Float>) / std::pow(t, 1.5);
}
template <typename Float = double>
Float meridianRadius(const Float &phi) {
    Float Rn;
    meridianRadius(Rn, phi);
    return Rn;
}

//! === GEOCENTRICRADIUS ===
/// @brief      Calculates the geocentric radius relative to user latitude
/// @param phi      Latitude [rad]
/// @param R_es_e   Geocentric radius [m]
/// @returns    Earth's geocentric radius at Latitude
template <typename Float = double>
void geocentricRadius(Float &R_es_e, const Float &phi) {
    Float sin_phi2 = std::sin(phi);
    sin_phi2 *= sin_phi2;
    Float cos_phi = std::cos(phi);
    Float t = 1.0 - WGS84_E2<Float> * sin_phi2;
    Float Re = WGS84_R0<Float> / std::sqrt(t);
    Float o_e2 = 1.0 - WGS84_E2<Float>;
    R_es_e = Re * std::sqrt(cos_phi * cos_phi + o_e2 * o_e2 * sin_phi2);
}
template <typename Float = double>
Float geocentricRadius(const Float &phi) {
    Float R_es_e;
    geocentricRadius(R_es_e, phi);
    return R_es_e;
}

//! === TRANSANDMERRADII ===
/// @brief      Calculates the {Transverse, Meridian} radii relative to user latitude
/// @param phi      Latitude [rad]
/// @param radii    vector containing {Transverse, Meridian} radii [m]
/// @returns    Earth's {Transverse, Meridian} radii at Latitude
template <typename Float = double>
void transAndMerRadii(Vec2<Float> &radii, const Float &phi) {
    Float sin_phi = std::sin(phi);
    Float t = 1.0 - WGS84_E2<Float> * sin_phi * sin_phi;

    radii(0) = WGS84_R0<Float> / std::sqrt(t);
    radii(1) = WGS84_R0<Float> * (1.0 - WGS84_E2<Float>) / std::pow(t, 1.5);
}
template <typename Float = double>
Vec2<Float> transAndMerRadii(const Float &phi) {
    Vec2<Float> radii;
    transAndMerRadii(radii, phi);
    return radii;
}

//! === RADIIOFCURVATURE ===
/// @brief      Calculates the {Transverse, Meridian, Geocentric} radii relative to user latitude
/// @param phi      Latitude [rad]
/// @param radii    vector containing {Transverse, Meridian, Geocentric} radii [m]
/// @returns    Earth's {Transverse, Meridian, Geocentric} radii at Latitude
template <typename Float = double>
void radiiOfCurvature(Vec3<Float> &radii, const Float &phi) {
    Float sin_phi2 = std::sin(phi);
    sin_phi2 *= sin_phi2;
    Float cos_phi = std::cos(phi);
    Float t = 1.0 - WGS84_E2<Float> * sin_phi2;
    Float o_e2 = 1.0 - WGS84_E2<Float>;

    radii(0) = WGS84_R0<Float> / std::sqrt(t);
    radii(1) = WGS84_R0<Float> * o_e2 / std::pow(t, 1.5);
    radii(2) = radii(0) * std::sqrt(cos_phi * cos_phi + o_e2 * o_e2 * sin_phi2);
}
template <typename Float = double>
Vec3<Float> radiiOfCurvature(const Float &phi) {
    Vec3<Float> radii;
    radiiOfCurvature(radii, phi);
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
template <typename Float = double>
void earthRate(Vec3<Float> &w_ie_n, const Float &phi, const std::string frame = "ned") {
    if (frame == "ned") {
        w_ie_n(0) = WGS84_OMEGA<Float> * std::cos(phi);
        w_ie_n(1) = 0.0;
        w_ie_n(2) = WGS84_OMEGA<Float> * std::sin(phi);
    } else if (frame == "enu") {
        w_ie_n(0) = 0.0;
        w_ie_n(1) = WGS84_OMEGA<Float> * std::cos(phi);
        w_ie_n(2) = -WGS84_OMEGA<Float> * std::sin(phi);
    }
}
template <typename Float = double>
Vec3<Float> earthRate(const Float &phi, const std::string frame = "ned") {
    Vec3<Float> w_ie_n;
    earthRate(w_ie_n, phi, frame);
    return w_ie_n;
}
template <typename Float = double>
void earthRateSkew(Mat3x3<Float> &W_ie_n, const Float &phi, const std::string frame = "ned") {
    W_ie_n = skew(earthRate(phi, frame));
}
template <typename Float = double>
Mat3x3<Float> earthRateSkew(const Float &phi, const std::string frame = "ned") {
    Mat3x3<Float> W_ie_n;
    earthRateSkew(W_ie_n, phi, frame);
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
template <typename Float = double>
void transportRate(
        Vec3<Float> &w_en_n,
        const Vec3<Float> &lla,
        const Vec3<Float> &v_nb_e,
        const std::string frame = "ned") {
    Vec2<Float> radii = transAndMerRadii(lla(0));
    if (frame == "ned") {
        Float ve_Reh = v_nb_e(1) / (radii(0) + lla(2));  // ve / (Re + h)
        w_en_n(0) = ve_Reh;                              // ve / (Re + h)
        w_en_n(1) = -v_nb_e(0) / (radii(1) + lla(2));    // -vn / (Rn + h)
        w_en_n(2) = -ve_Reh * std::tan(lla(0));          // -ve * tan(phi) / (Re + h)
    } else if (frame == "enu") {
        Float ve_Reh = v_nb_e(0) / (radii(0) + lla(2));  // ve / (Re + h)
        w_en_n(0) = -v_nb_e(1) / (radii(1) + lla(2));    // -vn / (Rn + h)
        w_en_n(1) = ve_Reh;                              // ve / (Re + h)
        w_en_n(2) = ve_Reh * std::tan(lla(0));           // ve * tan(phi) / (Re + h)
    }
}
template <typename Float = double>
Vec3<Float> transportRate(
        const Vec3<Float> &lla, const Vec3<Float> &v_nb_e, const std::string frame = "ned") {
    Vec3<Float> w_en_n;
    transportRate(w_en_n, lla, v_nb_e, frame);
    return w_en_n;
}
template <typename Float = double>
void transportRateSkew(
        Mat3x3<Float> &W_en_n,
        const Vec3<Float> &lla,
        const Vec3<Float> &v_nb_e,
        const std::string frame = "ned") {
    W_en_n = skew(transportRate(lla, v_nb_e, frame));
}
template <typename Float = double>
Mat3x3<Float> transportRateSkew(
        const Vec3<Float> &lla, const Vec3<Float> &v_nb_e, const std::string frame = "ned") {
    Mat3x3<Float> W_en_n;
    transportRateSkew(W_en_n, lla, v_nb_e, frame);
    return W_en_n;
}

//! === CORIOLISRATE ===
/// @brief      Coriolis effect perceived in the "NAV" frame
/// @param lla      Latitude, Longitude, Height [rad, rad, m]
/// @param v_nb_e   size 3 velocity vector in the 'NAV' coordinate system
/// @param frame    string representing the NAV-frame to rotate into
/// @param coriolis size 3 coriolis effect
/// @returns    Coriolis effect
template <typename Float = double>
void coriolisRate(
        Vec3<Float> &coriolis,
        const Vec3<Float> &lla,
        const Vec3<Float> &v_nb_e,
        const std::string frame = "ned") {
    Vec3<Float> w_ie_n = earthRate(lla(0), frame);
    Vec3<Float> w_en_n = transportRate(lla, v_nb_e, frame);
    coriolis = skew(w_en_n + 2.0 * w_ie_n) * v_nb_e;
}
template <typename Float = double>
Vec3<Float> coriolisRate(
        const Vec3<Float> &lla, const Vec3<Float> &v_nb_e, const std::string frame = "ned") {
    Vec3<Float> coriolis;
    coriolisRate(coriolis, lla, v_nb_e, frame);
    return coriolis;
}

//* ===== Gravity ============================================================================== *//

//! === SOMIGLIANA ===
/// @brief      Calculates the somilgiana model reference gravity
/// @param phi  Latitude [rad]
/// @param g0   Somigliana gravity
/// @returns Somgiliana gravity
template <typename Float = double>
void somigliana(Float &g0, const Float &phi) {
    Float sin_phi2 = std::sin(phi);
    sin_phi2 *= sin_phi2;
    g0 = 9.7803253359 *
         ((1.0 + 0.001931853 * sin_phi2) / std::sqrt(1.0 - WGS84_E2<Float> * sin_phi2));
}
template <typename Float = double>
Float somigliana(const Float &phi) {
    Float g0;
    somigliana(g0, phi);
    return g0;
}

//! === LOCALGRAVITY ===
/// @brief      Calculates gravity in the Local/NAV (ENU or NED) frame
/// @param lla      Latitude, Longitude, Height [rad, rad, m]
/// @param frame    string representing the NAV-frame to rotate into
/// @param g        size 3 Local/NAV frame gravity vector
/// @returns    Local/NAV frame gravity
template <typename Float = double>
void localGravity(Vec3<Float> &g, const Vec3<Float> &lla, const std::string frame = "ned") {
    Float sin_phi2 = std::sin(lla(0));
    sin_phi2 *= sin_phi2;
    Float g0 = 9.7803253359 *
               ((1.0 + 0.001931853 * sin_phi2) / std::sqrt(1.0 - WGS84_E2<Float> * sin_phi2));
    Float R02 = WGS84_R0<Float> * WGS84_R0<Float>;
    Float OMEGA2 = WGS84_OMEGA<Float> * WGS84_OMEGA<Float>;
    Float h2 = lla(2) * lla(2);
    if (frame == "ned") {
        // clang-format off
        g(0) = -8.08e-9 * lla(2) * std::sin(2.0 * lla(0));
        g(1) = 0.0;
        g(2) = g0 * (1.0 - 
                    (2.0 / WGS84_R0<Float>) *
                    (1.0 + WGS84_F<Float> * (1.0 - 2.0 * sin_phi2) + (OMEGA2 * R02 * WGS84_RP<Float> / GM<Float>)) *
                    lla(2) + (3.0 * h2 / R02));
        // clang-format on
    } else if (frame == "enu") {
        // clang-format off
        g(0) = 0.0;
        g(1) = -8.08e-9 * lla(2) * std::sin(2.0 * lla(0));
        g(2) = -g0 * (1.0 - 
                     (2.0 / WGS84_R0<Float>) *
                     (1.0 + WGS84_F<Float> * (1.0 - 2.0 * sin_phi2) + (OMEGA2 * R02 * WGS84_RP<Float> / GM<Float>)) *
                     lla(2) + (3.0 * h2 / R02));
        // clang-format on
    }
}
template <typename Float = double>
Vec3<Float> localGravity(const Vec3<Float> &lla, const std::string frame = "ned") {
    Vec3<Float> g;
    localGravity(g, lla, frame);
    return g;
}

//! === ECEFGRAVITY ===
/// @brief      Calculates gravity in the Earth-Centered-Earth-Fixed frame
/// @param xyz      ECEF position [m]
/// @param g        size 3 ECEF frame gravity vector
/// @param gamma    size 3 ECEF frame gravitational acceleration
/// @returns    ECEF frame gravity
template <typename Float = double>
void ecefGravity(Vec3<Float> &g, Vec3<Float> &gamma, const Vec3<Float> &xyz) {
    Float mag_r = xyz.norm();
    if (mag_r == 0.0) {
        g << 0.0, 0.0, 0.0;
        gamma << 0.0, 0.0, 0.0;
    } else {
        Float zeta = 5.0 * std::pow(xyz(2) / mag_r, 2.0);
        Float omega2 = WGS84_OMEGA<Float> * WGS84_OMEGA<Float>;
        Vec3<Float> v{1.0 - zeta, 1.0 - zeta, 3.0 - zeta};

        gamma = -GM<Float> / std::pow(mag_r, 3.0) *
                (xyz + 1.5 * J2<Float> * std::pow(WGS84_R0<Float> / mag_r, 2.0) * v * xyz);
        v << xyz(0) * omega2, xyz(1) * omega2, 0.0;
        g = gamma + v;
    }
}
template <typename Float = double>
void ecefGravity(Vec3<Float> &g, const Vec3<Float> &xyz) {
    Vec3<Float> gamma;
    ecefGravity(g, gamma, xyz);
}
template <typename Float = double>
Vec3<Float> ecefGravity(const Vec3<Float> &xyz) {
    Vec3<Float> g;
    ecefGravity(g, xyz);
    return g;
}

//* ===== 3d Ranging =========================================================================== *//

//! === CALCRANGE ===
/// @brief      Computes the range from user to satellite
/// @param sv_xyz   Satellite position [m]
/// @param user_xyz User position [m]
/// @param r        Calculated range [m]
/// @returns    Range
template <typename Float = double>
void calcRange(Float &r, const Vec3<Float> &sv_xyz, const Vec3<Float> &user_xyz) {
    r = (user_xyz - sv_xyz).norm();
}
template <typename Float = double>
Float calcRange(const Vec3<Float> &sv_xyz, const Vec3<Float> &user_xyz) {
    Float r;
    calcRange(r, sv_xyz, user_xyz);
    return r;
}

//! === CALCUNITVEC ===
/// @brief      Computes the unit vector from user to satellite
/// @param sv_xyz   Satellite position [m]
/// @param user_xyz User position [m]
/// @param u        Calculated unit vector [m]
/// @returns    Unit vector
template <typename Float = double>
void calcUnitVec(Vec3<Float> &u, const Vec3<Float> &sv_xyz, const Vec3<Float> &user_xyz) {
    Vec3<Float> dr = user_xyz - sv_xyz;
    u = dr / dr.norm();
}
template <typename Float = double>
Vec3<Float> calcUnitVec(const Vec3<Float> &sv_xyz, const Vec3<Float> &user_xyz) {
    Vec3<Float> u;
    calcUnitVec(u, sv_xyz, user_xyz);
    return u;
}

//! === CALCRANGEANDUNITVEC ===
/// @brief      Computes the range and unit vector from user to satellite
/// @param sv_xyz   Satellite position [m]
/// @param user_xyz User position [m]
/// @param r        Calculated range [m]
/// @param u        Calculated unit vector [m]
/// @returns    Range and unit vector
template <typename Float = double>
void calcRangeAndUnitVec(
        Float &r, Vec3<Float> &u, const Vec3<Float> &sv_xyz, const Vec3<Float> &user_xyz) {
    Vec3<Float> dr = user_xyz - sv_xyz;
    r = dr.norm();
    u = dr / r;
}

//! === CALCRANGERATE ===
/// @brief      Computes the range-rate between the user and satellite
/// @param u        Unit vector [m]
/// @param sv_vel   Satellite velocity [m/s]
/// @param user_vel User velocity [m/s]
/// @param rr       Calculated range-rate [m/s]
template <typename Float = double>
void calcRangeRate(
        Float &rr, const Vec3<Float> &u, const Vec3<Float> &sv_vel, const Vec3<Float> &user_vel) {
    Vec3<Float> dv = user_vel - sv_vel;
    rr = u(0) * dv(0) + u(1) * dv(1) + u(2) * dv(2);
}
template <typename Float = double>
Float calcRangeRate(const Vec3<Float> &u, const Vec3<Float> &sv_vel, const Vec3<Float> &user_vel) {
    Float rr;
    calcRangeRate(rr, u, sv_vel, user_vel);
    return rr;
}

}  // namespace navtools

#endif