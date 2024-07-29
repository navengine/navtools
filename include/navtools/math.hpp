/**
|========================================== math.hpp ==============================================|
|                                                                                                  |
|   @file     include/navtools/math.hpp                                                            |
|   @brief    Common mathematical operations.                                                      |
|   @ref      Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems          |
|               - (2013) Paul D. Groves                                                            |
|   @date     July 2024                                                                            |
|                                                                                                  |
|==================================================================================================|
*/

// TODO: move conversions

#ifndef NAVTOOLS_MATH_HPP
#define NAVTOOLS_MATH_HPP

#include <Eigen/Dense>
#include <cmath>

#include "navtools/constants.hpp"

namespace navtools {

//* ===== Skew-Symmmetric Matrix =============================================================== *//

//! === SKEW ===
/// @brief      Converts vector into its skew symmetric form
/// @param v    3x1 input vector
/// @param M    3x3 skew symmetric matrix
/// @returns    3x3 skew symmetric form of v
template <typename Float = double>
void Skew(Mat3x3<Float> &M, const Vec3<Float> &v) {
    M << 0.0, -v(2), v(1), v(2), 0.0, -v(0), -v(1), v(0), 0.0;
}
template <typename Float = double>
Mat3x3<Float> Skew(const Vec3<Float> &v) {
    Mat3x3<Float> M;
    Skew(M, v);
    return M;
}

//! === DESKEW ===
/// @brief      Converts skew symmetric form into its vector
/// @param M    3x3 input skew symmetric matrix
/// @param v    3x1 vector
/// @returns    3x1 vector form of M
template <typename Float = double>
void DeSkew(Vec3<Float> &v, const Mat3x3<Float> &M) {
    v << M(2, 1), M(0, 2), M(1, 0);
}
template <typename Float = double>
Vec3<Float> DeSkew(const Mat3x3<Float> &M) {
    Vec3<Float> v;
    DeSkew(v, M);
    return v;
}

//* ===== Sorting ============================================================================== *//

//! === MINPOSITIVEIDX ===
/// @brief find the index of the minimum value in vector greater than 0
template <typename Float = double>
size_t MinPositiveIdx(const std::vector<Float> &arr, const Float &set_point = 0.0) {
    size_t result = 0;
    Float tmp;
    Float max_val = std::numeric_limits<Float>::max();
    for (size_t i = 0; i < arr.size(); i++) {
        tmp = arr[i] - set_point;
        if ((0 < tmp) && (tmp < max_val)) {
            max_val = tmp;
            result = i;
        }
    }
    return result;
}

//* ===== Modulus Operations ===================================================================
//*//
//! === CIRCFMOD ===
/// @brief      Modulus of floating point number
/// @param x    user input and output
/// @param y    value to take modulus about
/// @returns    modulus of number
template <typename Float>
constexpr void CircMod(Float &x, const Float y) {
    x -= std::floor(x / y) * y;
}
// template <typename Float>
// constexpr Float CircMod(Float x, const Float y) {
//     return x - std::floor(x / y) * y;
// }

//* ===== Pi Wrapping ==========================================================================
//*//

//! === CIRCMOD2PI ===
/// @brief      Wraps angles from [0, 2*pi]
/// @param x    Inputs containing angles in radians
/// @returns    Wrapped/normalized angles [radians]
template <typename Float = double>
void CircMod2Pi(Float &x) {
    CircMod(x, TWO_PI<Float>);
}

//! === WRAPPITOPI ===
/// @brief      Wraps angles from [-pi, pi]
/// @param x    Inputs containing angles in radians
/// @returns    Wrapped/normalized angles [radians]
template <typename Float = double>
void WrapPiToPi(Float &x) {
    CircMod2Pi(x);
    if (x > PI<Float>) {
        x -= TWO_PI<Float>;
    }
}

//! === WRAPEULERANGLES ===
/// @brief      Auto wrap euler angles depending on pitch angle
/// @param x    Euler angles [radians]
/// @returns    Correctly wrapped Euler angles
template <typename Float = double>
void WrapEulerAngles(Vec3<Float> &x) {
    if (x(1) > HALF_PI<Float>) {
        x(0) += PI<Float>;
        x(1) = PI<Float> - x(1);
        x(2) += PI<Float>;
    } else if (x(1) < -HALF_PI<Float>) {
        x(0) += PI<Float>;
        x(1) = -PI<Float> - x(1);
        x(2) += PI<Float>;
    }
    WrapPiToPi(x(0));
    WrapPiToPi(x(2));
}

//! === RAD2DEG ===
/// @brief      Convert radians to degrees
/// @param  x   radians
/// @returns    degrees
template <typename Float>
void rad2deg(Float &x) {
    x *= RAD2DEG<Float>;
}
template <typename Float>
Float rad2deg(const Float x) {
    return x * RAD2DEG<Float>;
}

//! === DEG2RAD ===
/// @brief      Convert degrees to radians
/// @param  x   degrees
/// @returns    radiasn
template <typename Float>
void deg2rad(Float &x) {
    x *= DEG2RAD<Float>;
}
template <typename Float>
Float deg2rad(const Float x) {
    return x * DEG2RAD<Float>;
}

//* ===== Matrix Products/Normalization ========================================================
//*//

//! === QUATMAT ===
/// @brief      convert quaternion into its 4x4 matrix view
template <typename Float>
Mat4x4<Float> quatmat(const Vec4<Float> &q) {
    Mat4x4<Float> Q{
            {q(0), -q(1), -q(2), -q(3)},
            {q(1), q(0), -q(3), q(2)},
            {q(2), q(3), q(0), -q(1)},
            {q(3), -q(2), q(1), q(0)}};
    return Q;
}

//! === QUATDOT ===
/// @brief      quaternion product
template <typename Float = double>
Vec4<Float> quatdot(const Vec4<Float> &p, const Vec4<Float> &q) {
    // p o q
    Vec4<Float> r{
            p(0) * q(0) - p(1) * q(1) - p(2) * q(2) - p(3) * q(3),
            p(0) * q(1) + p(1) * q(0) + p(2) * q(3) - p(3) * q(2),
            p(0) * q(2) - p(1) * q(3) + p(2) * q(0) + p(3) * q(1),
            p(0) * q(3) + p(1) * q(2) - p(2) * q(1) + p(3) * q(0)};
    return r;
}
template <typename Float = double>
Vec4<Float> quatdot(const Vec3<Float> &a, const Vec4<Float> &q) {
    // a o q
    Vec4<Float> r{
            -a(0) * q(1) - a(1) * q(2) - a(2) * q(3),
            a(0) * q(0) + a(1) * q(3) - a(2) * q(2),
            -a(0) * q(3) + a(1) * q(0) + a(2) * q(1),
            a(0) * q(2) - a(1) * q(1) + a(2) * q(0)};
    return r;
}
template <typename Float = double>
Vec4<Float> quatdot(const Vec4<Float> &p, const Vec3<Float> &a) {
    // p o a
    Vec4<Float> r{
            -p(1) * a(0) - p(2) * a(1) - p(3) * a(2),
            p(0) * a(0) + p(2) * a(2) - p(3) * a(1),
            p(0) * a(1) - p(1) * a(2) + p(3) * a(0),
            p(0) * a(2) + p(1) * a(1) - p(2) * a(0)};
    return r;
}

//! === QUATCONJ/QUATINV ===
/// @brief      quaternion conjugate/inverse
template <typename Float = double>
void quatconj(Vec4<Float> &q) {
    q(1) *= -1.0;
    q(2) *= -1.0;
    q(3) *= -1.0;
}
template <typename Float = double>
void quatinv(Vec4<Float> &q) {
    q(1) *= -1.0;
    q(2) *= -1.0;
    q(3) *= -1.0;
}

//! === QUATNORM ===
/// @brief      quaternion normalization
template <typename Float = double>
void quatnorm(Vec4<Float> &q) {
    q /= q.norm();
}

//! === DCMNORM ===
/// @brief      DCM normalization
template <typename Float = double>
void dcmnorm(Mat3x3<Float> &R) {
    // Groves 5.79, 5.80
    Vec3<Float> c1 = R.col(0);
    Vec3<Float> c2 = R.col(1);
    Vec3<Float> c3 = R.col(2);

    Float c1c2 = R(0, 0) * R(0, 1) + R(1, 0) * R(1, 1) + R(2, 0) * R(2, 1);
    Float c1c3 = R(0, 0) * R(0, 2) + R(1, 0) * R(1, 2) + R(2, 0) * R(2, 2);
    Float c2c3 = R(0, 1) * R(0, 2) + R(1, 1) * R(1, 2) + R(2, 1) * R(2, 2);

    c1 -= 0.5 * (c1c2 * c2 + c1c3 * c3);
    c2 -= 0.5 * (c1c2 * c1 + c2c3 * c3);
    c3 -= 0.5 * (c1c3 * c1 + c2c3 * c2);
    c1 /= c1.norm();
    c2 /= c2.norm();
    c3 /= c3.norm();

    R.col(0) = c1;
    R.col(1) = c2;
    R.col(2) = c3;
}

//* ===== Matrix Exponential ===================================================================
//*//

//! === RODRIGUES ===
/// @brief      Rodrigues formula for the approximation of a matrix exponential
/// @param vec          size 3 vector
/// @param vec_norm     2-norm of vec
/// @returns    matrix exponential
template <typename Float = double>
Mat3x3<Float> Rodrigues(const Vec3<Float> &vec) {
    Float vec_norm = vec.norm();
    Mat3x3<Float> skew_sym = Skew(vec / vec_norm);
    return Eigen::Matrix<Float, 3, 3>::Identity() + (std::sin(vec_norm) * skew_sym) +
           ((1.0 - std::cos(vec_norm)) * skew_sym * skew_sym);
}
template <typename Float = double>
Mat3x3<Float> Rodrigues(const Vec3<Float> &vec, const Float &vec_norm) {
    Mat3x3<Float> skew_sym = Skew(vec / vec_norm);
    return Eigen::Matrix<Float, 3, 3>::Identity() + (std::sin(vec_norm) * skew_sym) +
           ((1.0 - std::cos(vec_norm)) * skew_sym * skew_sym);
}

//! === RODRIGUES4 ===
/// @brief      Rodrigues formula for the 4-th order approximation of a matrix exponential
/// @param vec          size 3 vector
/// @param vec_norm     2-norm of vec
/// @returns    matrix exponential
template <typename Float = double>
Mat3x3<Float> Rodrigues4(const Vec3<Float> &vec) {
    Float vec_norm = vec.norm();
    Mat3x3<Float> skew_sym = Skew(vec);
    Float norm_squared = vec_norm * vec_norm;
    return Eigen::Matrix<Float, 3, 3>::Identity() + ((1.0 - (norm_squared / 6.0)) * skew_sym) +
           ((0.5 - (norm_squared / 24.0)) * skew_sym * skew_sym);
}
template <typename Float = double>
Mat3x3<Float> Rodrigues4(const Vec3<Float> &vec, const Float &vec_norm) {
    Mat3x3<Float> skew_sym = Skew(vec);
    Float norm_squared = vec_norm * vec_norm;
    return Eigen::Matrix<Float, 3, 3>::Identity() + ((1.0 - (norm_squared / 6.0)) * skew_sym) +
           ((0.5 - (norm_squared / 24.0)) * skew_sym * skew_sym);
}

//! === VEC2EXPM ===
/// @brief      Converts vector into its matrix exponential approximation form
/// @param vec  size 3 vector
/// @returns    3x3 matrix exponential
template <typename Float = double>
Mat3x3<Float> vec2expm(const Vec3<Float> &vec) {
    Float vec_norm = vec.norm();
    if (vec_norm < 0.02) {
        return Rodrigues4(vec, vec_norm);
    } else {
        return Rodrigues(vec, vec_norm);
    }
}

//! === EXPM2VEC ===
/// @brief      Converts matrix exponential approximation into its vector form
/// @param mat  3x3 matrix exponential
/// @returns    size 3 vector
template <typename Float = double>
Vec3<Float> expm2vec(const Mat3x3<Float> &mat) {
    Float phi = std::acos((mat.trace() - 1.0) / 2.0);
    if (phi == 0.0) {
        return Eigen::Vector<Float, 3>::Zero();
    }
    return phi * DeSkew(mat - mat.transpose()) / (2.0 * std::sin(phi));
}

//* ===== Signal to Noise ======================================================================
//*//
// https://insidegnss.com/wp-content/uploads/2018/01/novdec10-Solutions.pdf

//! === WATT2DB ===
/// @brief      Convert unit of power to decibels
/// @param w    Input power
/// @param db   Output decibels
/// @returns    Decibels
template <typename Float = double>
void watt2db(Float &db, const Float &w) {
    db = 10.0 * std::log10(w);
}
template <typename Float = double>
Float watt2db(const Float &w) {
    Float db;
    watt2db(db, w);
    return db;
}

//! === DB2WATT ===
/// @brief      Convert decibels to unit of power
/// @param db   Input decibels
/// @param w    Output power
/// @returns    unit of power
template <typename Float = double>
void db2watt(Float &w, const Float &db) {
    w = std::pow(10.0, db / 10.0);
}
template <typename Float = double>
Float db2watt(const Float &db) {
    Float w;
    db2watt(w, db);
    return w;
}

//! === VOLT2DB ===
/// @brief      Convert unit of sqrt-power (volts) to decibels
/// @param v    Input sqrt-power
/// @param dB   Output decibels
/// @returns    Decibels
template <typename Float = double>
void volt2db(Float &db, const Float &v) {
    db = 20.0 * std::log10(v);
}
template <typename Float = double>
Float volt2db(const Float &v) {
    Float db;
    volt2db(db, v);
    return db;
}

//! === DB2VOLT ===
/// @brief      Convert decibels to unit of sqrt-power
/// @param db   Input decibels
/// @param v    Output sqrt-power
/// @returns    unit of sqrt-power
template <typename Float = double>
void db2volt(Float &v, const Float db) {
    v = std::pow(10.0, db / 20.0);
}
template <typename Float = double>
Float db2volt(const Float &db) {
    Float v;
    db2volt(v, db);
    return v;
}

//! === CN02SNR ===
/// @brief      Convert Carrier-to-Noise Ratio into raw Signal-to-Noise Ratio
/// @param cn0      Carrier-to-noise density ratio [dB/Hz]
/// @param fe_bw    Receiver front-end bandwidth [Hz]
/// @param T        Additional noise figure/loses from the given temperature [K]
/// @param eta      Additional noise figure/loses [dB]
/// @param snr      Signal to noise ratio [dB]
template <typename Float = double>
void cn02snr(
        Float &snr,
        const Float &cn0,
        const Float &fe_bw,
        const Float &T = 290.0,
        const Float &eta = 0.0) {
    snr = cn0 - watt2db(fe_bw) - watt2db(1.0 + T / 290.0) - eta;
}
template <typename Float = double>
Float cn02snr(
        const Float &cn0, const Float &fe_bw, const Float &T = 290.0, const Float &eta = 0.0) {
    Float snr;
    cn02snr(snr, cn0, fe_bw, T, eta);
    return snr;
}

//! === SNR2CN0 ===
/// @brief      Convert Carrier-to-Noise Ratio into raw Signal-to-Noise Ratio
/// @param snr      Signal to noise ratio [dB]
/// @param fe_bw    Receiver front-end bandwidth [Hz]
/// @param T        Additional noise figure/loses from the given temperature [K]
/// @param eta      Additional noise figure/loses [dB]
/// @param cn0      Carrier-to-noise density ratio [dB/Hz]
template <typename Float = double>
void snr2cn0(
        Float &cn0,
        const Float &snr,
        const Float &fe_bw,
        const Float &T = 290.0,
        const Float &eta = 0.0) {
    cn0 = snr + watt2db(fe_bw) + watt2db(1.0 + T / 290.0) + eta;
}
template <typename Float = double>
Float snr2cn0(
        const Float &snr, const Float &fe_bw, const Float &T = 290.0, const Float &eta = 0.0) {
    Float cn0;
    snr2cn0(cn0, snr, fe_bw, T, eta);
    return cn0;
}

}  // namespace navtools

#endif
