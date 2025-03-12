/**
 * *attitude.hpp*
 *
 * =======  ========================================================================================
 * @file    include/navtools/attitude.hpp
 * @brief   Attitude representations and conversions between them.
 * @author  Daniel Sturdivant, Blake Baker
 * @ref     Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems
 *            - (2013) Paul D. Groves
 * @date    March 2025
 * =======  ========================================================================================
 */

#ifndef NAVTOOLS_ATTITUDE_HPP
#define NAVTOOLS_ATTITUDE_HPP

#include <Eigen/Dense>

#include "navtools/constants.hpp"

namespace navtools {

//! === EULER2QUAT ===
/// @brief      Converts euler angles (roll-pitch-yaw) to corresponding BODY-to-NAV quaternion
/// @param e        size 3 RPY euler angles [radians]
/// @param frame    string representing the NAV-frame to rotate into
/// @param q        size 4 quaternion rotation
/// @returns    4x1 ZYX quaternion
template <typename T = double>
void euler2quat(
    Eigen::Ref<Eigen::Vector<T, 4>> q,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &e,
    const bool IsNed = true) {
  Eigen::Vector<T, 3> x = e / 2.0;
  T Sr = std::sin(x(0));
  T Sp = std::sin(x(1));
  T Cr = std::cos(x(0));
  T Cp = std::cos(x(1));
  T a = Sp * Cr;
  T b = Cp * Sr;
  T c = Sp * Sr;
  T d = Cp * Cr;
  if (IsNed) {
    // ned frame
    T Cy = std::cos(x(2));
    T Sy = std::sin(x(2));
    q << Cy * d + Sy * c, Cy * b - Sy * a, Cy * a + Sy * b, Sy * d - Cy * c;
  } else {
    // enu frame
    T Cy_pio4 = std::cos(x(2) + PI<T> / 4.0);
    T Sy_pio4 = std::sin(x(2) + PI<T> / 4.0);
    q << -a * Cy_pio4 - b * Sy_pio4, -c * Cy_pio4 + d * Sy_pio4, c * Sy_pio4 + d * Cy_pio4,
        a * Sy_pio4 - b * Cy_pio4;
  }
}
template <typename T = double>
Eigen::Vector<T, 4> euler2quat(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &e, const bool IsNed = true) {
  Eigen::Vector<T, 4> q;
  euler2quat<T>(q, e, IsNed);
  return q;
}

//! === EULER2DCM ===
/// @brief      Converts euler angles (roll-pitch-yaw) to corresponding BODY-to-NAV DCM
/// @param e        size 3 RPY euler angles [radians]
/// @param frame    string representing the NAV-frame to rotate into
/// @param R        3x3 DCM
/// @returns    3x3 ZYX DCM
template <typename T = double>
void euler2dcm(
    Eigen::Ref<Eigen::Matrix<T, 3, 3>> R,
    const Eigen::Ref<const Eigen::Vector<T, 3>> &e,
    const bool IsNed = true) {
  T Sr = std::sin(e(0));
  T Sp = std::sin(e(1));
  T Sy = std::sin(e(2));
  T Cr = std::cos(e(0));
  T Cp = std::cos(e(1));
  T Cy = std::cos(e(2));
  if (IsNed) {
    // clang-format off
        R << Cp * Cy, Sr * Sp * Cy - Cr * Sy, Cr * Sp * Cy + Sr * Sy,
             Cp * Sy, Sr * Sp * Sy + Cr * Cy, Cr * Sp * Sy - Cy * Sr,
                 -Sp,                Cp * Sr,                Cp * Cr;
    // clang-format on
  } else {
    // clang-format off
        R << Cp * Sy, Sr * Sp * Sy + Cr * Cy, Cr * Sp * Sy - Sr * Cy,
             Cp * Cy, Sr * Sp * Cy - Cr * Sy, Cr * Sp * Cy + Sr * Sy,
                  Sp,               -Sr * Cp,               -Cr * Cp;
    // clang-format on
  }
}
template <typename T = double>
Eigen::Matrix<T, 3, 3> euler2dcm(
    const Eigen::Ref<const Eigen::Vector<T, 3>> &e, const bool IsNed = true) {
  Eigen::Matrix<T, 3, 3> R;
  euler2dcm<T>(R, e, IsNed);
  return R;
}

//! === QUAT2EULER ===
/// @brief      Converts BODY-to-NAV quaternion to corresponding euler angles (roll-pitch-yaw)
/// @param q        size 4 quaternion rotation
/// @param frame    string representing the NAV-frame to rotate into
/// @param e        size 3 RPY euler angles [radians]
/// @returns    3x1 RPY euler angles [radians]
template <typename T = double>
void quat2euler(
    Eigen::Ref<Eigen::Vector<T, 3>> e,
    const Eigen::Ref<const Eigen::Vector<T, 4>> &q,
    const bool IsNed = true) {
  T w = q(0);
  T x = q(1);
  T y = q(2);
  T z = q(3);
  T x2 = x * x;
  T y2 = y * y;
  T z2 = z * z;
  T xw = x * w;
  T xy = x * y;
  T xz = x * z;
  T yw = y * w;
  T yz = y * z;
  T zw = z * w;
  if (IsNed) {
    e(0) = std::atan2(2.0 * (xw + yz), 1.0 - 2.0 * (x2 + y2));
    e(1) = std::asin(2.0 * (yw - xz));
    e(2) = std::atan2(2.0 * (zw + xy), 1.0 - 2.0 * (y2 + z2));
  } else {
    e(0) = PI<T> + std::atan2(2.0 * (xw + yz), 1.0 - 2.0 * (x2 + y2));
    e(1) = -std::asin(2.0 * (yw - xz));
    e(2) = HALF_PI<T> - std::atan2(2.0 * (zw + xy), 1.0 - 2.0 * (y2 + z2));
  }
}
template <typename T = double>
Eigen::Vector<T, 3> quat2euler(
    const Eigen::Ref<const Eigen::Vector<T, 4>> &q, const bool IsNed = true) {
  Eigen::Vector<T, 3> e;
  quat2euler<T>(e, q, IsNed);
  return e;
}

//! === QUAT2DCM ===
/// @brief      Converts BODY-to-NAV quaternion to corresponding BODY-to-NAV DCM
/// @param q        size 4 quaternion rotation
/// @param R        3x3 DCM
/// @returns    3x3 ZYX DCM
template <typename T = double>
void quat2dcm(
    Eigen::Ref<Eigen::Matrix<T, 3, 3>> R, const Eigen::Ref<const Eigen::Vector<T, 4>> &q) {
  T w = q(0);
  T x = q(1);
  T y = q(2);
  T z = q(3);
  T w2 = w * w;
  T x2 = x * x;
  T y2 = y * y;
  T z2 = z * z;
  T xw = x * w;
  T xy = x * y;
  T xz = x * z;
  T yw = y * w;
  T yz = y * z;
  T zw = z * w;
  // clang-format off
  R << w2 + x2 - y2 - z2,   2.0 * (xy - zw),   2.0 * (yw + xz),
         2.0 * (zw + xy), w2 - x2 + y2 - z2,   2.0 * (yz - xw),
         2.0 * (xz - yw),   2.0 * (yz + xw), w2 - x2 - y2 + z2;
  // R = {{w*w + x*x - y*y - z*z,       2.0*(x*y + z*w),       2.0*(x*z - y*w)},
  //      {      2.0*(x*y - z*w), w*w - x*x + y*y - z*z,       2.0*(y*z + x*w)},
  //      {      2.0*(x*z + y*w),       2.0*(y*z - x*w), w*w - x*x - y*y + z*z}};
  // clang-format on
}
template <typename T = double>
Eigen::Matrix<T, 3, 3> quat2dcm(const Eigen::Ref<const Eigen::Vector<T, 4>> &q) {
  Eigen::Matrix<T, 3, 3> R;
  quat2dcm<T>(R, q);
  return R;
}

//! === DCM2EULER ===
/// @brief      Converts BODY-to-NAV DCM to corresponding euler angles (roll-pitch-yaw)
/// @param R        3x3 DCM
/// @param frame    string representing the NAV-frame to rotate into
/// @param e        size 3 RPY euler angles [radians]
/// @returns    3x1 RPY euler angles [radians]
template <typename T = double>
void dcm2euler(
    Eigen::Ref<Eigen::Vector<T, 3>> e,
    const Eigen::Ref<const Eigen::Matrix<T, 3, 3>> &R,
    const bool IsNed = true) {
  if (IsNed) {
    e(0) = std::atan2(R(2, 1), R(2, 2));
    e(1) = -std::asin(R(2, 0));
    e(2) = std::atan2(R(1, 0), R(0, 0));
  } else {
    e(0) = std::atan2(-R(2, 1), -R(2, 2));
    e(1) = std::asin(R(2, 0));
    e(2) = std::atan2(R(0, 0), R(1, 0));
  }
}
template <typename T = double>
Eigen::Vector<T, 3> dcm2euler(
    const Eigen::Ref<const Eigen::Matrix<T, 3, 3>> &R, const bool IsNed = true) {
  Eigen::Vector<T, 3> e;
  dcm2euler<T>(e, R, IsNed);
  return e;
}

//! === DCM2QUAT ===
/// @brief      Converts BODY-to-NAV DCM to corresponding BODY-to-NAV quaternion
/// @param R        3x3 DCM
/// @param q        size 4 quaternion rotation
/// @returns    4x1 ZYX quaternion
template <typename T = double>
void dcm2quat(
    Eigen::Ref<Eigen::Vector<T, 4>> q, const Eigen::Ref<const Eigen::Matrix<T, 3, 3>> &R) {
  T q_w = std::sqrt(1.0 + R(0, 0) + R(1, 1) + R(2, 2)) / 2.0;
  if (q_w > 0.01) {
    T q_w_4 = 4.0 * q_w;
    q(0) = q_w;
    q(1) = (R(2, 1) - R(1, 2)) / q_w_4;
    q(2) = (R(0, 2) - R(2, 0)) / q_w_4;
    q(3) = (R(1, 0) - R(0, 1)) / q_w_4;
  } else {
    q(1) = std::sqrt(1.0 + R(0, 0) - R(1, 1) - R(2, 2)) / 2.0;
    T q_x_4 = 4.0 * q(1);
    q(0) = (R(2, 1) - R(1, 2)) / q_x_4;
    q(2) = (R(0, 1) + R(1, 0)) / q_x_4;
    q(3) = (R(0, 2) + R(2, 0)) / q_x_4;
  }
}
template <typename T = double>
Eigen::Vector<T, 4> dcm2quat(const Eigen::Ref<const Eigen::Matrix<T, 3, 3>> &R) {
  Eigen::Vector<T, 4> q;
  dcm2quat<T>(q, R);
  return q;
}

//* ===== Single Axes Rotation Matrices ======================================================== *//

//! === ROTX ===
/// @brief      Converts euler angle about x-axis into DCM rotation about x-axis
/// @param C    3x3 DCM
/// @param x    euler angle ]radians]
/// @returns    3x3 x-axis DCM rotation
template <typename T = double>
void RotX(Eigen::Ref<Eigen::Matrix<T, 3, 3>> C, const T &x) {
  T sx = std::sin(x);
  T cx = std::cos(x);
  C << 1.0, 0.0, 0.0, 0.0, cx, -sx, 0.0, sx, cx;
}
template <typename T = double>
Eigen::Matrix<T, 3, 3> RotX(const T &x) {
  Eigen::Matrix<T, 3, 3> C;
  RotX<T>(C, x);
  return C;
}

//! === ROTY ===
/// @brief      Converts euler angle about y-axis into DCM rotation about y-axis
/// @param C    3x3 DCM
/// @param y    euler angle ]radians]
/// @returns    3x3 y-axis DCM rotation
template <typename T = double>
void RotY(Eigen::Ref<Eigen::Matrix<T, 3, 3>> C, const T &y) {
  T sy = std::sin(y);
  T cy = std::cos(y);
  C << cy, 0.0, sy, 0.0, 1.0, 0.0, -sy, 0.0, cy;
}
template <typename T = double>
Eigen::Matrix<T, 3, 3> RotY(const T &y) {
  Eigen::Matrix<T, 3, 3> C;
  RotY<T>(C, y);
  return C;
}

//! === ROTZ ===
/// @brief      Converts euler angle about z-axis into DCM rotation about z-axis
/// @param C    3x3 DCM
/// @param z    euler angle ]radians]
/// @returns    3x3 z-axis DCM rotation
template <typename T = double>
void RotZ(Eigen::Ref<Eigen::Matrix<T, 3, 3>> C, const T &z) {
  T sz = std::sin(z);
  T cz = std::cos(z);
  C << cz, -sz, 0.0, sz, cz, 0.0, 0.0, 0.0, 1.0;
}
template <typename T = double>
Eigen::Matrix<T, 3, 3> RotZ(const T &z) {
  Eigen::Matrix<T, 3, 3> C;
  RotZ<T>(C, z);
  return C;
}

}  // namespace navtools

#endif
