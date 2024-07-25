/**
|========================================= attitude.hpp ===========================================|
|                                                                                                  |
|   @file     include/navtools/attitude.hpp                                                        |
|   @brief    Attitude representations and conversions between them.                               |
|   @ref      Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems          |
|               - (2013) Paul D. Groves                                                            |
|   @date     July 2024                                                                            |
|                                                                                                  |
|==================================================================================================|
*/

#ifndef NAVTOOLS_ATTITUDE_HPP
#define NAVTOOLS_ATTITUDE_HPP

#include <Eigen/Dense>

#include "navtools/constants.hpp"
#include "navtools/types.hpp"

namespace navtools {

//! === EULER2QUAT ===
/// @brief      Converts euler angles (roll-pitch-yaw) to corresponding BODY-to-NAV quaternion
/// @param e        size 3 RPY euler angles [radians]
/// @param frame    string representing the NAV-frame to rotate into
/// @param q        size 4 quaternion rotation
/// @returns    4x1 ZYX quaternion
template <typename Float = double>
void euler2quat(Vec4<Float> &q, const Vec3<Float> &e, const std::string frame = "ned") {
    Vec3<Float> x = e / 2.0;
    Float Sr = std::sin(x(0));
    Float Sp = std::sin(x(1));
    Float Cr = std::cos(x(0));
    Float Cp = std::cos(x(1));
    Float a = Sp * Cr;
    Float b = Cp * Sr;
    Float c = Sp * Sr;
    Float d = Cp * Cr;
    if (frame == "ned") {
        Float Cy = std::cos(x(2));
        Float Sy = std::sin(x(2));
        q << Cy * d + Sy * c, Cy * b - Sy * a, Cy * a + Sy * b, Sy * d - Cy * c;
    } else if (frame == "enu") {
        Float Cy_pio4 = std::cos(x(2) + PI<Float> / 4.0);
        Float Sy_pio4 = std::sin(x(2) + PI<Float> / 4.0);
        q << -a * Cy_pio4 - b * Sy_pio4, -c * Cy_pio4 + d * Sy_pio4, c * Sy_pio4 + d * Cy_pio4,
                a * Sy_pio4 - b * Cy_pio4;
    }
}
template <typename Float = double>
Vec4<Float> euler2quat(const Vec3<Float> &e, const std::string frame = "ned") {
    Vec4<Float> q;
    euler2quat(q, e, frame);
    return q;
}

//! === EULER2DCM ===
/// @brief      Converts euler angles (roll-pitch-yaw) to corresponding BODY-to-NAV DCM
/// @param e        size 3 RPY euler angles [radians]
/// @param frame    string representing the NAV-frame to rotate into
/// @param R        3x3 DCM
/// @returns    3x3 ZYX DCM
template <typename Float = double>
void euler2dcm(Mat3x3<Float> &R, const Vec3<Float> &e, const std::string frame = "ned") {
    Float Sr = std::sin(e(0));
    Float Sp = std::sin(e(1));
    Float Sy = std::sin(e(2));
    Float Cr = std::cos(e(0));
    Float Cp = std::cos(e(1));
    Float Cy = std::cos(e(2));
    if (frame == "ned") {
        // clang-format off
        R << Cp * Cy, Sr * Sp * Cy - Cr * Sy, Cr * Sp * Cy + Sr * Sy,
             Cp * Sy, Sr * Sp * Sy + Cr * Cy, Cr * Sp * Sy - Cy * Sr,
                 -Sp,                Cp * Sr,                Cp * Cr;
        // clang-format on
    } else if (frame == "enu") {
        // clang-format off
        R << Cp * Sy, Sr * Sp * Sy + Cr * Cy, Cr * Sp * Sy - Sr * Cy,
             Cp * Cy, Sr * Sp * Cy - Cr * Sy, Cr * Sp * Cy + Sr * Sy,
                  Sp,               -Sr * Cp,               -Cr * Cp;
        // clang-format on
    }
}
template <typename Float = double>
Mat3x3<Float> euler2dcm(const Vec3<Float> &e, const std::string frame = "ned") {
    Mat3x3<Float> R;
    euler2dcm(R, e, frame);
    return R;
}

//! === QUAT2EULER ===
/// @brief      Converts BODY-to-NAV quaternion to corresponding euler angles (roll-pitch-yaw)
/// @param q        size 4 quaternion rotation
/// @param frame    string representing the NAV-frame to rotate into
/// @param e        size 3 RPY euler angles [radians]
/// @returns    3x1 RPY euler angles [radians]
template <typename Float = double>
void quat2euler(Vec3<Float> &e, const Vec4<Float> &q, const std::string frame = "ned") {
    Float w = q(0);
    Float x = q(1);
    Float y = q(2);
    Float z = q(3);
    if (frame == "ned") {
        e(0) = std::atan2(2.0 * (w * x + y * z), 1.0 - 2.0 * (x * x + y * y));
        e(1) = std::asin(2.0 * (w * y - x * z));
        e(2) = std::atan2(2.0 * (w * z + x * y), 1.0 - 2.0 * (y * y + z * z));
    } else if (frame == "enu") {
        e(0) = PI<Float> + std::atan2(2.0 * (w * x + y * z), 1.0 - 2.0 * (x * x + y * y));
        e(1) = -std::asin(2.0 * (w * y - x * z));
        e(2) = HALF_PI<Float> - std::atan2(2.0 * (w * z + x * y), 1.0 - 2.0 * (y * y + z * z));
    }
}
template <typename Float = double>
Vec3<Float> quat2euler(const Vec4<Float> &q, const std::string frame = "ned") {
    Vec3<Float> e;
    quat2euler(e, q, frame);
    return e;
}

//! === QUAT2DCM ===
/// @brief      Converts BODY-to-NAV quaternion to corresponding BODY-to-NAV DCM
/// @param q        size 4 quaternion rotation
/// @param R        3x3 DCM
/// @returns    3x3 ZYX DCM
template <typename Float = double>
void quat2dcm(Mat3x3<Float> &R, const Vec4<Float> &q) {
    // clang-format off
    Float w = q(0);
    Float x = q(1);
    Float y = q(2);
    Float z = q(3);
    R << w * w + x * x - y * y - z * z,         2.0 * (x * y - w * z),         2.0 * (w * y + x * z),
                 2.0 * (w * z + x * y), w * w - x * x + y * y - z * z,         2.0 * (y * z - w * x),
                 2.0 * (x * z - w * y),         2.0 * (y * z + w * x), w * w - x * x - y * y + z * z;
    // R = {{w*w + x*x - y*y - z*z,       2.0*(x*y + z*w),       2.0*(x*z - y*w)},
    //      {      2.0*(x*y - z*w), w*w - x*x + y*y - z*z,       2.0*(y*z + x*w)},
    //      {      2.0*(x*z + y*w),       2.0*(y*z - x*w), w*w - x*x - y*y + z*z}};
    // clang-format on
}
template <typename Float = double>
Mat3x3<Float> quat2dcm(const Vec4<Float> &q) {
    Mat3x3<Float> R;
    quat2dcm(R, q);
    return R;
}

//! === DCM2EULER ===
/// @brief      Converts BODY-to-NAV DCM to corresponding euler angles (roll-pitch-yaw)
/// @param R        3x3 DCM
/// @param frame    string representing the NAV-frame to rotate into
/// @param e        size 3 RPY euler angles [radians]
/// @returns    3x1 RPY euler angles [radians]
template <typename Float = double>
void dcm2euler(Vec3<Float> &e, const Mat3x3<Float> &R, const std::string frame = "ned") {
    if (frame == "ned") {
        e(0) = std::atan2(R(2, 1), R(2, 2));
        e(1) = -std::asin(R(2, 0));
        e(2) = std::atan2(R(1, 0), R(0, 0));
    } else if (frame == "enu") {
        e(0) = std::atan2(-R(2, 1), -R(2, 2));
        e(1) = std::asin(R(2, 0));
        e(2) = std::atan2(R(0, 0), R(1, 0));
    }
}
template <typename Float = double>
Vec3<Float> dcm2euler(const Mat3x3<Float> &R, const std::string frame = "ned") {
    Vec3<Float> e;
    dcm2euler(e, R, frame);
    return e;
}

//! === DCM2QUAT ===
/// @brief      Converts BODY-to-NAV DCM to corresponding BODY-to-NAV quaternion
/// @param R        3x3 DCM
/// @param q        size 4 quaternion rotation
/// @returns    4x1 ZYX quaternion
template <typename Float = double>
void dcm2quat(Vec4<Float> &q, const Mat3x3<Float> &R) {
    Float q_w = std::sqrt(1.0 + R(0, 0) + R(1, 1) + R(2, 2)) / 2.0;
    if (q_w > 0.01) {
        Float q_w_4 = 4.0 * q_w;
        q(0) = q_w;
        q(1) = (R(2, 1) - R(1, 2)) / q_w_4;
        q(2) = (R(0, 2) - R(2, 0)) / q_w_4;
        q(3) = (R(1, 0) - R(0, 1)) / q_w_4;
    } else {
        q(1) = std::sqrt(1.0 + R(0, 0) - R(1, 1) - R(2, 2)) / 2.0;
        Float q_x_4 = 4.0 * q(1);
        q(0) = (R(2, 1) - R(1, 2)) / q_x_4;
        q(2) = (R(0, 1) + R(1, 0)) / q_x_4;
        q(3) = (R(0, 2) + R(2, 0)) / q_x_4;
    }
}
template <typename Float = double>
Vec4<Float> dcm2quat(const Mat3x3<Float> &R) {
    Vec4<Float> q;
    dcm2quat(q, R);
    return q;
}

//* ===== Single Axes Rotation Matrices ======================================================== *//

//! === ROTX ===
/// @brief      Converts euler angle about x-axis into DCM rotation about x-axis
/// @param C    3x3 DCM
/// @param x    euler angle ]radians]
/// @returns    3x3 x-axis DCM rotation
template <typename Float = double>
void rotX(Mat3x3<Float> &C, const Float &x) {
    Float sx = std::sin(x);
    Float cx = std::cos(x);
    C << 1.0, 0.0, 0.0, 0.0, cx, -sx, 0.0, sx, cx;
}
template <typename Float = double>
Mat3x3<Float> rotX(const Float &x) {
    Mat3x3<Float> C;
    rotX(C, x);
    return C;
}

//! === ROTY ===
/// @brief      Converts euler angle about y-axis into DCM rotation about y-axis
/// @param C    3x3 DCM
/// @param y    euler angle ]radians]
/// @returns    3x3 y-axis DCM rotation
template <typename Float = double>
void rotY(Mat3x3<Float> &C, const Float &y) {
    Float sy = std::sin(y);
    Float cy = std::cos(y);
    C = {{cy, 0.0, sy}, {0.0, 1.0, 0.0}, {-sy, 0.0, cy}};
}
template <typename Float = double>
Mat3x3<Float> rotY(const Float &y) {
    Mat3x3<Float> C;
    rotY(C, y);
    return C;
}

//! === ROTZ ===
/// @brief      Converts euler angle about z-axis into DCM rotation about z-axis
/// @param C    3x3 DCM
/// @param z    euler angle ]radians]
/// @returns    3x3 z-axis DCM rotation
template <typename Float = double>
void rotZ(Mat3x3<Float> &C, const Float &z) {
    Float sz = std::sin(z);
    Float cz = std::cos(z);
    C = {{cz, -sz, 0.0}, {sz, cz, 0.0}, {0.0, 0.0, 1.0}};
}
template <typename Float = double>
Mat3x3<Float> rotZ(const Float &z) {
    Mat3x3<Float> C;
    rotZ(C, z);
    return C;
}

}  // namespace navtools

#endif
