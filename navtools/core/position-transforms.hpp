#ifndef NAVTOOLS_CORE_POSITION_TRANSFORMS_HPP
#define NAVTOOLS_CORE_POSITION_TRANSFORMS_HPP

#include "navtools/core/constants.hpp"
#include "navtools/core/coordinate-dcms.hpp"
#include "navtools/core/macros.hpp"
#include <Eigen/Dense>

namespace nt {

//! ============================================================================================ !//

/**
 * *=== ecef2lla ===*
 * @brief converts Earth-centered-earth-fixed frame to Latitude-Longitude-Altitude position
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param xyz 3d vector of cartesian ECEF position (m,m,m)
 * @param lla 3d vector of ellipsoidal LLA position (rad,rad,m)
 */
template <typename Derived1, typename Derived2>
inline void ecef2lla(const Eigen::DenseBase<Derived1> &xyz, Eigen::DenseBase<Derived2> &lla) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, xyz, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, lla, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  using Scalar = typename Derived1::Scalar;

  const auto &in = xyz.derived();
  const Scalar &x = in(0);
  const Scalar &y = in(1);
  const Scalar &z = in(2);

  Scalar sign_z = std::copysign(1.0, z);
  Scalar sqrt_1_e2 = std::sqrt(1.0 - WGS84_E2<Scalar>);
  Scalar beta = std::sqrt(x * x + y * y);  // (Groves C.18)
  Scalar a = sqrt_1_e2 * std::abs(z);
  Scalar b = WGS84_E2<Scalar> * WGS84_A<Scalar>;
  Scalar E = (a - b) / beta;             // (Groves C.29)
  Scalar F = (a + b) / beta;             // (Groves C.30)
  Scalar P = 4.0 / 3.0 * (E * F + 1.0);  // (Groves C.31)
  Scalar Q = 2.0 * (E * E - F * F);      // (Groves C.32)
  Scalar D = P * P * P + Q * Q;          // (Groves C.33)
  Scalar sqrt_D = std::sqrt(D);
  Scalar V = std::pow(sqrt_D - Q, 1.0 / 3.0) - std::pow(sqrt_D + Q, 1.0 / 3.0);  // (Groves C.34)
  Scalar G = 0.5 * (std::sqrt(E * E + V) + E);                                   // (Groves C.35)
  Scalar T = std::sqrt(G * G + ((F - V * G) / (2.0 * G - E))) - G;               // (Groves C.36)

  // (Groves C.37) &  (Groves C.38)
  auto &out = lla.derived();
  out(0) = sign_z * std::atan((1.0 - T * T) / (2.0 * T * sqrt_1_e2));
  out(1) = std::atan2(y, x);
  out(2) = (beta - WGS84_A<Scalar> * T) * std::cos(lla(0)) +
           (z - sign_z * WGS84_A<Scalar> * sqrt_1_e2) * std::sin(lla(0));
}

/**
 * *=== ecef2lla ===*
 * @brief converts Earth-centered-earth-fixed frame to Latitude-Longitude-Altitude position
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param xyz 3d vector of cartesian ECEF position (m,m,m)
 * @returns the 3d vector of ellipsoidal LLA position (rad,rad,m)
 */
template <typename Derived>
inline Eigen::Vector3<typename Derived::Scalar> ecef2lla(const Eigen::DenseBase<Derived> &xyz) {
  Eigen::Vector3<typename Derived::Scalar> lla;
  ecef2lla(xyz, lla);
  return lla;
}

/**
 * *=== lla2ecef ===*
 * @brief converts Latitude-Longitude-Altitude to Earth-centered-earth-fixed frame position
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param lla 3d vector of ellipsoidal LLA position (rad,rad,m)
 * @param xyz 3d vector of cartesian ECEF position (m,m,m)
 */
template <typename Derived1, typename Derived2>
inline void lla2ecef(const Eigen::DenseBase<Derived1> &lla, Eigen::DenseBase<Derived2> &xyz) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, lla, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, xyz, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  using Scalar = typename Derived1::Scalar;

  const auto &in = lla.derived();
  auto &out = xyz.derived();
  const Scalar slat = std::sin(in(0));
  const Scalar clat = std::cos(in(0));
  const Scalar slong = std::sin(in(1));
  const Scalar clong = std::cos(in(1));
  const Scalar &h = in(2);

  Scalar Re = WGS84_A<Scalar> / std::sqrt(1.0 - WGS84_E2<Scalar> * slat * slat);
  out << (h + Re) * clat * clong,                  // x
      (h + Re) * clat * slong,                     // y
      (h + Re * (1.0 - WGS84_E2<Scalar>)) * slat;  // z
}

/**
 * *=== lla2ecef ===*
 * @brief converts Latitude-Longitude-Altitude to Earth-centered-earth-fixed frame position
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param lla 3d vector of ellipsoidal LLA position (rad,rad,m)
 * @returns the 3d vector of cartesian ECEF position (m,m,m)
 */
template <typename Derived>
inline Eigen::Vector3<typename Derived::Scalar> lla2ecef(const Eigen::DenseBase<Derived> &lla) {
  Eigen::Vector3<typename Derived::Scalar> xyz;
  lla2ecef(lla, xyz);
  return xyz;
}

//! ============================================================================================ !//

/**
 * *=== ecef2ned ===*
 * @brief converts Earth-centered-earth-fixed frame to North-East-Down position
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param xyz  3d vector of cartesian ECEF position (m,m,m)
 * @param lla0 3d vector of reference ellipsoidal LLA position (deg,deg,m)
 * @param ned  3d vector of cartesian NED local position (m,m,m)
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void ecef2ned(
    const Eigen::DenseBase<Derived1> &xyz,
    const Eigen::DenseBase<Derived2> &lla0,
    Eigen::DenseBase<Derived3> &ned) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, xyz, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, ned, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, lla0, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  ASSERT_EIGEN_TYPE(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix3<Scalar> C_e_n;
  Eigen::Vector3<Scalar> xyz0;
  ecef2nedDcm(lla0, C_e_n);
  lla2ecef(lla0, xyz0);
  ned.derived() = C_e_n * (xyz.derived() - xyz0);
}

/**
 * *=== ecef2ned ===*
 * @brief converts Earth-centered-earth-fixed frame to North-East-Down position
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param xyz  3d vector of cartesian ECEF position (m,m,m)
 * @param lla0 3d vector of reference ellipsoidal LLA position (deg,deg,m)
 * @returns the 3d vector of cartesian NED local position (m,m,m)
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> ecef2ned(
    const Eigen::DenseBase<Derived1> &xyz, const Eigen::DenseBase<Derived2> &lla0) {
  Eigen::Vector3<typename Derived1::Scalar> ned;
  ecef2ned(xyz, lla0, ned);
  return ned;
}

/**
 * *=== ecef2enu ===*
 * @brief converts Earth-centered-earth-fixed frame to East-North-Up position
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param xyz  3d vector of cartesian ECEF position (m,m,m)
 * @param lla0 3d vector of reference ellipsoidal LLA position (deg,deg,m)
 * @param enu  3d vector of cartesian ENU local position (m,m,m)
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void ecef2enu(
    const Eigen::DenseBase<Derived1> &xyz,
    const Eigen::DenseBase<Derived2> &lla0,
    Eigen::DenseBase<Derived3> &enu) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, xyz, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, enu, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, lla0, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  ASSERT_EIGEN_TYPE(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix3<Scalar> C_e_n;
  Eigen::Vector3<Scalar> xyz0;
  ecef2enuDcm(lla0, C_e_n);
  lla2ecef(lla0, xyz0);
  enu.derived() = C_e_n * (xyz.derived() - xyz0);
}

/**
 * *=== ecef2enu ===*
 * @brief converts Earth-centered-earth-fixed frame to East-North-Up position
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param xyz  3d vector of cartesian ECEF position (m,m,m)
 * @param lla0 3d vector of reference ellipsoidal LLA position (deg,deg,m)
 * @returns the 3d vector of cartesian ENU local position (m,m,m)
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> ecef2enu(
    const Eigen::DenseBase<Derived1> &xyz, const Eigen::DenseBase<Derived2> &lla0) {
  Eigen::Vector3<typename Derived1::Scalar> enu;
  ecef2enu(xyz, lla0, enu);
  return enu;
}

/**
 * *=== ned2ecef ===*
 * @brief converts North-East-Down to Earth-Centered-Earth-Fixed position coordinates
 * @param ned  3x1 NED position [m]
 * @param lla0 3x1 Reference LLA position [rad, rad, m]
 * @param xyz  3x1 ECEF position [m]
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void ned2ecef(
    const Eigen::DenseBase<Derived1> &ned,
    const Eigen::DenseBase<Derived2> &lla0,
    Eigen::DenseBase<Derived3> &xyz) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, ned, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, xyz, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  ASSERT_EIGEN_TYPE(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix<Scalar, 3, 3> C_n_e;
  ned2ecefDcm(lla0, C_n_e);
  lla2ecef(lla0, xyz);
  xyz.derived() += C_n_e * ned.derived();
}

/**
 * *=== ned2ecef ===*
 * @brief converts North-East-Down to Earth-Centered-Earth-Fixed position coordinates
 * @param ned  3x1 NED position [m]
 * @param lla0 3x1 Reference LLA position [rad, rad, m]
 * @param xyz  3x1 ECEF position [m]
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> ned2ecef(
    const Eigen::DenseBase<Derived1> &ned, const Eigen::DenseBase<Derived2> &lla0) {
  Eigen::Vector3<typename Derived1::Scalar> xyz;
  ned2ecef(ned, lla0, xyz);
  return xyz;
}

/**
 * *=== enu2ecef ===*
 * @brief converts East-North-Up to Earth-Centered-Earth-Fixed position coordinates
 * @param enu  3x1 ENU position [m]
 * @param lla0 3x1 Reference LLA position [rad, rad, m]
 * @param xyz  3x1 ECEF position [m]
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void enu2ecef(
    const Eigen::DenseBase<Derived1> &enu,
    const Eigen::DenseBase<Derived2> &lla0,
    Eigen::DenseBase<Derived3> &xyz) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, enu, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, xyz, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  ASSERT_EIGEN_TYPE(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix<Scalar, 3, 3> C_n_e;
  enu2ecefDcm(lla0, C_n_e);
  lla2ecef(lla0, xyz);
  xyz.derived() += C_n_e * enu.derived();
}

/**
 * *=== enu2ecef ===*
 * @brief converts East-North-Up to Earth-Centered-Earth-Fixed position coordinates
 * @param enu  3x1 ENU position [m]
 * @param lla0 3x1 Reference LLA position [rad, rad, m]
 * @param xyz  3x1 ECEF position [m]
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> enu2ecef(
    const Eigen::DenseBase<Derived1> &enu, const Eigen::DenseBase<Derived2> &lla0) {
  Eigen::Vector3<typename Derived1::Scalar> xyz;
  enu2ecef(enu, lla0, xyz);
  return xyz;
}

//! ============================================================================================ !//

/**
 * *=== ecef2eci ===*
 * @brief converts Earth-Centered-Earth-Fixed to Earth-Centered-Inertial position coordinates
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param dt  time elapsed between frames (s)
 * @param xyz 3x1 ECEF position (m,m,m)
 * @param eci 3x1 ECI position (m,m,m)
 */
template <typename Derived1, typename Derived2>
inline void ecef2eci(
    const typename Derived1::Scalar dt,
    const Eigen::DenseBase<Derived1> &xyz,
    Eigen::DenseBase<Derived2> &eci) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, xyz, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, eci, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix<Scalar, 3, 3> C_e_i;
  ecef2eciDcm(dt, C_e_i);
  eci.derived() = C_e_i * xyz.derived();
}

/**
 * *=== ecef2eci ===*
 * @brief converts Earth-Centered-Earth-Fixed to Earth-Centered-Inertial position coordinates
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param dt  time elapsed between frames (s)
 * @param xyz 3x1 ECEF position (m,m,m)
 * @returns the 3x1 ECI position (m,m,m)
 */
template <typename Derived>
inline Eigen::Vector3<typename Derived::Scalar> ecef2eci(
    const typename Derived::Scalar dt, const Eigen::DenseBase<Derived> &xyz) {
  Eigen::Vector3<typename Derived::Scalar> eci;
  ecef2eci(dt, xyz, eci);
  return eci;
}

/**
 * *=== eci2ecef ===*
 * @brief converts Earth-Centered-Inertial to Earth-Centered-Earth-Fixed position coordinates
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param dt  time elapsed between frames (s)
 * @param eci 3x1 ECI position (m,m,m)
 * @param xyz 3x1 ECEF position (m,m,m)
 */
template <typename Derived1, typename Derived2>
inline void eci2ecef(
    const typename Derived1::Scalar dt,
    const Eigen::DenseBase<Derived1> &eci,
    Eigen::DenseBase<Derived2> &xyz) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, eci, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, xyz, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix<Scalar, 3, 3> C_i_e;
  eci2ecefDcm(dt, C_i_e);
  xyz.derived() = C_i_e * eci.derived();
}

/**
 * *=== eci2ecef ===*
 * @brief converts Earth-Centered-Inertial to Earth-Centered-Earth-Fixed position coordinates
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param dt  time elapsed between frames (s)
 * @param eci 3x1 ECI position (m,m,m)
 * @returns the 3x1 ECEF position (m,m,m)
 */
template <typename Derived>
inline Eigen::Vector3<typename Derived::Scalar> eci2ecef(
    const typename Derived::Scalar dt, const Eigen::DenseBase<Derived> &eci) {
  Eigen::Vector3<typename Derived::Scalar> xyz;
  eci2ecef(dt, eci, xyz);
  return xyz;
}

//! ============================================================================================ !//

/**
 * *=== lla2ned ===*
 * @brief converts Latitude-Longitude-Altitude to North-East-Down frame position
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param lla  3d vector of ellipsoidal LLA position (rad,rad,m)
 * @param lla0 3d vector of reference ellipsoidal LLA position (rad,rad,m)
 * @param ned  3d vector of cartesian NED position (m,m,m)
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void lla2ned(
    const Eigen::DenseBase<Derived1> &lla,
    const Eigen::DenseBase<Derived2> &lla0,
    Eigen::DenseBase<Derived3> &ned) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, ned, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  ASSERT_EIGEN_TYPE(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Vector3<Scalar> xyz;
  lla2ecef(lla, xyz);
  ecef2ned(xyz, lla0, ned);
}

/**
 * *=== lla2ned ===*
 * @brief converts Latitude-Longitude-Altitude to North-East-Down frame position
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param lla  3d vector of ellipsoidal LLA position (rad,rad,m)
 * @param lla0 3d vector of reference ellipsoidal LLA position (rad,rad,m)
 * @returns the 3d vector of cartesian NED position (m,m,m)
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> lla2ned(
    const Eigen::DenseBase<Derived1> &lla, const Eigen::DenseBase<Derived2> &lla0) {
  Eigen::Vector3<typename Derived1::Scalar> ned;
  lla2ned(lla, lla0, ned);
  return ned;
}

/**
 * *=== lla2enu ===*
 * @brief converts Latitude-Longitude-Altitude to East-North-Up frame position
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param lla  3d vector of ellipsoidal LLA position (rad,rad,m)
 * @param lla0 3d vector of reference ellipsoidal LLA position (rad,rad,m)
 * @param enu  3d vector of cartesian ENU position (m,m,m)
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void lla2enu(
    const Eigen::DenseBase<Derived1> &lla,
    const Eigen::DenseBase<Derived2> &lla0,
    Eigen::DenseBase<Derived3> &enu) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, enu, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  ASSERT_EIGEN_TYPE(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Vector3<Scalar> xyz;
  lla2ecef(lla, xyz);
  ecef2enu(xyz, lla0, enu);
}

/**
 * *=== lla2enu ===*
 * @brief converts Latitude-Longitude-Altitude to North-East-Down frame position
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param lla  3d vector of ellipsoidal LLA position (rad,rad,m)
 * @param lla0 3d vector of reference ellipsoidal LLA position (rad,rad,m)
 * @returns the 3d vector of cartesian ENU position (m,m,m)
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> lla2enu(
    const Eigen::DenseBase<Derived1> &lla, const Eigen::DenseBase<Derived2> &lla0) {
  Eigen::Vector3<typename Derived1::Scalar> enu;
  lla2enu(lla, lla0, enu);
  return enu;
}

/**
 * *=== ned2lla ===*
 * @brief converts North-East-Down to Latitude-Longitude-Height position coordinates
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param ned  3x1 NED position [m]
 * @param lla0 3x1 Reference LLA position [rad, rad, m]
 * @param lla  3x1 LLA position [m]
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void ned2lla(
    const Eigen::DenseBase<Derived1> &ned,
    const Eigen::DenseBase<Derived2> &lla0,
    Eigen::DenseBase<Derived3> &lla) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, ned, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, lla, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  ASSERT_EIGEN_TYPE(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Vector3<Scalar> xyz;
  ned2ecef(ned, lla0, xyz);
  ecef2lla(xyz, lla);
}

/**
 * *=== ned2lla ===*
 * @brief converts North-East-Down to Latitude-Longitude-Height position coordinates
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param ned  3x1 NED position [m]
 * @param lla0 3x1 Reference LLA position [rad, rad, m]
 * @returns the 3x1 LLA position [m]
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> ned2lla(
    const Eigen::DenseBase<Derived1> &ned, const Eigen::DenseBase<Derived2> &lla0) {
  Eigen::Vector3<typename Derived1::Scalar> lla;
  ned2lla(ned, lla0, lla);
  return lla;
}

/**
 * *=== enu2lla ===*
 * @brief converts East-North-Up to Latitude-Longitude-Height position coordinates
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param enu  3x1 ENU position [m]
 * @param lla0 3x1 Reference LLA position [rad, rad, m]
 * @param lla  3x1 LLA position [m]
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void enu2lla(
    const Eigen::DenseBase<Derived1> &enu,
    const Eigen::DenseBase<Derived2> &lla0,
    Eigen::DenseBase<Derived3> &lla) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, enu, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, lla, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  ASSERT_EIGEN_TYPE(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Vector3<Scalar> xyz;
  enu2ecef(enu, lla0, xyz);
  ecef2lla(xyz, lla);
}

/**
 * *=== enu2lla ===*
 * @brief converts East-North-Up to Latitude-Longitude-Height position coordinates
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param ned  3x1 ENU position [m]
 * @param lla0 3x1 Reference LLA position [rad, rad, m]
 * @returns the 3x1 LLA position [m]
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> enu2lla(
    const Eigen::DenseBase<Derived1> &enu, const Eigen::DenseBase<Derived2> &lla0) {
  Eigen::Vector3<typename Derived1::Scalar> lla;
  enu2lla(enu, lla0, lla);
  return lla;
}

//! ============================================================================================ !//

/**
 * *=== lla2eci ===*
 * @brief converts Latitude-Longitude-Height to Earth-Centered-Inertial position coordinates
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param dt  time elapsed between frames (s)
 * @param lla 3x1 Geodetic Latitude, Longitude, Height (rad, rad, m)
 * @param eci 3x1 ECI position (m)
 */
template <typename Derived1, typename Derived2>
inline void lla2eci(
    const typename Derived1::Scalar dt,
    const Eigen::DenseBase<Derived1> &lla,
    Eigen::DenseBase<Derived2> &eci) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, lla, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, eci, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  using Scalar = typename Derived1::Scalar;

  Eigen::Vector3<Scalar> xyz;
  lla2ecef(lla, xyz);
  Eigen::Matrix<Scalar, 3, 3> C_e_i;
  ecef2eciDcm(dt, C_e_i);
  eci.derived() = C_e_i * xyz;
}

/**
 * *=== lla2eci ===*
 * @brief converts Latitude-Longitude-Height to Earth-Centered-Inertial position coordinates
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param dt  time elapsed between frames (s)
 * @param lla 3x1 Geodetic Latitude, Longitude, Height (rad, rad, m)
 * @returns the 3x1 ECI position (m)
 */
template <typename Derived>
inline Eigen::Vector3<typename Derived::Scalar> lla2eci(
    const typename Derived::Scalar dt, const Eigen::DenseBase<Derived> &lla) {
  Eigen::Vector3<typename Derived::Scalar> eci;
  lla2eci(dt, lla, eci);
  return eci;
}

/**
 * *=== eci2lla ===*
 * @brief converts Earth-Centered-Inertial to Latitude-Longitude-Height position coordinates
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param dt   time elapsed between frames (s)
 * @param eci  3x1 ECI position (m,m,m)
 * @param lla  3x1 LLA position (rad,rad,m)
 */
template <typename Derived1, typename Derived2>
inline void eci2lla(
    const typename Derived1::Scalar dt,
    const Eigen::DenseBase<Derived1> &eci,
    Eigen::DenseBase<Derived2> &lla) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, eci, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived1, lla, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  using Scalar = typename Derived1::Scalar;

  Eigen::Vector3<Scalar> xyz;
  eci2ecef(dt, eci, xyz);
  ecef2lla(xyz, lla);
}

/**
 * *=== eci2lla ===*
 * @brief converts Earth-Centered-Inertial to Latitude-Longitude-Height position coordinates
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param dt  time elapsed between frames (s)
 * @param eci 3x1 ECI position (m,m,m)
 * @returns the 3x1 LLA position (rad,rad,m)
 */
template <typename Derived>
inline Eigen::Vector3<typename Derived::Scalar> eci2lla(
    const typename Derived::Scalar dt, const Eigen::DenseBase<Derived> &eci) {
  Eigen::Vector3<typename Derived::Scalar> lla;
  eci2lla(dt, eci, lla);
  return lla;
}

/**
 * *=== eci2ned ===*
 * @brief converts Earth-Centered-Inertial to North-East-Down position coordinates
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param dt   time elapsed between frames (s)
 * @param eci  3x1 ECI position (m,m,m)
 * @param lla0 3x1 Reference Geodetic Latitude, Longitude, Height (rad,rad,m)
 * @param ned  3x1 NED position (m,m,m)
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void eci2ned(
    const typename Derived1::Scalar dt,
    const Eigen::DenseBase<Derived1> &eci,
    const Eigen::DenseBase<Derived2> &lla0,
    Eigen::DenseBase<Derived3> &ned) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, eci, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, ned, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  ASSERT_EIGEN_TYPE(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Vector3<Scalar> xyz;
  eci2ecef(dt, eci, xyz);
  ecef2ned(xyz, lla0, ned);
}

/**
 * *=== eci2ned ===*
 * @brief converts Earth-Centered-Inertial to North-East-Down position coordinates
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param dt   time elapsed between frames (s)
 * @param eci  3x1 ECI position (m,m,m)
 * @param lla0 3x1 Reference Geodetic Latitude, Longitude, Height (rad,rad,m)
 * @returns the 3x1 NED position (m,m,m)
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> eci2ned(
    const typename Derived1::Scalar dt,
    const Eigen::DenseBase<Derived1> &eci,
    const Eigen::DenseBase<Derived2> &lla0) {
  Eigen::Vector3<typename Derived1::Scalar> ned;
  eci2ned(dt, eci, lla0, ned);
  return ned;
}

/**
 * *=== ned2eci ===*
 * @brief converts North-East-Down to Earth-Centered-Inertial position coordinates
 * @param dt   time elapsed between frames [s]
 * @param ned  3x1 NED position [m]
 * @param lla0 3x1 Reference LLA position [rad, rad, m]
 * @param eci  3x1 ECI position [m]
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void ned2eci(
    const typename Derived1::Scalar dt,
    const Eigen::DenseBase<Derived1> &ned,
    const Eigen::DenseBase<Derived2> &lla0,
    Eigen::DenseBase<Derived3> &eci) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, ned, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, eci, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  ASSERT_EIGEN_TYPE(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Vector3<Scalar> xyz;
  ned2ecef(ned, lla0, xyz);
  Eigen::Matrix<Scalar, 3, 3> C_e_i;
  ecef2eciDcm(dt, C_e_i);
  eci.derived() = C_e_i * xyz;
}

/**
 * *=== ned2eci ===*
 * @brief converts North-East-Down to Earth-Centered-Inertial position coordinates
 * @param dt   time elapsed between frames [s]
 * @param ned  3x1 NED position [m]
 * @param lla0 3x1 Reference LLA position [rad, rad, m]
 * @returns the 3x1 ECI position [m]
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> ned2eci(
    const typename Derived1::Scalar dt,
    const Eigen::DenseBase<Derived1> &ned,
    const Eigen::DenseBase<Derived2> &lla0) {
  Eigen::Vector3<typename Derived1::Scalar> eci;
  ned2eci(dt, ned, lla0, eci);
  return eci;
}

/**
 * *=== eci2enu ===*
 * @brief converts Earth-Centered-Inertial to East-North-Up position coordinates
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param dt   time elapsed between frames (s)
 * @param eci  3x1 ECI position (m,m,m)
 * @param lla0 3x1 Reference Geodetic Latitude, Longitude, Height (rad,rad,m)
 * @param enu  3x1 ENU position (m,m,m)
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void eci2enu(
    const typename Derived1::Scalar dt,
    const Eigen::DenseBase<Derived1> &eci,
    const Eigen::DenseBase<Derived2> &lla0,
    Eigen::DenseBase<Derived3> &enu) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, eci, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, enu, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  ASSERT_EIGEN_TYPE(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Vector3<Scalar> xyz;
  eci2ecef(dt, eci, xyz);
  ecef2enu(xyz, lla0, enu);
}

/**
 * *=== eci2enu ===*
 * @brief converts Earth-Centered-Inertial to East-North-Up position coordinates
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param dt   time elapsed between frames (s)
 * @param eci  3x1 ECI position (m,m,m)
 * @param lla0 3x1 Reference Geodetic Latitude, Longitude, Height (rad,rad,m)
 * @returns the 3x1 NED position (m,m,m)
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> eci2enu(
    const typename Derived1::Scalar dt,
    const Eigen::DenseBase<Derived1> &eci,
    const Eigen::DenseBase<Derived2> &lla0) {
  Eigen::Vector3<typename Derived1::Scalar> enu;
  eci2enu(dt, eci, lla0, enu);
  return enu;
}

/**
 * *=== enu2eci ===*
 * @brief converts East-North-Up to Earth-Centered-Inertial position coordinates
 * @param dt   time elapsed between frames [s]
 * @param enu  3x1 ENU position [m]
 * @param lla0 3x1 Reference LLA position [rad, rad, m]
 * @param eci  3x1 ECI position [m]
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void enu2eci(
    const typename Derived1::Scalar dt,
    const Eigen::DenseBase<Derived1> &enu,
    const Eigen::DenseBase<Derived2> &lla0,
    Eigen::DenseBase<Derived3> &eci) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, enu, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, eci, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  ASSERT_EIGEN_TYPE(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Vector3<Scalar> xyz;
  enu2ecef(enu, lla0, xyz);
  Eigen::Matrix<Scalar, 3, 3> C_e_i;
  ecef2eciDcm(dt, C_e_i);
  eci.derived() = C_e_i * xyz;
}

/**
 * *=== enu2eci ===*
 * @brief converts East-North-Up to Earth-Centered-Inertial position coordinates
 * @param dt   time elapsed between frames [s]
 * @param enu  3x1 ENU position [m]
 * @param lla0 3x1 Reference LLA position [rad, rad, m]
 * @returns the 3x1 ECI position [m]
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> enu2eci(
    const typename Derived1::Scalar dt,
    const Eigen::DenseBase<Derived1> &enu,
    const Eigen::DenseBase<Derived2> &lla0) {
  Eigen::Vector3<typename Derived1::Scalar> eci;
  enu2eci(dt, enu, lla0, eci);
  return eci;
}

//! ============================================================================================ !//

/**
 * *=== ned2enu ===*
 * @brief converts North-East-Down to East-North-Up position coordinates
 * @param ned  3x1 NED position [m]
 * @param enu  3x1 ENU position [m]
 */
template <typename Derived1, typename Derived2>
inline void ned2enu(const Eigen::DenseBase<Derived1> &ned, Eigen::DenseBase<Derived2> &enu) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, ned, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, enu, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix<Scalar, 3, 3> R;
  ned2enuDcm(R);
  enu.derived() = R * ned.derived();
}

/**
 * *=== ned2enu ===*
 * @brief converts North-East-Down to East-North-Up position coordinates
 * @param ned  3x1 NED position [m]
 * @param enu  3x1 ENU position [m]
 */
template <typename Derived>
inline Eigen::Vector3<typename Derived::Scalar> ned2enu(const Eigen::DenseBase<Derived> &ned) {
  Eigen::Vector3<typename Derived::Scalar> enu;
  ned2enu(ned, enu);
  return enu;
}

/**
 * *=== enu2ned ===*
 * @brief converts East-North-Up to North-East-Dnw  position coordinates
 * @param enu  3x1 ENU position [m]
 * @param ned  3x1 ENU position [m]
 */
template <typename Derived1, typename Derived2>
inline void enu2ned(const Eigen::DenseBase<Derived1> &enu, Eigen::DenseBase<Derived2> &ned) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, ned, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, enu, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix<Scalar, 3, 3> R;
  enu2nedDcm(R);
  ned.derived() = R * enu.derived();
}

/**
 * *=== enu2ned ===*
 * @brief converts East-North-Up to North-East-Down position coordinates
 * @param enu  3x1 ENU position [m]
 * @param ned  3x1 NED position [m]
 */
template <typename Derived>
inline Eigen::Vector3<typename Derived::Scalar> enu2ned(const Eigen::DenseBase<Derived> &enu) {
  Eigen::Vector3<typename Derived::Scalar> ned;
  enu2ned(enu, ned);
  return ned;
}

}  // namespace nt

#endif