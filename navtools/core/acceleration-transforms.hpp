#ifndef NAVTOOLS_CORE_ACCELERATION_TRANSFORMS
#define NAVTOOLS_CORE_ACCELERATION_TRANSFORMS

#include "navtools/core/constants.hpp"
#include "navtools/core/coordinate-dcms.hpp"
#include "navtools/core/macros.hpp"
#include "navtools/core/position-transforms.hpp"
#include <Eigen/Dense>

namespace nt {

/**
 * *=== ecef2neda ===*
 * @brief Converts ECEF acceleration to NED acceleration
 * @param a_xyz Acceleration vector in the ECEF frame (m/s^2)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @param a_ned The acceleration vector in the NED frame (m/s^2)
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void ecef2neda(
    const Eigen::DenseBase<Derived1>& a_xyz,
    const Eigen::DenseBase<Derived2>& lla0,
    Eigen::DenseBase<Derived3>& a_ned) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, a_xyz, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, a_ned, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  ASSERT_EIGEN_TYPE(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix3<Scalar> C_e_n;
  ecef2nedDcm(lla0, C_e_n);
  a_ned.derived() = C_e_n * a_xyz.derived();
}

/**
 * *=== ecef2neda ===*
 * @brief Converts ECEF acceleration to NED acceleration
 * @param a_xyz Acceleration vector in the ECEF frame (m/s^2)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @returns The acceleration vector in the NED frame (m/s^2)
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> ecef2neda(
    const Eigen::DenseBase<Derived1>& a_xyz, const Eigen::DenseBase<Derived2>& lla0) {
  Eigen::Vector3<typename Derived1::Scalar> a_ned;
  ecef2neda(a_xyz, lla0, a_ned);
  return a_ned;
}

/**
 * *=== ecef2enua ===*
 * @brief Converts ECEF acceleration to ENU acceleration
 * @param a_xyz Acceleration vector in the ECEF frame (m/s^2)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @param a_enu The acceleration vector in the ENU frame (m/s^2)
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void ecef2enua(
    const Eigen::DenseBase<Derived1>& a_xyz,
    const Eigen::DenseBase<Derived2>& lla0,
    Eigen::DenseBase<Derived3>& a_enu) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, a_xyz, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, a_enu, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  ASSERT_EIGEN_TYPE(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix3<Scalar> C_e_n;
  ecef2enuDcm(lla0, C_e_n);
  a_enu.derived() = C_e_n * a_xyz.derived();
}

/**
 * *=== ecef2enua ===*
 * @brief Converts ECEF acceleration to ENU acceleration
 * @param a_xyz Acceleration vector in the ECEF frame (m/s^2)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @returns The acceleration vector in the ENU frame (m/s^2)
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> ecef2enua(
    const Eigen::DenseBase<Derived1>& a_xyz, const Eigen::DenseBase<Derived2>& lla0) {
  Eigen::Vector3<typename Derived1::Scalar> a_enu;
  ecef2enua(a_xyz, lla0, a_enu);
  return a_enu;
}

/**
 * *=== ned2ecefa ===*
 * @brief Converts NED acceleration to ECEF acceleration
 * @param a_ned Acceleration vector in the NED frame (m/s^2)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @param a_xyz The acceleration vector in the ECEF frame (m/s^2)
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void ned2ecefa(
    const Eigen::DenseBase<Derived1>& a_ned,
    const Eigen::DenseBase<Derived2>& lla0,
    Eigen::DenseBase<Derived3>& a_xyz) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, a_ned, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, a_xyz, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  ASSERT_EIGEN_TYPE(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix3<Scalar> C_n_e;
  ned2ecefDcm(lla0, C_n_e);
  a_xyz.derived() = C_n_e * a_ned.derived();
}

/**
 * *=== ned2ecefa ===*
 * @brief Converts NED acceleration to ECEF acceleration
 * @param a_ned Acceleration vector in the NED frame (m/s^2)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @returns The acceleration vector in the ECEF frame (m/s^2)
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> ned2ecefa(
    const Eigen::DenseBase<Derived1>& a_ned, const Eigen::DenseBase<Derived2>& lla0) {
  Eigen::Vector3<typename Derived1::Scalar> a_xyz;
  ned2ecefa(a_ned, lla0, a_xyz);
  return a_xyz;
}

/**
 * *=== enu2ecefa ===*
 * @brief Converts ENU acceleration to ECEF acceleration
 * @param a_enu Acceleration vector in the ENU frame (m/s^2)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @param a_xyz The acceleration vector in the ECEF frame (m/s^2)
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void enu2ecefa(
    const Eigen::DenseBase<Derived1>& a_enu,
    const Eigen::DenseBase<Derived2>& lla0,
    Eigen::DenseBase<Derived3>& a_xyz) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, a_enu, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, a_xyz, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  ASSERT_EIGEN_TYPE(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix3<Scalar> C_n_e;
  enu2ecefDcm(lla0, C_n_e);
  a_xyz.derived() = C_n_e * a_enu.derived();
}

/**
 * *=== enu2ecefa ===*
 * @brief Converts ENU acceleration to ECEF acceleration
 * @param a_enu Acceleration vector in the ENU frame (m/s^2)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @returns The acceleration vector in the ECEF frame (m/s^2)
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> enu2ecefa(
    const Eigen::DenseBase<Derived1>& a_enu, const Eigen::DenseBase<Derived2>& lla0) {
  Eigen::Vector3<typename Derived1::Scalar> a_xyz;
  enu2ecefa(a_enu, lla0, a_xyz);
  return a_xyz;
}

/**
 * *=== ecef2ecia ===*
 * @brief Converts ECEF acceleration to ECI acceleration
 * @param dt    Elapsed time (s)
 * @param a_xyz Acceleration vector in the ECEF frame (m/s^2)
 * @param v_xyz Velocity vector in the ECEF frame (m/s)
 * @param r_xyz Position vector in the ECEF frame (m)
 * @param a_eci The acceleration vector in the ECI frame (m/s^2)
 */
template <typename Derived1, typename Derived2, typename Derived3, typename Derived4>
inline void ecef2ecia(
    const typename Derived1::Scalar dt,
    const Eigen::DenseBase<Derived1>& a_xyz,
    const Eigen::DenseBase<Derived2>& v_xyz,
    const Eigen::DenseBase<Derived3>& r_xyz,
    Eigen::DenseBase<Derived4>& a_eci) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, a_xyz, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, v_xyz, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, r_xyz, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived4, a_eci, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  ASSERT_EIGEN_TYPE(Derived1, Derived3);
  ASSERT_EIGEN_TYPE(Derived1, Derived4);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix3<Scalar> C_e_i;
  ecef2eciDcm(dt, C_e_i);
  Eigen::Vector3<Scalar> a_coriolis = 2.0 * OMEGA_ECEF<Scalar>.cross(v_xyz.derived());
  Eigen::Vector3<Scalar> a_centripetal =
      OMEGA_ECEF<Scalar>.cross(OMEGA_ECEF<Scalar>.cross(r_xyz.derived()));
  a_eci.derived() = C_e_i * (a_xyz.derived() + a_coriolis + a_centripetal);
}

/**
 * *=== ecef2ecia ===*
 * @brief Converts ECEF acceleration to ECI acceleration
 * @param dt    Elapsed time (s)
 * @param a_xyz Acceleration vector in the ECEF frame (m/s^2)
 * @param v_xyz Velocity vector in the ECEF frame (m/s)
 * @param r_xyz Position vector in the ECEF frame (m)
 * @returns The acceleration vector in the ECI frame (m/s^2)
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline Eigen::Vector3<typename Derived1::Scalar> ecef2ecia(
    const typename Derived1::Scalar dt,
    const Eigen::DenseBase<Derived1>& a_xyz,
    const Eigen::DenseBase<Derived2>& v_xyz,
    const Eigen::DenseBase<Derived3>& r_xyz) {
  Eigen::Vector3<typename Derived1::Scalar> a_eci;
  ecef2ecia(dt, a_xyz, v_xyz, r_xyz, a_eci);
  return a_eci;
}

/**
 * *=== eci2ecefa ===*
 * @brief Converts ECI acceleration to ECEF acceleration
 * @param dt    Time elapsed (s)
 * @param a_eci Acceleration vector in the ECI frame (m/s^2)
 * @param v_eci Velocity vector in the ECI frame (m/s)
 * @param r_eci Position vector in the ECI frame (m)
 * @param a_xyz The acceleration vector in the ECEF frame (m/s^2)
 */
template <typename Derived1, typename Derived2, typename Derived3, typename Derived4>
inline void eci2ecefa(
    const typename Derived1::Scalar dt,
    const Eigen::DenseBase<Derived1>& a_eci,
    const Eigen::DenseBase<Derived2>& v_eci,
    const Eigen::DenseBase<Derived3>& r_eci,
    Eigen::DenseBase<Derived4>& a_xyz) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, a_eci, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, v_eci, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, r_eci, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived4, a_xyz, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  ASSERT_EIGEN_TYPE(Derived1, Derived3);
  ASSERT_EIGEN_TYPE(Derived1, Derived4);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix3<Scalar> C_i_e;
  eci2ecefDcm(dt, C_i_e);
  Eigen::Vector3<Scalar> a_coriolis = 2.0 * OMEGA_ECEF<Scalar>.cross(v_eci.derived());
  Eigen::Vector3<Scalar> a_centripetal =
      OMEGA_ECEF<Scalar>.cross(OMEGA_ECEF<Scalar>.cross(r_eci.derived()));
  a_xyz.derived() = C_i_e * (a_eci.derived() - a_coriolis + a_centripetal);
}

/**
 * *=== eci2ecefa ===*
 * @brief Converts ECI acceleration to ECEF acceleration
 * @param dt    Time elapsed (s)
 * @param a_eci Acceleration vector in the ECI frame (m/s^2)
 * @param v_eci Velocity vector in the ECI frame (m/s)
 * @param r_eci Position vector in the ECI frame (m)
 * @returns The acceleration vector in the ECEF frame (m/s^2)
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline Eigen::Vector3<typename Derived1::Scalar> eci2ecefa(
    const typename Derived1::Scalar dt,
    const Eigen::DenseBase<Derived1>& a_eci,
    const Eigen::DenseBase<Derived2>& v_eci,
    const Eigen::DenseBase<Derived3>& r_eci) {
  Eigen::Vector3<typename Derived1::Scalar> a_xyz;
  eci2ecefa(dt, a_eci, v_eci, r_eci, a_xyz);
  return a_xyz;
}

/**
 * *=== eci2neda ===*
 * @brief Converts ECI acceleration to NED acceleration
 * @param dt    Time elapsed (s)
 * @param a_eci Acceleration vector in the ECI frame (m/s^2)
 * @param v_eci Velocity vector in the ECI frame (m/s)
 * @param r_eci Position vector in the ECI frame (m)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @param a_ned The acceleration vector in the NED frame (m/s^2)
 */
template <
    typename Derived1,
    typename Derived2,
    typename Derived3,
    typename Derived4,
    typename Derived5>
inline void eci2neda(
    const typename Derived1::Scalar dt,
    const Eigen::DenseBase<Derived1>& a_eci,
    const Eigen::DenseBase<Derived2>& v_eci,
    const Eigen::DenseBase<Derived3>& r_eci,
    const Eigen::DenseBase<Derived4>& lla0,
    Eigen::DenseBase<Derived5>& a_ned) {
  using Scalar = typename Derived1::Scalar;

  Eigen::Vector3<Scalar> a_xyz;
  eci2ecefa(dt, a_eci, v_eci, r_eci, a_xyz);
  ecef2neda(a_xyz, lla0, a_ned);
}

/**
 * *=== eci2neda ===*
 * @brief Converts ECI acceleration to NED acceleration
 * @param dt    Time elapsed (s)
 * @param a_eci Acceleration vector in the ECI frame (m/s^2)
 * @param v_eci Velocity vector in the ECI frame (m/s)
 * @param r_eci Position vector in the ECI frame (m)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @returns The acceleration vector in the NED frame (m/s^2)
 */
template <typename Derived1, typename Derived2, typename Derived3, typename Derived4>
inline Eigen::Vector3<typename Derived1::Scalar> eci2neda(
    const typename Derived1::Scalar dt,
    const Eigen::DenseBase<Derived1>& a_eci,
    const Eigen::DenseBase<Derived2>& v_eci,
    const Eigen::DenseBase<Derived3>& r_eci,
    const Eigen::DenseBase<Derived4>& lla0) {
  Eigen::Vector3<typename Derived1::Scalar> a_ned;
  eci2neda(dt, a_eci, v_eci, r_eci, lla0, a_ned);
  return a_ned;
}

/**
 * *=== eci2enua ===*
 * @brief Converts ECI acceleration to ENU acceleration
 * @param dt    Time elapsed (s)
 * @param a_eci Acceleration vector in the ECI frame (m/s^2)
 * @param v_eci Velocity vector in the ECI frame (m/s)
 * @param r_eci Position vector in the ECI frame (m)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @param a_enu The acceleration vector in the ENU frame (m/s^2)
 */
template <
    typename Derived1,
    typename Derived2,
    typename Derived3,
    typename Derived4,
    typename Derived5>
inline void eci2enua(
    const typename Derived1::Scalar dt,
    const Eigen::DenseBase<Derived1>& a_eci,
    const Eigen::DenseBase<Derived2>& v_eci,
    const Eigen::DenseBase<Derived3>& r_eci,
    const Eigen::DenseBase<Derived4>& lla0,
    Eigen::DenseBase<Derived5>& a_enu) {
  using Scalar = typename Derived1::Scalar;

  Eigen::Vector3<Scalar> a_xyz;
  eci2ecefa(dt, a_eci, v_eci, r_eci, a_xyz);
  ecef2enua(a_xyz, lla0, a_enu);
}

/**
 * *=== eci2enua ===*
 * @brief Converts ECI acceleration to ENU acceleration
 * @param dt    Time elapsed (s)
 * @param a_eci Acceleration vector in the ECI frame (m/s^2)
 * @param v_eci Velocity vector in the ECI frame (m/s)
 * @param r_eci Position vector in the ECI frame (m)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @returns The acceleration vector in the ENU frame (m/s^2)
 */
template <typename Derived1, typename Derived2, typename Derived3, typename Derived4>
inline Eigen::Vector3<typename Derived1::Scalar> eci2enua(
    const typename Derived1::Scalar dt,
    const Eigen::DenseBase<Derived1>& a_eci,
    const Eigen::DenseBase<Derived2>& v_eci,
    const Eigen::DenseBase<Derived3>& r_eci,
    const Eigen::DenseBase<Derived4>& lla0) {
  Eigen::Vector3<typename Derived1::Scalar> a_enu;
  eci2enua(dt, a_eci, v_eci, r_eci, lla0, a_enu);
  return a_enu;
}

/**
 * *=== ned2ecia ===*
 * @brief Converts NED acceleration to ECI acceleration.
 * @param dt    Time elapsed (s)
 * @param a_ned Acceleration vector in the NED frame (m/s^2)
 * @param v_ned Velocity vector in the NED frame (m/s)
 * @param r_ned Position vector in the NED frame (m)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @param a_eci The acceleration vector in the ECI frame (m/s^2)
 */
template <
    typename Derived1,
    typename Derived2,
    typename Derived3,
    typename Derived4,
    typename Derived5>
inline void ned2ecia(
    const typename Derived1::Scalar dt,
    const Eigen::DenseBase<Derived1>& a_ned,
    const Eigen::DenseBase<Derived2>& v_ned,
    const Eigen::DenseBase<Derived3>& r_ned,
    const Eigen::DenseBase<Derived4>& lla0,
    Eigen::DenseBase<Derived5>& a_eci) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, a_ned, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, v_ned, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, r_ned, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived4, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived5, a_eci, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  ASSERT_EIGEN_TYPE(Derived1, Derived3);
  ASSERT_EIGEN_TYPE(Derived1, Derived4);
  ASSERT_EIGEN_TYPE(Derived1, Derived5);
  using Scalar = typename Derived1::Scalar;

  Eigen::Vector3<Scalar> a_xyz, v_xyz, r_xyz, r_xyz0;
  Eigen::Matrix3<Scalar> C_n_e;

  ned2ecefDcm(lla0, C_n_e);
  lla2ecef(lla0, r_xyz0);
  r_xyz = r_xyz0 + C_n_e * r_ned.derived();
  v_xyz = C_n_e * v_ned.derived();
  a_xyz = C_n_e * a_ned.derived();
  ecef2ecia(dt, a_xyz, v_xyz, r_xyz, a_eci);
}

/**
 * *=== ned2ecia ===*
 * @brief Converts NED acceleration to ECI acceleration.
 * @param dt    Time elapsed (s)
 * @param a_ned Acceleration vector in the NED frame (m/s^2)
 * @param v_ned Velocity vector in the NED frame (m/s)
 * @param r_ned Position vector in the NED frame (m)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @returns The acceleration vector in the ECI frame (m/s^2)
 */
template <typename Derived1, typename Derived2, typename Derived3, typename Derived4>
inline Eigen::Vector3<typename Derived1::Scalar> ned2ecia(
    const typename Derived1::Scalar dt,
    const Eigen::DenseBase<Derived1>& a_ned,
    const Eigen::DenseBase<Derived2>& v_ned,
    const Eigen::DenseBase<Derived3>& r_ned,
    const Eigen::DenseBase<Derived4>& lla0) {
  Eigen::Vector3<typename Derived1::Scalar> a_eci;
  ned2ecia(dt, a_ned, v_ned, r_ned, lla0, a_eci);
  return a_eci;
}

/**
 * *=== enu2ecia ===*
 * @brief Converts ENU acceleration to ECI acceleration.
 * @param dt    Time elapsed (s)
 * @param a_enu Acceleration vector in the ENU frame (m/s^2)
 * @param v_enu Velocity vector in the ENU frame (m/s)
 * @param r_enu Position vector in the ENU frame (m)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @param a_eci The acceleration vector in the ECI frame (m/s^2)
 */
template <
    typename Derived1,
    typename Derived2,
    typename Derived3,
    typename Derived4,
    typename Derived5>
inline void enu2ecia(
    const typename Derived1::Scalar dt,
    const Eigen::DenseBase<Derived1>& a_enu,
    const Eigen::DenseBase<Derived2>& v_enu,
    const Eigen::DenseBase<Derived3>& r_enu,
    const Eigen::DenseBase<Derived4>& lla0,
    Eigen::DenseBase<Derived5>& a_eci) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, a_enu, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, v_enu, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, r_enu, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived4, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived5, a_eci, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  ASSERT_EIGEN_TYPE(Derived1, Derived3);
  ASSERT_EIGEN_TYPE(Derived1, Derived4);
  ASSERT_EIGEN_TYPE(Derived1, Derived5);
  using Scalar = typename Derived1::Scalar;

  Eigen::Vector3<Scalar> a_xyz, v_xyz, r_xyz, r_xyz0;
  Eigen::Matrix3<Scalar> C_n_e;

  enu2ecefDcm(lla0, C_n_e);
  lla2ecef(lla0, r_xyz0);
  r_xyz = r_xyz0 + C_n_e * r_enu.derived();
  v_xyz = C_n_e * v_enu.derived();
  a_xyz = C_n_e * a_enu.derived();
  ecef2ecia(dt, a_xyz, v_xyz, r_xyz, a_eci);
}

/**
 * *=== enu2ecia ===*
 * @brief Converts ENU acceleration to ECI acceleration.
 * @param dt    Time elapsed (s)
 * @param a_enu Acceleration vector in the ENU frame (m/s^2)
 * @param v_enu Velocity vector in the ENU frame (m/s)
 * @param r_enu Position vector in the ENU frame (m)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @returns The acceleration vector in the ECI frame (m/s^2)
 */
template <typename Derived1, typename Derived2, typename Derived3, typename Derived4>
inline Eigen::Vector3<typename Derived1::Scalar> enu2ecia(
    const typename Derived1::Scalar dt,
    const Eigen::DenseBase<Derived1>& a_enu,
    const Eigen::DenseBase<Derived2>& v_enu,
    const Eigen::DenseBase<Derived3>& r_enu,
    const Eigen::DenseBase<Derived4>& lla0) {
  Eigen::Vector3<typename Derived1::Scalar> a_eci;
  enu2ecia(dt, a_enu, v_enu, r_enu, lla0, a_eci);
  return a_eci;
}

};  // namespace nt

#endif