#ifndef GLYPH_CORE_ANGULAR_RATE_TRANSFORMS
#define GLYPH_CORE_ANGULAR_RATE_TRANSFORMS

#include "navtools/core/coordinate-dcms.hpp"
#include "navtools/core/macros.hpp"
#include <Eigen/Dense>

namespace nt {

/**
 * *=== ecef2nedw ===*
 * @brief Converts ECEF angular velocity to NED angular velocity
 * @param w_xyz Angular velocity in the ECEF frame (rad/s)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @param w_ned The NED angular velocity w.r.t the ECI frame, expressed in the NED frame
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void ecef2nedw(
    const Eigen::DenseBase<Derived1>& w_xyz,
    const Eigen::DenseBase<Derived2>& lla0,
    Eigen::DenseBase<Derived3>& w_ned) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, w_xyz, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, w_ned, 3);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived2);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix3<Scalar> C_e_n;
  ecef2nedDcm(lla0, C_e_n);
  w_ned.derived() = C_e_n * w_xyz.derived();
}

/**
 * *=== ecef2nedw ===*
 * @brief Converts ECEF angular velocity to NED angular velocity
 * @param w_xyz Angular velocity in the ECEF frame (rad/s)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @returns The NED angular velocity w.r.t the ECI frame, expressed in the NED frame
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> ecef2nedw(
    const Eigen::DenseBase<Derived1>& w_xyz, const Eigen::DenseBase<Derived2>& lla0) {
  Eigen::Vector3<typename Derived1::Scalar> w_ned;
  ecef2nedw(w_xyz, lla0, w_ned);
  return w_ned;
}

/**
 * *=== ecef2enuw ===*
 * @brief Converts ECEF angular velocity to ENU angular velocity
 * @param w_xyz Angular velocity in the ECEF frame (rad/s)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @param w_enu The ENU angular velocity w.r.t the ECI frame, expressed in the ENU frame
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void ecef2enuw(
    const Eigen::DenseBase<Derived1>& w_xyz,
    const Eigen::DenseBase<Derived2>& lla0,
    Eigen::DenseBase<Derived3>& w_enu) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, w_xyz, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, w_enu, 3);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived2);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix3<Scalar> C_e_n;
  ecef2enuDcm(lla0, C_e_n);
  w_enu.derived() = C_e_n * w_xyz.derived();
}

/**
 * *=== ecef2enuw ===*
 * @brief Converts ECEF angular velocity to ENU angular velocity
 * @param w_xyz Angular velocity in the ECEF frame (rad/s)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @returns The ENU angular velocity w.r.t the ECI frame, expressed in the ENU frame
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> ecef2enuw(
    const Eigen::DenseBase<Derived1>& w_xyz, const Eigen::DenseBase<Derived2>& lla0) {
  Eigen::Vector3<typename Derived1::Scalar> w_enu;
  ecef2enuw(w_xyz, lla0, w_enu);
  return w_enu;
}

/**
 * *=== ned2ecefw ===*
 * @brief Converts NED angular velocity to ECEF angular velocity
 * @param w_ned The NED angular
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @param w_xyz Earth's angular velocity in the ECEF frame (rad/s)
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void ned2ecefw(
    const Eigen::DenseBase<Derived1>& w_ned,
    const Eigen::DenseBase<Derived2>& lla0,
    Eigen::DenseBase<Derived3>& w_xyz) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, w_xyz, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, w_ned, 3);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived2);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix3<Scalar> C_n_e;
  ned2ecefDcm(lla0, C_n_e);
  w_xyz.derived() = C_n_e * w_ned.derived();
}

/**
 * *=== ned2ecefw ===*
 * @brief Converts NED angular velocity to ECEF angular velocity
 * @param w_ned The NED angular velocity
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @returns The Earth's angular velocity in the ECEF frame (rad/s)
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> ned2ecefw(
    const Eigen::DenseBase<Derived1>& w_ned, const Eigen::DenseBase<Derived2>& lla0) {
  Eigen::Vector3<typename Derived1::Scalar> w_xyz;
  ned2ecefw(w_ned, lla0, w_xyz);
  return w_xyz;
}

/**
 * *=== enu2ecefw ===*
 * @brief Converts ENU angular velocity to ECEF angular velocity
 * @param w_enu The ENU angular velocity
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @param w_xyz Earth's angular velocity in the ECEF frame (rad/s)
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void enu2ecefw(
    const Eigen::DenseBase<Derived1>& w_enu,
    const Eigen::DenseBase<Derived2>& lla0,
    Eigen::DenseBase<Derived3>& w_xyz) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, w_xyz, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, w_enu, 3);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived2);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix3<Scalar> C_n_e;
  enu2ecefDcm(lla0, C_n_e);
  w_xyz.derived() = C_n_e * w_enu.derived();
}

/**
 * *=== enu2ecefw ===*
 * @brief Converts ENU angular velocity to ECEF angular velocity
 * @param w_enu The ENU angular velocity
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @returns The Earth's angular velocity in the ECEF frame (rad/s)
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> enu2ecefw(
    const Eigen::DenseBase<Derived1>& w_enu, const Eigen::DenseBase<Derived2>& lla0) {
  Eigen::Vector3<typename Derived1::Scalar> w_xyz;
  enu2ecefw(w_enu, lla0, w_xyz);
  return w_xyz;
}

/**
 * *=== ecef2eciw ===*
 * @brief Converts ECEF angular velocity to ECI angular velocity
 * @param dt    Time elapsed (s)
 * @param w_xyz Angular velocity in the ECEF frame (rad/s)
 * @param w_eci The ECI angular velocity w.r.t the ECI frame, expressed in the ECI frame (rad/s)
 */
template <typename Derived1, typename Derived2>
inline void ecef2eciw(
    const typename Derived1::Scalar& dt,
    const Eigen::DenseBase<Derived1>& w_xyz,
    Eigen::DenseBase<Derived2>& w_eci) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, w_xyz, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, w_eci, 3);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived2);
  using Scalar = typename Derived2::Scalar;

  Eigen::Matrix3<Scalar> C_e_i;
  ecef2eciDcm(dt, C_e_i);
  w_eci.derived() = C_e_i * (w_xyz.derived() + OMEGA_ECEF<Scalar>);
}

/**
 * *=== ecef2eciw ===*
 * @brief Converts ECEF angular velocity to ECI angular velocity
 * @param dt    Time elapsed (s)
 * @param w_xyz Earth's angular velocity in the ECEF frame (rad/s)
 * @returns The ECI angular velocity w.r.t the ECI frame, expressed in the ECI frame (rad/s)
 */
template <typename Derived>
inline Eigen::Vector3<typename Derived::Scalar> ecef2eciw(
    const typename Derived::Scalar& dt, const Eigen::DenseBase<Derived>& w_xyz) {
  Eigen::Vector3<typename Derived::Scalar> w_eci;
  ecef2eciw(dt, w_xyz, w_eci);
  return w_eci;
}

/**
 * *=== eci2ecefw ===*
 * @brief Converts ECI angular velocity to ECEF angular velocity
 * @param dt    Time elapsed (s)
 * @param w_eci The ECI angular velocity w.r.t the ECI frame, expressed in the ECI frame (rad/s)
 * @param w_xyz Earth's angular velocity in the ECEF frame (rad/s)
 */
template <typename Derived1, typename Derived2>
inline void eci2ecefw(
    const typename Derived1::Scalar& dt,
    const Eigen::DenseBase<Derived1>& w_eci,
    Eigen::DenseBase<Derived2>& w_xyz) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, w_eci, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, w_xyz, 3);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived2);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix3<Scalar> C_i_e;
  eci2ecefDcm(dt, C_i_e);
  w_xyz.derived() = C_i_e * (w_eci.derived() - OMEGA_ECEF<Scalar>);
}

/**
 * *=== eci2ecefw ===*
 * @brief Converts ECI angular velocity to ECEF angular velocity
 * @param dt    Time elapsed (s)
 * @param w_eci The ECI angular velocity w.r.t the ECI frame, expressed in the ECI frame (rad/s)
 * @returns Earth's angular velocity in the ECEF frame (rad/s)
 */
template <typename Derived>
inline Eigen::Vector3<typename Derived::Scalar> eci2ecefw(
    const typename Derived::Scalar& dt, const Eigen::DenseBase<Derived>& w_eci) {
  Eigen::Vector3<typename Derived::Scalar> w_xyz;
  eci2ecefw(dt, w_eci, w_xyz);
  return w_xyz;
}

/**
 * *=== eci2nedw ===*
 * @brief Converts ECI angular velocity to NED angular velocity
 * @param dt    Time elapsed (s)
 * @param w_eci The ECI angular velocity w.r.t the ECI frame, expressed in the ECI frame (rad/s)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @param w_ned The Earth's angular velocity w.r.t the ECI frame, expressed in the NED frame
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void eci2nedw(
    const typename Derived1::Scalar& dt,
    const Eigen::DenseBase<Derived1>& w_eci,
    const Eigen::DenseBase<Derived2>& lla0,
    Eigen::DenseBase<Derived3>& w_ned) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, w_eci, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, w_ned, 3);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived2);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Vector3<Scalar> w_xyz;
  eci2ecefw(dt, w_eci, w_xyz);
  ecef2nedw(w_xyz, lla0, w_ned);
}

/**
 * *=== eci2nedw ===*
 * @brief Converts ECI angular velocity to NED angular velocity
 * @param dt    Time elapsed (s)
 * @param w_eci The ECI angular velocity w.r.t the ECI frame, expressed in the ECI frame (rad/s)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @returns The Earth's angular velocity w.r.t the ECI frame, expressed in the NED frame
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> eci2nedw(
    const typename Derived1::Scalar& dt,
    const Eigen::DenseBase<Derived1>& w_eci,
    const Eigen::DenseBase<Derived2>& lla0) {
  Eigen::Vector3<typename Derived1::Scalar> w_ned;
  eci2nedw(dt, w_eci, lla0, w_ned);
  return w_ned;
}

/**
 * *=== eci2enuw ===*
 * @brief Converts ECI angular velocity to ENU angular velocity
 * @param dt    Time elapsed (s)
 * @param w_eci The ECI angular velocity w.r.t the ECI frame, expressed in the ECI frame (rad/s)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @param w_enu The Earth's angular velocity w.r.t the ECI frame, expressed in the ENU frame
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void eci2enuw(
    const typename Derived1::Scalar& dt,
    const Eigen::DenseBase<Derived1>& w_eci,
    const Eigen::DenseBase<Derived2>& lla0,
    Eigen::DenseBase<Derived3>& w_enu) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, w_eci, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, w_enu, 3);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived2);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Vector3<Scalar> w_xyz;
  eci2ecefw(dt, w_eci, w_xyz);
  ecef2enuw(w_xyz, lla0, w_enu);
}

/**
 * *=== eci2enuw ===*
 * @brief Converts ECI angular velocity to ENU angular velocity
 * @param dt    Time elapsed (s)
 * @param w_eci The ECI angular velocity w.r.t the ECI frame, expressed in the ECI frame (rad/s)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @returns The Earth's angular velocity w.r.t the ECI frame, expressed in the ENU frame
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> eci2enuw(
    const typename Derived1::Scalar& dt,
    const Eigen::DenseBase<Derived1>& w_eci,
    const Eigen::DenseBase<Derived2>& lla0) {
  Eigen::Vector3<typename Derived1::Scalar> w_enu;
  eci2enuw(dt, w_eci, lla0, w_enu);
  return w_enu;
}

/**
 * *=== ned2eciw ===*
 * @brief Converts NED angular velocity to ECI angular velocity
 * @param dt    Time elapsed (s)
 * @param w_ned The Earth's angular velocity w.r.t the ECI frame, expressed in the NED frame
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @param w_eci The ECI angular velocity w.r.t the ECI frame, expressed in the ECI frame (rad/s)
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void ned2eciw(
    const typename Derived1::Scalar& dt,
    const Eigen::DenseBase<Derived1>& w_ned,
    const Eigen::DenseBase<Derived2>& lla0,
    Eigen::DenseBase<Derived3>& w_eci) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, w_ned, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, w_eci, 3);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived2);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Vector3<Scalar> w_xyz;
  ned2ecefw(w_ned, lla0, w_xyz);
  ecef2eciw(dt, w_xyz, w_eci);
}

/**
 * *=== ned2eciw ===*
 * @brief Converts NED angular velocity to ECI angular velocity
 * @param dt    Time elapsed (s)
 * @param w_ned The Earth's angular velocity w.r.t the ECI frame, expressed in the NED frame
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @returns The ECI angular velocity w.r.t the ECI frame, expressed in the ECI frame (rad/s)
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> ned2eciw(
    const typename Derived1::Scalar& dt,
    const Eigen::DenseBase<Derived1>& w_ned,
    const Eigen::DenseBase<Derived2>& lla0) {
  Eigen::Vector3<typename Derived1::Scalar> w_eci;
  ned2eciw(dt, w_ned, lla0, w_eci);
  return w_eci;
}

/**
 * *=== enu2eciw ===*
 * @brief Converts ENU angular velocity to ECI angular velocity
 * @param dt    Time elapsed (s)
 * @param w_enu The Earth's angular velocity w.r.t the ECI frame, expressed in the ENU frame
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @param w_eci The ECI angular velocity w.r.t the ECI frame, expressed in the ECI frame (rad/s)
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void enu2eciw(
    const typename Derived1::Scalar& dt,
    const Eigen::DenseBase<Derived1>& w_enu,
    const Eigen::DenseBase<Derived2>& lla0,
    Eigen::DenseBase<Derived3>& w_eci) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, w_enu, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, w_eci, 3);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived2);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Vector3<Scalar> w_xyz;
  enu2ecefw(w_enu, lla0, w_xyz);
  ecef2eciw(dt, w_xyz, w_eci);
}

/**
 * *=== enu2eciw ===*
 * @brief Converts ENU angular velocity to ECI angular velocity
 * @param dt    Time elapsed (s)
 * @param w_enu The Earth's angular velocity w.r.t the ECI frame, expressed in the ENU frame
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @returns The ECI angular velocity w.r.t the ECI frame, expressed in the ECI frame (rad/s)
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> enu2eciw(
    const typename Derived1::Scalar& dt,
    const Eigen::DenseBase<Derived1>& w_enu,
    const Eigen::DenseBase<Derived2>& lla0) {
  Eigen::Vector3<typename Derived1::Scalar> w_eci;
  enu2eciw(dt, w_enu, lla0, w_eci);
  return w_eci;
}

};  // namespace nt

#endif
