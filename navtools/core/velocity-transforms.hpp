#ifndef NAVTOOLS_CORE_VELOCITY_TRANSFORMS
#define NAVTOOLS_CORE_VELOCITY_TRANSFORMS

#include "navtools/core/constants.hpp"
#include "navtools/core/coordinate-dcms.hpp"
#include "navtools/core/macros.hpp"
#include "navtools/core/position-transforms.hpp"
#include <Eigen/Dense>

namespace nt {

/**
 * *=== ecef2nedv ===*
 * @brief Converts ECEF velocity to NED velocity
 * @param v_ecef Velocity vector in the ECEF frame (m/s)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @param v_ned The velocity vector in the NED frame (m/s)
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void ecef2nedv(
    const Eigen::DenseBase<Derived1>& v_ecef,
    const Eigen::DenseBase<Derived2>& lla0,
    Eigen::DenseBase<Derived3>& v_ned)
{
  ASSERT_EIGEN_VEC_SIZE(Derived1, v_ecef, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, v_ned, 3);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived2);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix3<Scalar> C_e_n;
  ecef2nedDcm(lla0, C_e_n);
  v_ned.derived() = C_e_n * v_ecef.derived();
}

/**
 * *=== ecef2nedv ===*
 * @brief Converts ECEF velocity to NED velocity
 * @param v_ecef Velocity vector in the ECEF frame (m/s)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @returns The velocity vector in the NED frame (m/s)
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> ecef2nedv(
    const Eigen::DenseBase<Derived1>& v_ecef, const Eigen::DenseBase<Derived2>& lla0)
{
  Eigen::Vector3<typename Derived1::Scalar> v_ned;
  ecef2nedv(v_ecef, lla0, v_ned);
  return v_ned;
}

/**
 * *=== ecef2enuv ===*
 * @brief Converts ECEF velocity to ENU velocity
 * @param v_ecef Velocity vector in the ECEF frame (m/s)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @param v_enu The velocity vector in the ENU frame (m/s)
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void ecef2enuv(
    const Eigen::DenseBase<Derived1>& v_ecef,
    const Eigen::DenseBase<Derived2>& lla0,
    Eigen::DenseBase<Derived3>& v_enu)
{
  ASSERT_EIGEN_VEC_SIZE(Derived1, v_ecef, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, v_enu, 3);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived2);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix3<Scalar> C_e_n;
  ecef2enuDcm(lla0, C_e_n);
  v_enu.derived() = C_e_n * v_ecef.derived();
}

/**
 * *=== ecef2enuv ===*
 * @brief Converts ECEF velocity to ENU velocity
 * @param v_ecef Velocity vector in the ECEF frame (m/s)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @returns The velocity vector in the ENU frame (m/s)
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> ecef2enuv(
    const Eigen::DenseBase<Derived1>& v_ecef, const Eigen::DenseBase<Derived2>& lla0)
{
  Eigen::Vector3<typename Derived1::Scalar> v_enu;
  ecef2enuv(v_ecef, lla0, v_enu);
  return v_enu;
}

/**
 * *=== ned2ecefv ===*
 * @brief Converts NED velocity to ECEF velocity.
 * @param v_ned The velocity vector in the NED frame (m/s)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @param v_ecef Velocity vector in the ECEF frame (m/s)
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void ned2ecefv(
    const Eigen::DenseBase<Derived2>& lla0,
    const Eigen::DenseBase<Derived1>& v_ned,
    Eigen::DenseBase<Derived3>& v_ecef)
{
  ASSERT_EIGEN_VEC_SIZE(Derived1, v_ned, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, v_ecef, 3);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived2);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix3<Scalar> C_n_e;
  ned2ecefDcm(lla0, C_n_e);
  v_ecef.derived() = C_n_e * v_ned.derived();
}

/**
 * *=== ned2ecefv ===*
 * @brief Converts NED velocity to ECEF velocity.
 * @param v_ned The velocity vector in the NED frame (m/s)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @returns The Velocity vector in the ECEF frame (m/s)
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> ned2ecefv(
    const Eigen::DenseBase<Derived2>& lla0, const Eigen::DenseBase<Derived1>& v_ned)
{
  Eigen::Vector3<typename Derived1::Scalar> v_ecef;
  ned2ecefv(lla0, v_ned, v_ecef);
  return v_ecef;
}

/**
 * *=== enu2ecefv ===*
 * @brief Converts ENU velocity to ECEF velocity.
 * @param v_enu The velocity vector in the ENU frame (m/s)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @param v_ecef Velocity vector in the ECEF frame (m/s)
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void enu2ecefv(
    const Eigen::DenseBase<Derived2>& lla0,
    const Eigen::DenseBase<Derived1>& v_enu,
    Eigen::DenseBase<Derived3>& v_ecef)
{
  ASSERT_EIGEN_VEC_SIZE(Derived1, v_enu, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, v_ecef, 3);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived2);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix3<Scalar> C_n_e;
  enu2ecefDcm(lla0, C_n_e);
  v_ecef.derived() = C_n_e * v_enu.derived();
}

/**
 * *=== enu2ecefv ===*
 * @brief Converts ENU velocity to ECEF velocity.
 * @param v_enu The velocity vector in the ENU frame (m/s)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @returns The Velocity vector in the ECEF frame (m/s)
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> enu2ecefv(
    const Eigen::DenseBase<Derived2>& lla0, const Eigen::DenseBase<Derived1>& v_enu)
{
  Eigen::Vector3<typename Derived1::Scalar> v_ecef;
  enu2ecefv(lla0, v_enu, v_ecef);
  return v_ecef;
}

/**
 * *=== ecef2eciv ===*
 * @brief Converts ECEF velocity to ECI velocity
 * @param r_ecef Position vector in the ECEF frame (m)
 * @param v_ecef Velocity vector in the ECEF frame (m/s)
 * @param dt Elapsed time in seconds
 * @param v_eci The velocity vector in the ECI frame (m/s)
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void ecef2eciv(
    const typename Derived1::Scalar& dt,
    const Eigen::DenseBase<Derived2>& r_ecef,
    const Eigen::DenseBase<Derived1>& v_ecef,
    Eigen::DenseBase<Derived3>& v_eci)
{
  ASSERT_EIGEN_VEC_SIZE(Derived1, v_ecef, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, r_ecef, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, v_eci, 3);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived2);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix3<Scalar> C_e_i;
  ecef2eciDcm(dt, C_e_i);
  v_eci.derived() = C_e_i * (v_ecef.derived() + OMEGA_ECEF<Scalar>.cross(r_ecef.derived()));
}

/**
 * *=== ecef2eciv ===*
 * @brief Converts ECEF velocity to ECI velocity
 * @param dt Elapsed time in seconds
 * @param v_ecef Velocity vector in the ECEF frame (m/s)
 * @param r_ecef Position vector in the ECEF frame (m)
 * @returns The velocity vector in the ECI frame (m/s)
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> ecef2eciv(
    const typename Derived1::Scalar& dt,
    const Eigen::DenseBase<Derived2>& r_ecef,
    const Eigen::DenseBase<Derived1>& v_ecef)
{
  Eigen::Vector3<typename Derived1::Scalar> v_eci;
  ecef2eciv(dt, r_ecef, v_ecef, v_eci);
  return v_eci;
}

/**
 * *=== eci2ecefv ===*
 * @brief Converts ECI velocity to ECEF velocity
 * @param dt Elapsed time in seconds
 * @param v_eci Velocity vector in the ECI frame (m/s)
 * @param r_eci Position vector in the ECI frame (m)
 * @param v_ecef The velocity vector in the ECEF frame (m/s)
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline void eci2ecefv(
    const typename Derived1::Scalar& dt,
    const Eigen::DenseBase<Derived2>& r_eci,
    const Eigen::DenseBase<Derived1>& v_eci,
    Eigen::DenseBase<Derived3>& v_ecef)
{
  ASSERT_EIGEN_VEC_SIZE(Derived1, v_eci, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, r_eci, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, v_ecef, 3);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived2);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived3);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix3<Scalar> C_i_e;
  eci2ecefDcm(dt, C_i_e);
  v_ecef.derived() = C_i_e * (v_eci.derived() - OMEGA_ECEF<Scalar>.cross(r_eci.derived()));
}

/**
 * *=== eci2ecefv ===*
 * @brief Converts ECI velocity to ECEF velocity
 * @param dt Elapsed time in seconds
 * @param v_eci Velocity vector in the ECI frame (m/s)
 * @returns The velocity vector in the ECEF frame (m/s)
 */
template <typename Derived1, typename Derived2>
inline Eigen::Vector3<typename Derived1::Scalar> eci2ecefv(
    const typename Derived1::Scalar& dt,
    const Eigen::DenseBase<Derived2>& r_eci,
    const Eigen::DenseBase<Derived1>& v_eci)
{
  Eigen::Vector3<typename Derived1::Scalar> v_ecef;
  eci2ecefv(dt, r_eci, v_eci, v_ecef);
  return v_ecef;
}

/**
 * *=== eci2nedv ===*
 * @brief Converts ECI velocity to NED velocity
 * @param dt Elapsed time in seconds
 * @param v_eci Velocity vector in the ECI frame (m/s)
 * @param r_eci Position vector in the ECI frame (m)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @param v_ned The velocity vector in the NED frame (m/s)
 */
template <typename Derived1, typename Derived2, typename Derived3, typename Derived4>
inline void eci2nedv(
    const typename Derived1::Scalar& dt,
    const Eigen::DenseBase<Derived2>& r_eci,
    const Eigen::DenseBase<Derived1>& v_eci,
    const Eigen::DenseBase<Derived3>& lla0,
    Eigen::DenseBase<Derived4>& v_ned)
{
  ASSERT_EIGEN_VEC_SIZE(Derived1, v_eci, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived1, r_eci, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived4, v_ned, 3);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived2);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived3);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived4);
  using Scalar = typename Derived1::Scalar;

  Eigen::Vector3<Scalar> v_ecef;
  eci2ecefv(dt, r_eci, v_eci, v_ecef);
  ecef2nedv(v_ecef, lla0, v_ned);
}

/**
 * *=== eci2nedv ===*
 * @brief Converts ECI velocity to NED velocity
 * @param dt Elapsed time in seconds
 * @param v_eci Velocity vector in the ECI frame (m/s)
 * @param r_eci Position vector in the ECI frame (m)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @returns The velocity vector in the NED frame (m/s)
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline Eigen::Vector3<typename Derived1::Scalar> eci2nedv(
    const typename Derived1::Scalar& dt,
    const Eigen::DenseBase<Derived2>& r_eci,
    const Eigen::DenseBase<Derived1>& v_eci,
    const Eigen::DenseBase<Derived3>& lla0)
{
  Eigen::Vector3<typename Derived1::Scalar> v_ned;
  eci2nedv(dt, r_eci, v_eci, lla0, v_ned);
  return v_ned;
}

/**
 * *=== eci2enuv ===*
 * @brief Converts ECI velocity to ENU velocity
 * @param dt Elapsed time in seconds
 * @param v_eci Velocity vector in the ECI frame (m/s)
 * @param r_eci Position vector in the ECI frame (m)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @param v_enu The velocity vector in the ENU frame (m/s)
 */
template <typename Derived1, typename Derived2, typename Derived3, typename Derived4>
inline void eci2enuv(
    const typename Derived1::Scalar& dt,
    const Eigen::DenseBase<Derived2>& r_eci,
    const Eigen::DenseBase<Derived1>& v_eci,
    const Eigen::DenseBase<Derived3>& lla0,
    Eigen::DenseBase<Derived4>& v_enu)
{
  ASSERT_EIGEN_VEC_SIZE(Derived1, v_eci, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived1, r_eci, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived4, v_enu, 3);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived2);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived3);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived4);
  using Scalar = typename Derived1::Scalar;

  Eigen::Vector3<Scalar> v_ecef;
  eci2ecefv(dt, r_eci, v_eci, v_ecef);
  ecef2enuv(v_ecef, lla0, v_enu);
}

/**
 * *=== eci2enuv ===*
 * @brief Converts ECI velocity to ENU velocity
 * @param dt Elapsed time in seconds
 * @param v_eci Velocity vector in the ECI frame (m/s)
 * @param r_eci Position vector in the ECI frame (m)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @returns The velocity vector in the ENU frame (m/s)
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline Eigen::Vector3<typename Derived1::Scalar> eci2enuv(
    const typename Derived1::Scalar& dt,
    const Eigen::DenseBase<Derived2>& r_eci,
    const Eigen::DenseBase<Derived1>& v_eci,
    const Eigen::DenseBase<Derived3>& lla0)
{
  Eigen::Vector3<typename Derived1::Scalar> v_enu;
  eci2enuv(dt, r_eci, v_eci, lla0, v_enu);
  return v_enu;
}

/**
 * *=== ned2eciv ===*
 * @brief Converts NED velocity to ECI velocity
 * @param dt Elapsed time in seconds
 * @param v_ned Velocity vector in the NED frame (m/s)
 * @param r_ned Position vector in the NED frame (m)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @param v_eci The velocity vector in the ECI frame (m/s)
 */
template <typename Derived1, typename Derived2, typename Derived3, typename Derived4>
inline void ned2eciv(
    const typename Derived1::Scalar& dt,
    const Eigen::DenseBase<Derived2>& r_ned,
    const Eigen::DenseBase<Derived1>& v_ned,
    const Eigen::DenseBase<Derived3>& lla0,
    Eigen::DenseBase<Derived4>& v_eci)
{
  ASSERT_EIGEN_VEC_SIZE(Derived1, v_ned, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, r_ned, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived4, v_eci, 3);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived2);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived3);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived4);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix3<Scalar> C_n_e;
  ned2ecefDcm(lla0, C_n_e);
  Eigen::Vector3<Scalar> v_xyz = C_n_e * v_ned.derived();
  Eigen::Vector3<Scalar> r_xyz0 = lla2ecef(lla0);
  Eigen::Vector3<Scalar> r_xyz = r_xyz0 + C_n_e * r_ned.derived();
  ecef2eciv(dt, r_xyz, v_xyz, v_eci);
}

/**
 * *=== ned2eciv ===*
 * @brief Converts NED velocity to ECI velocity
 * @param dt Elapsed time in seconds
 * @param v_ned Velocity vector in the NED frame (m/s)
 * @param r_ned Position vector in the NED frame (m)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @returns The velocity vector in the ECI frame (m/s)
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline Eigen::Vector3<typename Derived1::Scalar> ned2eciv(
    const typename Derived1::Scalar& dt,
    const Eigen::DenseBase<Derived2>& r_ned,
    const Eigen::DenseBase<Derived1>& v_ned,
    const Eigen::DenseBase<Derived3>& lla0)
{
  Eigen::Vector3<typename Derived1::Scalar> v_eci;
  ned2eciv(dt, r_ned, v_ned, lla0, v_eci);
  return v_eci;
}

/**
 * *=== enu2eciv ===*
 * @brief Converts ENU velocity to ECI velocity
 * @param dt Elapsed time in seconds
 * @param v_enu Velocity vector in the ENU frame (m/s)
 * @param r_enu Position vector in the ENU frame (m)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @param v_eci The velocity vector in the ECI frame (m/s)
 */
template <typename Derived1, typename Derived2, typename Derived3, typename Derived4>
inline void enu2eciv(
    const typename Derived1::Scalar& dt,
    const Eigen::DenseBase<Derived2>& r_enu,
    const Eigen::DenseBase<Derived1>& v_enu,
    const Eigen::DenseBase<Derived3>& lla0,
    Eigen::DenseBase<Derived4>& v_eci)
{
  ASSERT_EIGEN_VEC_SIZE(Derived1, v_enu, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, r_enu, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived3, lla0, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived4, v_eci, 3);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived2);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived3);
  ASSERT_EIGEN_SAME_SCALAR(Derived1, Derived4);
  using Scalar = typename Derived1::Scalar;

  Eigen::Matrix3<Scalar> C_n_e;
  enu2ecefDcm(lla0, C_n_e);
  Eigen::Vector3<Scalar> v_xyz = C_n_e * v_enu.derived();
  Eigen::Vector3<Scalar> r_xyz0 = lla2ecef(lla0);
  Eigen::Vector3<Scalar> r_xyz = r_xyz0 + C_n_e * r_enu.derived();
  ecef2eciv(dt, r_xyz, v_xyz, v_eci);
}

/**
 * *=== enu2eciv ===*
 * @brief Converts ENU velocity to ECI velocity
 * @param dt Elapsed time in seconds
 * @param v_enu Velocity vector in the ENU frame (m/s)
 * @param r_ned Position vector in the NED frame (m)
 * @param lla0  Reference point Latitude, Longitude, Altitude (rad,rad,m)
 * @returns The velocity vector in the ECI frame (m/s)
 */
template <typename Derived1, typename Derived2, typename Derived3>
inline Eigen::Vector3<typename Derived1::Scalar> enu2eciv(
    const typename Derived1::Scalar& dt,
    const Eigen::DenseBase<Derived2>& r_enu,
    const Eigen::DenseBase<Derived1>& v_enu,
    const Eigen::DenseBase<Derived3>& lla0)
{
  Eigen::Vector3<typename Derived1::Scalar> v_eci;
  enu2eciv(dt, r_enu, v_enu, lla0, v_eci);
  return v_eci;
}

};  // namespace nt

#endif
