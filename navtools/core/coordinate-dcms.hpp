#ifndef NAVTOOLS_CORE_COORDINATE_DCMS
#define NAVTOOLS_CORE_COORDINATE_DCMS

#include "navtools/core/constants.hpp"
#include "navtools/core/macros.hpp"
#include <Eigen/Dense>
#include <cmath>

namespace nt {

/**
 * *=== eci2ecefDcm ===*
 * @brief Earth-Centered-Inertial to Earth-Centered-Earth-Fixed rotation matrix
 * @tparam Derived The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param dt time elapsed between frames (s)
 * @param R  3x3 ECI->ECEF rotation matrix
 */
template <typename Derived>
inline void eci2ecefDcm(const typename Derived::Scalar dt, Eigen::DenseBase<Derived> &R) {
  ASSERT_EIGEN_MAT_SIZE(Derived, R, 3, 3);
  using Scalar = typename Derived::Scalar;

  Scalar omega_dt = WGS84_OMEGA<Scalar> * dt;
  Scalar swt = std::sin(omega_dt);
  Scalar cwt = std::cos(omega_dt);
  // clang-format off
  R.derived() <<  cwt, swt, 0.0, 
                 -swt, cwt, 0.0, 
                  0.0, 0.0, 1.0;
  // clang-format on
}

/**
 * *=== eci2ecefDcm ===*
 * @brief Earth-Centered-Inertial to Earth-Centered-Earth-Fixed rotation matrix
 * @tparam Derived The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param dt time elapsed between frames (s)
 * @returns the 3x3 ECI->ECEF rotation matrix
 */
template <typename Derived>
inline Derived eci2ecefDcm(const typename Derived::Scalar dt) {
  Derived R;
  eci2ecefDcm(dt, R);
  return R;
}

/**
 * *=== eci2nedDcm ===*
 * @brief Earth-Centered-Inertial to North-East-Down rotation matrix
 * @tparam DerivedVec The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param dt   time elapsed between frames (s)
 * @param lla0 3x1 reference Geodetic Latitude, Longitude, Height (rad, rad, m)
 * @param R    3x3 ECI->NED rotation matrix
 */
template <typename DerivedVec, typename DerivedMat>
inline void eci2nedDcm(
    const typename DerivedVec::Scalar dt,
    const Eigen::DenseBase<DerivedVec> &lla0,
    Eigen::DenseBase<DerivedMat> &R) {
  ASSERT_EIGEN_VEC_SIZE(DerivedVec, lla0, 3);
  ASSERT_EIGEN_MAT_SIZE(DerivedMat, R, 3, 3);
  ASSERT_EIGEN_TYPE(DerivedVec, DerivedMat);
  using Scalar = typename DerivedVec::Scalar;

  const auto &ref = lla0.derived();
  Scalar omega_dt = WGS84_OMEGA<Scalar> * dt;
  Scalar slat = std::sin(ref(0));
  Scalar clat = std::cos(ref(0));
  Scalar slongwt = std::sin(ref(1) + omega_dt);
  Scalar clongwt = std::cos(ref(1) + omega_dt);
  // clang-format off
  R.derived() << -slat * clongwt, -slat * slongwt,  clat,
                        -slongwt,         clongwt,   0.0,
                 -clat * clongwt, -clat * slongwt, -slat;
  // clang-format on
}

/**
 * *=== eci2nedDcm ===*
 * @brief Earth-Centered-Inertial to North-East-Down rotation matrix
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @tparam DerivedVec The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param dt   time elapsed between frames (s)
 * @param lla0 3x1 reference Geodetic Latitude, Longitude, Height (rad, rad, m)
 * @return the 3x3 ECI->NED rotation matrix
 */
template <typename DerivedMat, typename DerivedVec>
inline DerivedMat eci2nedDcm(
    const typename DerivedVec::Scalar dt, const Eigen::DenseBase<DerivedVec> &lla0) {
  DerivedMat R;
  eci2nedDcm(dt, lla0.derived(), R);
  return R;
}

/**
 * *=== eci2enuDcm ===*
 * @brief Earth-Centered-Inertial to East-North-Up rotation matrix
 * @tparam DerivedVec The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param lla0 3x1 reference Geodetic Latitude, Longitude, Height (rad, rad, m)
 * @param dt   time elapsed between frames (s)
 * @param R    3x3 ECI->ENU rotation matrix
 */
template <typename DerivedVec, typename DerivedMat>
inline void eci2enuDcm(
    const typename DerivedVec::Scalar dt,
    const Eigen::DenseBase<DerivedVec> &lla0,
    Eigen::DenseBase<DerivedMat> &R) {
  ASSERT_EIGEN_VEC_SIZE(DerivedVec, lla0, 3);
  ASSERT_EIGEN_MAT_SIZE(DerivedMat, R, 3, 3);
  ASSERT_EIGEN_TYPE(DerivedVec, DerivedMat);
  using Scalar = typename DerivedVec::Scalar;

  const auto &ref = lla0.derived();
  Scalar omega_dt = WGS84_OMEGA<Scalar> * dt;
  Scalar slat = std::sin(ref(0));
  Scalar clat = std::cos(ref(0));
  Scalar slongwt = std::sin(ref(1) + omega_dt);
  Scalar clongwt = std::cos(ref(1) + omega_dt);
  // clang-format off
  R.derived() <<        -slongwt,         clongwt,  0.0, 
                 -slat * clongwt, -slat * slongwt, clat, 
                  clat * clongwt,  clat * slongwt, slat;
  // clang-format on
}

/**
 * *=== eci2enuDcm ===*
 * @brief Earth-Centered-Inertial to East-North-Up rotation matrix
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @tparam DerivedVec The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param lla0 3x1 reference Geodetic Latitude, Longitude, Height (rad, rad, m)
 * @param dt   time elapsed between frames (s)
 * @returns the 3x3 ECI->ENU rotation matrix
 */
template <typename DerivedMat, typename DerivedVec>
inline DerivedMat eci2enuDcm(
    const typename DerivedVec::Scalar dt, const Eigen::DenseBase<DerivedVec> &lla0) {
  DerivedMat R;
  eci2enuDcm(dt, lla0.derived(), R);
  return R;
}

/**
 * *=== ecef2eciDcm ===*
 * @brief Earth-Centered-Earth-Fixed to Earth-Centered-Inertial rotation matrix
 * @tparam Derived The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param dt time elapsed between frames (s)
 * @param R  3x3 ECEF->ECI rotation matrix
 */
template <typename Derived>
inline void ecef2eciDcm(const typename Derived::Scalar dt, Eigen::DenseBase<Derived> &R) {
  ASSERT_EIGEN_MAT_SIZE(Derived, R, 3, 3);
  using Scalar = typename Derived::Scalar;

  Scalar omega_dt = WGS84_OMEGA<Scalar> * dt;
  Scalar swt = std::sin(omega_dt);
  Scalar cwt = std::cos(omega_dt);
  // clang-format off
  R.derived() << cwt, -swt, 0.0, 
                 swt,  cwt, 0.0, 
                 0.0,  0.0, 1.0;
  // clang-format on
}

/**
 * *=== ecef2eciDcm ===*
 * @brief Earth-Centered-Earth-Fixed to Earth-Centered-Inertial rotation matrix
 * @tparam Derived The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param dt time elapsed between frames (s)
 * @returns the 3x3 ECEF->ECI rotation matrix
 */
template <typename Derived>
inline Derived ecef2eciDcm(const typename Derived::Scalar dt) {
  Derived R;
  ecef2eciDcm(dt, R);
  return R;
}

/**
 * *=== ecef2nedDcm ===*
 * @brief Earth-Centered-Earth-Fixed to North-East-Down rotation matrix
 * @tparam DerivedVec The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param lla0 3x1 reference Geodetic Latitude, Longitude, Height (rad, rad, m)
 * @param R    3x3 ECEF->NED rotation matrix
 */
template <typename DerivedVec, typename DerivedMat>
inline void ecef2nedDcm(const Eigen::DenseBase<DerivedVec> &lla0, Eigen::DenseBase<DerivedMat> &R) {
  ASSERT_EIGEN_VEC_SIZE(DerivedVec, lla0, 3);
  ASSERT_EIGEN_MAT_SIZE(DerivedMat, R, 3, 3);
  ASSERT_EIGEN_TYPE(DerivedVec, DerivedMat);
  using Scalar = typename DerivedVec::Scalar;

  const auto &ref = lla0.derived();
  Scalar slat = std::sin(ref(0));
  Scalar clat = std::cos(ref(0));
  Scalar slong = std::sin(ref(1));
  Scalar clong = std::cos(ref(1));

  // clang-format off
  R.derived() << -slat * clong, -slat * slong,  clat,
                        -slong,         clong,   0.0,
                 -clat * clong, -clat * slong, -slat;
  // clang-format on
}

/**
 * *=== ecef2nedDcm ===*
 * @brief Earth-Centered-Earth-Fixed to North-East-Down rotation matrix
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @tparam DerivedVec The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param lla0 3x1 reference Geodetic Latitude, Longitude, Height (rad, rad, m)
 * @returns the 3x3 ECEF->NED rotation matrix
 */
template <typename DerivedMat, typename DerivedVec>
inline DerivedMat ecef2nedDcm(const Eigen::DenseBase<DerivedVec> &lla0) {
  DerivedMat R;
  ecef2nedDcm(lla0.derived(), R);
  return R;
}

/**
 * *=== ecef2enuDcm ===*
 * @brief Earth-Centered-Earth-Fixed to East-North-Up rotation matrix
 * @tparam DerivedVec The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param lla0 3x1 reference Geodetic Latitude, Longitude, Height (rad, rad, m)
 * @param R    3x3 ECEF->ENU rotation matrix
 */
template <typename DerivedVec, typename DerivedMat>
inline void ecef2enuDcm(const Eigen::DenseBase<DerivedVec> &lla0, Eigen::DenseBase<DerivedMat> &R) {
  ASSERT_EIGEN_VEC_SIZE(DerivedVec, lla0, 3);
  ASSERT_EIGEN_MAT_SIZE(DerivedMat, R, 3, 3);
  ASSERT_EIGEN_TYPE(DerivedVec, DerivedMat);
  using Scalar = typename DerivedVec::Scalar;

  const auto &ref = lla0.derived();
  Scalar slat = std::sin(ref(0));
  Scalar clat = std::cos(ref(0));
  Scalar slong = std::sin(ref(1));
  Scalar clong = std::cos(ref(1));

  // clang-format off
  R.derived() <<        -slong,        clong,   0.0,
                 -slat * clong, -slat * slong, clat,
                  clat * clong,  clat * slong, slat;
  // clang-format on
}

/**
 * *=== ecef2enuDcm ===*
 * @brief Earth-Centered-Earth-Fixed to East-North-Up rotation matrix
 * @tparam DerivedVec The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param lla0 3x1 reference Geodetic Latitude, Longitude, Height (rad, rad, m)
 * @param R    3x3 ECEF->ENU rotation matrix
 */
template <typename DerivedMat, typename DerivedVec>
inline DerivedMat ecef2enuDcm(const Eigen::DenseBase<DerivedVec> &lla0) {
  DerivedMat R;
  ecef2enuDcm(lla0.derived(), R);
  return R;
}

/**
 * *=== ned2eciDcm ===*
 * @brief North-East-Down to Earth-Centered-Inertial rotation matrix
 * @tparam DerivedVec The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param lla0 3x1 reference Geodetic Latitude, Longitude, Height (rad, rad, m)
 * @param dt   time elapsed between frames (s)
 * @param R    3x3 NED->ECI rotation matrix
 */
template <typename DerivedVec, typename DerivedMat>
inline void ned2eciDcm(
    const typename DerivedVec::Scalar dt,
    const Eigen::DenseBase<DerivedVec> &lla0,
    Eigen::DenseBase<DerivedMat> &R) {
  ASSERT_EIGEN_VEC_SIZE(DerivedVec, lla0, 3);
  ASSERT_EIGEN_MAT_SIZE(DerivedMat, R, 3, 3);
  ASSERT_EIGEN_TYPE(DerivedVec, DerivedMat);
  using Scalar = typename DerivedVec::Scalar;

  const auto &ref = lla0.derived();
  Scalar omega_dt = WGS84_OMEGA<Scalar> * dt;
  Scalar slat = std::sin(ref(0));
  Scalar clat = std::cos(ref(0));
  Scalar slongwt = std::sin(ref(1) + omega_dt);
  Scalar clongwt = std::cos(ref(1) + omega_dt);
  // clang-format off
  R.derived() << -slat * clongwt, -slongwt, -clat * clongwt,
                 -slat * slongwt,  clongwt, -clat * slongwt,
                            clat,      0.0,           -slat;
  // clang-format on
}

/**
 * *=== ned2eciDcm ===*
 * @brief North-East-Down to Earth-Centered-Inertial rotation matrix
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @tparam DerivedVec The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param lla0 3x1 reference Geodetic Latitude, Longitude, Height (rad, rad, m)
 * @param dt  time elapsed between frames (s)
 * @returns the 3x3 NED->ECI rotation matrix
 */
template <typename DerivedMat, typename DerivedVec>
inline DerivedMat ned2eciDcm(
    const typename DerivedVec::Scalar dt, const Eigen::DenseBase<DerivedVec> &lla0) {
  DerivedMat R;
  ned2eciDcm(dt, lla0.derived(), R);
  return R;
}

/**
 * *=== ned2ecefDcm ===*
 * @brief North-East-Down to Earth-Centered-Earth-Fixed rotation matrix
 * @tparam DerivedVec The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param lla0 3x1 reference Geodetic Latitude, Longitude, Height (rad, rad, m)
 * @param R    3x3 NED->ECEF rotation matrix
 */
template <typename DerivedVec, typename DerivedMat>
inline void ned2ecefDcm(const Eigen::DenseBase<DerivedVec> &lla0, Eigen::DenseBase<DerivedMat> &R) {
  ASSERT_EIGEN_MAT_SIZE(DerivedMat, R, 3, 3);
  ASSERT_EIGEN_VEC_SIZE(DerivedVec, lla0, 3);
  ASSERT_EIGEN_TYPE(DerivedVec, DerivedMat);
  using Scalar = typename DerivedVec::Scalar;

  const auto &ref = lla0.derived();
  Scalar slat = std::sin(ref(0));
  Scalar clat = std::cos(ref(0));
  Scalar slong = std::sin(ref(1));
  Scalar clong = std::cos(ref(1));
  // clang-format off
  R.derived() << -slat * clong, -slong, -clat * clong,
                 -slat * slong,  clong, -clat * slong,
                          clat,    0.0,         -slat;
  // clang-format on
}

/**
 * *=== ned2ecefDcm ===*
 * @brief North-East-Down to Earth-Centered-Earth-Fixed rotation matrix
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @tparam DerivedVec The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param lla 3x1 Geodetic Latitude, Longitude, Height (rad, rad, m)
 * @returns the 3x3 NED->ECEF rotation matrix
 */
template <typename DerivedMat, typename DerivedVec>
inline DerivedMat ned2ecefDcm(const Eigen::DenseBase<DerivedVec> &lla0) {
  DerivedMat R;
  ned2ecefDcm(lla0.derived(), R);
  return R;
}

/**
 * *=== ned2enuDcm ===*
 * @brief North-East-Down to East-North-Up rotation matrix
 * @tparam Derived The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param R 3x3 NED->ENU rotation matrix
 */
template <typename Derived>
inline void ned2enuDcm(Eigen::DenseBase<Derived> &R) {
  ASSERT_EIGEN_MAT_SIZE(Derived, R, 3, 3);
  R.derived() << 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0;
}

/**
 * *=== ned2enuDcm ===*
 * @brief North-East-Down to East-North-Up rotation matrix
 * @tparam Derived The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @returns 3x3 NED->ENU rotation matrix
 */
template <typename Derived>
inline Derived ned2enuDcm() {
  Derived R;
  ned2enuDcm(R);
  return R;
}

/**
 * *=== enu2eciDcm ===*
 * @brief East-North-Up to Earth-Centered-Inertial rotation matrix
 * @tparam DerivedVec The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param lla0 3x1 reference Geodetic Latitude, Longitude, Height (rad, rad, m)
 * @param dt   time elapsed between frames (s)
 * @param R    3x3 ENU->ECI rotation matrix
 */
template <typename DerivedVec, typename DerivedMat>
inline void enu2eciDcm(
    const typename DerivedVec::Scalar dt,
    const Eigen::DenseBase<DerivedVec> &lla0,
    Eigen::DenseBase<DerivedMat> &R) {
  ASSERT_EIGEN_VEC_SIZE(DerivedVec, lla0, 3);
  ASSERT_EIGEN_MAT_SIZE(DerivedMat, R, 3, 3);
  ASSERT_EIGEN_TYPE(DerivedVec, DerivedMat);
  using Scalar = typename DerivedVec::Scalar;

  const auto &ref = lla0.derived();
  Scalar omega_dt = WGS84_OMEGA<Scalar> * dt;
  Scalar slat = std::sin(ref(0));
  Scalar clat = std::cos(ref(0));
  Scalar slongwt = std::sin(ref(1) + omega_dt);
  Scalar clongwt = std::cos(ref(1) + omega_dt);
  // clang-format off
  R.derived() << -slongwt, -slat * clongwt, clat * clongwt,
                  clongwt, -slat * slongwt, clat * slongwt,
                      0.0,            clat,           slat;
  // clang-format on
}

/**
 * *=== enu2eciDcm ===*
 * @brief East-North-Up to Earth-Centered-Inertial rotation matrix
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @tparam DerivedVec The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param lla0 3x1 reference Geodetic Latitude, Longitude, Height (rad, rad, m)
 * @param dt  time elapsed between frames (s)
 * @returns the 3x3 ENU->ECI rotation matrix
 */
template <typename DerivedMat, typename DerivedVec>
inline DerivedMat enu2eciDcm(
    const typename DerivedVec::Scalar dt, const Eigen::DenseBase<DerivedVec> &lla0) {
  DerivedMat R;
  enud2eciDcm(dt, lla0.derived(), R);
  return R;
}

/**
 * *=== enu2ecefDcm ===*
 * @brief East-North-Up to Earth-Centered-Earth-Fixed rotation matrix
 * @tparam DerivedVec The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param lla0 3x1 reference Geodetic Latitude, Longitude, Height (rad, rad, m)
 * @param R    3x3 ENU->ECEF rotation matrix
 */
template <typename DerivedVec, typename DerivedMat>
inline void enu2ecefDcm(const Eigen::DenseBase<DerivedVec> &lla0, Eigen::DenseBase<DerivedMat> &R) {
  ASSERT_EIGEN_MAT_SIZE(DerivedMat, R, 3, 3);
  ASSERT_EIGEN_VEC_SIZE(DerivedVec, lla0, 3);
  ASSERT_EIGEN_TYPE(DerivedVec, DerivedMat);
  using Scalar = typename DerivedVec::Scalar;

  const auto &ref = lla0.derived();
  Scalar slat = std::sin(ref(0));
  Scalar clat = std::cos(ref(0));
  Scalar slong = std::sin(ref(1));
  Scalar clong = std::cos(ref(1));
  // clang-format off
  R.derived() << -slong, -clong * slat, clong * clat,
                  clong, -slong * slat, slong * clat,
                    0.0,          clat,         slat;
  // clang-format on
}

/**
 * *=== enu2ecefDcm ===*
 * @brief East-North-Up to Earth-Centered-Earth-Fixed rotation matrix
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @tparam DerivedVec The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param lla 3x1 Geodetic Latitude, Longitude, Height (rad, rad, m)
 * @returns the 3x3 ENU->ECEF rotation matrix
 */
template <typename DerivedMat, typename DerivedVec>
inline DerivedMat enu2ecefDcm(const Eigen::DenseBase<DerivedVec> &lla0) {
  DerivedMat R;
  enu2ecefDcm(lla0.derived(), R);
  return R;
}

/**
 * *=== enu2nedDcm ===*
 * @brief East-North-Up to North-East-Down rotation matrix
 * @tparam Derived The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param R 3x3 ENU->NED rotation matrix
 */
template <typename Derived>
inline void enu2nedDcm(Eigen::DenseBase<Derived> &R) {
  ASSERT_EIGEN_MAT_SIZE(Derived, R, 3, 3);
  R.derived() << 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0;
}

/**
 * *=== ned2enuDcm ===*
 * @brief East-North-Up to North-East-Down rotation matrix
 * @tparam Derived The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @returns 3x3 ENU->NED rotation matrix
 */
template <typename Derived>
inline Derived enu2nedDcm() {
  Derived R;
  ned2enuDcm(R);
  return R;
}

}  // namespace nt

#endif