#ifndef NAVTOOLS_CORE_WGS84_MODELS_HPP
#define NAVTOOLS_CORE_WGS84_MODELS_HPP

#include "navtools/core/constants.hpp"
#include "navtools/core/macros.hpp"
#include <Eigen/Dense>

namespace nt {

/**
 * *=== TransverseRadius ===*
 * @brief Calculates the transverse radius relative to user latitude
 * @tparam T A floating point type (e.g. double, float)
 * @param phi Latitude (rad)
 * @returns Earth's transverse radius at "phi"
 */
template <typename T = double>
inline T TransverseRadius(const T phi) {
  T sin_phi = std::sin(phi);
  T t = 1.0 - WGS84_E2<T> * sin_phi * sin_phi;
  return WGS84_A<T> / std::sqrt(t);
}

/**
 * *=== MeridianRadius ===*
 * @brief Calculates the meridian radius relative to user latitude
 * @tparam T A floating point type (e.g. double, float)
 * @param phi Latitude (rad)
 * @returns Earth's meridian radius at "phi"
 */
template <typename T = double>
inline T MeridianRadius(const T phi) {
  T sin_phi = std::sin(phi);
  T t = 1.0 - WGS84_E2<T> * sin_phi * sin_phi;
  return WGS84_A<T> * (1.0 - WGS84_E2<T>) / std::pow(t, 1.5);
}

/**
 * *=== GeocentricRadius ===
 * @brief Calculates the geocentric radius relative to user latitude
 * @tparam T A floating point type (e.g. double, float)
 * @param phi Latitude (rad)
 * @returns Earth's geocentric radius at "phi"
 */
template <typename T = double>
inline T GeocentricRadius(const T phi) {
  T sin_phi2 = std::sin(phi);
  sin_phi2 *= sin_phi2;
  T cos_phi = std::cos(phi);
  T t = 1.0 - WGS84_E2<T> * sin_phi2;
  T Re = WGS84_A<T> / std::sqrt(t);
  T o_e2 = 1.0 - WGS84_E2<T>;
  return Re * std::sqrt(cos_phi * cos_phi + o_e2 * o_e2 * sin_phi2);
}

/**
 * *=== EarthRadii2 ===
 * @brief Calculates the {Transverse, Meridian} radii relative to user latitude
 * @tparam T A floating point type (e.g. double, float)
 * @param phi Latitude (rad)
 * @returns Earth's {Transverse, Meridian} radii at "phi"
 */
template <typename T = double>
inline std::pair<T, T> EarthRadii2(const T phi) {
  T sin_phi = std::sin(phi);
  T t = 1.0 - WGS84_E2<T> * sin_phi * sin_phi;

  return std::make_pair(
      WGS84_A<T> / std::sqrt(t), WGS84_A<T> * (1.0 - WGS84_E2<T>) / std::pow(t, 1.5));
}

/**
 * *=== EarthRadii3 ===
 * @brief Calculates the {Transverse, Meridian, Geocentric} radii relative to user latitude
 * @tparam T A floating point type (e.g. double, float)
 * @param phi Latitude (rad)
 * @returns Earth's {Transverse, Meridian, Geocentric} radii at Latitude
 */
template <typename T = double>
inline std::tuple<T, T, T> EarthRadii3(const T phi) {
  T sin_phi2 = std::sin(phi);
  sin_phi2 *= sin_phi2;
  T cos_phi = std::cos(phi);
  T t = 1.0 - WGS84_E2<T> * sin_phi2;
  T o_e2 = 1.0 - WGS84_E2<T>;

  T Re = WGS84_A<T> / std::sqrt(t);
  return std::make_tuple(
      Re,
      WGS84_A<T> * o_e2 / std::pow(t, 1.5),
      Re * std::sqrt(cos_phi * cos_phi + o_e2 * o_e2 * sin_phi2));
}

/**
 * *=== EarthRate ===
 * @brief Rotation rate of the earth relative to the specified local navigation frame
 * @tparam isNed Specify the local navigation frame as NED (true) or ENU (false)
 * @tparam T A floating point type (e.g. double, float)
 * @param phi Latitude (rad)
 * @returns Earth's rotation in the local navigation frame
 */
template <bool isNed = true, typename T = double>
inline auto EarthRate(const T phi) {
  if constexpr (isNed) {
    return Eigen::Vector3<T>(WGS84_OMEGA<T> * std::cos(phi), 0.0, WGS84_OMEGA<T> * std::sin(phi));
  } else {
    return Eigen::Vector3<T>(0.0, WGS84_OMEGA<T> * std::cos(phi), -WGS84_OMEGA<T> * std::sin(phi));
  }
}

/**
 * *=== TransportRate ===*
 * @brief Transport rate of the 'ECEF' frame relative to the local navigation frame
 * @tparam isNed Specify the local navigation frame as NED (true) or ENU (false)
 * @tparam Derived1 An Eigen size 3 object (i.e Vector3d, RowVector3d, Array3f)
 * @tparam Derived2 An Eigen size 3 object (i.e Vector3d, RowVector3d, Array3f)
 * @param lla    Latitude, Longitude, Height (rad, rad, m)
 * @param v_nb_e size 3 velocity vector in the local navigation coordinate system
 * @returns Transport rate in the local navigation frame
 */
template <bool isNed = true, typename Derived1, typename Derived2>
inline auto TransportRate(
    const Eigen::DenseBase<Derived1> &lla, const Eigen::DenseBase<Derived2> &v_nb_e) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, lla, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, v_nb_e, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  using Scalar = typename Derived1::Scalar;

  const Scalar phi = lla.derived()(0), h = lla.derived()(2);
  const auto [Re, Rn] = EarthRadii2(lla(0));
  if constexpr (isNed) {
    const Scalar vn = v_nb_e.derived()(0), ve = v_nb_e.derived()(1);
    Scalar ve_Reh = ve / (Re + h);
    return Eigen::Vector3<Scalar>(ve_Reh, -vn / (Rn + h), -ve_Reh * std::tan(phi));
  } else {
    const Scalar ve = v_nb_e.derived()(0), vn = v_nb_e.derived()(1);
    Scalar ve_Reh = ve / (Re + h);
    return Eigen::Vector3<Scalar>(-vn / (Rn + h), ve_Reh, ve_Reh * std::tan(phi));
  }
}

/**
 * *=== CoriolisRate ===
 * @brief Coriolis effect perceived in the local navigation frame
 * @tparam isNed Specify the local navigation frame as NED (true) or ENU (false)
 * @tparam Derived1 An Eigen size 3 object (i.e Vector3d, RowVector3d, Array3f)
 * @tparam Derived2 An Eigen size 3 object (i.e Vector3d, RowVector3d, Array3f)
 * @param lla Latitude, Longitude, Height (rad, rad, m)
 * @param v_nb_e size 3 velocity vector in the local navigation coordinate system
 * @returns Coriolis effect
 */
template <bool isNed = true, typename Derived1, typename Derived2>
inline auto CoriolisRate(
    const Eigen::DenseBase<Derived1> &lla, const Eigen::DenseBase<Derived2> &v_nb_e) {
  ASSERT_EIGEN_VEC_SIZE(Derived1, lla, 3);
  ASSERT_EIGEN_VEC_SIZE(Derived2, v_nb_e, 3);
  ASSERT_EIGEN_TYPE(Derived1, Derived2);
  using Scalar = typename Derived1::Scalar;

  Eigen::Vector3<Scalar> w_ie_n = EarthRate<isNed>(lla.derived()(0));
  Eigen::Vector3<Scalar> w_en_n = TransportRate<isNed>(lla, v_nb_e);
  return (w_en_n + 2.0 * w_ie_n).cross(v_nb_e);
}

/**
 * *=== LocalGravity ===
 * @brief Calculates gravity in the local navigation (ENU or NED) frame
 * @tparam isNed Specify the local navigation frame as NED (true) or ENU (false)
 * @tparam Derived An Eigen size 3 object (i.e Vector3d, RowVector3d, Array3f)
 * @param lla Latitude, Longitude, Height (rad, rad, m)
 * @returns Local navigation frame gravity
 */
template <bool isNed = true, typename Derived>
inline auto LocalGravity(const Eigen::DenseBase<Derived> &lla) {
  ASSERT_EIGEN_VEC_SIZE(Derived, lla, 3);
  using Scalar = typename Derived::Scalar;

  const Scalar phi = lla.derived()(0), h = lla.derived()(2);
  Scalar sin_phi2 = std::sin(phi);
  sin_phi2 *= sin_phi2;
  const Scalar g0 = WGS84_GRAVITY<Scalar> * ((1.0 + WGS84_SOMGLIANA<Scalar> * sin_phi2) /
                                             std::sqrt(1.0 - WGS84_E2<Scalar> * sin_phi2));
  const Scalar R02 = WGS84_A<Scalar> * WGS84_A<Scalar>;
  const Scalar OMEGA2 = WGS84_OMEGA<Scalar> * WGS84_OMEGA<Scalar>;
  const Scalar h2 = h * h;
  if (isNed) {
    return Eigen::Vector3<Scalar>(
        -8.08e-9 * h * std::sin(2.0 * phi),
        0.0,
        g0 * (1.0 - (1.0 + WGS84_F<Scalar> * (1.0 - 2.0 * sin_phi2) +
                     h * (OMEGA2 * R02 * WGS84_B<Scalar> / WGS84_GM<Scalar>)) *
                        (2.0 / WGS84_A<Scalar>)+(3.0 * h2 / R02)));
  } else {
    return Eigen::Vector3<Scalar>(
        0.0,
        -8.08e-9 * h * std::sin(2.0 * phi),
        -g0 * (1.0 - (1.0 + WGS84_F<Scalar> * (1.0 - 2.0 * sin_phi2) +
                      h * (OMEGA2 * R02 * WGS84_B<Scalar> / WGS84_GM<Scalar>)) *
                         (2.0 / WGS84_A<Scalar>)+(3.0 * h2 / R02)));
  }
}

/**
 * *=== EcefGravitation ===*
 * @brief Calculates gravitational acceleration in the Earth-Centered-Earth-Fixed frame
 * @tparam Derived An Eigen size 3 object (i.e Vector3d, RowVector3d, Array3f)
 * @param r_eb_e ECEF position (m,m,m)
 * @returns ECEF frame gravitational acceleration
 */
template <typename Derived>
inline auto EcefGravitation(const Eigen::DenseBase<Derived> &r_eb_e) {
  ASSERT_EIGEN_VEC_SIZE(Derived, r_eb_e, 3);
  using Scalar = typename Derived::Scalar;

  Scalar mag_r = r_eb_e.norm();
  if (mag_r == 0) return Eigen::Vector3<Scalar>::Zero();

  const Scalar x = r_eb_e.derived(0), y = r_eb_e.derived()(1), z = r_eb_e.derived()(2);
  Scalar c0 = 5.0 * std::pow(z / mag_r, 2);
  Scalar c1 = -WGS84_GM<Scalar> / std::pow(mag_r, 3);
  Scalar c2 = 1.5 * WGS84_J2<Scalar> * std::pow(WGS84_A<Scalar> / mag_r, 2);
  return c1 * (r_eb_e.derived().array() +
               c2 * Eigen::Array3<Scalar>((1.0 - c0) * x, (1.0 - c0) * y, (3.0 - c0) * z));
}

/**
 * *=== EcefGravity ===*
 * @brief Calculates gravity in the Earth-Centered-Earth-Fixed frame
 * @tparam Derived An Eigen size 3 object (i.e Vector3d, RowVector3d, Array3f)
 * @param r_eb_e ECEF position (m,m,m)
 * @returns ECEF frame gravity
 */
template <typename Derived>
inline auto EcefGravity(const Eigen::DenseBase<Derived> &r_eb_e) {
  using Scalar = typename Derived::Scalar;

  const Scalar w2 = WGS84_OMEGA<Scalar> * WGS84_OMEGA<Scalar>;
  const Scalar x = r_eb_e.derived(0), y = r_eb_e.derived()(1);
  return EcefGravitation(r_eb_e) + Eigen::Vector3<Scalar>(w2 * x, w2 * y, 0);
}

/**
 * *=== EciGravitation ===*
 * @brief Calculates gravitational acceleration in the Earth-Centered-Inertial frame
 * @tparam Derived An Eigen size 3 object (i.e Vector3d, RowVector3d, Array3f)
 * @param r_ib_i ECI position (m,m,m)
 * @returns ECI frame gravitational acceleration
 */
template <typename Derived>
inline auto EciGravitation(const Eigen::DenseBase<Derived> &r_ib_i) {
  ASSERT_EIGEN_VEC_SIZE(Derived, r_ib_i, 3);
  using Scalar = typename Derived::Scalar;

  const Scalar mag_r = r_ib_i.norm();
  if (mag_r == 0) return Eigen::Vector3<Scalar>::Zero();

  const Scalar x = r_ib_i.derived(0), y = r_ib_i.derived()(1), z = r_ib_i.derived()(2);
  Scalar c0 = 5.0 * std::pow(z / mag_r, 2);
  Scalar c1 = -WGS84_GM<Scalar> / std::pow(mag_r, 3);
  Scalar c2 = 1.5 * WGS84_J2<Scalar> * std::pow(WGS84_A<Scalar> / mag_r, 2);
  return c1 * (r_ib_i.derived().array() +
               c2 * Eigen::Array3<Scalar>((1.0 - c0) * x, (1.0 - c0) * y, (3.0 - c0) * z));
}

}  // namespace nt

#endif