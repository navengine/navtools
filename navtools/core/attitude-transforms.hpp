#ifndef NAVTOOLS_CORE_ATTITUDE_TRANSFORMS_HPP
#define NAVTOOLS_CORE_ATTITUDE_TRANSFORMS_HPP

#include "navtools/core/macros.hpp"
#include <Eigen/Dense>
#include <tuple>

namespace nt {

enum RotationAxis { X, Y, Z };
enum RotationOrder { XYZ, XZY, YXZ, YZX, ZXY, ZYX, XYX, XZX, YXY, YZY, ZXZ, ZYZ };

//! ============================================================================================ !//

/**
 * *=== angle2dcm ===*
 * @brief Builds the rotation matrix defined by "angle" and "axis"
 * @tparam axis Desired rotation axis (X, Y, or Z)
 * @tparam Derived The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param angle the rotation angle (rad)
 * @param R     the ouput rotation matrix
 */
template <RotationAxis axis, typename Derived>
inline void angle2dcm(const typename Derived::Scalar &angle, Eigen::DenseBase<Derived> &R) {
  ASSERT_EIGEN_MAT_SIZE(Derived, R, 3, 3);
  using Scalar = typename Derived::Scalar;

  Scalar cang = std::cos(angle);
  Scalar sang = std::sin(angle);

  if constexpr (axis == RotationAxis::X) {
    R.derived() << 1.0, 0.0, 0.0, 0.0, cang, sang, 0.0, -sang, cang;
  } else if constexpr (axis == RotationAxis::Y) {
    R.derived() << cang, 0.0, -sang, 0.0, 1.0, 0.0, sang, 0.0, cang;
  } else if constexpr (axis == RotationAxis::Z) {
    R.derived() << cang, sang, 0.0, -sang, cang, 0.0, 0.0, 0.0, 1.0;
  } else {
    throw std::runtime_error("Error: Invalid axis provided!");
  }
}

/**
 * *=== angle2dcm ===*
 * @brief Builds the rotation matrix defined by "angle" and "axis"
 * @tparam axis Desired rotation axis (X, Y, or Z)
 * @tparam Derived The type of the output Eigen dense matrix (e.g., Matrix3d...)
 * @param angle the rotation angle (rad)
 * @returns the ouput rotation matrix
 */
template <RotationAxis axis, typename Derived>
inline Derived angle2dcm(const typename Derived::Scalar &angle) {
  Derived R;
  angle2dcm<axis>(angle, R);
  return R;
}

/**
 * *=== angle2quat ===*
 * @brief Builds the quaternion defined by "angle" and "axis"
 * @tparam axis Desired rotation axis (X, Y, or Z)
 * @tparam Derived The type of the output Eigen dense quaternion (e.g., Vector4d, quat.coeffs())
 * @param angle the rotation angle (rad)
 * @param q     the ouput quaternion
 */
template <RotationAxis axis, typename Derived>
inline void angle2quat(const typename Derived::Scalar angle, Eigen::DenseBase<Derived> &q) {
  EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived, 4);
  using Scalar = typename Derived::Scalar;

  auto &quat = q.derived();
  Scalar sin_half = std::sin(angle / 2.0);
  Scalar cos_half = std::cos(angle / 2.0);
  if constexpr (axis == RotationAxis::X) {
    quat(0) = sin_half;
    quat(1) = 0.0;
    quat(2) = 0.0;
    quat(3) = cos_half;
  } else if constexpr (axis == RotationAxis::Y) {
    quat(0) = 0.0;
    quat(1) = sin_half;
    quat(2) = 0.0;
    quat(3) = cos_half;
  } else if constexpr (axis == RotationAxis::Z) {
    quat(0) = 0.0;
    quat(1) = 0.0;
    quat(2) = sin_half;
    quat(3) = cos_half;
  } else {
    throw std::runtime_error("Error: Invalid axis provided!");
  }
}

/**
 * *=== angle2quat ===*
 * @brief Builds the quaternion defined by "angle" and "axis"
 * @tparam axis Desired rotation axis (X, Y, or Z)
 * @tparam Derived The type of the output Eigen dense quaternion (e.g., Vector4d, quat.coeffs())
 * @param angle the rotation angle (rad)
 * @returns the ouput quaternion
 */
template <RotationAxis axis, typename Derived>
inline Derived angle2quat(const typename Derived::Scalar &angle) {
  Derived q;
  angle2quat<axis>(angle, q);
  return q;
}

//! ============================================================================================ !//

/**
 * *=== angles2dcm ===*
 * @brief Builds the rotation matrix defined by "roll", "pitch", and "yaw"
 * @tparam order Desired rotation order (e.g., XYZ...)
 * @tparam Derived The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param roll  roll angle (rad)
 * @param pitch pitch angle (rad)
 * @param yaw   yaw angle (rad)
 * @param R     the ouput rotation matrix
 */
template <RotationOrder order, typename Derived>
inline void angles2dcm(
    const typename Derived::Scalar &roll,
    const typename Derived::Scalar &pitch,
    const typename Derived::Scalar &yaw,
    Eigen::DenseBase<Derived> &R) {
  Derived R1, R2, R3;
  ASSERT_EIGEN_MAT_SIZE(Derived, R, 3, 3);
  ASSERT_EIGEN_MAT_SIZE(Derived, R1, 3, 3);
  ASSERT_EIGEN_MAT_SIZE(Derived, R2, 3, 3);
  ASSERT_EIGEN_MAT_SIZE(Derived, R3, 3, 3);

  if constexpr (order == RotationOrder::XYZ) {
    angle2dcm<RotationAxis::X>(roll, R1);
    angle2dcm<RotationAxis::Y>(pitch, R2);
    angle2dcm<RotationAxis::Z>(yaw, R3);
    R.derived() = R3 * R2 * R1;
  } else if constexpr (order == RotationOrder::XZY) {
    angle2dcm<RotationAxis::X>(roll, R1);
    angle2dcm<RotationAxis::Y>(pitch, R2);
    angle2dcm<RotationAxis::Z>(yaw, R3);
    R.derived() = R2 * R3 * R1;
  } else if constexpr (order == RotationOrder::YXZ) {
    angle2dcm<RotationAxis::X>(roll, R1);
    angle2dcm<RotationAxis::Y>(pitch, R2);
    angle2dcm<RotationAxis::Z>(yaw, R3);
    R.derived() = R3 * R1 * R2;
  } else if constexpr (order == RotationOrder::YZX) {
    angle2dcm<RotationAxis::X>(roll, R1);
    angle2dcm<RotationAxis::Y>(pitch, R2);
    angle2dcm<RotationAxis::Z>(yaw, R3);
    R.derived() = R1 * R3 * R2;
  } else if constexpr (order == RotationOrder::ZXY) {
    angle2dcm<RotationAxis::X>(roll, R1);
    angle2dcm<RotationAxis::Y>(pitch, R2);
    angle2dcm<RotationAxis::Z>(yaw, R3);
    R.derived() = R2 * R1 * R3;
  } else if constexpr (order == RotationOrder::ZYX) {
    angle2dcm<RotationAxis::X>(roll, R1);
    angle2dcm<RotationAxis::Y>(pitch, R2);
    angle2dcm<RotationAxis::Z>(yaw, R3);
    R.derived() = R1 * R2 * R3;
  } else if constexpr (order == RotationOrder::XYX) {
    angle2dcm<RotationAxis::X>(roll, R1);
    angle2dcm<RotationAxis::Y>(pitch, R2);
    angle2dcm<RotationAxis::X>(yaw, R3);
    R.derived() = R1 * R2 * R3;
  } else if constexpr (order == RotationOrder::XZX) {
    angle2dcm<RotationAxis::X>(roll, R1);
    angle2dcm<RotationAxis::Z>(pitch, R2);
    angle2dcm<RotationAxis::X>(yaw, R3);
    R.derived() = R1 * R2 * R3;
  } else if constexpr (order == RotationOrder::YXY) {
    angle2dcm<RotationAxis::Y>(roll, R1);
    angle2dcm<RotationAxis::X>(pitch, R2);
    angle2dcm<RotationAxis::Y>(yaw, R3);
    R.derived() = R1 * R2 * R3;
  } else if constexpr (order == RotationOrder::YZY) {
    angle2dcm<RotationAxis::Y>(roll, R1);
    angle2dcm<RotationAxis::Z>(pitch, R2);
    angle2dcm<RotationAxis::Y>(yaw, R3);
    R.derived() = R1 * R2 * R3;
  } else if constexpr (order == RotationOrder::ZXZ) {
    angle2dcm<RotationAxis::Z>(roll, R1);
    angle2dcm<RotationAxis::X>(pitch, R2);
    angle2dcm<RotationAxis::Z>(yaw, R3);
    R.derived() = R1 * R2 * R3;
  } else if constexpr (order == RotationOrder::ZYZ) {
    angle2dcm<RotationAxis::Z>(roll, R1);
    angle2dcm<RotationAxis::Y>(pitch, R2);
    angle2dcm<RotationAxis::Z>(yaw, R3);
    R.derived() = R1 * R2 * R3;
  } else {
    throw std::runtime_error("Error: Invalid order provided!");
  }
}

/**
 * *=== angles2dcm ===*
 * @brief Builds the rotation matrix defined by "roll", "pitch", and "yaw"
 * @tparam order Desired rotation order (e.g., XYZ...)
 * @tparam Derived The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param roll  roll angle (rad)
 * @param pitch pitch angle (rad)
 * @param yaw   yaw angle (rad)
 * @returns the ouput rotation matrix
 */
template <RotationOrder order, typename Derived>
inline Derived angles2dcm(
    const typename Derived::Scalar &roll,
    const typename Derived::Scalar &pitch,
    const typename Derived::Scalar &yaw) {
  Derived R;
  angles2dcm<order>(roll, pitch, yaw, R);
  return R;
}

/**
 * *=== angles2euler ===*
 * @brief Converts roll, pitch, and yaw angles to euler angle vector
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param roll  roll angle (rad)
 * @param pitch pitch angle (rad)
 * @param yaw   yaw angle (rad)
 * @param e     euler angle vector
 */
template <typename Derived>
inline void angles2euler(
    const typename Derived::Scalar &roll,
    const typename Derived::Scalar &pitch,
    const typename Derived::Scalar &yaw,
    Eigen::DenseBase<Derived> &e) {
  ASSERT_EIGEN_VEC_SIZE(Derived, e, 3);
  e.derived() << roll, pitch, yaw;
}

/**
 * *=== angles2euler ===*
 * @brief Converts roll, pitch, and yaw angles to euler angle vector
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param roll  roll angle (rad)
 * @param pitch pitch angle (rad)
 * @param yaw   yaw angle (rad)
 * @param e     euler angle vector
 */
template <typename Derived>
inline Derived angles2euler(
    const typename Derived::Scalar &roll,
    const typename Derived::Scalar &pitch,
    const typename Derived::Scalar &yaw) {
  Derived e;
  angles2euler(roll, pitch, yaw, e);
  return e;
}

/**
 * *=== angles2quat ===*
 * @brief Converts roll, pitch, and yaw angles to quaternion
 * @tparam order Desired rotation order (e.g., XYZ...)
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector4d, Array4f)
 * @param roll  roll angle (rad)
 * @param pitch pitch angle (rad)
 * @param yaw   yaw angle (rad)
 * @param q     the ouput coefficients for the quaternion (e.g. quat.coeffs())
 */
template <RotationOrder order, typename Derived>
inline void angles2quat(
    const typename Derived::Scalar &roll,
    const typename Derived::Scalar &pitch,
    const typename Derived::Scalar &yaw,
    Eigen::DenseBase<Derived> &q) {
  ASSERT_EIGEN_VEC_SIZE(Derived, q, 4);
  using Scalar = typename Derived::Scalar;

  Scalar s1 = std::sin(0.5 * roll), s2 = std::sin(0.5 * pitch), s3 = std::sin(0.5 * yaw);
  Scalar c1 = std::cos(0.5 * roll), c2 = std::cos(0.5 * pitch), c3 = std::cos(0.5 * yaw);
  auto &quat = q.derived();  // order is [x, y, z, w]

  if constexpr (order == RotationOrder::XYZ) {
    quat(3) = c1 * c2 * c3 - s1 * s2 * s3;
    quat(0) = c1 * s2 * s3 + s1 * c2 * c3;
    quat(1) = c1 * s2 * c3 - s1 * c2 * s3;
    quat(2) = c1 * c2 * s3 + s1 * s2 * c3;
  } else if constexpr (order == RotationOrder::XZY) {
    quat(3) = c1 * c2 * c3 + s1 * s2 * s3;
    quat(0) = s1 * c2 * c3 - c1 * s2 * s3;
    quat(1) = c1 * s2 * c3 - s1 * c2 * s3;
    quat(2) = c1 * c2 * s3 + s1 * s2 * c3;
  } else if constexpr (order == RotationOrder::YXZ) {
    quat(3) = c1 * c2 * c3 + s1 * s2 * s3;
    quat(0) = s1 * c2 * c3 + c1 * s2 * s3;
    quat(1) = c1 * s2 * c3 - s1 * c2 * s3;
    quat(2) = c1 * c2 * s3 - s1 * s2 * c3;
  } else if constexpr (order == RotationOrder::YZX) {
    quat(3) = c1 * c2 * c3 - s1 * s2 * s3;
    quat(0) = s1 * c2 * c3 + c1 * s2 * s3;
    quat(1) = s1 * c2 * s3 + c1 * s2 * c3;
    quat(2) = c1 * c2 * s3 - s1 * s2 * c3;
  } else if constexpr (order == RotationOrder::ZXY) {
    quat(3) = c1 * c2 * c3 - s1 * s2 * s3;
    quat(0) = s1 * c2 * c3 - c1 * s2 * s3;
    quat(1) = c1 * s2 * c3 + s1 * c2 * s3;
    quat(2) = s1 * s2 * c3 + c1 * c2 * s3;
  } else if constexpr (order == RotationOrder::ZYX) {
    quat(3) = c1 * c2 * c3 + s1 * s2 * s3;
    quat(0) = s1 * c2 * c3 - c1 * s2 * s3;
    quat(1) = c1 * s2 * c3 + s1 * c2 * s3;
    quat(2) = c1 * c2 * s3 - s1 * s2 * c3;
  } else if constexpr (order == RotationOrder::XYX) {
    quat(3) = c1 * c2 * c3 - s1 * c2 * s3;
    quat(0) = c1 * c2 * s3 + s1 * c2 * c3;
    quat(1) = c1 * s2 * c3 + s1 * s2 * s3;
    quat(2) = c1 * s2 * s3 - s1 * s2 * c3;
  } else if constexpr (order == RotationOrder::XZX) {
    quat(3) = c1 * c2 * c3 - s1 * c2 * s3;
    quat(0) = c1 * c2 * s3 + s1 * c2 * c3;
    quat(1) = s1 * s2 * c3 - c1 * s2 * s3;
    quat(2) = c1 * s2 * c3 + s1 * s2 * s3;
  } else if constexpr (order == RotationOrder::YXY) {
    quat(3) = c1 * c2 * c3 - s1 * c2 * s3;
    quat(0) = c1 * s2 * c3 + s1 * s2 * s3;
    quat(1) = c1 * c2 * s3 + s1 * c2 * c3;
    quat(2) = s1 * s2 * c3 - c1 * s2 * s3;
  } else if constexpr (order == RotationOrder::YZY) {
    quat(3) = c1 * c2 * c3 - s1 * c2 * s3;
    quat(0) = c1 * s2 * s3 - s1 * s2 * c3;
    quat(1) = s1 * c2 * c3 + c1 * c2 * s3;
    quat(2) = c1 * s2 * c3 + s1 * s2 * s3;
  } else if constexpr (order == RotationOrder::ZXZ) {
    quat(3) = c1 * c2 * c3 - s1 * c2 * s3;
    quat(0) = c1 * s2 * c3 + s1 * s2 * s3;
    quat(1) = c1 * s2 * s3 - s1 * s2 * c3;
    quat(2) = s1 * c2 * c3 + c1 * c2 * s3;
  } else if constexpr (order == RotationOrder::ZYZ) {
    quat(3) = c1 * c2 * c3 - s1 * c2 * s3;
    quat(0) = s1 * s2 * c3 - c1 * s2 * s3;
    quat(1) = c1 * s2 * c3 + s1 * s2 * s3;
    quat(2) = c1 * c2 * s3 + s1 * c2 * c3;
  } else {
    throw std::runtime_error("Error: Invalid order provided!");
  }
}

/**
 * *=== angles2quat ===*
 * @brief Converts roll, pitch, and yaw angles to quaternion
 * @tparam order Desired rotation order (e.g., XYZ...)
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector4d, Array4f)
 * @param roll  roll angle (rad)
 * @param pitch pitch angle (rad)
 * @param yaw   yaw angle (rad)
 * @returns the ouput coefficients for the quaternion (e.g. quat.coeffs())
 */
template <RotationOrder order, typename Derived>
inline Derived angles2quat(
    const typename Derived::Scalar &roll,
    const typename Derived::Scalar &pitch,
    const typename Derived::Scalar &yaw) {
  Derived q;
  ASSERT_EIGEN_VEC_SIZE(Derived, q, 4);
  angles2quat<order>(roll, pitch, yaw, q);
  return q;
}

//! ============================================================================================ !//

/**
 * *=== euler2angles ===*
 * @brief Converts euler angle vector to roll, pitch, and yaw angles
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param e     euler angle vector
 * @param roll  roll angle (rad)
 * @param pitch pitch angle (rad)
 * @param yaw   yaw angle (rad)
 */
template <typename Derived>
inline void euler2angles(
    const Eigen::DenseBase<Derived> &e,
    typename Derived::Scalar &roll,
    typename Derived::Scalar &pitch,
    typename Derived::Scalar &yaw) {
  ASSERT_EIGEN_VEC_SIZE(Derived, e, 3);

  const auto &eul = e.derived();
  roll = eul(0);
  pitch = eul(1);
  yaw = eul(2);
}

/**
 * *=== euler2angles ===*
 * @brief Converts euler angle vector to roll, pitch, and yaw angles
 * @tparam Derived The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param e euler angle vector
 * @returns a tuple containing the roll, pitch, and yaw angles (rad)
 */
template <typename Derived>
inline std::tuple<typename Derived::Scalar> euler2angles(const Eigen::DenseBase<Derived> &e) {
  ASSERT_EIGEN_VEC_SIZE(Derived, e, 3);

  const auto &eul = e.derived();
  return std::make_tuple(eul(0), eul(1), eul(2));
}

/**
 * *=== euler2dcm ===*
 * @brief Builds the rotation matrix defined by "roll", "pitch", and "yaw" euler angle vector
 * @tparam order Desired rotation order (e.g., XYZ...)
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @tparam DerivedVec The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param e euler angles vector (assumes "roll", "pitch", "yaw" ordered vector) (rad)
 * @returns the ouput rotation matrix
 */
template <RotationOrder order, typename DerivedMat, typename DerivedVec>
inline DerivedMat euler2dcm(const Eigen::DenseBase<DerivedVec> &e) {
  ASSERT_EIGEN_VEC_SIZE(DerivedVec, e, 3);

  const auto &eul = e.derived();
  return angles2dcm<order, DerivedMat>(eul(0), eul(1), eul(2));
}

/**
 * *=== euler2dcm ===*
 * @brief Builds the rotation matrix defined by "roll", "pitch", and "yaw"
 * @tparam order Desired rotation order (e.g., XYZ...)
 * @tparam DerivedVec The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param e euler angles vector (assumes "roll", "pitch", "yaw" ordered vector) (rad)
 * @param R the ouput rotation matrix
 */
template <RotationOrder order, typename DerivedVec, typename DerivedMat>
inline void euler2dcm(const Eigen::DenseBase<DerivedVec> &e, Eigen::DenseBase<DerivedMat> &R) {
  ASSERT_EIGEN_VEC_SIZE(DerivedVec, e, 3);
  ASSERT_EIGEN_SAME_SCALAR(DerivedMat, DerivedVec);

  const auto &eul = e.derived();
  angles2dcm<order>(eul(0), eul(1), eul(2), R);
}

/**
 * *=== euler2quat ===*
 * @brief Builds the quaternion defined by "roll", "pitch", and "yaw" euler angle vector
 * @tparam order Desired rotation order (e.g., XYZ...)
 * @tparam DerivedQuat The type of the input Eigen dense quaternion (e.g., Vector4d, quat.coeffs())
 * @tparam DerivedVec The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @param e euler angles vector (assumes "roll", "pitch", "yaw" ordered vector) (rad)
 * @returns the ouput quaternion
 */
template <RotationOrder order, typename DerivedQuat, typename DerivedVec>
inline DerivedQuat euler2quat(const Eigen::DenseBase<DerivedVec> &e) {
  ASSERT_EIGEN_VEC_SIZE(DerivedVec, e, 3);

  const auto &eul = e.derived();
  return angles2quat<order, DerivedQuat>(eul(0), eul(1), eul(2));
}

/**
 * *=== euler2quat ===*
 * @brief Builds the quaternion defined by "roll", "pitch", and "yaw"
 * @tparam axis Desired rotation axis (e.g., XYZ...)
 * @tparam DerivedVec The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @tparam DerivedQuat The type of the input Eigen dense quaternion (e.g., Vector4d, quat.coeffs())
 * @param e euler angles vector (assumes "roll", "pitch", "yaw" ordered vector) (rad)
 * @param q the ouput quaternion
 */
template <RotationOrder order, typename DerivedVec, typename DerivedQuat>
inline void euler2quat(const Eigen::DenseBase<DerivedVec> &e, Eigen::DenseBase<DerivedQuat> &q) {
  ASSERT_EIGEN_VEC_SIZE(DerivedVec, e, 3);
  ASSERT_EIGEN_VEC_SIZE(DerivedQuat, q, 4);
  ASSERT_EIGEN_SAME_SCALAR(DerivedVec, DerivedQuat);

  auto &eul = e.derived();
  angles2quat<order>(eul(0), eul(1), eul(2), q);
}

//! ============================================================================================ !//

/**
 * *=== dcm2angles ===*
 * @brief Extracts "roll", "pitch", and "yaw" from the defined rotation matrix
 * @tparam order Desired rotation order (e.g., XYZ...)
 * @tparam Derived The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param R     the input rotation matrix
 * @param roll  roll angle (rad)
 * @param pitch pitch angle (rad)
 * @param yaw   yaw angle (rad)
 */
template <RotationOrder order, typename Derived>
inline void dcm2angles(
    const Eigen::DenseBase<Derived> &R,
    typename Derived::Scalar &roll,
    typename Derived::Scalar &pitch,
    typename Derived::Scalar &yaw) {
  ASSERT_EIGEN_MAT_SIZE(Derived, R, 3, 3);

  const auto &C = R.derived();
  if constexpr (order == RotationOrder::XYZ) {
    roll = std::atan2(-C(2, 1), C(2, 2));
    pitch = std::asin(C(2, 0));
    yaw = std::atan2(-C(1, 0), C(0, 0));
  } else if constexpr (order == RotationOrder::XZY) {
    roll = std::atan2(C(1, 2), C(1, 1));
    pitch = std::atan2(C(2, 0), C(0, 0));
    yaw = std::asin(-C(1, 0));
  } else if constexpr (order == RotationOrder::YXZ) {
    roll = std::asin(-C(2, 1));
    pitch = std::atan2(C(2, 0), C(2, 2));
    yaw = std::atan2(C(0, 1), C(1, 1));
  } else if constexpr (order == RotationOrder::YZX) {
    roll = std::atan2(-C(2, 1), C(1, 1));
    pitch = std::atan2(-C(0, 2), C(0, 0));
    yaw = std::asin(C(0, 1));
  } else if constexpr (order == RotationOrder::ZXY) {
    roll = std::asin(C(1, 2));
    pitch = std::atan2(-C(0, 2), C(2, 2));
    yaw = std::atan2(-C(1, 0), C(1, 1));
  } else if constexpr (order == RotationOrder::ZYX) {
    roll = std::atan2(C(1, 2), C(2, 2));
    pitch = std::asin(-C(0, 2));
    yaw = std::atan2(C(0, 1), C(0, 0));
  } else if constexpr (order == RotationOrder::XYX) {
    roll = std::atan2(C(1, 0), C(2, 0));
    pitch = std::acos(C(0, 0));
    yaw = std::atan2(C(0, 1), -C(0, 2));
  } else if constexpr (order == RotationOrder::XZX) {
    roll = std::atan2(C(2, 0), -C(1, 0));
    pitch = std::acos(C(0, 0));
    yaw = std::atan2(C(0, 2), C(0, 1));
  } else if constexpr (order == RotationOrder::YXY) {
    roll = std::atan2(C(0, 1), -C(2, 1));
    pitch = std::acos(C(1, 1));
    yaw = std::atan2(C(1, 0), C(1, 2));
  } else if constexpr (order == RotationOrder::YZY) {
    roll = std::atan2(C(2, 1), C(0, 1));
    pitch = std::acos(C(1, 1));
    yaw = std::atan2(C(1, 2), -C(1, 0));
  } else if constexpr (order == RotationOrder::ZXZ) {
    roll = std::atan2(C(0, 2), C(1, 2));
    pitch = std::acos(C(2, 2));
    yaw = std::atan2(C(2, 0), -C(2, 1));
  } else if constexpr (order == RotationOrder::ZYZ) {
    roll = std::atan2(C(1, 2), -C(0, 2));
    pitch = std::acos(C(2, 2));
    yaw = std::atan2(C(2, 1), C(2, 0));
  } else {
    throw std::runtime_error("Error: Invalid order provided!");
  }
}

/**
 * *=== dcm2angles ===*
 * @brief Extracts "roll", "pitch", and "yaw" from the defined rotation matrix
 * @tparam order Desired rotation order (e.g., XYZ...)
 * @tparam Derived The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param R     the input rotation matrix
 * @param roll  roll angle (rad)
 * @param pitch pitch angle (rad)
 * @param yaw   yaw angle (rad)
 */
template <RotationOrder order, typename Derived>
inline std::tuple<typename Derived::Scalar> dcm2angles(const Eigen::DenseBase<Derived> &R) {
  using Scalar = typename Derived::Scalar;
  Scalar roll, pitch, yaw;
  dcm2angles<order>(R, roll, pitch, yaw);
  return std::make_tuple(roll, pitch, yaw);
}

/**
 * *=== dcm2euler ===*
 * @brief Extracts euler angle vector from the defined rotation matrix
 * @tparam order Desired rotation order (e.g., XYZ...)
 * @tparam DerivedVec The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param R the input rotation matrix
 * @param e the output euler angle vector
 */
template <RotationOrder order, typename DerivedVec, typename DerivedMat>
inline void dcm2euler(const Eigen::DenseBase<DerivedMat> &R, Eigen::DenseBase<DerivedVec> &e) {
  ASSERT_EIGEN_VEC_SIZE(DerivedVec, e, 3);
  ASSERT_EIGEN_SAME_SCALAR(DerivedMat, DerivedVec);
  auto &eul = e.derived();
  dcm2angles<order>(R, eul(0), eul(1), eul(2));
}

/**
 * *=== dcm2euler ===*
 * @brief Extracts euler angle vector from the defined rotation matrix
 * @tparam order Desired rotation order (e.g., XYZ...)
 * @tparam DerivedVec The type of the input Eigen dense vector (e.g., Vector3d, Array3f)
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param R the input rotation matrix
 * @param e the output euler angle vector
 */
template <RotationOrder order, typename DerivedVec, typename DerivedMat>
inline DerivedVec dcm2euler(const Eigen::DenseBase<DerivedMat> &R) {
  DerivedVec e;
  dcm2euler<order>(R, e);
  return e;
}

/**
 * *=== dcm2quat ===*
 * @brief Converts a rotation matrix to its quaternion form
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @tparam DerivedQuat The type of the input Eigen dense quaternion (e.g., Vector4d, quat.coeffs())
 * @param R the input rotation matrix
 * @param q the ouput quaternion
 */
template <typename DerivedMat, typename DerivedQuat>
inline void dcm2quat(const Eigen::DenseBase<DerivedMat> &R, Eigen::DenseBase<DerivedQuat> &q) {
  ASSERT_EIGEN_MAT_SIZE(DerivedMat, R, 3, 3);
  ASSERT_EIGEN_VEC_SIZE(DerivedQuat, q, 4);
  ASSERT_EIGEN_SAME_SCALAR(DerivedMat, DerivedQuat);
  using Scalar = typename DerivedMat::Scalar;

  const auto C = R.derived();
  auto &quat = q.derived();
  Scalar trace = C(0, 0) + C(1, 1) + C(2, 2);
  if (trace > 0.0) {
    Scalar S = 2 * std::sqrt(1.0 + trace);
    quat(0) = (C(1, 2) - C(2, 1)) / S;
    quat(1) = (C(2, 0) - C(0, 2)) / S;
    quat(2) = (C(0, 1) - C(1, 0)) / S;
    quat(3) = 0.25 * S;
  } else if ((C(0, 0) > C(1, 1)) && (C(0, 0) > C(2, 2))) {
    Scalar S = 2 * std::sqrt(1.0 + C(0, 0) - C(1, 1) - C(2, 2));
    quat(0) = 0.25 * S;
    quat(1) = (C(0, 1) + C(1, 0)) / S;
    quat(2) = (C(2, 0) + C(0, 2)) / S;
    quat(2) = (C(1, 2) - C(2, 1)) / S;
  } else if (C(1, 1) > C(2, 2)) {
    Scalar S = 2 * std::sqrt(1.0 - C(0, 0) + C(1, 1) - C(2, 2));
    quat(0) = (C(1, 0) + C(0, 1)) / S;
    quat(1) = 0.25 * S;
    quat(2) = (C(1, 2) + C(2, 1)) / S;
    quat(3) = (C(2, 0) - C(0, 2)) / S;
  } else {
    Scalar S = 2.0 * std::sqrt(1.0 - C(0, 0) - C(1, 1) + C(2, 2));
    quat(0) = (C(2, 0) + C(0, 2)) / S;
    quat(1) = (C(1, 2) + C(2, 1)) / S;
    quat(2) = 0.25 * S;
    quat(3) = (C(0, 1) - C(1, 0)) / S;
  }
}

/**
 * *=== dcm2quat ===*
 * @brief Converts a rotation matrix to its quaternion form
 * @tparam DerivedQuat The type of the input Eigen dense quaternion (e.g., Vector4d, quat.coeffs())
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param R the input rotation matrix
 * @returns the ouput quaternion
 */
template <typename DerivedQuat, typename DerivedMat>
inline DerivedQuat dcm2quat(const Eigen::DenseBase<DerivedMat> &R) {
  DerivedQuat q;
  dcm2quat(R, q);
  return q;
}

//! ============================================================================================ !//

/**
 * *=== quat2angles ===*
 * @brief Covert quaternion into individual euler angles "roll", "pitch", and "yaw"
 * @tparam order Desired rotation order (e.g., XYZ...)
 * @tparam Derived The type of the input Eigen dense quaternion (e.g., Vector4d, quat.coeffs())
 * @param q     the input quaternion
 * @param roll  roll angle (rad)
 * @param pitch pitch angle (rad)
 * @param yaw   yaw angle (rad)
 */
template <RotationOrder order, typename Derived>
inline void quat2angles(
    const Eigen::DenseBase<Derived> &q,
    typename Derived::Scalar &roll,
    typename Derived::Scalar &pitch,
    typename Derived::Scalar &yaw) {
  ASSERT_EIGEN_VEC_SIZE(Derived, q, 4);
  using Scalar = typename Derived::Scalar;

  const auto &quat = q.derived();
  Scalar qx = quat(0), qy = quat(1), qz = quat(2), qw = quat(3);
  Scalar qx2 = qx * qx;
  Scalar qy2 = qy * qy;
  Scalar qz2 = qz * qz;
  Scalar qw2 = qw * qw;
  if constexpr (order == RotationOrder::XYZ) {
    roll = std::atan2(2 * (qw * qx - qy * qz), (qw2 - qx2 - qy2 + qz2));
    pitch = std::asin(2 * (qw * qy + qx * qz));
    yaw = std::atan2(2 * (qw * qz - qx * qy), (qw2 + qx2 - qy2 - qz2));
  } else if constexpr (order == RotationOrder::XZY) {
    roll = std::atan2(2 * (qy * qz + qw * qx), (qw2 - qx2 + qy2 - qz2));
    pitch = std::atan2(2 * (qx * qz + qw * qy), (qw2 + qx2 - qy2 - qz2));
    yaw = std::asin(2 * (qw * qz - qx * qy));
  } else if constexpr (order == RotationOrder::YXZ) {
    roll = std::asin(2 * (qw * qx - qy * qz));
    pitch = std::atan2(2 * (qx * qz + qw * qy), (qw2 - qx2 - qy2 + qz2));
    yaw = std::atan2(2 * (qx * qy + qw * qz), (qw2 - qx2 + qy2 - qz2));
  } else if constexpr (order == RotationOrder::YZX) {
    roll = std::atan2(2 * (qw * qx - qy * qz), (qw2 - qx2 + qy2 - qz2));
    pitch = std::atan2(2 * (qw * qy - qx * qz), (qw2 + qx2 - qy2 - qz2));
    yaw = std::asin(2 * (qx * qy + qw * qz));
  } else if constexpr (order == RotationOrder::ZXY) {
    roll = std::asin(2 * (qy * qz + qw * qx));
    pitch = std::atan2(2 * (qw * qy - qx * qz), (qw2 - qx2 - qy2 + qz2));
    yaw = std::atan2(2 * (qw * qz - qx * qy), (qw2 - qx2 + qy2 - qz2));
  } else if constexpr (order == RotationOrder::ZYX) {
    roll = std::atan2(2 * (qy * qz + qw * qx), (qw2 - qx2 - qy2 + qz2));
    pitch = std::asin(2 * (qw * qy - qx * qz));
    yaw = std::atan2(2 * (qx * qy + qw * qz), (qw2 + qx2 - qy2 - qz2));
  } else if constexpr (order == RotationOrder::XYX) {
    roll = std::atan2(-2 * (qw * qz - qx * qy), 2 * (qx * qz + qw * qy));
    pitch = std::acos(qw2 + qx2 - qy2 - qz2);
    yaw = std::atan2(2 * (qx * qy + qw * qz), -2 * (qx * qz - qw * qy));
  } else if constexpr (order == RotationOrder::XZX) {
    roll = std::atan2(2 * (qw * qy + qx * qz), -2 * (qx * qy - qw * qz));
    pitch = std::acos(qw2 + qx2 - qy2 - qz2);
    yaw = std::atan2(2 * (qx * qz - qw * qy), 2 * (qx * qy + qw * qz));
  } else if constexpr (order == RotationOrder::YXY) {
    roll = std::atan2(2 * (qw * qz + qx * qy), -2 * (qy * qz - qw * qx));
    pitch = std::acos(qw2 - qx2 + qy2 - qz2);
    yaw = std::atan2(-2 * (qw * qz - qx * qy), 2 * (qy * qz + qw * qx));
  } else if constexpr (order == RotationOrder::YZY) {
    roll = std::atan2(-2 * (qw * qx - qy * qz), 2 * (qy * qx + qw * qz));
    pitch = std::acos(qw2 - qx2 + qy2 - qz2);
    yaw = std::atan2(2 * (qw * qx + qy * qz), -2 * (qy * qx - qw * qz));
  } else if constexpr (order == RotationOrder::ZXZ) {
    roll = std::atan2(-2 * (qw * qy - qx * qz), 2 * (qy * qz + qw * qx));
    pitch = std::acos(qw2 - qx2 - qy2 + qz2);
    yaw = std::atan2(2 * (qw * qy + qx * qz), -2 * (qy * qz - qw * qx));
  } else if constexpr (order == RotationOrder::ZYZ) {
    roll = std::atan2(2 * (qw * qx + qy * qz), -2 * (qx * qz - qw * qy));
    pitch = std::acos(qw2 - qx2 - qy2 + qz2);
    yaw = std::atan2(-2 * (qw * qx - qy * qz), 2 * (qx * qz + qw * qy));
  } else {
    throw std::runtime_error("Error: Invalid rotation order provided!");
  }
}

/**
 * *=== quat2angles ===*
 * @brief Covert quaternion into individual euler angles "roll", "pitch", and "yaw"
 * @tparam order Desired rotation order (e.g., XYZ...)
 * @tparam Derived The type of the input Eigen dense quaternion (e.g., Vector4d, quat.coeffs())
 * @param q the input quaternion
 * @returns tuple of individual euler angles
 */
template <RotationOrder order, typename Derived>
inline std::tuple<typename Derived::Scalar> quat2angles(const Eigen::DenseBase<Derived> &q) {
  using Scalar = typename Derived::Scalar;
  Scalar roll, pitch, yaw;
  quat2angles<order>(q, roll, pitch, yaw);
  return std::make_tuple(roll, pitch, yaw);
}

/**
 * *=== quat2euler ===*
 * @brief Covert quaternion into the euler angles vector of "roll", "pitch", and "yaw"
 * @tparam order Desired rotation order (e.g., XYZ...)
 * @tparam DerivedQuat The type of the input Eigen dense quaternion (e.g., Vector4d, quat.coeffs())
 * @param q the input quaternion
 * @param e the output euler angle vector
 */
template <RotationOrder order, typename DerivedQuat, typename DerivedVec>
inline void quat2euler(const Eigen::DenseBase<DerivedQuat> &q, Eigen::DenseBase<DerivedVec> &e) {
  ASSERT_EIGEN_VEC_SIZE(DerivedVec, e, 3);
  ASSERT_EIGEN_SAME_SCALAR(DerivedVec, DerivedQuat);

  auto &eul = e.derived();
  quat2angles<order>(q, eul(0), eul(1), eul(2));
}

/**
 * *=== quat2euler ===*
 * @brief Covert quaternion into the euler angles vector of "roll", "pitch", and "yaw"
 * @tparam order Desired rotation order (e.g., XYZ...)
 * @tparam DerivedQuat The type of the input Eigen dense quaternion (e.g., Vector4d, quat.coeffs())
 * @param q the input quaternion
 * @returns the output euler angle vector
 */
template <RotationOrder order, typename DerivedVec, typename DerivedQuat>
inline DerivedVec quat2euler(const Eigen::DenseBase<DerivedQuat> &q) {
  DerivedVec e;
  quat2euler<order>(q, e);
  return e;
}

/**
 * *=== quat2dcm ===*
 * @brief Converts a quaternion to its rotation matrix form
 * @tparam DerivedQuat The type of the input Eigen dense quaternion (e.g., Vector4d, quat.coeffs())
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param q the input quaternion
 * @param R           the ouput rotation matrix
 */
template <typename DerivedQuat, typename DerivedMat>
inline void quat2dcm(const Eigen::DenseBase<DerivedQuat> &q, Eigen::DenseBase<DerivedMat> &R) {
  ASSERT_EIGEN_VEC_SIZE(DerivedQuat, q, 4);
  ASSERT_EIGEN_MAT_SIZE(DerivedMat, R, 3, 3);
  ASSERT_EIGEN_SAME_SCALAR(DerivedQuat, DerivedMat);
  using Scalar = typename DerivedMat::Scalar;

  const auto &quat = q.derived() / q.derived().norm();
  Scalar qx = quat(0), qy = quat(1), qz = quat(2), qw = quat(3);
  Scalar w2 = qw * qw;
  Scalar x2 = qx * qx;
  Scalar y2 = qy * qy;
  Scalar z2 = qz * qz;
  Scalar wx = qw * qx;
  Scalar wy = qw * qy;
  Scalar wz = qw * qz;
  Scalar xy = qx * qy;
  Scalar xz = qx * qz;
  Scalar yz = qy * qz;
  // clang-format off
  R.derived() << w2 + x2 - y2 - z2,   2.0 * (xy + wz),   2.0 * (xz - wy),
                   2.0 * (xy - wz), w2 - x2 + y2 - z2,   2.0 * (yz + wx),
                   2.0 * (xz + wy),   2.0 * (yz - wx), w2 - x2 - y2 + z2;
  // clang-format on
}

/**
 * *=== quat2dcm ===*
 * @brief Converts a quaternion to its rotation matrix form
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @tparam DerivedQuat The type of the input Eigen dense quaternion (e.g., Vector4d, quat.coeffs())
 * @param q the input quaternion
 * @returns the ouput rotation matrix
 */
template <typename DerivedMat, typename DerivedQuat>
inline DerivedMat quat2dcm(const Eigen::DenseBase<DerivedQuat> &q) {
  DerivedMat R;
  quat2dcm(q, R);
  return R;
}

}  // namespace nt

#endif
