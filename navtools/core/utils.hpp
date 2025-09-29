#ifndef NAVTOOLS_CORE_UTILS_HPP
#define NAVTOOLS_CORE_UTILS_HPP

#include "navtools/core/macros.hpp"
#include <Eigen/Dense>

namespace nt {

/**
 * *=== CircMod ===*
 * @brief Modulus of floating point number
 * @tparam T A floating point type (e.g. double, float)
 * @param x user input and output
 * @param y value to take modulus about
 * @returns modulus of number
 */
template <typename T>
inline constexpr T CircMod(T x, const T y) {
  return x - std::floor(x / y) * y;
}

/**
 * *=== Skew ===*
 * @brief Creates a 3x3 skew-symmetric matrix from a 3x1 vector
 * @tparam DerivedVec The type of the input Eigen dense vector (e.g., Vector3d, Array3ff)
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param v The 3-component input vector. It can be a column vector, row vector, etc
 * @param M 3x3 skew-symmetric matrix of the same scalar type as the input
 */
template <typename DerivedVec, typename DerivedMat>
inline void Skew(const Eigen::DenseBase<DerivedVec>& v, Eigen::DenseBase<DerivedMat>& M) {
  ASSERT_EIGEN_VEC_SIZE(DerivedVec, v, 3);
  ASSERT_EIGEN_MAT_SIZE(DerivedMat, M, 3, 3);
  ASSERT_EIGEN_SCALAR_TYPE(DerivedVec, DerivedMat);

  const auto& vec = v.derived();
  // clang-format off
  M.derived() <<    0.0, -vec(2), vec(1), 
                 vec(2),    0.0, -vec(0), 
                -vec(1), vec(0),     0.0;
  // clang-format on
}

/**
 * *=== Skew ===*
 * @brief Creates a 3x3 skew-symmetric matrix from a 3x1 vector
 * @tparam DerivedVec The type of the input Eigen dense vector (e.g., Vector3d, Array3ff)
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param v The 3-component input vector. It can be a column vector, row vector, etc
 * @return A 3x3 skew-symmetric matrix of the same scalar type as the input
 */
template <typename DerivedMat, typename DerivedVec>
inline DerivedMat Skew(const Eigen::DenseBase<DerivedVec>& v) {
  DerivedMat M;
  Skew(v, M);
  return M;
}

/**
 * *=== Deskew ===*
 * @brief Converts a 3x3 skew-symmetric matrix or array to a 3-component vector or array
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @tparam DerivedVec The type of the input Eigen dense vector (e.g., Vector3d, Array3ff)
 * @param M The 3x3 skew-symmetric input matrix
 * @param v The 3-component vector/array corresponding to the input
 */
template <typename DerivedMat, typename DerivedVec>
inline void Deskew(const Eigen::DenseBase<DerivedMat>& M, Eigen::DenseBase<DerivedVec>& v) {
  ASSERT_EIGEN_VEC_SIZE(DerivedVec, v, 3);
  ASSERT_EIGEN_MAT_SIZE(DerivedMat, M, 3, 3);
  ASSERT_EIGEN_SCALAR_TYPE(DerivedVec, DerivedMat);

  const auto& Mat = M.derived();
  v << (Mat(2, 1) - Mat(1, 2)) / 2.0, (Mat(0, 2) - Mat(2, 0)) / 2.0, (Mat(1, 0) - Mat(0, 1)) / 2.0;
}

/**
 * *=== Deskew ===*
 * @brief Converts a 3x3 skew-symmetric matrix or array to a 3-component vector or array
 * @tparam DerivedVec The type of the input Eigen dense vector (e.g., Vector3d, Array3ff)
 * @tparam DerivedMat The type of the input Eigen dense matrix (e.g., Matrix3d...)
 * @param M The 3x3 skew-symmetric input matrix
 * @return The 3-component vector/array corresponding to the input
 */
template <typename DerivedVec, typename DerivedMat>
inline DerivedVec Deskew(const Eigen::DenseBase<DerivedMat>& M) {
  DerivedVec v;
  Deskew(M, v);
  return v;
}

/**
 * *=== Rodrigues ===*
 * @brief Rodrigues formula for the approximation of a matrix exponential
 * @tparam Derived An Eigen size 3 object (i.e Vector3d, RowVector3d, Array3f)
 * @param vec      size 3 vector
 * @param vec_norm 2-norm of vec
 * @returns matrix exponential
 */
template <typename Derived>
inline auto Rodrigues(
    const Eigen::DenseBase<Derived>& vec, const typename Derived::Scalar& vec_norm) {
  typedef typename Derived::Scalar Scalar;
  ASSERT_EIGEN_VEC_SIZE(Derived, vec, 3);

  Eigen::Matrix<Scalar, 3, 3> skew_sym = Skew(vec.derived() / vec_norm);
  return Eigen::Matrix<Scalar, 3, 3>::Identity() + (std::sin(vec_norm) * skew_sym) +
         ((1.0 - std::cos(vec_norm)) * skew_sym * skew_sym);
}

/**
 * *=== Rodrigues ===*
 * @brief Rodrigues formula for the approximation of a matrix exponential
 * @tparam Derived An Eigen size 3 object (i.e Vector3d, RowVector3d, Array3f)
 * @param vec size 3 vector
 * @returns matrix exponential
 */
template <typename Derived>
inline auto Rodrigues(const Eigen::DenseBase<Derived>& vec) {
  typedef typename Derived::Scalar Scalar;
  ASSERT_EIGEN_VEC_SIZE(Derived, vec, 3);

  Scalar vec_norm = vec.derived().norm();
  return Rodrigues(vec, vec_norm);
}

/**
 * *=== Rodrigues4 ===*
 * @brief Rodrigues formula for the 4-th order approximation of a matrix exponential
 * @tparam Derived An Eigen size 3 object (i.e Vector3d, RowVector3d, Array3f)
 * @param vec      size 3 vector
 * @param vec_norm 2-norm of vec
 * @returns matrix exponential
 */
template <typename Derived>
inline auto Rodrigues4(
    const Eigen::DenseBase<Derived>& vec, const typename Derived::Scalar& vec_norm) {
  typedef typename Derived::Scalar Scalar;
  ASSERT_EIGEN_VEC_SIZE(Derived, vec, 3);

  Eigen::Matrix3<Scalar> skew_sym = Skew(vec);
  Scalar norm_squared = vec_norm * vec_norm;
  return Eigen::Matrix3<Scalar>::Identity() + ((1.0 - (norm_squared / 6.0)) * skew_sym) +
         ((0.5 - (norm_squared / 24.0)) * skew_sym * skew_sym);
}

/**
 * *=== Rodrigues4 ===*
 * @brief Rodrigues formula for the 4-th order approximation of a matrix exponential
 * @tparam Derived An Eigen size 3 object (i.e Vector3d, RowVector3d, Array3f)
 * @param vec size 3 vector
 * @returns matrix exponential
 */
template <typename Derived>
inline auto Rodrigues4(const Eigen::DenseBase<Derived>& vec) {
  typedef typename Derived::Scalar Scalar;
  ASSERT_EIGEN_VEC_SIZE(Derived, vec, 3);

  Scalar vec_norm = vec.derived().norm();
  return Rodrigues4(vec, vec_norm);
}

/**
 * *=== scalar2expm ===
 * @brief Converts scalar value into 2x2 matrix corresponding to rotating by the value in radians
 *        (positive CCW)
 * @tparam T A floating point type (e.g. double, float)
 * @param scalar rotation angle as some kind of floating point type
 * @returns 2x2 rotation matrix
 */
template <typename T = double>
inline auto scalar2expm(const T& scalar) {
  T cs = std::cos(scalar);
  T ss = std::sin(scalar);
  return Eigen::Matrix2<T>({{cs, -ss}, {ss, cs}});
}

/**
 * *=== vec2expm ===*
 * @brief Converts vector into its matrix exponential approximation form
 * @tparam Derived An Eigen size 3 object (i.e Vector3d, RowVector3d, Array3f)
 * @param vec size 3 vector
 * @returns 3x3 matrix exponential
 */
template <typename Derived>
inline auto vec2expm(const Eigen::DenseBase<Derived>& vec) {
  typedef typename Derived::Scalar Scalar;
  static constexpr int RACT = Derived::RowsAtCompileTime;
  if constexpr (RACT == Eigen::Dynamic) {
    assert(vec.rows() == 1 || vec.rows() == 3);
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatX;
    if (vec.rows() == 1) {
      return MatX(scalar2expm<Scalar>(vec(0)));
    } else {
      Scalar vec_norm = vec.derived().norm();
      if (vec_norm < 0.02) {
        return MatX(Rodrigues4<Derived>(vec, vec_norm));
      } else {
        return MatX(Rodrigues<Derived>(vec, vec_norm));
      }
    }
  } else {
    static_assert(RACT == 1 || RACT == 3);
    if constexpr (RACT == 1) {
      return scalar2expm<Scalar>(vec(0));
    } else {
      Scalar vec_norm = vec.derived().norm();
      if (vec_norm < 0.02) {
        return Rodrigues4(vec, vec_norm);
      } else {
        return Rodrigues(vec, vec_norm);
      }
    }
  }
}

/**
 * *=== expm2vec ===*
 * @brief Converts matrix exponential approximation into its vector form
 * @tparam Derived An Eigen size 3x3 object (i.e Matrix3d, Array3f3)
 * @param mat 3x3 matrix exponential
 * @returns size 3 vector
 */
template <typename Derived>
inline auto expm2vec(const Eigen::DenseBase<Derived>& mat) {
  ASSERT_EIGEN_MAT_SIZE(Derived, mat, 3, 3);
  using Scalar = typename Derived::Scalar;

  Scalar phi = std::acos((mat.trace() - 1.0) / 2.0);
  if (phi == 0.0) {
    return Eigen::Vector3<Scalar>::Zero();
  }
  return phi * deskew(mat - mat.transpose()) / (2.0 * std::sin(phi));
}

}  // namespace nt

#endif