#ifndef KINEMATICS_INCLUDE_ROTATIONS_HPP
#define KINEMATICS_INCLUDE_ROTATIONS_HPP

#include <cmath>
#include <numbers>

#include <Eigen/Core>

#include "earth_info.hpp"

/*
E - ECEF
I - ECI
L - local tangent plane (defined by LLA)
*/

template<typename T>
constexpr T Deg2Rad(T input)
{
  return input * std::numbers::pi_v<T> / 180.0;
}

template<typename T>
constexpr T Rad2Deg(T input)
{
  return input * 180.0 / std::numbers::pi_v<T>;
}

template<typename Derived>
Eigen::Matrix<typename Derived::Scalar,3,3> SkewSymmetrize(const Eigen::MatrixBase<Derived>& vec)
{
  static_assert((Derived::ColsAtCompileTime == 1) || (Derived::ColsAtCompileTime == Eigen::Dynamic));
  if constexpr (Derived::ColsAtCompileTime == Eigen::Dynamic) {
    assert(vec.cols() == 1);
  }
  static_assert((Derived::RowsAtCompileTime == 3) || (Derived::RowsAtCompileTime == Eigen::Dynamic));
  if constexpr (Derived::RowsAtCompileTime == Eigen::Dynamic) {
    assert(vec.rows() == 3);
  }

  Eigen::Matrix<typename Derived::Scalar,3,3> result = Eigen::Matrix<typename Derived::Scalar,3,3>::Zero();
  result(2,1) = vec(0);
  result(1,2) = -vec(0);
  result(2,0) = -vec(1);
  result(0,2) = vec(1);
  result(1,0) = vec(2);
  result(0,1) = -vec(2);
  return result;
}


template<class Derived>
auto InverseSkewSymmetrize(const Eigen::MatrixBase<Derived>& mat)
{
  static_assert((Derived::ColsAtCompileTime == 3) || (Derived::ColsAtCompileTime == Eigen::Dynamic));
  if constexpr (Derived::ColsAtCompileTime == Eigen::Dynamic) {
    assert(mat.cols() == 3);
  }
  static_assert((Derived::RowsAtCompileTime == 3) || (Derived::RowsAtCompileTime == Eigen::Dynamic));
  if constexpr (Derived::RowsAtCompileTime == Eigen::Dynamic) {
    assert(mat.rows() == 3);
  }
  assert(mat.transpose() == -mat);

  Eigen::Vector<typename Derived::Scalar,3> result;
  result(0) = mat(2,1);
  result(1) = mat(0,2);
  result(2) = mat(1,0);
  return result;
}


template<typename Derived>
Eigen::Matrix<typename Derived::Scalar,3,3> RodriguesFormula(const Eigen::MatrixBase<Derived>& vec,
  const typename Derived::Scalar& vec_norm)
{
  static_assert((Derived::ColsAtCompileTime == 1) || (Derived::ColsAtCompileTime == Eigen::Dynamic));
  if constexpr (Derived::ColsAtCompileTime == Eigen::Dynamic) {
    assert(vec.cols() == 1);
  }
  static_assert((Derived::RowsAtCompileTime == 3) || (Derived::RowsAtCompileTime == Eigen::Dynamic));
  if constexpr (Derived::RowsAtCompileTime == Eigen::Dynamic) {
    assert(vec.rows() == 3);
  }
  assert(vec_norm > 0);

  typedef typename Derived::Scalar Scalar;

  Eigen::Matrix<Scalar,3,3> skew_sym = SkewSymmetrize(vec / vec_norm);
  return Eigen::Matrix<Scalar,3,3>::Identity()
        + ( std::sin(vec_norm) * skew_sym )
        + ( (1. - std::cos(vec_norm)) * skew_sym * skew_sym );
}

template<typename Derived>
Eigen::Matrix<typename Derived::Scalar,3,3> RodriguesFormula(const Eigen::MatrixBase<Derived>& vec)
{
  return RodriguesFormula(vec, vec.norm());
}


// Rodriguez formula 4th order approximation - for use when angles are small
template<typename Derived>
Eigen::Matrix<typename Derived::Scalar,3,3> RodriguesFormula4(const Eigen::MatrixBase<Derived>& vec,
  const typename Derived::Scalar& vec_norm)
{
  static_assert((Derived::ColsAtCompileTime == 1) || (Derived::ColsAtCompileTime == Eigen::Dynamic));
  if constexpr (Derived::ColsAtCompileTime == Eigen::Dynamic) {
    assert(vec.cols() == 1);
  }
  static_assert((Derived::RowsAtCompileTime == 3) || (Derived::RowsAtCompileTime == Eigen::Dynamic));
  if constexpr (Derived::RowsAtCompileTime == Eigen::Dynamic) {
    assert(vec.rows() == 3);
  }

  typedef typename Derived::Scalar Scalar;

  Eigen::Matrix<Scalar,3,3> skew_sym = SkewSymmetrize(vec);
  Scalar norm_squared = vec_norm * vec_norm;
  return Eigen::Matrix<Scalar,3,3>::Identity()
        + ( (1.0 - (norm_squared / 6.0)) * skew_sym )
        + ( (0.5 - (norm_squared / 24.0)) * skew_sym * skew_sym );
}

template<typename Derived>
Eigen::Matrix<typename Derived::Scalar,3,3> RodriguesFormula4(const Eigen::MatrixBase<Derived>& vec)
{
  return RodriguesFormula4(vec, vec.norm());
}


template<typename Derived>
auto VecExp(const Eigen::MatrixBase<Derived>& vec)
{
  typedef typename Derived::Scalar Scalar;

  constexpr Scalar threshold = 0.02;
  Scalar norm = vec.norm();
  if (norm < threshold) {
    return RodriguesFormula4(vec,norm);
  } else {
    return RodriguesFormula(vec,norm);
  }
}


template<class Derived>
void VecLog(Eigen::Vector<typename Derived::Scalar,3>& result, const Eigen::MatrixBase<Derived>& mat)
{
  static_assert((Derived::ColsAtCompileTime == 3) || (Derived::ColsAtCompileTime == Eigen::Dynamic));
  if constexpr (Derived::ColsAtCompileTime == Eigen::Dynamic) {
    assert(mat.cols() == 3);
  }
  static_assert((Derived::RowsAtCompileTime == 3) || (Derived::RowsAtCompileTime == Eigen::Dynamic));
  if constexpr (Derived::RowsAtCompileTime == Eigen::Dynamic) {
    assert(mat.rows() == 3);
  }

  typedef typename Derived::Scalar Scalar;

  Scalar phi = std::acos((mat.trace() - 1.0) / 2.0);
  if (phi == 0.0) {
    result = Eigen::Vector<Scalar,3>::Zero();
    return;
  }
  result = phi * InverseSkewSymmetrize(mat - mat.transpose()) / (2.0 * std::sin(phi));
}


template<class Derived>
auto VecLog(const Eigen::MatrixBase<Derived>& mat)
{
  Eigen::Vector<typename Derived::Scalar,3> result;
  VecLog(result, mat);
  return result;
}


using Eigen::Matrix, Eigen::Vector;

// latitude, longitude, altitude
template<class Derived>
void EtoL(Eigen::Matrix<typename Derived::Scalar,3,3>& result, 
        const Eigen::DenseBase<Derived>& lla)
{
  typedef typename Derived::Scalar Scalar;

  Scalar slat = std::sin(lla(0));
  Scalar slon = std::sin(lla(1));
  Scalar clat = std::cos(lla(0));
  Scalar clon = std::cos(lla(1));

  result(0,0) = -slat * clon;
  result(0,1) = -slat * slon;
  result(0,2) = clat;
  result(1,0) = -slon;
  result(1,1) = clon;
  result(1,2) = 0.;
  result(2,0) = -clat * clon;
  result(2,1) = -clat * slon;
  result(2,2) = -slat;
}


template<class Derived>
auto EtoL(const Eigen::DenseBase<Derived>& lla)
{
  Matrix<typename Derived::Scalar,3,3> result;
  EtoL(result,lla);
  return result;
}

template<class Derived>
void LtoE(Eigen::Matrix<typename Derived::Scalar,3,3>& result, 
        const Eigen::DenseBase<Derived>& lla)
{
  typedef typename Derived::Scalar Scalar;

  Scalar slat = std::sin(lla(0));
  Scalar slon = std::sin(lla(1));
  Scalar clat = std::cos(lla(0));
  Scalar clon = std::cos(lla(1));

  result(0,0) = -slat * clon;
  result(1,0) = -slat * slon;
  result(2,0) = clat;
  result(0,1) = -slon;
  result(1,1) = clon;
  result(2,1) = 0.;
  result(0,2) = -clat * clon;
  result(1,2) = -clat * slon;
  result(2,2) = -slat;
}

template<class Derived>
auto LtoE(const Eigen::DenseBase<Derived>& lla)
{
  Matrix<typename Derived::Scalar,3,3> result;
  LtoE(result,lla);
  return result;
}


template<typename Scalar>
void ItoE(Matrix<Scalar,3,3> result, Scalar time, Scalar last_alignment_time)
{
  Scalar angle = WGS84::EarthRate<Scalar> * (time - last_alignment_time);
  result(2,0) = result(2,1) = result(0,2) = result(1,2) = 0.0;
  result(0,0) = std::cos(angle);
  result(1,1) = result(0,0);
  result(0,1) = std::sin(angle);
  result(1,0) = -result(0,1);
  result(2,2) = 1.0;
}

template<typename Scalar>
auto ItoE(Scalar time, Scalar last_alignment_time)
{
  Matrix<Scalar,3,3> result;
  ItoE(result, time, last_alignment_time);
  return result;
}

template<typename Scalar>
void EtoI(Matrix<Scalar,3,3> result, Scalar time, Scalar last_alignment_time)
{
  Scalar angle = WGS84::EarthRate<Scalar> * (time - last_alignment_time);
  result(2,0) = result(2,1) = result(0,2) = result(1,2) = 0.0;
  result(0,0) = std::cos(angle);
  result(1,1) = result(0,0);
  result(1,0) = std::sin(angle);
  result(0,1) = -result(0,1);
  result(2,2) = 1.0;
}

template<typename Scalar>
auto EtoI(Scalar time, Scalar last_alignment_time)
{
  Matrix<Scalar,3,3> result;
  EtoI(result, time, last_alignment_time);
  return result;
}

#endif