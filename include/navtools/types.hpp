/**
|========================================== typing.hpp ============================================|
|                                                                                                  |
|   @file     include/navtools/types.hpp                                                           |
|   @brief    Useful type definitions.                                                             |
|   @date     July 2024                                                                            |
|                                                                                                  |
|==================================================================================================|
*/

#ifndef NAVTOOLS_TYPES_HPP
#define NAVTOOLS_TYPES_HPP

#include <cassert>
#include <type_traits>

#include <Eigen/Dense>

#define DEFINE_FP_CONSTANT(name, value) \
    template <typename T = double>      \
    inline constexpr T name = static_cast<T>(value)

namespace navtools {

//* ===== Linalg Library Calls ================================================================= *//

template <typename T = double>
using Mat2x2 = Eigen::Matrix<T, 2, 2>;  // 2x2 matrix of Ts
template <typename T = double>
using Mat3x3 = Eigen::Matrix<T, 3, 3>;  // 3x3 matrix of Ts
template <typename T = double>
using Mat4x4 = Eigen::Matrix<T, 4, 4>;  // 4x4 matrix of Ts

template <typename T = double>
using Mat2 = Eigen::Matrix<T, 2, 2>;  // 2x2 matrix of Ts
template <typename T = double>
using Mat3 = Eigen::Matrix<T, 3, 3>;  // 3x3 matrix of Ts
template <typename T = double>
using Mat4 = Eigen::Matrix<T, 4, 4>;  // 4x4 matrix of Ts

template <typename T = double>
using Vec2 = Eigen::Vector<T, 2>;  // size 2 vector of Ts
template <typename T = double>
using Vec3 = Eigen::Vector<T, 3>;  // size 3 vector of Ts
template <typename T = double>
using Vec4 = Eigen::Vector<T, 4>;  // size 4 vector of Ts

using Eigen::DenseBase;

// Using macros instead of functions so that the location in code is traceable when assertions fail
#define ASSERT_EIGEN_NUM_ROWS(T,obj,Rows)                 \
  if constexpr (T::RowsAtCompileTime != Eigen::Dynamic) { \
    static_assert(T::RowsAtCompileTime == Rows,           \
                  "Inappropriate number of rows");        \
  }                                                       \
  else {                                                  \
    assert(obj.rows() == Rows);                           \
  }

#define ASSERT_EIGEN_NUM_COLS(T,obj,Cols)                 \
  if constexpr (T::ColsAtCompileTime != Eigen::Dynamic) { \
    static_assert(T::ColsAtCompileTime == Cols,           \
                  "Inappropriate number of columns");     \
  }                                                       \
  else {                                                  \
    assert(obj.cols() == Cols);                           \
  }

#define ASSERT_EIGEN_OBJ_SIZE(T,obj,Rows,Cols) \
        ASSERT_EIGEN_NUM_ROWS(T,obj,Rows) \
        ASSERT_EIGEN_NUM_COLS(T,obj,Cols)

#define ASSERT_EIGEN_SCALAR_TYPE(Derived, T)                \
  static_assert(std::is_same_v<typename Derived::Scalar,T>, \
  "Eigen scalar type constraint not upheld")


// template<int Rows, typename Derived>
// void AssertNumRows(const Eigen::EigenBase<Derived>& obj)
// {
//   if constexpr (Derived::RowsAtCompileTime != Eigen::Dynamic) {
//     static_assert(Derived::RowsAtCompileTime == Rows,
//                   "Inappropriate number of rows");
//   }
//   else {
//     assert(obj.rows() == Rows);
//   }
// }
// 
// template<int Cols, typename Derived>
// void AssertNumCols(const Eigen::EigenBase<Derived>& obj)
// {
//   if constexpr (Derived::ColsAtCompileTime != Eigen::Dynamic) {
//     static_assert(Derived::ColsAtCompileTime == Cols,
//                   "Inappropriate number of columns");
//   }
//   else {
//     assert(obj.cols() == Cols);
//   }
// }
// 
// template<int Rows, int Cols, typename Derived>
// void AssertDimensions(const Eigen::EigenBase<Derived>& obj)
// {
//   AssertNumRows<Rows>(obj);
//   AssertNumCols<Cols>(obj);
// }

}  // namespace navtools

#endif
