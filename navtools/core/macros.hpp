#ifndef NAVTOOLS_CORE_MACROS_HPP
#define NAVTOOLS_CORE_MACROS_HPP

#include <Eigen/Dense>

#define NEW_FP_CONST(name, value) \
  template <typename T = double>  \
  inline constexpr T name = static_cast<T>(value)

#define NEW_EIGEN_CONST(name, rows, cols, ...)                   \
  template <typename T = double>                                 \
  inline static const Eigen::Matrix<T, rows, cols> name = []() { \
    Eigen::Matrix<T, rows, cols> m;                              \
    m << __VA_ARGS__;                                            \
    return m;                                                    \
  }();

#define ASSERT_EIGEN_NUM_ROWS(Derived, obj, Rows)                                      \
  if constexpr (Derived::RowsAtCompileTime != Eigen::Dynamic) {                        \
    static_assert(Derived::RowsAtCompileTime == Rows, "Inappropriate number of rows"); \
  } else {                                                                             \
    assert(obj.rows() == Rows);                                                        \
  }

#define ASSERT_EIGEN_NUM_COLS(Derived, obj, Cols)                                         \
  if constexpr (Derived::ColsAtCompileTime != Eigen::Dynamic) {                           \
    static_assert(Derived::ColsAtCompileTime == Cols, "Inappropriate number of columns"); \
  } else {                                                                                \
    assert(obj.cols() == Cols);                                                           \
  }

#define ASSERT_EIGEN_VEC_SIZE(Derived, obj, Size)                                   \
  static_assert(                                                                    \
      (Derived::RowsAtCompileTime == 1 || Derived::ColsAtCompileTime == 1),         \
      "Eigen object is not vector type.");                                          \
  if constexpr (                                                                    \
      Derived::RowsAtCompileTime != Eigen::Dynamic &&                               \
      Derived::ColsAtCompileTime != Eigen::Dynamic) {                               \
    static_assert(                                                                  \
        (Derived::RowsAtCompileTime == Size || Derived::ColsAtCompileTime == Size), \
        "Inappropriate vector size");                                               \
  } else {                                                                          \
    assert(obj.size() == Size);                                                     \
  }

#define ASSERT_EIGEN_MAT_SIZE(Derived, obj, Rows, Cols) \
  ASSERT_EIGEN_NUM_ROWS(Derived, obj, Rows)             \
  ASSERT_EIGEN_NUM_COLS(Derived, obj, Cols)

#define ASSERT_EIGEN_TYPE(Derived1, Derived2)                               \
  static_assert(                                                            \
      std::is_same_v<typename Derived1::Scalar, typename Derived2::Scalar>, \
      "Eigen scalar type constraint not upheld")

#define ASSERT_EIGEN_SCALAR_TYPE(Derived, Scalar) \
  static_assert(                                  \
      std::is_same_v<typename Derived::Scalar, Scalar>, "Eigen scalar type constraint not upheld")

#endif