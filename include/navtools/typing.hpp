#ifndef NAVTOOLS_TYPING_HPP
#define NAVTOOLS_TYPING_HPP

#include <Eigen/Dense>

namespace navtools {

//* ===== Linalg Library Calls ================================================================= *//

template <typename T = double>          //
using Mat4x4 = Eigen::Matrix<T, 4, 4>;  // 3x3 matrix of Ts
template <typename T = double>          //
using Mat3x3 = Eigen::Matrix<T, 3, 3>;  // 3x3 matrix of Ts
template <typename T = double>          //
using Vec2 = Eigen::Vector<T, 2>;       // size 2 vector of Ts
template <typename T = double>          //
using Vec3 = Eigen::Vector<T, 3>;       // size 3 vector of Ts
template <typename T = double>          //
using Vec4 = Eigen::Vector<T, 4>;       // size 4 vector of Ts

}  // namespace navtools

#endif
