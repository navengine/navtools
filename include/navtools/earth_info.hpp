#ifndef NAVENGINE_INCLUDE_NAVTOOLS_EARTH_INFO_HPP 
#define NAVENGINE_INCLUDE_NAVTOOLS_EARTH_INFO_HPP

#include <cmath>

#include <Eigen/Dense>

#include "constants.hpp"

namespace WGS84 {

NAV_DEFINE_FP_CONSTANT(EquatorialRadius, 6378137.0); // meters
NAV_DEFINE_FP_CONSTANT(PolarRadius, 6356752.31425);
NAV_DEFINE_FP_CONSTANT(f, 1.0 / 298.257223563);
NAV_DEFINE_FP_CONSTANT(e, 0.0818191908425);
NAV_DEFINE_FP_CONSTANT(e_squared, e<Scalar> * e<Scalar>);
NAV_DEFINE_FP_CONSTANT(EarthRate, 7.292115e-5);

template<typename FloatType>
FloatType R_E(FloatType latitude)
{
  FloatType slat = std::sin(latitude);
  return EquatorialRadius<FloatType> / std::sqrt(1. - (e_squared<FloatType> * slat * slat));
}

template<typename FloatType>
FloatType R_L(FloatType latitude)
{
  FloatType slat = std::sin(latitude);
  return EquatorialRadius<FloatType> * (1. - e_squared<FloatType>)
          * std::pow(1. - (e_squared<FloatType> * slat * slat), -1.5);
}

template<typename FloatType>
void LLAtoECEF(Eigen::Vector<FloatType,3>& result, FloatType latitude, FloatType longitude, FloatType altitude)
{
  FloatType eastward_radius = R_E<FloatType>(latitude); 
  FloatType xy_mag = (eastward_radius + altitude) * std::cos(latitude);
  result(0) = xy_mag * std::cos(longitude);
  result(1) = xy_mag * std::sin(longitude);
  result(2) = (((1 - e_squared<FloatType>) * eastward_radius) + altitude) * std::sin(latitude);
}

template<typename FloatType>
Eigen::Vector<FloatType,3> LLAtoECEF(FloatType latitude, FloatType longitude, FloatType altitude)
{
  Eigen::Vector<FloatType,3> result;
  LLAtoECEF(result, latitude, longitude, altitude);
  return result;
}

template<typename FloatType>
Eigen::Vector<FloatType,3> LLAtoECEF(Eigen::Vector<FloatType,3>& lla)
{
  return LLAtoECEF(lla(0),lla(1),lla(2));
}


// TODO add gravity model

} // namespace WGS84


// namespace GRS80 {

// template<typename FloatType>
// inline constexpr FloatType EquatorialRadius = 6378137.0; // meters 

// template<typename FloatType>
// inline constexpr FloatType PolarRadius = 6356752.31414; // meters 

// template<typename FloatType>
// inline constexpr FloatType f = 1.0 / 298.257222101; // flattening

// template<typename FloatType>
// inline constexpr FloatType e = 0.0818191910428; // eccentricity

// template<typename FloatType>
// inline constexpr FloatType e_squared = e<FloatType> * e<FloatType>; // eccentricity squared

// template<typename FloatType>
// FloatType R_E(FloatType latitude)
// {
//   FloatType slat = std::sin(latitude);
//   return EquatorialRadius<FloatType> / std::sqrt(1. - (e_squared<FloatType> * slat * slat));
// }

// template<typename FloatType>
// FloatType R_L(FloatType latitude)
// {
//   FloatType slat = std::sin(latitude);
//   return EquatorialRadius<FloatType> * (1. - e_squared<FloatType>)
//           * std::pow(1. - (e_squared<FloatType> * slat * slat), -1.5);
// }

// template<typename FloatType>
// void LLAtoECEF(Eigen::Vector<FloatType,3>& result, FloatType latitude, FloatType longitude, FloatType altitude)
// {
//   FloatType eastward_radius = R_E<FloatType>(latitude); 
//   FloatType xy_mag = (eastward_radius + altitude) * std::cos(latitude);
//   result(0) = xy_mag * std::cos(longitude);
//   result(1) = xy_mag * std::sin(longitude);
//   result(2) = (((1 - e_squared<FloatType>) * eastward_radius) + altitude) * std::sin(latitude);
// }

// template<typename FloatType>
// Eigen::Vector<FloatType,3> LLAtoECEF(FloatType latitude, FloatType longitude, FloatType altitude)
// {
//   Eigen::Vector<FloatType,3> result;
//   LLAtoECEF(result, latitude, longitude, altitude);
//   return result;
// }

// } // namespace GRS80


#endif
