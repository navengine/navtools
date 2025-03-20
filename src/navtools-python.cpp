/**
 * *navtools-python.cpp*
 *
 * =======  ========================================================================================
 * @file    src/navtools-python.cpp
 * @brief   PyBind11 wrapper for using navtools in python!
 * @date    March 2025
 * =======  ========================================================================================
 */

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include "navtools/attitude.hpp"
#include "navtools/binary-ops.hpp"
#include "navtools/constants.hpp"
#include "navtools/frames.hpp"
#include "navtools/math.hpp"
#include "navtools/models.hpp"

namespace py = pybind11;
using namespace navtools;

PYBIND11_MODULE(_navtools_core, h) {
  h.doc() =
      R"pbdoc(
      NavTools
      ========
      
      Common functions and definitions used for navigation. When using numpy arrays, ensure 
      that they are saved in column-wise contiguous memory (e.g. set order='F')!

      Contains the following submodules:

        1. `attitude`
        2. `binaryops`
        3. `frames`
        4. `math`
        5. `models`
      )pbdoc";

  h.attr("__version__") = "1.0.0";

  //! === Constants ================================================================================
  h.attr("PI") = PI<double>;
  h.attr("HALF_PI") = HALF_PI<double>;
  h.attr("TWO_PI") = TWO_PI<double>;
  h.attr("PI_SQU") = PI_SQU<double>;
  h.attr("SQRT_PI") = SQRT_PI<double>;
  h.attr("LIGHT_SPEED") = LIGHT_SPEED<double>;
  h.attr("BOLTZMANN") = BOLTZMANN<double>;
  h.attr("GAUSS_TO_TESLA") = GAUSS_TO_TESLA<double>;
  h.attr("METERS_TO_FOOT") = METERS_PER_FOOT<double>;
  h.attr("MINUTES_PER_DAY") = MINUTES_PER_DAY<double>;
  h.attr("GRAVITY") = GRAVITY<double>;
  h.attr("RE") = RE<double>;
  h.attr("J2") = J2<double>;
  h.attr("J3") = J3<double>;
  h.attr("J4") = J4<double>;
  h.attr("F") = F<double>;
  h.attr("WGS84_MU") = WGS84_MU<double>;
  h.attr("WGS84_R0") = WGS84_R0<double>;
  h.attr("WGS84_RP") = WGS84_RP<double>;
  h.attr("WGS84_E") = WGS84_E<double>;
  h.attr("WGS84_E2") = WGS84_E2<double>;
  h.attr("WGS84_F") = WGS84_F<double>;
  h.attr("WGS84_OMEGA") = WGS84_OMEGA<double>;
  h.attr("WGS84_OMEGA_VEC") = WGS84_OMEGA_VEC<double>;
  h.attr("WGS84_OMEGA_SKEW") = WGS84_OMEGA_SKEW<double>;
  h.attr("RAD2DEG") = RAD2DEG<double>;
  h.attr("DEG2RAD") = DEG2RAD<double>;
  h.def("deg2rad", [](double x) { return DEG2RAD<double> * x; });
  h.def("rad2deg", [](double x) { return RAD2DEG<double> * x; });

  //! === Attitude submodule =======================================================================
  py::module_ att = h.def_submodule("attitude", R"pbdoc(
        Attitude
        ========
        
        Attitude representations and conversions between them.)pbdoc");

  // euler2quat
  att.def(
      "euler2quat",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector4d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const bool>(&euler2quat<double>),
      py::arg("q"),
      py::arg("e"),
      py::arg("IsNed") = true,
      R"pbdoc(
      euler2quat
      ==========
      
      Converts euler angles (roll-pitch-yaw) to corresponding BODY-to-NAV quaternion

      Parameters
      ----------

      q : np.ndarray

          size 4 NAV quaternion

      e : np.ndarray

          size 3 RPY euler angles [rad]

      IsNed : bool

          Is desired frame NED, default is True
      )pbdoc");
  att.def(
      "euler2quat",
      py::overload_cast<const Eigen::Ref<const Eigen::Vector3d> &, const bool>(&euler2quat<double>),
      py::arg("q"),
      py::arg("IsNed") = true,
      R"pbdoc(
      euler2quat
      ==========
      
      Converts euler angles (roll-pitch-yaw) to corresponding BODY-to-NAV quaternion

      Parameters
      ----------

      q : np.ndarray

          size 4 NAV quaternion

      IsNed : bool

          Is desired frame NED, default is True 

      Returns
      -------

      e : np.ndarray

          size 3 RPY euler angles [rad]
      )pbdoc");

  // euler2dcm
  att.def(
      "euler2dcm",
      py::overload_cast<
          Eigen::Ref<Eigen::Matrix3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const bool>(&euler2dcm<double>),
      py::arg("C"),
      py::arg("e"),
      py::arg("IsNed") = true,
      R"pbdoc(
      euler2dcm
      =========
      
      Converts euler angles (roll-pitch-yaw) to corresponding BODY-to-NAV DCM
    
      Parameters
      ----------

      C : np.ndarray

          size 3x3 NAV DCM (ZYX)

      e : np.ndarray

          size 3 RPY euler angles [rad]

      IsNed : bool

          Is desired frame NED, default is True
      )pbdoc");
  att.def(
      "euler2dcm",
      py::overload_cast<const Eigen::Ref<const Eigen::Vector3d> &, const bool>(&euler2dcm<double>),
      py::arg("e"),
      py::arg("IsNed") = true,
      R"pbdoc(
      euler2dcm
      =========

      Converts euler angles (roll-pitch-yaw) to corresponding BODY-to-NAV DCM
    
      Parameters
      ----------

      e : np.ndarray

          size 3 RPY euler angles [rad]

      IsNed : bool

          Is desired frame NED, default is True 

      Returns
      -------

      C : np.ndarray

          size 3x3 NAV DCM (ZYX)
      )pbdoc");

  // quat2euler
  att.def(
      "quat2euler",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector4d> &,
          const bool>(&quat2euler<double>),
      py::arg("e"),
      py::arg("q"),
      py::arg("IsNed") = true,
      R"pbdoc(
      quat2euler
      ==========

      Converts BODY-to-NAV quaternion to corresponding euler angles (roll-pitch-yaw)
        
      Parameters
      ----------

      e : np.ndarray

          size 3 RPY euler angles [rad]

      q : np.ndarray

          size 4 NAV quaternion

      IsNed : bool

          Is desired frame NED, default is True 
      )pbdoc");
  att.def(
      "quat2euler",
      py::overload_cast<const Eigen::Ref<const Eigen::Vector4d> &, const bool>(&quat2euler<double>),
      py::arg("q"),
      py::arg("IsNed") = true,
      R"pbdoc(
      quat2euler
      ==========
      
      Converts BODY-to-NAV quaternion to corresponding euler angles (roll-pitch-yaw)
        
      Parameters
      ----------

      q : np.ndarray

          size 4 NAV quaternion

      IsNed : bool

          Is desired frame NED, default is True 

      Returns
      -------

      e : np.ndarray

          size 3 RPY euler angles [rad]
      )pbdoc");

  // quat2dcm
  att.def(
      "quat2dcm",
      py::overload_cast<Eigen::Ref<Eigen::Matrix3d>, const Eigen::Ref<const Eigen::Vector4d> &>(
          &quat2dcm<double>),
      py::arg("q"),
      py::arg("C"),
      R"pbdoc(
      quat2dcm
      ========
      
      Converts BODY-to-NAV quaternion to corresponding BODY-to-NAV DCM
        
      Parameters
      ----------

      q : np.ndarray

          size 4 NAV quaternion

      C : np.ndarray

          size 3x3 NAV DCM (ZYX)
      )pbdoc");
  att.def(
      "quat2dcm",
      py::overload_cast<const Eigen::Ref<const Eigen::Vector4d> &>(&quat2dcm<double>),
      py::arg("q"),
      R"pbdoc(
      quat2dcm
      ========
      
      Converts BODY-to-NAV quaternion to corresponding BODY-to-NAV DCM
        
      Parameters
      ----------

      q : np.ndarray

          size 4 NAV quaternion

      Returns
      -------

      C : np.ndarray

          size 3x3 NAV DCM (ZYX)
      )pbdoc");

  // dcm2euler
  att.def(
      "dcm2euler",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Matrix3d> &,
          const bool>(&dcm2euler<double>),
      py::arg("e"),
      py::arg("C"),
      py::arg("IsNed") = true,
      R"pbdoc(
      dcm2euler
      =========
      
      Converts BODY-to-NAV DCM to corresponding euler angles (roll-pitch-yaw)
    
      Parameters
      ----------

      e : np.ndarray

          size 3 RPY euler angles [rad]

      C : np.ndarray

          size 3x3 NAV DCM (ZYX)

      IsNed : bool

          Is desired frame NED, default is True 
      )pbdoc");
  att.def(
      "dcm2euler",
      py::overload_cast<const Eigen::Ref<const Eigen::Matrix3d> &, const bool>(&dcm2euler<double>),
      py::arg("C"),
      py::arg("IsNed") = true,
      R"pbdoc(
      dcm2euler
      =========

      Converts BODY-to-NAV DCM to corresponding euler angles (roll-pitch-yaw)
        
      Parameters
      ----------

      C : np.ndarray

          size 3x3 NAV DCM (ZYX)

      IsNed : bool

          Is desired frame NED, default is True 

      Returns
      -------

      e : np.ndarray

          size 3 RPY euler angles [rad]
      )pbdoc");

  // dcm2quat
  att.def(
      "dcm2quat",
      py::overload_cast<Eigen::Ref<Eigen::Vector4d>, const Eigen::Ref<const Eigen::Matrix3d> &>(
          &dcm2quat<double>),
      py::arg("q"),
      py::arg("C"),
      R"pbdoc(
        dcm2quat
        ========
  
        Converts BODY-to-NAV DCM to corresponding BODY-to-NAV quaternion
          
        Parameters
        ----------

        q : np.ndarray
  
            size 4 NAV quaternion
  
        C : np.ndarray
  
            size 3x3 NAV DCM (ZYX)
        )pbdoc");
  att.def(
      "dcm2quat",
      py::overload_cast<const Eigen::Ref<const Eigen::Matrix3d> &>(&dcm2quat<double>),
      py::arg("C"),
      R"pbdoc(
      dcm2quat
      ========

      Converts BODY-to-NAV DCM to corresponding BODY-to-NAV quaternion
        
      Parameters
      ----------

      C : np.ndarray

          size 3x3 NAV DCM (ZYX)

      Returns
      -------

      q : np.ndarray

          size 4 NAV quaternion
      )pbdoc");

  // RotX
  att.def(
      "RotX",
      py::overload_cast<Eigen::Ref<Eigen::Matrix3d>, const double &>(&RotX<double>),
      py::arg("C"),
      py::arg("x"),
      R"pbdoc(
      RotX
      ====

      Converts euler angle about x-axis into DCM rotation about x-axis

      Parameters
      ----------

      C : np.ndarray

          3x3 x-axis DCM rotation

      x : double

          euler angle [rad]
      )pbdoc");
  att.def(
      "RotX",
      py::overload_cast<const double &>(&RotX<double>),
      py::arg("x"),
      R"pbdoc(
      RotX
      ====

      Converts euler angle about x-axis into DCM rotation about xy-axis

      Parameters
      ----------

      x : double

          euler angle [rad]

      Returns
      -------

          3x3 x-axis DCM rotation
      )pbdoc");

  // RotY
  att.def(
      "RotY",
      py::overload_cast<Eigen::Ref<Eigen::Matrix3d>, const double &>(&RotY<double>),
      py::arg("C"),
      py::arg("y"),
      R"pbdoc(
      RotY
      ====

      Converts euler angle about y-axis into DCM rotation about y-axis

      Parameters
      ----------

      C : np.ndarray

          3x3 y-axis DCM rotation

      y : double

          euler angle [rad]
      )pbdoc");
  att.def(
      "RotY",
      py::overload_cast<const double &>(&RotY<double>),
      py::arg("y"),
      R"pbdoc(
      RotY
      ====

      Converts euler angle about y-axis into DCM rotation about y-axis

      Parameters
      ----------

      y : double

          euler angle [rad]

      Returns
      -------

          3x3 y-axis DCM rotation
      )pbdoc");

  // RotZ
  att.def(
      "RotZ",
      py::overload_cast<Eigen::Ref<Eigen::Matrix3d>, const double &>(&RotZ<double>),
      py::arg("C"),
      py::arg("z"),
      R"pbdoc(
      RotZ
      ====

      Converts euler angle about z-axis into DCM rotation about z-axis

      Parameters
      ----------

      C : np.ndarray

          3x3 z-axis DCM rotation

      z : double

          euler angle [rad]
      )pbdoc");
  att.def(
      "RotZ",
      py::overload_cast<const double &>(&RotZ<double>),
      py::arg("z"),
      R"pbdoc(
      RotZ
      ====

      Converts euler angle about z-axis into DCM rotation about z-axis

      Parameters
      ----------

      z : double

          euler angle [rad]

      Returns
      -------

          3x3 z-axis DCM rotation
      )pbdoc");

  //! === Binary-Ops submodule =====================================================================
  py::module_ bin = h.def_submodule("binaryops", R"pbdoc(
        Binary-Ops
        ==========
        
        Useful binary operations.)pbdoc");

  // SetBit
  bin.def(
      "SetBit",
      &SetBit<false>,
      py::arg("x"),
      py::arg("n"),
      R"pbdoc(
      SetBit
      ======

      Set a data bit to 1
    
      Parameters
      ----------

      x : uint32

          Number to modify

      n : uint8

          Position of bit to set (Position 0 is MSB and 31 is LSB by default)
      )pbdoc");

  // UnsetBit
  bin.def(
      "UnsetBit",
      &UnsetBit<false>,
      py::arg("x"),
      py::arg("n"),
      R"pbdoc(
      UnsetBit
      ========
      
      Set a data bit to 0

      Parameters
      ----------

      x : uint32

          Number to modify

      n : uint8

          Position of bit to set (Position 0 is MSB and 31 is LSB by default)
      )pbdoc");

  // SetBitTo
  bin.def(
      "SetBitTo",
      py::overload_cast<uint32_t &, const uint8_t, const bool>(&SetBitTo<false>),
      py::arg("x"),
      py::arg("n"),
      py::arg("b"),
      R"pbdoc(
      SetBitTo
      ========

      Set a data bit to 1

      Parameters
      ----------

      x : uint32

          Number to modify

      n : uint8

          Position of bit to set (Position 0 is MSB and 31 is LSB by default)

      b : bool

          New bit value
      )pbdoc");

  // GetBit
  bin.def(
      "GetBit",
      py::overload_cast<const uint32_t &, const uint8_t>(&GetBit<false>),
      py::arg("x"),
      py::arg("n"),
      R"pbdoc(
      GetBit
      ========

      Check the value of a bit

      Parameters
      ----------

      x : uint32

          Number to modify

      n : uint8

          Position of bit to set (Position 0 is MSB and 31 is LSB by default)
      )pbdoc");

  // GetBits
  bin.def(
      "GetBits",
      &GetBits<false>,
      py::arg("x"),
      py::arg("b"),
      py::arg("n"),
      R"pbdoc(
      GetBits
      =========

      Check the value of multiple bits in series

      Parameters
      ----------

      x : uint32

          Number to modify

      b : unit8

          Position of first bit to check (Position 0 is MSB and 31 is LSB by default)

      n : uint8

          Position of last bit to set (Position 0 is MSB and 31 is LSB by default)
      )pbdoc");

  // MultiXor
  bin.def(
      "MultiXor",
      py::overload_cast<const uint32_t &, const uint8_t[], const int>(&MultiXor<false>),
      py::arg("x"),
      py::arg("n"),
      py::arg("Size"),
      R"pbdoc(
      MultiXor
      ========

      Perform XOR operation on multiple bits

      Parameters
      ----------

      x : uint32

          Number to modify

      n : np.ndarray

          Positions of bits to check (Position 0 is MSB and 31 is LSB by default)

      Size: unit8

          Size of the input array
      )pbdoc");

  // TwosComp
  bin.def(
      "TwosComp",
      &TwosComp,
      py::arg("x"),
      py::arg("n"),
      R"pbdoc(
      TwosComp
      ========

      Two's compliment to create a signed integer

      Parameters
      ----------

      x : uint32

          Unsigned integer to use

      n : uint8

          Number of bits in the integer

      Returns
      -------

      y : int32

          Signed integer result
      )pbdoc");

  //! === Frames submodule =========================================================================
  py::module_ frm = h.def_submodule("frames", R"pbdoc(
        Frames
        ======
        
        Common coordinate frame transformations.)pbdoc");

  // eci2ecefDcm
  frm.def(
      "eci2ecefDcm",
      py::overload_cast<Eigen::Ref<Eigen::Matrix3d>, const double &>(&eci2ecefDcm<double>),
      py::arg("C_i_e"),
      py::arg("dt"),
      R"pbdoc(
      eci2ecefDcm
      ===========

      Earth-Centered-Inertial to Earth-Centered-Earth-Fixed direction cosine matrix

      Parameters
      ----------

      C_i_e : np.ndarray

          3x3 ECI->ECEF direction cosine matrix

      dt : double 
          time elapsed between frames [s]
      )pbdoc");
  frm.def(
      "eci2ecefDcm",
      py::overload_cast<const double &>(&eci2ecefDcm<double>),
      py::arg("dt"),
      R"pbdoc(
      eci2ecefDcm
      ===========

      Earth-Centered-Inertial to Earth-Centered-Earth-Fixed direction cosine matrix

      Parameters
      ----------

      dt : double 

          time elapsed between frames [s]

      Returns
      -------

      C_i_e : np.ndarray

          3x3 ECI->ECEF direction cosine matrix
      )pbdoc");

  // eci2nedDcm
  frm.def(
      "eci2nedDcm",
      py::overload_cast<
          Eigen::Ref<Eigen::Matrix3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2nedDcm<double>),
      py::arg("C_i_n"),
      py::arg("lla"),
      py::arg("dt"),
      R"pbdoc(
      ecef2nedDcm
      ===========

      Earth-Centered-Inertial to North-East-Down direction cosine matrix

      Parameters
      ----------

      C_i_n : np.ndarray

          3x3 ECI->NED direction cosine matrix

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]

      dt : double 

          time elapsed between frames [s]
      )pbdoc");
  frm.def(
      "eci2nedDcm",
      py::overload_cast<const Eigen::Ref<const Eigen::Vector3d> &, const double &>(
          &eci2nedDcm<double>),
      py::arg("lla"),
      py::arg("dt"),
      R"pbdoc(
      eci2nedDcm
      ==========

      Earth-Centered-Inertial to North-East-Down direction cosine matrix

      Parameters
      ----------

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]

      dt : double 

          time elapsed between frames [s]

      Returns
      -------

      C_i_n : np.ndarray

          3x3 ECI->NED direction cosine matrix
      )pbdoc");

  // eci2enuDcm
  frm.def(
      "eci2enuDcm",
      py::overload_cast<
          Eigen::Ref<Eigen::Matrix3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2enuDcm<double>),
      py::arg("C_i_n"),
      py::arg("lla"),
      py::arg("dt"),
      R"pbdoc(
      eci2enuDcm
      ==========

      Earth-Centered-Inertial to East-North-Up direction cosine matrix

      Parameters
      ----------

      C_i_n : np.ndarray

          3x3 ECI->ENU direction cosine matrix

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]

      dt : double 

          time elapsed between frames [s]
      )pbdoc");
  frm.def(
      "eci2enuDcm",
      py::overload_cast<const Eigen::Ref<const Eigen::Vector3d> &, const double &>(
          &eci2enuDcm<double>),
      py::arg("lla"),
      py::arg("dt"),
      R"pbdoc(
      eci2enuDcm
      ==========

      Earth-Centered-Inertial to East-North-Up direction cosine matrix

      Parameters
      ----------

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]
      dt : double 

          time elapsed between frames [s]

      Returns
      -------

      C_i_n : np.ndarray

          3x3 ECI->ENU direction cosine matrix
      )pbdoc");

  // ecef2eciDcm
  frm.def(
      "ecef2eciDcm",
      py::overload_cast<Eigen::Ref<Eigen::Matrix3d>, const double &>(&ecef2eciDcm<double>),
      py::arg("C_e_i"),
      py::arg("dt"),
      R"pbdoc(
      ecef2eciDcm
      ===========

      Earth-Centered-Earth-Fixed to Earth-Centered-Inertial direction cosine matrix

      Parameters
      ----------

      C_e_i : np.ndarray

          3x3 ECEF->ECI direction cosine matrix

      dt : double 

          time elapsed between frames [s]
      )pbdoc");
  frm.def(
      "ecef2eciDcm",
      py::overload_cast<const double &>(&ecef2eciDcm<double>),
      py::arg("dt"),
      R"pbdoc(
      ecef2eciDcm
      ===========

      Earth-Centered-Earth-Fixed to Earth-Centered-Inertial direction cosine matrix

      Parameters
      ----------

      dt : double 

          time elapsed between frames [s]

      Returns
      -------

      C_e_i : np.ndarray

          3x3 ECEF->ECI direction cosine matrix
      )pbdoc");

  // ecef2nedDcm
  frm.def(
      "ecef2nedDcm",
      py::overload_cast<Eigen::Ref<Eigen::Matrix3d>, const Eigen::Ref<const Eigen::Vector3d> &>(
          &ecef2nedDcm<double>),
      py::arg("C_e_n"),
      py::arg("lla"),
      R"pbdoc(
      ecef2nedDcm
      ===========

      Earth-Centered-Earth-Fixed to North-East-Down direction cosine matrix

      Parameters
      ----------

      C_e_n : np.ndarray

          3x3 ECEF->NED direction cosine matrix

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]
      )pbdoc");
  frm.def(
      "ecef2nedDcm",
      py::overload_cast<const Eigen::Ref<const Eigen::Vector3d> &>(&ecef2nedDcm<double>),
      py::arg("lla"),
      R"pbdoc(
      ecef2nedDcm
      ===========

      Earth-Centered-Earth-Fixed to North-East-Down direction cosine matrix

      Parameters
      ----------

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]

      Returns
      -------

      C_e_n : np.ndarray

          3x3 ECEF->NED direction cosine matrix
      )pbdoc");

  // ecef2enuDcm
  frm.def(
      "ecef2enuDcm",
      py::overload_cast<Eigen::Ref<Eigen::Matrix3d>, const Eigen::Ref<const Eigen::Vector3d> &>(
          &ecef2enuDcm<double>),
      py::arg("C_e_n"),
      py::arg("lla"),
      R"pbdoc(
      ecef2enuDcm
      ===========

      Earth-Centered-Earth-Fixed to East-North-Up direction cosine matrix

      Parameters
      ----------

      C_e_n : np.ndarray

          3x3 ECEF->ENU direction cosine matrix

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]
      )pbdoc");
  frm.def(
      "ecef2enuDcm",
      py::overload_cast<const Eigen::Ref<const Eigen::Vector3d> &>(&ecef2enuDcm<double>),
      py::arg("lla"),
      R"pbdoc(
      ecef2enuDcm
      ===========

      Earth-Centered-Earth-Fixed to East-North-Up direction cosine matrix

      Parameters
      ----------

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]

      Returns
      -------

      C_e_n : np.ndarray

          3x3 ECEF->ENU direction cosine matrix
      )pbdoc");

  // ned2eciDcm
  frm.def(
      "ned2eciDcm",
      py::overload_cast<
          Eigen::Ref<Eigen::Matrix3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&ned2eciDcm<double>),
      py::arg("C_n_i"),
      py::arg("lla"),
      py::arg("dt"),
      R"pbdoc(
      ned2ecuDcm
      ==========

      North-East-Down to Earth-Centered-Inertial direction cosine matrix

      Parameters
      ----------

      C_n_i : np.ndarray

          3x3 NED->ECI direction cosine matrix

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]

      dt : double 

          time elapsed between frames [s]
      )pbdoc");
  frm.def(
      "ned2eciDcm",
      py::overload_cast<const Eigen::Ref<const Eigen::Vector3d> &, const double &>(
          &ned2eciDcm<double>),
      py::arg("lla"),
      py::arg("dt"),
      R"pbdoc(
      ned2eciDcm
      ==========

      North-East-Down to Earth-Centered-Inertial direction cosine matrix

      Parameters
      ----------

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]

      dt : double 

          time elapsed between frames [s]

      Returns
      -------

      C_n_i : np.ndarray

          3x3 NED->ECI direction cosine matrix
      )pbdoc");

  // ned2ecefDcm
  frm.def(
      "ned2ecefDcm",
      py::overload_cast<Eigen::Ref<Eigen::Matrix3d>, const Eigen::Ref<const Eigen::Vector3d> &>(
          &ned2ecefDcm<double>),
      py::arg("C_e_n"),
      py::arg("lla"),
      R"pbdoc(
      ned2ecefDcm
      ===========

      North-East-Down to Earth-Centered-Earth-Fixed direction cosine matrix

      Parameters
      ----------

      C_e_n : np.ndarray

          3x3 NED->ECEF direction cosine matrix

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]
      )pbdoc");
  frm.def(
      "ned2ecefDcm",
      py::overload_cast<const Eigen::Ref<const Eigen::Vector3d> &>(&ned2ecefDcm<double>),
      py::arg("lla"),
      R"pbdoc(
      ned2ecefDcm
      ===========

      North-East-Down to Earth-Centered-Earth-Fixed direction cosine matrix

      Parameters
      ----------

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]

      Returns
      -------

      C_e_n : np.ndarray

          3x3 NED->ECEF direction cosine matrix
      )pbdoc");

  // ned2enuDcm
  frm.def(
      "ned2enuDcm",
      py::overload_cast<Eigen::Ref<Eigen::Matrix3d>>(&ned2enuDcm<double>),
      py::arg("C_n_n"),
      R"pbdoc(
      ned2enuDcm
      ==========

      North-East-Down to East-North-Up direction cosine matrix

      Parameters
      ----------

      C_n_n : np.ndarray

          3x3 NED->ENU direction cosine matrix
      )pbdoc");
  frm.def(
      "ned2enuDcm",
      py::overload_cast<>(&ned2enuDcm<double>),
      R"pbdoc(
      ned2enuDcm
      ==========

      North-East-Down to East-North-Up direction cosine matrix

      Returns
      -------

      C_n_n : np.ndarray

          3x3 NED->ENU direction cosine matrix
      )pbdoc");

  // enu2eciDcm
  frm.def(
      "enu2eciDcm",
      py::overload_cast<
          Eigen::Ref<Eigen::Matrix3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&enu2eciDcm<double>),
      py::arg("C_n_i"),
      py::arg("lla"),
      py::arg("dt"),
      R"pbdoc(
      enu2eciDcm
      ==========

      East-North-Up to Earth-Centered-Inertial direction cosine matrix

      Parameters
      ----------

      C_n_i : np.ndarray

          3x3 ENU->ECI direction cosine matrix

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]

      dt : double 

          time elapsed between frames [s]
      )pbdoc");
  frm.def(
      "enu2eciDcm",
      py::overload_cast<const Eigen::Ref<const Eigen::Vector3d> &, const double &>(
          &enu2eciDcm<double>),
      py::arg("lla"),
      py::arg("dt"),
      R"pbdoc(
      enu2eciDcm
      ==========

      East-North-Up to Earth-Centered-Inertial direction cosine matrix

      Parameters
      ----------

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]

      dt : double 

          time elapsed between frames [s]

      Returns
      -------

      C_n_i : np.ndarray

          3x3 ENU->ECI direction cosine matrix
      )pbdoc");

  // enu2ecefDcm
  frm.def(
      "enu2ecefDcm",
      py::overload_cast<Eigen::Ref<Eigen::Matrix3d>, const Eigen::Ref<const Eigen::Vector3d> &>(
          &enu2ecefDcm<double>),
      py::arg("C_e_n"),
      py::arg("lla"),
      R"pbdoc(
      enu2ecefDcm
      ===========

      East-North-Up to Earth-Centered-Earth-Fixed direction cosine matrix

      Parameters
      ----------

      C_e_n : np.ndarray

          3x3 ENU->ECEF direction cosine matrix

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]
      )pbdoc");
  frm.def(
      "enu2ecefDcm",
      py::overload_cast<const Eigen::Ref<const Eigen::Vector3d> &>(&enu2ecefDcm<double>),
      py::arg("lla"),
      R"pbdoc(
      enu2ecefDcm
      ===========

      East-North-Up to Earth-Centered-Earth-Fixed direction cosine matrix

      Parameters
      ----------

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]

      Returns
      -------

      C_e_n : np.ndarray

          3x3 ENU->ECEF direction cosine matrix
      )pbdoc");

  // enu2nedDcm
  frm.def(
      "enu2nedDcm",
      py::overload_cast<Eigen::Ref<Eigen::Matrix3d>>(&enu2nedDcm<double>),
      py::arg("C_n_n"),
      R"pbdoc(
      enu2nedDcm
      ==========

      East-North-Up to North-East-Down direction cosine matrix

      Parameters
      ----------

      C_n_n : np.ndarray

          3x3 ENU->NED direction cosine matrix
      )pbdoc");
  frm.def(
      "enu2nedDcm",
      py::overload_cast<>(&enu2nedDcm<double>),
      R"pbdoc(
      enu2nedDcm
      ==========

      East-North-Up to North-East-Down direction cosine matrix

      Returns
      -------

      C_n_n : np.ndarray

          3x3 ENU->NED direction cosine matrix
      )pbdoc");

  //*-----------------------------------------------------------------------------------------------

  // lla2eci
  frm.def(
      "lla2eci",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&lla2eci<double>),
      py::arg("eci"),
      py::arg("lla"),
      py::arg("dt"),
      R"pbdoc(
      lla2eci
      =======

      Latitude-Longitude-Height to Earth-Centered-Inertial position coordinates

      Parameters
      ----------

      eci : np.ndarray

          3x1 ECI Position [m]

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]

      dt : double 

          time elapsed between frames [s]
      )pbdoc");
  frm.def(
      "lla2eci",
      py::overload_cast<const Eigen::Ref<const Eigen::Vector3d> &, const double &>(
          &lla2eci<double>),
      py::arg("lla"),
      py::arg("dt"),
      R"pbdoc(
      lla2eci
      =======

      Latitude-Longitude-Height to Earth-Centered-Inertial position coordinates

      Parameters
      ----------

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]

      dt : double 

          time elapsed between frames [s]

      Returns
      -------

      eci : np.ndarray

          3x1 ECI Position [m]
      )pbdoc");

  // lla2ecef
  frm.def(
      "lla2ecef",
      py::overload_cast<Eigen::Ref<Eigen::Vector3d>, const Eigen::Ref<const Eigen::Vector3d> &>(
          &lla2ecef<double>),
      py::arg("xyz"),
      py::arg("lla"),
      R"pbdoc(
      lla2ecef
      ========

      Latitude-Longitude-Height to Earth-Centered-Earth-Fixed position coordinates

      Parameters
      ----------

      xyz : np.ndarray

          3x1 ECEF Position [m]

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]
      )pbdoc");
  frm.def(
      "lla2ecef",
      py::overload_cast<const Eigen::Ref<const Eigen::Vector3d> &>(&lla2ecef<double>),
      py::arg("lla"),
      R"pbdoc(
      lla2ecef
      ========

      Latitude-Longitude-Height to Earth-Centered-Earth-Fixed position coordinates

      Parameters
      ----------

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]

      Returns
      -------

      xyz : np.ndarray

          3x1 ECEF Position [m]
      )pbdoc");

  // lla2ned
  frm.def(
      "lla2ned",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&lla2ned<double>),
      py::arg("ned"),
      py::arg("lla"),
      py::arg("lla0"),
      R"pbdoc(
      lla2ned
      =======

      Latitude-Longitude-Height to North-East-Down position coordinates

      Parameters
      ----------

      ned : np.ndarray

          3x1 NED Position [m]

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height [rad, rad, m]
      )pbdoc");
  frm.def(
      "lla2ned",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&lla2ned<double>),
      py::arg("lla"),
      py::arg("lla0"),
      R"pbdoc(
      lla2ned
      =======

      Latitude-Longitude-Height to North-East-Down position coordinates

      Parameters
      ----------

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height [rad, rad, m]

      Returns
      -------

      ned : np.ndarray

          3x1 NED Position [m]
      )pbdoc");

  // lla2enu
  frm.def(
      "lla2enu",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&lla2enu<double>),
      py::arg("enu"),
      py::arg("lla"),
      py::arg("lla0"),
      R"pbdoc(
      lla2enu
      =======

      Latitude-Longitude-Height to East-North-Up position coordinates

      Parameters
      ----------

      enu : np.ndarray

          3x1 ENU Position [m]

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height [rad, rad, m]
      )pbdoc");
  frm.def(
      "lla2enu",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&lla2enu<double>),
      py::arg("lla"),
      py::arg("lla0"),
      R"pbdoc(
      lla2enu
      =======

      Latitude-Longitude-Height to East-North-Up position coordinates

      Parameters
      ----------

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height [rad, rad, m]

      Returns
      -------

      enu : np.ndarray

          3x1 ENU Position [m]
      )pbdoc");

  // lla2aer
  frm.def(
      "lla2aer",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&lla2aer<double>),
      py::arg("aer"),
      py::arg("llaR"),
      py::arg("llaT"),
      R"pbdoc(
      lla2aer
      =======

      Latitude-Longitude-Height to Azimuth-Elevation-Range position coordinates

      Parameters
      ----------

      aer : np.ndarray

          3x1 AER position [rad, rad, m]

      llaR : np.ndarray

          3x1 Reference Latitude, Longitude, Height [rad, rad, m]

      llaT : np.ndarray

          3x1 Target Latitude, Longitude, Height [rad, rad, m]
      )pbdoc");
  frm.def(
      "lla2aer",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&lla2aer<double>),
      py::arg("llaR"),
      py::arg("llaT"),
      R"pbdoc(
      lla2aer
      =======

      Latitude-Longitude-Height to Azimuth-Elevation-Range position coordinates

      Parameters
      ----------

      llaR : np.ndarray

          3x1 Reference Latitude, Longitude, Height [rad, rad, m]

      llaT : np.ndarray

          3x1 Target Latitude, Longitude, Height [rad, rad, m]

      Returns
      -------

      aer : np.ndarray

          3x1 AER position [rad, rad, m]
      )pbdoc");

  // eci2ecef
  frm.def(
      "eci2ecef",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2ecef<double>),
      py::arg("eci"),
      py::arg("xyz"),
      py::arg("dt"),
      R"pbdoc(
      eci2ecef
      ========

      Earth-Centered-Inertial to Earth-Centered-Earth-Fixed position coordinates

      Parameters
      ----------

      eci : np.ndarray

          3x1 ECI position [m]

      xyz : np.ndarray

          3x1 ECEF position [m]

      dt : double

          time elapsed between frames [s]
      )pbdoc");
  frm.def(
      "eci2ecef",
      py::overload_cast<const Eigen::Ref<const Eigen::Vector3d> &, const double &>(
          &eci2ecef<double>),
      py::arg("xyz"),
      py::arg("dt"),
      R"pbdoc(
      eci2ecef
      ========

      Earth-Centered-Inertial to Earth-Centered-Earth-Fixed position coordinates

      Parameters
      ----------

      xyz : np.ndarray

          3x1 ECEF position [m]

      dt : double

          time elapsed between frames [s]

      Returns
      -------

      eci : np.ndarray

          3x1 ECI position [m]
      )pbdoc");

  // eci2lla
  frm.def(
      "eci2lla",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2lla<double>),
      py::arg("lla"),
      py::arg("eci"),
      py::arg("dt"),
      R"pbdoc(
      eci2lla
      =======

      Earth-Centered-Inertial to Latitude-Longitude-Height position coordinates

      Parameters
      ----------

      lla : np.ndarray
          3x1 Latitude, Longitude, Height [rad, rad, m]

      eci : np.ndarray

          3x1 ECI Position [m]

      dt : double 

          time elapsed between frames [s]
      )pbdoc");
  frm.def(
      "eci2lla",
      py::overload_cast<const Eigen::Ref<const Eigen::Vector3d> &, const double &>(
          &eci2lla<double>),
      py::arg("eci"),
      py::arg("dt"),
      R"pbdoc(
      eci2lla
      =======

      Earth-Centered-Inertial to Latitude-Longitude-Height position coordinates

      Parameters
      ----------

      eci : np.ndarray

          3x1 ECI Position [m]

      dt : double 

          time elapsed between frames [s]

      Returns
      -------

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]
      )pbdoc");

  // eci2ned
  frm.def(
      "eci2ned",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2ned<double>),
      py::arg("ned"),
      py::arg("eci"),
      py::arg("lla0"),
      py::arg("dt"),
      R"pbdoc(
      eci2ned
      =======

      Earth-Centered-Inertial to North-East-Down position coordinates

      Parameters
      ----------

      ned : np.ndarray

          3x1 NED Position [m]

      eci : np.ndarray

          3x1 ECI Position [m]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height [rad, rad, m]

      dt : double 

          time elapsed between frames [s]
      )pbdoc");
  frm.def(
      "eci2ned",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2ned<double>),
      py::arg("eci"),
      py::arg("lla0"),
      py::arg("dt"),
      R"pbdoc(
      eci2ned
      =======

      Earth-Centered-Inertial to North-East-Down position coordinates

      Parameters
      ----------

      eci : np.ndarray

          3x1 ECI Position [m]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height [rad, rad, m]

      dt : double 

          time elapsed between frames [s]

      Returns
      -------

      ned : np.ndarray

          3x1 NED Position [m]
      )pbdoc");

  // eci2enu
  frm.def(
      "eci2enu",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2enu<double>),
      py::arg("enu"),
      py::arg("eci"),
      py::arg("lla0"),
      py::arg("dt"),
      R"pbdoc(
      eci2enu
      =======

      Earth-Centered-Inertial to East-North-Up position coordinates

      Parameters
      ----------

      enu : np.ndarray

          3x1 ENU Position [m]

      eci : np.ndarray

          3x1 ECI Position [m]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height [rad, rad, m]

      dt : double 

          time elapsed between frames [s]
      )pbdoc");
  frm.def(
      "eci2enu",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2enu<double>),
      py::arg("eci"),
      py::arg("lla0"),
      py::arg("dt"),
      R"pbdoc(
      eci2enu
      =======

      Earth-Centered-Inertial to East-North-Up position coordinates

      Parameters
      ----------

      eci : np.ndarray

          3x1 ECI Position [m]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height [rad, rad, m]

      dt : double 

          time elapsed between frames [s]

      Returns
      -------

      enu : np.ndarray

          3x1 ENU Position [m]
      )pbdoc");

  // eci2aer
  frm.def(
      "eci2aer",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2aer<double>),
      py::arg("aer"),
      py::arg("eciR"),
      py::arg("eciT"),
      py::arg("dt"),
      R"pbdoc(
      eci2aer
      =======

      Earth-Centered-Inertial to Azimuth-Elevation-Range position coordinates

      Parameters
      ----------

      aer : np.ndarray

          3x1 AER position [rad, rad, m]

      eciR : np.ndarray

          3x1 Reference ECI Position [m]

      eciT : np.ndarray

          3x1 Target ECI Position [m]

      dt : double 

          time elapsed between frames [s]
      )pbdoc");
  frm.def(
      "eci2aer",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2aer<double>),
      py::arg("eciR"),
      py::arg("eciT"),
      py::arg("dt"),
      R"pbdoc(
      eci2aer
      =======

      Earth-Centered-Inertial to Azimuth-Elevation-Range position coordinates

      Parameters
      ----------

      eciR : np.ndarray

          3x1 Reference ECI Position [m]

      eciT : np.ndarray

          3x1 Target ECI Position [m]

      dt : double 

          time elapsed between frames [s]

      Returns
      -------

      aer : np.ndarray

          3x1 AER position [rad, rad, m]
      )pbdoc");

  // ecef2eci
  frm.def(
      "ecef2eci",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&ecef2eci<double>),
      py::arg("xyz"),
      py::arg("eci"),
      py::arg("dt"),
      R"pbdoc(
      ecef2eci
      ========

      Earth-Centered-Earth-Fixed to Earth-Centered-Inertial position coordinates

      Parameters
      ----------

      xyz : np.ndarray

          3x1 ECEF position [m]

      eci : np.ndarray

          3x1 ECI position [m]

      dt : double

          time elapsed between frames [s]
      )pbdoc");
  frm.def(
      "ecef2eci",
      py::overload_cast<const Eigen::Ref<const Eigen::Vector3d> &, const double &>(
          &ecef2eci<double>),
      py::arg("eci"),
      py::arg("dt"),
      R"pbdoc(
      ecef2eci
      ========
      
      Earth-Centered-Earth-Fixed to Earth-Centered-Inertial position coordinates

      Parameters
      ----------

      eci : np.ndarray
          3x1 ECI position [m]

      dt : double

          time elapsed between frames [s]

      Returns
      -------

      xyz : np.ndarray

          3x1 ECEF position [m]
      )pbdoc");

  // ecef2lla
  frm.def(
      "ecef2lla",
      py::overload_cast<Eigen::Ref<Eigen::Vector3d>, const Eigen::Ref<const Eigen::Vector3d> &>(
          &ecef2lla<double>),
      py::arg("lla"),
      py::arg("xyz"),
      R"pbdoc(
      ecef2lla
      ========
      
      Earth-Centered-Earth-Fixed to Latitude-Longitude-Height position coordinates

      Parameters
      ----------

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]

      xyz : np.ndarray

          3x1 ECEF Position [m]
      )pbdoc");
  frm.def(
      "ecef2lla",
      py::overload_cast<const Eigen::Ref<const Eigen::Vector3d> &>(&ecef2lla<double>),
      py::arg("xyz"),
      R"pbdoc(
      ecef2lla
      ========

      Earth-Centered-Earth-Fixed to Latitude-Longitude-Height position coordinates

      Parameters
      ----------

      xyz : np.ndarray

          3x1 ECEF Position [m]

      Returns
      -------

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]
      )pbdoc");

  // ecef2ned
  frm.def(
      "ecef2ned",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ecef2ned<double>),
      py::arg("ned"),
      py::arg("xyz"),
      py::arg("lla0"),
      R"pbdoc(
      ecef2ned
      ========

      Earth-Centered-Earth-Fixed to North-East-Down position coordinates

      Parameters
      ----------

      xyz : np.ndarray

          3x1 ECEF Position [m]

      ned : np.ndarray

          3x1 NED Position [m]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height Position [rad, rad, m]
      )pbdoc");
  frm.def(
      "ecef2ned",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ecef2ned<double>),
      py::arg("xyz"),
      py::arg("lla0"),
      R"pbdoc(
      ecef2ned
      ========

      Earth-Centered-Earth-Fixed to North-East-Down position coordinates

      Parameters
      ----------

      xyz : np.ndarray

          3x1 ECEF Position [m]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height Position [rad, rad, m]

      Returns
      -------

      ned : np.ndarray

          3x1 NED Position [m]
      )pbdoc");

  // ecef2enu
  frm.def(
      "ecef2enu",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ecef2enu<double>),
      py::arg("enu"),
      py::arg("xyz"),
      py::arg("lla0"),
      R"pbdoc(
      ecef2enu
      ========

      Earth-Centered-Earth-Fixed to East-North-Up position coordinates

      Parameters
      ----------

      enu : np.ndarray

          3x1 ENU Position [m]

      xyz : np.ndarray

          3x1 ECEF Position [m]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height Position [rad, rad, m]
      )pbdoc");
  frm.def(
      "ecef2enu",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ecef2enu<double>),
      py::arg("xyz"),
      py::arg("lla0"),
      R"pbdoc(
      ecef2enu
      ========

      Earth-Centered-Earth-Fixed to East-North-Up position coordinates

      Parameters
      ----------

      xyz : np.ndarray

          3x1 ECEF Position [m]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height Position [rad, rad, m]

      Returns
      -------

      enu : np.ndarray

          3x1 ENU Position [m]
      )pbdoc");

  // ecef2aer
  frm.def(
      "ecef2aer",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ecef2aer<double>),
      py::arg("aer"),
      py::arg("xyzR"),
      py::arg("xyzT"),
      R"pbdoc(
      ecef2aer
      ========

      Earth-Centered-Earth-Fixed to Azimuth-Elevation-Range position coordinates

      Parameters
      ----------

      aer : np.ndarray

          3x1 AER position [rad, rad, m]

      xyzR : np.ndarray

          3x1 Reference ECEF Position [m]

      xyzT : np.ndarray

          3x1 Target ECEF Position [m]
      )pbdoc");
  frm.def(
      "ecef2aer",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ecef2aer<double>),
      py::arg("xyzR"),
      py::arg("xyzT"),
      R"pbdoc(
      ecef2aer
      ========

      Earth-Centered-Earth-Fixed to Azimuth-Elevation-Range position coordinates

      Parameters
      ----------

      xyzR : np.ndarray

          3x1 Reference ECEF Position [m]

      xyzT : np.ndarray

          3x1 Target ECEF Position [m]

      Returns
      -------

      aer : np.ndarray

          3x1 AER position [rad, rad, m]
      )pbdoc");

  // ned2eci
  frm.def(
      "ned2eci",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&ned2eci<double>),
      py::arg("eci"),
      py::arg("ned"),
      py::arg("lla0"),
      py::arg("dt"),
      R"pbdoc(
      ned2eci
      =======

      North-East-Down to Earth-Centered-Inertial position coordinates

      Parameters
      ----------

      eci : np.ndarray

          3x1 ECI Position [m]

      ned : np.ndarray

          3x1 NED Position [m]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height [rad, rad, m]

      dt : double 

          time elapsed between frames [s]
      )pbdoc");
  frm.def(
      "ned2eci",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&ned2eci<double>),
      py::arg("ned"),
      py::arg("lla0"),
      py::arg("dt"),
      R"pbdoc(
      ned2eci
      =======

      North-East-Down to Earth-Centered-Inertial position coordinates

      Parameters
      ----------

      ned : np.ndarray

          3x1 NED Position [m]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height [rad, rad, m]

      dt : double 

          time elapsed between frames [s]

      Returns
      -------

      eci : np.ndarray

          3x1 ECI Position [m]
      )pbdoc");

  // ned2ecef
  frm.def(
      "ned2ecef",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ned2ecef<double>),
      py::arg("xyz"),
      py::arg("ned"),
      py::arg("lla0"),
      R"pbdoc(
      ned2ecef
      ========

      North-East-Down to Earth-Centered-Earth-Fixed position coordinates

      Parameters
      ----------

      xyz : np.ndarray

          3x1 ECEF Position [m]

      ned : np.ndarray

          3x1 NED Position [m]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height Position [rad, rad, m]
      )pbdoc");
  frm.def(
      "ned2ecef",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ned2ecef<double>),
      py::arg("ned"),
      py::arg("lla0"),
      R"pbdoc(
      ned2ecef
      ========

      North-East-Down to Earth-Centered-Earth-Fixed position coordinates

      Parameters
      ----------

      ned : np.ndarray

          3x1 NED Position [m]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height Position [rad, rad, m]

      Returns
      -------

      xyz : np.ndarray

          3x1 ECEF Position [m]
      )pbdoc");

  // ned2lla
  frm.def(
      "ned2lla",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ned2lla<double>),
      py::arg("lla"),
      py::arg("ned"),
      py::arg("lla0"),
      R"pbdoc(
      ned2lla
      =======

      North-East-Down to Latitude-Longitude-Height position coordinates

      Parameters
      ----------

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]

      ned : np.ndarray

          3x1 NED Position [m]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height [rad, rad, m]
      )pbdoc");
  frm.def(
      "ned2lla",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ned2lla<double>),
      py::arg("ned"),
      py::arg("lla0"),
      R"pbdoc(
      ned2lla
      =======

      North-East-Down to Latitude-Longitude-Height position coordinates

      Parameters
      ----------

      ned : np.ndarray

          3x1 NED Position [m]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height [rad, rad, m]

      Returns
      -------

      lla : np.ndarray
      
          3x1 Latitude, Longitude, Height [rad, rad, m]
      )pbdoc");

  // ned2aer
  frm.def(
      "ned2aer",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ned2aer<double>),
      py::arg("aer"),
      py::arg("nedR"),
      py::arg("nedT"),
      R"pbdoc(
      ned2aer
      =======

      North-East-Down to Azimuth-Elevation-Range position coordinates

      Parameters
      ----------

      aer : np.ndarray

          3x1 AER Position [rad, rad, m]

      nedR : np.ndarray

          3x1 Reference NED Position [m]

      nedT : np.ndarray

          3x1 Target NED position [m]
      )pbdoc");
  frm.def(
      "ned2aer",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ned2aer<double>),
      py::arg("nedR"),
      py::arg("nedT"),
      R"pbdoc(
      ned2aer
      =======

      North-East-Down to Azimuth-Elevation-Range position coordinates

      Parameters
      ----------

      nedR : np.ndarray

          3x1 Reference NED Position [m]

      nedT : np.ndarray

          3x1 Target NED position [m]

      Returns
      -------

      aer : np.ndarray

          3x1 AER Position [rad, rad, m]
      )pbdoc");

  // enu2eci
  frm.def(
      "enu2eci",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&enu2eci<double>),
      py::arg("eci"),
      py::arg("enu"),
      py::arg("lla0"),
      py::arg("dt"),
      R"pbdoc(
      enu2eci
      =======

      East-North-Up to Earth-Centered-Inertial position coordinates

      Parameters
      ----------

      eci : np.ndarray

          3x1 ECI Position [m]

      enu : np.ndarray

          3x1 ENU Position [m]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height [rad, rad, m]

      dt : double 

          time elapsed between frames [s]
      )pbdoc");
  frm.def(
      "enu2eci",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&enu2eci<double>),
      py::arg("enu"),
      py::arg("lla0"),
      py::arg("dt"),
      R"pbdoc(
      enu2eci
      =======

      East-North-Up to Earth-Centered-Inertial position coordinates

      Parameters
      ----------

      enu : np.ndarray

          3x1 ENU Position [m]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height [rad, rad, m]

      dt : double 

          time elapsed between frames [s]

      Returns
      -------

      eci : np.ndarray

          3x1 ECI Position [m]
      )pbdoc");

  // enu2ecef
  frm.def(
      "enu2ecef",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&enu2ecef<double>),
      py::arg("xyz"),
      py::arg("enu"),
      py::arg("lla0"),
      R"pbdoc(
      enu2ecef
      ========

      East-North-Up to Earth-Centered-Earth-Fixed position coordinates

      Parameters
      ----------

      xyz : np.ndarray

          3x1 ECEF Position [m]

      enu : np.ndarray

          3x1 ENU Position [m]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height Position [rad, rad, m]
      )pbdoc");
  frm.def(
      "enu2ecef",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&enu2ecef<double>),
      py::arg("enu"),
      py::arg("lla0"),
      R"pbdoc(
      enu2ecef
      ========

      East-North-Up to Earth-Centered-Earth-Fixed position coordinates

      Parameters
      ----------

      enu : np.ndarray

          3x1 ENU Position [m]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height Position [rad, rad, m]

      Returns
      -------

      xyz : np.ndarray

          3x1 ECEF Position [m]
      )pbdoc");

  // enu2lla
  frm.def(
      "enu2lla",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&enu2lla<double>),
      py::arg("lla"),
      py::arg("enu"),
      py::arg("lla0"),
      R"pbdoc(
      enu2lla
      =======

      East-North-Up to Latitude-Longitude-Height position coordinates

      Parameters
      ----------

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]

      enu : np.ndarray

          3x1 ENU Position [m]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height [rad, rad, m]
      )pbdoc");
  frm.def(
      "enu2lla",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&enu2lla<double>),
      py::arg("enu"),
      py::arg("lla0"),
      R"pbdoc(
      enu2lla
      =======

      East-North-Up to Latitude-Longitude-Height position coordinates

      Parameters
      ----------

      enu : np.ndarray

          3x1 ENU Position [m]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height [rad, rad, m]

      Returns
      -------

      lla : np.ndarray

          3x1 Latitude, Longitude, Height [rad, rad, m]
      )pbdoc");

  // enu2aer
  frm.def(
      "enu2aer",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&enu2aer<double>),
      py::arg("aer"),
      py::arg("enuR"),
      py::arg("enuT"),
      R"pbdoc(
      enu2aer
      =======

      East-North-Up to Azimuth-Elevation-Range position coordinates

      Parameters
      ----------

      aer : np.ndarray

          3x1 AER Position [rad, rad, m]

      enuR : np.ndarray

          3x1 Reference ENU Position [m]

      enuT : np.ndarray

          3x1 Target ENU position [m]
      )pbdoc");
  frm.def(
      "enu2aer",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&enu2aer<double>),
      py::arg("enuR"),
      py::arg("enuT"),
      R"pbdoc(
      enu2aer
      =======

      East-North-Up to Azimuth-Elevation-Range position coordinates

      Parameters
      ----------

      enuR : np.ndarray

          3x1 Reference ENU Position [m]

      enuT : np.ndarray

          3x1 Target ENU position [m]

      Returns
      -------

      aer : np.ndarray

          3x1 AER Position [rad, rad, m]
      )pbdoc");

  //*-----------------------------------------------------------------------------------------------

  // eci2ecefv
  frm.def(
      "eci2ecefv",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2ecefv<double>),
      py::arg("v_eb_e"),
      py::arg("r_ib_i"),
      py::arg("v_ib_i"),
      py::arg("dt"),
      R"pbdoc(
      eci2ecefv
      =========

      Converts Earth-Centered-Inertial to Earth-Centered-Earth-Fixed velocity

      Parameters
      ----------

      v_eb_e : np.ndarray

          3x1 ECEF velocity [m/s]

      r_ib_i : np.ndarray

          3x1 ECI position [m]

      v_ib_i : np.ndarray

          3x1 ECI velocity [m/s]

      dt : double

          time elapsed between frames [s]
      )pbdoc");
  frm.def(
      "eci2ecefv",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2ecefv<double>),
      py::arg("r_ib_i"),
      py::arg("v_ib_i"),
      py::arg("dt"),
      R"pbdoc(
      eci2ecefv
      =========

      Converts Earth-Centered-Inertial to Earth-Centered-Earth-Fixed velocity

      Parameters
      ----------

      r_ib_i : np.ndarray

          3x1 ECI position [m]

      v_ib_i : np.ndarray

          3x1 ECI velocity [m/s]

      dt : double

          time elapsed between frames [s]

      Returns
      -------

      v_eb_e : np.ndarray

          3x1 ECEF velocity [m/s]
      )pbdoc");

  // eci2nedv
  frm.def(
      "eci2nedv",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2nedv<double>),
      py::arg("v_bn_e"),
      py::arg("r_ib_i"),
      py::arg("v_ib_i"),
      py::arg("lla0"),
      py::arg("dt"),
      R"pbdoc(
      eci2nedv
      ========

      Converts Earth-Centered-Inertial to North-East-Down velocity

      Parameters
      ----------

      v_bn_e : np.ndarray

          3x1 NED velocity [m/s]

      r_ib_i : np.ndarray

          3x1 ECI position [m]

      v_ib_i : np.ndarray

          3x1 ECI velocity [m/s]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height [rad, rad, m]

      dt : double

          time elapsed between frames [s]
      )pbdoc");
  frm.def(
      "eci2nedv",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2nedv<double>),
      py::arg("r_ib_i"),
      py::arg("v_ib_i"),
      py::arg("lla0"),
      py::arg("dt"),
      R"pbdoc(
      eci2nedv
      ========

      Converts Earth-Centered-Inertial to North-East-Down velocity

      Parameters
      ----------

      r_ib_i : np.ndarray

          3x1 ECI position [m]

      v_ib_i : np.ndarray

          3x1 ECI velocity [m/s]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height [rad, rad, m]

      dt : double

          time elapsed between frames [s]

      Returns
      -------

      v_bn_e : np.ndarray
  3x1 NED velocity [m/s]
      )pbdoc");

  // eci2enuv
  frm.def(
      "eci2enuv",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2enuv<double>),
      py::arg("v_bn_e"),
      py::arg("r_ib_i"),
      py::arg("v_ib_i"),
      py::arg("lla0"),
      py::arg("dt"),
      R"pbdoc(
      eci2enuv
      ========

      Converts Earth-Centered-Inertial to East-North-Up velocity

      Parameters
      ----------

      v_bn_e : np.ndarray

          3x1 ENU velocity [m/s]

      r_ib_i : np.ndarray

          3x1 ECI position [m]

      v_ib_i : np.ndarray

          3x1 ECI velocity [m/s]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height [rad, rad, m]

      dt : double

          time elapsed between frames [s]
      )pbdoc");
  frm.def(
      "eci2enuv",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2enuv<double>),
      py::arg("r_ib_i"),
      py::arg("v_ib_i"),
      py::arg("lla0"),
      py::arg("dt"),
      R"pbdoc(
      eci2enuv
      ========

      Converts Earth-Centered-Inertial to East-North-Up velocity

      Parameters
      ----------

      r_ib_i : np.ndarray

          3x1 ECI position [m]

      v_ib_i : np.ndarray

          3x1 ECI velocity [m/s]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height [rad, rad, m]

      dt : double

          time elapsed between frames [s]

      Returns
      -------

      v_bn_e : np.ndarray

          3x1 ENU velocity [m/s]
      )pbdoc");

  // ecef2eciv
  frm.def(
      "ecef2eciv",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&ecef2eciv<double>),
      py::arg("v_ib_i"),
      py::arg("r_eb_e"),
      py::arg("v_eb_e"),
      py::arg("dt"),
      R"pbdoc(
      ecef2eciv
      =========

      Converts Earth-Centered-Earth-Fixed to Earth-Centered-Inertial velocity

      Parameters
      ----------

      v_ib_i : np.ndarray

          3x1 ECI velocity [m/s]

      r_eb_e : np.ndarray

          3x1 ECEF position [m]

      v_eb_e : np.ndarray

          3x1 ECEF velocity [m/s]

      dt : double

          time elapsed between frames [s]
      )pbdoc");
  frm.def(
      "ecef2eciv",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&ecef2eciv<double>),
      py::arg("r_eb_e"),
      py::arg("v_eb_e"),
      py::arg("dt"),
      R"pbdoc(
      ecef2eciv
      =========

      Converts Earth-Centered-Earth-Fixed to Earth-Centered-Inertial velocity

      Parameters
      ----------

      r_eb_e : np.ndarray

          3x1 ECEF position [m]

      v_eb_e : np.ndarray

          3x1 ECEF velocity [m/s]

      dt : double

          time elapsed between frames [s]

      Returns
      -------

      v_ib_i : np.ndarray

          3x1 ECI velocity [m/s]
      )pbdoc");

  // ecef2nedv
  frm.def(
      "ecef2nedv",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ecef2nedv<double>),
      py::arg("v_nb_e"),
      py::arg("r_eb_e"),
      py::arg("v_eb_e"),
      R"pbdoc(
      ecef2nedv
      =========

      Converts Earth-Centered-Earth-Fixed to North-East-Down velocity

      Parameters
      ----------

      v_nb_e : np.ndarray

          3x1 NED velocity [m/s]

      r_eb_e : np.ndarray

          3x1 ECEF position [m]

      v_eb_e : np.ndarray

          3x1 ECEF velocity [m/s]
      )pbdoc");
  frm.def(
      "ecef2nedv",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ecef2nedv<double>),
      py::arg("r_eb_e"),
      py::arg("v_eb_e"),
      R"pbdoc(
      ecef2nedv
      =========

      Converts Earth-Centered-Earth-Fixed to North-East-Down velocity

      Parameters
      ----------

      r_eb_e : np.ndarray

          3x1 ECEF position [m]

      v_eb_e : np.ndarray

          3x1 ECEF velocity [m/s]

      Returns
      -------

      v_nb_e : np.ndarray

          3x1 NED velocity [m/s]
      )pbdoc");

  // ecef2enuv
  frm.def(
      "ecef2enuv",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ecef2enuv<double>),
      py::arg("v_nb_e"),
      py::arg("r_eb_e"),
      py::arg("v_eb_e"),
      R"pbdoc(
      ecef2enuv
      =========
      Converts Earth-Centered-Earth-Fixed to East-North-Up velocity

      Parameters
      ----------

      v_nb_e : np.ndarray

        3x1 ENU velocity [m/s]

      r_eb_e : np.ndarray

          3x1 ECEF position [m]

      v_eb_e : np.ndarray

          3x1 ECEF velocity [m/s]
      )pbdoc");
  frm.def(
      "ecef2enuv",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ecef2enuv<double>),
      py::arg("r_eb_e"),
      py::arg("v_eb_e"),
      R"pbdoc(
      ecef2enuv
      =========

      Converts Earth-Centered-Earth-Fixed to East-North-Up velocity

      Parameters
      ----------

      r_eb_e : np.ndarray

          3x1 ECEF position [m]

      v_eb_e : np.ndarray

          3x1 ECEF velocity [m/s]

      Returns
      -------

      v_nb_e : np.ndarray

          3x1 ENU velocity [m/s]
      )pbdoc");

  // ned2eciv
  frm.def(
      "ned2eciv",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&ned2eciv<double>),
      py::arg("v_in_i"),
      py::arg("r_nb_e"),
      py::arg("v_nb_e"),
      py::arg("lla0"),
      py::arg("dt"),
      R"pbdoc(
      ned2eciv
      ========

      Converts North-East-Down to Earth-Centered-Inertial velocity

      Parameters
      ----------

      v_ib_i : np.ndarray

          3x1 ECI velocity [m/s]

      r_nb_e : np.ndarray

          3x1 NED position [m]

      v_nb_e : np.ndarray

          3x1 NED velocity [m/s]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height [rad, rad, m]

      dt : double

          time elapsed between frames [s]
      )pbdoc");
  frm.def(
      "ned2eciv",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&ned2eciv<double>),
      py::arg("r_nb_e"),
      py::arg("v_nb_e"),
      py::arg("lla0"),
      py::arg("dt"),
      R"pbdoc(
      ned2eciv
      ========

      Converts North-East-Down to Earth-Centered-Inertial velocity

      Parameters
      ----------

      r_nb_e : np.ndarray

          3x1 NED position [m]

      v_nb_e : np.ndarray

          3x1 NED velocity [m/s]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height [rad, rad, m]

      dt : double

          time elapsed between frames [s]

      Returns
      -------

      v_ib_i : np.ndarray

          3x1 ECI velocity [m/s]
      )pbdoc");

  // ned2ecefv
  frm.def(
      "ned2ecefv",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ned2ecefv<double>),
      py::arg("v_eb_e"),
      py::arg("v_nb_e"),
      py::arg("lla0"),
      R"pbdoc(
      ned2ecefv
      =========

      Converts North-East-Down to Earth-Centered-Earth-Fixed velocity

      Parameters
      ----------

      v_eb_e : np.ndarray

          3x1 ECEF velocity [m/s]

      v_nb_e : np.ndarray

          3x1 NED velocity [m/s]

      v_nb_e : np.ndarray

          3x1 Latitude, Longitude, Altitude [rad,rad,m]
      )pbdoc");
  frm.def(
      "ned2ecefv",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ned2ecefv<double>),
      py::arg("v_nb_e"),
      py::arg("lla0"),
      R"pbdoc(
      ned2ecefv
      =========

      Converts North-East-Down to Earth-Centered-Earth-Fixed velocity

      Parameters
      ----------

      v_nb_e : np.ndarray

          3x1 NED velocity [m/s]

      v_nb_e : np.ndarray

          3x1 Latitude, Longitude, Altitude [rad,rad,m]

      Returns
      -------

      v_eb_e : np.ndarray

          3x1 ECEF velocity [m/s]
      )pbdoc");

  // enu2eciv
  frm.def(
      "enu2eciv",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&enu2eciv<double>),
      py::arg("v_in_i"),
      py::arg("r_nb_e"),
      py::arg("v_nb_e"),
      py::arg("lla0"),
      py::arg("dt"),
      R"pbdoc(
      enu2eciv
      ========

      Converts East-North-Up to Earth-Centered-Inertial velocity

      Parameters
      ----------

      v_ib_i : np.ndarray

          3x1 ECI velocity [m/s]

      r_nb_e : np.ndarray

          3x1 ENU position [m]

      v_nb_e : np.ndarray

          3x1 ENU velocity [m/s]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height [rad, rad, m]

      dt : double

          time elapsed between frames [s]
      )pbdoc");
  frm.def(
      "enu2eciv",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&enu2eciv<double>),
      py::arg("r_nb_e"),
      py::arg("v_nb_e"),
      py::arg("lla0"),
      py::arg("dt"),
      R"pbdoc(
      enu2eciv
      ========

      Converts East-North-Up to Earth-Centered-Inertial velocity

      Parameters
      ----------

      r_nb_e : np.ndarray

          3x1 ENU position [m]

      v_nb_e : np.ndarray

          3x1 ENU velocity [m/s]

      lla0 : np.ndarray

          3x1 Reference Latitude, Longitude, Height [rad, rad, m]
      dt : double

          time elapsed between frames [s]

      Returns
      -------

      v_ib_i : np.ndarray

          3x1 ECI velocity [m/s]
      )pbdoc");

  // enu2ecefv
  frm.def(
      "enu2ecefv",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&enu2ecefv<double>),
      py::arg("v_eb_e"),
      py::arg("v_nb_e"),
      py::arg("lla0"),
      R"pbdoc(
      enu2ecefv
      =========

      Converts East-North-Up to Earth-Centered-Earth-Fixed velocity

      Parameters
      ----------

      v_eb_e : np.ndarray

          3x1 ECEF velocity [m/s]

      v_nb_e : np.ndarray

          3x1 ENU velocity [m/s]

      lla0 : np.ndarray

          3x1 Latitude, Longitude, Altitude [rad,rad,m]
      )pbdoc");
  frm.def(
      "enu2ecefv",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&enu2ecefv<double>),
      py::arg("v_nb_e"),
      py::arg("lla0"),
      R"pbdoc(
      enu2ecefv
      =========

      Converts East-North-Up to Earth-Centered-Earth-Fixed velocity

      Parameters
      ----------

      v_nb_e : np.ndarray

          3x1 ENU velocity [m/s]

      lla0 : np.ndarray

          3x1 Latitude, Longitude, Altitude [rad,rad,m]

      Returns
      -------

      v_eb_e : np.ndarray

          3x1 ECEF velocity [m/s]
      )pbdoc");

  // eci2ecefw
  frm.def(
      "eci2ecefw",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2ecefw<double>));
  frm.def(
      "eci2ecefw",
      py::overload_cast<const Eigen::Ref<const Eigen::Vector3d> &, const double &>(
          &eci2ecefw<double>));

  // eci2nedw
  frm.def(
      "eci2nedw",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2nedw<double>));
  frm.def(
      "eci2nedw",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2nedw<double>));

  // eci2enuw
  frm.def(
      "eci2enuw",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2enuw<double>));
  frm.def(
      "eci2enuw",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2enuw<double>));

  // ecef2eciw
  frm.def(
      "ecef2eciw",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&ecef2eciw<double>));
  frm.def(
      "ecef2eciw",
      py::overload_cast<const Eigen::Ref<const Eigen::Vector3d> &, const double &>(
          &ecef2eciw<double>));

  // ecef2nedw
  frm.def(
      "ecef2nedw",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ecef2nedw<double>));
  frm.def(
      "ecef2nedw",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ecef2nedw<double>));

  // ecef2enuw
  frm.def(
      "ecef2enuw",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ecef2enuw<double>));
  frm.def(
      "ecef2enuw",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ecef2enuw<double>));

  // ned2eciw
  frm.def(
      "ned2eciw",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&ned2eciw<double>));
  frm.def(
      "ned2eciw",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&ned2eciw<double>));

  // ned2ecefw
  frm.def(
      "ned2ecefw",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ned2ecefw<double>));
  frm.def(
      "ned2ecefw",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ned2ecefw<double>));

  // enu2eciw
  frm.def(
      "ned2eciw",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&ned2eciw<double>));
  frm.def(
      "ned2eciw",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&ned2eciw<double>));

  // enu2ecefw
  frm.def(
      "enu2ecefw",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&enu2ecefw<double>));
  frm.def(
      "enu2ecefw",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&enu2ecefw<double>));

  // eci2ecefa
  frm.def(
      "eci2ecefa",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2ecefa<double>));
  frm.def(
      "eci2ecefa",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2ecefa<double>));

  // eci2neda
  frm.def(
      "eci2neda",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2neda<double>));
  frm.def(
      "eci2neda",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2neda<double>));

  // eci2enua
  frm.def(
      "eci2enua",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2enua<double>));
  frm.def(
      "eci2enua",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&eci2enua<double>));

  // ecef2ecia
  frm.def(
      "ecef2ecia",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&ecef2ecia<double>));
  frm.def(
      "ecef2ecia",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&ecef2ecia<double>));

  // ecef2neda
  frm.def(
      "ecef2neda",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ecef2neda<double>));
  frm.def(
      "ecef2neda",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ecef2neda<double>));

  // ecef2enua
  frm.def(
      "ecef2enua",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ecef2enua<double>));
  frm.def(
      "ecef2enua",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ecef2enua<double>));

  // ned2ecia
  frm.def(
      "ned2ecia",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&ned2ecia<double>));
  frm.def(
      "ned2ecia",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&ned2ecia<double>));

  // ned2ecefa
  frm.def(
      "ned2ecefa",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ned2ecefa<double>));
  frm.def(
      "ned2ecefa",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&ned2ecefa<double>));

  // enu2ecia
  frm.def(
      "enu2ecia",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&enu2ecia<double>));
  frm.def(
      "enu2ecia",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const double &>(&enu2ecia<double>));

  // enu2ecefa
  frm.def(
      "enu2ecefa",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&enu2ecefa<double>));
  frm.def(
      "enu2ecefa",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &>(&enu2ecefa<double>));

  //! === Math submodule ===========================================================================
  py::module_ math = h.def_submodule("math", R"pbdoc(
        Math
        ====
        
        Common mathematical operations.)pbdoc");

  // Skew
  math.def(
      "Skew",
      [](Eigen::Matrix3d &M, const Eigen::Vector3d &v) {
        Skew<Eigen::Matrix3d, Eigen::Vector3d>(M, v);
      },
      py::arg("R"),
      py::arg("v"),
      R"pbdoc(
      Skew
      ====
      
      Converts vector into its skew symmetric form

      Parameters
      ----------

      R : np.ndarray

          3x3 skew symmetric matrix

      v : np.ndarray

          3x1 vector
      )pbdoc");
  math.def(
      "Skew",
      [](const Eigen::Vector3d &v) { return Skew(v); },
      py::arg("v"),
      R"pbdoc(
        Skew
        ====

        Converts vector into its skew symmetric form

        Parameters
        ----------

        v : np.ndarray

            3x1 vector

        Returns
        -------

        R : np.ndarray

            3x3 skew symmetric matrix
        )pbdoc");

  // DeSkew
  math.def(
      "DeSkew",
      py::overload_cast<Eigen::Ref<Eigen::Vector3d>, const Eigen::Ref<const Eigen::Matrix3d> &>(
          &DeSkew<double>),
      py::arg("v"),
      py::arg("R"),
      R"pbdoc(
      DeSkew
      ======

      Converts skew symmetric matrix into its vector form

      Parameters
      ----------

      v : np.ndarray

          3x1 vector

      R : np.ndarray

          3x3 skew symmetric matrix
      )pbdoc");
  math.def(
      "DeSkew",
      py::overload_cast<const Eigen::Ref<const Eigen::Matrix3d> &>(&DeSkew<double>),
      py::arg("R"),
      R"pbdoc(
      DeSkew
      ======

      Converts skew symmetric matrix into its vector form

      Parameters
      ----------

      R : np.ndarray

          3x3 skew symmetric matrix

      Returns
      -------

      v : np.ndarray

          3x1 vector
      )pbdoc");

  // CircMod
  math.def(
      "CircMod",
      &CircMod<double>,
      py::arg("x"),
      py::arg("y"),
      R"pbdoc(
      CircMod
      =======

      Modulus of floating point number

      Parameters
      ----------

      x : double

          user input and output

      y : double

          value to take modulus about
      )pbdoc");

  // WrapTo2Pi
  math.def(
      "WrapTo2Pi",
      &CircMod2Pi<double>,
      py::arg("x"),
      R"pbdoc(
      WrapTo2Pi
      =========

      Wraps angles from [0, 2*pi]

      Parameters
      ----------

      x : double

          user input and output [rad]
      )pbdoc");

  // WrapToPi
  math.def(
      "WrapToPi",
      &WrapPiToPi<double>,
      py::arg("x"),
      R"pbdoc(
      WrapToPi
      ========

      Wraps angles from [-pi, pi]

      Parameters
      ----------

      x : double

          user input and output [rad]
      )pbdoc");

  // WrapEulerAngles
  math.def(
      "WrapEulerAngles",
      &WrapEulerAngles<double>,
      py::arg("e"),
      R"pbdoc(
      WrapEulerAngles
      ===============

      Auto wrap euler angles depending on pitch angle

      Parameters
      ----------

      e : np.ndarray

          3x1 user input and output euler angles [rad]
      )pbdoc");

  // quatmat
  math.def(
      "quatmat",
      &quatmat<double>,
      py::arg("q"),
      R"pbdoc(
      quatmat
      =======

      Converts quaternion into its 4x4 matrix view

      Parameters
      ----------

      q : np.ndarray

          4x1 quaternion

      Returns
      -------

      Q : np.ndarray

          4x4 quaternion matrix view
      )pbdoc");

  // quatdot
  math.def(
      "quatdot",
      &quatdot<double>,
      py::arg("p"),
      py::arg("q"),
      R"pbdoc(
      quatdot
      =======

      Quaternion product (r = p . q)

      Parameters
      ----------

      p : np.ndarray

          4x1 quaternion

      q : np.ndarray

          4x1 quaternion

      Returns
      -------

      r : np.ndarray

          4x1 quaternion product
      )pbdoc");

  // quatconj
  math.def(
      "quatconj",
      &quatconj<double>,
      py::arg("q"),
      R"pbdoc(
      quatconj
      ========

      Conjugates/Inverts input quaternion

      Parameters
      ----------

      q : np.ndarray

          4x1 quaternion
      )pbdoc");

  // quatinv
  math.def(
      "quatinv",
      &quatinv<double>,
      py::arg("q"),
      R"pbdoc(
      quatinv
      =======

      Inverts/Conjugates input quaternion

      Parameters
      ----------

      q : np.ndarray

          4x1 quaternion
      )pbdoc");

  // quatnorm
  math.def(
      "quatnorm",
      &quatnorm<double>,
      py::arg("q"),
      R"pbdoc(
      quatnorm
      ========

      Normalizes input quaternion

      Parameters
      ----------

      q : np.ndarray

          4x1 quaternion
      )pbdoc");

  // dcmnorm
  math.def(
      "dcmnorm",
      &dcmnorm<double>,
      py::arg("R"),
      R"pbdoc(
      dcmnorm
      =======

      Normalizes input DCM

      Parameters
      ----------

      R : np.ndarray

          3x3 DCM
      )pbdoc");

  // Rodrigues
  math.def(
      "Rodrigues",
      [](const Eigen::Vector3d &v) { return Rodrigues(v); },
      py::arg("v"),
      R"pbdoc(
      Rodrigues
      =========

      Rodrigues formula for converting a 3x1 vector into its matrix exponential form

      Parameters
      ----------

      v : np.ndarray

          3x1 vector

      Returns
      -------

      R : np.ndarray

          3x3 matrix exponential
      )pbdoc");

  // Rodrigues4
  math.def(
      "Rodrigues4",
      [](const Eigen::Vector3d &v) { return Rodrigues4(v); },
      py::arg("v"),
      R"pbdoc(
      Rodrigues4
      ==========

      4th order approximation of Rodrigues formula

      Parameters
      ----------

      v : np.ndarray

          3x1 vector

      Returns
      -------

      R : np.ndarray

          3x3 matrix exponential
      )pbdoc");

  // scalar2expm
  math.def(
      "scalar2expm",
      &scalar2expm<double>,
      py::arg("x"),
      R"pbdoc(
      scalar2expm
      ===========

      Converts scalar value into 2x2 matrix by rotating the value in radians (positive CCW)

      Parameters
      ----------

      x : double

          Rotation angle

      Returns
      -------

      R : np.ndarray

          2x2 matrix exponential
      )pbdoc");

  // vec2expm
  math.def(
      "vec2expm",
      [](const Eigen::VectorXd &v) { return vec2expm(v); },
      py::arg("v"),
      R"pbdoc(
      vec2expm
      ========

      Converts vector into its matrix exponential form

      Parameters
      ----------

      v : np.ndarray

          3x1 vector

      Returns
      -------

      R : np.ndarray

          3x3 matrix exponential
      )pbdoc");

  // expm2vec
  math.def(
      "expm2vec",
      &expm2vec<double>,
      py::arg("R"),
      R"pbdoc(
      expm2vec
      ========

      Converts matrix exponential into its vector form

      Parameters
      ----------

      R : np.ndarray

          3x3 matrix exponential

      Returns
      -------

      v : np.ndarray

          3x1 vector
      )pbdoc");

  // pow2db
  math.def(
      "pow2db",
      py::overload_cast<const double &>(&watt2db<double>),
      py::arg("x"),
      R"pbdoc(
      pow2db
      ======
      
      Convert unit of power to decibels

      Parameters
      ----------

      x : double

          Input power

      Returns
      -------

      y : double

          Output dB
      )pbdoc");

  // db2pow
  math.def(
      "db2pow",
      py::overload_cast<const double &>(&db2watt<double>),
      py::arg("x"),
      R"pbdoc(
      db2pow
      ======

      Convert decibels to unit of power

      Parameters
      ----------

      x : double

          Input dB

      Returns
      -------

      y : double

          Output power
      )pbdoc");

  // mag2db
  math.def(
      "mag2db",
      py::overload_cast<const double &>(&volt2db<double>),
      py::arg("x"),
      R"pbdoc(
      mag2db
      ======

      Convert unit of magnitude to decibels

      Parameters
      ----------

      x : double

          Input magnitude

      Returns
      -------

      y : double

          Output dB
      )pbdoc");

  // db2mag
  math.def(
      "db2mag",
      py::overload_cast<const double &>(&db2volt<double>),
      py::arg("x"),
      R"pbdoc(
      db2mag
      ======

      Convert decibels to unit of magnitude

      Parameters
      ----------

      x : double

          Input dB

      Returns
      -------

      y : double

          Output magnitude
      )pbdoc");

  //! === Earth-Models submodule ===================================================================
  py::module_ mod = h.def_submodule("models", R"pbdoc(
        Models
        ======
        
        Simple Earth models commonly used in navigation equations.)pbdoc");

  // TransverseRadius
  mod.def(
      "TransverseRadius",
      py::overload_cast<double &, const double &>(&TransverseRadius<double>),
      py::arg("Re"),
      py::arg("phi"),
      R"pbdoc(
      TransverseRadius
      ================

      Calculates the transverse radius relative to user latitude

      Parameters
      ----------

      Re : double

          Transverse radius [m]

      phi : double

          Latitude [rad]
      )pbdoc");
  mod.def(
      "TransverseRadius",
      py::overload_cast<const double &>(&TransverseRadius<double>),
      py::arg("phi"),
      R"pbdoc(
      TransverseRadius
      ================

      Calculates the transverse radius relative to user latitude

      Parameters
      ----------

      phi : double

          Latitude [rad]

      Returns
      -------

      Re : double

          Transverse radius [m]
  
      )pbdoc");

  // MeridianRadius
  mod.def(
      "MeridianRadius",
      py::overload_cast<double &, const double &>(&MeridianRadius<double>),
      py::arg("Rn"),
      py::arg("phi"),
      R"pbdoc(
      MeridianRadius
      ==============

      Calculates the meridian radius relative to user latitude

      Parameters
      ----------

      Rn : double

          Meridian radius [m]

      phi : double

          Latitude [rad]
      )pbdoc");
  mod.def(
      "MeridianRadius",
      py::overload_cast<const double &>(&MeridianRadius<double>),
      py::arg("phi"),
      R"pbdoc(
      MeridianRadius
      ==============

      Calculates the meridian radius relative to user latitude

      Parameters
      ----------

      phi : double

          Latitude [rad]

      Returns
      -------

      Rn : double

          Meridian radius [m]
      )pbdoc");

  // GeocentricRadius
  mod.def(
      "GeocentricRadius",
      py::overload_cast<double &, const double &>(&GeocentricRadius<double>),
      py::arg("R_es_e"),
      py::arg("phi"),
      R"pbdoc(
      GeocentricRadius
      ================

      Calculates the geocentric radius relative to user latitude

      Parameters
      ----------

      R_es_e : double

          Geocentric radius [m]

      phi : double

          Latitude [rad]
      )pbdoc");
  mod.def(
      "GeocentricRadius",
      py::overload_cast<const double &>(&GeocentricRadius<double>),
      py::arg("phi"),
      R"pbdoc(
      GeocentricRadius
      ================

      Calculates the geocentric radius relative to user latitude

      Parameters
      ----------

      phi : double

          Latitude [rad]

      Returns
      -------

      R_es_e : double

          Geocentric radius [m]
      )pbdoc");

  // TransAndMerRadii
  mod.def(
      "TransAndMerRadii",
      py::overload_cast<double &, double &, const double &>(&TransAndMerRadii<double>),
      py::arg("Re"),
      py::arg("Rn"),
      py::arg("phi"),
      R"pbdoc(
      TransAndMerRadii
      ================

      Calculates the {Transverse, Meridian} radii relative to user latitude

      Parameters
      ----------

      Re : double
          Transverse radius [m]

      Rn : double

          Meridian radius [m]

      phi : double

          Latitude [rad]
      )pbdoc");

  // RadiiOfCurvature
  mod.def(
      "RadiiOfCurvature",
      py::overload_cast<double &, double &, double &, const double &>(&RadiiOfCurvature<double>),
      py::arg("Re"),
      py::arg("Rn"),
      py::arg("R_es_e"),
      py::arg("phi"),
      R"pbdoc(
      RadiiOfCurvature
      ================

      Calculates the {Transverse, Meridian, Geocentric} radii relative to user latitude

      Parameters
      ----------

      Re : double

          Transverse radius [m]

      Rn : double

          Meridian radius [m]

      R_es_e : double

          Geocentric radius [m]

      phi : double

          Latitude [rad]
      )pbdoc");

  // EarthRate
  mod.def(
      "EarthRate",
      py::overload_cast<Eigen::Ref<Eigen::Vector3d>, const double &, const bool>(
          &EarthRate<double>),
      py::arg("w_ie_n"),
      py::arg("phi"),
      py::arg("IsNed") = true,
      R"pbdoc(
      EarthRate
      =========

      Rotation rate of the earth relative to the 'NAV' frame

      Parameters
      ----------

      w_ie_n : np.ndarray

          size 3 vector of earth's rotation in the 'NAV' frame

      phi : double

          Latitude [rad]

      IsNed : bool

          Is desired frame NED, default is True 
      )pbdoc");
  mod.def(
      "EarthRate",
      py::overload_cast<const double &, const bool>(&EarthRate<double>),
      py::arg("phi"),
      py::arg("IsNed") = true,
      R"pbdoc(
      EarthRate
      =========

      Rotation rate of the earth relative to the 'NAV' frame

      Parameters
      ----------

      phi : double

          Latitude [rad]

      IsNed : bool

          Is desired frame NED, default is True 

      Returns
      -------

      w_ie_n : np.ndarray

          size 3 vector of earth's rotation in the 'NAV' frame
      )pbdoc");

  // TransportRate
  mod.def(
      "TransportRate",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const bool>(&TransportRate<double>),
      py::arg("w_en_n"),
      py::arg("phi"),
      py::arg("v_nb_e"),
      py::arg("IsNed") = true,
      R"pbdoc(
      TransportRate
      =============

      Transport rate of the 'ECEF' frame relative to the 'NAV' frame

      Parameters
      ----------

      w_en_n : np.ndarray

          size 3 vector of earth's rotation in the 'NAV' frame

      lla : np.ndarray

          size 3 vector of Latitude, Longitude, Height [rad, rad, m]

      v_nb_e : np.ndarray

          size 3 velocity vector in the 'NAV' coordinate system

      IsNed : bool

          Is desired frame NED, default is True 
      )pbdoc");
  mod.def(
      "EarthRate",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const bool>(&TransportRate<double>),
      py::arg("phi"),
      py::arg("v_nb_e"),
      py::arg("IsNed") = true,
      R"pbdoc(
      EarthRate
      =========

      Transport rate of the 'ECEF' frame relative to the 'NAV' frame

      Parameters
      ----------

      lla : np.ndarray

          size 3 vector of Latitude, Longitude, Height [rad, rad, m]

      v_nb_e : np.ndarray

          size 3 velocity vector in the 'NAV' coordinate system

      IsNed : bool

          Is desired frame NED, default is True 

      Returns
      -------

      w_en_n : np.ndarray

          size 3 vector of earth's rotation in the 'NAV' frame
      )pbdoc");

  // CoriolisRate
  mod.def(
      "CoriolisRate",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const bool>(&CoriolisRate<double>),
      py::arg("coriolis"),
      py::arg("lla"),
      py::arg("v_nb_e"),
      py::arg("IsNed") = true,
      R"pbdoc(
      CoriolisRate
      ============

      Coriolis effect perceived in the "NAV" frame

      Parameters
      ----------

      coriolis : np.ndarray

          size 3 coriolis effect

      lla : np.ndarray

          size 3 vector of Latitude, Longitude, Height [rad, rad, m]

      v_nb_e : np.ndarray

          size 3 velocity vector in the 'NAV' coordinate system

      IsNed : bool

          Is desired frame NED, default is True 
      )pbdoc");
  mod.def(
      "CoriolisRate",
      py::overload_cast<
          const Eigen::Ref<const Eigen::Vector3d> &,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const bool>(&CoriolisRate<double>),
      py::arg("lla"),
      py::arg("v_nb_e"),
      py::arg("IsNed") = true,
      R"pbdoc(
      CoriolisRate
      ============

      Coriolis effect perceived in the "NAV" frame

      Parameters
      ----------

      lla : np.ndarray

          size 3 vector of Latitude, Longitude, Height [rad, rad, m]
      
      v_nb_e : np.ndarray
  
          size 3 velocity vector in the 'NAV' coordinate system

      IsNed : bool

          Is desired frame NED, default is True 

      Returns
      -------

      coriolis : np.ndarray

          size 3 coriolis effect
      )pbdoc");

  // Somigliana
  mod.def(
      "Somigliana",
      py::overload_cast<double &, const double &>(&Somigliana<double>),
      py::arg("g0"),
      py::arg("phi"),
      R"pbdoc(
      Somigliana
      ==========

      Calculates the somilgiana model reference gravity

      Parameters
      ----------

      g0 : double

          Somigliana gravitation

      phi : double

          Latitude [rad]
      )pbdoc");
  mod.def(
      "Somigliana",
      py::overload_cast<const double &>(&Somigliana<double>),
      py::arg("phi"),
      R"pbdoc(
      Somigliana
      ==========

      Calculates the somilgiana model reference gravity

      Parameters
      ----------

      phi : double

          Latitude [rad]

      Returns
      -------

      g0 : double

          Somigliana gravitation
      )pbdoc");

  // LocalGravity
  mod.def(
      "LocalGravity",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &,
          const bool>(&LocalGravity<double>),
      py::arg("g"),
      py::arg("lla"),
      py::arg("IsNed") = true,
      R"pbdoc(
      LocalGravity
      ============

      Calculates gravity in the Local/NAV (ENU or NED) frame

      Parameters
      ----------

      g : np.ndarray

          size 3 Local/NAV frame gravity vector

      lla : np.ndarray

          size 3 vector of Latitude, Longitude, Height [rad, rad, m]

      IsNed : bool

          Is desired frame NED, default is True 
      )pbdoc");
  mod.def(
      "LocalGravity",
      py::overload_cast<const Eigen::Ref<const Eigen::Vector3d> &, const bool>(
          &LocalGravity<double>),
      py::arg("lla"),
      py::arg("IsNed") = true,
      R"pbdoc(
      LocalGravity
      ============

      Calculates gravity in the Local/NAV (ENU or NED) frame

      Parameters
      ----------

      lla : np.ndarray

          size 3 vector of Latitude, Longitude, Height [rad, rad, m]

      IsNed : bool

          Is desired frame NED, default is True 

      Returns
      -------

      g : np.ndarray

          size 3 Local/NAV frame gravity vector
      )pbdoc");

  // EcefGravity
  mod.def(
      "EcefGravity",
      py::overload_cast<
          Eigen::Ref<Eigen::Vector3d>,
          Eigen::Ref<Eigen::Vector3d>,
          const Eigen::Ref<const Eigen::Vector3d> &>(&EcefGravity<double>),
      py::arg("g"),
      py::arg("gamma"),
      py::arg("xyz"),
      R"pbdoc(
      EcefGravity
      ===========

      Calculates gravity in the Earth-Centered-Earth-Fixed frame

      Parameters
      ----------

      g : np.ndarray

          size 3 ECEF frame gravity vector

      gamma : np.ndarray

          size 3 ECEF frame gravitation

      xyz : np.ndarray

          size 3 ECEF position [m]
      )pbdoc");
}
