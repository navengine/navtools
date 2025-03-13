"""

        Math
        ====
        
        Common mathematical operations.
"""
from __future__ import annotations
import numpy
import typing
__all__ = ['CircMod', 'DeSkew', 'Rodrigues', 'Rodrigues4', 'Skew', 'WrapEulerAngles', 'WrapTo2Pi', 'WrapToPi', 'db2mag', 'db2pow', 'dcmnorm', 'expm2vec', 'mag2db', 'pow2db', 'quatconj', 'quatdot', 'quatinv', 'quatmat', 'quatnorm', 'scalar2expm', 'vec2expm']
def CircMod(x: float, y: float) -> None:
    """
          CircMod
          =======
    
          Modulus of floating point number
    
          Parameters
          ----------
    
          x : double
    
              user input and output
    
          y : double
    
              value to take modulus about
    """
@typing.overload
def DeSkew(v: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], R: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.f_contiguous]) -> None:
    """
          DeSkew
          ======
    
          Converts skew symmetric matrix into its vector form
    
          Parameters
          ----------
    
          v : np.ndarray
    
              3x1 vector
    
          R : np.ndarray
    
              3x3 skew symmetric matrix
    """
@typing.overload
def DeSkew(R: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.f_contiguous]) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
def Rodrigues(v: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 3]]:
    """
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
    """
def Rodrigues4(v: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 3]]:
    """
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
    """
@typing.overload
def Skew(R: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.writeable, numpy.ndarray.flags.f_contiguous], v: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    """
          Skew
          ====
          
          Converts vector into its skew symmetric form
    
          Parameters
          ----------
    
          R : np.ndarray
    
              3x3 skew symmetric matrix
    
          v : np.ndarray
    
              3x1 vector
    """
@typing.overload
def Skew(v: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 3]]:
    """
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
    """
def WrapEulerAngles(e: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable]) -> None:
    """
          WrapEulerAngles
          ===============
    
          Auto wrap euler angles depending on pitch angle
    
          Parameters
          ----------
    
          e : np.ndarray
    
              3x1 user input and output euler angles [rad]
    """
def WrapTo2Pi(x: float) -> None:
    """
          WrapTo2Pi
          =========
    
          Wraps angles from [0, 2*pi]
    
          Parameters
          ----------
    
          x : double
    
              user input and output [rad]
    """
def WrapToPi(x: float) -> None:
    """
          WrapToPi
          ========
    
          Wraps angles from [-pi, pi]
    
          Parameters
          ----------
    
          x : double
    
              user input and output [rad]
    """
def db2mag(x: float) -> float:
    """
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
    """
def db2pow(x: float) -> float:
    """
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
    """
def dcmnorm(R: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.writeable, numpy.ndarray.flags.f_contiguous]) -> None:
    """
          dcmnorm
          =======
    
          Normalizes input DCM
    
          Parameters
          ----------
    
          R : np.ndarray
    
              3x3 DCM
    """
def expm2vec(R: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.f_contiguous]) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
def mag2db(x: float) -> float:
    """
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
    """
def pow2db(x: float) -> float:
    """
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
    """
def quatconj(q: numpy.ndarray[numpy.float64[4, 1], numpy.ndarray.flags.writeable]) -> None:
    """
          quatconj
          ========
    
          Conjugates/Inverts input quaternion
    
          Parameters
          ----------
    
          q : np.ndarray
    
              4x1 quaternion
    """
def quatdot(p: numpy.ndarray[numpy.float64[4, 1]], q: numpy.ndarray[numpy.float64[4, 1]]) -> numpy.ndarray[numpy.float64[4, 1]]:
    """
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
    """
def quatinv(q: numpy.ndarray[numpy.float64[4, 1], numpy.ndarray.flags.writeable]) -> None:
    """
          quatinv
          =======
    
          Inverts/Conjugates input quaternion
    
          Parameters
          ----------
    
          q : np.ndarray
    
              4x1 quaternion
    """
def quatmat(q: numpy.ndarray[numpy.float64[4, 1]]) -> numpy.ndarray[numpy.float64[4, 4]]:
    """
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
    """
def quatnorm(q: numpy.ndarray[numpy.float64[4, 1], numpy.ndarray.flags.writeable]) -> None:
    """
          quatnorm
          ========
    
          Normalizes input quaternion
    
          Parameters
          ----------
    
          q : np.ndarray
    
              4x1 quaternion
    """
def scalar2expm(x: float) -> numpy.ndarray[numpy.float64[2, 2]]:
    """
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
    """
def vec2expm(v: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 3]]:
    """
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
    """
