"""

        Attitude
        ========
        
        Attitude representations and conversions between them.
"""
from __future__ import annotations
import numpy
import typing
__all__ = ['RotX', 'RotY', 'RotZ', 'dcm2euler', 'dcm2quat', 'euler2dcm', 'euler2quat', 'quat2dcm', 'quat2euler']
@typing.overload
def RotX(C: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.writeable, numpy.ndarray.flags.f_contiguous], x: float) -> None:
    """
          RotX
          ====
    
          Converts euler angle about x-axis into DCM rotation about x-axis
    
          Parameters
          ----------
    
          C : np.ndarray
    
              3x3 x-axis DCM rotation
    
          x : double
    
              euler angle [rad]
    """
@typing.overload
def RotX(x: float) -> numpy.ndarray[numpy.float64[3, 3]]:
    """
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
    """
@typing.overload
def RotY(C: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.writeable, numpy.ndarray.flags.f_contiguous], y: float) -> None:
    """
          RotY
          ====
    
          Converts euler angle about y-axis into DCM rotation about y-axis
    
          Parameters
          ----------
    
          C : np.ndarray
    
              3x3 y-axis DCM rotation
    
          y : double
    
              euler angle [rad]
    """
@typing.overload
def RotY(y: float) -> numpy.ndarray[numpy.float64[3, 3]]:
    """
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
    """
@typing.overload
def RotZ(C: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.writeable, numpy.ndarray.flags.f_contiguous], z: float) -> None:
    """
          RotZ
          ====
    
          Converts euler angle about z-axis into DCM rotation about z-axis
    
          Parameters
          ----------
    
          C : np.ndarray
    
              3x3 z-axis DCM rotation
    
          z : double
    
              euler angle [rad]
    """
@typing.overload
def RotZ(z: float) -> numpy.ndarray[numpy.float64[3, 3]]:
    """
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
    """
@typing.overload
def dcm2euler(e: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], C: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.f_contiguous], IsNed: bool) -> None:
    """
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
    """
@typing.overload
def dcm2euler(C: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.f_contiguous], IsNed: bool) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def dcm2quat(q: numpy.ndarray[numpy.float64[4, 1], numpy.ndarray.flags.writeable], C: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.f_contiguous]) -> None:
    """
            dcm2quat
            ========
      
            Converts BODY-to-NAV DCM to corresponding BODY-to-NAV quaternion
              
            Parameters
            ----------
    
            q : np.ndarray
      
                size 4 NAV quaternion
      
            C : np.ndarray
      
                size 3x3 NAV DCM (ZYX)
    """
@typing.overload
def dcm2quat(C: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.f_contiguous]) -> numpy.ndarray[numpy.float64[4, 1]]:
    """
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
    """
@typing.overload
def euler2dcm(C: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.writeable, numpy.ndarray.flags.f_contiguous], e: numpy.ndarray[numpy.float64[3, 1]], IsNed: bool) -> None:
    """
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
    """
@typing.overload
def euler2dcm(e: numpy.ndarray[numpy.float64[3, 1]], IsNed: bool) -> numpy.ndarray[numpy.float64[3, 3]]:
    """
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
    """
@typing.overload
def euler2quat(q: numpy.ndarray[numpy.float64[4, 1], numpy.ndarray.flags.writeable], e: numpy.ndarray[numpy.float64[3, 1]], IsNed: bool) -> None:
    """
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
    """
@typing.overload
def euler2quat(q: numpy.ndarray[numpy.float64[3, 1]], IsNed: bool) -> numpy.ndarray[numpy.float64[4, 1]]:
    """
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
    """
@typing.overload
def quat2dcm(q: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.writeable, numpy.ndarray.flags.f_contiguous], C: numpy.ndarray[numpy.float64[4, 1]]) -> None:
    """
          quat2dcm
          ========
          
          Converts BODY-to-NAV quaternion to corresponding BODY-to-NAV DCM
            
          Parameters
          ----------
    
          q : np.ndarray
    
              size 4 NAV quaternion
    
          C : np.ndarray
    
              size 3x3 NAV DCM (ZYX)
    """
@typing.overload
def quat2dcm(q: numpy.ndarray[numpy.float64[4, 1]]) -> numpy.ndarray[numpy.float64[3, 3]]:
    """
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
    """
@typing.overload
def quat2euler(e: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], q: numpy.ndarray[numpy.float64[4, 1]], IsNed: bool) -> None:
    """
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
    """
@typing.overload
def quat2euler(q: numpy.ndarray[numpy.float64[4, 1]], IsNed: bool) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
