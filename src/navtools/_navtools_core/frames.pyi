"""

        Frames
        ======
        
        Common coordinate frame transformations.
"""
from __future__ import annotations
import numpy
import typing
__all__ = ['ecef2aer', 'ecef2eci', 'ecef2eciDcm', 'ecef2ecia', 'ecef2eciv', 'ecef2eciw', 'ecef2enu', 'ecef2enuDcm', 'ecef2enua', 'ecef2enuv', 'ecef2enuw', 'ecef2lla', 'ecef2ned', 'ecef2nedDcm', 'ecef2neda', 'ecef2nedv', 'ecef2nedw', 'eci2aer', 'eci2ecef', 'eci2ecefDcm', 'eci2ecefa', 'eci2ecefv', 'eci2ecefw', 'eci2enu', 'eci2enuDcm', 'eci2enua', 'eci2enuv', 'eci2enuw', 'eci2lla', 'eci2ned', 'eci2nedDcm', 'eci2neda', 'eci2nedv', 'eci2nedw', 'enu2aer', 'enu2ecef', 'enu2ecefDcm', 'enu2ecefa', 'enu2ecefv', 'enu2ecefw', 'enu2eci', 'enu2eciDcm', 'enu2ecia', 'enu2eciv', 'enu2lla', 'enu2nedDcm', 'lla2aer', 'lla2ecef', 'lla2eci', 'lla2enu', 'lla2ned', 'ned2aer', 'ned2ecef', 'ned2ecefDcm', 'ned2ecefa', 'ned2ecefv', 'ned2ecefw', 'ned2eci', 'ned2eciDcm', 'ned2ecia', 'ned2eciv', 'ned2eciw', 'ned2enuDcm', 'ned2lla']
@typing.overload
def ecef2aer(aer: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], xyzR: numpy.ndarray[numpy.float64[3, 1]], xyzT: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    """
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
    """
@typing.overload
def ecef2aer(xyzR: numpy.ndarray[numpy.float64[3, 1]], xyzT: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def ecef2eci(xyz: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], eci: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> None:
    """
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
    """
@typing.overload
def ecef2eci(eci: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def ecef2eciDcm(C_e_i: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.writeable, numpy.ndarray.flags.f_contiguous], dt: float) -> None:
    """
          ecef2eciDcm
          ===========
    
          Earth-Centered-Earth-Fixed to Earth-Centered-Inertial direction cosine matrix
    
          Parameters
          ----------
    
          C_e_i : np.ndarray
    
              3x3 ECEF->ECI direction cosine matrix
    
          dt : double 
    
              time elapsed between frames [s]
    """
@typing.overload
def ecef2eciDcm(dt: float) -> numpy.ndarray[numpy.float64[3, 3]]:
    """
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
    """
@typing.overload
def ecef2ecia(arg0: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]], arg3: numpy.ndarray[numpy.float64[3, 1]], arg4: float) -> None:
    ...
@typing.overload
def ecef2ecia(arg0: numpy.ndarray[numpy.float64[3, 1]], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]], arg3: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def ecef2eciv(v_ib_i: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], r_eb_e: numpy.ndarray[numpy.float64[3, 1]], v_eb_e: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> None:
    """
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
    """
@typing.overload
def ecef2eciv(r_eb_e: numpy.ndarray[numpy.float64[3, 1]], v_eb_e: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def ecef2eciw(arg0: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: float) -> None:
    ...
@typing.overload
def ecef2eciw(arg0: numpy.ndarray[numpy.float64[3, 1]], arg1: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def ecef2enu(enu: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], xyz: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    """
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
    """
@typing.overload
def ecef2enu(xyz: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def ecef2enuDcm(C_e_n: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.writeable, numpy.ndarray.flags.f_contiguous], lla: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    """
          ecef2enuDcm
          ===========
    
          Earth-Centered-Earth-Fixed to East-North-Up direction cosine matrix
    
          Parameters
          ----------
    
          C_e_n : np.ndarray
    
              3x3 ECEF->ENU direction cosine matrix
    
          lla : np.ndarray
    
              3x1 Latitude, Longitude, Height [rad, rad, m]
    """
@typing.overload
def ecef2enuDcm(lla: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 3]]:
    """
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
    """
@typing.overload
def ecef2enua(arg0: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    ...
@typing.overload
def ecef2enua(arg0: numpy.ndarray[numpy.float64[3, 1]], arg1: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def ecef2enuv(v_nb_e: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], r_eb_e: numpy.ndarray[numpy.float64[3, 1]], v_eb_e: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    """
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
    """
@typing.overload
def ecef2enuv(r_eb_e: numpy.ndarray[numpy.float64[3, 1]], v_eb_e: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def ecef2enuw(arg0: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    ...
@typing.overload
def ecef2enuw(arg0: numpy.ndarray[numpy.float64[3, 1]], arg1: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def ecef2lla(lla: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], xyz: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    """
          ecef2lla
          ========
          
          Earth-Centered-Earth-Fixed to Latitude-Longitude-Height position coordinates
    
          Parameters
          ----------
    
          lla : np.ndarray
    
              3x1 Latitude, Longitude, Height [rad, rad, m]
    
          xyz : np.ndarray
    
              3x1 ECEF Position [m]
    """
@typing.overload
def ecef2lla(xyz: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def ecef2ned(ned: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], xyz: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    """
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
    """
@typing.overload
def ecef2ned(xyz: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def ecef2nedDcm(C_e_n: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.writeable, numpy.ndarray.flags.f_contiguous], lla: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    """
          ecef2nedDcm
          ===========
    
          Earth-Centered-Earth-Fixed to North-East-Down direction cosine matrix
    
          Parameters
          ----------
    
          C_e_n : np.ndarray
    
              3x3 ECEF->NED direction cosine matrix
    
          lla : np.ndarray
    
              3x1 Latitude, Longitude, Height [rad, rad, m]
    """
@typing.overload
def ecef2nedDcm(lla: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 3]]:
    """
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
    """
@typing.overload
def ecef2neda(arg0: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    ...
@typing.overload
def ecef2neda(arg0: numpy.ndarray[numpy.float64[3, 1]], arg1: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def ecef2nedv(v_nb_e: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], r_eb_e: numpy.ndarray[numpy.float64[3, 1]], v_eb_e: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    """
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
    """
@typing.overload
def ecef2nedv(r_eb_e: numpy.ndarray[numpy.float64[3, 1]], v_eb_e: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def ecef2nedw(arg0: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    ...
@typing.overload
def ecef2nedw(arg0: numpy.ndarray[numpy.float64[3, 1]], arg1: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def eci2aer(aer: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], eciR: numpy.ndarray[numpy.float64[3, 1]], eciT: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> None:
    """
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
    """
@typing.overload
def eci2aer(eciR: numpy.ndarray[numpy.float64[3, 1]], eciT: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def eci2ecef(eci: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], xyz: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> None:
    """
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
    """
@typing.overload
def eci2ecef(xyz: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def eci2ecefDcm(C_i_e: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.writeable, numpy.ndarray.flags.f_contiguous], dt: float) -> None:
    """
          eci2ecefDcm
          ===========
    
          Earth-Centered-Inertial to Earth-Centered-Earth-Fixed direction cosine matrix
    
          Parameters
          ----------
    
          C_i_e : np.ndarray
    
              3x3 ECI->ECEF direction cosine matrix
    
          dt : double 
              time elapsed between frames [s]
    """
@typing.overload
def eci2ecefDcm(dt: float) -> numpy.ndarray[numpy.float64[3, 3]]:
    """
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
    """
@typing.overload
def eci2ecefa(arg0: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]], arg3: numpy.ndarray[numpy.float64[3, 1]], arg4: float) -> None:
    ...
@typing.overload
def eci2ecefa(arg0: numpy.ndarray[numpy.float64[3, 1]], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]], arg3: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def eci2ecefv(v_eb_e: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], r_ib_i: numpy.ndarray[numpy.float64[3, 1]], v_ib_i: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> None:
    """
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
    """
@typing.overload
def eci2ecefv(r_ib_i: numpy.ndarray[numpy.float64[3, 1]], v_ib_i: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def eci2ecefw(arg0: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: float) -> None:
    ...
@typing.overload
def eci2ecefw(arg0: numpy.ndarray[numpy.float64[3, 1]], arg1: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def eci2enu(enu: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], eci: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> None:
    """
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
    """
@typing.overload
def eci2enu(eci: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def eci2enuDcm(C_i_n: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.writeable, numpy.ndarray.flags.f_contiguous], lla: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> None:
    """
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
    """
@typing.overload
def eci2enuDcm(lla: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> numpy.ndarray[numpy.float64[3, 3]]:
    """
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
    """
@typing.overload
def eci2enua(arg0: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]], arg3: numpy.ndarray[numpy.float64[3, 1]], arg4: numpy.ndarray[numpy.float64[3, 1]], arg5: float) -> None:
    ...
@typing.overload
def eci2enua(arg0: numpy.ndarray[numpy.float64[3, 1]], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]], arg3: numpy.ndarray[numpy.float64[3, 1]], arg4: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def eci2enuv(v_bn_e: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], r_ib_i: numpy.ndarray[numpy.float64[3, 1]], v_ib_i: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> None:
    """
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
    """
@typing.overload
def eci2enuv(r_ib_i: numpy.ndarray[numpy.float64[3, 1]], v_ib_i: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def eci2enuw(arg0: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]], arg3: numpy.ndarray[numpy.float64[3, 1]], arg4: numpy.ndarray[numpy.float64[3, 1]], arg5: float) -> None:
    ...
@typing.overload
def eci2enuw(arg0: numpy.ndarray[numpy.float64[3, 1]], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]], arg3: numpy.ndarray[numpy.float64[3, 1]], arg4: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def eci2lla(lla: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], eci: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> None:
    """
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
    """
@typing.overload
def eci2lla(eci: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def eci2ned(ned: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], eci: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> None:
    """
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
    """
@typing.overload
def eci2ned(eci: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def eci2nedDcm(C_i_n: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.writeable, numpy.ndarray.flags.f_contiguous], lla: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> None:
    """
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
    """
@typing.overload
def eci2nedDcm(lla: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> numpy.ndarray[numpy.float64[3, 3]]:
    """
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
    """
@typing.overload
def eci2neda(arg0: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]], arg3: numpy.ndarray[numpy.float64[3, 1]], arg4: numpy.ndarray[numpy.float64[3, 1]], arg5: float) -> None:
    ...
@typing.overload
def eci2neda(arg0: numpy.ndarray[numpy.float64[3, 1]], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]], arg3: numpy.ndarray[numpy.float64[3, 1]], arg4: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def eci2nedv(v_bn_e: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], r_ib_i: numpy.ndarray[numpy.float64[3, 1]], v_ib_i: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> None:
    """
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
    """
@typing.overload
def eci2nedv(r_ib_i: numpy.ndarray[numpy.float64[3, 1]], v_ib_i: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def eci2nedw(arg0: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]], arg3: numpy.ndarray[numpy.float64[3, 1]], arg4: numpy.ndarray[numpy.float64[3, 1]], arg5: float) -> None:
    ...
@typing.overload
def eci2nedw(arg0: numpy.ndarray[numpy.float64[3, 1]], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]], arg3: numpy.ndarray[numpy.float64[3, 1]], arg4: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def enu2aer(aer: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], enuR: numpy.ndarray[numpy.float64[3, 1]], enuT: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    """
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
    """
@typing.overload
def enu2aer(enuR: numpy.ndarray[numpy.float64[3, 1]], enuT: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def enu2ecef(xyz: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], enu: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    """
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
    """
@typing.overload
def enu2ecef(enu: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def enu2ecefDcm(C_e_n: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.writeable, numpy.ndarray.flags.f_contiguous], lla: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    """
          enu2ecefDcm
          ===========
    
          East-North-Up to Earth-Centered-Earth-Fixed direction cosine matrix
    
          Parameters
          ----------
    
          C_e_n : np.ndarray
    
              3x3 ENU->ECEF direction cosine matrix
    
          lla : np.ndarray
    
              3x1 Latitude, Longitude, Height [rad, rad, m]
    """
@typing.overload
def enu2ecefDcm(lla: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 3]]:
    """
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
    """
@typing.overload
def enu2ecefa(arg0: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    ...
@typing.overload
def enu2ecefa(arg0: numpy.ndarray[numpy.float64[3, 1]], arg1: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def enu2ecefv(v_eb_e: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], r_nb_e: numpy.ndarray[numpy.float64[3, 1]], v_nb_e: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    """
          enu2ecefv
          =========
    
          Converts East-North-Up to Earth-Centered-Earth-Fixed velocity
    
          Parameters
          ----------
    
          v_eb_e : np.ndarray
    
              3x1 ECEF velocity [m/s]
    
          r_nb_e : np.ndarray
    
              3x1 ENU position [m]
    
          v_nb_e : np.ndarray
    
              3x1 ENU velocity [m/s]
    """
@typing.overload
def enu2ecefv(r_nb_e: numpy.ndarray[numpy.float64[3, 1]], v_nb_e: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
          enu2ecefv
          =========
    
          Converts East-North-Up to Earth-Centered-Earth-Fixed velocity
    
          Parameters
          ----------
    
          r_nb_e : np.ndarray
    
              3x1 ENU position [m]
    
          v_nb_e : np.ndarray
    
              3x1 ENU velocity [m/s]
    
          Returns
          -------
    
          v_eb_e : np.ndarray
    
              3x1 ECEF velocity [m/s]
    """
@typing.overload
def enu2ecefw(arg0: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    ...
@typing.overload
def enu2ecefw(arg0: numpy.ndarray[numpy.float64[3, 1]], arg1: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def enu2eci(eci: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], enu: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> None:
    """
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
    """
@typing.overload
def enu2eci(enu: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def enu2eciDcm(C_n_i: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.writeable, numpy.ndarray.flags.f_contiguous], lla: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> None:
    """
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
    """
@typing.overload
def enu2eciDcm(lla: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> numpy.ndarray[numpy.float64[3, 3]]:
    """
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
    """
@typing.overload
def enu2ecia(arg0: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]], arg3: numpy.ndarray[numpy.float64[3, 1]], arg4: numpy.ndarray[numpy.float64[3, 1]], arg5: float) -> None:
    ...
@typing.overload
def enu2ecia(arg0: numpy.ndarray[numpy.float64[3, 1]], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]], arg3: numpy.ndarray[numpy.float64[3, 1]], arg4: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def enu2eciv(v_in_i: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], r_nb_e: numpy.ndarray[numpy.float64[3, 1]], v_nb_e: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> None:
    """
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
    """
@typing.overload
def enu2eciv(r_nb_e: numpy.ndarray[numpy.float64[3, 1]], v_nb_e: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def enu2lla(lla: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], enu: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    """
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
    """
@typing.overload
def enu2lla(enu: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def enu2nedDcm(C_n_n: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.writeable, numpy.ndarray.flags.f_contiguous]) -> None:
    """
          enu2nedDcm
          ==========
    
          East-North-Up to North-East-Down direction cosine matrix
    
          Parameters
          ----------
    
          C_n_n : np.ndarray
    
              3x3 ENU->NED direction cosine matrix
    """
@typing.overload
def enu2nedDcm() -> numpy.ndarray[numpy.float64[3, 3]]:
    """
          enu2nedDcm
          ==========
    
          East-North-Up to North-East-Down direction cosine matrix
    
          Returns
          -------
    
          C_n_n : np.ndarray
    
              3x3 ENU->NED direction cosine matrix
    """
@typing.overload
def lla2aer(aer: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], llaR: numpy.ndarray[numpy.float64[3, 1]], llaT: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    """
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
    """
@typing.overload
def lla2aer(llaR: numpy.ndarray[numpy.float64[3, 1]], llaT: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def lla2ecef(xyz: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], lla: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    """
          lla2ecef
          ========
    
          Latitude-Longitude-Height to Earth-Centered-Earth-Fixed position coordinates
    
          Parameters
          ----------
    
          xyz : np.ndarray
    
              3x1 ECEF Position [m]
    
          lla : np.ndarray
    
              3x1 Latitude, Longitude, Height [rad, rad, m]
    """
@typing.overload
def lla2ecef(lla: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def lla2eci(eci: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], lla: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> None:
    """
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
    """
@typing.overload
def lla2eci(lla: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def lla2enu(enu: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], lla: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    """
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
    """
@typing.overload
def lla2enu(lla: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def lla2ned(ned: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], lla: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    """
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
    """
@typing.overload
def lla2ned(lla: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def ned2aer(aer: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], nedR: numpy.ndarray[numpy.float64[3, 1]], nedT: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    """
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
    """
@typing.overload
def ned2aer(nedR: numpy.ndarray[numpy.float64[3, 1]], nedT: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def ned2ecef(xyz: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], ned: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    """
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
    """
@typing.overload
def ned2ecef(ned: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def ned2ecefDcm(C_e_n: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.writeable, numpy.ndarray.flags.f_contiguous], lla: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    """
          ned2ecefDcm
          ===========
    
          North-East-Down to Earth-Centered-Earth-Fixed direction cosine matrix
    
          Parameters
          ----------
    
          C_e_n : np.ndarray
    
              3x3 NED->ECEF direction cosine matrix
    
          lla : np.ndarray
    
              3x1 Latitude, Longitude, Height [rad, rad, m]
    """
@typing.overload
def ned2ecefDcm(lla: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 3]]:
    """
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
    """
@typing.overload
def ned2ecefa(arg0: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    ...
@typing.overload
def ned2ecefa(arg0: numpy.ndarray[numpy.float64[3, 1]], arg1: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def ned2ecefv(v_eb_e: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], r_nb_e: numpy.ndarray[numpy.float64[3, 1]], v_nb_e: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    """
          ned2ecefv
          =========
    
          Converts North-East-Down to Earth-Centered-Earth-Fixed velocity
    
          Parameters
          ----------
    
          v_eb_e : np.ndarray
    
              3x1 ECEF velocity [m/s]
    
          r_nb_e : np.ndarray
    
              3x1 NED position [m]
    
          v_nb_e : np.ndarray
    
              3x1 NED velocity [m/s]
    """
@typing.overload
def ned2ecefv(r_nb_e: numpy.ndarray[numpy.float64[3, 1]], v_nb_e: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
          ned2ecefv
          =========
    
          Converts North-East-Down to Earth-Centered-Earth-Fixed velocity
    
          Parameters
          ----------
    
          r_nb_e : np.ndarray
    
              3x1 NED position [m]
    
          v_nb_e : np.ndarray
    
              3x1 NED velocity [m/s]
    
          Returns
          -------
    
          v_eb_e : np.ndarray
    
              3x1 ECEF velocity [m/s]
    """
@typing.overload
def ned2ecefw(arg0: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    ...
@typing.overload
def ned2ecefw(arg0: numpy.ndarray[numpy.float64[3, 1]], arg1: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def ned2eci(eci: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], ned: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> None:
    """
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
    """
@typing.overload
def ned2eci(ned: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def ned2eciDcm(C_n_i: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.writeable, numpy.ndarray.flags.f_contiguous], lla: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> None:
    """
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
    """
@typing.overload
def ned2eciDcm(lla: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> numpy.ndarray[numpy.float64[3, 3]]:
    """
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
    """
@typing.overload
def ned2ecia(arg0: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]], arg3: numpy.ndarray[numpy.float64[3, 1]], arg4: numpy.ndarray[numpy.float64[3, 1]], arg5: float) -> None:
    ...
@typing.overload
def ned2ecia(arg0: numpy.ndarray[numpy.float64[3, 1]], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]], arg3: numpy.ndarray[numpy.float64[3, 1]], arg4: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def ned2eciv(v_in_i: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], r_nb_e: numpy.ndarray[numpy.float64[3, 1]], v_nb_e: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> None:
    """
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
    """
@typing.overload
def ned2eciv(r_nb_e: numpy.ndarray[numpy.float64[3, 1]], v_nb_e: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]], dt: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def ned2eciw(arg0: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]], arg3: numpy.ndarray[numpy.float64[3, 1]], arg4: float) -> None:
    ...
@typing.overload
def ned2eciw(arg0: numpy.ndarray[numpy.float64[3, 1]], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]], arg3: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def ned2eciw(arg0: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]], arg3: numpy.ndarray[numpy.float64[3, 1]], arg4: float) -> None:
    ...
@typing.overload
def ned2eciw(arg0: numpy.ndarray[numpy.float64[3, 1]], arg1: numpy.ndarray[numpy.float64[3, 1]], arg2: numpy.ndarray[numpy.float64[3, 1]], arg3: float) -> numpy.ndarray[numpy.float64[3, 1]]:
    ...
@typing.overload
def ned2enuDcm(C_n_n: numpy.ndarray[numpy.float64[3, 3], numpy.ndarray.flags.writeable, numpy.ndarray.flags.f_contiguous]) -> None:
    """
          ned2enuDcm
          ==========
    
          North-East-Down to East-North-Up direction cosine matrix
    
          Parameters
          ----------
    
          C_n_n : np.ndarray
    
              3x3 NED->ENU direction cosine matrix
    """
@typing.overload
def ned2enuDcm() -> numpy.ndarray[numpy.float64[3, 3]]:
    """
          ned2enuDcm
          ==========
    
          North-East-Down to East-North-Up direction cosine matrix
    
          Returns
          -------
    
          C_n_n : np.ndarray
    
              3x3 NED->ENU direction cosine matrix
    """
@typing.overload
def ned2lla(lla: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], ned: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    """
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
    """
@typing.overload
def ned2lla(ned: numpy.ndarray[numpy.float64[3, 1]], lla0: numpy.ndarray[numpy.float64[3, 1]]) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
