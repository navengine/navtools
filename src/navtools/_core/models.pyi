"""
Simple Earth models commonly used in navigation equations.
"""
from __future__ import annotations
import numpy
import typing
__all__ = ['CoriolisRate', 'EarthRate', 'EcefGravity', 'GeocentricRadius', 'LocalGravity', 'MeridianRadius', 'RadiiOfCurvature', 'Somigliana', 'TransAndMerRadii', 'TransportRate', 'TransverseRadius']
@typing.overload
def CoriolisRate(coriolis: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], lla: numpy.ndarray[numpy.float64[3, 1]], v_nb_e: numpy.ndarray[numpy.float64[3, 1]], IsNed: bool) -> None:
    """
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
    """
@typing.overload
def CoriolisRate(lla: numpy.ndarray[numpy.float64[3, 1]], v_nb_e: numpy.ndarray[numpy.float64[3, 1]], IsNed: bool) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def EarthRate(w_ie_n: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], phi: float, IsNed: bool) -> None:
    """
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
    """
@typing.overload
def EarthRate(phi: float, IsNed: bool) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def EarthRate(phi: numpy.ndarray[numpy.float64[3, 1]], v_nb_e: numpy.ndarray[numpy.float64[3, 1]], IsNed: bool) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
def EcefGravity(g: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], gamma: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], xyz: numpy.ndarray[numpy.float64[3, 1]]) -> None:
    """
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
    """
@typing.overload
def GeocentricRadius(R_es_e: float, phi: float) -> None:
    """
          GeocentricRadius
          ================
    
          Calculates the geocentric radius relative to user latitude
    
          Parameters
          ----------
    
          R_es_e : double
    
              Geocentric radius [m]
    
          phi : double
    
              Latitude [rad]
    """
@typing.overload
def GeocentricRadius(phi: float) -> float:
    """
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
    """
@typing.overload
def LocalGravity(g: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], lla: numpy.ndarray[numpy.float64[3, 1]], IsNed: bool) -> None:
    """
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
    """
@typing.overload
def LocalGravity(lla: numpy.ndarray[numpy.float64[3, 1]], IsNed: bool) -> numpy.ndarray[numpy.float64[3, 1]]:
    """
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
    """
@typing.overload
def MeridianRadius(Rn: float, phi: float) -> None:
    """
          MeridianRadius
          ==============
    
          Calculates the meridian radius relative to user latitude
    
          Parameters
          ----------
    
          Rn : double
    
              Meridian radius [m]
    
          phi : double
    
              Latitude [rad]
    """
@typing.overload
def MeridianRadius(phi: float) -> float:
    """
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
    """
def RadiiOfCurvature(Re: float, Rn: float, R_es_e: float, phi: float) -> None:
    """
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
    """
@typing.overload
def Somigliana(g0: float, phi: float) -> None:
    """
          Somigliana
          ==========
    
          Calculates the somilgiana model reference gravity
    
          Parameters
          ----------
    
          g0 : double
    
              Somigliana gravitation
    
          phi : double
    
              Latitude [rad]
    """
@typing.overload
def Somigliana(phi: float) -> float:
    """
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
    """
def TransAndMerRadii(Re: float, Rn: float, phi: float) -> None:
    """
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
    """
def TransportRate(w_en_n: numpy.ndarray[numpy.float64[3, 1], numpy.ndarray.flags.writeable], phi: numpy.ndarray[numpy.float64[3, 1]], v_nb_e: numpy.ndarray[numpy.float64[3, 1]], IsNed: bool) -> None:
    """
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
    """
@typing.overload
def TransverseRadius(Re: float, phi: float) -> None:
    """
          TransverseRadius
          ================
    
          Calculates the transverse radius relative to user latitude
    
          Parameters
          ----------
    
          Re : double
    
              Transverse radius [m]
    
          phi : double
    
              Latitude [rad]
    """
@typing.overload
def TransverseRadius(phi: float) -> float:
    """
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
    """
