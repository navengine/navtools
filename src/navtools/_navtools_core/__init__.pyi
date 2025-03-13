"""

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
      
"""
from __future__ import annotations
import numpy
from . import attitude
from . import binaryops
from . import frames
from . import math
from . import models
__all__ = ['BOLTZMANN', 'DEG2RAD', 'F', 'GAUSS_TO_TESLA', 'GRAVITY', 'HALF_PI', 'J2', 'J3', 'J4', 'LIGHT_SPEED', 'METERS_TO_FOOT', 'MINUTES_PER_DAY', 'PI', 'PI_SQU', 'RAD2DEG', 'RE', 'SQRT_PI', 'TWO_PI', 'WGS84_E', 'WGS84_E2', 'WGS84_F', 'WGS84_MU', 'WGS84_OMEGA', 'WGS84_OMEGA_SKEW', 'WGS84_OMEGA_VEC', 'WGS84_R0', 'WGS84_RP', 'attitude', 'binaryops', 'deg2rad', 'frames', 'math', 'models', 'rad2deg']
def deg2rad(arg0: float) -> float:
    ...
def rad2deg(arg0: float) -> float:
    ...
BOLTZMANN: float = 1.38e-23
DEG2RAD: float = 0.017453292519943295
F: float = -4.442807633e-10
GAUSS_TO_TESLA: float = 0.0001
GRAVITY: float = 9.80665
HALF_PI: float = 1.5707963267948966
J2: float = 0.0010826269
J3: float = -2.5323e-06
J4: float = -1.6204e-06
LIGHT_SPEED: float = 299792458.0
METERS_TO_FOOT: float = 0.3048
MINUTES_PER_DAY: float = 1440.0
PI: float = 3.141592653589793
PI_SQU: float = 9.869604401089358
RAD2DEG: float = 57.29577951308232
RE: float = 6378136.3
SQRT_PI: float = 1.772453850905516
TWO_PI: float = 6.283185307179586
WGS84_E: float = 0.0818191908429654
WGS84_E2: float = 0.00669437999019758
WGS84_F: float = 0.00335281066477569
WGS84_MU: float = 398600500000000.0
WGS84_OMEGA: float = 7.2921151467e-05
WGS84_OMEGA_SKEW: numpy.ndarray  # value = array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],...
WGS84_OMEGA_VEC: numpy.ndarray  # value = array([0.00000000e+00, 0.00000000e+00, 7.29211515e-05])
WGS84_R0: float = 6378137.0
WGS84_RP: float = 6356752.314245
__version__: str = '1.0.0'
