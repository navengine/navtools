"""

Binary-Ops
==========

Useful binary operations.
"""

from __future__ import annotations

__all__ = ["GetBit", "GetBits", "UnsetBit", "MultiXor", "SetBit", "SetBitTo", "TwosComp"]

def GetBit(x: int, n: int) -> bool:
    """
    GetBit
    ========

    Check the value of a bit

    Parameters
    ----------

    x : uint32

        Number to modify

    n : uint8

        Position of bit to set (Position 0 is MSB and 31 is LSB by default)
    """

def GetBits(x: int, b: int, n: int) -> int:
    """
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
    """

def UnsetBit(x: int, n: int) -> None:
    """
    UnsetBit
    ========

    Set a data bit to 0

    Parameters
    ----------

    x : uint32

        Number to modify

    n : uint8

        Position of bit to set (Position 0 is MSB and 31 is LSB by default)
    """

def MultiXor(x: int, n: int, Size: int) -> bool:
    """
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
    """

def SetBit(x: int, n: int) -> None:
    """
    SetBit
    ======

    Set a data bit to 1

    Parameters
    ----------

    x : uint32

        Number to modify

    n : uint8

        Position of bit to set (Position 0 is MSB and 31 is LSB by default)
    """

def SetBitTo(x: int, n: int, b: bool) -> None:
    """
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
    """

def TwosComp(x: int, n: int) -> float:
    """
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
    """
