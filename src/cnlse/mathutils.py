"""Small math utilities ported from the original MATLAB simulator.

Original MATLAB files credited to Ryan Maher Hamze, Krayem Ahmad Ajjawi,
and Meriam Hawwa.
"""

from __future__ import annotations

import numpy as np


def nmod(a, n: int):
    """Reduce integer values into the MATLAB-style inclusive range 1..n.

    This ports `Nmod.m`, where `nmod(0, 8) == 8` and `nmod(9, 8) == 1`.
    Scalars return Python integers; arrays return NumPy arrays.
    """

    if n <= 0:
        raise ValueError("n must be positive")

    result = np.mod(np.asarray(a) - 1 - n, n) + 1
    if np.isscalar(a):
        return int(result)
    return result


def fastshift(x, shift: int):
    """Circularly shift a vector or matrix along the first dimension.

    This ports `Fastshift.m`. For vectors, it behaves like `circshift(x,n)`.
    For matrices, rows are shifted and columns are preserved.
    """

    if isinstance(shift, (list, tuple, np.ndarray)):
        arr = np.asarray(shift)
        if arr.size != 1:
            raise ValueError("shift must be a single value")
        shift = int(arr.item())

    values = np.asarray(x)
    return np.roll(values, int(shift), axis=0)


def fastexp(x):
    """Return exp(1j*x) using cosine and sine, ported from `Fastexp.m`."""

    values = np.asarray(x)
    result = np.cos(values) + 1j * np.sin(values)
    if np.isscalar(x):
        return complex(result)
    return result
