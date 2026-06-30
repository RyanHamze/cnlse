"""Electrical driving signal generation.

This module ports the first practical path of the original MATLAB
`Electricsource.m`: OOK modulation with ideal or cosine-roll pulses.

Original MATLAB implementation credited to Ryan Maher Hamze, Krayem Ahmad
Ajjawi, and Meriam Hawwa.
"""

from __future__ import annotations

import numpy as np


def pulse_shape(kind: str, samples_per_symbol: int, duty: float = 1.0, roll: float = 0.2) -> np.ndarray:
    """Create a two-symbol elementary pulse.

    Parameters mirror the MATLAB `pulseshape` helper. The result has length
    `2 * samples_per_symbol` because roll-off can spread into the neighboring
    symbol interval.
    """

    if samples_per_symbol <= 0:
        raise ValueError("samples_per_symbol must be positive")
    if not 0 <= duty <= 1:
        raise ValueError("duty must be between 0 and 1")
    if not 0 <= roll <= 1:
        raise ValueError("roll must be between 0 and 1")

    kind = kind.lower()
    nt = samples_per_symbol
    pulse = np.zeros(2 * nt, dtype=float)

    if kind == "dirac":
        pulse[nt] = 1.0
        return pulse

    if kind == "idpulse":
        half = round(0.5 * duty * nt)
        if half == 0:
            pulse[nt] = 1.0
        else:
            pulse[nt - half : nt + half] = 1.0
        return pulse

    if kind != "cosroll":
        raise NotImplementedError(f"pulse kind is not implemented yet: {kind}")

    flat = round(0.5 * (1 - roll) * duty * nt)
    right = int(duty * nt - flat - 1)
    if flat > 0:
        pulse[nt : nt + flat] = 1.0

    transition = np.arange(flat, max(flat, right))
    half_period = duty * nt - 2 * flat
    if half_period != 0 and transition.size:
        pulse[nt + transition] = 0.5 * (1 + np.cos(np.pi / half_period * (transition - flat + 0.5)))

    pulse[:nt] = pulse[nt: 2 * nt][::-1]
    return pulse


def ook_signal(
    bits: list[int] | np.ndarray,
    samples_per_symbol: int,
    pulse: str = "cosroll",
    duty: float = 1.0,
    roll: float = 0.2,
    shape_power: bool = False,
) -> np.ndarray:
    """Create an OOK electrical drive signal.

    This is the Python equivalent of the common MATLAB call:

    `electricsource(pat, 'ook', symbrate, 'cosroll', duty, roll)`
    """

    bits_array = np.asarray(bits, dtype=float)
    if bits_array.ndim != 1:
        raise ValueError("bits must be a one-dimensional sequence")
    if np.any((bits_array != 0) & (bits_array != 1)):
        raise ValueError("OOK bits must contain only 0 and 1")

    nt = samples_per_symbol
    n_symbols = len(bits_array)
    nfft = n_symbols * nt
    elementary = pulse_shape(pulse, nt, duty=duty, roll=roll)
    signal = np.zeros(nfft, dtype=float)

    # Match MATLAB's cyclic construction: the first pulse wraps around the
    # periodic sequence boundary, then each following pulse is added normally.
    signal[nfft - nt : nfft] = bits_array[0] * elementary[:nt]
    signal[:nt] = bits_array[0] * elementary[nt : 2 * nt]

    for symbol_index in range(1, n_symbols):
        start = (symbol_index - 1) * nt
        end = (symbol_index + 1) * nt
        signal[start:end] += bits_array[symbol_index] * elementary

    if shape_power:
        signal = np.sqrt(np.abs(signal))

    return signal
