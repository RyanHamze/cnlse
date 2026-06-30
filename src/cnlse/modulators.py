"""Optical modulators."""

from __future__ import annotations

import math

import numpy as np

from .mathutils import fastexp


def mach_zehnder(
    optical_field,
    drive_signal,
    extinction_ratio_db: float = math.inf,
    bias: float = 0.0,
    amplitude: float = 1.0,
    nochirp: bool = False,
):
    """Modulate an optical field with a Mach-Zehnder interferometer.

    Ports the transfer equation from original `Mz_Modulator.m`.
    """

    ein = np.asarray(optical_field, dtype=complex)
    modsig = np.asarray(drive_signal, dtype=float)

    if ein.shape[0] != modsig.shape[0]:
        raise ValueError("optical_field and drive_signal must have the same sample length")

    if math.isinf(extinction_ratio_db):
        exr_linear = 0.0
    else:
        exr_linear = 10 ** (-extinction_ratio_db / 10)

    gamma = (1 - math.sqrt(exr_linear)) / (math.sqrt(exr_linear) + 1)

    if nochirp:
        phi_upper = modsig * math.pi / (1 + gamma**2) * amplitude + bias
        phi_lower = modsig * -gamma**2 * math.pi / (1 + gamma**2) * amplitude + bias
    else:
        phi_upper = modsig * math.pi / 2 * amplitude + bias
        phi_lower = modsig * -math.pi / 2 * amplitude + bias

    if ein.ndim == 2 and modsig.ndim == 1:
        phi_upper = phi_upper[:, np.newaxis]
        phi_lower = phi_lower[:, np.newaxis]

    return 1j * ein * (fastexp(phi_lower) - gamma * fastexp(phi_upper)) / (1 + gamma)
