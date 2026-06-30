"""Laser source helpers."""

from __future__ import annotations

import numpy as np


def channel_wavelengths(center_nm: float, spacing_nm: float, channels: int) -> np.ndarray:
    """Return WDM channel wavelengths around a center wavelength."""

    if channels <= 0:
        raise ValueError("channels must be positive")
    return np.asarray([center_nm + spacing_nm * (chan - (channels + 1) / 2) for chan in range(1, channels + 1)])


def laser_source(
    peak_power_mw: float | list[float] | np.ndarray,
    samples: int,
    channels: int = 1,
    wavelength_nm: float | list[float] | np.ndarray = 1550.0,
    spacing_nm: float = 0.4,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Create a constant optical field for one or more laser channels.

    Returns `(field, wavelengths_nm, powers_mw)`.
    """

    if samples <= 0:
        raise ValueError("samples must be positive")
    if channels <= 0:
        raise ValueError("channels must be positive")

    powers = np.asarray(peak_power_mw, dtype=float)
    if powers.ndim == 0:
        powers = np.full(channels, float(powers))
    if len(powers) != channels:
        raise ValueError("peak_power_mw must be scalar or have one value per channel")

    wavelengths = np.asarray(wavelength_nm, dtype=float)
    if wavelengths.ndim == 0:
        wavelengths = channel_wavelengths(float(wavelengths), spacing_nm, channels)
    if len(wavelengths) != channels:
        raise ValueError("wavelength_nm must be scalar or have one value per channel")

    field = np.ones((samples, channels), dtype=complex) * np.sqrt(powers)
    if channels == 1:
        field = field[:, 0]

    return field, wavelengths, powers
