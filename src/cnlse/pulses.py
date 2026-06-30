"""Pulse-shaping helpers for optical modulation examples."""

from __future__ import annotations

import numpy as np


def ook_rectangular(bits: list[int], samples_per_symbol: int) -> np.ndarray:
    """Create a rectangular on-off-keying power envelope."""

    if samples_per_symbol <= 0:
        raise ValueError("samples_per_symbol must be positive")
    return np.repeat(np.asarray(bits, dtype=float), samples_per_symbol)


def gaussian_smooth(values: np.ndarray, sigma_samples: float) -> np.ndarray:
    """Smooth a 1D waveform with a lightweight Gaussian kernel."""

    if sigma_samples <= 0:
        return values.copy()

    radius = max(1, int(4 * sigma_samples))
    x = np.arange(-radius, radius + 1)
    kernel = np.exp(-0.5 * (x / sigma_samples) ** 2)
    kernel /= kernel.sum()
    return np.convolve(values, kernel, mode="same")
