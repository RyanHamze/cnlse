"""Simulation state.

The original MATLAB project used global `GSTATE` and `CONSTANTS`
structures. Python keeps the same ideas in explicit dataclasses so examples
and tests can pass state around without hidden globals.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class SimulationState:
    """Container for common optical simulation state."""

    symbols: int
    samples_per_symbol: int
    channels: int
    symbol_rate_gbaud: float | None = None
    normalized_frequency: np.ndarray = field(init=False)
    field_x: np.ndarray | None = None
    field_y: np.ndarray | None = None
    field_x_tx: np.ndarray | None = None
    field_y_tx: np.ndarray | None = None
    delay: np.ndarray | None = None
    dispersion: np.ndarray | None = None
    wavelengths_nm: np.ndarray | None = None
    powers_mw: np.ndarray | None = None

    def __post_init__(self) -> None:
        if self.symbols <= 0:
            raise ValueError("symbols must be positive")
        if self.samples_per_symbol <= 0:
            raise ValueError("samples_per_symbol must be positive")
        if self.channels <= 0:
            raise ValueError("channels must be positive")

        step = 1 / self.symbols
        values = np.arange(-self.samples_per_symbol / 2, self.samples_per_symbol / 2, step)
        self.normalized_frequency = np.fft.fftshift(values)

    @property
    def samples(self) -> int:
        """Total FFT/sample count."""

        return self.symbols * self.samples_per_symbol


def create_state(symbols: int, samples_per_symbol: int, channels: int) -> SimulationState:
    """Create a new simulation state, equivalent to the core of `Reset_All.m`."""

    return SimulationState(symbols=symbols, samples_per_symbol=samples_per_symbol, channels=channels)
