"""Optical field creation and multiplexing."""

from __future__ import annotations

import numpy as np

from .constants import C
from .mathutils import fastshift
from .state import SimulationState


def _as_channel_matrix(values, samples: int, channels: int) -> np.ndarray:
    array = np.asarray(values, dtype=complex)
    if array.ndim == 1:
        array = array[:, np.newaxis]
    if array.shape != (samples, channels):
        raise ValueError(f"field must have shape {(samples, channels)} or ({samples},)")
    return array


def create_field(
    state: SimulationState,
    field_x,
    field_y=None,
    mode: str = "unique",
    delay=None,
    power_mode: str = "peak",
    rng: np.random.Generator | None = None,
) -> SimulationState:
    """Create the optical field carried by ``state``.

    Parameters follow the original `Create_Field.m` concepts:

    - `mode="sepfields"` keeps one column per channel.
    - `mode="unique"` multiplexes all channels into one combined field.
    - `power_mode="average"` rescales fields so `state.powers_mw` is treated
      as desired average power.
    """

    mode = mode.lower()
    if mode not in {"unique", "sepfields"}:
        raise ValueError("mode must be 'unique' or 'sepfields'")

    sigx = _as_channel_matrix(field_x, state.samples, state.channels)
    sigy = None if field_y is None else _as_channel_matrix(field_y, state.samples, state.channels)
    npol = 2 if sigy is not None else 1

    if state.powers_mw is not None and power_mode.lower() == "average":
        if sigy is None:
            average_energy = np.mean(np.abs(sigx) ** 2, axis=0)
        else:
            average_energy = np.mean(np.abs(sigx) ** 2 + np.abs(sigy) ** 2, axis=0)

        scale = np.sqrt(state.powers_mw / average_energy)
        sigx = sigx * scale[np.newaxis, :]
        if sigy is not None:
            sigy = sigy * scale[np.newaxis, :]
        state.powers_mw = state.powers_mw * state.powers_mw / average_energy

    if delay is None:
        tau = np.zeros((npol, state.channels), dtype=int)
    elif isinstance(delay, str) and delay.lower() == "rand":
        rng = rng or np.random.default_rng()
        tau = np.rint(rng.random((npol, state.channels)) * state.samples_per_symbol).astype(int)
    else:
        tau = np.rint(np.asarray(delay, dtype=float) * state.samples_per_symbol).astype(int)
        if tau.shape != (npol, state.channels):
            raise ValueError(f"delay must have shape {(npol, state.channels)}")

    for channel in range(state.channels):
        sigx[:, channel] = fastshift(sigx[:, channel], int(tau[0, channel]))
        if sigy is not None:
            sigy[:, channel] = fastshift(sigy[:, channel], int(tau[1, channel]))

    state.delay = tau
    state.dispersion = np.zeros((npol, state.channels), dtype=float)

    if mode == "sepfields":
        state.field_x = sigx.copy()
        state.field_x_tx = sigx.copy()
        if sigy is not None:
            state.field_y = sigy.copy()
            state.field_y_tx = sigy.copy()
        return state

    if state.wavelengths_nm is None:
        if state.channels == 1:
            state.wavelengths_nm = np.array([1550.0])
        else:
            raise ValueError("state.wavelengths_nm is required for unique multi-channel fields")
    if state.symbol_rate_gbaud is None:
        raise ValueError("state.symbol_rate_gbaud is required for unique fields")

    lamt = state.wavelengths_nm
    maxl = np.max(lamt)
    minl = np.min(lamt)
    lamc = 2 * maxl * minl / (maxl + minl)
    deltafn = C * (1 / lamc - 1 / lamt)
    minfreq = state.normalized_frequency[1] - state.normalized_frequency[0]
    ndfn = np.rint(deltafn / state.symbol_rate_gbaud / minfreq).astype(int)

    field_x = np.zeros(state.samples, dtype=complex)
    zfieldx = np.fft.fft(sigx, axis=0)
    for channel in range(state.channels):
        field_x += fastshift(zfieldx[:, channel], int(-ndfn[channel]))
    state.field_x = np.fft.ifft(field_x)
    state.field_x_tx = state.field_x.copy()

    if sigy is not None:
        field_y = np.zeros(state.samples, dtype=complex)
        zfieldy = np.fft.fft(sigy, axis=0)
        for channel in range(state.channels):
            field_y += fastshift(zfieldy[:, channel], int(-ndfn[channel]))
        state.field_y = np.fft.ifft(field_y)
        state.field_y_tx = state.field_y.copy()

    return state
