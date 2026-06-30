"""Power measurement helpers."""

from __future__ import annotations

import numpy as np

from .constants import C
from .mathutils import fastshift
from .state import SimulationState


def average_power(state: SimulationState, channel: int = 1, normalized: bool = False) -> float:
    """Evaluate average power per symbol for a channel.

    This ports the no-filter path of original `Avg_Power.m`.
    """

    if state.field_x is None:
        raise ValueError("state.field_x has not been created")
    if channel < 1 or channel > state.channels:
        raise ValueError("channel does not exist")

    nfft = len(state.normalized_frequency)
    field_x = state.field_x
    if field_x.ndim == 1:
        nfc = 1
    else:
        nfc = field_x.shape[1]

    minfreq = state.normalized_frequency[1] - state.normalized_frequency[0]
    ich = channel - 1

    if nfc != state.channels:
        if state.wavelengths_nm is None or state.symbol_rate_gbaud is None:
            raise ValueError("state wavelengths and symbol rate are required for unique fields")
        maxl = np.max(state.wavelengths_nm)
        minl = np.min(state.wavelengths_nm)
        lamc = 2 * maxl * minl / (maxl + minl)
        deltafn = C * (1 / lamc - 1 / state.wavelengths_nm[ich])
        ndfn = int(round(deltafn / state.symbol_rate_gbaud / minfreq))
        nch = 0
    else:
        ndfn = 0
        nch = ich

    if channel == 1 or nfc == state.channels:
        ndfnl = int(nfft / 2 + ndfn)
    else:
        deltafn = C * (1 / lamc - 1 / state.wavelengths_nm[ich - 1])
        ndfnl_prev = round(deltafn / state.symbol_rate_gbaud / minfreq)
        ndfnl = int(round((ndfn - ndfnl_prev) * 0.5))

    if channel == state.channels or nfc == state.channels:
        ndfnr = int(nfft / 2 - ndfn)
    else:
        deltafn = C * (1 / lamc - 1 / state.wavelengths_nm[ich + 1])
        ndfnr_next = round(deltafn / state.symbol_rate_gbaud / minfreq)
        ndfnr = int(round((ndfnr_next - ndfn) * 0.5))

    nrm = 1.0
    if normalized:
        if state.powers_mw is None:
            raise ValueError("state.powers_mw is required for normalized power")
        nrm = state.powers_mw[ich]

    if field_x.ndim == 1:
        spectrum = np.fft.fft(field_x)
    else:
        spectrum = np.fft.fft(field_x[:, nch])
    spectrum = fastshift(spectrum, ndfn)

    energy = np.sum(np.abs(spectrum[:ndfnl]) ** 2) + np.sum(np.abs(spectrum[nfft - ndfnr :]) ** 2)
    power = energy / nrm / nfft**2

    if state.field_y is not None:
        field_y = state.field_y
        if field_y.ndim == 1:
            spectrum_y = np.fft.fft(field_y)
        else:
            spectrum_y = np.fft.fft(field_y[:, nch])
        spectrum_y = fastshift(spectrum_y, ndfn)
        energy_y = np.sum(np.abs(spectrum_y[:ndfnl]) ** 2) + np.sum(np.abs(spectrum_y[nfft - ndfnr :]) ** 2)
        power += energy_y / nrm / nfft**2

    return float(np.real(power))
