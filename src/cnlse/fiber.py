"""Scalar optical fiber propagation.

This module ports the first practical scalar path from the original MATLAB
`Fiber.m`: attenuation, chromatic dispersion, self-phase modulation, and a
basic split-step Fourier method (SSFM).

Original MATLAB implementation credited to Ryan Maher Hamze, Krayem Ahmad
Ajjawi, and Meriam Hawwa.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import exp, isinf, log, pi

import numpy as np

from .constants import C
from .mathutils import fastexp
from .state import SimulationState


@dataclass
class FiberParameters:
    """Fiber parameters matching the original MATLAB structure names."""

    length_m: float
    alpha_db_per_km: float = 0.0
    effective_area_um2: float = 80.0
    nonlinear_index_m2_w: float = 2.7e-20
    wavelength_nm: float = 1550.0
    dispersion_ps_nm_km: float = 17.0
    slope_ps_nm2_km: float = 0.0
    dz_max_m: float | None = None
    dphi_max_rad: float = 3e-3


@dataclass
class FiberResult:
    """Summary from a fiber propagation run."""

    first_step_m: float
    cycles: int
    alpha_linear: float
    beta: np.ndarray
    gamma: np.ndarray


def alpha_db_to_linear(alpha_db_per_km: float) -> float:
    """Convert attenuation from dB/km to 1/m."""

    return log(10) * 1e-4 * alpha_db_per_km


def effective_length(length_m: float, alpha_linear: float) -> float:
    """Return attenuation effective length."""

    if alpha_linear == 0:
        return length_m
    return (1 - exp(-alpha_linear * length_m)) / alpha_linear


def parse_fiber_flag(flag: str) -> tuple[bool, bool, bool, bool]:
    """Return booleans for GVD, PMD, SPM, XPM."""

    valid = {"----", "g---", "--s-", "g-s-", "g-sx", "--sx"}
    if flag.lower() not in valid:
        raise NotImplementedError(f"fiber flag is not implemented yet: {flag}")
    f = flag.lower()
    return f[0] == "g", f[1] == "p", f[2] == "s", f[3] == "x"


def beta_operator(state: SimulationState, fiber: FiberParameters, use_gvd: bool) -> tuple[np.ndarray, np.ndarray]:
    """Build beta(omega) and gamma arrays for scalar propagation."""

    if state.wavelengths_nm is None:
        state.wavelengths_nm = np.array([fiber.wavelength_nm], dtype=float)
    if state.symbol_rate_gbaud is None:
        raise ValueError("state.symbol_rate_gbaud is required for fiber propagation")

    field_x = state.field_x
    if field_x is None:
        raise ValueError("state.field_x has not been created")
    nfc = 1 if field_x.ndim == 1 else field_x.shape[1]

    lambda_ref = fiber.wavelength_nm
    b20 = -lambda_ref**2 / (2 * pi * C) * fiber.dispersion_ps_nm_km * 1e-6
    b30 = (lambda_ref / (2 * pi * C)) ** 2 * (
        2 * lambda_ref * fiber.dispersion_ps_nm_km + lambda_ref**2 * fiber.slope_ps_nm2_km
    ) * 1e-6
    b30 = b30 if use_gvd else 0.0

    maxl = np.max(state.wavelengths_nm)
    minl = np.min(state.wavelengths_nm)
    lamc = 2 * maxl * minl / (maxl + minl)

    if nfc == 1:
        beta1 = np.array([0.0])
        domega_i0 = np.array([2 * pi * C * (1 / lamc - 1 / lambda_ref)])
        gamma = np.array([2 * pi * fiber.nonlinear_index_m2_w / (lamc * fiber.effective_area_um2) * 1e18])
    else:
        domega_i0_all = 2 * pi * C * (1 / state.wavelengths_nm - 1 / lambda_ref)
        domega_ic = 2 * pi * C * (1 / state.wavelengths_nm - 1 / lamc)
        domega_c0 = 2 * pi * C * (1 / lamc - 1 / lambda_ref)
        beta1 = b20 * domega_ic + 0.5 * b30 * (domega_i0_all**2 - domega_c0**2)
        domega_i0 = domega_i0_all
        gamma = 2 * pi * fiber.nonlinear_index_m2_w / (state.wavelengths_nm * fiber.effective_area_um2) * 1e18

    beta2 = (b20 + b30 * domega_i0) if use_gvd else np.zeros_like(domega_i0)
    omega = 2 * pi * state.symbol_rate_gbaud * state.normalized_frequency
    beta = np.zeros((state.samples, nfc), dtype=float)
    for channel in range(nfc):
        beta[:, channel] = omega * beta1[channel] + 0.5 * omega**2 * beta2[channel] + omega**3 * b30 / 6

    if nfc == 1:
        beta = beta[:, 0]

    return beta, gamma


def linear_step(field, beta_times_dz):
    """Apply the linear frequency-domain fiber step."""

    return np.fft.ifft(np.fft.fft(field, axis=0) * fastexp(-beta_times_dz), axis=0)


def nonlinear_step(field, gamma, dz_m: float, alpha_linear: float, spm: bool = True, xpm: bool = False):
    """Apply scalar nonlinear phase rotation."""

    if not spm and not xpm:
        return field

    leff = effective_length(dz_m, alpha_linear)
    values = np.asarray(field, dtype=complex)
    power = np.abs(values) ** 2

    if values.ndim == 1:
        nonlinear_power = power if spm else np.zeros_like(power)
        return values * fastexp(-gamma[0] * nonlinear_power * leff)

    if xpm:
        row_sum = np.sum(power, axis=1)[:, np.newaxis]
        nonlinear_power = 2 * row_sum - power if spm else 2 * (row_sum - power)
    else:
        nonlinear_power = power if spm else np.zeros_like(power)

    return values * fastexp(-gamma[np.newaxis, :] * nonlinear_power * leff)


def next_step(dz_max_m: float, dphi_max_rad: float, gamma, alpha_linear: float, field) -> float:
    """Choose an SSFM step from maximum nonlinear phase rotation."""

    if isinf(dphi_max_rad):
        return dz_max_m

    power_max = np.max(np.abs(field) ** 2, axis=0)
    phase_rate = float(np.max(gamma * power_max))
    if phase_rate <= 0:
        return dz_max_m

    leff = dphi_max_rad / phase_rate
    dl = alpha_linear * leff
    if dl >= 1:
        return dz_max_m
    step = leff if alpha_linear == 0 else -1 / alpha_linear * log(1 - dl)
    return min(step, dz_max_m)


def propagate_scalar(state: SimulationState, fiber: FiberParameters, flag: str = "g-s-") -> FiberResult:
    """Propagate ``state.field_x`` through a scalar fiber."""

    use_gvd, use_pmd, use_spm, use_xpm = parse_fiber_flag(flag)
    if use_pmd:
        raise NotImplementedError("PMD/CNLSE propagation will be added after scalar NLSE support")
    if state.field_x is None:
        raise ValueError("state.field_x has not been created")

    length = fiber.length_m
    dz_max = min(fiber.dz_max_m or length, length)
    alpha_linear = alpha_db_to_linear(fiber.alpha_db_per_km)
    beta, gamma = beta_operator(state, fiber, use_gvd)

    field = np.asarray(state.field_x, dtype=complex)
    dphi = fiber.dphi_max_rad if use_spm or use_xpm else float("inf")
    dz = next_step(dz_max, dphi, gamma, alpha_linear, field)
    first_step = dz
    zdone = 0.0
    cycles = 0

    while zdone < length - 1e-12:
        step = min(dz, length - zdone)
        field = nonlinear_step(field, gamma, step, alpha_linear, spm=use_spm, xpm=use_xpm)
        field = linear_step(field, beta * step)
        field = field * exp(-0.5 * alpha_linear * step)
        zdone += step
        cycles += 1
        if zdone < length:
            dz = next_step(dz_max, dphi, gamma, alpha_linear, field)

    state.field_x = field

    if state.wavelengths_nm is not None and state.dispersion is not None:
        dch = fiber.dispersion_ps_nm_km + fiber.slope_ps_nm2_km * (state.wavelengths_nm - fiber.wavelength_nm)
        state.dispersion = state.dispersion + (dch * length * 1e-3 if use_gvd else 0)

    return FiberResult(
        first_step_m=float(first_step),
        cycles=cycles,
        alpha_linear=float(alpha_linear),
        beta=beta,
        gamma=gamma,
    )
