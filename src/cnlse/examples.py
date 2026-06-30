"""Runnable examples for the cnlse package."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .electric import ook_signal
from .fields import create_field
from .fiber import FiberParameters, propagate_scalar
from .laser import laser_source
from .link import phi_to_power
from .modulators import mach_zehnder
from .patterns import de_bruijn_binary
from .power import average_power
from .plotting import plot_power_traces
from .state import create_state


def run_ook_100km_compensated(output: str | Path = "ook_100km_compensated.png") -> Path:
    """Create the first OOK fiber-compensation demonstration plot."""

    bits = de_bruijn_binary(3)[:12]
    samples_per_symbol = 120
    state = create_state(symbols=len(bits), samples_per_symbol=samples_per_symbol, channels=1)
    state.symbol_rate_gbaud = 10
    time = np.linspace(0, len(bits), len(bits) * samples_per_symbol, endpoint=False)

    drive = ook_signal(bits, samples_per_symbol, pulse="cosroll", duty=1.0, roll=0.2)
    gamma = 2 * np.pi * 2.7e-20 / (1550 * 80) * 1e21
    peak_power = phi_to_power(0.4 * np.pi, 1e5, 0.2, gamma, 0, 1)
    laser, wavelengths, powers = laser_source(peak_power, samples=len(drive), channels=1, wavelength_nm=1550)
    optical_field = mach_zehnder(laser, drive, extinction_ratio_db=13)
    state.wavelengths_nm = wavelengths
    state.powers_mw = powers
    create_field(state, optical_field, mode="unique")

    input_power = np.abs(state.field_x) ** 2
    tx_fiber = FiberParameters(
        length_m=1e5,
        alpha_db_per_km=0.2,
        effective_area_um2=80,
        nonlinear_index_m2_w=2.7e-20,
        wavelength_nm=1550,
        dispersion_ps_nm_km=17,
        slope_ps_nm2_km=0,
        dphi_max_rad=3e-3,
        dz_max_m=2e4,
    )
    comp_fiber = FiberParameters(
        length_m=(-tx_fiber.dispersion_ps_nm_km * tx_fiber.length_m / 1e3) / -100 * 1e3,
        alpha_db_per_km=0.6,
        effective_area_um2=20,
        nonlinear_index_m2_w=2.7e-20,
        wavelength_nm=1550,
        dispersion_ps_nm_km=-100,
        slope_ps_nm2_km=0,
        dphi_max_rad=3e-3,
        dz_max_m=2e4,
    )

    uncomp_state = create_state(symbols=len(bits), samples_per_symbol=samples_per_symbol, channels=1)
    uncomp_state.symbol_rate_gbaud = state.symbol_rate_gbaud
    uncomp_state.wavelengths_nm = wavelengths.copy()
    uncomp_state.powers_mw = powers.copy()
    create_field(uncomp_state, optical_field, mode="unique")
    propagate_scalar(uncomp_state, tx_fiber, flag="g-s-")
    uncompensated = np.abs(uncomp_state.field_x) ** 2

    comp_state = create_state(symbols=len(bits), samples_per_symbol=samples_per_symbol, channels=1)
    comp_state.symbol_rate_gbaud = state.symbol_rate_gbaud
    comp_state.wavelengths_nm = wavelengths.copy()
    comp_state.powers_mw = powers.copy()
    create_field(comp_state, optical_field, mode="unique")
    propagate_scalar(comp_state, tx_fiber, flag="g-s-")
    propagate_scalar(comp_state, comp_fiber, flag="g-s-")
    compensated = np.abs(comp_state.field_x) ** 2
    measured_power = average_power(state, 1)

    traces = [
        ("Input OOK optical pulse train", input_power),
        ("Output after 100 km fiber without compensation", uncompensated),
        ("Output after dispersion-compensating fiber", compensated),
    ]
    return plot_power_traces(
        time,
        traces,
        output,
        title="OOK Fiber Simulation With Dispersion Compensation",
        subtitle=f"Measured launch average power: {measured_power:.3f} mW",
    )


def main_example() -> None:
    """Run the packaged example from the command line."""

    output = run_ook_100km_compensated()
    print(f"Wrote {output}")
