"""Streamlit interface for the cnlse demo tool."""

from __future__ import annotations

from io import BytesIO

import numpy as np

from .electric import ook_signal
from .fields import create_field
from .fiber import FiberParameters, propagate_scalar
from .laser import laser_source
from .link import phi_to_power
from .modulators import mach_zehnder
from .patterns import de_bruijn_binary, random_bits
from .state import create_state


def _require_streamlit():
    try:
        import matplotlib.pyplot as plt
        import streamlit as st
    except Exception as exc:  # pragma: no cover - exercised by real app install
        raise RuntimeError('Install app dependencies with: python -m pip install -e ".[app]"') from exc
    return st, plt


def build_simulation(
    symbols: int,
    samples_per_symbol: int,
    symbol_rate_gbaud: float,
    pattern_type: str,
    seed: int,
    phi_rad: float,
    extinction_ratio_db: float,
    wavelength_nm: float,
    tx_length_km: float,
    tx_loss_db_km: float,
    tx_dispersion: float,
    comp_loss_db_km: float,
    comp_dispersion: float,
    nonlinear: bool,
):
    """Run the OOK compensation simulation and return time/traces."""

    if pattern_type == "De Bruijn":
        bits = de_bruijn_binary(max(1, int(np.ceil(np.log2(symbols)))))
        bits = (bits * ((symbols // len(bits)) + 1))[:symbols]
    else:
        bits = random_bits(symbols, seed=seed)

    state = create_state(symbols=symbols, samples_per_symbol=samples_per_symbol, channels=1)
    state.symbol_rate_gbaud = symbol_rate_gbaud
    time = np.linspace(0, symbols, symbols * samples_per_symbol, endpoint=False)

    drive = ook_signal(bits, samples_per_symbol, pulse="cosroll", duty=1.0, roll=0.2)
    gamma = 2 * np.pi * 2.7e-20 / (wavelength_nm * 80) * 1e21
    peak_power = phi_to_power(phi_rad, tx_length_km * 1000, tx_loss_db_km, gamma, 0, 1)
    laser, wavelengths, powers = laser_source(peak_power, samples=len(drive), channels=1, wavelength_nm=wavelength_nm)
    optical_field = mach_zehnder(laser, drive, extinction_ratio_db=extinction_ratio_db)

    tx_fiber = FiberParameters(
        length_m=tx_length_km * 1000,
        alpha_db_per_km=tx_loss_db_km,
        effective_area_um2=80,
        nonlinear_index_m2_w=2.7e-20,
        wavelength_nm=wavelength_nm,
        dispersion_ps_nm_km=tx_dispersion,
        slope_ps_nm2_km=0,
        dphi_max_rad=3e-3,
        dz_max_m=20_000,
    )
    comp_length_m = (-tx_fiber.dispersion_ps_nm_km * tx_fiber.length_m / 1000) / comp_dispersion * 1000
    comp_fiber = FiberParameters(
        length_m=max(0, comp_length_m),
        alpha_db_per_km=comp_loss_db_km,
        effective_area_um2=20,
        nonlinear_index_m2_w=2.7e-20,
        wavelength_nm=wavelength_nm,
        dispersion_ps_nm_km=comp_dispersion,
        slope_ps_nm2_km=0,
        dphi_max_rad=3e-3,
        dz_max_m=20_000,
    )
    flag = "g-s-" if nonlinear else "g---"

    def fresh_state():
        sim = create_state(symbols=symbols, samples_per_symbol=samples_per_symbol, channels=1)
        sim.symbol_rate_gbaud = symbol_rate_gbaud
        sim.wavelengths_nm = wavelengths.copy()
        sim.powers_mw = powers.copy()
        create_field(sim, optical_field, mode="unique")
        return sim

    input_power = np.abs(optical_field) ** 2
    uncomp = fresh_state()
    tx_result = propagate_scalar(uncomp, tx_fiber, flag=flag)
    uncomp_power = np.abs(uncomp.field_x) ** 2

    comp = fresh_state()
    propagate_scalar(comp, tx_fiber, flag=flag)
    comp_result = propagate_scalar(comp, comp_fiber, flag=flag)
    comp_power = np.abs(comp.field_x) ** 2

    return {
        "bits": bits,
        "time": time,
        "input_power": input_power,
        "uncompensated_power": uncomp_power,
        "compensated_power": comp_power,
        "peak_power": peak_power,
        "tx_result": tx_result,
        "comp_result": comp_result,
        "comp_length_km": comp_fiber.length_m / 1000,
    }


def main() -> None:  # pragma: no cover - run by streamlit
    st, plt = _require_streamlit()

    st.set_page_config(page_title="CNLSE Fiber Optics Simulator", layout="wide")
    st.title("Fiber Optics Modulation Using Coupled Nonlinear Schrodinger Equation")
    st.caption("Open-source educational simulator based on the original 2009 MATLAB senior project.")

    with st.sidebar:
        st.header("Signal")
        symbols = st.slider("Number of bits", 4, 128, 16, help="How many digital symbols to simulate.")
        samples_per_symbol = st.slider("Detail per bit", 16, 256, 96, step=16, help="More samples makes smoother plots but runs slower.")
        symbol_rate_gbaud = st.slider("Symbol rate", 1.0, 40.0, 10.0, step=1.0, help="Transmission speed in billions of symbols per second.")
        pattern_type = st.selectbox("Bit pattern", ["De Bruijn", "Random"], help="De Bruijn gives compact transitions; Random creates seeded random bits.")
        seed = st.number_input("Random seed", value=7, step=1, help="Used only for random bit patterns.")

        st.header("Source And Modulator")
        wavelength_nm = st.slider("Laser wavelength", 1260.0, 1625.0, 1550.0, step=1.0, help="Optical carrier wavelength in nanometers.")
        phi_rad = st.slider("Nonlinear phase target", 0.01, float(np.pi), 0.4 * float(np.pi), help="Target cumulative nonlinear phase used to choose launch power.")
        extinction_ratio_db = st.slider("Extinction ratio", 3.0, 30.0, 13.0, help="Mach-Zehnder on/off contrast in dB.")

        st.header("Transmit Fiber")
        tx_length_km = st.slider("Fiber length", 1.0, 200.0, 100.0, step=1.0, help="Transmit fiber length in kilometers.")
        tx_loss_db_km = st.slider("Fiber loss", 0.0, 1.0, 0.2, step=0.05, help="Power loss per kilometer.")
        tx_dispersion = st.slider("Fiber dispersion", -50.0, 50.0, 17.0, step=1.0, help="Chromatic dispersion in ps/nm/km.")

        st.header("Compensating Fiber")
        comp_loss_db_km = st.slider("Compensating fiber loss", 0.0, 2.0, 0.6, step=0.05, help="Power loss in the compensation fiber.")
        comp_dispersion = st.slider("Compensating dispersion", -200.0, -1.0, -100.0, step=1.0, help="Negative dispersion used to cancel transmit fiber dispersion.")
        nonlinear = st.checkbox("Include nonlinear SPM", value=True, help="Turns self-phase modulation on or off.")

    data = build_simulation(
        symbols=symbols,
        samples_per_symbol=samples_per_symbol,
        symbol_rate_gbaud=symbol_rate_gbaud,
        pattern_type=pattern_type,
        seed=int(seed),
        phi_rad=phi_rad,
        extinction_ratio_db=extinction_ratio_db,
        wavelength_nm=wavelength_nm,
        tx_length_km=tx_length_km,
        tx_loss_db_km=tx_loss_db_km,
        tx_dispersion=tx_dispersion,
        comp_loss_db_km=comp_loss_db_km,
        comp_dispersion=comp_dispersion,
        nonlinear=nonlinear,
    )

    metric_cols = st.columns(4)
    metric_cols[0].metric("Launch peak power", f"{data['peak_power']:.3f} mW")
    metric_cols[1].metric("Comp fiber length", f"{data['comp_length_km']:.2f} km")
    metric_cols[2].metric("TX SSFM steps", data["tx_result"].cycles)
    metric_cols[3].metric("Comp SSFM steps", data["comp_result"].cycles)

    fig, axes = plt.subplots(3, 1, figsize=(12, 7), sharex=True)
    traces = [
        ("Input OOK optical pulse train", data["input_power"]),
        ("Output after transmit fiber", data["uncompensated_power"]),
        ("Output after compensation", data["compensated_power"]),
    ]
    for axis, (label, values) in zip(axes, traces):
        axis.plot(data["time"], values, linewidth=2)
        axis.fill_between(data["time"], 0, values, alpha=0.13)
        axis.set_title(label, loc="left")
        axis.set_ylabel("Power")
        axis.grid(True, alpha=0.25)
    axes[-1].set_xlabel("Time / symbol slots")
    fig.tight_layout()
    st.pyplot(fig, clear_figure=False)

    buffer = BytesIO()
    fig.savefig(buffer, format="png", dpi=160)
    st.download_button("Export Plot", buffer.getvalue(), file_name="cnlse_simulation.png", mime="image/png")

    with st.expander("Bit pattern"):
        st.code("".join(str(bit) for bit in data["bits"]))


if __name__ == "__main__":
    main()
