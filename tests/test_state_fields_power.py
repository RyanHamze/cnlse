import numpy as np

from cnlse.electric import ook_signal
from cnlse.fields import create_field
from cnlse.laser import laser_source
from cnlse.modulators import mach_zehnder
from cnlse.power import average_power
from cnlse.state import create_state


def test_create_state_frequency_grid():
    state = create_state(4, 8, 1)
    assert state.samples == 32
    assert len(state.normalized_frequency) == 32
    assert state.normalized_frequency[1] - state.normalized_frequency[0] == 0.25


def test_create_field_sepfields_and_average_power():
    state = create_state(4, 8, 1)
    state.symbol_rate_gbaud = 10
    field = np.ones(state.samples, dtype=complex) * 2
    state.wavelengths_nm = np.array([1550.0])
    state.powers_mw = np.array([4.0])

    create_field(state, field, mode="sepfields")

    assert state.field_x.shape == (32, 1)
    assert state.delay.tolist() == [[0]]
    assert average_power(state, 1) == 4.0
    assert average_power(state, 1, normalized=True) == 1.0


def test_create_field_unique_single_channel_chain():
    state = create_state(8, 32, 1)
    state.symbol_rate_gbaud = 10

    drive = ook_signal([1, 0, 1, 1, 0, 0, 1, 0], state.samples_per_symbol, pulse="cosroll")
    laser, wavelengths, powers = laser_source(2.0, state.samples, channels=1, wavelength_nm=1550)
    optical = mach_zehnder(laser, drive, extinction_ratio_db=13)

    state.wavelengths_nm = wavelengths
    state.powers_mw = powers
    create_field(state, optical, mode="unique")

    assert state.field_x.shape == (state.samples,)
    assert average_power(state, 1) > 0
