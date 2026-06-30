import numpy as np

from cnlse.fields import create_field
from cnlse.fiber import FiberParameters, alpha_db_to_linear, effective_length, propagate_scalar
from cnlse.laser import laser_source
from cnlse.power import average_power
from cnlse.state import create_state


def test_attenuation_conversion_and_effective_length():
    alpha = alpha_db_to_linear(0.2)
    assert alpha > 0
    assert effective_length(100_000, 0) == 100_000
    assert effective_length(100_000, alpha) < 100_000


def test_attenuation_only_reduces_field_power():
    state = create_state(8, 16, 1)
    state.symbol_rate_gbaud = 10
    laser, wavelengths, powers = laser_source(1.0, state.samples, channels=1, wavelength_nm=1550)
    state.wavelengths_nm = wavelengths
    state.powers_mw = powers
    create_field(state, laser, mode="unique")

    before = average_power(state, 1)
    result = propagate_scalar(state, FiberParameters(length_m=100_000, alpha_db_per_km=0.2), flag="----")
    after = average_power(state, 1)

    assert result.cycles == 1
    assert after < before
    np.testing.assert_allclose(after / before, 10 ** (-0.2 * 100 / 10), rtol=1e-6)


def test_linear_dispersion_preserves_power_when_lossless():
    state = create_state(8, 32, 1)
    state.symbol_rate_gbaud = 10
    rng = np.random.default_rng(7)
    field = rng.normal(size=state.samples) + 1j * rng.normal(size=state.samples)
    state.wavelengths_nm = np.array([1550.0])
    state.powers_mw = np.array([1.0])
    create_field(state, field, mode="unique")

    before = average_power(state, 1)
    propagate_scalar(
        state,
        FiberParameters(length_m=10_000, alpha_db_per_km=0.0, dispersion_ps_nm_km=17.0),
        flag="g---",
    )
    after = average_power(state, 1)

    np.testing.assert_allclose(after, before, rtol=1e-10)


def test_spm_changes_phase_but_preserves_power_when_lossless():
    state = create_state(8, 16, 1)
    state.symbol_rate_gbaud = 10
    laser, wavelengths, powers = laser_source(2.0, state.samples, channels=1, wavelength_nm=1550)
    state.wavelengths_nm = wavelengths
    state.powers_mw = powers
    create_field(state, laser, mode="unique")

    before_power = average_power(state, 1)
    before_field = state.field_x.copy()
    propagate_scalar(
        state,
        FiberParameters(length_m=1_000, alpha_db_per_km=0.0, dphi_max_rad=0.01),
        flag="--s-",
    )

    np.testing.assert_allclose(average_power(state, 1), before_power, rtol=1e-10)
    assert np.max(np.abs(state.field_x - before_field)) > 0
