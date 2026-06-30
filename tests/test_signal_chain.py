import math

import numpy as np

from cnlse.electric import ook_signal, pulse_shape
from cnlse.laser import channel_wavelengths, laser_source
from cnlse.link import phi_to_power
from cnlse.modulators import mach_zehnder
from cnlse.patterns import de_bruijn_binary


def test_pulse_shape_lengths():
    assert len(pulse_shape("cosroll", 64, duty=1.0, roll=0.2)) == 128
    assert len(pulse_shape("idpulse", 64, duty=1.0)) == 128
    assert len(pulse_shape("dirac", 64)) == 128


def test_ook_signal_shape_and_levels():
    signal = ook_signal([1, 0, 1], 16, pulse="idpulse", duty=1.0)
    assert signal.shape == (48,)
    assert signal.max() == 1.0
    assert signal.min() == 0.0


def test_laser_source_single_channel():
    field, wavelengths, powers = laser_source(2.5, samples=10, channels=1, wavelength_nm=1550)
    assert field.shape == (10,)
    assert wavelengths.tolist() == [1550.0]
    assert powers.tolist() == [2.5]
    np.testing.assert_allclose(np.abs(field) ** 2, np.full(10, 2.5))


def test_channel_wavelengths():
    np.testing.assert_allclose(channel_wavelengths(1550, 0.8, 4), np.array([1548.8, 1549.6, 1550.4, 1551.2]))


def test_mach_zehnder_runs_on_ex01_like_chain():
    bits = de_bruijn_binary(3)
    drive = ook_signal(bits, 32, pulse="cosroll", duty=1.0, roll=0.2)

    gamma = 2 * math.pi * 2.7e-20 / (1550 * 80) * 1e21
    peak_power = phi_to_power(0.4 * math.pi, 1e5, 0.2, gamma, 0, 1)
    laser, _, _ = laser_source(peak_power, samples=len(drive), channels=1, wavelength_nm=1550)

    modulated = mach_zehnder(laser, drive, extinction_ratio_db=13)
    assert modulated.shape == drive.shape
    assert np.iscomplexobj(modulated)
    assert np.max(np.abs(modulated) ** 2) > np.min(np.abs(modulated) ** 2)
