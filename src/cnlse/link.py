"""Optical link helper equations."""

from __future__ import annotations

import numpy as np


def phi_to_power(phi, length, alpha_db, gamma, gain_db, spans: int = 1):
    """Convert cumulative nonlinear phase rotation into transmit power.

    This ports the core equation from the original MATLAB `Phi2pow.m`.

    Parameters
    ----------
    phi:
        Desired cumulative nonlinear phase rotation in radians.
    length:
        Fiber length or lengths in meters.
    alpha_db:
        Attenuation in dB/km for each fiber.
    gamma:
        Nonlinear coefficient in 1/W/m for each fiber.
    gain_db:
        Net gain in dB before each fiber.
    spans:
        Number of repeated spans.

    Returns
    -------
    float
        Transmit power in mW.
    """

    if spans <= 0:
        raise ValueError("spans must be positive")

    length = np.atleast_1d(np.asarray(length, dtype=float))
    alpha_db = np.atleast_1d(np.asarray(alpha_db, dtype=float))
    gamma = np.atleast_1d(np.asarray(gamma, dtype=float))
    gain_db = np.atleast_1d(np.asarray(gain_db, dtype=float))

    if not (len(length) == len(alpha_db) == len(gamma) == len(gain_db)):
        raise ValueError("length, alpha_db, gamma, and gain_db must have the same length")

    alpha_linear = np.log(10) * 1e-4 * alpha_db
    effective_length = np.where(
        alpha_linear == 0,
        length,
        (1 - np.exp(-alpha_linear * length)) / alpha_linear,
    )

    loss_db = -alpha_db * length * 1e-3
    prior_loss = np.concatenate(([0.0], loss_db[:-1]))
    net_gain = prior_loss + gain_db
    cumulative_gain = 10 ** (np.cumsum(net_gain) * 0.1)

    denominator = np.sum(effective_length * gamma * cumulative_gain * spans)
    if denominator == 0:
        raise ValueError("nonlinear phase denominator is zero")

    return float(phi / denominator * 1e3)
