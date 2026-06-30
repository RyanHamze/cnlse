import math

import numpy as np

from cnlse.link import phi_to_power
from cnlse.mathutils import fastexp, fastshift, nmod


def test_nmod_scalar_examples():
    assert nmod(-2, 8) == 6
    assert nmod(-1, 8) == 7
    assert nmod(0, 8) == 8
    assert nmod(1, 8) == 1
    assert nmod(9, 8) == 1


def test_nmod_array():
    np.testing.assert_array_equal(nmod(np.array([-2, -1, 0, 1, 9]), 8), np.array([6, 7, 8, 1, 1]))


def test_fastshift_vector():
    values = np.arange(11)
    np.testing.assert_array_equal(fastshift(values, 2), np.array([9, 10, 0, 1, 2, 3, 4, 5, 6, 7, 8]))
    np.testing.assert_array_equal(fastshift(values, -2), np.array([2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 1]))


def test_fastshift_matrix_rows():
    values = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    np.testing.assert_array_equal(fastshift(values, -1), np.array([[4, 5, 6], [7, 8, 9], [1, 2, 3]]))


def test_fastexp():
    assert fastexp(0) == complex(1, 0)
    assert abs(fastexp(math.pi) + 1) < 1e-12


def test_phi_to_power_single_span_matches_formula_shape():
    phi = 0.4 * math.pi
    length = 1e5
    alpha = 0.2
    gamma = 2 * math.pi * 2.7e-20 / (1550 * 80) * 1e21
    value = phi_to_power(phi, length, alpha, gamma, 0, 1)
    assert value > 0
    assert value < 100
