"""Bit-pattern helpers inspired by the original MATLAB Pattern.m module."""

from __future__ import annotations

import random


def random_bits(length: int, seed: int | None = None) -> list[int]:
    """Return a reproducible random binary sequence."""

    if length <= 0:
        raise ValueError("length must be positive")

    rng = random.Random(seed)
    return [rng.randrange(2) for _ in range(length)]


def de_bruijn_binary(order: int) -> list[int]:
    """Return a binary De Bruijn sequence of subsequence length ``order``.

    The original MATLAB examples use De Bruijn patterns for compact test
    sequences. This implementation returns one cyclic sequence as a list of
    0/1 integers.
    """

    if order <= 0:
        raise ValueError("order must be positive")

    alphabet_size = 2
    sequence: list[int] = []
    a = [0] * (alphabet_size * order)

    def db(t: int, p: int) -> None:
        if t > order:
            if order % p == 0:
                sequence.extend(a[1 : p + 1])
            return

        a[t] = a[t - p]
        db(t + 1, p)
        for value in range(a[t - p] + 1, alphabet_size):
            a[t] = value
            db(t + 1, t)

    db(1, 1)
    return sequence
