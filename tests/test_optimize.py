"""Stage Z10, part 2: design optimisation (T26)."""

import math

import numpy as np
import pytest

from bindoptix.classify import TrajectoryLibrary
from bindoptix.optimize import (
    brute_force_arrangement,
    information_per_cell,
    optimize_arrangement,
)


def synthetic_library(n: int, s: int, seed: int) -> TrajectoryLibrary:
    """Cells that slide over a random subset of spots in increasing spot order."""
    g = np.random.default_rng(seed)
    entry = np.full((n, s), math.nan)
    time = np.zeros((n, s))
    off = np.full((n, s), math.nan)
    for c in range(n):
        first = int(g.integers(0, s))
        t = g.uniform(0, 100)
        for j in range(first, s):
            entry[c, j] = t
            time[c, j] = g.uniform(20, 200)
            off[c, j] = 0.0
            t += time[c, j] + g.uniform(0, 50)
    return TrajectoryLibrary(entry, time, off)


def test_information_limits() -> None:
    # Perfect: every cell captured on its strain's own spot → log2(number of strains).
    assert information_per_cell(np.eye(4)) == pytest.approx(2.0, rel=1e-12)
    # Nothing captured, or identical profiles → no information.
    assert information_per_cell(np.zeros((3, 5))) == pytest.approx(0.0, abs=1e-15)
    assert information_per_cell(np.tile([0.1, 0.2, 0.05], (3, 1))) == pytest.approx(0.0, abs=1e-12)
    # Half of the cells captured, perfectly sorted: 1 bit · 0.5 for two strains.
    assert information_per_cell(0.5 * np.eye(2)) == pytest.approx(0.5, rel=1e-12)


@pytest.mark.parametrize("seed", [1, 2, 3, 4])
def test_t26_search_finds_brute_force_optimum(seed: int) -> None:
    # Known optimum by enumeration (5 spots, 4 aptamers with one repeated → 60 arrangements).
    lib = synthetic_library(400, 5, seed)
    g = np.random.default_rng(seed + 10)
    rates = 0.003 * g.lognormal(0, 1.0, (4, 4))
    np.fill_diagonal(rates, 0.02)
    aptamers = (0, 1, 2, 3, 0)
    exact = brute_force_arrangement(lib, rates, aptamers)
    found = optimize_arrangement(lib, rates, aptamers, seed, restarts=6)
    assert found.information_bits == pytest.approx(exact.information_bits, rel=1e-9)
    assert exact.evaluations == 60
