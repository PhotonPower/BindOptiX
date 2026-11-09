"""Strain classification from spot counts (stage Z10, part 1).

Trajectory library: cells are simulated once per layout and protocol with all capture rates 0
(``particles.run(record_visits=True)``). For every cell and spot it holds the first entry time,
the time spent on the spot, and the floor time off the spots before that entry. Captures are
independent events along a path, so for any capture rates k_j of a strain and non-specific
adhesion k_ns the probability that the cell ends on spot j is

    P_j = exp(−k_ns T_off,j) · exp(−Σ_{i entered before j} k_i τ_i) · (1 − exp(−k_j τ_j)),

exact for paths that enter each spot once (one-way coupling, no revisits). The mean over the
library is the capture profile of the strain (per injected cell).

Classifier: the captured counts n_j of one sample follow a multinomial distribution given their
total; the strain with the largest Σ_j n_j log(p_s,j / Σ_k p_s,k) wins (the unknown
concentration drops out). The binding matrix (strain × aptamer capture rates) and the
arrangement of aptamers on the spots are inputs (D3, D13, D15); all values P until M5.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from itertools import product

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.stats import multinomial

from bindoptix.chipsim import ChannelSpots
from bindoptix.particles import RunResult

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class TrajectoryLibrary:
    """Spot visits of simulated cells, columns are grid spots.

    Attributes:
        entry_s, time_s, off_floor_s: (N, S); entry NaN if the spot was never entered.
    """

    entry_s: FloatArray
    time_s: FloatArray
    off_floor_s: FloatArray

    @classmethod
    def from_runs(
        cls, runs: Sequence[RunResult], spot_maps: Sequence[ChannelSpots], n_grid: int
    ) -> TrajectoryLibrary:
        rows_e, rows_t, rows_o = [], [], []
        for r, cs in zip(runs, spot_maps, strict=True):
            n = r.position_um.shape[0]
            if r.visit_entry_s.shape != (n, len(cs.spots)):
                raise ValueError("runs must be recorded with record_visits=True")
            e = np.full((n, n_grid), math.nan)
            t = np.zeros((n, n_grid))
            o = np.full((n, n_grid), math.nan)
            for local, g in enumerate(cs.grid_index):
                e[:, g] = r.visit_entry_s[:, local]
                t[:, g] = r.visit_time_s[:, local]
                o[:, g] = r.visit_off_floor_s[:, local]
            rows_e.append(e)
            rows_t.append(t)
            rows_o.append(o)
        return cls(np.vstack(rows_e), np.vstack(rows_t), np.vstack(rows_o))

    @property
    def n_cells(self) -> int:
        return int(self.entry_s.shape[0])

    def capture_probabilities(self, rates_s: ArrayLike, k_ns_s: float = 0.0) -> FloatArray:
        """P_j per cell (N, S) for capture rates per grid spot and non-specific adhesion."""
        k = np.asarray(rates_s, dtype=np.float64)
        n, s = self.entry_s.shape
        if k.shape != (s,) or np.any(k < 0) or k_ns_s < 0:
            raise ValueError("one non-negative rate per grid spot and k_ns >= 0 required")
        entered = ~np.isnan(self.entry_s)
        order = np.argsort(np.where(entered, self.entry_s, np.inf), axis=1)
        rows = np.arange(n)[:, None]
        tau = self.time_s[rows, order]
        ent = entered[rows, order]
        off = np.nan_to_num(self.off_floor_s[rows, order], nan=0.0)
        kk = k[order]
        with np.errstate(invalid="ignore"):
            lost_ns = np.where(off > 0.0, k_ns_s * off, 0.0)
        exposure = np.where(ent, kk * tau, 0.0)
        before = np.cumsum(exposure, axis=1) - exposure
        p_sorted = np.where(ent, np.exp(-lost_ns - before) * -np.expm1(-exposure), 0.0)
        out = np.zeros_like(p_sorted)
        out[rows, order] = p_sorted
        return out

    def profile(self, rates_s: ArrayLike, k_ns_s: float = 0.0) -> FloatArray:
        """Expected captures per grid spot per injected cell (S,)."""
        return np.asarray(self.capture_probabilities(rates_s, k_ns_s).mean(axis=0))


@dataclass(frozen=True)
class BindingMatrix:
    """Capture rates [1/s] of strain s on aptamer a, and the aptamer on each grid spot."""

    rates_s: FloatArray  # (strains, aptamers)
    strains: tuple[str, ...]
    arrangement: tuple[int, ...]  # aptamer index per grid spot

    def __post_init__(self) -> None:
        r = np.asarray(self.rates_s, dtype=np.float64)
        if r.ndim != 2 or r.shape[0] != len(self.strains) or np.any(r < 0):
            raise ValueError("rates_s must be (strains, aptamers) and non-negative")
        if any(not 0 <= a < r.shape[1] for a in self.arrangement):
            raise ValueError("arrangement refers to unknown aptamers")
        object.__setattr__(self, "rates_s", r)

    def spot_rates(self, strain: int) -> FloatArray:
        return np.asarray(self.rates_s[strain, list(self.arrangement)])

    def profiles(self, lib: TrajectoryLibrary, k_ns_s: float = 0.0) -> FloatArray:
        """Capture profiles (strains, S) per injected cell."""
        return np.vstack(
            [lib.profile(self.spot_rates(s), k_ns_s) for s in range(len(self.strains))]
        )


def classify(counts: ArrayLike, profiles: ArrayLike, min_total: int = 1) -> int:
    """Strain index with the largest multinomial log-likelihood; −1 below ``min_total``."""
    n = np.asarray(counts, dtype=np.float64)
    p = np.asarray(profiles, dtype=np.float64)
    if n.sum() < min_total:
        return -1
    q = p / np.maximum(p.sum(axis=1, keepdims=True), 1e-300)
    with np.errstate(divide="ignore"):
        ll = np.where(n[None, :] > 0, n[None, :] * np.log(q), 0.0).sum(axis=1)
    return int(np.argmax(ll))


def confusion_matrix(
    profiles: ArrayLike, n_cells: int, trials: int, seed: int, min_total: int = 1
) -> FloatArray:
    """Monte Carlo: rows true strain, columns decision (last column: no decision)."""
    p = np.asarray(profiles, dtype=np.float64)
    ns, s = p.shape
    g = np.random.default_rng(seed)
    out = np.zeros((ns, ns + 1))
    for a in range(ns):
        probs = np.append(p[a], max(0.0, 1.0 - p[a].sum()))
        draws = g.multinomial(n_cells, probs, size=trials)[:, :s]
        for d in draws:
            out[a, classify(d, p, min_total)] += 1
    return out / trials


def confusion_given_total(profiles: ArrayLike, total: int) -> FloatArray:
    """Exact confusion when exactly ``total`` cells are captured (enumeration, small S)."""
    p = np.asarray(profiles, dtype=np.float64)
    ns, s = p.shape
    q = p / p.sum(axis=1, keepdims=True)
    out = np.zeros((ns, ns))
    for counts in product(range(total + 1), repeat=s):
        if sum(counts) != total:
            continue
        dec = classify(counts, p)
        for a in range(ns):
            out[a, dec] += multinomial.pmf(counts, total, q[a])
    return out
