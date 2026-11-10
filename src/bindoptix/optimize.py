"""Design optimisation (stage Z10, part 2).

Objective: information per injected cell, the mutual information between the strain (uniform
prior) and the outcome of one cell (captured on spot j, or not captured), in bits. It combines
how many cells are caught and how well their spots tell the strains apart, and it is exact and
smooth in the profiles (no Monte Carlo noise), so local search is deterministic for a seed.

- :func:`optimize_arrangement`: aptamers on the grid spots by pair swaps (hill climbing with
  random restarts); the multiset of aptamers is kept. Uses a trajectory library, so one
  evaluation costs only the profiles.
- :func:`brute_force_arrangement`: exhaustive search for small problems (tests).
- Continuous parameters (flow rate, channel height) need new libraries per value; they are
  swept by the caller (``experiments/z10_optimierung.py``).

D11 named BoTorch; this first version uses no extra dependency.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from itertools import permutations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from bindoptix.classify import BindingMatrix, TrajectoryLibrary

FloatArray = NDArray[np.float64]


def information_per_cell(profiles: ArrayLike) -> float:
    """I(strain; outcome) in bits for one injected cell, uniform prior over the strains."""
    p = np.asarray(profiles, dtype=np.float64)
    ns = p.shape[0]
    joint = np.column_stack([p, np.clip(1.0 - p.sum(axis=1), 0.0, 1.0)]) / ns
    pj = joint.sum(axis=0)
    ps = joint.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = joint * np.log2(joint / (ps[:, None] * pj[None, :]))
    return float(np.nansum(np.where(joint > 0, terms, 0.0)))


def evaluate(
    lib: TrajectoryLibrary, rates: FloatArray, arrangement: Sequence[int], k_ns: float = 0.0
) -> float:
    bm = BindingMatrix(rates, tuple(str(i) for i in range(rates.shape[0])), tuple(arrangement))
    return information_per_cell(bm.profiles(lib, k_ns))


@dataclass(frozen=True)
class ArrangementResult:
    arrangement: tuple[int, ...]
    information_bits: float
    evaluations: int


def optimize_arrangement(
    lib: TrajectoryLibrary,
    rates: ArrayLike,
    initial: Sequence[int],
    seed: int,
    restarts: int = 4,
    max_passes: int = 20,
    k_ns: float = 0.0,
) -> ArrangementResult:
    """Pair-swap hill climbing from the initial arrangement and random shuffles of it."""
    r = np.asarray(rates, dtype=np.float64)
    g = np.random.default_rng(seed)
    best = tuple(initial)
    best_val = evaluate(lib, r, best, k_ns)
    n_eval = 1
    s = len(best)
    pairs = np.array([(a, b) for a in range(s) for b in range(a + 1, s)])
    for k in range(restarts + 1):
        cur = list(best) if k == 0 else list(g.permutation(np.asarray(initial)))
        val = evaluate(lib, r, cur, k_ns)
        n_eval += 1
        for _ in range(max_passes):
            improved = False
            for i, j in g.permutation(pairs):
                if cur[i] == cur[j]:
                    continue
                cur[i], cur[j] = cur[j], cur[i]
                v = evaluate(lib, r, cur, k_ns)
                n_eval += 1
                if v > val + 1e-12:
                    val, improved = v, True
                else:
                    cur[i], cur[j] = cur[j], cur[i]
            if not improved:
                break
        if val > best_val + 1e-12:
            best, best_val = tuple(int(a) for a in cur), val
    return ArrangementResult(tuple(int(a) for a in best), best_val, n_eval)


def brute_force_arrangement(
    lib: TrajectoryLibrary, rates: ArrayLike, aptamers: Sequence[int], k_ns: float = 0.0
) -> ArrangementResult:
    r = np.asarray(rates, dtype=np.float64)
    best: tuple[int, ...] = tuple(aptamers)
    best_val = -1.0
    n = 0
    for perm in set(permutations(aptamers)):
        v = evaluate(lib, r, perm, k_ns)
        n += 1
        if v > best_val + 1e-12:
            best, best_val = tuple(perm), v
    return ArrangementResult(best, best_val, n)
