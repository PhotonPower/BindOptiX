"""Counter-based random numbers: value = f(seed, particle id, step, lane).

The numbers do not depend on how particles are split into batches or in which order they are
processed (test case T11). SplitMix64 finaliser chained over the four keys; uniform numbers from
the upper 53 bits, normal numbers by Box–Muller.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

_GOLDEN = np.uint64(0x9E3779B97F4A7C15)
_M1 = np.uint64(0xBF58476D1CE4E5B9)
_M2 = np.uint64(0x94D049BB133111EB)


def _mix(x: NDArray[np.uint64]) -> NDArray[np.uint64]:
    z = x + _GOLDEN
    z = (z ^ (z >> np.uint64(30))) * _M1
    z = (z ^ (z >> np.uint64(27))) * _M2
    return z ^ (z >> np.uint64(31))


def _hash(seed: int, pid: ArrayLike, step: int, lane: int) -> NDArray[np.uint64]:
    with np.errstate(over="ignore"):
        p = np.asarray(pid, dtype=np.uint64)
        h = _mix(np.full(p.shape, np.uint64(seed & 0xFFFFFFFFFFFFFFFF), dtype=np.uint64))
        h = _mix(h ^ p)
        h = _mix(h ^ np.uint64(step))
        return _mix(h ^ np.uint64(lane))


def uniform(seed: int, pid: ArrayLike, step: int, lane: int) -> NDArray[np.float64]:
    """Uniform numbers in [0, 1)."""
    h = _hash(seed, pid, step, lane)
    return np.asarray((h >> np.uint64(11)).astype(np.float64) * (1.0 / 9007199254740992.0))


def normal(seed: int, pid: ArrayLike, step: int, lane: int) -> NDArray[np.float64]:
    """Standard normal numbers (Box–Muller over two hashed uniforms)."""
    u1 = 1.0 - uniform(seed, pid, step, 2 * lane)
    u2 = uniform(seed, pid, step, 2 * lane + 1)
    return np.asarray(np.sqrt(-2.0 * np.log(u1)) * np.cos(2.0 * np.pi * u2))
