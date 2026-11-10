"""Stage Z11: CUDA backend against the CPU core and the reference (T40).

Runs only where the core was built with BINDOPTIX_ENABLE_CUDA=ON and a GPU is present; with
BINDOPTIX_REQUIRE_CUDA=1 a missing CUDA backend fails instead of skipping (large machine).
Device math (exp, log, sin, cos) rounds differently from libm; with --fmad=false the paths agree
to about 1e-9 µm over the test runs, discrete outcomes identically.
"""

from __future__ import annotations

import math
import os
import time

import numpy as np
import pytest

from bindoptix import fast
from bindoptix.cells import CellModel
from bindoptix.flow import DuctFlow
from bindoptix.particles import ChannelModel, Spot, inlet_positions, run

if not fast.cuda_available():
    if os.environ.get("BINDOPTIX_REQUIRE_CUDA") == "1":
        raise RuntimeError("BINDOPTIX_REQUIRE_CUDA=1 but the CUDA backend is not available")
    pytest.skip("CUDA backend not built or no GPU", allow_module_level=True)

SEED = 20261119
SPOTS = tuple(
    Spot(np.array([[x - 100, -100], [x + 100, -100], [x + 100, 100], [x - 100, 100]], float), k)
    for x, k in ((800.0, 2e-3), (1600.0, math.inf), (2400.0, 0.0))
)
MODELS = {
    "slide": ChannelModel(
        DuctFlow(200.0, 50.0, 5e5 / 6), CellModel(), spots=SPOTS, x_out_um=3000.0
    ),
    "wall_swim": ChannelModel(
        DuctFlow(200.0, 100.0, 5e4),
        CellModel(swim_speed_um_s=15.0),
        spots=SPOTS,
        wall_hydrodynamics=True,
        contact_gap_um=0.02,
        nonspecific_rate_s=1e-3,
    ),
}


@pytest.mark.parametrize("name", list(MODELS))
def test_t40_cuda_equals_cpu_and_reference(name: str) -> None:
    m = MODELS[name]
    start = inlet_positions(m.flow, 2000, SEED)
    cpu = fast.run(m, start, 1800.0, 5.0, SEED, record_visits=True, backend="cpu")
    gpu = fast.run(m, start, 1800.0, 5.0, SEED, record_visits=True, backend="cuda")
    np.testing.assert_allclose(gpu.position_um, cpu.position_um, rtol=0, atol=1e-8)
    np.testing.assert_array_equal(gpu.state, cpu.state)
    np.testing.assert_array_equal(gpu.spot, cpu.spot)
    np.testing.assert_allclose(gpu.visit_time_s, cpu.visit_time_s, atol=1e-9)
    ref = run(m, start[:300], 1800.0, 5.0, SEED)
    np.testing.assert_allclose(gpu.position_um[:300], ref.position_um, rtol=0, atol=1e-8)
    np.testing.assert_array_equal(gpu.state[:300], ref.state)


def test_t40_cuda_throughput() -> None:
    # Smoke test of the scale the GPU is meant for: 10⁶ cells, 100 steps.
    m = MODELS["slide"]
    start = inlet_positions(m.flow, 1_000_000, SEED)
    t = time.perf_counter()
    r = fast.run(m, start, 500.0, 5.0, SEED, backend="cuda")
    assert r.position_um.shape == (1_000_000, 3)
    assert time.perf_counter() - t < 120.0
