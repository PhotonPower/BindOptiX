"""Stage Z11: compiled core against the Python reference (T39, docs/testfaelle.md).

Runs where bindoptix._core is built (CI: yes); with BINDOPTIX_REQUIRE_CORE=1 a missing core
fails instead of skipping. Agreement: same random numbers (bit-exact integers), positions to
1e-9 µm, states, spots and event times identical.
"""

from __future__ import annotations

import math
from dataclasses import replace
from typing import Any

import numpy as np
import pytest

from bindoptix import fast, rng
from bindoptix.cells import CellModel
from bindoptix.flow import DuctFlow
from bindoptix.particles import ChannelModel, RunResult, Spot, inlet_positions, run

if not fast.CORE_AVAILABLE:
    pytest.skip("bindoptix._core not built", allow_module_level=True)

SEED = 20261118
CORE: Any = getattr(fast, "_core")  # noqa: B009
SPOTS = tuple(
    Spot(np.array([[x - 100, -100], [x + 100, -100], [x + 100, 100], [x - 100, 100]], float), k)
    for x, k in ((800.0, 2e-3), (1600.0, math.inf), (2400.0, 0.0))
)


def same(a: RunResult, b: RunResult, atol: float = 1e-9) -> None:
    np.testing.assert_allclose(a.position_um, b.position_um, rtol=0, atol=atol)
    np.testing.assert_array_equal(a.state, b.state)
    np.testing.assert_array_equal(a.spot, b.spot)
    np.testing.assert_array_equal(np.isnan(a.event_time_s), np.isnan(b.event_time_s))
    np.testing.assert_allclose(a.event_time_s, b.event_time_s, rtol=1e-12, equal_nan=True)
    np.testing.assert_allclose(a.orientation, b.orientation, atol=1e-12)
    assert a.steps == b.steps and a.t_s == pytest.approx(b.t_s)


def test_t39_random_numbers_bit_exact() -> None:
    ids = np.arange(5000, dtype=np.uint64) * 7919 + 3
    for lane in (0, 10, 22, 101):
        np.testing.assert_array_equal(
            CORE.uniform(SEED, ids, 12345, lane), rng.uniform(SEED, ids, 12345, lane)
        )
    np.testing.assert_allclose(
        CORE.normal(SEED, ids, 7, 2), rng.normal(SEED, ids, 7, 2), rtol=1e-14, atol=1e-15
    )


CASES = {
    "flow_sediment": ChannelModel(
        DuctFlow(200.0, 100.0, 5e5 / 6), CellModel(), spots=SPOTS, x_out_um=3000.0
    ),
    "slide_nonspecific": ChannelModel(
        DuctFlow(200.0, 50.0, 5e5 / 6),
        CellModel(),
        spots=SPOTS,
        nonspecific_rate_s=1e-3,
        x_out_um=3000.0,
    ),
    "stick": ChannelModel(
        DuctFlow(200.0, 100.0, 2.5e4), CellModel(), spots=SPOTS, nonspecific_rate_s=math.inf
    ),
    "wall_hydrodynamics": ChannelModel(
        DuctFlow(200.0, 100.0, 5e4),
        CellModel(),
        spots=SPOTS,
        wall_hydrodynamics=True,
        contact_gap_um=0.02,
    ),
    "swimming": ChannelModel(
        DuctFlow(200.0, 100.0, 5e4), CellModel(swim_speed_um_s=15.0), spots=SPOTS
    ),
    "reflect_tilted": ChannelModel(
        DuctFlow(200.0, 100.0, 0.0), CellModel(), floor="reflect", gravity_dir=(0.3, -0.2, -1.0)
    ),
}


@pytest.mark.parametrize("name", list(CASES))
def test_t39_core_equals_reference(name: str) -> None:
    m = CASES[name]
    start = inlet_positions(DuctFlow(200.0, m.flow.height_um, 1e5), 300, SEED)
    a = run(m, start, 1800.0, 5.0, SEED, snapshot_times_s=(600.0, 1200.0))
    b = fast.run(m, start, 1800.0, 5.0, SEED, snapshot_times_s=(600.0, 1200.0), backend="cpu")
    same(a, b)
    for t in (600.0, 1200.0):
        np.testing.assert_allclose(a.snapshots_um[t], b.snapshots_um[t], atol=1e-9)


def test_t39_waiting_resume_visits_and_rate_table() -> None:
    m = replace(
        CASES["slide_nonspecific"],
        nonspecific_rate_s=0.0,
        spot_rates_fn=lambda t: np.array([1e-3 * (t < 900.0), 0.0, 5e-3]),
    )
    start = inlet_positions(m.flow, 200, SEED)
    t_in = np.linspace(0.0, 600.0, 200)
    a1 = run(m, start, 700.0, 5.0, SEED, start_times_s=t_in, record_visits=True)
    b1 = fast.run(m, start, 700.0, 5.0, SEED, start_times_s=t_in, record_visits=True, backend="cpu")
    same(a1, b1)
    a2 = run(m, None, 1100.0, 5.0, SEED, resume=a1, record_visits=True)
    b2 = fast.run(m, None, 1100.0, 5.0, SEED, resume=b1, record_visits=True, backend="cpu")
    same(a2, b2)
    for name in ("visit_entry_s", "visit_time_s", "visit_off_floor_s", "off_floor_s"):
        np.testing.assert_allclose(getattr(a2, name), getattr(b2, name), atol=1e-9, equal_nan=True)


def test_backend_selection() -> None:
    assert fast._resolve("auto") in ("cpu", "cuda")
    assert fast._resolve("python") == "python"
    with pytest.raises(ValueError):
        fast._resolve("gpu")
