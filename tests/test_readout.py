"""Stage Z7: protocols, synthetic frames, readout (T21, T22, T37)."""

import math

import numpy as np
import pytest
from numpy.typing import NDArray

from bindoptix.cells import CellModel
from bindoptix.continuum import ChannelSection, ReactiveRegion, Species
from bindoptix.flow import DuctFlow
from bindoptix.imaging import FrameGeometry, fixed_pattern, render
from bindoptix.particles import FREE, WAITING, ChannelModel, inlet_positions, run
from bindoptix.protocol import Phase, injection_schedule, run_continuum, run_particles
from bindoptix.readout import Tracker, counts_per_spot, detect
from bindoptix.signal import Camera

SEED = 20261116
CAM = Camera(2.44, 0.8, 1.5, 0.2, 1.0, 1.2)
GEOM = FrameGeometry(0.0, 0.0, 240, 240, 2.44)
SPOT = np.array([[50.0, 50.0], [535.0, 50.0], [535.0, 535.0], [50.0, 535.0]])


def scene(seed: int) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Cells on the spot (≥ 25 µm apart), fixed pattern, transient debris positions."""
    g = np.random.default_rng(seed)
    cells: list[NDArray[np.float64]] = []
    while len(cells) < 25:
        p = g.uniform(80.0, 505.0, 2)
        if all(np.hypot(*(p - q)) > 25.0 for q in cells):
            cells.append(p)
    pattern = fixed_pattern(GEOM, 0.3, 10.0, seed)
    debris = g.uniform(80.0, 505.0, (10, 2))
    return np.array(cells), pattern, debris


def frames(seed: int) -> tuple[list[NDArray[np.float64]], NDArray[np.float64], NDArray[np.float64]]:
    cells, pattern, debris = scene(seed)
    g = np.random.default_rng(seed + 1)
    areas = ((SPOT, 50.0),)
    base = render(GEOM, 1.2, areas=areas, area_pattern=pattern, uniform_photons_per_um2=2.0)
    out = [CAM.sample(base, g)]  # reference
    for k in range(1, 4):
        pts, ph = cells, np.full(len(cells), 800.0)
        if k == 2:
            pts, ph = np.vstack([cells, debris]), np.concatenate([ph, np.full(10, 1500.0)])
        img = base + render(GEOM, 1.2, pts, ph)
        out.append(CAM.sample(img, g))
    expected_bg = CAM.quantum_efficiency * base + CAM.dark_e_s * CAM.exposure_s
    return out, expected_bg, cells


def run_tracker(
    persistence: int, use_reference: bool, seed: int
) -> tuple[Tracker, NDArray[np.float64]]:
    imgs, bg, cells = frames(seed)
    ref = imgs[0]
    if not use_reference:
        # model background without the fixed pattern
        flat = render(GEOM, 1.2, areas=((SPOT, 50.0),), uniform_photons_per_um2=2.0)
        ref = CAM.quantum_efficiency * flat + CAM.dark_e_s * CAM.exposure_s
    tr = Tracker(match_um=4.0, persistence=persistence)
    for k in range(1, 4):
        tr.update(k, detect(imgs[k], ref, bg, GEOM, CAM, 0.05))
    return tr, cells


def matched(tr: Tracker, cells: NDArray[np.float64]) -> tuple[int, int]:
    hit, false = 0, 0
    for t in tr.confirmed:
        d = np.hypot(*(cells - t.xy_um).T).min()
        if d <= 3.0:
            hit += 1
        else:
            false += 1
    return hit, false


def test_t21_difference_and_persistence_find_cells() -> None:
    tr, cells = run_tracker(persistence=2, use_reference=True, seed=SEED)
    hit, false = matched(tr, cells)
    assert hit >= 24
    assert false == 0


def test_t21_without_reference_the_pattern_gives_false_alarms() -> None:
    tr, cells = run_tracker(persistence=2, use_reference=False, seed=SEED)
    _, false = matched(tr, cells)
    assert false > 10


def test_t21_without_persistence_debris_counts() -> None:
    tr, cells = run_tracker(persistence=1, use_reference=True, seed=SEED)
    _, false = matched(tr, cells)
    assert false >= 8


def test_t21_counts_per_spot() -> None:
    tr, cells = run_tracker(persistence=2, use_reference=True, seed=SEED)
    left = np.array([[50.0, 50.0], [292.5, 50.0], [292.5, 535.0], [50.0, 535.0]])
    right = left + np.array([242.5, 0.0])
    c = counts_per_spot(tr.confirmed, [left, right], 4)
    n_left = int(np.sum(cells[:, 0] < 292.5))
    assert c[0].sum() == 0
    assert abs(int(c[3, 0]) - n_left) <= 1
    assert int(c[3].sum()) == len(tr.confirmed)


def test_render_conserves_photons() -> None:
    img = render(GEOM, 1.2, [[300.0, 300.0]], [1000.0])
    assert img.sum() == pytest.approx(1000.0, rel=1e-8)  # window ±6σ
    img2 = render(GEOM, 1.2, areas=((SPOT, 2.0),))
    assert img2.sum() == pytest.approx(2.0 * 485.0**2, rel=0.02)


# --- T22: washing ------------------------------------------------------------------------------


def test_t22_wash_releases_reversible_binding() -> None:
    # Fast transport, slow binding, low coverage: no rebinding; after the wash Γ ∝ exp(−k_off t).
    k_off = 5e-3

    def make(q: float) -> ChannelSection:
        reg = ReactiveRegion(50.0, 150.0, 1e-2, (5.0,), (k_off,))
        return ChannelSection(
            np.linspace(0.0, 200.0, 41),
            np.linspace(0.0, 20.0, 11),
            lambda z: q * 6.0 * z / 20.0 * (1 - z / 20.0) / 20.0,
            [Species("target", 500.0)],
            [reg],
        )

    phases = [Phase("load", 300.0, 4000.0, (1e-3,)), Phase("wash", 300.0, 8000.0, (0.0,))]
    ends, times, cov = run_continuum(make, phases, 1.0)
    g = np.array([c[0, 0] for c in cov])
    i_wash = int(np.argmin(np.abs(times - 300.0)))
    rel = g[i_wash:] / g[i_wash]
    ref = np.exp(-k_off * (times[i_wash:] - 300.0))
    # implicit Euler and the wash-out of the channel (~1 s) bias the decay by < 2 %
    np.testing.assert_allclose(rel[30:], ref[30:], rtol=0.02)


# --- T37: protocol mechanics -------------------------------------------------------------------


def channel(q: float) -> ChannelModel:
    return ChannelModel(DuctFlow(200.0, 100.0, q), CellModel(), x_out_um=5000.0)


def test_t37_resume_equals_single_run() -> None:
    m = channel(2.5e4)
    start = inlet_positions(m.flow, 100, SEED)
    one = run(m, start, 1200.0, 5.0, SEED)
    a = run(m, start, 600.0, 5.0, SEED)
    b = run(m, None, 600.0, 5.0, SEED, resume=a)
    np.testing.assert_array_equal(b.position_um, one.position_um)
    np.testing.assert_array_equal(b.state, one.state)
    assert b.t_s == pytest.approx(1200.0)


def test_t37_injection_and_waiting() -> None:
    phases = [Phase("load", 100.0, 2.5e4, cells_per_s=0.2), Phase("incubate", 50.0, 0.0)]
    t_in = injection_schedule(phases)
    assert t_in.size == 20
    np.testing.assert_allclose(t_in[:2], [2.5, 7.5])
    m = channel(2.5e4)
    start = inlet_positions(m.flow, 3, SEED)
    r = run(m, start, 10.0, 1.0, SEED, start_times_s=[0.0, 5.0, 50.0])
    assert r.state[2] == WAITING
    np.testing.assert_array_equal(r.position_um[2], start[2])
    assert r.state[0] != WAITING and r.state[1] != WAITING


def test_t37_incubation_without_flow() -> None:
    # During a phase with Q = 0 cells only sediment and diffuse: no net drift along x.
    phases = [Phase("load", 200.0, 2.5e4, cells_per_s=1.0), Phase("incubate", 600.0, 0.0)]
    res = run_particles(channel, phases, 2.0, SEED)
    a, b = res
    moving = (a.state == FREE) & (b.state == FREE)
    dx = b.position_um[moving, 0] - a.position_um[moving, 0]
    d = CellModel().diffusion_um2_s
    assert abs(dx.mean()) < 4 * math.sqrt(2 * d * 600.0 / max(moving.sum(), 1)) + 0.5
    assert len(res) == 2


def test_t21_counts_with_margin_at_spot_edge() -> None:
    from bindoptix.readout import Track

    poly = np.array([[0.0, 0.0], [200.0, 0.0], [200.0, 200.0], [0.0, 200.0]])
    edge_cell = Track(np.array([-1.5, 50.0]), 1, 2, 2)  # fitted just upstream of the edge
    assert counts_per_spot([edge_cell], [poly], 3)[2, 0] == 0
    assert counts_per_spot([edge_cell], [poly], 3, margin_um=3.0)[2, 0] == 1
