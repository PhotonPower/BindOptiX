"""Stage Z10, part 1: trajectory library and classification (T25, T38)."""

import math

import numpy as np
import pytest

from bindoptix.cells import CellModel
from bindoptix.chip import ChipLayout
from bindoptix.chipsim import channel_model
from bindoptix.classify import (
    BindingMatrix,
    TrajectoryLibrary,
    classify,
    confusion_given_total,
    confusion_matrix,
)
from bindoptix.frames import Mount
from bindoptix.particles import CAPTURED, STUCK, inlet_positions, run

SEED = 20261117
THETA = math.atan(1 / 5)


def library_and_direct(
    rates_grid: np.ndarray,  # type: ignore[type-arg]
    k_ns: float,
    n: int = 1500,
) -> tuple[TrajectoryLibrary, object, object]:
    layout = ChipLayout(grid_angle_rad=THETA)
    q = 5e6 / 60.0
    lib_model, cs = channel_model(layout, 2, CellModel(), q, Mount.sensor_floor(), 0.0)
    start = inlet_positions(lib_model.flow, n, SEED)
    lib_run = run(lib_model, start, 3 * 3600.0, 10.0, SEED, record_visits=True)
    lib = TrajectoryLibrary.from_runs([lib_run], [cs], 25)
    local = [float(rates_grid[g]) for g in cs.grid_index]
    model, _ = channel_model(
        layout, 2, CellModel(), q, Mount.sensor_floor(), 0.0, nonspecific_rate_s=k_ns
    )
    from dataclasses import replace

    from bindoptix.particles import Spot

    spots = tuple(Spot(s.polygon_um, r) for s, r in zip(model.spots, local, strict=True))
    direct = run(replace(model, spots=spots), start, 3 * 3600.0, 10.0, SEED)
    return lib, direct, cs


@pytest.mark.parametrize(("rates", "k_ns"), [({11: 2e-3, 13: 1e-2}, 0.0), ({10: 5e-3}, 2e-4)])
def test_t38_library_predicts_direct_capture(rates: dict[int, float], k_ns: float) -> None:
    # Same random streams: the library paths equal the direct paths until capture; the surrogate
    # probabilities must explain the direct captures per spot within the Poisson-binomial spread.
    r = np.zeros(25)
    for g, k in rates.items():
        r[g] = k
    lib, direct, cs = library_and_direct(r, k_ns)
    p = lib.capture_probabilities(r, k_ns)
    got = direct.state == CAPTURED  # type: ignore[attr-defined]
    for local, g in enumerate(cs.grid_index):  # type: ignore[attr-defined]
        expected = p[:, g].sum()
        var = (p[:, g] * (1 - p[:, g])).sum()
        observed = int(np.sum(got & (direct.spot == local)))  # type: ignore[attr-defined]
        assert abs(observed - expected) <= 4 * math.sqrt(var) + 1.0
    if k_ns > 0:
        stuck = direct.count(STUCK)  # type: ignore[attr-defined]
        assert stuck > 0


def test_t25_exact_confusion_one_cell() -> None:
    # With one captured cell the decision is argmax_s q_s,j; independent closed form.
    prof = np.array([[0.3, 0.1, 0.05], [0.1, 0.2, 0.05], [0.05, 0.05, 0.2]])
    q = prof / prof.sum(axis=1, keepdims=True)
    dec = np.argmax(q, axis=0)
    ref = np.zeros((3, 3))
    for a in range(3):
        for j in range(3):
            ref[a, dec[j]] += q[a, j]
    np.testing.assert_allclose(confusion_given_total(prof, 1), ref, atol=1e-14)


def test_t25_monte_carlo_confusion_matches_exact() -> None:
    prof = np.array([[0.02, 0.01, 0.005], [0.01, 0.02, 0.005], [0.005, 0.005, 0.02]])
    exact = confusion_given_total(prof, 3)
    # condition the Monte Carlo on exactly three captures by rejection
    g = np.random.default_rng(SEED)
    counts = np.zeros((3, 3))
    for a in range(3):
        probs = np.append(prof[a] / prof[a].sum(), 0.0)
        draws = g.multinomial(3, probs[:3], size=20000)
        for d in draws:
            counts[a, classify(d, prof)] += 1
    np.testing.assert_allclose(counts / 20000, exact, atol=0.015)
    mc = confusion_matrix(prof, 200, 2000, SEED)
    assert mc.shape == (3, 4)
    np.testing.assert_allclose(mc.sum(axis=1), 1.0)


def test_classify_ignores_concentration() -> None:
    prof = np.array([[0.2, 0.1], [0.05, 0.1]])
    assert classify([20, 10], prof) == 0
    assert classify([2, 1], prof) == 0  # same shape, ten times fewer cells
    assert classify([0, 0], prof) == -1


def test_binding_matrix_spot_rates() -> None:
    bm = BindingMatrix(np.array([[1.0, 0.1], [0.2, 2.0]]), ("A", "B"), (0, 1, 1, 0))
    np.testing.assert_allclose(bm.spot_rates(1), [0.2, 2.0, 2.0, 0.2])
    with pytest.raises(ValueError):
        BindingMatrix(np.array([[1.0]]), ("A",), (0, 1))
