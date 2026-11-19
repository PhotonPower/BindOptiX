"""T43: port Q from BioOptoMD and the mixture B · Q (DP7 of the ADR draft).

The ensemble of a spot is Σ_s f_s Q_s over the receptor states s, weighted per sample (ADR draft,
section 2). Tracker output of BioOptoMD (PHYSICS_MODEL.md 6.1): (n_replicas, n_frames, 7) =
(R_F [nm], μ̂, d [nm]) with d = z(R_F) − z_surface.
"""

import numpy as np
import pytest

from bindoptix.binding import Occupancy
from bindoptix.ports import EmitterResponseTable, RegularGrid
from bindoptix.signal import (
    EmitterEnsemble,
    Fluorophore,
    detected_rate_s,
    mix_ensembles,
    spot_rates_s,
    uniform_field_map,
)

DYE = Fluorophore(250_000.0, 0.27, 0.640, 0.670, bleach_yield=1e-6)
FIELD = uniform_field_map(50.0, DYE.excitation_um)


def height_table() -> EmitterResponseTable:
    """Quenching and collection that change with the height (values arbitrary, not optics)."""
    z = np.array([0.001, 0.005, 0.02])
    eye = np.eye(3)
    a_tot = np.stack([k * eye for k in (40.0, 4.0, 1.2)])[None, None]
    a_rad = np.stack([k * eye for k in (0.5, 0.8, 1.0)])[None, None]
    a_det = np.stack([[0.3 * k * eye, 0.6 * k * eye] for k in (0.5, 0.8, 1.0)])[None, None]
    grid = RegularGrid((np.zeros(1), np.zeros(1), z))
    return EmitterResponseTable(grid, 0.670, 1.333, ("top", "bottom"), a_tot, a_rad, a_det)


TABLE = height_table()


def tracker_series() -> np.ndarray:  # type: ignore[type-arg]
    rng = np.random.default_rng(20261119)
    s = np.empty((2, 5, 7))
    s[..., 0:2] = rng.uniform(10.0, 12.0, (2, 5, 2))  # R_F x, y [nm]
    s[..., 2] = rng.uniform(30.0, 40.0, (2, 5))  # R_F z in the box, not the height
    mu = rng.normal(size=(2, 5, 3))
    s[..., 3:6] = mu / np.linalg.norm(mu, axis=-1, keepdims=True)
    s[..., 6] = rng.uniform(2.0, 15.0, (2, 5))  # d [nm]
    return s


def test_t43_tracker_output_becomes_an_ensemble() -> None:
    s = tracker_series()
    q = EmitterEnsemble.from_tracker(s, anchor_xy_nm=(11.0, 11.0))
    flat = s.reshape(-1, 7)
    np.testing.assert_allclose(q.offsets_um[:, :2], (flat[:, :2] - 11.0) * 1e-3, rtol=1e-15)
    np.testing.assert_allclose(q.offsets_um[:, 2], flat[:, 6] * 1e-3, rtol=1e-15)
    np.testing.assert_allclose(q.mu_abs, flat[:, 3:6], rtol=1e-12)
    np.testing.assert_allclose(np.asarray(q.weights), 0.1)
    centred = EmitterEnsemble.from_tracker(s)
    np.testing.assert_allclose(centred.offsets_um[:, :2].mean(axis=0), 0.0, atol=1e-15)
    with pytest.raises(ValueError):
        EmitterEnsemble.from_tracker(s[..., :6])


def ensembles() -> dict[str, EmitterEnsemble]:
    s = tracker_series()
    return {
        "free": EmitterEnsemble.from_tracker(s[0], brightness=0.01),
        "target": EmitterEnsemble.from_tracker(s[1]),
        "albumin": EmitterEnsemble.isotropic(0.003, 7, 1, brightness=0.02),
    }


def test_t43_dp7_without_analyte_only_the_free_state() -> None:
    q = ensembles()
    mix = mix_ensembles(list(q.values()), [1.0, 0.0, 0.0])
    free = q["free"]
    np.testing.assert_array_equal(mix.offsets_um, free.offsets_um)
    np.testing.assert_array_equal(mix.mu_abs, free.mu_abs)
    np.testing.assert_allclose(np.asarray(mix.weights), np.asarray(free.weights), rtol=1e-15)
    assert mix.brightness == pytest.approx(free.brightness, rel=1e-15)


@pytest.mark.parametrize("f", [(0.7, 0.2, 0.1), (0.0, 1.0, 0.0), (0.2, 0.5, 0.3)])
def test_t43_dp7_mixture_weights_samples_linearly(f: tuple[float, float, float]) -> None:
    # Σ_s f_s Q_s weighted per sample: the signal is Σ_s f_s · signal(Q_s); weights sum to 1.
    q = list(ensembles().values())
    mix = mix_ensembles(q, f)
    assert mix.weights is not None and mix.weights.sum() == pytest.approx(1.0, abs=1e-15)
    got = detected_rate_s([0, 0, 0], mix, FIELD, TABLE, DYE, 3.0)
    want = sum(
        fs * detected_rate_s([0, 0, 0], qs, FIELD, TABLE, DYE, 3.0)
        for fs, qs in zip(f, q, strict=True)
    )
    np.testing.assert_allclose(got, want, rtol=1e-12)
    with pytest.raises(ValueError):
        mix_ensembles(q, [0.5, 0.2, 0.2])


def test_t43_spot_rates_follow_port_b() -> None:
    q = ensembles()
    fr = np.array([[[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]], [[0.6, 0.3, 0.1], [0.2, 0.7, 0.1]]])
    occ = Occupancy(np.array([0.0, 600.0]), (3, 4), ("free", "target", "albumin"), fr)
    anchors = np.array([[0.0, 0.0, 0.0], [800.0, 0.0, 0.0]])
    rates = spot_rates_s(occ, q, anchors, FIELD, TABLE, DYE, 3.0)
    assert rates.shape == (2, 2, 2)
    per_state = {s: detected_rate_s([0, 0, 0], e, FIELD, TABLE, DYE, 3.0) for s, e in q.items()}
    for t in range(2):
        for k in range(2):
            want = sum(fr[t, k, j] * per_state[s] for j, s in enumerate(occ.states))
            np.testing.assert_allclose(rates[t, k], want, rtol=1e-12)
    with pytest.raises(ValueError, match="albumin"):
        spot_rates_s(
            occ, {"free": q["free"], "target": q["target"]}, anchors, FIELD, TABLE, DYE, 3.0
        )
