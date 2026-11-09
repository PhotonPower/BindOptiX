"""T34: hydrodynamic hindrance near the floor (docs/testfaelle.md)."""

import math

import numpy as np
import pytest
from scipy.integrate import quad

from bindoptix import wall
from bindoptix.cells import CellModel
from bindoptix.flow import DuctFlow
from bindoptix.particles import ChannelModel, run

SEED = 20261113


def brenner_perp(gap: float, a: float, terms: int = 400) -> float:
    """Independent reference: Brenner's exact series for λ⊥ (brenner1961)."""
    alpha = math.acosh((a + gap) / a)
    s = 0.0
    for n in range(1, terms + 1):
        if (2 * n + 1) * alpha > 600.0:  # remaining terms below 1e-200
            break
        num = 2.0 * math.sinh((2 * n + 1) * alpha) + (2 * n + 1) * math.sinh(2 * alpha)
        den = 4.0 * math.sinh((n + 0.5) * alpha) ** 2 - (2 * n + 1) ** 2 * math.sinh(alpha) ** 2
        s += n * (n + 1) / ((2 * n - 1) * (2 * n + 3)) * (num / den - 1.0)
    return 1.0 / (4.0 / 3.0 * math.sinh(alpha) * s)


@pytest.mark.parametrize("gap_over_a", [0.05, 0.1, 0.3, 1.0, 3.0, 10.0])
def test_t34_perp_against_brenner(gap_over_a: float) -> None:
    a = 0.45
    approx = float(wall.lambda_perp(gap_over_a * a, a))
    assert approx == pytest.approx(brenner_perp(gap_over_a * a, a), rel=0.02)


def test_t34_perp_limits_and_derivative() -> None:
    a = 0.5
    assert float(wall.lambda_perp(1e-6, a)) == pytest.approx(1e-6 / a, rel=1e-4)
    assert float(wall.lambda_perp(1e6, a)) == pytest.approx(1.0, abs=1e-5)
    h = np.array([0.01, 0.2, 1.0, 5.0])
    eps = 1e-7
    fd = (wall.lambda_perp(h + eps, a) - wall.lambda_perp(h - eps, a)) / (2 * eps)
    np.testing.assert_allclose(wall.dlambda_perp(h, a), fd, rtol=1e-6)


def test_t34_parallel_faxen() -> None:
    # Hand evaluation of Faxén's series at r = a/(a + h) = 1/2.
    ref = 1 - 9 / 32 + 1 / 64 - 45 / 4096 - 1 / 512
    assert float(wall.lambda_par(1.0, 1.0)) == pytest.approx(ref, rel=1e-14)
    assert float(wall.lambda_par(1e6, 1.0)) == pytest.approx(1.0, abs=1e-5)


def test_t34_boltzmann_equilibrium_with_hindrance() -> None:
    # With position-dependent D and the Itô drift the equilibrium is still c ∝ exp(−z v/D₀);
    # ℓ = D₀/v ≈ 1.6 µm, so the hindered layer (a few a) dominates the profile.
    cell = CellModel(width_um=1.0, length_um=1.0, delta_rho_kg_m3=500.0)
    d, v = cell.diffusion_um2_s, cell.sedimentation_speed_um_s()
    ell, h = d / v, 30.0
    m = ChannelModel(DuctFlow(1e9, h, 0.0), cell, floor="reflect", wall_hydrodynamics=True)
    n = 2500
    g = np.random.default_rng(SEED)
    start = np.column_stack([np.zeros(n), np.zeros(n), g.exponential(ell, n)])
    r = run(m, start, 150.0, 0.03, SEED)
    z = r.position_um[:, 2]
    mean_ref = ell - h / math.expm1(h / ell)
    var_ref = ell**2 - h**2 * math.exp(h / ell) / math.expm1(h / ell) ** 2
    assert abs(z.mean() - mean_ref) < 4 * math.sqrt(var_ref / n) + 0.02 * ell


def test_t34_hindered_approach_time() -> None:
    # Deterministic sedimentation: t = ∫ dh / (v λ⊥(h)) from the start gap to the contact gap.
    cell = CellModel()
    a, v = cell.equivalent_radius_um, cell.sedimentation_speed_um_s()
    zc, z0 = 0.02, 5.0
    m = ChannelModel(
        DuctFlow(1e9, 1e9, 0.0),
        cell,
        nonspecific_rate_s=math.inf,
        wall_hydrodynamics=True,
        contact_gap_um=zc,
        brownian=False,
    )
    t_ref = quad(lambda hh: 1.0 / (v * float(wall.lambda_perp(hh, a))), zc, z0, limit=200)[0]
    r = run(m, np.array([[0.0, 0.0, z0]]), 2.0 * t_ref, 0.5, SEED)
    assert r.event_time_s[0] == pytest.approx(t_ref, rel=0.01)
    assert t_ref > 1.5 * (z0 - zc) / v  # hindrance slows the last micrometres noticeably


def test_t34_requires_contact_gap() -> None:
    with pytest.raises(ValueError, match="contact_gap"):
        ChannelModel(DuctFlow(1.0, 1.0, 0.0), CellModel(), wall_hydrodynamics=True)
