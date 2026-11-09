"""Stage Z4: convection–diffusion with surface reaction (T12–T14, docs/testfaelle.md)."""

import math

import numpy as np
import pytest
from numpy.typing import NDArray

from bindoptix.continuum import ChannelSection, ReactiveRegion, Species, graded_edges
from bindoptix.flow import DuctFlow

H, U_MEAN = 100.0, 1000.0
SHEAR = 6.0 * U_MEAN / H


def plane_poiseuille(z: NDArray[np.float64]) -> NDArray[np.float64]:
    return np.asarray(6.0 * U_MEAN * z / H * (1.0 - z / H), dtype=np.float64)


def x_grid(x0: float, length: float, tail: float) -> NDArray[np.float64]:
    """Edges 0 … x0 + L + tail, refined towards the leading edge x0 from both sides."""
    up = x0 - graded_edges(x0, 30, 0.5)[::-1]
    on = x0 + graded_edges(length, 120, 0.2)
    down = x0 + length + graded_edges(tail, 20, 2.0)
    return np.asarray(np.unique(np.concatenate([up, on, down])), dtype=np.float64)


def leveque_flux(c0: float, d: float, length: float) -> float:
    """Total flux per width into a perfectly absorbing strip (Lévêque; squires2008: 0.81)."""
    k = 3.0 / (2.0 * math.gamma(4.0 / 3.0) * 9.0 ** (1.0 / 3.0))
    return float(k * c0 * (d * d * SHEAR * length * length) ** (1.0 / 3.0))


def section(
    d: float, k_on: float, sites: float, k_off: float = 0.0, linear_shear: bool = False
) -> ChannelSection:
    x0, length = 100.0, 500.0
    reg = ReactiveRegion(x0, x0 + length, sites, (k_on,), (k_off,))
    return ChannelSection(
        x_grid(x0, length, 200.0),
        graded_edges(H, 60, 0.05),
        (lambda z: SHEAR * z) if linear_shear else plane_poiseuille,
        [Species("protein", d)],
        [reg],
    )


def test_graded_edges() -> None:
    e = graded_edges(100.0, 40, 0.1)
    assert e[0] == 0.0
    assert e[-1] == pytest.approx(100.0)
    assert np.diff(e)[0] == pytest.approx(0.1, rel=1e-9)
    assert np.all(np.diff(np.diff(e)) > 0)


def test_t12_leveque_transport_limited_flux() -> None:
    # Lévêque's assumption: linear shear profile; boundary layer (D L/γ̇)^(1/3) ≈ 7.5 µm.
    d, c0 = 50.0, 1.0
    s = section(d, math.inf, 1.0, linear_shear=True)
    j = s.surface_flux(s.steady_first_order(0, math.inf, c0)).sum()
    assert j == pytest.approx(leveque_flux(c0, d, 500.0), rel=0.01)


def test_t12_parabolic_profile_slightly_below_leveque() -> None:
    # The parabola is slower than its wall tangent at the edge of the layer: about 2 % less.
    d, c0 = 50.0, 1.0
    s = section(d, math.inf, 1.0)
    j = s.surface_flux(s.steady_first_order(0, math.inf, c0)).sum()
    ratio = j / leveque_flux(c0, d, 500.0)
    assert 0.96 < ratio < 1.0
    assert j < U_MEAN * H * c0


@pytest.mark.parametrize(("k", "limit"), [(1e-4, "reaction"), (1e3, "transport")])
def test_t13_damkoehler_limits(k: float, limit: str) -> None:
    d, c0, length = 50.0, 1.0, 500.0
    s = section(d, math.inf, 1.0)
    c = s.steady_first_order(0, k, c0)
    j = s.surface_flux(c).sum()
    j_r = k * c0 * length
    j_d = leveque_flux(c0, d, length)
    if limit == "reaction":
        assert j == pytest.approx(j_r, rel=0.02)
    else:
        assert j == pytest.approx(j_d, rel=0.04)


def test_t13_damkoehler_crossover() -> None:
    # Resistances in series (squires2008): J = J_R / (1 + Da), Da = J_R / J_D; approximate.
    d, c0, length = 50.0, 1.0, 500.0
    j_d = leveque_flux(c0, d, length)
    k = j_d / (c0 * length)  # Da = 1
    s = section(d, math.inf, 1.0)
    j = s.surface_flux(s.steady_first_order(0, k, c0)).sum()
    assert j == pytest.approx(0.5 * k * c0 * length, rel=0.15)


def test_t14_mass_conservation_two_species() -> None:
    reg = ReactiveRegion(100.0, 300.0, 0.05, (2e-3, 5e-4), (1e-3, 1e-2))
    s = ChannelSection(
        np.linspace(0.0, 600.0, 61),
        graded_edges(50.0, 20, 0.2),
        lambda z: 6.0 * 20.0 * z / 50.0 * (1 - z / 50.0),
        [Species("target", 40.0), Species("albumin", 60.0, site_weight=0.5)],
        [reg],
    )
    st = s.initial_state()
    for n in range(200):
        c_in = (2e-3, 5e-3) if n < 120 else (0.0, 0.0)
        st = s.step(st, 0.5, c_in)
    dissolved, bound = s.mass(st)
    np.testing.assert_allclose(st.inflow - st.outflow, dissolved + bound, rtol=1e-10)
    assert np.all(bound > 0)
    assert np.all(st.gamma >= -1e-15)
    used = st.gamma[0] + 0.5 * st.gamma[1]
    assert used.max() <= 0.05 * (1 + 1e-12)


def test_reaction_limited_langmuir_kinetics() -> None:
    # Fast transport, slow binding: Γ(t) = Γ_eq (1 − exp(−(k_on c₀ + k_off) t)) after the fill.
    d, c0, k_on, k_off, gmax = 500.0, 1e-3, 5.0, 2e-3, 1e-2
    reg = ReactiveRegion(50.0, 150.0, gmax, (k_on,), (k_off,))
    s = ChannelSection(
        np.linspace(0.0, 200.0, 41),
        np.linspace(0.0, 20.0, 11),
        lambda z: 6.0 * 200.0 * z / 20.0 * (1 - z / 20.0),
        [Species("target", d)],
        [reg],
    )
    st = s.initial_state()
    dt, n = 1.0, 600
    for _ in range(n):
        st = s.step(st, dt, (c0,))
    rate = k_on * c0 + k_off
    g_eq = gmax * k_on * c0 / rate
    m = s.region_of == 0
    ref = g_eq * (1 - math.exp(-rate * n * dt))
    # implicit Euler with dt·rate = 7·10⁻³ lags by about dt/2
    assert float(st.gamma[0, m].mean()) == pytest.approx(ref, rel=0.01)


def test_width_averaged_profile() -> None:
    f = DuctFlow(200.0, 100.0, 2.5e5)
    z = np.linspace(0, 100, 4001)
    q = np.trapezoid(f.width_averaged_velocity_um_s(z), z) * 200.0
    assert q == pytest.approx(2.5e5, rel=1e-5)


def test_t14_long_steps_respect_site_capacity() -> None:
    # Irreversible, fast binding with steps far beyond the filling time: sub-steps keep
    # Σ α Γ ≤ Γ_max and the balance exact.
    reg = ReactiveRegion(100.0, 300.0, 1e3, (5.0, 1.0), (0.0, 0.0))
    s = ChannelSection(
        np.linspace(0.0, 400.0, 41),
        graded_edges(50.0, 12, 0.5),
        lambda z: 6.0 * 20.0 * z / 50.0 * (1 - z / 50.0),
        [Species("target", 40.0), Species("albumin", 60.0, site_weight=2.0)],
        [reg],
    )
    st = s.initial_state()
    for _ in range(5):
        st = s.step(st, 600.0, (1.0, 5.0))
    used = st.gamma[0] + 2.0 * st.gamma[1]
    assert used.max() <= 1e3 * (1 + 1e-9)
    assert used[s.region_of == 0].min() > 0.99e3
    dissolved, bound = s.mass(st)
    np.testing.assert_allclose(st.inflow - st.outflow, dissolved + bound, rtol=1e-10)
