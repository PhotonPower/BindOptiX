"""Z4, part 2: 3D channel volume (T35, docs/testfaelle.md)."""

import math

import numpy as np
import pytest
from numpy.typing import NDArray

from bindoptix.continuum import ChannelSection, ReactiveRegion, Species, graded_edges
from bindoptix.continuum3d import ChannelVolume, FloorRegion
from bindoptix.flow import DuctFlow

W, H = 200.0, 100.0


def rect(x0: float, x1: float, y0: float, y1: float) -> NDArray[np.float64]:
    return np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]], dtype=np.float64)


def y_edges(n: int) -> NDArray[np.float64]:
    """Symmetric edges on [−w/2, w/2], refined towards the side walls."""
    half = graded_edges(W / 2, n // 2, 1.0)
    return np.asarray(np.concatenate([-half[::-1], half[1:]]), dtype=np.float64)


def x_edges(x0: float, length: float) -> NDArray[np.float64]:
    up = x0 - graded_edges(x0, 15, 1.0)[::-1]
    on = x0 + graded_edges(length, 50, 0.5)
    down = x0 + length + graded_edges(150.0, 8, 5.0)
    return np.asarray(np.unique(np.concatenate([up, on, down])), dtype=np.float64)


def test_t35_laterally_uniform_equals_2d() -> None:
    # u independent of y and a full-width spot: the 3D solution is the 2D solution times w.
    def u_z(z: NDArray[np.float64]) -> NDArray[np.float64]:
        return np.asarray(6.0 * 50.0 * z / H * (1 - z / H), dtype=np.float64)

    xe, ze = x_edges(100.0, 300.0), graded_edges(H, 24, 0.3)
    sp2 = Species("protein", 50.0)
    two = ChannelSection(xe, ze, u_z, [sp2], [ReactiveRegion(100.0, 400.0, 1.0, (2.0,), (0.0,))])
    three = ChannelVolume(
        xe,
        np.linspace(-W / 2, W / 2, 5),
        ze,
        lambda y, z: u_z(z),
        [sp2],
        [FloorRegion(rect(100.0, 400.0, -W, W), 1.0, (2.0,), (0.0,))],
    )
    j2 = two.surface_flux(two.steady_first_order(0, 2.0, 1.0)).sum() * W
    j3 = three.surface_flux(three.steady_first_order(0, 2.0, 1.0)).sum()
    assert j3 == pytest.approx(j2, rel=1e-10)


def test_t35_local_leveque_in_duct() -> None:
    # Absorbing full-width spot in the duct: thin layer → local Lévêque flux with the local
    # wall shear, J ≈ 0.8075 c₀ (D² L²)^(1/3) ∫ γ̇(y)^(1/3) dy (lateral diffusion neglected).
    d, c0, length = 50.0, 1.0, 300.0
    flow = DuctFlow(W, H, W * H * 1000.0)
    vol = ChannelVolume(
        x_edges(100.0, length),
        y_edges(20),
        graded_edges(H, 36, 0.05),
        lambda y, z: flow.velocity_um_s(y, z),
        [Species("protein", d)],
        [FloorRegion(rect(100.0, 100.0 + length, -W, W), 1.0, (math.inf,), (0.0,))],
    )
    j = vol.surface_flux(vol.steady_first_order(0, math.inf, c0)).sum()
    yy = np.linspace(-W / 2, W / 2, 4001)
    g13 = np.trapezoid(flow.floor_shear_rate_s(yy) ** (1 / 3), yy)
    k = 3.0 / (2.0 * math.gamma(4.0 / 3.0) * 9.0 ** (1.0 / 3.0))
    ref = k * c0 * (d * d * length * length) ** (1 / 3) * g13
    assert j == pytest.approx(ref, rel=0.04)


def test_t35_mass_conservation_and_symmetry() -> None:
    flow = DuctFlow(W, H, W * H * 20.0)
    sp_t, sp_a = Species("target", 40.0), Species("albumin", 60.0, site_weight=0.5)
    vol = ChannelVolume(
        np.linspace(0.0, 500.0, 41),
        y_edges(12),
        graded_edges(H, 14, 0.5),
        lambda y, z: flow.velocity_um_s(y, z),
        [sp_t, sp_a],
        [FloorRegion(rect(150.0, 350.0, -60.0, 60.0), 0.05, (2e-3, 5e-4), (1e-3, 1e-2))],
    )
    st = vol.initial_state()
    for n in range(40):
        st = vol.step(st, 2.0, (2e-3, 5e-3) if n < 25 else (0.0, 0.0))
    dissolved, bound = vol.mass(st)
    np.testing.assert_allclose(st.inflow - st.outflow, dissolved + bound, rtol=1e-9)
    assert np.all(bound > 0)
    np.testing.assert_allclose(st.gamma, st.gamma[:, :, ::-1], rtol=1e-8, atol=1e-14)
    assert (st.gamma[0] + 0.5 * st.gamma[1]).max() <= 0.05 * (1 + 1e-12)


def test_t35_chip_channel_volume() -> None:
    from bindoptix.chip import ChipLayout
    from bindoptix.chipsim import channel_volume

    layout = ChipLayout(grid_angle_rad=math.atan(1 / 5))
    kon = [[1.0]] * 25
    koff = [[0.0]] * 25
    vol, cs = channel_volume(layout, 2, 5e6 / 60, [Species("p", 50.0)], 1.0, kon, koff)
    assert len(vol.regions) == 5
    assert cs.grid_index == (10, 11, 12, 13, 14)
    for j in range(5):
        reactive = vol.floor_area[vol.region_of == j].sum()
        # spots are axis-aligned in the channel frame and the x grid contains their edges;
        # the y grid (16.7 µm) resolves the 200 µm width exactly (full-width spots)
        assert reactive == pytest.approx(200.0 * 200.0, rel=1e-9)
