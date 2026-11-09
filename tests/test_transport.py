"""Stage Z3: cells, channel flow and particles (T7–T11, T27, T28, T31; docs/testfaelle.md)."""

import math

import numpy as np
import pytest
from scipy.integrate import quad

from bindoptix import rng
from bindoptix.cells import CellModel
from bindoptix.flow import DuctFlow
from bindoptix.particles import (
    CAPTURED,
    STUCK,
    ChannelModel,
    Spot,
    inlet_positions,
    run,
)

SEED = 20261111


def sphere(radius_um: float, **kw: float) -> CellModel:
    return CellModel(width_um=2 * radius_um, length_um=2 * radius_um, **kw)


def big_box(q: float = 0.0) -> DuctFlow:
    return DuctFlow(width_um=1e9, height_um=1e9, flow_rate_um3_s=q)


# --- T27: cell transport properties ----------------------------------------------------------


def test_t27_sphere_limit_stokes_einstein() -> None:
    c = sphere(0.5, delta_rho_kg_m3=100.0)
    kt = 1.380649e-23 * 293.15
    d_ref = kt / (6 * math.pi * 1.002e-3 * 0.5e-6) * 1e12
    v_ref = 2 / 9 * 100.0 * 9.80665 * (0.5e-6) ** 2 / 1.002e-3 * 1e6
    assert c.diffusion_um2_s == pytest.approx(d_ref, rel=1e-12)
    assert c.sedimentation_speed_um_s() == pytest.approx(v_ref, rel=1e-12)


@pytest.mark.parametrize("ratio", [1.2, 2.0, 5.0, 20.0])
def test_t27_perrin_orientation_average(ratio: float) -> None:
    # Independent closed form (Perrin): f/f₀ = √(1 − p²) / (p^(2/3) ln((1 + √(1 − p²))/p)),
    # p = b/a, f₀ friction of the sphere of equal volume, for the averaged mobility.
    c = CellModel(width_um=1.0, length_um=ratio)
    p = 1.0 / ratio
    s = math.sqrt(1 - p * p)
    f_ratio = s / (p ** (2 / 3) * math.log((1 + s) / p))
    f0 = 6 * math.pi * 1.002e-3 * c.equivalent_radius_um * 1e-6
    assert 1.0 / c.mobility_s_kg == pytest.approx(f0 * f_ratio, rel=1e-12)


def test_t27_legionella_nominal_values() -> None:
    # Values of docs/abschaetzung-einfang.md
    c = CellModel()
    assert c.sedimentation_speed_um_s() == pytest.approx(0.0255, rel=0.01)
    assert c.diffusion_um2_s == pytest.approx(0.401, rel=0.01)


# --- flow ------------------------------------------------------------------------------------


def test_flow_rate_from_series_matches_q() -> None:
    f = DuctFlow(200.0, 100.0, 2.5e5)
    y_int = quad(lambda y: float(f.column_flux_um2_s([y])[0]), -100.0, 100.0, epsrel=1e-12)[0]
    assert y_int == pytest.approx(2.5e5, rel=1e-6)
    zz = np.linspace(0, 100, 2001)
    q0 = np.trapezoid(f.velocity_um_s(np.zeros_like(zz), zz), zz)
    assert q0 == pytest.approx(float(f.column_flux_um2_s([0.0])[0]), rel=1e-6)


def test_flow_wide_channel_limit_is_plane_poiseuille() -> None:
    # Far from the side walls: u = G z (h − z)/2, ∂u/∂z(0) = G h/2; Q → h³wG/12.
    f = DuctFlow(1e5, 100.0, 1e5 * 100.0 * 10.0)
    g = f.pressure_gradient
    z = np.array([25.0, 50.0, 75.0])
    np.testing.assert_allclose(f.velocity_um_s(np.zeros(3), z), g * z * (100 - z) / 2, rtol=1e-6)
    assert float(f.floor_shear_rate_s([0.0])[0]) == pytest.approx(g * 50.0, rel=1e-6)
    # side-wall layers of width ~h reduce Q by about 0.63 h/w
    assert 1 - f.shape_factor == pytest.approx(0.630 * 100.0 / 1e5, rel=0.01)


def test_floor_shear_against_finite_difference() -> None:
    # Reference: finite difference of a long series (5000 terms, tail of Σ n⁻² ≈ 5·10⁻⁵).
    f = DuctFlow(200.0, 100.0, 2.5e5)
    ref_flow = DuctFlow(200.0, 100.0, 2.5e5, n_terms=5000)
    y = np.array([-80.0, -30.0, 0.0, 55.0])
    dz = 1e-4
    fd = ref_flow.velocity_um_s(y, np.full(4, dz)) / dz
    np.testing.assert_allclose(f.floor_shear_rate_s(y), fd, rtol=2e-4)


def test_flow_vanishes_on_walls() -> None:
    f = DuctFlow(200.0, 100.0, 1e5)
    assert np.abs(f.velocity_um_s([100.0, -100.0, 0.0], [50.0, 50.0, 0.0])).max() < 1e-9


# --- T11: counter-based random numbers and determinism ---------------------------------------


def test_t11_random_numbers_independent_of_batching() -> None:
    ids = np.arange(1000)
    whole = rng.normal(SEED, ids, 7, 2)
    parts = np.concatenate([rng.normal(SEED, ids[:333], 7, 2), rng.normal(SEED, ids[333:], 7, 2)])
    np.testing.assert_array_equal(whole, parts)
    big = rng.normal(SEED, np.arange(200_000), 1, 0)
    assert abs(big.mean()) < 4 / math.sqrt(2e5)
    assert big.var() == pytest.approx(1.0, abs=4 * math.sqrt(2 / 2e5))


def test_t11_run_independent_of_batching() -> None:
    m = ChannelModel(DuctFlow(200.0, 100.0, 2e4), CellModel(swim_speed_um_s=5.0))
    start = inlet_positions(m.flow, 60, SEED)
    a = run(m, start, 20.0, 0.1, SEED)
    b1 = run(m, start[:25], 20.0, 0.1, SEED, ids=np.arange(25))
    b2 = run(m, start[25:], 20.0, 0.1, SEED, ids=np.arange(25, 60))
    # Random streams are bit-identical; floating-point results agree to rounding (NumPy uses
    # different SIMD paths for different array lengths), discrete states exactly.
    np.testing.assert_allclose(
        a.position_um, np.vstack([b1.position_um, b2.position_um]), rtol=1e-12, atol=1e-12
    )
    np.testing.assert_array_equal(a.state, np.concatenate([b1.state, b2.state]))
    c = run(m, start, 20.0, 0.1, SEED)
    np.testing.assert_array_equal(a.position_um, c.position_um)


# --- T7: free diffusion ----------------------------------------------------------------------


def test_t7_mean_square_displacement() -> None:
    m = ChannelModel(big_box(), sphere(0.5, delta_rho_kg_m3=0.0), floor="reflect")
    n, t = 20_000, 10.0
    start = np.tile([0.0, 0.0, 5e8], (n, 1))
    r = run(m, start, t, 0.5, SEED)
    msd = np.sum((r.position_um - start) ** 2, axis=1)
    d = m.cell.diffusion_um2_s
    se = math.sqrt(2 / 3) * 6 * d * t / math.sqrt(n)  # var of |Δx|² = 2/3 (6Dt)²
    assert abs(msd.mean() - 6 * d * t) < 4 * se


# --- T8: sedimentation -----------------------------------------------------------------------


def test_t8_drift_and_mason_weaver_equilibrium() -> None:
    # Test cell with larger Δρ so that the equilibrium length ℓ = D/v_s fits into 60 µm.
    cell = sphere(0.5, delta_rho_kg_m3=500.0)
    d, v = cell.diffusion_um2_s, cell.sedimentation_speed_um_s()
    ell, h = d / v, 60.0
    m = ChannelModel(DuctFlow(1e9, h, 0.0), cell, floor="reflect")
    n = 4000
    g = np.random.default_rng(SEED)
    start = np.column_stack([np.zeros(n), np.zeros(n), g.uniform(0, h, n)])
    t_eq = 8 * h / v
    r = run(m, start, t_eq, 0.5, SEED)
    z = r.position_um[:, 2]
    mean_ref = ell - h / math.expm1(h / ell)
    var_ref = ell**2 - h**2 * math.exp(h / ell) / math.expm1(h / ell) ** 2
    assert abs(z.mean() - mean_ref) < 4 * math.sqrt(var_ref / n)


def test_t8_mean_sedimentation_speed() -> None:
    cell = CellModel()
    m = ChannelModel(big_box(), cell, floor="reflect")
    n, t = 20_000, 200.0
    start = np.tile([0.0, 0.0, 5e8], (n, 1))
    r = run(m, start, t, 1.0, SEED)
    dz = r.position_um[:, 2] - start[:, 2]
    se = math.sqrt(2 * cell.diffusion_um2_s * t / n)
    assert abs(dz.mean() + cell.sedimentation_speed_um_s() * t) < 4 * se


# --- T9: run and tumble ----------------------------------------------------------------------


def test_t9_run_and_tumble_effective_diffusion() -> None:
    cell = sphere(0.5, delta_rho_kg_m3=0.0, swim_speed_um_s=20.0, run_time_s=1.0, tumble_cos=0.33)
    m = ChannelModel(big_box(), cell, floor="reflect", brownian=False)
    n = 4000
    start = np.tile([0.0, 0.0, 5e8], (n, 1))
    r = run(m, start, 60.0, 0.02, SEED, snapshot_times_s=(30.0, 60.0))
    msd = {t: np.sum((r.snapshots_um[t] - start) ** 2, axis=1).mean() for t in (30.0, 60.0)}
    slope = (msd[60.0] - msd[30.0]) / 30.0 / 6.0
    d_ref = cell.effective_diffusion_um2_s
    # statistical spread of the slope estimate ~ 6 %; discretisation of the runs small
    assert slope == pytest.approx(d_ref, rel=0.12)


# --- T10: first passage to an absorbing floor ------------------------------------------------


def test_t10_smoluchowski_flux() -> None:
    m = ChannelModel(
        DuctFlow(1e9, 1e9, 0.0), sphere(0.5, delta_rho_kg_m3=0.0), nonspecific_rate_s=math.inf
    )
    d = m.cell.diffusion_um2_s
    n, height, t = 40_000, 40.0, 100.0
    g = np.random.default_rng(SEED)
    start = np.column_stack([np.zeros(n), np.zeros(n), g.uniform(0, height, n)])
    r = run(m, start, t, 0.25, SEED)
    frac = r.count(STUCK) / n
    ref = 2 * math.sqrt(d * t / math.pi) / height
    assert abs(frac - ref) < 4 * math.sqrt(ref * (1 - ref) / n) + 0.02 * ref


# --- T28: Hazen limit with the exact column fluxes -------------------------------------------


def test_t28_hazen_capture_fraction_in_duct() -> None:
    cell = CellModel()
    vs = cell.sedimentation_speed_um_s()
    flow = DuctFlow(200.0, 100.0, 2.5e4)  # 1.5 nL/min
    length = 4000.0
    m = ChannelModel(flow, cell, nonspecific_rate_s=math.inf, x_out_um=length, brownian=False)
    n = 6000
    start = inlet_positions(flow, n, SEED)
    t_max = 1.2 * 100.0 / vs
    r = run(m, start, t_max, 5.0, SEED)
    frac = r.count(STUCK) / n
    # Exact: column y is captured up to min(q(y), L v_s) of its flux q(y).
    ref = (
        quad(lambda y: min(float(flow.column_flux_um2_s([y])[0]), length * vs), -100.0, 100.0)[0]
        / flow.flow_rate_um3_s
    )
    assert 0.05 < ref < 0.95
    assert abs(frac - ref) < 4 * math.sqrt(ref * (1 - ref) / n) + 0.01


# --- T31: sliding on the floor --------------------------------------------------------------


def test_t31_sliding_speed_on_the_floor() -> None:
    # Cells resting on the floor move with u_wall = s · a · ∂u/∂z(0) (no Brownian motion).
    cell = CellModel()
    flow = DuctFlow(200.0, 100.0, 2.5e5)
    m = ChannelModel(flow, cell, nonspecific_rate_s=0.0, slip_factor=0.5, brownian=False)
    y = np.array([-90.0, -40.0, 0.0, 20.0, 70.0])
    start = np.column_stack([np.zeros(5), y, np.zeros(5)])
    t = 600.0
    r = run(m, start, t, 1.0, SEED)
    dz = 1e-4
    ref_flow = DuctFlow(200.0, 100.0, 2.5e5, n_terms=5000)
    shear = ref_flow.velocity_um_s(y, np.full(5, dz)) / dz
    ref = 0.5 * cell.equivalent_radius_um * shear * t
    np.testing.assert_allclose(r.position_um[:, 0], ref, rtol=3e-4)
    np.testing.assert_allclose(r.position_um[:, 2], 0.0)
    # u_wall relative to the mean speed: 6 s a / h in a wide channel, ~1.3 % here
    ratio = float(ref[2] / t / flow.mean_speed_um_s)
    assert 0.008 < ratio < 0.03


def test_t31_sliding_cell_captured_on_arrival() -> None:
    cell = CellModel()
    flow = DuctFlow(200.0, 100.0, 2.5e5)
    x_spot = 50.0
    spot = Spot(
        np.array(
            [[x_spot, -100.0], [x_spot + 200, -100.0], [x_spot + 200, 100.0], [x_spot, 100.0]]
        ),
        math.inf,
    )
    m = ChannelModel(flow, cell, nonspecific_rate_s=0.0, spots=(spot,), brownian=False)
    r = run(m, np.array([[0.0, 0.0, 0.0]]), 400.0, 0.5, SEED)
    u_wall = 0.5 * cell.equivalent_radius_um * float(flow.floor_shear_rate_s([0.0])[0])
    assert r.state[0] == CAPTURED
    assert r.spot[0] == 0
    assert r.event_time_s[0] == pytest.approx(x_spot / u_wall, abs=0.5)
