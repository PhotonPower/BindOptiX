"""Stage Z5: binding kinetics, cell capture, receptor states (T15, T16, T36)."""

import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from bindoptix.binding import AvidityCapture, CompetitiveLangmuir, Occupancy
from bindoptix.cells import CellModel
from bindoptix.chip import ChipLayout
from bindoptix.chipsim import CoupledRun, channel_model, run_channel_coupled
from bindoptix.continuum import Species
from bindoptix.frames import Mount
from bindoptix.particles import CAPTURED, inlet_positions, run

SEED = 20261114
THETA = math.atan(1 / 5)


# --- T15: competitive Langmuir ---------------------------------------------------------------


def test_t15_single_species_analytic() -> None:
    m = CompetitiveLangmuir(1e4, (2e-3,), (1e-3,), (1.0,))
    c, t = 0.5, 400.0
    rate = 2e-3 * c + 1e-3
    ref = 1e4 * 2e-3 * c / rate * (1 - math.exp(-rate * t))
    assert m.integrate([0.0], [c], t)[0] == pytest.approx(ref, rel=1e-12)


def test_t15_competitive_equilibrium_closed_form() -> None:
    m = CompetitiveLangmuir(1e4, (2e-3, 5e-4, 1e-3), (1e-3, 1e-2, 2e-4), (1.0, 0.5, 2.0))
    c = [0.3, 4.0, 0.05]
    eq = m.equilibrium(c)
    np.testing.assert_allclose(m.rhs(eq, c), 0.0, atol=1e-9)
    # independent: K_i c_i θ with θ from the site balance
    k = np.array([2e-3 * 0.3 / 1e-3, 5e-4 * 4.0 / 1e-2, 1e-3 * 0.05 / 2e-4])
    theta = 1e4 / (1 + np.array([1.0, 0.5, 2.0]) @ k)
    np.testing.assert_allclose(eq, k * theta, rtol=1e-13)
    np.testing.assert_allclose(m.integrate([0, 0, 0], c, 1e6), eq, rtol=1e-9)


def test_t15_competitive_transient_against_ode_solver() -> None:
    m = CompetitiveLangmuir(1e4, (2e-3, 5e-4), (1e-3, 1e-2), (1.0, 0.5))
    c, g0 = [0.3, 4.0], np.array([100.0, 2000.0])
    sol = solve_ivp(lambda t, g: m.rhs(g, c), (0.0, 300.0), g0, rtol=1e-11, atol=1e-9)
    np.testing.assert_allclose(m.integrate(g0, c, 300.0), sol.y[:, -1], rtol=1e-8)


def test_t15_irreversible_has_no_equilibrium() -> None:
    m = CompetitiveLangmuir(1e4, (2e-3,), (0.0,), (1.0,))
    with pytest.raises(ValueError, match="irreversible"):
        m.equilibrium([1.0])
    assert m.integrate([0.0], [1.0], 1e5)[0] == pytest.approx(1e4, rel=1e-9)


# --- T16: receptor states (port B) -------------------------------------------------------------


def occupancy(fr: np.ndarray) -> Occupancy:  # type: ignore[type-arg]
    return Occupancy(np.array([0.0, 60.0]), (3, 4), ("free", "albumin", "cell"), fr)


def test_t16_states_sum_to_one_and_round_trip() -> None:
    fr = np.array([[[1, 0, 0], [1, 0, 0]], [[0.7, 0.25, 0.05], [0.9, 0.1, 0.0]]], dtype=float)
    occ = Occupancy(
        np.array([0.0, 60.0]),
        (3, 4),
        ("free", "albumin", "cell"),
        fr,
        np.array([42.0]),
        np.array([[1.0, -2.0]]),
        np.array([3]),
        ("Lp-SG1",),
    )
    back = Occupancy.from_json(occ.to_json())
    np.testing.assert_array_equal(back.fractions, fr)
    assert back.cell_strain == ("Lp-SG1",)
    np.testing.assert_array_equal(back.cell_xy_um, [[1.0, -2.0]])


@pytest.mark.parametrize("bad", ["sum", "range", "first"])
def test_t16_invalid_states(bad: str) -> None:
    fr = np.zeros((2, 2, 3))
    fr[..., 0] = 1.0
    if bad == "sum":
        fr[1, 0, 1] = 0.2
        with pytest.raises(ValueError, match="sum"):
            occupancy(fr)
    elif bad == "range":
        fr[1, 0] = [1.2, -0.2, 0.0]
        with pytest.raises(ValueError, match=r"\[0, 1\]"):
            occupancy(fr)
    else:
        with pytest.raises(ValueError, match="free"):
            Occupancy(np.array([0.0, 1.0]), (3, 4), ("albumin", "free", "cell"), fr)


# --- T36: coupling fouling → capture -----------------------------------------------------------

COARSE = {"dx_um": 250.0, "ny": 6, "nz": 8, "first_dz_um": 2.0}


def test_t36_without_fouling_equals_constant_rate() -> None:
    layout = ChipLayout(grid_angle_rad=THETA)
    cap = AvidityCapture(bond_rate_um2_s=2e-4, contact_area_um2=0.5)
    gmax = 1e4
    q, n, dur = 5e6 / 60.0, 200, 3 * 3600.0
    res = run_channel_coupled(layout, 2, CellModel(), q, cap, gmax, n, dur, 10.0, 600.0, SEED)
    model, _ = channel_model(
        layout, 2, CellModel(), q, Mount.sensor_floor(), float(cap.capture_rate_s(gmax))
    )
    start = inlet_positions(model.flow, n, SEED + 2)
    ref = run(model, start, dur, 10.0, SEED, ids=np.arange(n) + 2_000_000_000)
    np.testing.assert_array_equal(res.particles.state, ref.state)
    np.testing.assert_array_equal(res.particles.spot, ref.spot)
    occ = res.occupancy
    assert occ.states == ("free", "cell")
    assert occ.cell_time_s.size == ref.count(CAPTURED)
    assert np.all(occ.fractions[-1, :, 1] >= 0.0)


def test_t36_fouling_reduces_capture() -> None:
    layout = ChipLayout(grid_angle_rad=THETA)
    cap = AvidityCapture(bond_rate_um2_s=2e-4)
    albumin = Species("albumin", 60.0, site_weight=1.0)

    def go(fouled: bool) -> CoupledRun:
        return run_channel_coupled(
            layout,
            2,
            CellModel(),
            5e6 / 60.0,
            cap,
            1e4,
            300,
            3 * 3600.0,
            10.0,
            300.0,
            SEED,
            fouling=(albumin,) if fouled else (),
            fouling_k_on_um3_s=(1.0,) if fouled else (),
            fouling_k_off_s=(0.0,) if fouled else (),
            fouling_c_in=(10.0,) if fouled else (),
            grid=COARSE,
        )

    clean, fouled = go(False), go(True)
    assert fouled.particles.count(CAPTURED) < clean.particles.count(CAPTURED)
    th = fouled.free_sites_per_um2
    times = fouled.times_s
    assert th[-1].max() < 1e-3 * 1e4  # all spots fouled at the end
    # no capture after the spots are practically full (free sites < 1 %)
    t_full = times[np.argmax(th.max(axis=1) < 100.0)]
    occ = fouled.occupancy
    assert np.all(occ.cell_time_s <= t_full + 300.0)
    fr = occ.fractions
    assert occ.states == ("free", "albumin", "cell")
    np.testing.assert_allclose(fr.sum(axis=2), 1.0, atol=1e-9)
    assert fr[-1, :, 1].min() > 0.99
    # free sites fall first on the upstream spot (fouling front moves downstream)
    first_drop = [int(np.argmax(th[:, j] < 5e3)) for j in range(5)]
    assert first_drop == sorted(first_drop)


def test_t36_affinity_row_of_binding_matrix() -> None:
    # A strain that binds only the aptamers of spots 11 and 13 is caught only there.
    layout = ChipLayout(grid_angle_rad=THETA)
    aff = np.zeros(25)
    aff[[11, 13]] = 1.0
    res = run_channel_coupled(
        layout,
        2,
        CellModel(),
        5e6 / 60.0,
        AvidityCapture(bond_rate_um2_s=2e-4),
        1e4,
        300,
        3 * 3600.0,
        10.0,
        600.0,
        SEED,
        strain="Lp-SG1",
        affinity_per_spot=aff,
        nonspecific_rate_s=0.0,
    )
    spots = set(res.occupancy.cell_spot.tolist())
    assert spots <= {11, 13}
    assert len(spots) > 0
    assert set(res.occupancy.cell_strain) == {"Lp-SG1"}
