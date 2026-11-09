"""Optical maps as data ports (Z2): T6 and the interpolation, evaluation and file format."""

import math

import numpy as np
import pytest
from numpy.typing import NDArray

from bindoptix.ports import (
    EmitterResponseTable,
    FieldMap,
    FieldMapStack,
    RegularGrid,
    orientation_average,
    quadratic,
    unit_vectors,
)
from bindoptix.results import decode_array, encode_array, load_json, to_json

SEED = 20261110


def icosahedron() -> NDArray[np.float64]:
    """Vertices of the icosahedron: a spherical 5-design, exact for quadratic forms."""
    g = (1.0 + math.sqrt(5.0)) / 2.0
    v = []
    for a in (-1.0, 1.0):
        for b in (-g, g):
            v += [(0.0, a, b), (a, b, 0.0), (b, 0.0, a)]
    return unit_vectors(v)


def random_psd(rng: np.random.Generator, shape: tuple[int, ...]) -> NDArray[np.float64]:
    b = rng.normal(size=(*shape, 3, 3))
    return np.asarray(b @ np.swapaxes(b, -1, -2), dtype=np.float64)


def make_table(rng: np.random.Generator, n_det: int = 2) -> EmitterResponseTable:
    """Consistent random table: A_det,d ≤ A_rad ≤ A_tot by construction."""
    grid = RegularGrid((np.array([0.0, 10.0]), np.array([0.0, 5.0, 10.0]), np.array([0.01, 0.1])))
    g = grid.shape
    det = random_psd(rng, (*g, n_det))
    rad = det.sum(axis=3) + random_psd(rng, g)
    tot = rad + random_psd(rng, g)
    return EmitterResponseTable(
        grid, 0.67, 1.333, tuple(f"d{k}" for k in range(n_det)), tot, rad, det, {"solver": "test"}
    )


# --- T6 ------------------------------------------------------------------------------------


def test_t6_orientation_average_is_trace_over_three() -> None:
    a = random_psd(np.random.default_rng(SEED), (20,))
    mu = icosahedron()
    avg = np.array([quadratic(np.repeat(ai[None], 12, axis=0), mu).mean() for ai in a])
    np.testing.assert_allclose(orientation_average(a), avg, rtol=1e-13)


def test_t6_valid_table_accepted() -> None:
    make_table(np.random.default_rng(SEED))


@pytest.mark.parametrize("which", ["asym", "indefinite", "rad_gt_tot", "det_gt_rad"])
def test_t6_invalid_tables_rejected(which: str) -> None:
    t = make_table(np.random.default_rng(SEED))
    tot, rad, det = t.a_tot.copy(), t.a_rad.copy(), t.a_det.copy()
    if which == "asym":
        tot[0, 0, 0, 0, 1] += 1.0
        match = "symmetric"
    elif which == "indefinite":
        det[0, 0, 0, 0] = np.diag([1.0, 1.0, -1.0])
        match = "semidefinite"
    elif which == "rad_gt_tot":
        rad = tot * 1.5
        match = "absorbed"
    else:
        det = det * 0 + rad[..., None, :, :]
        match = "missing"
    with pytest.raises(ValueError, match=match):
        EmitterResponseTable(t.grid, 0.67, 1.333, t.detectors, tot, rad, det)


# --- rates ---------------------------------------------------------------------------------


def free_dipole_table(share: float) -> EmitterResponseTable:
    grid = RegularGrid((np.array([0.0]), np.array([0.0]), np.array([0.0, 1.0])))
    eye = np.broadcast_to(np.eye(3), (1, 1, 2, 3, 3)).copy()
    return EmitterResponseTable(grid, 0.67, 1.333, ("cam",), eye, eye, share * eye[..., None, :, :])


@pytest.mark.parametrize("q0", [0.1, 0.5, 0.9, 1.0])
def test_free_dipole_keeps_intrinsic_yield(q0: float) -> None:
    t = free_dipole_table(0.3)
    p = np.array([[123.0, -45.0, 0.5]])
    mu = unit_vectors([[1.0, 2.0, 3.0]])
    assert t.quantum_yield(p, mu, q0)[0] == pytest.approx(q0, rel=1e-15)
    assert t.collection(p, mu)[0, 0] == pytest.approx(0.3, rel=1e-15)


def test_quantum_yield_from_rates() -> None:
    # Independent: q = γ_r / (γ_r + γ_abs + γ_nr⁰) with γ in units of γ_r⁰.
    grid = RegularGrid((np.array([0.0]), np.array([0.0]), np.array([0.0])))
    rad = np.diag([2.0, 0.5, 0.5]).reshape(1, 1, 1, 3, 3)
    tot = rad + np.diag([40.0, 10.0, 10.0]).reshape(1, 1, 1, 3, 3)
    t = EmitterResponseTable(grid, 0.67, 1.333, ("cam",), tot, rad, 0.2 * rad[..., None, :, :])
    q0 = 0.6
    g_r, g_abs, g_nr0 = 2.0, 40.0, (1 - q0) / q0
    q = t.quantum_yield(np.zeros((1, 3)), [[1.0, 0.0, 0.0]], q0)[0]
    assert q == pytest.approx(g_r / (g_r + g_abs + g_nr0), rel=1e-14)


# --- interpolation -------------------------------------------------------------------------


def test_interpolation_exact_for_linear_fields() -> None:
    grid = RegularGrid((np.linspace(0, 10, 4), np.linspace(-5, 5, 3), np.array([0.0, 0.3, 2.0])))
    x, y, z = np.meshgrid(*grid.axes_um, indexing="ij")
    f = 1.5 + 2.0 * x - 0.5 * y + 7.0 * z
    rng = np.random.default_rng(SEED)
    p = np.column_stack([rng.uniform(0, 10, 50), rng.uniform(-5, 5, 50), rng.uniform(0, 2, 50)])
    ref = 1.5 + 2.0 * p[:, 0] - 0.5 * p[:, 1] + 7.0 * p[:, 2]
    np.testing.assert_allclose(grid.interpolate(f, p), ref, rtol=1e-13)


def test_interpolation_keeps_psd_and_symmetry() -> None:
    t = make_table(np.random.default_rng(SEED))
    rng = np.random.default_rng(SEED + 1)
    p = np.column_stack(
        [rng.uniform(0, 10, 200), rng.uniform(0, 10, 200), rng.uniform(0.01, 0.1, 200)]
    )
    tot, rad, det = t.tensors_at(p)
    for a in (tot, rad, det, tot - rad):
        assert np.max(np.abs(a - np.swapaxes(a, -1, -2))) < 1e-12 * np.abs(a).max()
        assert np.linalg.eigvalsh(a).min() > -1e-9 * np.abs(a).max()


def test_invariant_axis_and_outside() -> None:
    t = free_dipole_table(0.3)
    t.tensors_at([[1e4, -1e4, 0.5]])  # x and y invariant
    with pytest.raises(ValueError, match="outside"):
        t.tensors_at([[0.0, 0.0, 2.0]])


# --- field map -----------------------------------------------------------------------------


def plane_wave_map(nz: int) -> FieldMap:
    lam, n = 0.6, 1.333
    k = 2 * math.pi * n / lam
    z = np.linspace(0.0, 1.0, nz)
    grid = RegularGrid((np.array([0.0]), np.array([0.0]), z))
    e = np.zeros((1, 1, nz, 3), dtype=np.complex128)
    e[0, 0, :, 0] = np.exp(1j * k * z)
    return FieldMap(grid, lam, e)


def test_field_map_plane_wave_converges() -> None:
    # Linear interpolation of exp(ikz): |E|² error O((k Δz)²); halving Δz quarters it.
    p = np.column_stack([np.zeros(101), np.zeros(101), np.linspace(0.0, 1.0, 101)])
    mu = np.tile([1.0, 0.0, 0.0], (101, 1))
    errs = [np.abs(plane_wave_map(n).excitation(p, mu) - 1.0).max() for n in (41, 81, 161)]
    assert errs[0] < 0.05
    assert errs[1] / errs[2] == pytest.approx(4.0, rel=0.1)


def test_field_map_scales_with_power_and_guards_nonlinear() -> None:
    m = plane_wave_map(41)
    p, mu = [[0.0, 0.0, 0.0]], [[1.0, 0.0, 0.0]]
    assert m.excitation(p, mu, 2.5)[0] == pytest.approx(2.5, rel=1e-14)
    nl = FieldMap(m.grid, m.wavelength_um, m.e, nonlinear=True, power_mw=3.0)
    with pytest.raises(ValueError, match="nonlinear"):
        nl.excitation(p, mu, 2.5)
    with pytest.raises(ValueError, match="needs the power"):
        FieldMap(m.grid, m.wavelength_um, m.e, nonlinear=True)


def test_nested_maps_use_finest() -> None:
    coarse = FieldMap(
        RegularGrid((np.linspace(-10, 10, 3), np.linspace(-10, 10, 3), np.linspace(0, 10, 3))),
        0.6,
        np.ones((3, 3, 3, 3), dtype=np.complex128),
    )
    fine = FieldMap(
        RegularGrid((np.linspace(-1, 1, 5), np.linspace(-1, 1, 5), np.linspace(0, 1, 5))),
        0.6,
        2.0 * np.ones((5, 5, 5, 3), dtype=np.complex128),
    )
    stack = FieldMapStack((coarse, fine))
    e = stack.e_at([[0.0, 0.0, 0.5], [5.0, 0.0, 5.0]])
    np.testing.assert_allclose(e[:, 0], [2.0, 1.0])
    with pytest.raises(ValueError, match="outside"):
        stack.e_at([[50.0, 0.0, 0.0]])


# --- file format ---------------------------------------------------------------------------


def test_table_json_round_trip_bitwise() -> None:
    t = make_table(np.random.default_rng(SEED))
    back = EmitterResponseTable.from_json(t.to_json())
    for name in ("a_tot", "a_rad", "a_det"):
        np.testing.assert_array_equal(getattr(back, name), getattr(t, name))
    assert back.detectors == t.detectors
    assert back.source == t.source


def test_field_map_json_round_trip_bitwise() -> None:
    m = plane_wave_map(11)
    q = np.asarray(np.linspace(0, 1, 11).reshape(1, 1, 11), dtype=np.float64)
    m2 = FieldMap(m.grid, m.wavelength_um, m.e, q_w_m3=q, source={"solver": "hp-FEM"})
    back = FieldMap.from_json(m2.to_json())
    np.testing.assert_array_equal(back.e, m2.e)
    assert back.q_w_m3 is not None
    np.testing.assert_array_equal(back.q_w_m3, q)
    assert back.delta_t_k is None


def test_codec_nonfinite_and_envelope() -> None:
    a = np.array([1.0, math.nan, math.inf, -math.inf, -0.0])
    enc = encode_array(a)
    assert enc["values"][1:4] == ["NaN", "Infinity", "-Infinity"]
    b = decode_array(enc)
    np.testing.assert_array_equal(np.isnan(b), np.isnan(a))
    assert math.copysign(1.0, b[4]) == -1.0
    text = to_json("FieldMap", {"x": 1})
    assert load_json(text, "FieldMap") == {"x": 1}
    with pytest.raises(ValueError, match="expected type"):
        load_json(text, "EmitterResponseTable")
