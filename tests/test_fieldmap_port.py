"""Port F metadata and normalisation (Z2, issue #23): T42."""

import json
import math

import numpy as np
import pytest
from numpy.typing import NDArray

from bindoptix.ports import FieldMap, FieldMapStack, RegularGrid
from bindoptix.results import encode_array, load_json, to_json

# Independent constants (CODATA 2018), not taken from bindoptix.units.
C0 = 299_792_458.0
EPS0 = 8.8541878128e-12
Z0 = 1.0 / (EPS0 * C0)
N_WATER = 1.333


def plane_wave_grid() -> RegularGrid:
    return RegularGrid((np.array([0.0]), np.array([0.0]), np.linspace(0.0, 1.0, 5)))


def unit_field(grid: RegularGrid) -> NDArray[np.complex128]:
    """|E| = 1 V/m along x, as hp-FEM computes the excitation."""
    e = np.zeros((*grid.shape, 3), dtype=np.complex128)
    e[..., 0] = 1.0
    return e


def hpfem_map(**kwargs: object) -> FieldMap:
    """A plane wave of |E₀| = 1 V/m from the solution side (k̂ = −z), normalised to I₀."""
    grid = plane_wave_grid()
    i0 = N_WATER / (2.0 * Z0)  # I₀ = n |E₀|² / (2 Z₀)
    args: dict[str, object] = {
        "irradiance_w_m2": i0,
        "k_hat": (0.0, 0.0, -1.0),
        "polarisation": (1.0, 0.0, 0.0),
    }
    args.update(kwargs)
    return FieldMap(grid, 0.64, unit_field(grid), **args)  # type: ignore[arg-type]


def irradiance(e: NDArray[np.complex128], n: float) -> NDArray[np.float64]:
    return np.asarray(0.5 * n / Z0 * np.sum(np.abs(e) ** 2, axis=-1))


# --- T42: normalisation to the incident irradiance -----------------------------------------


def test_t42_irradiance_normalised_map_needs_conversion() -> None:
    m = hpfem_map()
    assert not m.is_port_normalised
    with pytest.raises(ValueError, match="per_port_power"):
        m.e_at([[0.0, 0.0, 0.5]], 1.0)


def test_t42_conversion_to_port_power() -> None:
    # 1 mW uniformly on a disc of 5 mm diameter: 1e-3 / (π · 2.5e-3²) W/m² per mW.
    per_mw = 1e-3 / (math.pi * 2.5e-3**2)
    m = hpfem_map().per_port_power(per_mw, "P")
    assert m.is_port_normalised
    e = m.e_at([[0.0, 0.0, 0.5]], 10.0)
    assert irradiance(e, N_WATER)[0] == pytest.approx(10.0 * per_mw, rel=1e-12)
    assert m.source["port_irradiance_status"] == "P"
    assert float(m.source["port_irradiance_w_m2_per_mw"]) == per_mw
    assert float(m.source["converted_from_irradiance_w_m2"]) == N_WATER / (2.0 * Z0)
    assert m.k_hat is not None and m.polarisation is not None  # metadata kept


def test_t42_thermal_fields_scale_linearly() -> None:
    grid = plane_wave_grid()
    q = np.full(grid.shape, 4.0)
    dt = np.full(grid.shape, 0.5)
    m = hpfem_map(irradiance_w_m2=200.0, q_w_m3=q, delta_t_k=dt)
    c = m.per_port_power(50.0, "P")
    assert c.q_w_m3 is not None and c.delta_t_k is not None
    np.testing.assert_allclose(c.q_w_m3, 1.0, rtol=1e-15)  # 4 · 50/200
    np.testing.assert_allclose(c.delta_t_k, 0.125, rtol=1e-15)


@pytest.mark.parametrize(
    ("args", "match"),
    [
        ((50.0, "L"), "reference"),
        ((50.0, "X"), "status"),
        ((-1.0, "P"), "positive"),
        ((math.inf, "P"), "positive"),
    ],
)
def test_t42_conversion_arguments_checked(args: tuple[float, str], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        hpfem_map().per_port_power(*args)


def test_t42_port_map_cannot_be_converted_twice() -> None:
    m = hpfem_map().per_port_power(50.0, "P")
    with pytest.raises(ValueError, match="already"):
        m.per_port_power(50.0, "P")


@pytest.mark.parametrize("value", [0.0, -1.0, math.nan])
def test_t42_invalid_incident_irradiance(value: float) -> None:
    with pytest.raises(ValueError, match="irradiance_w_m2"):
        hpfem_map(irradiance_w_m2=value)


def test_t42_nonlinear_map_must_be_port_normalised() -> None:
    with pytest.raises(ValueError, match="nonlinear"):
        hpfem_map(nonlinear=True, power_mw=3.0)


# --- T42: direction and polarisation -------------------------------------------------------


def test_t42_direction_and_polarisation_stored() -> None:
    m = hpfem_map(k_hat=(0.6, 0.0, -0.8), polarisation=(0.8, 0.0, 0.6))
    assert m.k_hat is not None and m.polarisation is not None
    np.testing.assert_array_equal(m.k_hat, [0.6, 0.0, -0.8])
    assert m.polarisation.dtype == np.complex128
    # circular polarisation ⟂ −z
    s = 1.0 / math.sqrt(2.0)
    hpfem_map(polarisation=(s, 1j * s, 0.0))


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"k_hat": (0.0, 0.0, -2.0)}, "k_hat"),
        ({"k_hat": (0.0, 0.0)}, "k_hat"),
        ({"polarisation": (2.0, 0.0, 0.0)}, "polarisation"),
        ({"polarisation": (0.0, 0.0, 1.0)}, "perpendicular"),
    ],
)
def test_t42_invalid_direction_or_polarisation(kwargs: dict[str, object], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        hpfem_map(**kwargs)


# --- T42: medium per grid point ------------------------------------------------------------


def layered_map() -> FieldMap:
    """Gold below z = 0, water from z = 0 (interface nodes hold the solution-side limit)."""
    grid = RegularGrid((np.array([0.0, 1.0]), np.array([0.0]), np.array([-0.1, 0.0, 0.1, 0.2])))
    medium = np.zeros(grid.shape, dtype=np.int64)
    medium[:, :, 1:] = 1
    return FieldMap(grid, 0.64, unit_field(grid), medium=medium, media=("gold", "water"))


def test_t42_evaluation_inside_one_medium() -> None:
    m = layered_map()
    e = m.e_at([[0.5, 0.0, 0.0], [0.3, 0.0, 0.05], [1.0, 0.0, 0.2]], 1.0)
    np.testing.assert_allclose(e[:, 0], 1.0)
    assert m.medium_at([[0.5, 0.0, 0.15], [0.5, 0.0, -0.1]]).tolist() == ["water", "gold"]


def test_t42_evaluation_across_media_rejected() -> None:
    m = layered_map()
    with pytest.raises(ValueError, match="media"):
        m.e_at([[0.5, 0.0, -0.02]], 1.0)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"media": ("gold",)}, "medium"),
        ({"medium": np.zeros((2, 1, 3), dtype=np.int64), "media": ("water",)}, "shape"),
        ({"medium": np.full((2, 1, 4), 2, dtype=np.int64), "media": ("a", "b")}, "index"),
        ({"medium": np.zeros((2, 1, 4)), "media": ("a",)}, "integer"),
        ({"medium": np.zeros((2, 1, 4), dtype=np.int64), "media": ("a", "a")}, "unique"),
    ],
)
def test_t42_invalid_medium(kwargs: dict[str, object], match: str) -> None:
    grid = RegularGrid((np.array([0.0, 1.0]), np.array([0.0]), np.array([-0.1, 0.0, 0.1, 0.2])))
    with pytest.raises(ValueError, match=match):
        FieldMap(grid, 0.64, unit_field(grid), **kwargs)  # type: ignore[arg-type]


# --- T42: nested maps ----------------------------------------------------------------------


def test_t42_stack_normalisation_consistent_and_converted() -> None:
    a = hpfem_map()
    b = FieldMap(
        RegularGrid((np.linspace(-1, 1, 3), np.linspace(-1, 1, 3), np.linspace(0, 2, 3))),
        0.64,
        np.ones((3, 3, 3, 3), dtype=np.complex128),
        irradiance_w_m2=a.irradiance_w_m2,
    )
    stack = FieldMapStack((a, b))
    with pytest.raises(ValueError, match="per_port_power"):
        stack.e_at([[0.0, 0.0, 0.5]], 1.0)
    conv = stack.per_port_power(50.0, "P")
    e = conv.e_at([[0.0, 0.0, 0.5], [0.5, 0.5, 1.5]], 1.0)
    assert irradiance(e[:1], N_WATER)[0] == pytest.approx(50.0, rel=1e-12)
    with pytest.raises(ValueError, match="normalisation"):
        FieldMapStack((a, a.per_port_power(50.0, "P")))


# --- T42: file format ----------------------------------------------------------------------


def test_t42_json_round_trip_with_metadata() -> None:
    grid = RegularGrid((np.array([0.0, 1.0]), np.array([0.0]), np.array([-0.1, 0.0, 0.1, 0.2])))
    medium = np.zeros(grid.shape, dtype=np.int64)
    medium[:, :, 1:] = 1
    s = 1.0 / math.sqrt(2.0)
    m = FieldMap(
        grid,
        0.64,
        unit_field(grid),
        irradiance_w_m2=1.5,
        k_hat=np.array([0.0, 0.0, -1.0]),
        polarisation=np.array([s, 1j * s, 0.0]),
        medium=medium,
        media=("gold", "water"),
        source={"phase_reference": "0,0,0"},
    )
    back = FieldMap.from_json(m.to_json())
    assert back.irradiance_w_m2 == 1.5
    assert back.k_hat is not None and back.polarisation is not None and back.medium is not None
    np.testing.assert_array_equal(back.k_hat, m.k_hat)
    np.testing.assert_array_equal(back.polarisation, m.polarisation)
    np.testing.assert_array_equal(back.medium, medium)
    assert back.media == ("gold", "water")
    assert back.source == m.source


def test_t42_old_files_without_metadata_load() -> None:
    grid = plane_wave_grid()
    old = {
        "grid": grid.to_data(),
        "wavelength_um": 0.64,
        "e": encode_array(unit_field(grid)),
        "nonlinear": False,
        "power_mw": None,
        "q_w_m3": None,
        "delta_t_k": None,
        "source": {"solver": "hp-FEM"},
    }
    m = FieldMap.from_json(to_json("FieldMap", old))
    assert m.is_port_normalised and m.k_hat is None and m.medium is None and m.media == ()


def test_t42_reader_takes_metadata_from_source() -> None:
    # Until the format has the keys, hp-FEM may write them as JSON strings into ``source``.
    grid = plane_wave_grid()
    src = {
        "solver": "hp-FEM",
        "irradiance_w_m2": "0.25",
        "k_hat": "[0.0, 0.0, -1.0]",
        "polarisation": "[[0.0, 1.0], [0.0, 0.0], [0.0, 0.0]]",
    }
    old = {
        "grid": grid.to_data(),
        "wavelength_um": 0.64,
        "e": encode_array(unit_field(grid)),
        "nonlinear": False,
        "power_mw": None,
        "q_w_m3": None,
        "delta_t_k": None,
        "source": src,
    }
    m = FieldMap.from_json(to_json("FieldMap", old))
    assert m.irradiance_w_m2 == 0.25
    assert m.k_hat is not None and m.polarisation is not None
    np.testing.assert_array_equal(m.k_hat, [0.0, 0.0, -1.0])
    np.testing.assert_array_equal(m.polarisation, [1j, 0.0, 0.0])
    data = load_json(m.to_json(), "FieldMap")
    assert data["irradiance_w_m2"] == 0.25  # written as proper keys from now on
    assert json.loads(data["source"]["k_hat"]) == [0.0, 0.0, -1.0]  # provenance kept


def test_t42_converted_map_stays_converted_after_reload() -> None:
    # hp-FEM wrote I₀ into ``source``; after the conversion a reload must not take it as I₀.
    m = hpfem_map(source={"irradiance_w_m2": "0.5"}).per_port_power(50.0, "P")
    assert "irradiance_w_m2" not in m.source
    back = FieldMap.from_json(m.to_json())
    assert back.is_port_normalised
    np.testing.assert_array_equal(back.e, m.e)


def test_t42_reader_takes_medium_from_source() -> None:
    grid = RegularGrid((np.array([0.0, 1.0]), np.array([0.0]), np.array([-0.1, 0.0, 0.1, 0.2])))
    medium = np.zeros(grid.shape, dtype=np.int64)
    medium[:, :, 1:] = 1
    old = {
        "grid": grid.to_data(),
        "wavelength_um": 0.64,
        "e": encode_array(unit_field(grid)),
        "nonlinear": False,
        "power_mw": None,
        "q_w_m3": None,
        "delta_t_k": None,
        "source": {"medium": json.dumps(medium.tolist()), "media": '["gold", "water"]'},
    }
    m = FieldMap.from_json(to_json("FieldMap", old))
    assert m.media == ("gold", "water")
    assert m.medium is not None
    np.testing.assert_array_equal(m.medium, medium)
