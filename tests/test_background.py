"""Background sources (Z6, part 2): T41."""

import math

import numpy as np
import pytest

from bindoptix.background import Background, CellAutofluorescence, SurfaceBackground
from bindoptix.imaging import FrameGeometry, render
from bindoptix.ports import EmitterResponseTable, RegularGrid
from bindoptix.signal import placeholder_response_table, uniform_field_map
from bindoptix.units import SPEED_OF_LIGHT, VACUUM_PERMITTIVITY

SQUARE = np.array([[100.0, 100.0], [300.0, 100.0], [300.0, 300.0], [100.0, 300.0]])


def make_background() -> Background:
    return Background(
        surfaces=(
            SurfaceBackground("gold PL", "spots", 0.01, "P"),
            SurfaceBackground("2PP lacquer", "field", 0.004, "P"),
            SurfaceBackground("glass", "field", 0.002, "L", "ref2000"),
        ),
        cells=CellAutofluorescence(3e-3, "P"),
    )


# --- T41: surface sources ------------------------------------------------------------------


def test_t41_surface_densities_are_linear_in_power_and_time() -> None:
    bg = make_background()
    field, spots = bg.surface_photons_per_um2(power_mw=10.0, exposure_s=2.0)
    # by hand: (0.004 + 0.002) · 10 · 2 and 0.01 · 10 · 2
    assert field == pytest.approx(0.12, rel=1e-15)
    assert spots == pytest.approx(0.2, rel=1e-15)


def test_t41_render_terms_put_spot_sources_on_spots_only() -> None:
    bg = make_background()
    geom = FrameGeometry(0.0, 0.0, 200, 200, 2.0)
    areas, uniform = bg.render_terms((SQUARE,), power_mw=10.0, exposure_s=1.0)
    img = render(geom, 1.0, areas=areas, uniform_photons_per_um2=uniform)
    # totals: field 0.06/µm² over 400 µm × 400 µm, spot 0.1/µm² over 200 µm × 200 µm
    assert img.sum() == pytest.approx(0.06 * 400**2 + 0.1 * 200**2, rel=1e-6)
    assert img[10, 10] == pytest.approx(0.06 * 4.0, rel=1e-12)
    assert img[100, 100] == pytest.approx((0.06 + 0.1) * 4.0, rel=1e-9)


# --- T41: autofluorescence of cells --------------------------------------------------------


def test_t41_cell_autofluorescence_uniform_field_and_collection() -> None:
    # Independent: low-yield emitter, detected = k · I · η with A_det = η · 1.
    irr_per_mw, power, eta, k, t = 50.0, 10.0, 0.3, 3e-3, 2.0
    field = uniform_field_map(irr_per_mw, 0.64)
    table = placeholder_response_table(eta, 0.67)
    bg = make_background()
    pts = np.array([[0.0, 0.0, 0.5], [120.0, -40.0, 1.0]])
    ph = bg.cell_photons(pts, field, table, power, t)
    assert ph.shape == (2,)
    np.testing.assert_allclose(ph, k * irr_per_mw * power * eta * t, rtol=1e-12)


def test_t41_cell_autofluorescence_averages_anisotropic_collection() -> None:
    # Randomly oriented emitters see the trace: η = (a + b + c) / 3 for A_det = diag(a, b, c).
    grid = RegularGrid((np.zeros(1), np.zeros(1), np.zeros(1)))
    a_det = np.diag([0.1, 0.2, 0.6]).reshape(1, 1, 1, 1, 3, 3)
    eye = np.eye(3).reshape(1, 1, 1, 3, 3)
    table = EmitterResponseTable(grid, 0.67, 1.333, ("camera",), 2 * eye, eye, a_det)
    field = uniform_field_map(1.0, 0.64)
    bg = Background(cells=CellAutofluorescence(1.0, "P"))
    ph = bg.cell_photons([[0.0, 0.0, 0.5]], field, table, 1.0, 1.0)
    assert ph[0] == pytest.approx(0.3, rel=1e-12)


def test_t41_cell_irradiance_from_field_amplitude() -> None:
    # I = ½ n c ε₀ |E|² computed here from the amplitude of the uniform map.
    field = uniform_field_map(80.0, 0.64, n_host=1.333)
    amp2 = float(np.sum(np.abs(field.e_at([[0.0, 0.0, 0.0]], 1.0)) ** 2))
    irr = 0.5 * 1.333 * SPEED_OF_LIGHT * VACUUM_PERMITTIVITY * amp2
    assert irr == pytest.approx(80.0, rel=1e-12)
    bg = Background(cells=CellAutofluorescence(1.0, "P"))
    ph = bg.cell_photons([[0.0, 0.0, 0.5]], field, placeholder_response_table(1.0, 0.67), 1.0, 1.0)
    assert ph[0] == pytest.approx(irr, rel=1e-12)


def test_t41_no_cells_and_no_sources() -> None:
    bg = Background()
    assert bg.surface_photons_per_um2(10.0, 1.0) == (0.0, 0.0)
    assert bg.render_terms((SQUARE,), 10.0, 1.0) == ((), 0.0)
    field, table = uniform_field_map(1.0, 0.64), placeholder_response_table(0.3, 0.67)
    assert bg.cell_photons([[0.0, 0.0, 0.5]], field, table, 1.0, 1.0)[0] == 0.0
    assert bg.cell_photons(np.zeros((0, 3)), field, table, 1.0, 1.0).shape == (0,)


# --- T41: status of the parameters (AGENTS.md rule 9) --------------------------------------


def test_t41_placeholders_are_reported() -> None:
    bg = make_background()
    assert bg.placeholders == ("gold PL", "2PP lacquer", "cell autofluorescence")
    assert (
        Background(surfaces=(SurfaceBackground("glass", "field", 1.0, "L", "x"),)).placeholders
        == ()
    )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"region": "channel"},
        {"status": "X"},
        {"photons_per_um2_s_per_mw": -1.0},
        {"photons_per_um2_s_per_mw": math.nan},
        {"status": "L", "reference": ""},
    ],
)
def test_t41_invalid_surface_source_rejected(kwargs: dict[str, object]) -> None:
    args: dict[str, object] = {
        "name": "s",
        "region": "field",
        "photons_per_um2_s_per_mw": 1.0,
        "status": "K",
        "reference": "calibration 1",
    }
    args.update(kwargs)
    with pytest.raises(ValueError):
        SurfaceBackground(**args)  # type: ignore[arg-type]


def test_t41_invalid_cell_source_rejected() -> None:
    with pytest.raises(ValueError):
        CellAutofluorescence(-1.0, "P")
    with pytest.raises(ValueError):
        CellAutofluorescence(1.0, "L")
