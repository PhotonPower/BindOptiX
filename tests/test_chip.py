"""T29: chip geometry against analytic values (docs/testfaelle.md)."""

import math
from pathlib import Path

import numpy as np
import nurbforge as nf
import pytest
from scipy.integrate import quad

from bindoptix.chip import ChipLayout


def test_default_layout_fits_and_counts() -> None:
    c = ChipLayout()
    assert c.spot_centers_um().shape == (25, 3)
    # farthest corner of the 5 × 5 array: √2 · (2 · 800 + 100) µm
    assert c.max_spot_radius_um == pytest.approx(math.sqrt(2) * 1700.0, rel=1e-15)
    assert c.max_spot_radius_um < 2500.0
    assert c.spot_area_total_um2 == pytest.approx(25 * 200.0**2, rel=1e-15)


def test_channel_floor_in_field_against_quadrature() -> None:
    # Independent reference: numerical integral of the chord length 2√(R² − y²).
    c = ChipLayout()
    r, w2 = c.field_radius_um, 0.5 * c.channel_width_um
    ref = sum(
        quad(lambda y: 2.0 * math.sqrt(r * r - y * y), y0 - w2, y0 + w2, epsabs=0, epsrel=1e-13)[0]
        for y0 in c.row_y_um()
    )
    assert c.channel_floor_in_field_um2() == pytest.approx(ref, rel=1e-11)


def test_channel_floor_full_circle_limit() -> None:
    # One channel as wide as the field covers the whole disc: πR².
    c = ChipLayout(
        field_diameter_um=200.0,
        n_rows=1,
        n_cols=1,
        spot_size_um=100.0,
        channel_width_um=200.0,
        channel_length_um=400.0,
    )
    assert c.channel_floor_in_field_um2() == pytest.approx(math.pi * 100.0**2, rel=1e-14)


def test_channel_shells_against_box_volume() -> None:
    c = ChipLayout()
    shells = c.channel_shells()
    assert len(shells) == 5
    for s in shells:
        assert s.is_closed
        assert s.is_consistently_oriented
        assert s.volume() == pytest.approx(10000.0 * 200.0 * 100.0, rel=1e-9)
    assert sum(s.volume() for s in shells) == pytest.approx(c.channel_volume_total_um3, rel=1e-9)


def test_spots_lie_in_their_channel() -> None:
    c = ChipLayout()
    shells = c.channel_shells()
    centres = c.spot_centers_um()
    for k in range(centres.shape[0]):
        row = k // c.n_cols
        probe = (float(centres[k, 0]), float(centres[k, 1]), 0.5 * c.channel_height_um)
        assert nf.point_in_shell(shells[row], probe)
        other = shells[(row + 1) % c.n_rows]
        assert not nf.point_in_shell(other, probe)


def test_spot_faces_area() -> None:
    c = ChipLayout()
    for f in c.spot_faces():
        assert nf.area(f.surface) == pytest.approx(200.0**2, rel=1e-12)


def test_field_face_area() -> None:
    c = ChipLayout()
    assert c.field_face().area() == pytest.approx(math.pi * 2500.0**2, rel=1e-6)


def test_step_round_trip(tmp_path: Path) -> None:
    c = ChipLayout()
    path = tmp_path / "chip.step"
    c.write_step(path)
    named = nf.read_step_named(str(path), unit="um")
    names = [str(item[0]) for item in named]
    assert names[:5] == [f"channel_{i}" for i in range(5)]
    vols = [item[1].volume() for item in named if str(item[0]).startswith("channel_")]
    np.testing.assert_allclose(vols, 10000.0 * 200.0 * 100.0, rtol=1e-9)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"channel_width_um": 250.0, "spot_size_um": 200.0}, "exceeds"),
        ({"spot_size_um": 220.0, "channel_width_um": 200.0}, "spot_size_um"),
        ({"pitch_x_um": 1000.0, "pitch_y_um": 1000.0}, "field of view"),
        ({"pitch_y_um": 150.0, "spot_size_um": 100.0, "channel_width_um": 150.0}, "touch"),
        ({"channel_length_um": 3000.0}, "shorter"),
        ({"channel_height_um": 0.0}, "positive"),
    ],
)
def test_invalid_layouts(kwargs: dict[str, float], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        ChipLayout(**kwargs)  # type: ignore[arg-type]
