"""Z3, part 2: particle runs on the chip layout (T33, docs/testfaelle.md)."""

import math

import numpy as np
import pytest

from bindoptix.cells import CellModel
from bindoptix.chip import Channel, ChipLayout
from bindoptix.chipsim import (
    channel_model,
    channel_spots,
    gravity_in_channel,
    run_chip,
    to_channel_frame,
    to_sensor_frame,
)
from bindoptix.flow import DuctFlow
from bindoptix.frames import Mount
from bindoptix.particles import CAPTURED, ChannelModel, inlet_positions, run

THETA = math.atan(1.0 / 5.0)
SEED = 20261112


def test_frame_round_trip() -> None:
    ch = Channel(offset_um=-350.0, angle_rad=0.7)
    rng = np.random.default_rng(SEED)
    p = rng.normal(scale=2000.0, size=(50, 3))
    np.testing.assert_allclose(to_sensor_frame(ch, to_channel_frame(ch, p)), p, atol=1e-9)


def test_spots_of_row_channels_in_channel_frame() -> None:
    # Independent: rows along the tilted grid → spots axis-aligned on the centre line at
    # x = L/2 + (j − 2) · 800 µm, each 200 × 200 µm.
    layout = ChipLayout(grid_angle_rad=THETA)
    for k in range(5):
        cs = channel_spots(layout, k)
        assert cs.grid_index == tuple(5 * k + j for j in range(5))
        for j, s in enumerate(cs.spots):
            c = 0.5 * (s.polygon_um.min(axis=0) + s.polygon_um.max(axis=0))
            np.testing.assert_allclose(c, [5000.0 + (j - 2) * 800.0, 0.0], atol=1e-9)
            span = s.polygon_um.max(axis=0) - s.polygon_um.min(axis=0)
            np.testing.assert_allclose(span, [200.0, 200.0], atol=1e-9)


def test_clipped_spot_area_matches_coverage() -> None:
    chans = (Channel(offset_um=-300.0), Channel(offset_um=150.0), Channel(offset_um=700.0))
    layout = ChipLayout(grid_angle_rad=THETA, channels=chans)
    for k in range(3):
        cs = channel_spots(layout, k)
        for s, f in zip(cs.spots, cs.fraction, strict=True):
            x, y = s.polygon_um[:, 0], s.polygon_um[:, 1]
            area = 0.5 * (x @ np.roll(y, -1) - y @ np.roll(x, -1))
            assert area == pytest.approx(f * 200.0**2, rel=1e-12)


def test_gravity_rotation() -> None:
    ch = Channel(offset_um=0.0, angle_rad=math.pi / 2)
    assert gravity_in_channel(ch, Mount.sensor_floor()) == pytest.approx((0.0, 0.0, -1.0))
    tilted = Mount(np.array([0.0, 1.0, -1.0]))  # gravity has a +y component in the sensor frame
    gx, gy, gz = gravity_in_channel(ch, tilted)
    # channel along +y: sensor +y is the flow direction
    assert (gx, gy, gz) == pytest.approx((math.sqrt(0.5), 0.0, -math.sqrt(0.5)))


def test_layout_run_equals_direct_run() -> None:
    layout = ChipLayout()
    cell = CellModel()
    q = 5e6 / 60.0
    model, _ = channel_model(layout, 2, cell, q, Mount.sensor_floor())
    direct = ChannelModel(DuctFlow(200.0, 100.0, q), cell, spots=model.spots, x_out_um=10000.0)
    start = inlet_positions(direct.flow, 200, SEED)
    a = run(model, start, 3600.0, 10.0, SEED)
    b = run(direct, start, 3600.0, 10.0, SEED)
    np.testing.assert_array_equal(a.position_um, b.position_um)
    np.testing.assert_array_equal(a.spot, b.spot)


def test_run_chip_counts_per_grid_spot() -> None:
    layout = ChipLayout(grid_angle_rad=THETA)
    res = run_chip(layout, CellModel(), 5e6 / 60.0, 300, 4 * 3600.0, 10.0, SEED)
    assert res.total_cells == 1500
    assert res.captured_per_spot.sum() == sum(r.count(CAPTURED) for r in res.runs)
    assert res.captured_per_spot.sum() > 0
    # captures of channel 3 land on the grid spots of row 3
    assert res.captured_per_spot[15:20].sum() == res.runs[3].count(CAPTURED)
