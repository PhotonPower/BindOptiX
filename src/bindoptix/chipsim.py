"""Particle runs on a chip layout (stage Z3, part 2).

Connects :class:`bindoptix.chip.ChipLayout` with :func:`bindoptix.particles.run`. Each channel gets
its own channel frame: x along the flow from the upstream end (x = 0) to the downstream end
(x = L), y along the left normal of the flow with y = 0 on the centre line, z above the floor.
For a sensor-frame point p (µm) and a channel with direction d, normal n, offset o and length L:

    x = p·d + L/2,   y = p·n − o,   z = p_z.

The parts of the spots that lie inside a channel become :class:`bindoptix.particles.Spot` in that
frame (convex clipping, orientation kept); their index is the grid index of the layout. Gravity
of a :class:`bindoptix.frames.Mount` is rotated into each channel frame.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from bindoptix.cells import CellModel
from bindoptix.chip import Channel, ChipLayout, _clip_convex, _polygon_area
from bindoptix.continuum import Species, graded_edges
from bindoptix.continuum3d import ChannelVolume, FloorRegion
from bindoptix.flow import DuctFlow
from bindoptix.frames import Mount
from bindoptix.particles import CAPTURED, ChannelModel, RunResult, Spot, inlet_positions, run

FloatArray = NDArray[np.float64]


def to_channel_frame(ch: Channel, p_um: ArrayLike) -> FloatArray:
    """Sensor-frame points (N, 2 or 3) → channel frame of ``ch`` [µm]."""
    p = np.atleast_2d(np.asarray(p_um, dtype=np.float64))
    x = p[:, :2] @ ch.direction + 0.5 * ch.length_um
    y = p[:, :2] @ ch.normal - ch.offset_um
    cols = [x, y] if p.shape[1] == 2 else [x, y, p[:, 2]]
    return np.column_stack(cols)


def to_sensor_frame(ch: Channel, q_um: ArrayLike) -> FloatArray:
    """Channel-frame points (N, 2 or 3) of ``ch`` → sensor frame [µm]."""
    q = np.atleast_2d(np.asarray(q_um, dtype=np.float64))
    s = q[:, 0] - 0.5 * ch.length_um
    t = q[:, 1] + ch.offset_um
    xy = s[:, None] * ch.direction + t[:, None] * ch.normal
    return np.column_stack([xy, q[:, 2]]) if q.shape[1] == 3 else xy


def gravity_in_channel(ch: Channel, mount: Mount) -> tuple[float, float, float]:
    """Unit gravity vector of the mount in the channel frame."""
    g = mount.gravity_direction
    return (float(g[:2] @ ch.direction), float(g[:2] @ ch.normal), float(g[2]))


@dataclass(frozen=True)
class ChannelSpots:
    """Spots of one channel in its frame, with their grid indices and covered area fractions."""

    spots: tuple[Spot, ...]
    grid_index: tuple[int, ...]
    fraction: tuple[float, ...]


def channel_spots(
    layout: ChipLayout, k: int, capture_rate_s: float | Sequence[float] = math.inf
) -> ChannelSpots:
    """Spots under channel k, clipped to the channel floor, in the channel frame.

    ``capture_rate_s`` is one rate for all spots or one per grid spot (length n_rows · n_cols).
    """
    ch = layout.all_channels()[k]
    n_spots = layout.n_rows * layout.n_cols
    rates = np.broadcast_to(np.asarray(capture_rate_s, dtype=np.float64), (n_spots,))
    floor = ch.floor_polygon_um()
    out, idx, frac = [], [], []
    for cov in layout.coverage()[k]:
        part = _clip_convex(layout.spot_polygon_um(cov.spot), floor)
        if part.shape[0] < 3 or _polygon_area(part) <= 0.0:
            continue
        out.append(Spot(to_channel_frame(ch, part), float(rates[cov.spot])))
        idx.append(cov.spot)
        frac.append(cov.fraction)
    return ChannelSpots(tuple(out), tuple(idx), tuple(frac))


def channel_model(
    layout: ChipLayout,
    k: int,
    cell: CellModel,
    flow_rate_um3_s: float,
    mount: Mount,
    capture_rate_s: float | Sequence[float] = math.inf,
    **kwargs: object,
) -> tuple[ChannelModel, ChannelSpots]:
    """ChannelModel for channel k of the layout; further keywords go to ChannelModel."""
    ch = layout.all_channels()[k]
    cs = channel_spots(layout, k, capture_rate_s)
    model = ChannelModel(
        flow=DuctFlow(ch.width_um, ch.height_um, flow_rate_um3_s),
        cell=cell,
        gravity_dir=gravity_in_channel(ch, mount),
        spots=cs.spots,
        x_out_um=ch.length_um,
        **kwargs,  # type: ignore[arg-type]
    )
    return model, cs


@dataclass
class ChipRun:
    """Result of a run over all channels: captures per grid spot and the raw runs."""

    captured_per_spot: NDArray[np.int64]
    cells_per_channel: NDArray[np.int64]
    runs: list[RunResult]
    spot_maps: list[ChannelSpots]

    @property
    def total_cells(self) -> int:
        return int(self.cells_per_channel.sum())


def run_chip(
    layout: ChipLayout,
    cell: CellModel,
    flow_rate_um3_s: float | Sequence[float],
    cells_per_channel: int | Sequence[int],
    duration_s: float,
    dt_s: float,
    seed: int,
    mount: Mount | None = None,
    capture_rate_s: float | Sequence[float] = math.inf,
    **kwargs: object,
) -> ChipRun:
    """Cells start flux-weighted at the inlet of every channel at t = 0.

    Particle ids are unique over the chip (channel k uses ids k·10⁹ + i); start positions use the
    seed ``seed + k``.
    """
    mnt = mount if mount is not None else Mount.sensor_floor()
    chans = layout.all_channels()
    nc = len(chans)
    q = np.broadcast_to(np.asarray(flow_rate_um3_s, dtype=np.float64), (nc,))
    n = np.broadcast_to(np.asarray(cells_per_channel, dtype=np.int64), (nc,))
    per_spot = np.zeros(layout.n_rows * layout.n_cols, dtype=np.int64)
    runs, maps = [], []
    for k in range(nc):
        model, cs = channel_model(layout, k, cell, float(q[k]), mnt, capture_rate_s, **kwargs)
        start = inlet_positions(model.flow, int(n[k]), seed + k, x_um=0.0)
        ids = np.arange(int(n[k]), dtype=np.int64) + k * 1_000_000_000
        r = run(model, start, duration_s, dt_s, seed, ids=ids)
        for local in r.spot[r.state == CAPTURED]:
            per_spot[cs.grid_index[int(local)]] += 1
        runs.append(r)
        maps.append(cs)
    return ChipRun(per_spot, np.array(n, dtype=np.int64), runs, maps)


def channel_volume(
    layout: ChipLayout,
    k: int,
    flow_rate_um3_s: float,
    species: Sequence[Species],
    sites_per_um2: float,
    k_on_um3_s: Sequence[Sequence[float]],
    k_off_s: Sequence[Sequence[float]],
    dx_um: float = 50.0,
    ny: int = 12,
    nz: int = 16,
    first_dz_um: float = 0.5,
) -> tuple[ChannelVolume, ChannelSpots]:
    """3D continuum model of channel k with its spots as reactive floor regions (Z4, part 2).

    ``k_on_um3_s[s][q]`` and ``k_off_s[s][q]`` per grid spot s and species q. The x grid is
    uniform with spacing ``dx_um`` plus the x extents of the spots; y uniform; z graded from the
    floor with first cell ``first_dz_um``.
    """
    ch = layout.all_channels()[k]
    cs = channel_spots(layout, k)
    flow = DuctFlow(ch.width_um, ch.height_um, flow_rate_um3_s)
    xs = [np.linspace(0.0, ch.length_um, int(math.ceil(ch.length_um / dx_um)) + 1)]
    regions = []
    for s, spot in zip(cs.grid_index, cs.spots, strict=True):
        xs.append(np.array([spot.polygon_um[:, 0].min(), spot.polygon_um[:, 0].max()]))
        regions.append(
            FloorRegion(spot.polygon_um, sites_per_um2, tuple(k_on_um3_s[s]), tuple(k_off_s[s]))
        )
    xe = np.unique(np.round(np.concatenate(xs), 9))
    ye = np.linspace(-0.5 * ch.width_um, 0.5 * ch.width_um, ny + 1)
    ze = graded_edges(ch.height_um, nz, first_dz_um)
    vol = ChannelVolume(xe, ye, ze, flow.velocity_um_s, species, regions)
    return vol, cs
