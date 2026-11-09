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
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from bindoptix.binding import AvidityCapture, Occupancy
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


@dataclass
class CoupledRun:
    """Particles with fouling-dependent capture on one channel, and port B."""

    particles: RunResult
    occupancy: Occupancy
    free_sites_per_um2: FloatArray  # (nt, S) spot means
    times_s: FloatArray


def run_channel_coupled(
    layout: ChipLayout,
    k: int,
    cell: CellModel,
    flow_rate_um3_s: float,
    capture: AvidityCapture,
    sites_per_um2: float,
    n_cells: int,
    duration_s: float,
    dt_particles_s: float,
    dt_continuum_s: float,
    seed: int,
    fouling: Sequence[Species] = (),
    fouling_k_on_um3_s: Sequence[float] = (),
    fouling_k_off_s: Sequence[float] = (),
    fouling_c_in: Sequence[float] = (),
    strain: str = "target",
    affinity_per_spot: ArrayLike | None = None,
    mount: Mount | None = None,
    grid: dict[str, Any] | None = None,
    **kwargs: object,
) -> CoupledRun:
    """Fouling species (continuum, Z4) occupy sites; cells (particles, Z3) are captured with
    k_cap = k_bond · θ(t) · A_contact on each spot (one-way coupling: cells occupy a negligible
    share of the sites). ``affinity_per_spot`` (one factor per grid spot, default 1) scales
    k_bond for this strain and the aptamer of each spot: one row of the binding matrix (D13).
    Returns the particle run and port B (states: free, one per fouling species, "cell").
    """
    mnt = mount if mount is not None else Mount.sensor_floor()
    ch = layout.all_channels()[k]
    cs = channel_spots(layout, k)
    n_spots = len(cs.spots)
    n_steps = int(math.ceil(duration_s / dt_continuum_s - 1e-9))
    times = np.arange(n_steps + 1, dtype=np.float64) * dt_continuum_s
    n_f = len(fouling)
    shares = np.zeros((n_steps + 1, n_spots, n_f))
    theta = np.full((n_steps + 1, n_spots), sites_per_um2)
    if n_f:
        kon = [list(fouling_k_on_um3_s)] * (layout.n_rows * layout.n_cols)
        koff = [list(fouling_k_off_s)] * (layout.n_rows * layout.n_cols)
        vol, _ = channel_volume(
            layout, k, flow_rate_um3_s, fouling, sites_per_um2, kon, koff, **(grid or {})
        )
        st = vol.initial_state()
        areas = [vol.floor_area[vol.region_of == j] for j in range(n_spots)]
        for n in range(1, n_steps + 1):
            st = vol.step(st, dt_continuum_s, tuple(fouling_c_in))
            for j in range(n_spots):
                m = vol.region_of == j
                a = areas[j]
                for q, sp_q in enumerate(fouling):
                    g = float(np.sum(st.gamma[q][m] * a) / a.sum())
                    shares[n, j, q] = sp_q.site_weight * g / sites_per_um2
                theta[n, j] = sites_per_um2 * max(0.0, 1.0 - shares[n, j].sum())

    n_grid = layout.n_rows * layout.n_cols
    aff_grid = np.ones(n_grid) if affinity_per_spot is None else np.asarray(affinity_per_spot)
    if aff_grid.shape != (n_grid,) or np.any(aff_grid < 0.0):
        raise ValueError("affinity_per_spot needs one non-negative factor per grid spot")
    aff = aff_grid[list(cs.grid_index)]

    def rates_fn(t: float) -> FloatArray:
        th = np.array([np.interp(t, times, theta[:, j]) for j in range(n_spots)])
        return np.asarray(aff * capture.capture_rate_s(th), dtype=np.float64)

    model, _ = channel_model(
        layout, k, cell, flow_rate_um3_s, mnt, 0.0, spot_rates_fn=rates_fn, **kwargs
    )
    start = inlet_positions(model.flow, n_cells, seed + k)
    ids = np.arange(n_cells, dtype=np.int64) + k * 1_000_000_000
    r = run(model, start, duration_s, dt_particles_s, seed, ids=ids)

    got = r.state == CAPTURED
    t_cap = r.event_time_s[got]
    local = r.spot[got]
    xy = to_sensor_frame(ch, r.position_um[got, :2])
    spot_area = np.array([abs(_polygon_area(s.polygon_um)) for s in cs.spots], dtype=np.float64)
    fractions = np.zeros((n_steps + 1, n_spots, 2 + n_f))
    for n, t in enumerate(times):
        counts = np.bincount(local[t_cap <= t], minlength=n_spots).astype(np.float64)
        cell_share = np.minimum(counts * capture.contact_area_um2 / spot_area, 1.0)
        foul = shares[n] * (1.0 - cell_share)[:, None]
        fractions[n, :, 1 : 1 + n_f] = foul
        fractions[n, :, -1] = cell_share
        fractions[n, :, 0] = np.maximum(1.0 - foul.sum(axis=1) - cell_share, 0.0)
    occ = Occupancy(
        times_s=times,
        spots=cs.grid_index,
        states=("free", *(s.name for s in fouling), "cell"),
        fractions=fractions,
        cell_time_s=t_cap,
        cell_xy_um=xy,
        cell_spot=np.array([cs.grid_index[int(s)] for s in local], dtype=np.int64),
        cell_strain=tuple(strain for _ in range(int(got.sum()))),
    )
    return CoupledRun(r, occ, theta, times)
