"""Chip geometry: field of view, flow channels and aptamer spots (stage Z1, docs/roadmap.md).

Sensor frame (docs/konventionen.md): µm, sensor surface and channel floor at z = 0, solution at
z > 0, centre of the field of view at the origin.

The spot array is a square n_rows × n_cols grid that may be tilted by ``grid_angle_rad`` against
the x axis (maintainer 2026-11-09: tilted 5 × 5 grid; final layout open, D15, D16). Channels are
straight strips with their own direction (flow along the channel, upstream first) and cross the
whole field of view. Without explicit channels, one channel runs along each grid row.

With tan θ = 1/n_cols every spot of the grid has its own position across a flow along x
(projections i·sin θ + j·cos θ · p are distinct, spacing p/√(1 + n²)); ``cross_flow_positions_um``
exposes this.

Values marked P are placeholders until the maintainer fixes them (AGENTS.md rule 9).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import nurbforge as nf
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]

MAX_CHANNEL_WIDTH_UM = 200.0
"""Upper limit of the channel width set by the maintainer (2026-11-09)."""


@dataclass(frozen=True)
class Channel:
    """Straight channel through the field of view.

    Attributes:
        offset_um: signed distance of the centre line from the field centre, measured along the
            left normal (−sin φ, cos φ) of the flow direction [µm].
        angle_rad: flow direction φ against +x.
        width_um: width, at most 200 µm.
        height_um: height above the floor z = 0 (P: 100 µm).
        length_um: length, centred on the foot point of the field centre (P: 10 mm).
    """

    offset_um: float
    angle_rad: float = 0.0
    width_um: float = 200.0
    height_um: float = 100.0
    length_um: float = 10000.0

    def __post_init__(self) -> None:
        for name in ("offset_um", "angle_rad"):
            if not math.isfinite(getattr(self, name)):
                raise ValueError(f"{name} must be finite")
        for name in ("width_um", "height_um", "length_um"):
            v = getattr(self, name)
            if not (math.isfinite(v) and v > 0.0):
                raise ValueError(f"{name} must be finite and positive")
        if self.width_um > MAX_CHANNEL_WIDTH_UM:
            raise ValueError(f"width_um exceeds {MAX_CHANNEL_WIDTH_UM} µm")

    @property
    def direction(self) -> FloatArray:
        return np.array([math.cos(self.angle_rad), math.sin(self.angle_rad)])

    @property
    def normal(self) -> FloatArray:
        return np.array([-math.sin(self.angle_rad), math.cos(self.angle_rad)])

    def floor_polygon_um(self) -> FloatArray:
        """Corners of the channel floor (4, 2), counter-clockwise [µm]."""
        c = self.offset_um * self.normal
        d = 0.5 * self.length_um * self.direction
        n = 0.5 * self.width_um * self.normal
        return np.array([c - d - n, c + d - n, c + d + n, c - d + n])


@dataclass(frozen=True)
class SpotCoverage:
    """A spot under a channel: grid index, covered area fraction, position along the flow."""

    spot: int
    fraction: float
    s_um: float


@dataclass(frozen=True)
class ChipLayout:
    """Tilted square spot grid in a circular field of view with straight channels.

    Attributes:
        field_diameter_um: diameter of the field of view of SIL and relay optics (D14: 5 mm).
        n_rows, n_cols: spot grid (5 × 5, maintainer 2026-11-09).
        spot_size_um: edge of a square spot, aligned with the grid (P: 200 µm).
        pitch_x_um, pitch_y_um: spot pitch along the grid axes (P: 800 µm).
        grid_angle_rad: tilt θ of the grid against x (P: 0; tan θ = 1/5 is a candidate).
        channel_width_um, channel_height_um, channel_length_um: used for the default channels
            (one per grid row) when ``channels`` is empty (P: 200 µm, 100 µm, 10 mm).
        channels: explicit channels; empty means one channel along each grid row.
    """

    field_diameter_um: float = 5000.0
    n_rows: int = 5
    n_cols: int = 5
    spot_size_um: float = 200.0
    pitch_x_um: float = 800.0
    pitch_y_um: float = 800.0
    grid_angle_rad: float = 0.0
    channel_width_um: float = 200.0
    channel_height_um: float = 100.0
    channel_length_um: float = 10000.0
    channels: tuple[Channel, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        for name in (
            "field_diameter_um",
            "spot_size_um",
            "pitch_x_um",
            "pitch_y_um",
            "channel_width_um",
            "channel_height_um",
            "channel_length_um",
        ):
            v = getattr(self, name)
            if not (math.isfinite(v) and v > 0.0):
                raise ValueError(f"{name} must be finite and positive")
        if not math.isfinite(self.grid_angle_rad):
            raise ValueError("grid_angle_rad must be finite")
        if self.n_rows < 1 or self.n_cols < 1:
            raise ValueError("n_rows and n_cols must be at least 1")
        if self.channel_width_um > MAX_CHANNEL_WIDTH_UM:
            raise ValueError(f"channel_width_um exceeds {MAX_CHANNEL_WIDTH_UM} µm")
        if self.spot_size_um > min(self.pitch_x_um, self.pitch_y_um):
            raise ValueError("spots overlap (spot_size_um larger than the pitch)")
        if self.max_spot_radius_um > self.field_radius_um:
            raise ValueError(
                f"spot array does not fit into the field of view: farthest spot corner at "
                f"{self.max_spot_radius_um:.1f} µm > {self.field_radius_um:.1f} µm"
            )
        if not self.channels:
            if self.spot_size_um > self.channel_width_um:
                raise ValueError("spot_size_um must not exceed channel_width_um")
            if self.n_rows > 1 and self.channel_width_um >= self.pitch_y_um:
                raise ValueError("neighbouring channels touch (channel_width_um >= pitch_y_um)")
        for k, ch in enumerate(self.all_channels()):
            if abs(ch.offset_um) >= self.field_radius_um:
                raise ValueError(f"channel {k} does not cross the field of view")
            half_chord = math.sqrt(self.field_radius_um**2 - ch.offset_um**2)
            if 0.5 * ch.length_um < half_chord + 0.5 * ch.width_um:
                raise ValueError(f"channel {k} is shorter than its chord through the field")
        self._check_channels_disjoint()

    # --- grid ------------------------------------------------------------------------------

    @property
    def field_radius_um(self) -> float:
        return 0.5 * self.field_diameter_um

    @property
    def _rot(self) -> FloatArray:
        c, s = math.cos(self.grid_angle_rad), math.sin(self.grid_angle_rad)
        return np.array([[c, -s], [s, c]])

    def _grid_local_um(self) -> FloatArray:
        i = np.arange(self.n_rows, dtype=np.float64) - 0.5 * (self.n_rows - 1)
        j = np.arange(self.n_cols, dtype=np.float64) - 0.5 * (self.n_cols - 1)
        u, v = np.meshgrid(j * self.pitch_x_um, i * self.pitch_y_um)
        return np.column_stack([u.ravel(), v.ravel()])

    def spot_centers_um(self) -> FloatArray:
        """Spot centres (n_rows · n_cols, 3) on z = 0, row by row in grid order [µm]."""
        xy = self._grid_local_um() @ self._rot.T
        return np.column_stack([xy, np.zeros(xy.shape[0])])

    def spot_polygon_um(self, k: int) -> FloatArray:
        """Corners of spot k (4, 2), counter-clockwise [µm]."""
        h = 0.5 * self.spot_size_um
        local = self._grid_local_um()[k] + np.array([[-h, -h], [h, -h], [h, h], [-h, h]])
        return np.asarray(local @ self._rot.T, dtype=np.float64)

    @property
    def max_spot_radius_um(self) -> float:
        """Largest distance of a spot corner from the field centre [µm] (rotation invariant)."""
        h = 0.5 * self.spot_size_um
        local = self._grid_local_um()
        return float(np.hypot(np.abs(local[:, 0]) + h, np.abs(local[:, 1]) + h).max())

    @property
    def spot_area_total_um2(self) -> float:
        return self.n_rows * self.n_cols * self.spot_size_um**2

    def cross_flow_positions_um(self, flow_angle_rad: float = 0.0) -> FloatArray:
        """Spot centres projected on the left normal of a flow direction [µm], grid order."""
        n = np.array([-math.sin(flow_angle_rad), math.cos(flow_angle_rad)])
        return np.asarray(self.spot_centers_um()[:, :2] @ n, dtype=np.float64)

    # --- channels --------------------------------------------------------------------------

    def all_channels(self) -> tuple[Channel, ...]:
        """Explicit channels, or one channel along each grid row."""
        if self.channels:
            return self.channels
        i = np.arange(self.n_rows, dtype=np.float64) - 0.5 * (self.n_rows - 1)
        return tuple(
            Channel(
                offset_um=float(y),
                angle_rad=self.grid_angle_rad,
                width_um=self.channel_width_um,
                height_um=self.channel_height_um,
                length_um=self.channel_length_um,
            )
            for y in i * self.pitch_y_um
        )

    def _check_channels_disjoint(self) -> None:
        chans = self.all_channels()
        for a in range(len(chans)):
            for b in range(a + 1, len(chans)):
                pa, pb = chans[a].floor_polygon_um(), chans[b].floor_polygon_um()
                if _polygon_area(_clip_convex(pa, pb)) > 1e-9 * chans[a].width_um ** 2:
                    raise ValueError(f"channels {a} and {b} overlap")

    def coverage(self) -> list[list[SpotCoverage]]:
        """Per channel: spots it covers (area fraction > 0), sorted upstream first."""
        out = []
        a_spot = self.spot_size_um**2
        centres = self.spot_centers_um()[:, :2]
        for ch in self.all_channels():
            floor = ch.floor_polygon_um()
            row = []
            for k in range(centres.shape[0]):
                frac = _polygon_area(_clip_convex(self.spot_polygon_um(k), floor)) / a_spot
                if frac > 1e-12:
                    row.append(SpotCoverage(k, frac, float(centres[k] @ ch.direction)))
            out.append(sorted(row, key=lambda c: c.s_um))
        return out

    def uncovered_spots(self) -> list[int]:
        """Spots not covered by any channel."""
        covered = {c.spot for row in self.coverage() for c in row}
        return [k for k in range(self.n_rows * self.n_cols) if k not in covered]

    def channel_floor_in_field_um2(self) -> float:
        """Channel floor inside the field of view, all channels [µm²] (exact).

        Area of the strip |y − y₀| ≤ w/2 inside the circle of radius R:
        F(y₂) − F(y₁) with F(y) = y√(R² − y²) + R² asin(y/R), y clipped to [−R, R];
        y₀ is the channel offset (the area does not depend on the direction).
        """
        r = self.field_radius_um

        def f(y: float) -> float:
            y = min(max(y, -r), r)
            return y * math.sqrt(max(r * r - y * y, 0.0)) + r * r * math.asin(y / r)

        return float(
            sum(
                f(ch.offset_um + 0.5 * ch.width_um) - f(ch.offset_um - 0.5 * ch.width_um)
                for ch in self.all_channels()
            )
        )

    def fluid_volume_in_field_um3(self) -> float:
        """Fluid volume above the channel floor inside the field of view [µm³]."""
        r = self.field_radius_um

        def f(y: float) -> float:
            y = min(max(y, -r), r)
            return y * math.sqrt(max(r * r - y * y, 0.0)) + r * r * math.asin(y / r)

        return float(
            sum(
                ch.height_um
                * (f(ch.offset_um + 0.5 * ch.width_um) - f(ch.offset_um - 0.5 * ch.width_um))
                for ch in self.all_channels()
            )
        )

    @property
    def channel_volume_total_um3(self) -> float:
        return float(sum(ch.length_um * ch.width_um * ch.height_um for ch in self.all_channels()))

    # --- B-rep via Nurbforge -----------------------------------------------------------------

    def channel_shells(self) -> list[Any]:
        """Closed Nurbforge shells of the fluid volume, one per channel (µm)."""
        shells = []
        for ch in self.all_channels():
            p = ch.floor_polygon_um()[0]
            ex = np.append(ch.length_um * ch.direction, 0.0)
            ey = np.append(ch.width_um * ch.normal, 0.0)
            shells.append(_box_shell(np.append(p, 0.0), ex, ey, ch.height_um))
        return shells

    def spot_faces(self) -> list[Any]:
        """Spots as planar faces on z = 0 (Nurbforge, µm), grid order."""
        faces = []
        s = self.spot_size_um
        ex = (s * math.cos(self.grid_angle_rad), s * math.sin(self.grid_angle_rad), 0.0)
        ey = (-s * math.sin(self.grid_angle_rad), s * math.cos(self.grid_angle_rad), 0.0)
        for k in range(self.n_rows * self.n_cols):
            p = self.spot_polygon_um(k)[0]
            faces.append(nf.TrimmedSurface(nf.make_rectangle((p[0], p[1], 0.0), ex, ey)))
        return faces

    def field_face(self) -> Any:
        """Field of view as a circular trimmed face on z = 0 (Nurbforge, µm)."""
        r = self.field_radius_um
        return nf.TrimmedSurface(
            nf.make_rectangle((-r, -r, 0.0), (2 * r, 0, 0), (0, 2 * r, 0)),
            nf.circle_loop(0.5, 0.5, 0.5),
        )

    def write_step(self, path: str | Path) -> None:
        """Channels (solids) and spots (surface models) as STEP AP214, unit µm."""
        channels = self.channel_shells()
        spots = [nf.Shell.sew([f]) for f in self.spot_faces()]
        names = [f"channel_{i}" for i in range(len(channels))]
        names += [f"spot_{i}_{j}" for i in range(self.n_rows) for j in range(self.n_cols)]
        nf.write_step(str(path), [*channels, *spots], unit="um", name="bindoptix-chip", names=names)


def _box_shell(origin: FloatArray, ex: FloatArray, ey: FloatArray, height: float) -> Any:
    """Closed, consistently oriented box shell (prism over the parallelogram origin, ex, ey)."""
    o = tuple(float(v) for v in origin)
    ez = (0.0, 0.0, height)
    exs, eys = tuple(float(v) for v in ex), tuple(float(v) for v in ey)
    top = (o[0], o[1], o[2] + height)
    far_x = tuple(o[i] + exs[i] for i in range(3))
    far_y = tuple(o[i] + eys[i] for i in range(3))
    faces = [
        nf.make_rectangle(o, exs, eys),
        nf.make_rectangle(top, exs, eys),
        nf.make_rectangle(o, exs, ez),
        nf.make_rectangle(far_y, exs, ez),
        nf.make_rectangle(o, eys, ez),
        nf.make_rectangle(far_x, eys, ez),
    ]
    return nf.Shell.sew([nf.TrimmedSurface(f) for f in faces]).oriented()


def _polygon_area(p: FloatArray) -> float:
    """Signed area of a polygon (n, 2), positive counter-clockwise (shoelace)."""
    if p.shape[0] < 3:
        return 0.0
    x, y = p[:, 0], p[:, 1]
    return float(0.5 * (x @ np.roll(y, -1) - y @ np.roll(x, -1)))


def _clip_convex(subject: FloatArray, clip: FloatArray) -> FloatArray:
    """Sutherland–Hodgman: subject polygon clipped by a convex counter-clockwise polygon."""
    out = [tuple(v) for v in subject]
    m = clip.shape[0]
    for e in range(m):
        a, b = clip[e], clip[(e + 1) % m]
        inp, out = out, []
        if not inp:
            break

        def side(p: tuple[float, ...], a: FloatArray = a, b: FloatArray = b) -> float:
            return float((b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]))

        prev = inp[-1]
        for cur in inp:
            sc, sp = side(cur), side(prev)
            if sc >= 0.0:
                if sp < 0.0:
                    t = sp / (sp - sc)
                    out.append(tuple(prev[i] + t * (cur[i] - prev[i]) for i in range(2)))
                out.append(cur)
            elif sp >= 0.0:
                t = sp / (sp - sc)
                out.append(tuple(prev[i] + t * (cur[i] - prev[i]) for i in range(2)))
            prev = cur
    return np.array(out, dtype=np.float64).reshape(-1, 2)
