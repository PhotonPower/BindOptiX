"""Chip geometry: field of view, flow channels and aptamer spots (stage Z1, docs/roadmap.md).

Sensor frame (docs/konventionen.md): µm, sensor surface and channel floor at z = 0, solution at
z > 0, centre of the field of view at the origin. Layout ``rows`` (D15 open): one straight
channel per spot row, flow along +x, the spots of a row lie in series on the channel floor.

Values marked P in the docstrings are placeholders until the maintainer fixes them (AGENTS.md
rule 9). The analytic quantities are exact; the B-rep bodies come from Nurbforge.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import nurbforge as nf
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]

MAX_CHANNEL_WIDTH_UM = 200.0
"""Upper limit of the channel width set by the maintainer (2026-11-09)."""


@dataclass(frozen=True)
class ChipLayout:
    """Spot array in a circular field of view with one channel per spot row.

    Attributes:
        field_diameter_um: diameter of the field of view of SIL and relay optics (D14: 5 mm).
        n_rows, n_cols: spot array (5 × 5, maintainer 2026-11-09).
        spot_size_um: edge of a square spot (P: 200 µm).
        pitch_x_um, pitch_y_um: spot pitch along and across the flow (P: 800 µm).
        channel_width_um: channel width, at most 200 µm (P: 200 µm).
        channel_height_um: channel height (P: 100 µm).
        channel_length_um: channel length, centred on the field (P: 10 mm).
    """

    field_diameter_um: float = 5000.0
    n_rows: int = 5
    n_cols: int = 5
    spot_size_um: float = 200.0
    pitch_x_um: float = 800.0
    pitch_y_um: float = 800.0
    channel_width_um: float = 200.0
    channel_height_um: float = 100.0
    channel_length_um: float = 10000.0

    def __post_init__(self) -> None:
        positive = {
            "field_diameter_um": self.field_diameter_um,
            "spot_size_um": self.spot_size_um,
            "pitch_x_um": self.pitch_x_um,
            "pitch_y_um": self.pitch_y_um,
            "channel_width_um": self.channel_width_um,
            "channel_height_um": self.channel_height_um,
            "channel_length_um": self.channel_length_um,
        }
        for name, value in positive.items():
            if not (math.isfinite(value) and value > 0.0):
                raise ValueError(f"{name} must be finite and positive")
        if self.n_rows < 1 or self.n_cols < 1:
            raise ValueError("n_rows and n_cols must be at least 1")
        if self.channel_width_um > MAX_CHANNEL_WIDTH_UM:
            raise ValueError(f"channel_width_um exceeds {MAX_CHANNEL_WIDTH_UM} µm")
        if self.spot_size_um > self.channel_width_um:
            raise ValueError("spot_size_um must not exceed channel_width_um")
        if self.spot_size_um > self.pitch_x_um:
            raise ValueError("spots overlap along the channel (spot_size_um > pitch_x_um)")
        if self.n_rows > 1 and self.channel_width_um >= self.pitch_y_um:
            raise ValueError("neighbouring channels touch (channel_width_um >= pitch_y_um)")
        if self.max_spot_radius_um > self.field_radius_um:
            raise ValueError(
                f"spot array does not fit into the field of view: farthest spot corner at "
                f"{self.max_spot_radius_um:.1f} µm > {self.field_radius_um:.1f} µm"
            )
        half_array_x = 0.5 * (self.n_cols - 1) * self.pitch_x_um + 0.5 * self.spot_size_um
        if half_array_x > 0.5 * self.channel_length_um:
            raise ValueError("channels are shorter than the spot rows")

    # --- analytic geometry ----------------------------------------------------------------

    @property
    def field_radius_um(self) -> float:
        return 0.5 * self.field_diameter_um

    def row_y_um(self) -> FloatArray:
        """Centre lines of the channels (= spot rows) along y [µm]."""
        i = np.arange(self.n_rows, dtype=np.float64)
        return (i - 0.5 * (self.n_rows - 1)) * self.pitch_y_um

    def col_x_um(self) -> FloatArray:
        """Spot centres along the flow (x) [µm]."""
        j = np.arange(self.n_cols, dtype=np.float64)
        return (j - 0.5 * (self.n_cols - 1)) * self.pitch_x_um

    def spot_centers_um(self) -> FloatArray:
        """Spot centres (n_rows · n_cols, 3) on z = 0, row by row, upstream first [µm]."""
        x, y = np.meshgrid(self.col_x_um(), self.row_y_um())
        return np.column_stack([x.ravel(), y.ravel(), np.zeros(x.size)])

    @property
    def max_spot_radius_um(self) -> float:
        """Largest distance of a spot corner from the field centre [µm]."""
        h = 0.5 * self.spot_size_um
        cx = np.abs(self.col_x_um()).max() + h
        cy = np.abs(self.row_y_um()).max() + h
        return float(math.hypot(cx, cy))

    @property
    def spot_area_total_um2(self) -> float:
        return self.n_rows * self.n_cols * self.spot_size_um**2

    @property
    def channel_volume_total_um3(self) -> float:
        return self.n_rows * self.channel_length_um * self.channel_width_um * self.channel_height_um

    def channel_floor_in_field_um2(self) -> float:
        """Channel floor inside the field of view, all channels [µm²] (exact).

        Area of the strip |y − y₀| ≤ w/2 inside the circle of radius R:
        F(y₂) − F(y₁) with F(y) = y√(R² − y²) + R² asin(y/R), y clipped to [−R, R].
        """
        r = self.field_radius_um

        def f(y: float) -> float:
            y = min(max(y, -r), r)
            return y * math.sqrt(max(r * r - y * y, 0.0)) + r * r * math.asin(y / r)

        w2 = 0.5 * self.channel_width_um
        return float(sum(f(y0 + w2) - f(y0 - w2) for y0 in self.row_y_um()))

    def fluid_volume_in_field_um3(self) -> float:
        """Fluid volume above the channel floor inside the field of view [µm³]."""
        return self.channel_floor_in_field_um2() * self.channel_height_um

    # --- B-rep via Nurbforge --------------------------------------------------------------

    def channel_shells(self) -> list[Any]:
        """Closed Nurbforge shells of the fluid volume, one per channel (µm)."""
        return [
            _box_shell(
                (-0.5 * self.channel_length_um, y0 - 0.5 * self.channel_width_um, 0.0),
                (self.channel_length_um, self.channel_width_um, self.channel_height_um),
            )
            for y0 in self.row_y_um()
        ]

    def spot_faces(self) -> list[Any]:
        """Spots as trimmed planar faces on z = 0 (Nurbforge, µm), order of spot_centers_um."""
        s = self.spot_size_um
        return [
            nf.TrimmedSurface(
                nf.make_rectangle((c[0] - s / 2, c[1] - s / 2, 0.0), (s, 0, 0), (0, s, 0))
            )
            for c in self.spot_centers_um()
        ]

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


def _box_shell(origin: tuple[float, float, float], size: tuple[float, float, float]) -> Any:
    """Closed, consistently oriented box shell from six planar faces."""
    x0, y0, z0 = origin
    lx, ly, lz = size
    ex, ey, ez = (lx, 0.0, 0.0), (0.0, ly, 0.0), (0.0, 0.0, lz)
    faces = [
        nf.make_rectangle((x0, y0, z0), ex, ey),
        nf.make_rectangle((x0, y0, z0 + lz), ex, ey),
        nf.make_rectangle((x0, y0, z0), ex, ez),
        nf.make_rectangle((x0, y0 + ly, z0), ex, ez),
        nf.make_rectangle((x0, y0, z0), ey, ez),
        nf.make_rectangle((x0 + lx, y0, z0), ey, ez),
    ]
    return nf.Shell.sew([nf.TrimmedSurface(f) for f in faces]).oriented()
