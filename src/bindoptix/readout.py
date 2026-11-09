"""Digital readout: difference images, aperture detection, persistence, counts (stage Z7).

Each frame is compared with a reference frame taken before cells arrive; fixed structures
(uneven spots, static particles) cancel. Square apertures of the camera's footprint are summed
on the difference image; with a Gaussian approximation of the noise of the difference (twice the
variance of one frame) a threshold gives the target false-alarm density. Local maxima above the
threshold are candidates; a candidate is confirmed when it is found again at the same place
(within ``match_um``) in ``persistence`` consecutive frames, which rejects transient particles.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray
from scipy.ndimage import maximum_filter, uniform_filter
from scipy.stats import norm

from bindoptix.imaging import FrameGeometry
from bindoptix.signal import Camera

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class Detection:
    """Candidates of one frame: positions in the sensor frame and aperture sums [e⁻]."""

    xy_um: FloatArray
    value_e: FloatArray


def detect(
    frame_e: FloatArray,
    reference_e: FloatArray,
    expected_bg_e: FloatArray,
    geom: FrameGeometry,
    camera: Camera,
    false_alarm_per_mm2: float,
) -> Detection:
    """Candidates in the difference ``frame − reference``.

    ``expected_bg_e``: expected electrons per pixel without cells (model or a mean of
    reference frames), used for the noise level of the difference.
    """
    side = int(round(math.sqrt(camera.footprint_pixels)))
    diff = frame_e - reference_e
    agg = uniform_filter(diff, side, mode="constant") * side * side
    var = 2.0 * (
        uniform_filter(expected_bg_e, side, mode="constant") * side * side
        + side * side * camera.read_noise_e**2
    )
    p = false_alarm_per_mm2 * side * side * camera.pixel_um**2 * 1e-6
    z = float(norm.isf(p))
    thr = z * np.sqrt(np.maximum(var, 1e-12))
    peak = (agg == maximum_filter(agg, size=side)) & (agg > thr)
    rows, cols = np.nonzero(peak)
    xy = geom.to_object(np.column_stack([cols, rows]))
    return Detection(xy, agg[rows, cols])


@dataclass
class Track:
    xy_um: FloatArray
    first_frame: int
    last_frame: int
    hits: int = 1


@dataclass
class Tracker:
    """Persistence criterion over consecutive frames."""

    match_um: float
    persistence: int
    tracks: list[Track] = field(default_factory=list)
    confirmed: list[Track] = field(default_factory=list)

    def update(self, frame: int, det: Detection) -> None:
        used = np.zeros(det.xy_um.shape[0], dtype=bool)
        alive: list[Track] = []
        for tr in self.tracks:
            if det.xy_um.shape[0]:
                d = np.hypot(*(det.xy_um - tr.xy_um).T)
                d[used] = np.inf
                k = int(np.argmin(d))
                if d[k] <= self.match_um and tr.last_frame == frame - 1:
                    used[k] = True
                    tr.last_frame, tr.hits = frame, tr.hits + 1
                    if tr.hits == self.persistence:
                        self.confirmed.append(tr)
                    alive.append(tr)
                    continue
            if tr.hits >= self.persistence:
                continue  # confirmed tracks stay confirmed
        for k in np.flatnonzero(~used):
            tr = Track(det.xy_um[k].copy(), frame, frame)
            if self.persistence == 1:
                self.confirmed.append(tr)
            alive.append(tr)
        self.tracks = alive


def counts_per_spot(
    tracks: Sequence[Track],
    spot_polygons: Sequence[FloatArray],
    n_frames: int,
    margin_um: float = 0.0,
) -> NDArray[np.int64]:
    """Confirmed tracks per spot and frame (cumulative from their first frame), (frames, S).

    A track belongs to the first spot whose convex polygon, grown by ``margin_um``, contains it;
    the margin catches cells captured at the leading edge of a spot whose fitted position lies
    just outside (localisation error of the order of the point image).
    """
    out = np.zeros((n_frames, len(spot_polygons)), dtype=np.int64)
    for tr in tracks:
        for s, poly in enumerate(spot_polygons):
            p = np.asarray(poly, dtype=np.float64)
            a, b = p, np.roll(p, -1, axis=0)
            edge = b - a
            length = np.hypot(edge[:, 0], edge[:, 1])
            ok = length > 1e-12  # clipped polygons may repeat vertices
            rel = tr.xy_um - a[ok]
            signed = (edge[ok, 0] * rel[:, 1] - edge[ok, 1] * rel[:, 0]) / length[ok]
            if np.all(signed >= -margin_um):
                out[tr.first_frame :, s] += 1
                break
    return out
