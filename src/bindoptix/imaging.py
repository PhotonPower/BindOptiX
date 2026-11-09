"""Synthetic camera frames of the field of view (stage Z7).

Expected photons per pixel from point sources (captured cells) and areal sources (rest light of
the spots, global background), imaged with a Gaussian point spread function: points are
integrated exactly over each pixel (erf), areal sources are rasterised at pixel centres and
blurred. A multiplicative fixed pattern (P) models the unevenness of the spots. Frames are in the
sensor frame (µm); pixel (i, j) covers x0 + [j, j+1)·p, y0 + [i, i+1)·p.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.ndimage import gaussian_filter
from scipy.special import erf

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class FrameGeometry:
    """Region of interest in the object plane."""

    x0_um: float
    y0_um: float
    nx: int
    ny: int
    pixel_um: float

    def __post_init__(self) -> None:
        if not (self.nx > 0 and self.ny > 0 and self.pixel_um > 0):
            raise ValueError("nx, ny and pixel_um must be positive")

    @property
    def area_mm2(self) -> float:
        return self.nx * self.ny * self.pixel_um**2 * 1e-6

    def centres_um(self) -> tuple[FloatArray, FloatArray]:
        x = self.x0_um + (np.arange(self.nx) + 0.5) * self.pixel_um
        y = self.y0_um + (np.arange(self.ny) + 0.5) * self.pixel_um
        return np.asarray(x), np.asarray(y)

    def to_pixel(self, xy_um: ArrayLike) -> FloatArray:
        """Sensor-frame points (N, 2) → fractional pixel coordinates (column, row)."""
        p = np.atleast_2d(np.asarray(xy_um, dtype=np.float64))
        return np.column_stack(
            [(p[:, 0] - self.x0_um) / self.pixel_um, (p[:, 1] - self.y0_um) / self.pixel_um]
        )

    def to_object(self, col_row: ArrayLike) -> FloatArray:
        c = np.atleast_2d(np.asarray(col_row, dtype=np.float64))
        return np.column_stack(
            [
                self.x0_um + (c[:, 0] + 0.5) * self.pixel_um,
                self.y0_um + (c[:, 1] + 0.5) * self.pixel_um,
            ]
        )


def _polygon_mask(geom: FrameGeometry, polygon_um: FloatArray) -> NDArray[np.bool_]:
    x, y = geom.centres_um()
    xx, yy = np.meshgrid(x, y)
    p = np.asarray(polygon_um, dtype=np.float64)
    inside = np.ones(xx.shape, dtype=bool)
    for i in range(p.shape[0]):
        a, b = p[i], p[(i + 1) % p.shape[0]]
        inside &= (b[0] - a[0]) * (yy - a[1]) - (b[1] - a[1]) * (xx - a[0]) >= 0.0
    return inside


def fixed_pattern(geom: FrameGeometry, rms: float, corr_um: float, seed: int) -> FloatArray:
    """Multiplicative unevenness, mean 1, relative rms ``rms``, correlation length ``corr_um``."""
    g = np.random.default_rng(seed)
    noise = g.normal(size=(geom.ny, geom.nx))
    s = max(corr_um / geom.pixel_um, 1e-9)
    sm = gaussian_filter(noise, s, mode="wrap")
    sm /= sm.std() if sm.std() > 0 else 1.0
    return np.asarray(np.clip(1.0 + rms * sm, 0.0, None))


def render(
    geom: FrameGeometry,
    psf_sigma_um: float,
    points_um: ArrayLike | None = None,
    point_photons: ArrayLike | None = None,
    areas: tuple[tuple[FloatArray, float], ...] = (),
    area_pattern: FloatArray | None = None,
    uniform_photons_per_um2: float = 0.0,
) -> FloatArray:
    """Expected photons per pixel (ny, nx).

    ``areas``: (convex counter-clockwise polygon, photons per µm²) pairs, multiplied by
    ``area_pattern`` (fixed pattern) where given.
    """
    img = np.zeros((geom.ny, geom.nx))
    pa = geom.pixel_um**2
    for poly, dens in areas:
        img += _polygon_mask(geom, poly) * dens * pa
    if area_pattern is not None:
        img *= area_pattern
    if areas:
        img = gaussian_filter(img, psf_sigma_um / geom.pixel_um, mode="constant")
    img += uniform_photons_per_um2 * pa
    pts = np.zeros((0, 2)) if points_um is None else np.asarray(points_um, float).reshape(-1, 2)
    ph = np.zeros(0) if point_photons is None else np.asarray(point_photons, float).reshape(-1)
    if pts.shape[0] != ph.size:
        raise ValueError("one photon number per point")
    xe = geom.x0_um + np.arange(geom.nx + 1) * geom.pixel_um
    ye = geom.y0_um + np.arange(geom.ny + 1) * geom.pixel_um
    r = 6.0 * psf_sigma_um
    s2 = math.sqrt(2.0) * psf_sigma_um
    for (px, py), n in zip(pts, ph, strict=True):
        i0 = max(int((px - r - geom.x0_um) // geom.pixel_um), 0)
        i1 = min(int((px + r - geom.x0_um) // geom.pixel_um) + 1, geom.nx)
        j0 = max(int((py - r - geom.y0_um) // geom.pixel_um), 0)
        j1 = min(int((py + r - geom.y0_um) // geom.pixel_um) + 1, geom.ny)
        if i0 >= i1 or j0 >= j1:
            continue
        fx = 0.5 * np.diff(erf((xe[i0 : i1 + 1] - px) / s2))
        fy = 0.5 * np.diff(erf((ye[j0 : j1 + 1] - py) / s2))
        img[j0:j1, i0:i1] += n * np.outer(fy, fx)
    return img
