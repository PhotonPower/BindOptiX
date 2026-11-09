"""Analytic pressure-driven flow in a rectangular channel (stage Z3, decision D7).

Channel frame: x along the flow, y across (−w/2 ≤ y ≤ w/2), z above the floor (0 ≤ z ≤ h), µm.
Fourier series of the Poiseuille problem in a rectangle (bruus2008):

    u(y, z) = (4h²G/π³) Σ_{n odd} n⁻³ [1 − cosh(nπy/h)/cosh(nπw/2h)] sin(nπz/h),
    Q = (h³wG/12) [1 − Σ_{n odd} 192 h/(π⁵ n⁵ w) · tanh(nπw/2h)],

with G = Δp/(ηL). The flow is given by Q; G follows from it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]


def _cosh_ratio(n_pi_over_h: FloatArray, y: FloatArray, half_w: float) -> FloatArray:
    """cosh(k y)/cosh(k w/2), overflow-free, broadcast over terms k and points y."""
    a = n_pi_over_h[:, None] * np.abs(y)[None, :]
    b = n_pi_over_h[:, None] * half_w
    return np.asarray(np.exp(a - b) * (1.0 + np.exp(-2.0 * a)) / (1.0 + np.exp(-2.0 * b)))


@dataclass(frozen=True)
class DuctFlow:
    """Poiseuille flow in a rectangular channel.

    Attributes:
        width_um, height_um: cross-section w × h.
        flow_rate_um3_s: volume flow Q [µm³/s] (1 nL/min = 1.667·10⁴ µm³/s).
        n_terms: number of odd Fourier terms.
    """

    width_um: float
    height_um: float
    flow_rate_um3_s: float
    n_terms: int = 60

    def __post_init__(self) -> None:
        for name in ("width_um", "height_um"):
            v = getattr(self, name)
            if not (math.isfinite(v) and v > 0.0):
                raise ValueError(f"{name} must be finite and positive")
        if not (math.isfinite(self.flow_rate_um3_s) and self.flow_rate_um3_s >= 0.0):
            raise ValueError("flow_rate_um3_s must be finite and non-negative")
        if self.n_terms < 1:
            raise ValueError("n_terms must be at least 1")

    @property
    def _n(self) -> FloatArray:
        return np.arange(1, 2 * self.n_terms, 2, dtype=np.float64)

    @property
    def shape_factor(self) -> float:
        """Q / (h³wG/12): 1 for wide channels, below 1 for finite aspect ratio."""
        n, h, w = self._n, self.height_um, self.width_um
        s = np.sum(192.0 * h / (math.pi**5 * n**5 * w) * np.tanh(n * math.pi * w / (2.0 * h)))
        return float(1.0 - s)

    @property
    def pressure_gradient(self) -> float:
        """G = Δp/(ηL) [1/(µm·s)] belonging to Q."""
        h, w = self.height_um, self.width_um
        return 12.0 * self.flow_rate_um3_s / (h**3 * w * self.shape_factor)

    @property
    def mean_speed_um_s(self) -> float:
        return self.flow_rate_um3_s / (self.width_um * self.height_um)

    def velocity_um_s(self, y_um: ArrayLike, z_um: ArrayLike) -> FloatArray:
        """u(y, z) [µm/s]; points outside the cross-section give 0."""
        y = np.asarray(y_um, dtype=np.float64).ravel()
        z = np.asarray(z_um, dtype=np.float64).ravel()
        h, half_w = self.height_um, 0.5 * self.width_um
        k = self._n * math.pi / h
        ratio = _cosh_ratio(k, np.clip(y, -half_w, half_w), half_w)
        terms = (1.0 / self._n**3)[:, None] * (1.0 - ratio) * np.sin(k[:, None] * z[None, :])
        u = 4.0 * h * h * self.pressure_gradient / math.pi**3 * terms.sum(axis=0)
        inside = (np.abs(y) <= half_w) & (z >= 0.0) & (z <= h)
        return np.asarray(np.where(inside, u, 0.0), dtype=np.float64)

    def floor_shear_rate_s(self, y_um: ArrayLike) -> FloatArray:
        """∂u/∂z at z = 0 [1/s].

        The part without the cosh ratio is summed in closed form, Σ_{n odd} n⁻² = π²/8, which
        gives G h/2; only the exponentially decaying wall terms are summed numerically.
        """
        y = np.asarray(y_um, dtype=np.float64).ravel()
        h, half_w = self.height_um, 0.5 * self.width_um
        g_p = self.pressure_gradient
        k = self._n * math.pi / h
        ratio = _cosh_ratio(k, np.clip(y, -half_w, half_w), half_w)
        wall = (1.0 / self._n**2)[:, None] * ratio
        g = 0.5 * g_p * h - 4.0 * h * g_p / math.pi**2 * wall.sum(axis=0)
        return np.asarray(np.where(np.abs(y) <= half_w, g, 0.0), dtype=np.float64)

    def column_flux_um2_s(self, y_um: ArrayLike) -> FloatArray:
        """q(y) = ∫₀ʰ u dz [µm²/s], exact termwise (∫ sin(nπz/h) dz = 2h/(nπ) for odd n)."""
        y = np.asarray(y_um, dtype=np.float64).ravel()
        h, half_w = self.height_um, 0.5 * self.width_um
        k = self._n * math.pi / h
        ratio = _cosh_ratio(k, np.clip(y, -half_w, half_w), half_w)
        terms = (1.0 / self._n**3 * 2.0 * h / (self._n * math.pi))[:, None] * (1.0 - ratio)
        q = 4.0 * h * h * self.pressure_gradient / math.pi**3 * terms.sum(axis=0)
        return np.asarray(np.where(np.abs(y) <= half_w, q, 0.0), dtype=np.float64)
