"""Hydrodynamic hindrance of a sphere near a plane wall (stage Z3, part 3).

Gap h between sphere surface and wall, radius a. Normal mobility: rational approximation of
Brenner's exact series (brenner1961) by Bevan and Prieve (bevan2000),

    λ⊥(h) = (6h² + 2ah) / (6h² + 9ah + 2a²),

λ⊥ → 1 far from the wall and λ⊥ ≈ h/a for h → 0 (lubrication). Parallel mobility: Faxén's
expansion in r = a/(a + h) (happel1983),

    λ∥ = 1 − 9/16 r + 1/8 r³ − 45/256 r⁴ − 1/16 r⁵,

accurate for h ≳ 0.1a; at contact it gives 0.32 instead of the logarithmic decay to 0 and is used
there only as a placeholder (P). Both multiply the free diffusion coefficient and the mobility
(sedimentation). For position-dependent diffusion the Itô drift D₀ dλ⊥/dh is added so that the
equilibrium stays the Boltzmann distribution (T34).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]


def lambda_perp(gap_um: ArrayLike, radius_um: float) -> FloatArray:
    """Normal mobility factor λ⊥ ∈ [0, 1) (bevan2000)."""
    h = np.maximum(np.asarray(gap_um, dtype=np.float64), 0.0)
    a = radius_um
    return np.asarray((6.0 * h * h + 2.0 * a * h) / (6.0 * h * h + 9.0 * a * h + 2.0 * a * a))


def dlambda_perp(gap_um: ArrayLike, radius_um: float) -> FloatArray:
    """dλ⊥/dh [1/µm]."""
    h = np.maximum(np.asarray(gap_um, dtype=np.float64), 0.0)
    a = radius_um
    num = 6.0 * h * h + 2.0 * a * h
    den = 6.0 * h * h + 9.0 * a * h + 2.0 * a * a
    return np.asarray(((12.0 * h + 2.0 * a) * den - num * (12.0 * h + 9.0 * a)) / (den * den))


def lambda_par(gap_um: ArrayLike, radius_um: float) -> FloatArray:
    """Parallel mobility factor λ∥ (Faxén), clipped to [0, 1]."""
    h = np.maximum(np.asarray(gap_um, dtype=np.float64), 0.0)
    r = radius_um / (radius_um + h)
    lam = 1.0 - 9.0 / 16.0 * r + r**3 / 8.0 - 45.0 / 256.0 * r**4 - r**5 / 16.0
    return np.asarray(np.clip(lam, 0.0, 1.0))
