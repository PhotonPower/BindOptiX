"""Independent reference: total power of an electric dipole above a planar stack (T4).

Test code only (AGENTS.md rule 4: no optics in the package). Two separately written forms of the
same classical result, normalised to the free dipole in the cover medium 1 (lossless, index
n1), dipole at height z above the top interface, layers (index, thickness) and a substrate:

* :func:`cps_rates`: Chance, Prock, Silbey (chance1978), Eqs. (2.17) and (2.29), integration
  variable u = k_∥/k₁, l_j = −i (ε_j/ε₁ − u²)^½, reflection coefficients in their convention
  R^∥ = (ε₁l₂ − ε₂l₁)/(ε₁l₂ + ε₂l₁) (2.15), R^⊥ = (l₁ − l₂)/(l₁ + l₂) (2.27), layers by
  (3.27b)/(3.28) applied recursively.
* :func:`lukosz_rates`: Lukosz and Kunz (lukosz1977), Eqs. (3.15) and (3.16) with the
  interchange (3.17) for the electric dipole, variable v = k_z,1/k₁ on the path 0 → 1 and
  i∞ → 0, Fresnel coefficients (2.34), (2.35) in the usual convention, k_z with Im ≥ 0 (2.27).

Returns (F_⊥, F_∥), the normalised total powers (decay rates for q = 1) of the perpendicular
and the parallel dipole.

Closed forms for the checks:

* :func:`mirror_rates`: perfect conductor from the field of the image dipole (perpendicular
  image parallel, parallel image antiparallel at distance 2z, dipole field as in
  jackson1999), F = 1 + Im(μ̂·E_image)/Im(μ̂·E_self) with Im(E_self) = k³p/(6πε).
* :func:`quasistatic_rates`: leading term for k₁z → 0 above an absorbing half-space from the
  electrostatic image β p, β = (ε₂ − ε₁)/(ε₂ + ε₁).
* :func:`lukosz_asymptotic`: Lukosz and Kunz (4.4)–(4.6) for a lossless half-space, k₁z ≫ 1.
"""

from __future__ import annotations

import cmath
import math
from collections.abc import Callable, Sequence
from itertools import pairwise

import numpy as np
from scipy.integrate import quad

Layer = tuple[complex, float]  # (refractive index, thickness in the units of the wavelength)


def _segments(lo: float, hi: float, width: float) -> list[tuple[float, float]]:
    n = max(1, math.ceil((hi - lo) / width))
    edges = [float(e) for e in np.linspace(lo, hi, n + 1)]
    return list(pairwise(edges))


def _integrate(f: Callable[[float], float], segments: Sequence[tuple[float, float]]) -> float:
    return float(sum(quad(f, a, b, limit=200, epsabs=1e-13, epsrel=1e-11)[0] for a, b in segments))


def _evanescent_extent(n1: float, layers: Sequence[Layer], substrate: complex, k1z: float) -> float:
    """Upper end of u (or s): past all modes and the decay of exp(−2 k₁z √(u²−1))."""
    n_max = max([abs(substrate)] + [abs(n) for n, _ in layers]) / n1
    return max(4.0 * n_max + 4.0, 40.0 / k1z)


# --- Chance, Prock, Silbey -------------------------------------------------------------------


def _cps_reflection(
    u: complex, n1: float, layers: Sequence[Layer], substrate: complex, k1: float
) -> tuple[complex, complex]:
    """(R^∥, R^⊥) seen from medium 1 in the CPS convention, layers from the top."""
    eps = [n1**2] + [complex(n) ** 2 for n, _ in layers] + [complex(substrate) ** 2]
    ell = [-1j * cmath.sqrt(e / eps[0] - u * u) for e in eps]
    ell[0] = -1j * cmath.sqrt(1.0 - u * u) if abs(u) <= 1.0 else cmath.sqrt(u * u - 1.0)
    for j in range(1, len(ell)):
        if ell[j].real < 0.0:
            ell[j] = -ell[j]
    r_par = (eps[-2] * ell[-1] - eps[-1] * ell[-2]) / (eps[-2] * ell[-1] + eps[-1] * ell[-2])
    r_perp = (ell[-2] - ell[-1]) / (ell[-2] + ell[-1])
    for j in range(len(layers), 0, -1):  # layer j between media j−1 (above) and j+1 (below)
        a_par = (eps[j - 1] * ell[j] - eps[j] * ell[j - 1]) / (
            eps[j - 1] * ell[j] + eps[j] * ell[j - 1]
        )
        a_perp = (ell[j - 1] - ell[j]) / (ell[j - 1] + ell[j])
        ph = cmath.exp(-2.0 * ell[j] * k1 * layers[j - 1][1])
        r_par = (a_par + r_par * ph) / (1.0 + a_par * r_par * ph)
        r_perp = (a_perp + r_perp * ph) / (1.0 + a_perp * r_perp * ph)
    return r_par, r_perp


def cps_rates(
    n1: float,
    layers: Sequence[Layer],
    substrate: complex,
    wavelength: float,
    z: float,
    reflection: Callable[[complex], tuple[complex, complex]] | None = None,
) -> tuple[float, float]:
    """(F_⊥, F_∥) after chance1978 (2.17), (2.29) with q = 1."""
    k1 = 2.0 * math.pi * n1 / wavelength
    dh = k1 * z

    def refl(u: float) -> tuple[complex, complex]:
        if reflection is not None:
            return reflection(u)
        return _cps_reflection(u, n1, layers, substrate, k1)

    # u = sin θ on [0, 1]: u³/l₁ du = i u³ dθ, u/l₁ du = i u dθ, e^{−2l₁d̂} = e^{2i d̂ cos θ}
    def perp_prop(th: float) -> float:
        u = math.sin(th)
        r_par, _ = refl(u)
        return (r_par * cmath.exp(2j * dh * math.cos(th)) * 1j * u**3).imag

    def par_prop(th: float) -> float:
        u = math.sin(th)
        r_par, r_perp = refl(u)
        val = ((1.0 - u * u) * r_par + r_perp) * cmath.exp(2j * dh * math.cos(th)) * 1j * u
        return val.imag

    # u = cosh t on [1, ∞): l₁ = sinh t, u³/l₁ du = u³ dt, u/l₁ du = u dt
    def perp_ev(t: float) -> float:
        u = math.cosh(t)
        r_par, _ = refl(u)
        return (r_par * u**3).imag * math.exp(-2.0 * dh * math.sinh(t))

    def par_ev(t: float) -> float:
        u = math.cosh(t)
        r_par, r_perp = refl(u)
        return (((1.0 - u * u) * r_par + r_perp) * u).imag * math.exp(-2.0 * dh * math.sinh(t))

    t_max = math.acosh(_evanescent_extent(n1, layers, substrate, dh))
    prop = _segments(0.0, 0.5 * math.pi, 0.05)
    ev = _segments(0.0, t_max, 0.01)
    perp = 1.0 - 1.5 * (_integrate(perp_prop, prop) + _integrate(perp_ev, ev))
    par = 1.0 + 0.75 * (_integrate(par_prop, prop) + _integrate(par_ev, ev))
    return perp, par


# --- Lukosz and Kunz -------------------------------------------------------------------------


def _fresnel_stack(
    kz: Sequence[complex], eps: Sequence[complex], thick: Sequence[float]
) -> tuple[complex, complex]:
    """(r_s, r_p) of the stack seen from medium 0, usual convention, Airy recursion."""
    r_s = (kz[-2] - kz[-1]) / (kz[-2] + kz[-1])
    r_p = (eps[-1] / eps[-2] * kz[-2] - kz[-1]) / (eps[-1] / eps[-2] * kz[-2] + kz[-1])
    for j in range(len(kz) - 2, 0, -1):
        s = (kz[j - 1] - kz[j]) / (kz[j - 1] + kz[j])
        p = (eps[j] / eps[j - 1] * kz[j - 1] - kz[j]) / (eps[j] / eps[j - 1] * kz[j - 1] + kz[j])
        ph = cmath.exp(2j * kz[j] * thick[j - 1])
        r_s = (s + r_s * ph) / (1.0 + s * r_s * ph)
        r_p = (p + r_p * ph) / (1.0 + p * r_p * ph)
    return r_s, r_p


def lukosz_rates(
    n1: float, layers: Sequence[Layer], substrate: complex, wavelength: float, z: float
) -> tuple[float, float]:
    """(F_⊥, F_∥) after lukosz1977 (3.15)–(3.17)."""
    k0 = 2.0 * math.pi / wavelength
    k1 = n1 * k0
    eps = [complex(n1**2)] + [complex(n) ** 2 for n, _ in layers] + [complex(substrate) ** 2]
    thick = [d for _, d in layers]

    def coeffs(v: complex) -> tuple[complex, complex]:
        kt2 = k1 * k1 * (1.0 - v * v)  # (k_x² + k_y²)
        kz = [k1 * v]
        for e in eps[1:]:
            w = cmath.sqrt(e * k0 * k0 - kt2)
            kz.append(w if w.imag >= 0.0 else -w)
        return _fresnel_stack(kz, eps, thick)

    def perp_kernel(v: complex) -> complex:
        _, r_p = coeffs(v)
        return r_p * (1.0 - v * v)

    def par_kernel(v: complex) -> complex:
        r_s, r_p = coeffs(v)
        return r_s - r_p

    def on_path(kernel: Callable[[complex], complex]) -> float:
        # ∫₀¹ f(v) e^{2ik₁zv} dv + ∫_{i∞}^0 f(v) e^{2ik₁zv} dv, the second with v = i s
        a = _integrate(
            lambda v: (kernel(v) * cmath.exp(2j * k1 * z * v)).real, _segments(0, 1, 0.05)
        )
        s_max = math.sqrt(_evanescent_extent(n1, layers, substrate, k1 * z) ** 2 - 1.0)
        b = _integrate(
            lambda s: (-1j * kernel(1j * s)).real * math.exp(-2.0 * k1 * z * s),
            _segments(0.0, s_max, 0.01),
        )
        return a + b

    perp = 1.0 + 1.5 * on_path(perp_kernel)
    par = 0.5 * (1.0 + perp) + 0.75 * on_path(par_kernel)
    return perp, par


# --- closed forms ----------------------------------------------------------------------------


def mirror_rates(k1z: float) -> tuple[float, float]:
    """Perfect conductor, image dipole at x = 2k₁z.

    Perpendicular: on-axis field of a parallel image, e^{ix}(2/x³ − 2i/x²) in units of
    p k³/(4πε); parallel: broadside field of an antiparallel image, −e^{ix}(1/x + i/x² − 1/x³).
    Divided by Im(E_self) = 2/3 in the same units.
    """
    x = 2.0 * k1z
    e = cmath.exp(1j * x)
    perp = 1.0 + 1.5 * (e * (2.0 / x**3 - 2j / x**2)).imag
    par = 1.0 - 1.5 * (e * (1.0 / x + 1j / x**2 - 1.0 / x**3)).imag
    return perp, par


def quasistatic_rates(eps1: float, eps2: complex, k1z: float) -> tuple[float, float]:
    """Leading term for k₁z → 0: image β p at distance 2z in the near field.

    Perpendicular: E = 2βp/(4πε₁D³) along the axis; parallel: E = βp/(4πε₁D³) (antiparallel image
    −βp, broadside field −p/(4πε₁D³)); divided by Im(E_self) = k³p/(6πε₁), D = 2z.
    """
    beta = (eps2 - eps1) / (eps2 + eps1)
    par = 3.0 * beta.imag / (16.0 * k1z**3)
    return 2.0 * par, par


def lukosz_asymptotic(n: float, k1z: float) -> tuple[float, float]:
    """lukosz1977 (4.4)–(4.6), electric dipole (upper signs), z̄ = 4πz/λ₁ = 2k₁z."""
    zb = 2.0 * k1z
    r = (n - 1.0) / (n + 1.0)
    perp = 1.0 + 3.0 * r * (math.sin(zb) / zb**3 - math.cos(zb) / zb**2)
    par = 0.5 * (1.0 + perp) - 3.0 * r * math.sin(zb) / zb
    return perp, par
