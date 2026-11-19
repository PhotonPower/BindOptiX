"""Independent reference: dipole above a planar stack (T4, T5), plane wave through it (T3).

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

Guided modes of lossless layers are poles on the real axis of u (on the real s axis of the
evanescent branch v = i s). Both forms therefore integrate below the real axis up to
U = max(n)/n₁ + 1, past every guided mode, and along it beyond: the causal limit of a
vanishing absorption, which moves the poles above the axis. On the real axis alone the power
into the guided modes is lost (found by hp-FEM, issue #23); the option ``contour="real"`` is
kept for absorbing stacks and for the check of that limit.

Closed forms for the checks:

* :func:`mirror_rates`: perfect conductor from the field of the image dipole (perpendicular
  image parallel, parallel image antiparallel at distance 2z, dipole field as in
  jackson1999), F = 1 + Im(μ̂·E_image)/Im(μ̂·E_self) with Im(E_self) = k³p/(6πε).
* :func:`quasistatic_rates`: leading term for k₁z → 0 above an absorbing half-space from the
  electrostatic image β p, β = (ε₂ − ε₁)/(ε₂ + ε₁).
* :func:`lukosz_asymptotic`: Lukosz and Kunz (4.4)–(4.6) for a lossless half-space, k₁z ≫ 1.

Radiation patterns (T5): :func:`substrate_pattern` (angular spectrum and transmission of the
stack), :func:`cover_pattern` (lukosz1979 (3.21), (3.22)), :func:`radiated_powers`, and the
closed forms :func:`lukosz_half_space_pattern` (lukosz1979 (3.7), (3.8), (3.11)–(3.15)).

Excitation field (T3): :func:`plane_wave_field` and :func:`plane_wave_rt` for a plane wave
incident from the substrate.
"""

from __future__ import annotations

import cmath
import math
import warnings
from collections.abc import Callable, Sequence
from itertools import pairwise

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import IntegrationWarning, quad

Layer = tuple[complex, float]  # (refractive index, thickness in the units of the wavelength)

PATH_DEPTH = 0.05  # depth δ of the integration path below the real u axis


def _segments(lo: float, hi: float, width: float) -> list[tuple[float, float]]:
    n = max(1, math.ceil((hi - lo) / width))
    edges = [float(e) for e in np.linspace(lo, hi, n + 1)]
    return list(pairwise(edges))


def _integrate(f: Callable[[float], float], segments: Sequence[tuple[float, float]]) -> float:
    """Sum of adaptive quadratures; a quadrature warning (e.g. a pole on the path) is an error."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", IntegrationWarning)
        return float(
            sum(quad(f, a, b, limit=200, epsabs=1e-13, epsrel=1e-11)[0] for a, b in segments)
        )


def _evanescent_extent(n1: float, layers: Sequence[Layer], substrate: complex, k1z: float) -> float:
    """Upper end of u (or s): past all modes and the decay of exp(−2 k₁z √(u²−1))."""
    n_max = max([abs(substrate)] + [abs(n) for n, _ in layers]) / n1
    return max(4.0 * n_max + 4.0, 40.0 / k1z)


def _path_end(n1: float, layers: Sequence[Layer], substrate: complex) -> float:
    """End U of the lowered path in u: beyond the largest index, so past all guided modes."""
    return max([abs(substrate)] + [abs(n) for n, _ in layers] + [n1]) / n1 + 1.0


def _lowered(s: float, end: float) -> tuple[complex, complex]:
    """Point and derivative of the path s − iδ sin(πs/S) below the real axis, s ∈ [0, S]."""
    w = math.pi / end
    return s - 1j * PATH_DEPTH * math.sin(w * s), 1.0 - 1j * PATH_DEPTH * w * math.cos(w * s)


def _l1(u: complex) -> complex:
    """l₁ = −i(1 − u²)^½, continuous on and below the real axis (l₁ = √(u²−1) for real u > 1)."""
    return -1j * cmath.sqrt(1.0 - u * u)


# --- Chance, Prock, Silbey -------------------------------------------------------------------


def _cps_reflection(
    u: complex, n1: float, layers: Sequence[Layer], substrate: complex, k1: float
) -> tuple[complex, complex]:
    """(R^∥, R^⊥) seen from medium 1 in the CPS convention, layers from the top."""
    eps = [n1**2] + [complex(n) ** 2 for n, _ in layers] + [complex(substrate) ** 2]
    ell = [-1j * cmath.sqrt(e / eps[0] - u * u) for e in eps]
    ell[0] = _l1(u)
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
    contour: str = "below",
) -> tuple[float, float]:
    """(F_⊥, F_∥) after chance1978 (2.17), (2.29) with q = 1.

    ``contour``: "below" (default) runs from u = 0 to U below the real axis and on along it,
    which includes the guided modes of lossless layers (poles on the real axis); "real" stays on
    the real axis and is valid only when every pole lies off it (absorbing media).
    """
    if contour not in ("below", "real"):
        raise ValueError("contour must be 'below' or 'real'")
    k1 = 2.0 * math.pi * n1 / wavelength
    dh = k1 * z

    def refl(u: complex) -> tuple[complex, complex]:
        if reflection is not None:
            return reflection(u)
        return _cps_reflection(u, n1, layers, substrate, k1)

    def kernels(u: complex) -> tuple[complex, complex]:
        """Integrands of (2.17) and (2.29) per du."""
        l1 = _l1(u)
        r_par, r_perp = refl(u)
        e = cmath.exp(-2.0 * l1 * dh)
        return r_par * e * u**3 / l1, ((1.0 - u * u) * r_par + r_perp) * e * u / l1

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
    if contour == "real":
        # u = sin θ on [0, 1]: u³/l₁ du = i u³ dθ, u/l₁ du = i u dθ, e^{−2l₁d̂} = e^{2i d̂ cos θ}
        def perp_head(x: float) -> float:  # x = θ
            u = math.sin(x)
            r_par, _ = refl(u)
            return (r_par * cmath.exp(2j * dh * math.cos(x)) * 1j * u**3).imag

        def par_head(x: float) -> float:
            u = math.sin(x)
            r_par, r_perp = refl(u)
            val = ((1.0 - u * u) * r_par + r_perp) * cmath.exp(2j * dh * math.cos(x)) * 1j * u
            return val.imag

        head = _segments(0.0, 0.5 * math.pi, 0.05)
        t_head = 0.0
    else:
        end = _path_end(n1, layers, substrate)

        def perp_head(x: float) -> float:  # x = s on the lowered path
            u, du = _lowered(x, end)
            return (kernels(u)[0] * du).imag

        def par_head(x: float) -> float:
            u, du = _lowered(x, end)
            return (kernels(u)[1] * du).imag

        head = _segments(0.0, end, 0.02)
        t_head = math.acosh(end)
    tail = _segments(t_head, max(t_max, t_head), 0.01)
    perp = 1.0 - 1.5 * (_integrate(perp_head, head) + _integrate(perp_ev, tail))
    par = 1.0 + 0.75 * (_integrate(par_head, head) + _integrate(par_ev, tail))
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
    n1: float,
    layers: Sequence[Layer],
    substrate: complex,
    wavelength: float,
    z: float,
    contour: str = "below",
) -> tuple[float, float]:
    """(F_⊥, F_∥) after lukosz1977 (3.15)–(3.17); ``contour`` as in :func:`cps_rates`.

    On the imaginary branch v = i s the guided modes are poles on the real s axis; "below" runs
    s on σ − iδ sin(πσ/S) (u below the real axis as in :func:`cps_rates`), then along the axis.
    """
    if contour not in ("below", "real"):
        raise ValueError("contour must be 'below' or 'real'")
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

    s_max = math.sqrt(_evanescent_extent(n1, layers, substrate, k1 * z) ** 2 - 1.0)
    s_end = math.sqrt(_path_end(n1, layers, substrate) ** 2 - 1.0) if contour == "below" else 0.0

    def on_path(kernel: Callable[[complex], complex]) -> float:
        # ∫₀¹ f(v) e^{2ik₁zv} dv + ∫_{i∞}^0 f(v) e^{2ik₁zv} dv, the second with v = i s
        a = _integrate(
            lambda v: (kernel(v) * cmath.exp(2j * k1 * z * v)).real, _segments(0, 1, 0.05)
        )

        def lowered(sig: float) -> float:
            s, ds = _lowered(sig, s_end)
            return (-1j * kernel(1j * s) * cmath.exp(-2.0 * k1 * z * s) * ds).real

        b = _integrate(lowered, _segments(0.0, s_end, 0.02)) if s_end > 0.0 else 0.0
        c = _integrate(
            lambda s: (-1j * kernel(1j * s)).real * math.exp(-2.0 * k1 * z * s),
            _segments(s_end, max(s_max, s_end), 0.01),
        )
        return a + b + c

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


# --- radiation patterns (T5) -----------------------------------------------------------------


def _stack_coefficients(
    kz: Sequence[complex], eps: Sequence[complex], thick: Sequence[float]
) -> tuple[complex, complex, complex, complex]:
    """(r_s, r_p, t_s, t_p) of the stack seen from medium 0, Airy recursion.

    r_p as in lukosz1977 (2.35) (magnetic-field amplitudes); t_s and t_p for the electric field
    of the wave transmitted into the last medium (t_p^E = (n₀/n_N) t_p^H).
    """

    def pair(i: int, j: int) -> tuple[complex, complex]:
        rs = (kz[i] - kz[j]) / (kz[i] + kz[j])
        rp = (eps[j] * kz[i] - eps[i] * kz[j]) / (eps[j] * kz[i] + eps[i] * kz[j])
        return rs, rp

    n = len(kz) - 1
    rs, rp = pair(n - 1, n)
    ts, tp = 1.0 + rs, 1.0 + rp
    for j in range(n - 1, 0, -1):
        s, p = pair(j - 1, j)
        ph = cmath.exp(1j * kz[j] * thick[j - 1])
        ds, dp = 1.0 + s * rs * ph * ph, 1.0 + p * rp * ph * ph
        ts, tp = (1.0 + s) * ts * ph / ds, (1.0 + p) * tp * ph / dp
        rs, rp = (s + rs * ph * ph) / ds, (p + rp * ph * ph) / dp
    return rs, rp, ts, tp * cmath.sqrt(eps[0]) / cmath.sqrt(eps[n])


def _kz_list(
    n1: float, layers: Sequence[Layer], substrate: complex, k0: float, kpar: complex
) -> tuple[list[complex], list[complex], list[float]]:
    eps = [complex(n1**2)] + [complex(n) ** 2 for n, _ in layers] + [complex(substrate) ** 2]
    kz = []
    for e in eps:
        w = cmath.sqrt(e * k0 * k0 - kpar * kpar)
        kz.append(w if w.imag >= 0.0 else -w)
    return kz, eps, [d for _, d in layers]


def substrate_pattern(
    n1: float, layers: Sequence[Layer], substrate: float, wavelength: float, z: float, alpha: float
) -> tuple[float, float, float]:
    """Radiation into a lossless substrate per solid angle, normalised as in lukosz1979 (2.8).

    Angle α from −z in the substrate. Returns the coefficients (a_⊥, a_s, a_p) of
    P_⊥ = a_⊥ (perpendicular dipole), P_∥ = a_s sin²φ + a_p cos²φ (dipole along x).

    From the angular spectrum of the dipole and the transmission of the stack: the downward wave
    k_∥ = k₃ sin α of the dipole has the s and p amplitudes ∝ (p̂·ŝ), (p̂·ê_p) e^{ik_z1 z} / k_z1
    with ê_p = (k_z1 cos φ, k_z1 sin φ, k_∥)/k₁; the far field by stationary phase gives
    P = (3/8π) (n₃³/n₁) k₀² cos²α |t|² |p̂·ê|² e^{−2 Im(k_z1) z} / |k_z1|². For a half-space this
    is lukosz1979 (3.7), (3.8) below and (3.11)–(3.15) above the critical angle.
    """
    k0 = 2.0 * math.pi / wavelength
    k1 = n1 * k0
    kpar = substrate * k0 * math.sin(alpha)
    kz, eps, thick = _kz_list(n1, layers, substrate, k0, kpar)
    _, _, ts, tp = _stack_coefficients(kz, eps, thick)
    kz1 = kz[0]
    common = (
        3.0
        / (8.0 * math.pi)
        * substrate**3
        / n1
        * k0**2
        * math.cos(alpha) ** 2
        * math.exp(-2.0 * kz1.imag * z)
        / abs(kz1) ** 2
    )
    a_perp = common * abs(tp) ** 2 * abs(kpar / k1) ** 2
    a_s = common * abs(ts) ** 2
    a_p = common * abs(tp) ** 2 * abs(kz1 / k1) ** 2
    return a_perp, a_s, a_p


def cover_pattern(
    n1: float,
    layers: Sequence[Layer],
    substrate: complex,
    wavelength: float,
    z: float,
    alpha: float,
) -> tuple[float, float, float]:
    """Radiation into the cover per solid angle, lukosz1979 (3.21), (3.22) for the electric dipole
    ((m) → (e), (s) ↔ (p)); α from +z. Returns (a_⊥, a_s, a_p) as in :func:`substrate_pattern`.
    """
    k0 = 2.0 * math.pi / wavelength
    kpar = n1 * k0 * math.sin(alpha)
    kz, eps, thick = _kz_list(n1, layers, substrate, k0, kpar)
    rs, rp, _, _ = _stack_coefficients(kz, eps, thick)
    ph = cmath.exp(2j * kz[0] * z)
    c = 3.0 / (8.0 * math.pi)
    a_perp = c * math.sin(alpha) ** 2 * (1.0 + abs(rp) ** 2 + 2.0 * (rp * ph).real)
    a_p = c * math.cos(alpha) ** 2 * (1.0 + abs(rp) ** 2 - 2.0 * (rp * ph).real)
    a_s = c * (1.0 + abs(rs) ** 2 + 2.0 * (rs * ph).real)
    return a_perp, a_s, a_p


def radiated_powers(
    n1: float, layers: Sequence[Layer], substrate: float, wavelength: float, z: float
) -> tuple[tuple[float, float], tuple[float, float]]:
    """((cover ⊥, cover ∥), (substrate ⊥, substrate ∥)), integrated over the half-spaces.

    φ integrals: 2π for ⊥, π for sin²φ and cos²φ; the cross term of (3.21) vanishes.
    """
    out = []
    for pattern, n_side in ((cover_pattern, n1), (substrate_pattern, substrate)):
        crit = math.asin(min(1.0, n1 / n_side)) if pattern is substrate_pattern else 0.0
        edges = sorted({0.0, crit, 0.5 * math.pi})
        segs = [s for a, b in pairwise(edges) for s in _segments(a, b, 0.005)]

        def perp(a: float, f: Callable[..., tuple[float, float, float]] = pattern) -> float:
            return 2.0 * math.pi * f(n1, layers, substrate, wavelength, z, a)[0] * math.sin(a)

        def par(a: float, f: Callable[..., tuple[float, float, float]] = pattern) -> float:
            v = f(n1, layers, substrate, wavelength, z, a)
            return math.pi * (v[1] + v[2]) * math.sin(a)

        out.append((_integrate(perp, segs), _integrate(par, segs)))
    return (out[0][0], out[0][1]), (out[1][0], out[1][1])


def lukosz_half_space_pattern(
    n1: float, n2: float, wavelength: float, z: float, alpha: float
) -> tuple[float, float, float]:
    """lukosz1979 (3.7), (3.8), (3.11)–(3.15): (a_⊥, a_s, a_p) into a dielectric half-space."""
    n = n2 / n1
    lam1 = wavelength / n1
    s2, c2 = math.sin(alpha), math.cos(alpha)
    if n * s2 <= 1.0:
        c1 = math.sqrt(1.0 - (n * s2) ** 2)
        g = 3.0 / (2.0 * math.pi) * n**3 * c2**2
        a_perp = g * (n * s2) ** 2 / (n * c1 + c2) ** 2  # (3.7) at Θ = 0, sin α₁ = n sin α₂
        a_s = g / (c1 + n * c2) ** 2
        a_p = g * c1**2 / (n * c1 + c2) ** 2
        return a_perp, a_s, a_p
    dz = lam1 / (2.0 * math.pi) / math.sqrt(n * n * s2 * s2 - 1.0)
    f_s = 3.0 / (2.0 * math.pi) * n**3 / (n * n - 1.0) * c2**2 * math.exp(-2.0 * z / dz)
    f_p = f_s / ((n * n + 1.0) * s2 * s2 - 1.0)
    return f_p * n * n * s2 * s2, f_s, f_p * (n * n * s2 * s2 - 1.0)


# --- excitation field of a plane wave (T3) ---------------------------------------------------

Z0 = 1.0 / (8.8541878128e-12 * 299792458.0)  # vacuum impedance [Ω], CODATA 2018 ε₀

Profile = Callable[[float], tuple[complex, complex, complex]]


def _profile(
    cover: float, layers: Sequence[Layer], substrate: float, k0: float, kpar: float, pol: str
) -> tuple[Profile, complex, complex]:
    """Field profile ψ(z) of one polarisation, wave incident from the substrate.

    ψ is E_s ("s") or Z₀H_s ("p"); ψ and ψ'/q are continuous (q = 1 for s, ε for p). The top
    interface is z = 0, the layers lie below. Starting from the transmitted wave ψ = e^{ik_z z}
    in the cover, the pair (ψ, ψ') is carried down through each layer in closed form.
    Returns ψ(z) as (ψ, ψ', ε) and the amplitudes (a, b) of the incident and the reflected wave
    in the substrate at its interface z_L: ψ = a e^{ik_z(z−z_L)} + b e^{−ik_z(z−z_L)}.
    """
    eps = [complex(cover) ** 2] + [complex(n) ** 2 for n, _ in layers] + [complex(substrate) ** 2]
    kz = []
    for e in eps:
        w = cmath.sqrt(e * k0 * k0 - kpar * kpar)
        kz.append(w if w.imag >= 0.0 else -w)
    q = [1.0 + 0j] * len(eps) if pol == "s" else eps
    tops = [0.0]  # upper boundary of medium j ≥ 1
    pairs = [(1.0 + 0j, 1j * kz[0])]  # (ψ, ψ') just below that boundary, medium j
    for j in range(1, len(eps)):
        psi, dpsi = pairs[-1]
        if j > 1:  # carry through layer j − 1
            k, d = kz[j - 1], layers[j - 2][1]
            psi, dpsi = (
                psi * cmath.cos(k * d) - dpsi / k * cmath.sin(k * d),
                dpsi * cmath.cos(k * d) + psi * k * cmath.sin(k * d),
            )
            tops.append(tops[-1] - d)
        pairs.append((psi, dpsi * q[j] / q[j - 1]))
    del pairs[0]
    zl = tops[-1]
    psi_l, dpsi_l = pairs[-1]
    a = 0.5 * (psi_l + dpsi_l / (1j * kz[-1]))
    b = 0.5 * (psi_l - dpsi_l / (1j * kz[-1]))

    def at(z: float) -> tuple[complex, complex, complex]:
        if z >= 0.0:
            v = cmath.exp(1j * kz[0] * z)
            return v, 1j * kz[0] * v, eps[0]
        if z <= zl:
            k = kz[-1]
            up, down = a * cmath.exp(1j * k * (z - zl)), b * cmath.exp(-1j * k * (z - zl))
            return up + down, 1j * k * (up - down), eps[-1]
        j = max(i for i, top in enumerate(tops) if top > z)  # z inside layer j + 1
        psi, dpsi = pairs[j]
        k, d = kz[j + 1], tops[j] - z
        return (
            psi * cmath.cos(k * d) - dpsi / k * cmath.sin(k * d),
            dpsi * cmath.cos(k * d) + psi * k * cmath.sin(k * d),
            eps[j + 1],
        )

    return at, a, b


def plane_wave_rt(
    cover: float,
    layers: Sequence[Layer],
    substrate: float,
    wavelength: float,
    theta: float,
    pol: str,
) -> tuple[complex, complex]:
    """(r, t) of the stack for a wave from the substrate at the angle θ (in the substrate).

    r = b/a at the substrate interface, t = ψ(0⁺)/a; for "p" both refer to the magnetic field
    (raether1988 (2.16), (2.25) with medium 0 = substrate, 2 = cover).
    """
    k0 = 2.0 * math.pi / wavelength
    _, a, b = _profile(cover, layers, substrate, k0, substrate * k0 * math.sin(theta), pol)
    return b / a, 1.0 / a


def plane_wave_field(
    cover: float,
    layers: Sequence[Layer],
    substrate: float,
    wavelength: float,
    k_hat: ArrayLike,
    polarisation: ArrayLike,
    irradiance: float,
    points: ArrayLike,
) -> NDArray[np.complex128]:
    """Electric field (N, 3) [V/m] of a plane wave from the substrate at the points (N, 3).

    The incident wave is E₀ ê e^{ik·r} (factor e^{−iωt}, zero phase at the origin, continued
    through the stack as if homogeneous), k = n_s k₀ k̂ with k̂_z > 0, irradiance
    I₀ = n_s |E₀|²/(2Z₀) [W/m²] perpendicular to the beam. Top interface at z = 0, layers
    (index, thickness) from the top below it; lengths in the units of the wavelength. ê is
    split into s (ŝ = ẑ × k̂, ŷ at normal incidence) and p (p̂ = ŝ × k̂); for p, E follows
    from Z₀H_s = ψ as E_u = −iψ'/(k₀ε), E_z = −k_u ψ/(k₀ε), û = ŝ × ẑ.
    """
    kh = np.asarray(k_hat, dtype=np.float64)
    if kh[2] <= 0.0:
        raise ValueError("the wave must come from the substrate (k_hat[2] > 0)")
    e_hat = np.asarray(polarisation, dtype=np.complex128)
    k0 = 2.0 * math.pi / wavelength
    zk = np.cross([0.0, 0.0, 1.0], kh)
    s_hat = zk / np.linalg.norm(zk) if np.linalg.norm(zk) > 1e-15 else np.array([0.0, 1.0, 0.0])
    p_hat = np.cross(s_hat, kh)
    u_hat = np.cross(s_hat, [0.0, 0.0, 1.0])
    kpar = substrate * k0 * float(kh @ u_hat)
    e0 = math.sqrt(2.0 * Z0 * irradiance / substrate)
    kz_s = substrate * k0 * float(kh[2])
    pts = np.atleast_2d(np.asarray(points, dtype=np.float64))
    out = np.zeros((pts.shape[0], 3), dtype=np.complex128)
    for pol, amp in (("s", e0 * (s_hat @ e_hat)), ("p", substrate * e0 * (p_hat @ e_hat))):
        if amp == 0.0:
            continue
        at, a, _ = _profile(cover, layers, substrate, k0, kpar, pol)
        zl = -sum(d for _, d in layers)
        scale = amp / (a * cmath.exp(-1j * kz_s * zl))  # incident amplitude at the origin
        for i in range(pts.shape[0]):
            psi, dpsi, eps = at(float(pts[i, 2]))
            ph = scale * cmath.exp(1j * kpar * float(u_hat[:2] @ pts[i, :2]))
            if pol == "s":
                out[i] += ph * psi * s_hat
            else:
                out[i] += (
                    ph * (-1j * dpsi * u_hat - kpar * psi * np.array([0.0, 0.0, 1.0])) / (k0 * eps)
                )
    return out
