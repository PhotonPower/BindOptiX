"""T3: excitation field of a plane wave at a gold film (Kretschmann, p polarisation).

The reference (tests/planar_reference.py) is checked against the homogeneous medium, against
the closed three-layer forms of raether1988 and against the energy balance; then it serves as
the target for a field map exported by a light solver (hp-FEM), if one is given in
BINDOPTIX_T3_MAP.
"""

import cmath
import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np
import pytest
from numpy.typing import NDArray
from planar_reference import Z0, plane_wave_field, plane_wave_rt

from bindoptix.ports import FieldMap

LAM = 0.67  # vacuum wavelength [µm]
N_WATER = 1.333
N_GLASS = 1.52
N_GOLD = cmath.sqrt(-16.0 + 1.1j)  # gold-like test permittivity, not a material value
I0 = 1.0e7  # irradiance [W/m²]


def _k_hat(theta: float, phi: float = 0.0) -> NDArray[np.float64]:
    return np.array(
        [math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta)]
    )


# --- the reference against closed forms ----------------------------------------------------


@pytest.mark.parametrize(("theta", "phi"), [(0.0, 0.0), (0.7, 0.4), (1.2, 2.5)])
def test_t3_homogeneous_medium_gives_the_incident_wave(theta: float, phi: float) -> None:
    # Without contrast the field is E₀ ê e^{ik·r} everywhere, I₀ = n|E₀|²/(2Z₀).
    kh = _k_hat(theta, phi)
    a = np.cross(kh, [0.3, -0.5, 0.8])
    b = np.cross(kh, a)
    e_hat = (a / np.linalg.norm(a) + 0.6j * b / np.linalg.norm(b)) / math.sqrt(1.36)
    pts = np.array([[0.1, -0.2, 0.3], [0.0, 0.0, -0.04], [-0.3, 0.2, -0.5]])
    e = plane_wave_field(N_GLASS, [(N_GLASS, 0.05)], N_GLASS, LAM, kh, e_hat, I0, pts)
    e0 = math.sqrt(2.0 * Z0 * I0 / N_GLASS)
    k = 2.0 * math.pi / LAM * N_GLASS * kh
    expect = e0 * e_hat[None, :] * np.exp(1j * pts @ k)[:, None]
    np.testing.assert_allclose(e, expect, rtol=1e-12, atol=1e-12 * e0)


def _raether_three_layer(theta: float, d: float) -> tuple[complex, complex]:
    """(r₀₁₂, t₀₁₂) for p after raether1988 (2.16), (2.17), (2.25), (2.26); 0 glass, 1 gold."""
    k0 = 2.0 * math.pi / LAM
    kx = N_GLASS * k0 * math.sin(theta)
    eps = [N_GLASS**2, N_GOLD**2, N_WATER**2]
    kz = [cmath.sqrt(e * k0 * k0 - kx * kx) for e in eps]
    kz = [k if k.imag >= 0.0 else -k for k in kz]

    def r(i: int, k: int) -> complex:
        return (kz[i] / eps[i] - kz[k] / eps[k]) / (kz[i] / eps[i] + kz[k] / eps[k])

    ph = cmath.exp(1j * kz[1] * d)
    den = 1.0 + r(0, 1) * r(1, 2) * ph * ph
    return (r(0, 1) + r(1, 2) * ph * ph) / den, (1.0 + r(0, 1)) * (1.0 + r(1, 2)) * ph / den


@pytest.mark.parametrize("theta_deg", [30.0, 62.0, 68.4, 75.0])  # critical angle 61.3°
@pytest.mark.parametrize("d", [0.0, 0.03, 0.05])
def test_t3_reflection_and_enhancement_equal_raether(theta_deg: float, d: float) -> None:
    th = math.radians(theta_deg)
    r, t = plane_wave_rt(N_WATER, [(N_GOLD, d)], N_GLASS, LAM, th, "p")
    r_ref, t_ref = _raether_three_layer(th, d)
    assert abs(r - r_ref) <= 1e-12 * abs(r_ref)
    assert abs(t - t_ref) <= 1e-12 * abs(t_ref)


@pytest.mark.parametrize("pol", ["s", "p"])
@pytest.mark.parametrize("theta_deg", [30.0, 68.4, 75.0])
def test_t3_reflected_transmitted_and_absorbed_power_add_up(pol: str, theta_deg: float) -> None:
    # Poynting's theorem for harmonic fields (jackson1999): R + T + A = 1 with A from
    # (ω ε₀/2) ∫ Im ε |E|² dz in the film, all per incident flux I₀ cos θ.
    th = math.radians(theta_deg)
    d = 0.05
    kh = _k_hat(th)
    e_hat = (
        np.array([0.0, 1.0, 0.0]) if pol == "s" else np.array([math.cos(th), 0.0, -math.sin(th)])
    )
    layers = [(N_GOLD, d)]
    k0 = 2.0 * math.pi / LAM
    x, w = np.polynomial.legendre.leggauss(80)
    zs = -0.5 * d * (x + 1.0)
    e = plane_wave_field(N_WATER, layers, N_GLASS, LAM, kh, e_hat, I0, np.c_[0 * zs, 0 * zs, zs])
    absorbed = 0.5 * k0 / Z0 * (N_GOLD**2).imag * 0.5 * d * float(w @ np.sum(abs(e) ** 2, axis=1))
    flux = I0 * math.cos(th)
    r, _ = plane_wave_rt(N_WATER, layers, N_GLASS, LAM, th, pol)
    kz_c = cmath.sqrt(N_WATER**2 - (N_GLASS * math.sin(th)) ** 2)
    e_c = plane_wave_field(N_WATER, layers, N_GLASS, LAM, kh, e_hat, I0, [[0.0, 0.0, 0.0]])[0]
    transmitted = kz_c.real * float(np.sum(abs(e_c) ** 2)) / (2.0 * Z0)
    total = abs(r) ** 2 + (transmitted + absorbed) / flux
    assert total == pytest.approx(1.0, abs=1e-10)


# --- the exported map ----------------------------------------------------------------------

MAP = os.environ.get("BINDOPTIX_T3_MAP")


@pytest.mark.skipif(not MAP, reason="no exported field map given (BINDOPTIX_T3_MAP)")
def test_t3_exported_map_matches_reference() -> None:
    """E of a planar field map from hp-FEM, normalised to I₀, against the reference.

    Stack in ``source["stack"]`` as for T4 (top interface at z = 0, layers below), incidence
    from the substrate (k̂_z > 0), ``irradiance_w_m2``, ``k_hat`` and ``polarisation`` set. The
    global phase is fitted (the point of zero phase is free); grid nodes within 1 pm of an
    interface are left out (E_z jumps there).
    """
    m = FieldMap.from_json(Path(str(MAP)).read_text(encoding="utf-8"))
    assert m.irradiance_w_m2 is not None and m.k_hat is not None and m.polarisation is not None
    stack = m.source["stack"]
    if "content_hash" in m.source:
        assert hashlib.sha256(stack.encode("utf-8")).hexdigest() == m.source["content_hash"]
    st = json.loads(stack)
    lam = float(st["wavelength_m"]) * 1e6
    assert lam == pytest.approx(m.wavelength_um, rel=1e-12)
    cover = float(st["cover"][0])
    layers = [(complex(re, im), d * 1e6) for re, im, d in st["layers"]]
    sub = complex(*st["substrate"])
    assert sub.imag == 0.0, "T3 needs a lossless substrate"
    gx, gy, gz = np.meshgrid(*m.grid.axes, indexing="ij")
    pts = np.c_[gx.ravel(), gy.ravel(), gz.ravel()]
    bounds = np.cumsum([0.0] + [-d for _, d in layers])
    keep = np.min(np.abs(pts[:, 2:3] - bounds[None, :]), axis=1) > 1e-6
    e = m.e.reshape(-1, 3)[keep]
    ref = plane_wave_field(
        cover, layers, sub.real, lam, m.k_hat, m.polarisation, m.irradiance_w_m2, pts[keep]
    )
    c = np.vdot(ref, e) / np.vdot(ref, ref)
    assert abs(c) == pytest.approx(1.0, rel=1e-6)
    np.testing.assert_allclose(e, c * ref, rtol=1e-6, atol=1e-6 * np.abs(ref).max())
