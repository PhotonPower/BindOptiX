"""T4: total power (decay rate) of a dipole above glass / gold / water against the distance.

The reference (tests/planar_reference.py) is checked first against closed forms and against
a second, separately written form; then it serves as the target for the A_tot of an emitter
response table exported by a light solver (hp-FEM), if one is given in BINDOPTIX_T4_TABLE.
"""

import cmath
import json
import math
import os
from pathlib import Path

import numpy as np
import pytest
from planar_reference import (
    cps_rates,
    lukosz_asymptotic,
    lukosz_rates,
    mirror_rates,
    quasistatic_rates,
)

from bindoptix.ports import EmitterResponseTable

LAM = 0.67  # vacuum wavelength [µm]
N_WATER = 1.333
N_GLASS = 1.52
EPS_TEST = -16.0 + 1.1j  # gold-like test permittivity, not a material value
K1 = 2.0 * math.pi * N_WATER / LAM


# --- the reference against closed forms ----------------------------------------------------


@pytest.mark.parametrize("k1z", [0.05, 0.5, 2.0, 7.0])
def test_t4_perfect_mirror_equals_image_dipole(k1z: float) -> None:
    ref = cps_rates(N_WATER, [], 1.0, LAM, k1z / K1, reflection=lambda u: (-1.0 + 0j, -1.0 + 0j))
    np.testing.assert_allclose(ref, mirror_rates(k1z), rtol=1e-10, atol=1e-12)


def test_t4_mirror_limits() -> None:
    # Close to a perfect mirror: perpendicular doubled, parallel cancelled; far: free dipole.
    assert mirror_rates(1e-3) == pytest.approx((2.0, 0.0), abs=1e-5)
    assert mirror_rates(1e4) == pytest.approx((1.0, 1.0), abs=2e-4)


def test_t4_quasistatic_limit_above_absorbing_half_space() -> None:
    # The relative deviation from the electrostatic image term vanishes for k₁z → 0, at least
    # quadratically (next terms of the near-field expansion).
    n2 = cmath.sqrt(EPS_TEST)
    err = []
    for k1z in (0.01, 0.005):
        ex = cps_rates(N_WATER, [], n2, LAM, k1z / K1)
        qs = quasistatic_rates(N_WATER**2, EPS_TEST, k1z)
        err.append(max(abs(e - q) / q for e, q in zip(ex, qs, strict=True)))
    assert err[0] < 2e-3
    assert err[1] < err[0] / 4.0


@pytest.mark.parametrize("z_over_lambda1", [1.0, 2.0, 4.0])
def test_t4_lossless_half_space_approaches_lukosz_asymptotics(z_over_lambda1: float) -> None:
    # lukosz1977 (4.5), (4.6) keep the terms up to z̄⁻¹ (parallel) and z̄⁻³ (perpendicular);
    # the deviation must stay below the first neglected order, 0.1/z̄².
    n = 1.5
    z = z_over_lambda1 * LAM / N_WATER
    ex = cps_rates(N_WATER, [], n * N_WATER, LAM, z)
    asym = lukosz_asymptotic(n, K1 * z)
    zbar = 2.0 * K1 * z
    for e, a in zip(ex, asym, strict=True):
        assert abs(e - a) * zbar**2 < 0.1


# --- two separately written forms ----------------------------------------------------------

STACKS = {
    "glass": ([], N_GLASS),
    "gold film on glass": ([(cmath.sqrt(EPS_TEST), 0.05)], N_GLASS),
}


@pytest.mark.parametrize("stack", sorted(STACKS))
@pytest.mark.parametrize("z", [0.005, 0.05, 0.3])
def test_t4_chance_prock_silbey_equals_lukosz_kunz(stack: str, z: float) -> None:
    layers, sub = STACKS[stack]
    cps = cps_rates(N_WATER, layers, sub, LAM, z)
    luk = lukosz_rates(N_WATER, layers, sub, LAM, z)
    np.testing.assert_allclose(cps, luk, rtol=1e-9)


def test_t4_film_limits() -> None:
    z = 0.02
    n_au = cmath.sqrt(EPS_TEST)
    thin = cps_rates(N_WATER, [(n_au, 0.0)], N_GLASS, LAM, z)
    np.testing.assert_allclose(thin, cps_rates(N_WATER, [], N_GLASS, LAM, z), rtol=1e-12)
    thick = cps_rates(N_WATER, [(n_au, 2.0)], N_GLASS, LAM, z)  # decay length ≈ 27 nm
    np.testing.assert_allclose(thick, cps_rates(N_WATER, [], n_au, LAM, z), rtol=1e-10)


# --- the exported table --------------------------------------------------------------------

TABLE = os.environ.get("BINDOPTIX_T4_TABLE")


@pytest.mark.skipif(not TABLE, reason="no exported table given (BINDOPTIX_T4_TABLE)")
def test_t4_exported_table_matches_reference() -> None:
    """A_tot of a planar emitter response table from hp-FEM against the reference.

    The table must name its stack in ``source["stack"]`` (JSON: ``cover`` [n, 0], ``layers``
    [[Re n, Im n, thickness m], …] from the top, ``substrate`` [Re n, Im n], ``wavelength_m``).
    """
    table = EmitterResponseTable.from_json(Path(str(TABLE)).read_text(encoding="utf-8"))
    assert "stack" in table.source, "the table must name its stack in source['stack']"
    st = json.loads(table.source["stack"])
    n1 = float(st["cover"][0])
    assert n1 == pytest.approx(table.n_host, rel=1e-12)
    lam = float(st["wavelength_m"]) * 1e6
    assert lam == pytest.approx(table.wavelength_um, rel=1e-12)
    layers = [(complex(re, im), d * 1e6) for re, im, d in st["layers"]]
    sub = complex(*st["substrate"])
    heights = table.grid.axes[2]
    picks = np.unique(np.linspace(0, heights.size - 1, min(heights.size, 8)).round().astype(int))
    for k in picks:
        a = table.a_tot[0, 0, k]
        perp, par = cps_rates(n1, layers, sub, lam, float(heights[k]))
        np.testing.assert_allclose(np.diag(a), [par, par, perp], rtol=1e-6)
        assert np.max(np.abs(a - np.diag(np.diag(a)))) <= 1e-12 * np.abs(a).max()
