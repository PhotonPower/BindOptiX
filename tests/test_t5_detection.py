"""T5: radiation into the substrate (supercritical part, SPCE) and the detection tensors.

Reference in tests/planar_reference.py; an exported table (hp-FEM) is compared when given in
BINDOPTIX_T5_TABLE.
"""

import cmath
import math
import os
from pathlib import Path

import numpy as np
import pytest
from planar_reference import (
    _kz_list,
    _stack_coefficients,
    cps_rates,
    lukosz_half_space_pattern,
    radiated_powers,
    substrate_pattern,
)

from bindoptix.ports import EmitterResponseTable

LAM = 0.67  # vacuum wavelength [µm]
N_WATER = 1.333
N_GLASS = 1.52
GOLD_FILM = [(cmath.sqrt(-16.0 + 1.1j), 0.05)]  # gold-like test permittivity, 50 nm


@pytest.mark.parametrize("z", [0.0, 0.05])
@pytest.mark.parametrize("alpha", [0.2, 0.9, 1.1, 1.4])  # critical angle 1.07 rad
def test_t5_substrate_pattern_equals_lukosz_closed_form(z: float, alpha: float) -> None:
    ref = lukosz_half_space_pattern(N_WATER, N_GLASS, LAM, z, alpha)
    np.testing.assert_allclose(substrate_pattern(N_WATER, [], N_GLASS, LAM, z, alpha), ref, 1e-12)


@pytest.mark.parametrize("z", [0.01, 0.1])
def test_t5_lossless_half_space_radiates_the_total_power(z: float) -> None:
    # No absorption, no guided modes: cover + substrate = A_tot of the T4 reference.
    (up_perp, up_par), (dn_perp, dn_par) = radiated_powers(N_WATER, [], N_GLASS, LAM, z)
    tot_perp, tot_par = cps_rates(N_WATER, [], N_GLASS, LAM, z)
    np.testing.assert_allclose([up_perp + dn_perp, up_par + dn_par], [tot_perp, tot_par], 1e-9)


def test_t5_spce_angle_is_the_spr_angle() -> None:
    # lakowicz2004: a fluorophore near the film emits p-polarised light into the glass at the
    # SPR angle (reflectance minimum, Kretschmann). The emission maximum must lie well inside
    # the resonance: within a quarter of the width of the reflectance dip.
    k0 = 2.0 * math.pi / LAM
    crit = math.asin(N_WATER / N_GLASS)
    alpha = np.linspace(crit + 1e-4, 0.5 * math.pi - 1e-3, 20001)
    emission = np.array(
        [substrate_pattern(N_WATER, GOLD_FILM, N_GLASS, LAM, 0.02, a)[0] for a in alpha]
    )

    def reflectance(a: float) -> float:
        kz, eps, th = _kz_list(N_GLASS, GOLD_FILM[::-1], N_WATER, k0, N_GLASS * k0 * math.sin(a))
        return abs(_stack_coefficients(kz, eps, th)[1]) ** 2

    refl = np.array([reflectance(a) for a in alpha])
    dip = alpha[refl < 0.5 * (refl.max() + refl.min())]
    width = dip.max() - dip.min()
    assert abs(alpha[emission.argmax()] - alpha[refl.argmin()]) < 0.25 * width


TABLE = os.environ.get("BINDOPTIX_T5_TABLE")


@pytest.mark.skipif(not TABLE, reason="no exported table given (BINDOPTIX_T5_TABLE)")
def test_t5_exported_table_matches_reference() -> None:
    """A_rad and A_det of a planar table from hp-FEM against the reference.

    Stack in ``source["stack"]`` as for T4; A_det is compared for the default regions
    ("solution", "substrate") = whole cover and whole substrate.
    """
    import json

    table = EmitterResponseTable.from_json(Path(str(TABLE)).read_text(encoding="utf-8"))
    st = json.loads(table.source["stack"])
    n1 = float(st["cover"][0])
    lam = float(st["wavelength_m"]) * 1e6
    layers = [(complex(re, im), d * 1e6) for re, im, d in st["layers"]]
    sub = complex(*st["substrate"])
    assert sub.imag == 0.0, "T5 needs a lossless substrate"
    heights = table.grid.axes[2]
    picks = np.unique(np.linspace(0, heights.size - 1, min(heights.size, 6)).round().astype(int))
    default = table.detectors == ("solution", "substrate")
    for k in picks:
        (up_perp, up_par), (dn_perp, dn_par) = radiated_powers(
            n1, layers, sub.real, lam, float(heights[k])
        )
        rad = table.a_rad[0, 0, k]
        expect = [up_par + dn_par, up_par + dn_par, up_perp + dn_perp]
        np.testing.assert_allclose(np.diag(rad), expect, rtol=1e-6)
        if default:
            for d, (perp, par) in enumerate(((up_perp, up_par), (dn_perp, dn_par))):
                np.testing.assert_allclose(
                    np.diag(table.a_det[0, 0, k, d]), [par, par, perp], rtol=1e-6
                )
