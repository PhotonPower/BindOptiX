"""Units and constants (docs/konventionen.md).

Internal units: length µm, time s, temperature K, power mW. Concentrations are particle counts
per µm³, surface coverages per µm². All conversions of the package go through this module.
"""

from __future__ import annotations

# Exact SI defining constants (SI brochure, 9th edition, 2019).
AVOGADRO: float = 6.02214076e23
"""Avogadro constant N_A [1/mol]."""
BOLTZMANN: float = 1.380649e-23
"""Boltzmann constant k_B [J/K]."""

# Conventional value (CGPM 1901), used only when no local value is given.
STANDARD_GRAVITY_M_S2: float = 9.80665
"""Standard acceleration of gravity [m/s²]."""

UM3_PER_LITRE: float = 1e15
UM3_PER_ML: float = 1e12
MM_PER_UM: float = 1e-3
UM_PER_MM: float = 1e3
UM_PER_M: float = 1e6


def molar_to_per_um3(c_molar: float) -> float:
    """Molar concentration [mol/L] → particles per µm³."""
    return c_molar * AVOGADRO / UM3_PER_LITRE


def per_um3_to_molar(c_per_um3: float) -> float:
    """Particles per µm³ → molar concentration [mol/L]."""
    return c_per_um3 * UM3_PER_LITRE / AVOGADRO


def cfu_per_ml_to_per_um3(c_cfu_per_ml: float) -> float:
    """Colony-forming units per mL → cells per µm³ (one CFU counted as one cell)."""
    return c_cfu_per_ml / UM3_PER_ML


def per_um3_to_cfu_per_ml(c_per_um3: float) -> float:
    """Cells per µm³ → colony-forming units per mL (one cell counted as one CFU)."""
    return c_per_um3 * UM3_PER_ML


def kon_molar_to_um3(k_on_per_molar_s: float) -> float:
    """Association rate constant [1/(M·s)] → [µm³/s] per binding site.

    The rate per site is k_on · c; with c in mol/L the product has unit 1/s, with c in 1/µm³ the
    constant has unit µm³/s, so k[µm³/s] = k[1/(M·s)] / (N_A / 10¹⁵ µm³ per L).
    """
    return k_on_per_molar_s * UM3_PER_LITRE / AVOGADRO


def kon_um3_to_molar(k_on_um3_per_s: float) -> float:
    """Association rate constant [µm³/s] per binding site → [1/(M·s)]."""
    return k_on_um3_per_s * AVOGADRO / UM3_PER_LITRE
