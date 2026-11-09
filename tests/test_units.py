"""Units and constants (docs/konventionen.md); reference values typed from the SI brochure."""

import math

import pytest

from bindoptix import units


def test_exact_si_constants() -> None:
    assert units.AVOGADRO == 6.02214076e23
    assert units.BOLTZMANN == 1.380649e-23
    assert units.STANDARD_GRAVITY_M_S2 == 9.80665


def test_one_nanomolar_in_particles_per_um3() -> None:
    # 1 nM = 1e-9 mol/L · 6.02214076e23 /mol / 1e15 µm³/L
    assert units.molar_to_per_um3(1e-9) == pytest.approx(0.602214076, rel=1e-15)


def test_one_cfu_per_ml_in_cells_per_um3() -> None:
    # 1 mL = 1 cm³ = (1e4 µm)³
    assert units.cfu_per_ml_to_per_um3(1.0) == pytest.approx(1.0 / (1e4) ** 3, rel=1e-15)


def test_regulatory_value_in_cells_per_um3() -> None:
    # TrinkwV technical action value 100 CFU / 100 mL
    assert units.cfu_per_ml_to_per_um3(100.0 / 100.0) == pytest.approx(1e-12, rel=1e-15)


def test_kon_conversion_gives_same_rate() -> None:
    # The binding rate per site k_on·c must not depend on the unit system.
    k_molar, c_molar = 1.0e6, 2.5e-9  # 1/(M·s), mol/L
    rate_molar = k_molar * c_molar
    rate_counts = units.kon_molar_to_um3(k_molar) * units.molar_to_per_um3(c_molar)
    assert rate_counts == pytest.approx(rate_molar, rel=1e-14)


@pytest.mark.parametrize("x", [1e-12, 3.7e-9, 1.0, 42.0])
def test_round_trips(x: float) -> None:
    assert units.per_um3_to_molar(units.molar_to_per_um3(x)) == pytest.approx(x, rel=1e-15)
    assert units.per_um3_to_cfu_per_ml(units.cfu_per_ml_to_per_um3(x)) == pytest.approx(
        x, rel=1e-15
    )
    assert units.kon_um3_to_molar(units.kon_molar_to_um3(x)) == pytest.approx(x, rel=1e-15)
    assert math.isclose(units.UM_PER_MM * units.MM_PER_UM, 1.0)
