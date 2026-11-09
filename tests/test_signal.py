"""Stage Z6: signal model (T17–T20, docs/testfaelle.md)."""

import math

import numpy as np
import pytest
from scipy.integrate import quad

from bindoptix.signal import (
    Camera,
    EmitterEnsemble,
    Fluorophore,
    detected_rate_s,
    excitation_rate_s,
    irradiance_w_m2,
    photons_with_bleaching,
    placeholder_response_table,
    spot_detection,
    surviving_fraction,
    tail_probability,
    uniform_field_map,
)
from bindoptix.units import extinction_to_cross_section_cm2

DYE = Fluorophore(250_000.0, 0.27, 0.640, 0.670, bleach_yield=1e-6)
SEED = 20261115


def three_axis_ensemble(height: float = 0.005, brightness: float = 1.0) -> EmitterEnsemble:
    """x, y, z orientations with equal weights: exact orientation average of quadratic forms."""
    return EmitterEnsemble(np.tile([0.0, 0.0, height], (3, 1)), np.eye(3), brightness=brightness)


def test_cross_section_literature_constant() -> None:
    # σ [cm²] = 3.82·10⁻²¹ · ε [1/(M cm)] (common textbook factor)
    assert extinction_to_cross_section_cm2(250_000.0) == pytest.approx(3.82e-21 * 250_000, rel=1e-3)


def test_uniform_field_irradiance() -> None:
    f = uniform_field_map(509.3, 0.640)
    e = f.e_at([[123.0, -4.0, 0.01]], 2.0)
    assert float(irradiance_w_m2(e, 1.333)[0]) == pytest.approx(2 * 509.3, rel=1e-12)


def test_t17_single_emitter_rate() -> None:
    # Freely rotating dye in the host medium: rate = σ I/(hν) · q₀ · η.
    irr, eta, p = 509.3, 0.3, 10.0
    f = uniform_field_map(irr, DYE.excitation_um)
    t = placeholder_response_table(eta, DYE.emission_um)
    r = detected_rate_s([0, 0, 0], three_axis_ensemble(), f, t, DYE, p)
    h, c = 6.62607015e-34, 299792458.0
    sigma = 3.82e-21 * 250_000 * 1e-4
    ref = sigma * irr * p / (h * c / 0.640e-6) * 0.27 * eta
    assert r[0] == pytest.approx(ref, rel=1e-3)


def test_t17_dipole_orientation_factor() -> None:
    # A dipole along the polarisation is excited three times as strongly as the average.
    f = uniform_field_map(100.0, 0.640, polarization=(1, 0, 0))
    e = f.e_at([[0, 0, 0]] * 3, 1.0)
    rates = excitation_rate_s(e, np.eye(3), DYE, 1.333)
    assert rates[1] == 0.0 and rates[2] == 0.0
    assert rates[0] == pytest.approx(3 * rates.mean(), rel=1e-14)


def test_t18_quenched_yield_from_tensors() -> None:
    # Placeholder quenching: A_tot = Q · I → q = q₀ / (Q q₀ + 1 − q₀) (rates in units of γ_r⁰).
    for big_q in (1.0, 5.0, 100.0):
        t = placeholder_response_table(0.3, DYE.emission_um, quench_total=big_q)
        q = t.quantum_yield([[0, 0, 0]], [[0, 0, 1]], DYE.quantum_yield)[0]
        assert q == pytest.approx(0.27 / (big_q * 0.27 + 1 - 0.27), rel=1e-14)


def test_t19_bleaching() -> None:
    r_exc, r_det = 2e4, 300.0
    k = DYE.bleach_yield * r_exc
    assert surviving_fraction(r_exc, DYE, 60.0) == pytest.approx(math.exp(-k * 60.0), rel=1e-14)
    ref = quad(lambda t: r_det * math.exp(-k * t), 10.0, 70.0)[0]
    assert photons_with_bleaching(r_det, r_exc, DYE, 10.0, 70.0) == pytest.approx(ref, rel=1e-10)


def test_t20_camera_statistics() -> None:
    cam = Camera(2.4, 0.8, 1.5, 0.2, 2.0, 1.0)
    g = np.random.default_rng(SEED)
    x = cam.sample(np.full(200_000, 50.0), g)
    mean = 0.8 * 50 + 0.4
    assert x.mean() == pytest.approx(mean, abs=4 * math.sqrt((mean + 2.25) / 2e5))
    assert x.var() == pytest.approx(mean + 2.25, rel=0.02)


@pytest.mark.parametrize(("mean", "sig", "thr"), [(20.0, 3.0, 30.0), (200.0, 6.0, 230.0)])
def test_t20_tail_probability_against_monte_carlo(mean: float, sig: float, thr: float) -> None:
    g = np.random.default_rng(SEED)
    n = 400_000
    x = g.poisson(mean, n) + g.normal(0, sig, n)
    emp = float(np.mean(x > thr))
    p = tail_probability(mean, sig, thr)
    assert p == pytest.approx(emp, abs=4 * math.sqrt(emp * (1 - emp) / n) + 1e-5)


def test_t20_detection_meets_false_alarm_target() -> None:
    cam = Camera(2.4, 0.8, 1.5, 0.2, 30.0, 1.2)
    d = spot_detection(cam, 960.0, 1.3 * 30.0, 0.05)
    assert d.false_alarms_per_mm2 <= 0.05 * (1 + 1e-6)
    assert d.false_alarms_per_mm2 > 0.01
    assert 0.0 < d.detection_probability <= 1.0
    weak = spot_detection(cam, 50.0, 1.3 * 30.0, 0.05)
    assert weak.detection_probability < d.detection_probability
