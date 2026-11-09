"""Signal model: excitation, emission, collection, bleaching, background (stage Z6).

Per emitter of an ensemble sample (BioOptoMD ensembles Q per receptor state; not averaged
tensors, ADR draft section 3):

    R_exc = σ · I_eff / (h ν_ex),   I_eff = ½ n c ε₀ · 3 |μ̂_abs · E|²,
    detected rate per region d = R_exc · q(r, μ̂_em) · η_d(r, μ̂_em) · B,

with σ from the molar extinction coefficient, E from the field map F (per mW, scaled with the
power), q and η_d from the emitter response table (M) and B a brightness factor. For a freely
rotating dye 3 |μ̂·Ê|² averages to 1 (I_eff = I). Bleaching: each excitation destroys the dye
with probability φ_b, so N(t) = N₀ exp(−φ_b R_exc t) (linear regime, no saturation, no
triplet blinking).

Until hp-FEM and CliffordBem export real maps (Z2), :func:`placeholder_response_table` and
:func:`uniform_field_map` provide maps without any optics computed in BindOptiX: a collection
efficiency and a uniform irradiance given as parameters (P). The quenching of "off" aptamers is
then a brightness factor of the off ensemble (P); with real tables it follows from the heights of
the ensemble samples above gold.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import ArrayLike, NDArray

from bindoptix.ports import EmitterResponseTable, FieldMap, FieldMapStack, RegularGrid
from bindoptix.units import (
    SPEED_OF_LIGHT,
    VACUUM_PERMITTIVITY,
    extinction_to_cross_section_cm2,
    photon_energy_j,
)

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class Fluorophore:
    """Dye properties.

    Attributes:
        epsilon_per_molar_cm: molar extinction coefficient at the excitation wavelength.
        quantum_yield: intrinsic q₀ in the host medium.
        excitation_um, emission_um: vacuum wavelengths.
        bleach_yield: φ_b, destructions per excitation (P: 10⁻⁶).
    """

    epsilon_per_molar_cm: float
    quantum_yield: float
    excitation_um: float
    emission_um: float
    bleach_yield: float = 1e-6

    def __post_init__(self) -> None:
        if not (self.epsilon_per_molar_cm > 0 and 0 < self.quantum_yield <= 1):
            raise ValueError("epsilon > 0 and 0 < quantum_yield <= 1 required")
        if not (self.excitation_um > 0 and self.emission_um > 0 and self.bleach_yield >= 0):
            raise ValueError("wavelengths > 0 and bleach_yield >= 0 required")

    @property
    def cross_section_m2(self) -> float:
        return extinction_to_cross_section_cm2(self.epsilon_per_molar_cm) * 1e-4


@dataclass(frozen=True)
class EmitterEnsemble:
    """Samples of one receptor state relative to its anchor on the surface.

    Attributes:
        offsets_um: (N, 3) position relative to the anchor; z is the height above the surface.
        mu_abs: (N, 3) absorption dipole orientations (normalised on construction).
        mu_em: (N, 3) emission dipole orientations; default = mu_abs.
        weights: (N,) sample weights, normalised to sum 1.
        brightness: factor B (P: quenching of off states while tables are placeholders).
    """

    offsets_um: FloatArray
    mu_abs: FloatArray
    mu_em: FloatArray | None = None
    weights: FloatArray | None = None
    brightness: float = 1.0

    def __post_init__(self) -> None:
        off = np.atleast_2d(np.asarray(self.offsets_um, dtype=np.float64))
        n = off.shape[0]

        def unit(v: ArrayLike) -> FloatArray:
            a = np.atleast_2d(np.asarray(v, dtype=np.float64))
            if a.shape != (n, 3):
                raise ValueError("orientations must have shape (N, 3)")
            nrm = np.linalg.norm(a, axis=1, keepdims=True)
            if np.any(nrm == 0):
                raise ValueError("zero orientation")
            return np.asarray(a / nrm)

        if off.shape[1] != 3:
            raise ValueError("offsets_um must have shape (N, 3)")
        mu_a = unit(self.mu_abs)
        mu_e = mu_a if self.mu_em is None else unit(self.mu_em)
        w = np.full(n, 1.0 / n) if self.weights is None else np.asarray(self.weights, float)
        if w.shape != (n,) or np.any(w < 0) or w.sum() <= 0:
            raise ValueError("weights must be non-negative with one entry per sample")
        if not self.brightness >= 0.0:
            raise ValueError("brightness must be non-negative")
        object.__setattr__(self, "offsets_um", off)
        object.__setattr__(self, "mu_abs", mu_a)
        object.__setattr__(self, "mu_em", mu_e)
        object.__setattr__(self, "weights", w / w.sum())

    @classmethod
    def isotropic(
        cls, height_um: float, n: int, seed: int, brightness: float = 1.0
    ) -> EmitterEnsemble:
        """Random orientations at one height (placeholder for a BioOptoMD ensemble)."""
        g = np.random.default_rng(seed)
        z = g.uniform(-1, 1, n)
        phi = g.uniform(0, 2 * np.pi, n)
        r = np.sqrt(1 - z * z)
        mu = np.column_stack([r * np.cos(phi), r * np.sin(phi), z])
        off = np.column_stack([np.zeros(n), np.zeros(n), np.full(n, height_um)])
        return cls(off, mu, brightness=brightness)


def irradiance_w_m2(e_v_m: ArrayLike, n_host: float) -> FloatArray:
    """I = ½ n c ε₀ |E|² for complex fields (…, 3)."""
    e = np.asarray(e_v_m)
    return np.asarray(
        0.5 * n_host * SPEED_OF_LIGHT * VACUUM_PERMITTIVITY * np.sum(np.abs(e) ** 2, axis=-1)
    )


def excitation_rate_s(
    e_v_m: ArrayLike, mu_abs: ArrayLike, dye: Fluorophore, n_host: float
) -> FloatArray:
    """R_exc = σ · ½ n c ε₀ · 3|μ̂·E|² / (h ν) [1/s] per sample."""
    e = np.atleast_2d(np.asarray(e_v_m))
    mu = np.atleast_2d(np.asarray(mu_abs, dtype=np.float64))
    proj = np.abs(np.einsum("ni,ni->n", mu, e)) ** 2
    i_eff = 0.5 * n_host * SPEED_OF_LIGHT * VACUUM_PERMITTIVITY * 3.0 * proj
    return np.asarray(dye.cross_section_m2 * i_eff / photon_energy_j(dye.excitation_um))


def detected_rate_s(
    anchor_um: ArrayLike,
    ensemble: EmitterEnsemble,
    excitation: FieldMap | FieldMapStack,
    table: EmitterResponseTable,
    dye: Fluorophore,
    power_mw: float,
) -> FloatArray:
    """Mean detected photon rate per emitter of the ensemble at an anchor, per region (D,)."""
    anchor = np.asarray(anchor_um, dtype=np.float64).reshape(1, 3)
    pts = anchor + ensemble.offsets_um
    e = excitation.e_at(pts, power_mw)
    r_exc = excitation_rate_s(e, ensemble.mu_abs, dye, table.n_host)
    mu_em = ensemble.mu_em if ensemble.mu_em is not None else ensemble.mu_abs
    q = table.quantum_yield(pts, mu_em, dye.quantum_yield)
    eta = table.collection(pts, mu_em)
    w = ensemble.weights if ensemble.weights is not None else np.full(pts.shape[0], 1.0)
    per_sample = (r_exc * q)[:, None] * eta
    return np.asarray(ensemble.brightness * (w @ per_sample))


def surviving_fraction(rate_exc_s: float, dye: Fluorophore, t_s: float) -> float:
    """N(t)/N₀ = exp(−φ_b R_exc t)."""
    return math.exp(-dye.bleach_yield * rate_exc_s * t_s)


def photons_with_bleaching(
    rate_det_s: float, rate_exc_s: float, dye: Fluorophore, t0: float, t1: float
) -> float:
    """Detected photons per initial emitter between t0 and t1 with bleaching."""
    k = dye.bleach_yield * rate_exc_s
    if k == 0.0:
        return rate_det_s * (t1 - t0)
    return rate_det_s / k * (math.exp(-k * t0) - math.exp(-k * t1))


# --- placeholder maps (no optics computed here; parameters P) ---------------------------------


def placeholder_response_table(
    collection: float,
    emission_um: float,
    n_host: float = 1.333,
    quench_total: float = 1.0,
) -> EmitterResponseTable:
    """A_rad = I, A_tot = quench_total · I (≥ 1), A_det = collection · I, position invariant (P)."""
    if not (0.0 <= collection <= 1.0 and quench_total >= 1.0):
        raise ValueError("0 <= collection <= 1 and quench_total >= 1 required")
    grid = RegularGrid((np.zeros(1), np.zeros(1), np.zeros(1)))
    eye = np.eye(3).reshape(1, 1, 1, 3, 3)
    return EmitterResponseTable(
        grid,
        emission_um,
        n_host,
        ("camera",),
        quench_total * eye,
        eye.copy(),
        collection * eye[..., None, :, :],
        {"solver": "placeholder", "status": "P"},
    )


def uniform_field_map(
    irradiance_w_m2_per_mw: float,
    excitation_um: float,
    n_host: float = 1.333,
    polarization: ArrayLike = (1.0, 0.0, 0.0),
) -> FieldMap:
    """Uniform field with the given irradiance per mW incident power (P)."""
    pol = np.asarray(polarization, dtype=np.complex128)
    pol = pol / np.linalg.norm(pol)
    amp = math.sqrt(2.0 * irradiance_w_m2_per_mw / (n_host * SPEED_OF_LIGHT * VACUUM_PERMITTIVITY))
    grid = RegularGrid((np.zeros(1), np.zeros(1), np.zeros(1)))
    return FieldMap(
        grid,
        excitation_um,
        (amp * pol).reshape(1, 1, 1, 3),
        source={"solver": "placeholder", "status": "P"},
    )


# --- camera ------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Camera:
    """Camera referred to the object plane.

    Attributes:
        pixel_um: pixel pitch in the object plane.
        quantum_efficiency, read_noise_e, dark_e_s, exposure_s: sensor (P).
        psf_sigma_um: Gaussian width of the image of a point in the object plane (P; from the
            K5 chain later).
    """

    pixel_um: float
    quantum_efficiency: float
    read_noise_e: float
    dark_e_s: float
    exposure_s: float
    psf_sigma_um: float

    def __post_init__(self) -> None:
        if not (self.pixel_um > 0 and 0 < self.quantum_efficiency <= 1 and self.exposure_s > 0):
            raise ValueError("pixel_um > 0, 0 < QE <= 1, exposure > 0 required")
        if not (self.read_noise_e >= 0 and self.dark_e_s >= 0 and self.psf_sigma_um > 0):
            raise ValueError("noise parameters must be non-negative, psf_sigma_um > 0")

    @property
    def footprint_pixels(self) -> int:
        """Pixels of the detection aperture: square of side ≈ 2·FWHM of the PSF."""
        side = max(1, int(math.ceil(2 * 2.3548 * self.psf_sigma_um / self.pixel_um)))
        return side * side

    @property
    def aperture_fraction(self) -> float:
        """Share of a point image inside the square aperture (Gaussian PSF)."""
        side = math.sqrt(self.footprint_pixels) * self.pixel_um
        return math.erf(side / (2 * math.sqrt(2) * self.psf_sigma_um)) ** 2

    def sample(self, expected_photons: ArrayLike, rng: np.random.Generator) -> FloatArray:
        """Electrons per pixel: Poisson(QE·photons + dark·t) + Gaussian read noise."""
        lam = self.quantum_efficiency * np.asarray(expected_photons, dtype=np.float64)
        lam = lam + self.dark_e_s * self.exposure_s
        return np.asarray(rng.poisson(lam) + rng.normal(0.0, self.read_noise_e, lam.shape))


def tail_probability(mean_e: float, read_sigma_e: float, threshold_e: float) -> float:
    """P(X > threshold) for X = Poisson(mean) + Normal(0, read_sigma²), summed exactly."""
    from scipy.stats import norm, poisson

    if read_sigma_e == 0.0:
        return float(poisson.sf(math.floor(threshold_e), mean_e))
    lo = max(0, int(mean_e - 12 * math.sqrt(mean_e + 1) - 10))
    hi = int(mean_e + 12 * math.sqrt(mean_e + 1) + 10 + 6 * read_sigma_e)
    k = np.arange(lo, hi + 1)
    p = float(np.sum(poisson.pmf(k, mean_e) * norm.sf(threshold_e - k, scale=read_sigma_e)))
    return min(max(p, 0.0), 1.0)


@dataclass(frozen=True)
class SpotDetection:
    """Detection of single cells in one aperture against a uniform background."""

    signal_e: float
    background_e: float
    read_sigma_e: float
    threshold_e: float
    detection_probability: float
    false_alarms_per_mm2: float
    snr: float = field(default=0.0)


def spot_detection(
    camera: Camera,
    cell_photons: float,
    background_photons_per_um2: float,
    false_alarm_target_per_mm2: float,
) -> SpotDetection:
    """Threshold for a false-alarm density target, then the detection probability of a cell.

    ``cell_photons``: detected photons of one cell in the exposure (all of its image);
    ``background_photons_per_um2``: detected photons per object-plane µm² in the exposure.
    Apertures tile the image, so the false-alarm density is P_fa / aperture area.
    """
    npix = camera.footprint_pixels
    area_um2 = npix * camera.pixel_um**2
    qe = camera.quantum_efficiency
    bg = qe * background_photons_per_um2 * area_um2 + npix * camera.dark_e_s * camera.exposure_s
    sig = qe * cell_photons * camera.aperture_fraction
    rs = math.sqrt(npix) * camera.read_noise_e
    p_target = false_alarm_target_per_mm2 * area_um2 * 1e-6
    lo, hi = bg, bg + 50 * math.sqrt(bg + rs * rs + 1) + 50
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if tail_probability(bg, rs, mid) > p_target:
            lo = mid
        else:
            hi = mid
    thr = hi
    p_fa = tail_probability(bg, rs, thr)
    p_det = tail_probability(bg + sig, rs, thr)
    snr = sig / math.sqrt(sig + bg + rs * rs) if sig > 0 else 0.0
    return SpotDetection(sig, bg, rs, thr, p_det, p_fa / (area_um2 * 1e-6), snr)
