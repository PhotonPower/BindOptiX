"""Background light without the analyte (stage Z6, part 2).

Two kinds of sources, both linear in the excitation power (one-photon excitation):

* Surface sources (photoluminescence of the gold, autofluorescence of the 2PP lacquer of the SIL
  and of the glass): detected photons per object-plane µm² and second per mW incident power,
  either over the whole field or on the spots only (gold under the spots only, D4 open). These
  numbers include the collection of the optics; they are what a measurement without sample
  gives (Z9) and stay parameters until Scattuccino or the K5 chain deliver them.
* Autofluorescence of cells: photons emitted per second into the host per W/m² local
  irradiance, as for a free emitter. The irradiance comes from the field map F, the collection
  from the emitter response table (M). With the rate model of the ADR draft (section 3), the
  detected rate of an emitter with orientation μ̂ is R_exc · γ_r⁰ μ̂ᵀA_det μ̂ / (γ_r⁰ μ̂ᵀA_tot μ̂
  + γ_nr⁰). Autofluorescence has a low intrinsic yield q₀; while q₀ · μ̂ᵀA_tot μ̂ ≪ 1 the
  denominator is γ_nr⁰ ≈ γ_r⁰/q₀ and the rate is R_emit · μ̂ᵀA_det μ̂ with R_emit = R_exc · q₀.
  Randomly oriented emitters in a cell average this to R_emit · tr(A_det)/3. Quenching close to
  the gold (large A_tot within a few 10 nm) is neglected; most of a cell is far from it.

Every source carries its status (L literature, K calibration, P placeholder; AGENTS.md rule 9);
:attr:`Background.placeholders` lists the placeholders that enter a result.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from bindoptix.ports import EmitterResponseTable, FieldMap, FieldMapStack, orientation_average
from bindoptix.signal import irradiance_w_m2

FloatArray = NDArray[np.float64]

REGIONS = ("field", "spots")
STATUSES = ("L", "K", "P")


def _check_status(status: str, reference: str) -> None:
    if status not in STATUSES:
        raise ValueError(f"status must be one of {STATUSES}")
    if status in ("L", "K") and not reference:
        raise ValueError("literature values and calibrations need a reference")


@dataclass(frozen=True)
class SurfaceBackground:
    """Areal background source in the object plane.

    Attributes:
        name: label, e.g. "gold PL".
        region: "field" (whole field of view) or "spots" (inside the spot polygons only).
        photons_per_um2_s_per_mw: detected photons per µm² and second per mW incident power
            (≥ 0, collection of the optics included).
        status: "L", "K" or "P".
        reference: key in docs/quellen.md or calibration record; required for L and K.
    """

    name: str
    region: str
    photons_per_um2_s_per_mw: float
    status: str
    reference: str = ""

    def __post_init__(self) -> None:
        if self.region not in REGIONS:
            raise ValueError(f"region must be one of {REGIONS}")
        v = self.photons_per_um2_s_per_mw
        if not (math.isfinite(v) and v >= 0.0):
            raise ValueError("photons_per_um2_s_per_mw must be finite and non-negative")
        _check_status(self.status, self.reference)


@dataclass(frozen=True)
class CellAutofluorescence:
    """Autofluorescence of one cell at the excitation wavelength.

    Attributes:
        photons_per_s_per_w_m2: photons emitted per second into the free host per W/m² local
            irradiance (≥ 0; absorption cross-section times intrinsic yield of all endogenous
            fluorophores of the cell).
        status: "L", "K" or "P".
        reference: key in docs/quellen.md or calibration record; required for L and K.
    """

    photons_per_s_per_w_m2: float
    status: str
    reference: str = ""

    def __post_init__(self) -> None:
        v = self.photons_per_s_per_w_m2
        if not (math.isfinite(v) and v >= 0.0):
            raise ValueError("photons_per_s_per_w_m2 must be finite and non-negative")
        _check_status(self.status, self.reference)


@dataclass(frozen=True)
class Background:
    """All background sources of a setup."""

    surfaces: tuple[SurfaceBackground, ...] = ()
    cells: CellAutofluorescence | None = None

    @property
    def placeholders(self) -> tuple[str, ...]:
        """Names of the sources with status P."""
        names = [s.name for s in self.surfaces if s.status == "P"]
        if self.cells is not None and self.cells.status == "P":
            names.append("cell autofluorescence")
        return tuple(names)

    def surface_photons_per_um2(self, power_mw: float, exposure_s: float) -> tuple[float, float]:
        """Detected photons per µm² in one exposure: (whole field, additionally on spots)."""
        dens = {r: 0.0 for r in REGIONS}
        for s in self.surfaces:
            dens[s.region] += s.photons_per_um2_s_per_mw * power_mw * exposure_s
        return dens["field"], dens["spots"]

    def render_terms(
        self, spot_polygons_um: tuple[FloatArray, ...], power_mw: float, exposure_s: float
    ) -> tuple[tuple[tuple[FloatArray, float], ...], float]:
        """Arguments ``areas`` and ``uniform_photons_per_um2`` of :func:`imaging.render`."""
        field, spots = self.surface_photons_per_um2(power_mw, exposure_s)
        areas = tuple((p, spots) for p in spot_polygons_um) if spots > 0.0 else ()
        return areas, field

    def cell_photons(
        self,
        points_um: ArrayLike,
        excitation: FieldMap | FieldMapStack,
        table: EmitterResponseTable,
        power_mw: float,
        exposure_s: float,
        detector: int = 0,
    ) -> FloatArray:
        """Detected autofluorescence photons per cell (N,) in one exposure at ``detector``.

        ``points_um``: cell centres (N, 3) in the sensor frame; the cell is treated as a point
        at its centre (the size of the cell against the variation of F and M is neglected).
        """
        pts = np.asarray(points_um, dtype=np.float64).reshape(-1, 3)
        if self.cells is None or pts.shape[0] == 0:
            return np.zeros(pts.shape[0])
        irr = irradiance_w_m2(excitation.e_at(pts, power_mw), table.n_host)
        _, _, det = table.tensors_at(pts)
        eta = orientation_average(det[:, detector])
        return np.asarray(self.cells.photons_per_s_per_w_m2 * irr * eta * exposure_s)
