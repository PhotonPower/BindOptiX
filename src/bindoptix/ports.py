"""Optical maps as data ports (stage Z2; ADR draft ``docs/adr-entwurf-datenports.md``).

BindOptiX computes no light. The maps come from the light solvers of the fusion (hp-FEM,
CliffordBem, the chain of pilot K5) and are read, checked, interpolated and evaluated here.

- :class:`EmitterResponseTable` (in M): per emitter position the real symmetric 3×3 tensors
  A_tot (total power), A_rad (radiated power) and A_det,d (power reaching detector region d), each
  normalised to the free dipole in the host medium at that position. For a unit orientation μ̂ the
  power is a = μ̂ᵀ A μ̂.
- :class:`FieldMap` (port F): complex E on a regular grid per 1 mW incident power, time factor
  e^(−iωt); optionally absorbed power density and temperature rise.

All positions are in the sensor frame (µm, docs/konventionen.md). Interpolation is trilinear;
it is a convex combination, so it keeps symmetric tensors symmetric and positive semidefinite.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from bindoptix.results import decode_array, encode_array, load_json, to_json

FloatArray = NDArray[np.float64]
ComplexArray = NDArray[np.complex128]

SYMMETRY_TOL = 1e-12
PSD_TOL = 1e-9


@dataclass(frozen=True)
class RegularGrid:
    """Axis-aligned grid in the sensor frame [µm]; an axis of length 1 means invariance along it.

    Attributes:
        axes_um: three strictly increasing coordinate arrays (x, y, z).
    """

    axes_um: tuple[ArrayLike, ArrayLike, ArrayLike]

    def __post_init__(self) -> None:
        axes = tuple(np.asarray(a, dtype=np.float64).ravel() for a in self.axes_um)
        if len(axes) != 3:
            raise ValueError("three axes required")
        for a in axes:
            if a.size < 1 or not np.isfinite(a).all():
                raise ValueError("axes must be non-empty and finite")
            if a.size > 1 and not np.all(np.diff(a) > 0.0):
                raise ValueError("axes must be strictly increasing")
        object.__setattr__(self, "axes_um", axes)

    @property
    def axes(self) -> tuple[FloatArray, FloatArray, FloatArray]:
        """The three axes as float arrays [µm]."""
        a = self.axes_um
        return (np.asarray(a[0]), np.asarray(a[1]), np.asarray(a[2]))

    @property
    def shape(self) -> tuple[int, int, int]:
        a = self.axes
        return (a[0].size, a[1].size, a[2].size)

    @property
    def cell_volume_um3(self) -> float:
        """Product of the mean spacings of the non-invariant axes (for nesting)."""
        v = 1.0
        for a in self.axes:
            if a.size > 1:
                v *= float(np.mean(np.diff(a)))
        return v

    def contains(self, points_um: ArrayLike) -> NDArray[np.bool_]:
        p = np.atleast_2d(np.asarray(points_um, dtype=np.float64))
        ok = np.ones(p.shape[0], dtype=bool)
        for k, a in enumerate(self.axes):
            if a.size > 1:
                ok &= (p[:, k] >= a[0]) & (p[:, k] <= a[-1])
        return ok

    def weights(self, points_um: ArrayLike) -> tuple[NDArray[np.intp], FloatArray]:
        """Flat corner indices (N, 8) and trilinear weights (N, 8) of the points."""
        p = np.atleast_2d(np.asarray(points_um, dtype=np.float64))
        if p.shape[-1] != 3:
            raise ValueError("points must have shape (N, 3)")
        if not self.contains(p).all():
            raise ValueError("points outside the grid")
        idx: list[NDArray[np.intp]] = []
        frac: list[FloatArray] = []
        for k, a in enumerate(self.axes):
            if a.size == 1:
                idx.append(np.zeros((p.shape[0], 2), dtype=np.intp))
                frac.append(np.zeros(p.shape[0]))
                continue
            i = np.clip(np.searchsorted(a, p[:, k], side="right") - 1, 0, a.size - 2)
            t = (p[:, k] - a[i]) / (a[i + 1] - a[i])
            idx.append(np.stack([i, i + 1], axis=1))
            frac.append(t)
        nx, ny, nz = self.shape
        corners, w = [], []
        for bx in (0, 1):
            for by in (0, 1):
                for bz in (0, 1):
                    ix, iy, iz = idx[0][:, bx], idx[1][:, by], idx[2][:, bz]
                    corners.append((ix * ny + iy) * nz + iz)
                    w.append(
                        (frac[0] if bx else 1.0 - frac[0])
                        * (frac[1] if by else 1.0 - frac[1])
                        * (frac[2] if bz else 1.0 - frac[2])
                    )
        return np.stack(corners, axis=1), np.stack(w, axis=1)

    def interpolate(self, values: np.ndarray[Any, Any], points_um: ArrayLike) -> Any:
        """Trilinear interpolation of ``values`` with shape (nx, ny, nz, …) at points (N, 3)."""
        if values.shape[:3] != self.shape:
            raise ValueError("values do not match the grid")
        flat = values.reshape((-1, *values.shape[3:]))
        c, w = self.weights(points_um)
        return np.einsum("nk,nk...->n...", w, flat[c])

    def to_data(self) -> dict[str, Any]:
        return {"axes_um": [encode_array(a) for a in self.axes]}

    @classmethod
    def from_data(cls, d: dict[str, Any]) -> RegularGrid:
        a = [decode_array(x) for x in d["axes_um"]]
        return cls((a[0], a[1], a[2]))


def quadratic(a: ArrayLike, mu: ArrayLike) -> FloatArray:
    """μ̂ᵀ A μ̂ for tensors (N, …, 3, 3) and orientations (N, 3) (unit vectors)."""
    m = np.asarray(mu, dtype=np.float64)
    t = np.asarray(a, dtype=np.float64)
    if t.ndim == 4:
        return np.asarray(np.einsum("ni,ndij,nj->nd", m, t, m), dtype=np.float64)
    return np.asarray(np.einsum("ni,nij,nj->n", m, t, m), dtype=np.float64)


def _check_psd(name: str, a: FloatArray) -> None:
    if not np.isfinite(a).all():
        raise ValueError(f"{name} contains non-finite values")
    scale = max(float(np.max(np.abs(a))), 1e-300)
    if np.max(np.abs(a - np.swapaxes(a, -1, -2))) > SYMMETRY_TOL * scale:
        raise ValueError(f"{name} is not symmetric")
    ev = np.linalg.eigvalsh(a)
    tr = np.abs(np.trace(a, axis1=-2, axis2=-1))[..., None]
    if np.any(ev < -PSD_TOL * np.maximum(tr, 1e-300)):
        raise ValueError(f"{name} is not positive semidefinite")


@dataclass(frozen=True)
class EmitterResponseTable:
    """Emitter response tensors on a grid (ADR draft, section 3).

    Attributes:
        grid: positions of the emitter [µm, sensor frame].
        wavelength_um: emission wavelength (vacuum) [µm].
        n_host: real refractive index of the host medium used for the normalisation.
        detectors: names of the detector regions d.
        a_tot, a_rad: shape (nx, ny, nz, 3, 3), normalised to the free dipole in the host.
        a_det: shape (nx, ny, nz, D, 3, 3), normalised the same way.
        source: provenance (solver, version, content hash, convergence).
    """

    grid: RegularGrid
    wavelength_um: float
    n_host: float
    detectors: tuple[str, ...]
    a_tot: FloatArray
    a_rad: FloatArray
    a_det: FloatArray
    source: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not (math.isfinite(self.wavelength_um) and self.wavelength_um > 0.0):
            raise ValueError("wavelength_um must be finite and positive")
        if not (math.isfinite(self.n_host) and self.n_host > 0.0):
            raise ValueError("n_host must be finite and positive")
        g = self.grid.shape
        tot = np.asarray(self.a_tot, dtype=np.float64)
        rad = np.asarray(self.a_rad, dtype=np.float64)
        det = np.asarray(self.a_det, dtype=np.float64)
        if tot.shape != (*g, 3, 3) or rad.shape != (*g, 3, 3):
            raise ValueError("a_tot and a_rad must have shape (nx, ny, nz, 3, 3)")
        if det.shape != (*g, len(self.detectors), 3, 3):
            raise ValueError("a_det must have shape (nx, ny, nz, D, 3, 3)")
        _check_psd("a_tot", tot)
        _check_psd("a_rad", rad)
        _check_psd("a_det", det)
        _check_psd("a_tot - a_rad (absorbed power)", tot - rad)
        _check_psd("a_rad - sum(a_det) (power missing the detectors)", rad - det.sum(axis=3))
        object.__setattr__(self, "a_tot", tot)
        object.__setattr__(self, "a_rad", rad)
        object.__setattr__(self, "a_det", det)

    def tensors_at(self, points_um: ArrayLike) -> tuple[FloatArray, FloatArray, FloatArray]:
        """Interpolated A_tot (N, 3, 3), A_rad (N, 3, 3), A_det (N, D, 3, 3)."""
        g = self.grid
        return (
            np.asarray(g.interpolate(self.a_tot, points_um), dtype=np.float64),
            np.asarray(g.interpolate(self.a_rad, points_um), dtype=np.float64),
            np.asarray(g.interpolate(self.a_det, points_um), dtype=np.float64),
        )

    def quantum_yield(self, points_um: ArrayLike, mu: ArrayLike, q0: float) -> FloatArray:
        """Modified quantum yield q = a_rad / (a_tot + (1 − q₀)/q₀) for intrinsic yield q₀.

        Rates normalised to the radiative rate of the free dipole; the intrinsic non-radiative
        rate is γ_nr⁰/γ_r⁰ = (1 − q₀)/q₀ (anger2006, novotny2012).
        """
        if not 0.0 < q0 <= 1.0:
            raise ValueError("q0 must lie in (0, 1]")
        tot, rad, _ = self.tensors_at(points_um)
        a_tot, a_rad = quadratic(tot, mu), quadratic(rad, mu)
        return np.asarray(a_rad / (a_tot + (1.0 - q0) / q0), dtype=np.float64)

    def collection(self, points_um: ArrayLike, mu: ArrayLike) -> FloatArray:
        """Share of the radiated power reaching each detector region, shape (N, D)."""
        _, rad, det = self.tensors_at(points_um)
        a_rad = quadratic(rad, mu)
        return np.asarray(quadratic(det, mu) / a_rad[:, None], dtype=np.float64)

    def to_json(self) -> str:
        return to_json(
            "EmitterResponseTable",
            {
                "grid": self.grid.to_data(),
                "wavelength_um": self.wavelength_um,
                "n_host": self.n_host,
                "detectors": list(self.detectors),
                "a_tot": encode_array(self.a_tot),
                "a_rad": encode_array(self.a_rad),
                "a_det": encode_array(self.a_det),
                "source": dict(self.source),
            },
        )

    @classmethod
    def from_json(cls, text: str) -> EmitterResponseTable:
        d = load_json(text, "EmitterResponseTable")
        return cls(
            grid=RegularGrid.from_data(d["grid"]),
            wavelength_um=float(d["wavelength_um"]),
            n_host=float(d["n_host"]),
            detectors=tuple(str(x) for x in d["detectors"]),
            a_tot=decode_array(d["a_tot"]),
            a_rad=decode_array(d["a_rad"]),
            a_det=decode_array(d["a_det"]),
            source={str(k): str(v) for k, v in d["source"].items()},
        )


@dataclass(frozen=True)
class FieldMap:
    """Port F: excitation field on a regular grid (ADR draft, section 1).

    Attributes:
        grid: positions [µm, sensor frame].
        wavelength_um: excitation wavelength (vacuum) [µm].
        e: complex field (nx, ny, nz, 3) [V/m] for 1 mW incident power, factor e^(−iωt).
        nonlinear: computed with temperature feedback; then valid only at ``power_mw``.
        power_mw: incident power the map was computed for (required if nonlinear).
        q_w_m3: absorbed power density per mW, shape (nx, ny, nz), optional.
        delta_t_k: temperature rise per mW, shape (nx, ny, nz), optional.
        source: provenance.
    """

    grid: RegularGrid
    wavelength_um: float
    e: ComplexArray
    nonlinear: bool = False
    power_mw: float | None = None
    q_w_m3: FloatArray | None = None
    delta_t_k: FloatArray | None = None
    source: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not (math.isfinite(self.wavelength_um) and self.wavelength_um > 0.0):
            raise ValueError("wavelength_um must be finite and positive")
        e = np.asarray(self.e, dtype=np.complex128)
        if e.shape != (*self.grid.shape, 3) or not np.isfinite(e).all():
            raise ValueError("e must be finite with shape (nx, ny, nz, 3)")
        object.__setattr__(self, "e", e)
        if self.nonlinear and not (self.power_mw is not None and self.power_mw > 0.0):
            raise ValueError("a nonlinear map needs the power it was computed for")
        for name in ("q_w_m3", "delta_t_k"):
            v = getattr(self, name)
            if v is not None:
                arr = np.asarray(v, dtype=np.float64)
                if arr.shape != self.grid.shape or not np.isfinite(arr).all():
                    raise ValueError(f"{name} must be finite with the grid shape")
                object.__setattr__(self, name, arr)

    def _scale(self, power_mw: float) -> float:
        if not (math.isfinite(power_mw) and power_mw >= 0.0):
            raise ValueError("power_mw must be finite and non-negative")
        if self.nonlinear and power_mw != self.power_mw:
            raise ValueError(f"nonlinear map valid only at {self.power_mw} mW")
        return power_mw

    def e_at(self, points_um: ArrayLike, power_mw: float = 1.0) -> ComplexArray:
        """Interpolated field (N, 3) at the given incident power [V/m]."""
        p = self._scale(power_mw)
        e = self.grid.interpolate(self.e, points_um)
        return np.asarray(e * math.sqrt(p), dtype=np.complex128)

    def excitation(self, points_um: ArrayLike, mu: ArrayLike, power_mw: float = 1.0) -> FloatArray:
        """|μ̂ · E|² at the points (N,) [V²/m²]; proportional to the excitation rate."""
        e = self.e_at(points_um, power_mw)
        m = np.asarray(mu, dtype=np.float64)
        return np.asarray(np.abs(np.einsum("ni,ni->n", m, e)) ** 2, dtype=np.float64)

    def to_json(self) -> str:
        return to_json(
            "FieldMap",
            {
                "grid": self.grid.to_data(),
                "wavelength_um": self.wavelength_um,
                "e": encode_array(self.e),
                "nonlinear": self.nonlinear,
                "power_mw": self.power_mw,
                "q_w_m3": None if self.q_w_m3 is None else encode_array(self.q_w_m3),
                "delta_t_k": None if self.delta_t_k is None else encode_array(self.delta_t_k),
                "source": dict(self.source),
            },
        )

    @classmethod
    def from_json(cls, text: str) -> FieldMap:
        d = load_json(text, "FieldMap")
        return cls(
            grid=RegularGrid.from_data(d["grid"]),
            wavelength_um=float(d["wavelength_um"]),
            e=decode_array(d["e"]),
            nonlinear=bool(d["nonlinear"]),
            power_mw=None if d["power_mw"] is None else float(d["power_mw"]),
            q_w_m3=None if d["q_w_m3"] is None else decode_array(d["q_w_m3"]),
            delta_t_k=None if d["delta_t_k"] is None else decode_array(d["delta_t_k"]),
            source={str(k): str(v) for k, v in d["source"].items()},
        )


@dataclass(frozen=True)
class FieldMapStack:
    """Nested field maps: each point is taken from the finest grid that contains it."""

    maps: tuple[FieldMap, ...]

    def __post_init__(self) -> None:
        if not self.maps:
            raise ValueError("at least one map required")
        lam = {m.wavelength_um for m in self.maps}
        if len(lam) != 1:
            raise ValueError("all maps must have the same wavelength")
        order = sorted(self.maps, key=lambda m: m.grid.cell_volume_um3)
        object.__setattr__(self, "maps", tuple(order))

    def select(self, points_um: ArrayLike) -> NDArray[np.intp]:
        """Index of the map used for each point (finest containing), −1 if none."""
        p = np.atleast_2d(np.asarray(points_um, dtype=np.float64))
        out = np.full(p.shape[0], -1, dtype=np.intp)
        for k, m in enumerate(self.maps):
            free = (out < 0) & m.grid.contains(p)
            out[free] = k
        return out

    def e_at(self, points_um: ArrayLike, power_mw: float = 1.0) -> ComplexArray:
        p = np.atleast_2d(np.asarray(points_um, dtype=np.float64))
        sel = self.select(p)
        if np.any(sel < 0):
            raise ValueError("points outside all maps")
        out = np.empty((p.shape[0], 3), dtype=np.complex128)
        for k, m in enumerate(self.maps):
            mask = sel == k
            if mask.any():
                out[mask] = m.e_at(p[mask], power_mw)
        return out


def orientation_average(a: ArrayLike) -> FloatArray:
    """Average of μ̂ᵀ A μ̂ over all orientations: tr(A)/3, for tensors (..., 3, 3)."""
    return np.asarray(np.trace(np.asarray(a), axis1=-2, axis2=-1) / 3.0, dtype=np.float64)


def unit_vectors(mu: Sequence[Sequence[float]] | ArrayLike) -> FloatArray:
    """Normalise orientations (N, 3) to unit length."""
    m = np.atleast_2d(np.asarray(mu, dtype=np.float64))
    n = np.linalg.norm(m, axis=1, keepdims=True)
    if np.any(n == 0.0):
        raise ValueError("zero orientation")
    return np.asarray(m / n, dtype=np.float64)
