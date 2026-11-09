"""Binding kinetics, cell capture and receptor states (stage Z5).

- Competitive Langmuir kinetics of dissolved species on one kind of site,
      dΓ_i/dt = k_on,i c_i θ − k_off,i Γ_i,   θ = Γ_max − Σ_j α_j Γ_j,
  for constant concentrations a linear system; :class:`CompetitiveLangmuir` integrates it
  exactly with the matrix exponential and gives the closed-form equilibrium
  Γ_i = (k_on,i c_i / k_off,i) θ with θ = Γ_max / (1 + Σ_j α_j k_on,j c_j / k_off,j).
- Capture of whole cells (:class:`AvidityCapture`): the first aptamer–cell bond within the
  contact patch captures the cell irreversibly (many bonds follow within seconds, avidity);
  the rate is k_cap = k_bond · θ · A_contact and falls with fouling of the spot (P).
- Port B (:class:`Occupancy`, ADR draft ``docs/adr-entwurf-datenports.md``): shares of the
  receptor states per spot over time and the list of captured cells.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import scipy.linalg
from numpy.typing import ArrayLike, NDArray

from bindoptix.results import decode_array, encode_array, load_json, to_json

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class CompetitiveLangmuir:
    """Competitive Langmuir kinetics of several species on Γ_max sites [1/µm²].

    Attributes:
        sites_per_um2: Γ_max.
        k_on_um3_s, k_off_s, weights: per species (α: sites blocked per bound molecule).
    """

    sites_per_um2: float
    k_on_um3_s: tuple[float, ...]
    k_off_s: tuple[float, ...]
    weights: tuple[float, ...]

    def __post_init__(self) -> None:
        n = len(self.k_on_um3_s)
        if not (len(self.k_off_s) == n == len(self.weights) and n > 0):
            raise ValueError("one k_on, k_off and weight per species")
        if not self.sites_per_um2 >= 0.0:
            raise ValueError("sites_per_um2 must be non-negative")
        for v in (*self.k_on_um3_s, *self.k_off_s, *self.weights):
            if not (math.isfinite(v) and v >= 0.0):
                raise ValueError("rates and weights must be finite and non-negative")

    def _system(self, c: Sequence[float]) -> tuple[FloatArray, FloatArray]:
        kc = np.asarray(self.k_on_um3_s) * np.asarray(c, dtype=np.float64)
        a = -np.outer(kc, np.asarray(self.weights)) - np.diag(self.k_off_s)
        return a, kc * self.sites_per_um2

    def free_sites(self, gamma: ArrayLike) -> float:
        return float(self.sites_per_um2 - np.asarray(self.weights) @ np.asarray(gamma))

    def rhs(self, gamma: ArrayLike, c: Sequence[float]) -> FloatArray:
        a, b = self._system(c)
        return np.asarray(a @ np.asarray(gamma, dtype=np.float64) + b)

    def equilibrium(self, c: Sequence[float]) -> FloatArray:
        """Closed form; requires k_off > 0 for every species with k_on c > 0."""
        kc = np.asarray(self.k_on_um3_s) * np.asarray(c, dtype=np.float64)
        koff = np.asarray(self.k_off_s, dtype=np.float64)
        if np.any((kc > 0) & (koff == 0)):
            raise ValueError("irreversible binding has no equilibrium below saturation")
        ratio = np.divide(kc, koff, out=np.zeros_like(kc), where=kc > 0)
        theta = self.sites_per_um2 / (1.0 + np.asarray(self.weights) @ ratio)
        return np.asarray(ratio * theta)

    def integrate(self, gamma0: ArrayLike, c: Sequence[float], t_s: float) -> FloatArray:
        """Γ(t) for constant concentrations, exact (augmented matrix exponential)."""
        a, b = self._system(c)
        n = a.shape[0]
        m = np.zeros((n + 1, n + 1))
        m[:n, :n], m[:n, n] = a, b
        e = scipy.linalg.expm(m * t_s)
        return np.asarray(e[:n, :n] @ np.asarray(gamma0, dtype=np.float64) + e[:n, n])


@dataclass(frozen=True)
class AvidityCapture:
    """Irreversible capture of a cell by its first bond in the contact patch (P).

    Attributes:
        bond_rate_um2_s: k_bond, formation rate of one bond per free site in the patch [µm²/s
            per site, i.e. rate = k_bond · θ · A].
        contact_area_um2: A_contact of a cell lying on the spot (P: 0.5 µm² for a 0.5 µm rod).
    """

    bond_rate_um2_s: float
    contact_area_um2: float = 0.5

    def __post_init__(self) -> None:
        if not (self.bond_rate_um2_s >= 0.0 and self.contact_area_um2 > 0.0):
            raise ValueError("bond_rate_um2_s ≥ 0 and contact_area_um2 > 0 required")

    def capture_rate_s(self, free_sites_per_um2: ArrayLike) -> FloatArray:
        theta = np.maximum(np.asarray(free_sites_per_um2, dtype=np.float64), 0.0)
        return np.asarray(self.bond_rate_um2_s * theta * self.contact_area_um2)


@dataclass(frozen=True)
class Occupancy:
    """Port B: receptor states per spot over time and captured cells.

    Attributes:
        times_s: (nt,).
        spots: grid indices of the spots (S,).
        states: names of the receptor states (first one "free").
        fractions: (nt, S, n_states), each row sums to 1.
        cell_time_s, cell_xy_um, cell_spot, cell_strain: captured cells (sensor frame).
    """

    times_s: FloatArray
    spots: tuple[int, ...]
    states: tuple[str, ...]
    fractions: FloatArray
    cell_time_s: FloatArray = field(default_factory=lambda: np.zeros(0))
    cell_xy_um: FloatArray = field(default_factory=lambda: np.zeros((0, 2)))
    cell_spot: NDArray[np.int64] = field(default_factory=lambda: np.zeros(0, dtype=np.int64))
    cell_strain: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        t = np.asarray(self.times_s, dtype=np.float64)
        f = np.asarray(self.fractions, dtype=np.float64)
        if not self.states or self.states[0] != "free":
            raise ValueError('the first state must be "free"')
        if f.shape != (t.size, len(self.spots), len(self.states)):
            raise ValueError("fractions must have shape (nt, S, n_states)")
        if np.any(f < -1e-12) or np.any(f > 1 + 1e-12):
            raise ValueError("fractions must lie in [0, 1]")
        if np.max(np.abs(f.sum(axis=2) - 1.0), initial=0.0) > 1e-9:
            raise ValueError("fractions of each spot must sum to 1")
        n = np.asarray(self.cell_time_s).size
        if (
            np.asarray(self.cell_xy_um).shape != (n, 2)
            or np.asarray(self.cell_spot).shape != (n,)
            or len(self.cell_strain) != n
        ):
            raise ValueError("cell event arrays must have one entry per cell")
        object.__setattr__(self, "times_s", t)
        object.__setattr__(self, "fractions", f)

    def to_json(self) -> str:
        return to_json(
            "Occupancy",
            {
                "times_s": encode_array(self.times_s),
                "spots": list(self.spots),
                "states": list(self.states),
                "fractions": encode_array(self.fractions),
                "cell_time_s": encode_array(np.asarray(self.cell_time_s, dtype=np.float64)),
                "cell_xy_um": encode_array(np.asarray(self.cell_xy_um, dtype=np.float64)),
                "cell_spot": encode_array(np.asarray(self.cell_spot, dtype=np.int64)),
                "cell_strain": list(self.cell_strain),
            },
        )

    @classmethod
    def from_json(cls, text: str) -> Occupancy:
        d: dict[str, Any] = load_json(text, "Occupancy")
        return cls(
            times_s=decode_array(d["times_s"]),
            spots=tuple(int(s) for s in d["spots"]),
            states=tuple(str(s) for s in d["states"]),
            fractions=decode_array(d["fractions"]),
            cell_time_s=decode_array(d["cell_time_s"]),
            cell_xy_um=decode_array(d["cell_xy_um"]),
            cell_spot=decode_array(d["cell_spot"]),
            cell_strain=tuple(str(s) for s in d["cell_strain"]),
        )
