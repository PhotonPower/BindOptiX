"""Convection–diffusion of dissolved species with surface reactions (stage Z4).

Two-dimensional finite-volume model of a channel section in the plane (x along the flow,
z above the floor), per unit width: the flow profile u(z) is given (e.g. the width-averaged
profile of :class:`bindoptix.flow.DuctFlow`), lateral variations across the channel are not
resolved (part 2). Species i with diffusion D_i enter at x = 0 with concentration c_in,i(t)
[1/µm³]; reactive floor regions (spots) carry binding sites Γ_max [1/µm²] and react by
competitive Langmuir kinetics

    dΓ_i/dt = k_on,ij c_i,wall (Γ_max,j − Σ_k α_k Γ_k) − k_off,ij Γ_i

(squires2008 for the transport–reaction coupling). Discretisation: upwind convection, central
diffusion, implicit Euler; the free-site factor is taken from the previous step, so each species
is a linear system per step and the scheme is exactly conservative (T14). The wall concentration
follows from the half-cell diffusion resistance, c_wall = c₀ − r Δz₀/(2D), which turns into the
Dirichlet condition c_wall = 0 for k_on → ∞ (Lévêque limit, T12).

Units (docs/konventionen.md): µm, s, particles per µm³ and µm², k_on in µm³/s.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]


def graded_edges(length_um: float, n: int, first_um: float) -> FloatArray:
    """Edges 0 … length with n cells, first cell ``first_um``, geometric growth (bisection)."""
    if not (0.0 < first_um * n <= length_um * (1 + 1e-12)) and n > 1:
        raise ValueError("first cell too large for a graded grid")
    if n == 1 or abs(first_um * n - length_um) < 1e-12 * length_um:
        return np.asarray(np.linspace(0.0, length_um, n + 1), dtype=np.float64)
    lo, hi = 1.0 + 1e-12, 2.0
    while first_um * (hi**n - 1) / (hi - 1) < length_um:
        hi *= 2.0
    for _ in range(200):
        r = 0.5 * (lo + hi)
        if first_um * (r**n - 1) / (r - 1) < length_um:
            lo = r
        else:
            hi = r
    r = 0.5 * (lo + hi)
    e = np.concatenate([[0.0], np.cumsum(first_um * r ** np.arange(n))])
    e[-1] = length_um
    return e


@dataclass(frozen=True)
class Species:
    """Dissolved species.

    Attributes:
        name: label.
        diffusion_um2_s: D.
        site_weight: α, sites blocked by one bound molecule (1 for 1:1 binding).
    """

    name: str
    diffusion_um2_s: float
    site_weight: float = 1.0

    def __post_init__(self) -> None:
        if not (math.isfinite(self.diffusion_um2_s) and self.diffusion_um2_s > 0.0):
            raise ValueError("diffusion_um2_s must be finite and positive")
        if not self.site_weight >= 0.0:
            raise ValueError("site_weight must be non-negative")


@dataclass(frozen=True)
class ReactiveRegion:
    """Floor interval [x0, x1] with binding sites and rate constants per species.

    Attributes:
        x0_um, x1_um: extent along the flow (should coincide with grid edges).
        sites_per_um2: Γ_max.
        k_on_um3_s: one rate per species (0: species does not bind here).
        k_off_s: one rate per species.
    """

    x0_um: float
    x1_um: float
    sites_per_um2: float
    k_on_um3_s: tuple[float, ...]
    k_off_s: tuple[float, ...]

    def __post_init__(self) -> None:
        if not self.x1_um > self.x0_um:
            raise ValueError("x1_um must exceed x0_um")
        if not self.sites_per_um2 >= 0.0:
            raise ValueError("sites_per_um2 must be non-negative")
        if len(self.k_on_um3_s) != len(self.k_off_s):
            raise ValueError("k_on and k_off need one entry per species")
        if any(k < 0.0 for k in (*self.k_on_um3_s, *self.k_off_s)):
            raise ValueError("rate constants must be non-negative")


@dataclass
class ContinuumState:
    """Concentrations (S, nx, nz) [1/µm³], coverages (S, nx) [1/µm²], time and mass ledger."""

    c: FloatArray
    gamma: FloatArray
    t_s: float = 0.0
    inflow: FloatArray = field(default_factory=lambda: np.zeros(0))
    outflow: FloatArray = field(default_factory=lambda: np.zeros(0))


class ChannelSection:
    """2D finite-volume channel section with reactive floor regions.

    Args:
        x_edges_um, z_edges_um: cell edges (strictly increasing, starting at 0).
        velocity: u(z) [µm/s] at the z cell centres, u ≥ 0.
        species: dissolved species.
        regions: reactive floor regions; a floor cell is reactive if its centre lies inside.
    """

    def __init__(
        self,
        x_edges_um: ArrayLike,
        z_edges_um: ArrayLike,
        velocity: Callable[[FloatArray], FloatArray],
        species: Sequence[Species],
        regions: Sequence[ReactiveRegion] = (),
    ) -> None:
        self.xe = np.asarray(x_edges_um, dtype=np.float64)
        self.ze = np.asarray(z_edges_um, dtype=np.float64)
        for e in (self.xe, self.ze):
            if e.ndim != 1 or e.size < 2 or not np.all(np.diff(e) > 0.0):
                raise ValueError("edges must be strictly increasing")
        self.dx = np.diff(self.xe)
        self.dz = np.diff(self.ze)
        self.xc = 0.5 * (self.xe[1:] + self.xe[:-1])
        self.zc = 0.5 * (self.ze[1:] + self.ze[:-1])
        self.nx, self.nz = self.dx.size, self.dz.size
        self.u = np.asarray(velocity(self.zc), dtype=np.float64)
        if self.u.shape != (self.nz,) or np.any(self.u < 0.0):
            raise ValueError("velocity must return u(z) >= 0 at the z cell centres")
        self.species = tuple(species)
        self.regions = tuple(regions)
        ns = len(self.species)
        for r in self.regions:
            if len(r.k_on_um3_s) != ns:
                raise ValueError("each region needs one rate per species")
        self.region_of = np.full(self.nx, -1, dtype=np.int64)
        for j, r in enumerate(self.regions):
            self.region_of[(self.xc > r.x0_um) & (self.xc < r.x1_um)] = j
        self._base = [self._transport_matrix(s.diffusion_um2_s) for s in self.species]

    def _idx(self, i: NDArray[np.int64], k: NDArray[np.int64]) -> NDArray[np.int64]:
        return np.asarray(i * self.nz + k, dtype=np.int64)

    def _transport_matrix(self, d: float) -> sp.csr_matrix:
        """Operator L with dV·dc/dt = L c (per unit width), no reaction, no inflow term."""
        nx, nz = self.nx, self.nz
        rows: list[NDArray[np.int64]] = []
        cols: list[NDArray[np.int64]] = []
        vals: list[FloatArray] = []

        def add(r: NDArray[np.int64], c: NDArray[np.int64], v: FloatArray) -> None:
            rows.append(r)
            cols.append(c)
            vals.append(v)

        ii, kk = np.meshgrid(np.arange(nx), np.arange(nz), indexing="ij")
        # x faces between i and i+1
        i0, k0 = ii[:-1].ravel(), kk[:-1].ravel()
        a, b = self._idx(i0, k0), self._idx(i0 + 1, k0)
        gx = d * self.dz[k0] / (self.xc[i0 + 1] - self.xc[i0])
        fx = self.u[k0] * self.dz[k0]  # upwind, u ≥ 0: carries c_a
        add(a, a, -gx - fx)
        add(a, b, gx)
        add(b, b, -gx)
        add(b, a, gx + fx)
        # outlet: convective outflow of the last cells
        il = np.full(nz, nx - 1)
        kl = np.arange(nz)
        add(self._idx(il, kl), self._idx(il, kl), -self.u * self.dz)
        # z faces between k and k+1
        i1, k1 = ii[:, :-1].ravel(), kk[:, :-1].ravel()
        a, b = self._idx(i1, k1), self._idx(i1, k1 + 1)
        gz = d * self.dx[i1] / (self.zc[k1 + 1] - self.zc[k1])
        add(a, a, -gz)
        add(a, b, gz)
        add(b, b, -gz)
        add(b, a, gz)
        n = nx * nz
        return sp.csr_matrix(
            (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)
        )

    @property
    def cell_area_um2(self) -> FloatArray:
        """Cell areas (nx, nz) per unit width."""
        return np.asarray(np.outer(self.dx, self.dz))

    def initial_state(self) -> ContinuumState:
        ns = len(self.species)
        return ContinuumState(
            np.zeros((ns, self.nx, self.nz)),
            np.zeros((ns, self.nx)),
            0.0,
            np.zeros(ns),
            np.zeros(ns),
        )

    def _free_sites(self, gamma: FloatArray) -> FloatArray:
        theta = np.zeros(self.nx)
        for j, r in enumerate(self.regions):
            m = self.region_of == j
            used = sum(s.site_weight * gamma[q, m] for q, s in enumerate(self.species))
            theta[m] = np.maximum(r.sites_per_um2 - used, 0.0)
        return theta

    def _rates(self, q: int, theta: FloatArray) -> tuple[FloatArray, FloatArray]:
        """β and γ per floor cell for species q (zero off the regions)."""
        d = self.species[q].diffusion_um2_s
        beta = np.zeros(self.nx)
        gam = np.zeros(self.nx)
        half = 0.5 * self.dz[0] / d
        for j, r in enumerate(self.regions):
            m = self.region_of == j
            kon, koff = r.k_on_um3_s[q], r.k_off_s[q]
            if math.isinf(kon):
                beta[m] = 1.0 / half if r.sites_per_um2 > 0 else 0.0
                gam[m] = 0.0
            else:
                den = 1.0 + kon * theta[m] * half
                beta[m] = kon * theta[m] / den
                gam[m] = koff / den
        return beta, gam

    def step(self, state: ContinuumState, dt_s: float, c_in: Sequence[float]) -> ContinuumState:
        """One implicit Euler step with inlet concentrations ``c_in`` (one per species)."""
        if len(c_in) != len(self.species):
            raise ValueError("one inlet concentration per species")
        theta = self._free_sites(state.gamma)
        vol = self.cell_area_um2.ravel()
        floor = self._idx(np.arange(self.nx), np.zeros(self.nx, dtype=np.int64))
        inlet = self._idx(np.zeros(self.nz, dtype=np.int64), np.arange(self.nz))
        c_new = np.empty_like(state.c)
        g_new = np.empty_like(state.gamma)
        inflow, outflow = state.inflow.copy(), state.outflow.copy()
        for q in range(len(self.species)):
            beta, gam = self._rates(q, theta)
            sink = beta / (1.0 + dt_s * gam) * self.dx
            src = gam * state.gamma[q] / (1.0 + dt_s * gam) * self.dx
            a = (
                sp.diags(vol / dt_s)
                - self._base[q]
                + sp.csr_matrix((sink, (floor, floor)), shape=self._base[q].shape)
            )
            rhs = vol / dt_s * state.c[q].ravel()
            rhs[floor] += src
            rhs[inlet] += self.u * self.dz * c_in[q]
            c = spla.spsolve(a.tocsc(), rhs)
            c_floor = c[floor]
            gq = (state.gamma[q] + dt_s * beta * c_floor) / (1.0 + dt_s * gam)
            c_new[q] = c.reshape(self.nx, self.nz)
            g_new[q] = gq
            inflow[q] += dt_s * float(np.sum(self.u * self.dz)) * c_in[q]
            outflow[q] += dt_s * float(np.sum(self.u * self.dz * c_new[q, -1, :]))
        return ContinuumState(c_new, g_new, state.t_s + dt_s, inflow, outflow)

    def mass(self, state: ContinuumState) -> tuple[FloatArray, FloatArray]:
        """Dissolved and bound amount per species (per unit width) [1/µm]."""
        dissolved = np.einsum("sxz,xz->s", state.c, self.cell_area_um2)
        bound = state.gamma @ self.dx
        return dissolved, bound

    def steady_first_order(self, q: int, k_surface_um_s: ArrayLike, c_in: float) -> FloatArray:
        """Steady concentration (nx, nz) for a linear surface sink k·c_wall on the regions.

        ``k_surface_um_s`` per region (k_on Γ_max for Γ ≪ Γ_max; ∞ for c_wall = 0).
        """
        d = self.species[q].diffusion_um2_s
        half = 0.5 * self.dz[0] / d
        ks = np.broadcast_to(np.asarray(k_surface_um_s, dtype=np.float64), (len(self.regions),))
        beta = np.zeros(self.nx)
        for j in range(len(self.regions)):
            m = self.region_of == j
            beta[m] = 1.0 / half if math.isinf(ks[j]) else ks[j] / (1.0 + ks[j] * half)
        floor = self._idx(np.arange(self.nx), np.zeros(self.nx, dtype=np.int64))
        inlet = self._idx(np.zeros(self.nz, dtype=np.int64), np.arange(self.nz))
        a = -self._base[q] + sp.csr_matrix(
            (beta * self.dx, (floor, floor)), shape=self._base[q].shape
        )
        rhs = np.zeros(self.nx * self.nz)
        rhs[inlet] = self.u * self.dz * c_in
        c = spla.spsolve(a.tocsc(), rhs).reshape(self.nx, self.nz)
        self._last_beta = beta
        return np.asarray(c, dtype=np.float64)

    def surface_flux(self, c: FloatArray) -> FloatArray:
        """Flux into the floor per cell (nx,) [1/(µm·s)] of the last steady solution."""
        return np.asarray(self._last_beta * c[:, 0] * self.dx, dtype=np.float64)
