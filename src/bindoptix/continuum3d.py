"""Convection–diffusion with surface reactions in the full channel volume (stage Z4, part 2).

Three-dimensional version of :mod:`bindoptix.continuum`: channel frame x along the flow,
−w/2 ≤ y ≤ w/2 across, 0 ≤ z ≤ h above the floor; flow u(y, z) ≥ 0 (e.g.
:meth:`bindoptix.flow.DuctFlow.velocity_um_s`); side walls and ceiling are impermeable. Reactive
floor regions are convex polygons (e.g. the clipped spots of :mod:`bindoptix.chipsim`); a floor
cell reacts if its centre lies inside. Kinetics, wall concentration and time stepping as in
:mod:`bindoptix.continuum` (competitive Langmuir, half-cell resistance, implicit Euler with the
free sites of the previous step, exactly conservative).

Each step solves one linear system per species with BiCGSTAB, preconditioned by the LU
factorisation of the transport operator with the reaction sink of the first step; the sink only
changes on the floor cells, so few iterations suffice.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from numpy.typing import ArrayLike, NDArray

from bindoptix.continuum import ContinuumState, Species

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class FloorRegion:
    """Convex counter-clockwise floor polygon (n, 2) [µm] with sites and rates per species."""

    polygon_um: FloatArray
    sites_per_um2: float
    k_on_um3_s: tuple[float, ...]
    k_off_s: tuple[float, ...]

    def __post_init__(self) -> None:
        p = np.asarray(self.polygon_um, dtype=np.float64)
        if p.ndim != 2 or p.shape[1] != 2 or p.shape[0] < 3:
            raise ValueError("polygon_um must have shape (n >= 3, 2)")
        x, y = p[:, 0], p[:, 1]
        if 0.5 * (x @ np.roll(y, -1) - y @ np.roll(x, -1)) <= 0.0:
            raise ValueError("polygon_um must be counter-clockwise")
        if not self.sites_per_um2 >= 0.0:
            raise ValueError("sites_per_um2 must be non-negative")
        if len(self.k_on_um3_s) != len(self.k_off_s):
            raise ValueError("k_on and k_off need one entry per species")
        if any(k < 0.0 for k in (*self.k_on_um3_s, *self.k_off_s)):
            raise ValueError("rate constants must be non-negative")
        object.__setattr__(self, "polygon_um", p)

    def contains(self, x: FloatArray, y: FloatArray) -> NDArray[np.bool_]:
        p = self.polygon_um
        inside = np.ones(x.shape, dtype=bool)
        for i in range(p.shape[0]):
            a, b = p[i], p[(i + 1) % p.shape[0]]
            inside &= (b[0] - a[0]) * (y - a[1]) - (b[1] - a[1]) * (x - a[0]) >= 0.0
        return inside


class ChannelVolume:
    """3D finite-volume channel with reactive floor polygons.

    Args:
        x_edges_um, y_edges_um, z_edges_um: cell edges (y from −w/2 to w/2, z from 0 to h).
        velocity: u(y, z) [µm/s] for arrays of y and z cell centres (same shape), u ≥ 0.
        species: dissolved species.
        regions: reactive floor polygons.
    """

    def __init__(
        self,
        x_edges_um: ArrayLike,
        y_edges_um: ArrayLike,
        z_edges_um: ArrayLike,
        velocity: Callable[[FloatArray, FloatArray], FloatArray],
        species: Sequence[Species],
        regions: Sequence[FloorRegion] = (),
    ) -> None:
        self.xe, self.ye, self.ze = (
            np.asarray(e, dtype=np.float64) for e in (x_edges_um, y_edges_um, z_edges_um)
        )
        for e in (self.xe, self.ye, self.ze):
            if e.ndim != 1 or e.size < 2 or not np.all(np.diff(e) > 0.0):
                raise ValueError("edges must be strictly increasing")
        self.dx, self.dy, self.dz = np.diff(self.xe), np.diff(self.ye), np.diff(self.ze)
        self.xc = 0.5 * (self.xe[1:] + self.xe[:-1])
        self.yc = 0.5 * (self.ye[1:] + self.ye[:-1])
        self.zc = 0.5 * (self.ze[1:] + self.ze[:-1])
        self.nx, self.ny, self.nz = self.dx.size, self.dy.size, self.dz.size
        yy, zz = np.meshgrid(self.yc, self.zc, indexing="ij")
        self.u = np.asarray(velocity(yy.ravel(), zz.ravel()), dtype=np.float64).reshape(
            self.ny, self.nz
        )
        if np.any(self.u < 0.0) or not np.isfinite(self.u).all():
            raise ValueError("velocity must be finite and non-negative")
        self.species = tuple(species)
        self.regions = tuple(regions)
        ns = len(self.species)
        if any(len(r.k_on_um3_s) != ns for r in self.regions):
            raise ValueError("each region needs one rate per species")
        fx, fy = np.meshgrid(self.xc, self.yc, indexing="ij")
        self.region_of = np.full((self.nx, self.ny), -1, dtype=np.int64)
        for j, r in enumerate(self.regions):
            self.region_of[(self.region_of < 0) & r.contains(fx, fy)] = j
        self.floor_area = np.outer(self.dx, self.dy)
        self._floor = self._idx(
            *np.meshgrid(np.arange(self.nx), np.arange(self.ny), 0, indexing="ij")
        )
        self._floor = self._floor.reshape(self.nx, self.ny)
        self._base = [self._transport_matrix(s.diffusion_um2_s) for s in self.species]
        self._precond: dict[tuple[int, float], spla.SuperLU] = {}

    def _idx(self, i: ArrayLike, j: ArrayLike, k: ArrayLike) -> NDArray[np.int64]:
        return np.asarray(
            (np.asarray(i) * self.ny + np.asarray(j)) * self.nz + np.asarray(k), dtype=np.int64
        )

    @property
    def n_cells(self) -> int:
        return self.nx * self.ny * self.nz

    @property
    def cell_volume_um3(self) -> FloatArray:
        return np.asarray(self.dx[:, None, None] * self.dy[None, :, None] * self.dz[None, None, :])

    def _transport_matrix(self, d: float) -> sp.csr_matrix:
        nx, ny, nz = self.nx, self.ny, self.nz
        rows: list[NDArray[np.int64]] = []
        cols: list[NDArray[np.int64]] = []
        vals: list[FloatArray] = []
        ii, jj, kk = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")

        def couple(a: NDArray[np.int64], b: NDArray[np.int64], g: FloatArray) -> None:
            rows.extend([a, a, b, b])
            cols.extend([a, b, b, a])
            vals.extend([-g, g, -g, g])

        # x faces: diffusion and upwind convection (u ≥ 0)
        i0, j0, k0 = ii[:-1].ravel(), jj[:-1].ravel(), kk[:-1].ravel()
        a, b = self._idx(i0, j0, k0), self._idx(i0 + 1, j0, k0)
        area = self.dy[j0] * self.dz[k0]
        couple(a, b, d * area / (self.xc[i0 + 1] - self.xc[i0]))
        f = self.u[j0, k0] * area
        rows.extend([a, b])
        cols.extend([a, a])
        vals.extend([-f, f])
        jl, kl = np.meshgrid(np.arange(ny), np.arange(nz), indexing="ij")
        last = self._idx(nx - 1, jl.ravel(), kl.ravel())
        rows.append(last)
        cols.append(last)
        vals.append(-self.u.ravel() * np.outer(self.dy, self.dz).ravel())
        # y faces
        i1, j1, k1 = ii[:, :-1].ravel(), jj[:, :-1].ravel(), kk[:, :-1].ravel()
        couple(
            self._idx(i1, j1, k1),
            self._idx(i1, j1 + 1, k1),
            d * self.dx[i1] * self.dz[k1] / (self.yc[j1 + 1] - self.yc[j1]),
        )
        # z faces
        i2, j2, k2 = ii[:, :, :-1].ravel(), jj[:, :, :-1].ravel(), kk[:, :, :-1].ravel()
        couple(
            self._idx(i2, j2, k2),
            self._idx(i2, j2, k2 + 1),
            d * self.dx[i2] * self.dy[j2] / (self.zc[k2 + 1] - self.zc[k2]),
        )
        n = self.n_cells
        return sp.csr_matrix(
            (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)
        )

    def initial_state(self) -> ContinuumState:
        ns = len(self.species)
        return ContinuumState(
            np.zeros((ns, self.nx, self.ny, self.nz)),
            np.zeros((ns, self.nx, self.ny)),
            0.0,
            np.zeros(ns),
            np.zeros(ns),
        )

    def _free_sites(self, gamma: FloatArray) -> FloatArray:
        theta = np.zeros((self.nx, self.ny))
        for j, r in enumerate(self.regions):
            m = self.region_of == j
            used = sum(s.site_weight * gamma[q][m] for q, s in enumerate(self.species))
            theta[m] = np.maximum(r.sites_per_um2 - used, 0.0)
        return theta

    def _rates(self, q: int, theta: FloatArray) -> tuple[FloatArray, FloatArray]:
        half = 0.5 * self.dz[0] / self.species[q].diffusion_um2_s
        beta = np.zeros((self.nx, self.ny))
        gam = np.zeros((self.nx, self.ny))
        for j, r in enumerate(self.regions):
            m = self.region_of == j
            kon, koff = r.k_on_um3_s[q], r.k_off_s[q]
            if math.isinf(kon):
                beta[m] = 1.0 / half if r.sites_per_um2 > 0 else 0.0
            else:
                den = 1.0 + kon * theta[m] * half
                beta[m] = kon * theta[m] / den
                gam[m] = koff / den
        return beta, gam

    def _solve(
        self, q: int, dt_s: float, a: sp.csr_matrix, rhs: FloatArray, x0: FloatArray
    ) -> FloatArray:
        key = (q, dt_s)
        if key not in self._precond:
            self._precond[key] = spla.splu(a.tocsc())
            return np.asarray(self._precond[key].solve(rhs), dtype=np.float64)
        lu = self._precond[key]
        m = spla.LinearOperator(a.shape, matvec=lambda v: lu.solve(v))
        x, info = spla.bicgstab(a, rhs, x0=x0, M=m, rtol=1e-13, atol=0.0, maxiter=200)
        if info != 0:
            x = spla.spsolve(a.tocsc(), rhs)
        return np.asarray(x, dtype=np.float64)

    def step(
        self, state: ContinuumState, dt_s: float, c_in: Sequence[float], max_halvings: int = 30
    ) -> ContinuumState:
        """Implicit Euler over ``dt_s`` with inlet concentrations ``c_in`` (one per species).

        The free sites are taken from the start of a step, so a long step can bind more than
        the free sites; such a step is repeated as two half steps (each conservative) until no
        site capacity is exceeded.
        """
        new = self._step_once(state, dt_s, c_in)
        if max_halvings <= 0 or not self._overfilled(new):
            return new
        half = self.step(state, 0.5 * dt_s, c_in, max_halvings - 1)
        return self.step(half, 0.5 * dt_s, c_in, max_halvings - 1)

    def _overfilled(self, state: ContinuumState) -> bool:
        for j, r in enumerate(self.regions):
            m = self.region_of == j
            used = sum(s.site_weight * state.gamma[q][m] for q, s in enumerate(self.species))
            if np.any(used > r.sites_per_um2 * (1.0 + 1e-9) + 1e-300):
                return True
        return False

    def _step_once(
        self, state: ContinuumState, dt_s: float, c_in: Sequence[float]
    ) -> ContinuumState:
        if len(c_in) != len(self.species):
            raise ValueError("one inlet concentration per species")
        theta = self._free_sites(state.gamma)
        vol = self.cell_volume_um3.ravel()
        inlet_j, inlet_k = np.meshgrid(np.arange(self.ny), np.arange(self.nz), indexing="ij")
        inlet = self._idx(0, inlet_j.ravel(), inlet_k.ravel())
        inflow_rate = (self.u * np.outer(self.dy, self.dz)).ravel()
        floor = self._floor.ravel()
        c_new = np.empty_like(state.c)
        g_new = np.empty_like(state.gamma)
        inflow, outflow = state.inflow.copy(), state.outflow.copy()
        for q in range(len(self.species)):
            beta, gam = self._rates(q, theta)
            sink = (beta / (1.0 + dt_s * gam) * self.floor_area).ravel()
            src = (gam * state.gamma[q] / (1.0 + dt_s * gam) * self.floor_area).ravel()
            a = (
                sp.diags(vol / dt_s)
                - self._base[q]
                + sp.csr_matrix((sink, (floor, floor)), shape=self._base[q].shape)
            ).tocsr()
            rhs = vol / dt_s * state.c[q].ravel()
            rhs[floor] += src
            rhs[inlet] += inflow_rate * c_in[q]
            c = self._solve(q, dt_s, a, rhs, state.c[q].ravel())
            cq = c.reshape(self.nx, self.ny, self.nz)
            c_new[q] = cq
            g_new[q] = (state.gamma[q] + dt_s * beta * cq[:, :, 0]) / (1.0 + dt_s * gam)
            inflow[q] += dt_s * float(inflow_rate.sum()) * c_in[q]
            outflow[q] += dt_s * float(np.sum(self.u * np.outer(self.dy, self.dz) * cq[-1]))
        return ContinuumState(c_new, g_new, state.t_s + dt_s, inflow, outflow)

    def mass(self, state: ContinuumState) -> tuple[FloatArray, FloatArray]:
        """Dissolved and bound amount per species (particles)."""
        dissolved = np.einsum("sxyz,xyz->s", state.c, self.cell_volume_um3)
        bound = np.einsum("sxy,xy->s", state.gamma, self.floor_area)
        return dissolved, bound

    def steady_first_order(self, q: int, k_surface_um_s: ArrayLike, c_in: float) -> FloatArray:
        """Steady concentration (nx, ny, nz) for a linear surface sink k·c_wall per region."""
        half = 0.5 * self.dz[0] / self.species[q].diffusion_um2_s
        ks = np.broadcast_to(np.asarray(k_surface_um_s, dtype=np.float64), (len(self.regions),))
        beta = np.zeros((self.nx, self.ny))
        for j in range(len(self.regions)):
            m = self.region_of == j
            beta[m] = 1.0 / half if math.isinf(ks[j]) else ks[j] / (1.0 + ks[j] * half)
        floor = self._floor.ravel()
        a = -self._base[q] + sp.csr_matrix(
            ((beta * self.floor_area).ravel(), (floor, floor)), shape=self._base[q].shape
        )
        inlet_j, inlet_k = np.meshgrid(np.arange(self.ny), np.arange(self.nz), indexing="ij")
        rhs = np.zeros(self.n_cells)
        rhs[self._idx(0, inlet_j.ravel(), inlet_k.ravel())] = (
            self.u * np.outer(self.dy, self.dz)
        ).ravel() * c_in
        c = spla.spsolve(a.tocsc(), rhs).reshape(self.nx, self.ny, self.nz)
        self._last_beta = beta
        return np.asarray(c, dtype=np.float64)

    def surface_flux(self, c: FloatArray) -> FloatArray:
        """Flux into the floor per floor cell (nx, ny) [1/s] of the last steady solution."""
        return np.asarray(self._last_beta * c[:, :, 0] * self.floor_area, dtype=np.float64)
