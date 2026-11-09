"""Lagrangian particles in a channel: Brownian motion, sedimentation, swimming, flow, capture.

Stage Z3 (docs/roadmap.md). Overdamped Langevin dynamics in the channel frame (x along the flow,
y across, z above the floor, µm, s), Euler–Maruyama step

    dx = (u(y, z) x̂ + v_s ĝ + v ê) dt + √(2D dt) ξ,

with u from :class:`bindoptix.flow.DuctFlow`, v_s ĝ the sedimentation velocity along gravity and
v ê the run of a swimming cell; tumbles with rate 1/τ turn ê by the fixed angle arccos α about a
random axis. Side walls and the ceiling reflect. The floor either reflects (equilibrium checks) or
is a contact surface: a crossing is detected also between steps with the Brownian-bridge
probability exp(−z₀z₁/(D dt)) (berg1993-type first passage, exact for pure diffusion).

On the floor a cell slides with the flow, u_wall = s · a · ∂u/∂z(z = 0) (s: slip factor, P,
goldman1967; a: equivalent radius), and diffuses laterally with D·f_floor (P). Over a spot k it is
captured with rate k_cap[k]; elsewhere it adheres non-specifically with rate k_ns (∞: cells stick
where they land, 0: they slide until they reach a spot). A swimming cell on the floor whose run
points upwards leaves it. Cells beyond ``x_out_um`` leave the channel.

Random numbers are counter-based (``bindoptix.rng``): results depend only on the seed and the
particle ids, not on batching (T11).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import ArrayLike, NDArray

from bindoptix import rng
from bindoptix.cells import CellModel
from bindoptix.flow import DuctFlow

FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]

FREE, FLOOR, CAPTURED, STUCK, EXITED = 0, 1, 2, 3, 4
STATE_NAMES = ("free", "floor", "captured", "stuck", "exited")


@dataclass(frozen=True)
class Spot:
    """Capture area on the floor: convex polygon (n, 2) in the channel frame, rate [1/s]."""

    polygon_um: FloatArray
    capture_rate_s: float

    def __post_init__(self) -> None:
        p = np.asarray(self.polygon_um, dtype=np.float64)
        if p.ndim != 2 or p.shape[1] != 2 or p.shape[0] < 3:
            raise ValueError("polygon_um must have shape (n >= 3, 2)")
        x, y = p[:, 0], p[:, 1]
        if 0.5 * (x @ np.roll(y, -1) - y @ np.roll(x, -1)) <= 0.0:
            raise ValueError("polygon_um must be counter-clockwise")
        if not self.capture_rate_s >= 0.0:
            raise ValueError("capture_rate_s must be non-negative (inf allowed)")
        object.__setattr__(self, "polygon_um", p)

    def contains(self, x: FloatArray, y: FloatArray) -> NDArray[np.bool_]:
        p = self.polygon_um
        inside = np.ones(x.shape, dtype=bool)
        for i in range(p.shape[0]):
            a, b = p[i], p[(i + 1) % p.shape[0]]
            inside &= (b[0] - a[0]) * (y - a[1]) - (b[1] - a[1]) * (x - a[0]) >= 0.0
        return inside


@dataclass(frozen=True)
class ChannelModel:
    """Everything a particle run needs.

    Attributes:
        flow: channel cross-section and flow rate.
        cell: transport properties of the cells.
        gravity_dir: unit vector of gravity in the channel frame (sensor as floor: (0, 0, −1)).
        floor: "contact" or "reflect".
        nonspecific_rate_s: adhesion rate off the spots on the floor (P; ∞ = stick on landing).
        slip_factor: s in u_wall = s · a · ∂u/∂z (P: 0.5).
        floor_diffusion_factor: lateral diffusion on the floor relative to D (P: 0.5).
        spots: capture areas.
        x_out_um: cells beyond this x leave the channel.
        brownian: switch off Brownian motion for deterministic checks.
    """

    flow: DuctFlow
    cell: CellModel
    gravity_dir: tuple[float, float, float] = (0.0, 0.0, -1.0)
    floor: str = "contact"
    nonspecific_rate_s: float = 0.0
    slip_factor: float = 0.5
    floor_diffusion_factor: float = 0.5
    spots: tuple[Spot, ...] = field(default_factory=tuple)
    x_out_um: float = math.inf
    brownian: bool = True

    def __post_init__(self) -> None:
        if self.floor not in ("contact", "reflect"):
            raise ValueError("floor must be 'contact' or 'reflect'")
        g = np.asarray(self.gravity_dir, dtype=np.float64)
        n = float(np.linalg.norm(g))
        if g.shape != (3,) or n == 0.0:
            raise ValueError("gravity_dir must be a non-zero 3-vector")
        object.__setattr__(self, "gravity_dir", tuple(float(v) for v in g / n))
        for name in ("nonspecific_rate_s", "slip_factor", "floor_diffusion_factor"):
            if not getattr(self, name) >= 0.0:
                raise ValueError(f"{name} must be non-negative")


@dataclass
class RunResult:
    """Final state of a run; ``spot`` is the capturing spot index or −1."""

    position_um: FloatArray
    state: IntArray
    spot: IntArray
    event_time_s: FloatArray
    snapshots_um: dict[float, FloatArray]

    def count(self, state: int) -> int:
        return int(np.count_nonzero(self.state == state))


def _random_unit(seed: int, pid: IntArray, step: int, lane: int) -> FloatArray:
    z = 2.0 * rng.uniform(seed, pid, step, lane) - 1.0
    phi = 2.0 * math.pi * rng.uniform(seed, pid, step, lane + 1)
    r = np.sqrt(1.0 - z * z)
    return np.column_stack([r * np.cos(phi), r * np.sin(phi), z])


def _turn(e: FloatArray, cos_a: float, phi: FloatArray) -> FloatArray:
    """Rotate unit vectors e by the polar angle arccos(cos_a) about a random azimuth phi."""
    helper = np.where(np.abs(e[:, 2:3]) < 0.9, [[0.0, 0.0, 1.0]], [[1.0, 0.0, 0.0]])
    p1 = np.cross(e, helper)
    p1 /= np.linalg.norm(p1, axis=1, keepdims=True)
    p2 = np.cross(e, p1)
    sin_a = math.sqrt(max(0.0, 1.0 - cos_a * cos_a))
    out = cos_a * e + sin_a * (np.cos(phi)[:, None] * p1 + np.sin(phi)[:, None] * p2)
    return np.asarray(out / np.linalg.norm(out, axis=1, keepdims=True))


def run(
    model: ChannelModel,
    start_um: ArrayLike,
    duration_s: float,
    dt_s: float,
    seed: int,
    ids: ArrayLike | None = None,
    snapshot_times_s: tuple[float, ...] = (),
) -> RunResult:
    """Simulate cells from ``start_um`` (N, 3) for ``duration_s`` with time step ``dt_s``.

    ``ids`` (default 0 … N−1) select the random streams; the same ids give the same paths in any
    batch. Snapshots record positions at the first step at or after each requested time.
    """
    x = np.array(start_um, dtype=np.float64).reshape(-1, 3)
    n = x.shape[0]
    pid = np.arange(n, dtype=np.int64) if ids is None else np.asarray(ids, dtype=np.int64)
    if pid.shape != (n,):
        raise ValueError("ids must have one entry per particle")
    if not (dt_s > 0.0 and duration_s >= 0.0):
        raise ValueError("dt_s must be positive and duration_s non-negative")
    fl, cell = model.flow, model.cell
    w2, h = 0.5 * fl.width_um, fl.height_um
    d = cell.diffusion_um2_s if model.brownian else 0.0
    sig = math.sqrt(2.0 * d * dt_s)
    sig_floor = math.sqrt(2.0 * d * model.floor_diffusion_factor * dt_s)
    v_sed = cell.sedimentation_speed_um_s() * np.asarray(model.gravity_dir)
    v_swim = cell.swim_speed_um_s
    p_tumble = -math.expm1(-dt_s / cell.run_time_s)
    a_eq = cell.equivalent_radius_um
    rates = np.array([s.capture_rate_s for s in model.spots], dtype=np.float64)

    state = np.full(n, FREE, dtype=np.int64)
    if model.floor == "contact":
        state[x[:, 2] <= 0.0] = FLOOR
        x[x[:, 2] <= 0.0, 2] = 0.0
    spot = np.full(n, -1, dtype=np.int64)
    t_event = np.full(n, math.nan)
    e = _random_unit(seed, pid, 0, 100) if v_swim > 0.0 else np.zeros((n, 3))
    snaps: dict[float, FloatArray] = {}
    pending = sorted(snapshot_times_s)
    n_steps = int(math.ceil(duration_s / dt_s - 1e-9))

    for step in range(1, n_steps + 1):
        t = step * dt_s
        free = np.flatnonzero(state == FREE)
        if free.size:
            p, ids_f = x[free], pid[free]
            u = fl.velocity_um_s(p[:, 1], p[:, 2]) if fl.flow_rate_um3_s > 0.0 else 0.0
            drift = np.zeros_like(p)
            drift[:, 0] += u
            drift += v_sed + v_swim * e[free]
            noise = np.column_stack([rng.normal(seed, ids_f, step, k) for k in range(3)])
            z0 = p[:, 2].copy()
            q = p + drift * dt_s + sig * noise
            q[:, 1] = np.where(q[:, 1] > w2, 2.0 * w2 - q[:, 1], q[:, 1])
            q[:, 1] = np.where(q[:, 1] < -w2, -2.0 * w2 - q[:, 1], q[:, 1])
            q[:, 2] = np.where(q[:, 2] > h, 2.0 * h - q[:, 2], q[:, 2])
            if model.floor == "reflect":
                q[:, 2] = np.abs(q[:, 2])
            else:
                hit = q[:, 2] <= 0.0
                if d > 0.0:
                    both = ~hit
                    pb = np.exp(-np.maximum(z0, 0.0) * np.maximum(q[:, 2], 0.0) / (d * dt_s))
                    hit |= both & (rng.uniform(seed, ids_f, step, 10) < pb)
                q[hit, 2] = 0.0
                state[free[hit]] = FLOOR
            x[free] = q

        on = np.flatnonzero(state == FLOOR)
        if on.size:
            p, ids_o = x[on], pid[on]
            slide = model.slip_factor * a_eq * fl.floor_shear_rate_s(p[:, 1])
            p[:, 0] += slide * dt_s + sig_floor * rng.normal(seed, ids_o, step, 20)
            p[:, 1] += sig_floor * rng.normal(seed, ids_o, step, 21)
            p[:, 1] = np.clip(p[:, 1], -w2, w2)
            x[on] = p
            r = np.full(on.size, model.nonspecific_rate_s)
            which = np.full(on.size, -1, dtype=np.int64)
            for k, s in enumerate(model.spots):
                inside = (which < 0) & s.contains(p[:, 0], p[:, 1])
                which[inside] = k
                r[inside] = rates[k]
            prob = np.where(np.isinf(r), 1.0, -np.expm1(-r * dt_s))
            caught = rng.uniform(seed, ids_o, step, 22) < prob
            spec = caught & (which >= 0)
            state[on[spec]] = CAPTURED
            spot[on[spec]] = which[spec]
            state[on[caught & (which < 0)]] = STUCK
            t_event[on[caught]] = t
            if v_swim > 0.0:
                lift = (~caught) & (e[on, 2] > 0.0)
                state[on[lift]] = FREE

        if v_swim > 0.0:
            mov = np.flatnonzero((state == FREE) | (state == FLOOR))
            if mov.size:
                tumble = rng.uniform(seed, pid[mov], step, 30) < p_tumble
                if tumble.any():
                    phi = 2.0 * math.pi * rng.uniform(seed, pid[mov[tumble]], step, 31)
                    e[mov[tumble]] = _turn(e[mov[tumble]], cell.tumble_cos, phi)

        out = ((state == FREE) | (state == FLOOR)) & (x[:, 0] > model.x_out_um)
        state[out] = EXITED
        t_event[out] = t
        while pending and t >= pending[0] - 1e-9 * dt_s:
            snaps[pending.pop(0)] = x.copy()

    return RunResult(x, state, spot, t_event, snaps)


def inlet_positions(flow: DuctFlow, n: int, seed: int, x_um: float = 0.0) -> FloatArray:
    """Flux-weighted start positions over the cross-section (rejection sampling)."""
    g = np.random.default_rng(seed)
    u_max = float(flow.velocity_um_s([0.0], [0.5 * flow.height_um])[0])
    out: list[FloatArray] = []
    have = 0
    while have < n:
        y = g.uniform(-0.5 * flow.width_um, 0.5 * flow.width_um, 4 * n)
        z = g.uniform(0.0, flow.height_um, 4 * n)
        keep = g.uniform(0.0, u_max * 1.0001, 4 * n) < flow.velocity_um_s(y, z)
        pts = np.column_stack([np.full(keep.sum(), x_um), y[keep], z[keep]])
        out.append(pts)
        have += pts.shape[0]
    return np.concatenate(out)[:n]
