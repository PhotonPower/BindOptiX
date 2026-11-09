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
from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import ArrayLike, NDArray

from bindoptix import rng, wall
from bindoptix.cells import CellModel
from bindoptix.flow import DuctFlow

FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]

FREE, FLOOR, CAPTURED, STUCK, EXITED, WAITING = 0, 1, 2, 3, 4, 5
STATE_NAMES = ("free", "floor", "captured", "stuck", "exited", "waiting")


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
        wall_hydrodynamics: hindered diffusion and sedimentation near the floor
            (``bindoptix.wall``); z is then the gap between cell surface and floor.
        contact_gap_um: gap at which a cell counts as touching the floor (P; required > 0 with
            wall hydrodynamics, because λ⊥ → 0 at the wall).
        spot_rates_fn: optional capture rates per spot as a function of time [1/s]; overrides
            the constant rates of the spots (coupling to fouling, Z5).
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
    wall_hydrodynamics: bool = False
    contact_gap_um: float = 0.0
    spot_rates_fn: Callable[[float], FloatArray] | None = None

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
        if not (math.isfinite(self.contact_gap_um) and self.contact_gap_um >= 0.0):
            raise ValueError("contact_gap_um must be finite and non-negative")
        if self.wall_hydrodynamics and self.floor == "contact" and self.contact_gap_um <= 0.0:
            raise ValueError("wall hydrodynamics need contact_gap_um > 0 (λ⊥ vanishes at the wall)")


@dataclass
class RunResult:
    """Final state of a run; ``spot`` is the capturing spot index or −1."""

    position_um: FloatArray
    state: IntArray
    spot: IntArray
    event_time_s: FloatArray
    snapshots_um: dict[float, FloatArray]
    orientation: FloatArray = field(default_factory=lambda: np.zeros((0, 3)))
    ids: IntArray = field(default_factory=lambda: np.zeros(0, dtype=np.int64))
    steps: int = 0
    t_s: float = 0.0
    start_times_s: FloatArray = field(default_factory=lambda: np.zeros(0))
    visit_entry_s: FloatArray = field(default_factory=lambda: np.zeros((0, 0)))
    visit_time_s: FloatArray = field(default_factory=lambda: np.zeros((0, 0)))
    visit_off_floor_s: FloatArray = field(default_factory=lambda: np.zeros((0, 0)))
    off_floor_s: FloatArray = field(default_factory=lambda: np.zeros(0))

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
    start_um: ArrayLike | None,
    duration_s: float,
    dt_s: float,
    seed: int,
    ids: ArrayLike | None = None,
    snapshot_times_s: tuple[float, ...] = (),
    resume: RunResult | None = None,
    start_times_s: ArrayLike | None = None,
    record_visits: bool = False,
) -> RunResult:
    """Simulate cells from ``start_um`` (N, 3) for ``duration_s`` with time step ``dt_s``.

    ``ids`` (default 0 … N−1) select the random streams; the same ids give the same paths in any
    batch. Snapshots record positions at the first step at or after each requested (absolute)
    time. ``resume`` continues a previous run (its positions, states, orientations, ids, clock
    and step counter; ``start_um`` and ``ids`` are then ignored), e.g. with another model for
    the next phase of a protocol. ``start_times_s`` lets cells enter later: they wait at their
    start position until then (continuous injection). ``record_visits`` records for every cell
    and spot the first entry time, the time spent on the spot and the floor time off the spots
    before that entry (trajectory library for Z10; meaningful with capture rates 0).
    """
    if not (dt_s > 0.0 and duration_s >= 0.0):
        raise ValueError("dt_s must be positive and duration_s non-negative")
    if resume is not None:
        x = resume.position_um.copy()
        n = x.shape[0]
        pid = resume.ids.copy()
        step0, t0 = resume.steps, resume.t_s
        t_start = resume.start_times_s.copy()
    else:
        if start_um is None:
            raise ValueError("start_um is required without resume")
        x = np.array(start_um, dtype=np.float64).reshape(-1, 3)
        n = x.shape[0]
        pid = np.arange(n, dtype=np.int64) if ids is None else np.asarray(ids, dtype=np.int64)
        step0, t0 = 0, 0.0
        t_start = (
            np.zeros(n) if start_times_s is None else np.asarray(start_times_s, dtype=np.float64)
        )
    if pid.shape != (n,):
        raise ValueError("ids must have one entry per particle")
    if t_start.shape != (n,):
        raise ValueError("start_times_s must have one entry per particle")
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

    if resume is not None:
        state, spot = resume.state.copy(), resume.spot.copy()
        t_event = resume.event_time_s.copy()
        e = resume.orientation.copy() if resume.orientation.shape == (n, 3) else np.zeros((n, 3))
        if v_swim > 0.0 and not np.any(e):
            e = _random_unit(seed, pid, 0, 100)
    else:
        state = np.where(t_start > 0.0, WAITING, FREE).astype(np.int64)
        if model.floor == "contact":
            on_floor = (state == FREE) & (x[:, 2] <= 0.0)
            state[on_floor] = FLOOR
            x[on_floor, 2] = 0.0
        spot = np.full(n, -1, dtype=np.int64)
        t_event = np.full(n, math.nan)
        e = _random_unit(seed, pid, 0, 100) if v_swim > 0.0 else np.zeros((n, 3))
    snaps: dict[float, FloatArray] = {}
    pending = sorted(snapshot_times_s)
    n_steps = int(math.ceil(duration_s / dt_s - 1e-9))
    ns = len(model.spots)
    if resume is not None and resume.visit_entry_s.shape == (n, ns):
        v_entry, v_time = resume.visit_entry_s.copy(), resume.visit_time_s.copy()
        v_off, off_floor = resume.visit_off_floor_s.copy(), resume.off_floor_s.copy()
    else:
        v_entry, v_time = np.full((n, ns), math.nan), np.zeros((n, ns))
        v_off, off_floor = np.full((n, ns), math.nan), np.zeros(n)

    for k_step in range(1, n_steps + 1):
        step = step0 + k_step
        t = t0 + k_step * dt_s
        wake = (state == WAITING) & (t_start <= t - dt_s + 1e-9 * dt_s)
        if wake.any():
            state[wake] = FREE
            if model.floor == "contact":
                down = wake & (x[:, 2] <= 0.0)
                state[down] = FLOOR
                x[down, 2] = 0.0
        free = np.flatnonzero(state == FREE)
        if free.size:
            p, ids_f = x[free], pid[free]
            u = fl.velocity_um_s(p[:, 1], p[:, 2]) if fl.flow_rate_um3_s > 0.0 else 0.0
            drift = np.zeros_like(p)
            drift[:, 0] += u
            z0 = p[:, 2].copy()
            if model.wall_hydrodynamics:
                lam_p = wall.lambda_perp(z0, a_eq)
                lam_s = wall.lambda_par(z0, a_eq)
                drift[:, :2] += v_sed[:2] * lam_s[:, None]
                drift[:, 2] += v_sed[2] * lam_p + d * wall.dlambda_perp(z0, a_eq)
                drift += v_swim * e[free]
                scale = np.column_stack([np.sqrt(lam_s), np.sqrt(lam_s), np.sqrt(lam_p)])
                d_perp = d * lam_p
            else:
                drift += v_sed + v_swim * e[free]
                scale = np.ones_like(p)
                d_perp = np.full(free.size, d)
            noise = np.column_stack([rng.normal(seed, ids_f, step, k) for k in range(3)])
            q = p + drift * dt_s + sig * scale * noise
            q[:, 1] = np.where(q[:, 1] > w2, 2.0 * w2 - q[:, 1], q[:, 1])
            q[:, 1] = np.where(q[:, 1] < -w2, -2.0 * w2 - q[:, 1], q[:, 1])
            q[:, 2] = np.where(q[:, 2] > h, 2.0 * h - q[:, 2], q[:, 2])
            if model.floor == "reflect":
                q[:, 2] = np.abs(q[:, 2])
            else:
                zc = model.contact_gap_um
                hit = q[:, 2] <= zc
                if d > 0.0:
                    both = ~hit
                    g0 = np.maximum(z0 - zc, 0.0)
                    g1 = np.maximum(q[:, 2] - zc, 0.0)
                    with np.errstate(divide="ignore", invalid="ignore"):
                        pb = np.where(d_perp > 0.0, np.exp(-g0 * g1 / (d_perp * dt_s)), 0.0)
                    hit |= both & (rng.uniform(seed, ids_f, step, 10) < pb)
                q[hit, 2] = zc
                state[free[hit]] = FLOOR
            x[free] = q

        if model.spot_rates_fn is not None:
            rates = np.asarray(model.spot_rates_fn(t), dtype=np.float64)
            if rates.shape != (len(model.spots),) or np.any(rates < 0.0):
                raise ValueError("spot_rates_fn must return one non-negative rate per spot")
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
            if record_visits:
                off_floor[on[which < 0]] += dt_s
                ins = np.flatnonzero(which >= 0)
                pi, si = on[ins], which[ins]
                first = np.isnan(v_entry[pi, si])
                v_entry[pi[first], si[first]] = t
                v_off[pi[first], si[first]] = off_floor[pi[first]]
                v_time[pi, si] += dt_s
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

    return RunResult(
        x,
        state,
        spot,
        t_event,
        snaps,
        e,
        pid,
        step0 + n_steps,
        t0 + n_steps * dt_s,
        t_start,
        v_entry if record_visits else np.zeros((0, 0)),
        v_time if record_visits else np.zeros((0, 0)),
        v_off if record_visits else np.zeros((0, 0)),
        off_floor if record_visits else np.zeros(0),
    )


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
