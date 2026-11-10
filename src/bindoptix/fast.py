"""Facade for the compiled core ``bindoptix._core`` (stage Z11, like BWPM's ``bwpm.fast``).

:func:`run` has the signature and the result of :func:`bindoptix.particles.run` and computes the
same paths (same equations, order of operations and counter-based random numbers; agreement up
to the rounding of exp/log/sin/cos). ``backend``: "auto" (CUDA if built and a device is
present, else the CPU core), "cpu", "cuda", or "python" (the reference). Without the core
(``BINDOPTIX_BUILD_CORE=OFF``) "auto" falls back to the reference. ``CORE_AVAILABLE`` and
:func:`cuda_available` tell what is active.

Differences to the reference: ``spot_rates_fn`` is evaluated once per step in Python and passed
as a table; the core runs in double precision.
"""

from __future__ import annotations

import math
import os

import numpy as np
from numpy.typing import ArrayLike, NDArray

from bindoptix import particles as _ref
from bindoptix.particles import (
    FLOOR,
    FREE,
    WAITING,
    ChannelModel,
    RunResult,
    _random_unit,
)

try:
    from bindoptix import _core  # type: ignore[attr-defined]

    CORE_AVAILABLE = True
except ImportError:  # pure Python installation
    _core = None
    CORE_AVAILABLE = False

if os.environ.get("BINDOPTIX_REQUIRE_CORE") == "1" and not CORE_AVAILABLE:
    raise ImportError("BINDOPTIX_REQUIRE_CORE=1 but bindoptix._core is not built")

FloatArray = NDArray[np.float64]


def cuda_available() -> bool:
    return bool(CORE_AVAILABLE and _core.cuda_built() and _core.cuda_available())


def _resolve(backend: str) -> str:
    if backend == "auto":
        if not CORE_AVAILABLE:
            return "python"
        return "cuda" if cuda_available() else "cpu"
    if backend in ("cpu", "cuda") and not CORE_AVAILABLE:
        raise RuntimeError("bindoptix._core is not built")
    if backend == "cuda" and not cuda_available():
        raise RuntimeError("CUDA backend not built or no device")
    if backend not in ("cpu", "cuda", "python"):
        raise ValueError("backend must be auto, cpu, cuda or python")
    return backend


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
    backend: str = "auto",
) -> RunResult:
    """Like :func:`bindoptix.particles.run`, on the compiled core."""
    which = _resolve(backend)
    if which == "python":
        return _ref.run(
            model,
            start_um,
            duration_s,
            dt_s,
            seed,
            ids,
            snapshot_times_s,
            resume,
            start_times_s,
            record_visits,
        )
    if not (dt_s > 0.0 and duration_s >= 0.0):
        raise ValueError("dt_s must be positive and duration_s non-negative")
    # --- initial state, exactly as the reference -------------------------------------------
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
        t_start = np.zeros(n) if start_times_s is None else np.asarray(start_times_s, float)
    if pid.shape != (n,) or t_start.shape != (n,):
        raise ValueError("ids and start_times_s need one entry per particle")
    fl, cell = model.flow, model.cell
    v_swim = cell.swim_speed_um_s
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
    ns = len(model.spots)
    if resume is not None and resume.visit_entry_s.shape == (n, ns):
        v_entry, v_time = resume.visit_entry_s.copy(), resume.visit_time_s.copy()
        v_off, off_floor = resume.visit_off_floor_s.copy(), resume.off_floor_s.copy()
    else:
        v_entry, v_time = np.full((n, ns), math.nan), np.zeros((n, ns))
        v_off, off_floor = np.full((n, ns), math.nan), np.zeros(n)
    n_steps = int(math.ceil(duration_s / dt_s - 1e-9))
    # --- tables: rates per step, snapshot steps, spot polygons ------------------------------
    times = t0 + dt_s * np.arange(1, n_steps + 1)
    if model.spot_rates_fn is not None:
        rates = np.array([model.spot_rates_fn(float(t)) for t in times], dtype=np.float64)
        if rates.shape != (n_steps, ns) or np.any(rates < 0.0):
            raise ValueError("spot_rates_fn must return one non-negative rate per spot")
    else:
        rates = np.tile([s.capture_rate_s for s in model.spots], (n_steps, 1)).astype(float)
    snap_t = sorted(snapshot_times_s)
    snap_k: list[int] = []
    keys: list[float] = []
    for ts in snap_t:
        ok = np.flatnonzero(times >= ts - 1e-9 * dt_s)
        if ok.size:
            snap_k.append(int(ok[0]) + 1)
            keys.append(ts)
    offsets = np.cumsum([0] + [s.polygon_um.shape[0] for s in model.spots]).astype(np.int32)
    poly = (np.vstack([s.polygon_um for s in model.spots]) if ns else np.zeros((0, 2))).astype(
        np.float64
    )
    snaps = np.zeros((len(snap_k), n, 3))
    params = {
        "width": fl.width_um,
        "height": fl.height_um,
        "pressure_gradient": fl.pressure_gradient,
        "n_terms": fl.n_terms,
        "has_flow": int(fl.flow_rate_um3_s > 0.0),
        "diffusion": cell.diffusion_um2_s,
        "radius": cell.equivalent_radius_um,
        "vsed_x": cell.sedimentation_speed_um_s() * model.gravity_dir[0],
        "vsed_y": cell.sedimentation_speed_um_s() * model.gravity_dir[1],
        "vsed_z": cell.sedimentation_speed_um_s() * model.gravity_dir[2],
        "swim_speed": v_swim,
        "p_tumble": -math.expm1(-dt_s / cell.run_time_s),
        "tumble_cos": cell.tumble_cos,
        "floor_contact": int(model.floor == "contact"),
        "k_ns": model.nonspecific_rate_s,
        "slip": model.slip_factor,
        "floor_diffusion_factor": model.floor_diffusion_factor,
        "x_out": model.x_out_um,
        "contact_gap": model.contact_gap_um,
        "brownian": int(model.brownian),
        "wall_hydrodynamics": int(model.wall_hydrodynamics),
        "n_spots": ns,
        "seed": seed & 0xFFFFFFFFFFFFFFFF,
        "dt": dt_s,
        "t0": t0,
        "n_steps": n_steps,
        "step0": step0,
        "record_visits": int(record_visits),
        "n_snap": len(snap_k),
    }
    _core.run(
        params,
        x,
        state,
        spot,
        t_event,
        e,
        t_start.astype(np.float64),
        pid.astype(np.uint64),
        v_entry,
        v_time,
        v_off,
        off_floor,
        snaps,
        offsets,
        poly,
        rates.reshape(-1),
        np.asarray(snap_k, dtype=np.int32),
        which,
    )
    return RunResult(
        x,
        state,
        spot,
        t_event,
        {k: snaps[j].copy() for j, k in enumerate(keys)},
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
