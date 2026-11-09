"""Assay protocols: phases of flow, injection and washing (stage Z7).

A protocol is a list of phases (load, incubate, wash, …). Each phase has a duration, a flow rate
(0 for incubation without flow), inlet concentrations of the dissolved species, and a rate of
cells entering the channel. Particles keep their state between phases (``particles.run`` with
``resume``; the flow, hence the model, changes per phase); the continuum keeps its fields (one
model per flow rate on the same grid; u scales with Q in Stokes flow).
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from bindoptix.continuum import ChannelSection, ContinuumState
from bindoptix.continuum3d import ChannelVolume
from bindoptix.particles import ChannelModel, RunResult, inlet_positions, run

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class Phase:
    """One protocol phase.

    Attributes:
        name: label ("load", "incubate", "wash", …).
        duration_s: length.
        flow_rate_um3_s: volume flow (0 allowed).
        c_in: inlet concentrations of the dissolved species [1/µm³].
        cells_per_s: cells entering the channel per second (evenly spaced).
    """

    name: str
    duration_s: float
    flow_rate_um3_s: float
    c_in: tuple[float, ...] = field(default_factory=tuple)
    cells_per_s: float = 0.0

    def __post_init__(self) -> None:
        if not (self.duration_s > 0 and self.flow_rate_um3_s >= 0 and self.cells_per_s >= 0):
            raise ValueError("duration > 0, flow rate >= 0, cells_per_s >= 0 required")


def phase_starts(phases: Sequence[Phase]) -> FloatArray:
    return np.concatenate([[0.0], np.cumsum([p.duration_s for p in phases])])


def injection_schedule(phases: Sequence[Phase]) -> FloatArray:
    """Entry times of the cells: evenly spaced within each phase at its rate."""
    starts = phase_starts(phases)
    out = []
    for p, t0 in zip(phases, starts[:-1], strict=True):
        n = int(math.floor(p.cells_per_s * p.duration_s + 1e-9))
        if n:
            out.append(t0 + (np.arange(n) + 0.5) / p.cells_per_s)
    return np.concatenate(out) if out else np.zeros(0)


def run_particles(
    make_model: Callable[[float], ChannelModel],
    phases: Sequence[Phase],
    dt_s: float,
    seed: int,
    id_offset: int = 0,
    record_visits: bool = False,
) -> list[RunResult]:
    """Cells enter at the inlet according to the schedule; one result per phase end.

    Start positions are flux-weighted for the flow of the phase in which a cell enters (uniform
    over the cross-section if that flow is 0).
    """
    t_in = injection_schedule(phases)
    n = t_in.size
    starts = phase_starts(phases)
    pos = np.zeros((n, 3))
    for k, p in enumerate(phases):
        sel = (t_in >= starts[k]) & (t_in < starts[k + 1])
        m = int(sel.sum())
        if not m:
            continue
        flow = make_model(max(p.flow_rate_um3_s, 1e-9)).flow
        pos[sel] = inlet_positions(flow, m, seed + 7919 * (k + 1))
    ids = np.arange(n, dtype=np.int64) + id_offset
    results: list[RunResult] = []
    prev: RunResult | None = None
    for p in phases:
        model = make_model(p.flow_rate_um3_s)
        if prev is None:
            prev = run(
                model,
                pos,
                p.duration_s,
                dt_s,
                seed,
                ids=ids,
                start_times_s=t_in,
                record_visits=record_visits,
            )
        else:
            prev = run(
                model, None, p.duration_s, dt_s, seed, resume=prev, record_visits=record_visits
            )
        results.append(prev)
    return results


def run_continuum(
    make_model: Callable[[float], ChannelSection | ChannelVolume],
    phases: Sequence[Phase],
    dt_s: float,
    record_every: int = 1,
) -> tuple[list[ContinuumState], FloatArray, list[FloatArray]]:
    """Continuum over the phases; returns states at the phase ends, record times and the
    region-mean coverages (S, regions) at those times."""
    models: dict[float, ChannelSection | ChannelVolume] = {}

    def get(q: float) -> ChannelSection | ChannelVolume:
        if q not in models:
            models[q] = make_model(q)
        return models[q]

    first = get(phases[0].flow_rate_um3_s)
    st = first.initial_state()
    ends: list[ContinuumState] = []
    times: list[float] = [0.0]
    cov: list[FloatArray] = [_region_means(first, st)]
    step = 0
    for p in phases:
        m = get(p.flow_rate_um3_s)
        n = int(math.ceil(p.duration_s / dt_s - 1e-9))
        for _ in range(n):
            st = m.step(st, p.duration_s / n, p.c_in)
            step += 1
            if step % record_every == 0:
                times.append(st.t_s)
                cov.append(_region_means(m, st))
        ends.append(st)
    return ends, np.asarray(times), cov


def _region_means(m: ChannelSection | ChannelVolume, st: ContinuumState) -> FloatArray:
    ns = st.gamma.shape[0]
    out = np.zeros((ns, len(m.regions)))
    if isinstance(m, ChannelVolume):
        area = m.floor_area
    else:
        area = m.dx
    for j in range(len(m.regions)):
        sel = m.region_of == j
        a = area[sel]
        for q in range(ns):
            out[q, j] = float(np.sum(st.gamma[q][sel] * a) / a.sum())
    return out
