"""Z10: optimisation of the aptamer arrangement and of flow rate and channel height (P).

Explorative script (not part of the tests); results in docs/ergebnisse-z10.md.
"""

import math

import numpy as np

from bindoptix.cells import CellModel
from bindoptix.chip import ChipLayout
from bindoptix.chipsim import channel_model
from bindoptix.classify import BindingMatrix, TrajectoryLibrary, confusion_matrix
from bindoptix.frames import Mount
from bindoptix.optimize import evaluate, optimize_arrangement
from bindoptix.protocol import Phase, run_particles

NL = 1e6 / 60.0


def library(layout: ChipLayout, q: float, n: int, rest_s: float = 4200.0) -> TrajectoryLibrary:
    load_s = 0.2e9 / q
    phases = [Phase("load", load_s, q, cells_per_s=n / load_s), Phase("flow", rest_s, q)]
    runs, maps = [], []
    for k in range(len(layout.all_channels())):

        def make(qq: float, k: int = k):  # type: ignore[no-untyped-def]
            return channel_model(layout, k, CellModel(), qq, Mount.sensor_floor(), 0.0)[0]

        runs.append(run_particles(make, phases, 10.0, 100 + k, k * 10**9, True)[-1])
        maps.append(channel_model(layout, k, CellModel(), q, Mount.sensor_floor())[1])
    return TrajectoryLibrary.from_runs(runs, maps, 25)


def accuracy(lib: TrajectoryLibrary, rates: np.ndarray, arr: tuple[int, ...], n: int) -> float:  # type: ignore[type-arg]
    bm = BindingMatrix(rates, tuple(str(i) for i in range(rates.shape[0])), arr)
    return float(np.mean(np.diag(confusion_matrix(bm.profiles(lib), n, 1000, 2))))


def arrangement_study() -> None:
    layout = ChipLayout(grid_angle_rad=math.atan(1 / 5))
    lib = library(layout, 5 * NL, 600)
    g = np.random.default_rng(1)
    rates = 0.1 * 1e-2 * g.lognormal(0.0, 1.0, (5, 5))
    np.fill_diagonal(rates, 1e-2)
    latin = tuple((j + r) % 5 for r in range(5) for j in range(5))
    blocks = tuple(r for r in range(5) for _ in range(5))  # one aptamer per row
    res = optimize_arrangement(lib, rates, latin, seed=3, restarts=3, max_passes=6)
    print("arrangement (5 classes, cross-reactivity 10 %, 5 nL/min)")
    for name, arr in (
        ("one aptamer per row", blocks),
        ("Latin square", latin),
        ("optimised", res.arrangement),
    ):
        info = evaluate(lib, rates, arr)
        a50, a200 = accuracy(lib, rates, arr, 50), accuracy(lib, rates, arr, 200)
        print(f"  {name:20s} info {1000 * info:6.1f} mbit/cell  acc 50: {a50:.2f}  200: {a200:.2f}")
    print("  optimised rows:", [res.arrangement[5 * r : 5 * r + 5] for r in range(5)])


def flow_height_study() -> None:
    print("capture fraction (any spot, k = 0.01/s), 0.2 µL per channel + 70 min flow")
    print("  h[µm]  Q[nL/min]  captured/cell")
    for h in (50.0, 100.0, 150.0):
        layout = ChipLayout(grid_angle_rad=math.atan(1 / 5), channel_height_um=h)
        for q in (2.0, 5.0, 10.0):
            lib = library(layout, q * NL, 400)
            cap = float(lib.profile(np.full(25, 1e-2)).sum())
            print(f"  {h:5.0f}  {q:9.1f}  {cap:13.3f}")


if __name__ == "__main__":
    arrangement_study()
    flow_height_study()
