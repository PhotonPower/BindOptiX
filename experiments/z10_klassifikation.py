"""Z10: classification of 25 strains vs flow rate and cross-reactivity (synthetic, P).

Explorative script (not part of the tests); results in docs/ergebnisse-z10.md.
"""

import math

import numpy as np

from bindoptix.cells import CellModel
from bindoptix.chip import ChipLayout
from bindoptix.chipsim import channel_model
from bindoptix.classify import BindingMatrix, TrajectoryLibrary, confusion_matrix
from bindoptix.frames import Mount
from bindoptix.protocol import Phase, run_particles


def library(layout: ChipLayout, q_load: float, n_per_channel: int) -> TrajectoryLibrary:
    nl = 1e6 / 60.0
    load_s = 0.2e9 / q_load  # 0.2 µL concentrate per channel
    phases = [
        Phase("load", load_s, q_load, cells_per_s=n_per_channel / load_s),
        Phase("incubate", 3600.0, 0.0),
        Phase("wash", 600.0, 20 * nl),
    ]
    runs, maps = [], []
    for k in range(len(layout.all_channels())):

        def make(q: float, k: int = k):  # type: ignore[no-untyped-def]
            return channel_model(layout, k, CellModel(), q, Mount.sensor_floor(), 0.0)[0]

        res = run_particles(make, phases, 10.0, 100 + k, id_offset=k * 10**9, record_visits=True)
        runs.append(res[-1])
        maps.append(channel_model(layout, k, CellModel(), q_load, Mount.sensor_floor())[1])
    return TrajectoryLibrary.from_runs(runs, maps, 25)


def binding(cross: float, seed: int) -> BindingMatrix:
    g = np.random.default_rng(seed)
    k_spec = 1e-2
    r = cross * k_spec * g.lognormal(0.0, 1.0, (25, 25))
    np.fill_diagonal(r, k_spec)
    return BindingMatrix(r, tuple(f"S{i}" for i in range(25)), tuple(range(25)))


def main() -> None:
    layout = ChipLayout(grid_angle_rad=math.atan(1 / 5))
    nl = 1e6 / 60.0
    print("Q[nL/min] cross  captured/cell  accuracy(50 cells)  accuracy(200 cells)")
    for q in (1.5, 5.0, 15.0):
        lib = library(layout, q * nl, 600)
        for cross in (0.01, 0.1, 0.3):
            bm = binding(cross, 1)
            prof = bm.profiles(lib)
            cap = float(np.mean(prof.sum(axis=1)))
            acc = []
            for n in (50, 200):
                cm = confusion_matrix(prof, n, 400, 2)
                acc.append(float(np.mean(np.diag(cm))))
            print(f"{q:8.1f}  {cross:5.2f}  {cap:13.3f}  {acc[0]:18.3f}  {acc[1]:19.3f}")


def latin_rows() -> tuple[int, ...]:
    """Five aptamers, every row carries all five in a cyclically shifted order."""
    return tuple((j + r) % 5 for r in range(5) for j in range(5))


def main_five_classes() -> None:
    layout = ChipLayout(grid_angle_rad=math.atan(1 / 5))
    nl = 1e6 / 60.0
    lib = library(layout, 5 * nl, 600)
    g = np.random.default_rng(1)
    print("5 classes, all aptamers in every row (Latin square), 5 nL/min")
    print("cross  captured/cell  accuracy(10)  accuracy(50)  accuracy(200)")
    for cross in (0.01, 0.1, 0.3):
        r = cross * 1e-2 * g.lognormal(0.0, 1.0, (5, 5))
        np.fill_diagonal(r, 1e-2)
        bm = BindingMatrix(r, tuple(f"C{i}" for i in range(5)), latin_rows())
        prof = bm.profiles(lib)
        cap = float(np.mean(prof.sum(axis=1)))
        acc = [float(np.mean(np.diag(confusion_matrix(prof, n, 1000, 2)))) for n in (10, 50, 200)]
        print(f"{cross:5.2f}  {cap:13.3f}  {acc[0]:12.3f}  {acc[1]:12.3f}  {acc[2]:13.3f}")


if __name__ == "__main__":
    main()
    main_five_classes()
