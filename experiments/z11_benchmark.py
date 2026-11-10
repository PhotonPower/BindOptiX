"""Z11: run time of the reference, the CPU core and (if present) the CUDA core.

python experiments/z11_benchmark.py [n_cells] [n_steps]
"""

import sys
import time

import numpy as np

from bindoptix import fast
from bindoptix.cells import CellModel
from bindoptix.flow import DuctFlow
from bindoptix.particles import ChannelModel, Spot, inlet_positions


def main() -> None:
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
    steps = int(sys.argv[2]) if len(sys.argv) > 2 else 720
    spots = tuple(
        Spot(
            np.array([[x - 100, -100], [x + 100, -100], [x + 100, 100], [x - 100, 100]], float),
            1e-2,
        )
        for x in (3400.0, 4200.0, 5000.0, 5800.0, 6600.0)
    )
    m = ChannelModel(DuctFlow(200.0, 50.0, 5e6 / 60), CellModel(), spots=spots, x_out_um=10000.0)
    start = inlet_positions(m.flow, n, 1)
    backends = ["python", "cpu"] + (["cuda"] if fast.cuda_available() else [])
    for b in backends:
        t = time.perf_counter()
        r = fast.run(m, start, steps * 10.0, 10.0, 7, backend=b)
        dt = time.perf_counter() - t
        print(
            f"{b:6s}  {n} cells × {steps} steps: {dt:8.2f} s  "
            f"({n * steps / dt / 1e6:7.2f} M cell-steps/s), captured {int(np.sum(r.state == 2))}"
        )


if __name__ == "__main__":
    main()
