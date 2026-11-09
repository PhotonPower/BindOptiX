"""Z3: captures per spot position along the flow (order bias, D15).

Explorative script (not part of the tests); results in docs/ergebnisse-z3.md.
"""

import math

import numpy as np

from bindoptix.cells import CellModel
from bindoptix.chip import ChipLayout
from bindoptix.chipsim import run_chip


def main() -> None:
    layout = ChipLayout(grid_angle_rad=math.atan(1.0 / 5.0))
    n = 1500
    for mode, k_ns in (("stick", math.inf), ("slide", 0.0)):
        res = run_chip(
            layout, CellModel(), 5e6 / 60.0, n, 6 * 3600.0, 10.0, 7, nonspecific_rate_s=k_ns
        )
        per_col = res.captured_per_spot.reshape(5, 5).sum(axis=0) / res.total_cells
        total = res.captured_per_spot.sum() / res.total_cells
        cols = "  ".join(f"{100 * v:4.1f}%" for v in per_col)
        print(f"{mode:5s}  total {100 * total:4.1f}%  per column upstream→downstream: {cols}")
        print("       per spot:", np.array2string(res.captured_per_spot.reshape(5, 5)))


if __name__ == "__main__":
    main()
