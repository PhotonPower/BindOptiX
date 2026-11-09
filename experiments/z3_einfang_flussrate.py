"""Z3: capture on the five spots of one channel vs. flow rate, cells stick or slide.

Explorative script (not part of the tests); results in docs/ergebnisse-z3.md.
"""

import math

import numpy as np

from bindoptix.cells import CellModel
from bindoptix.flow import DuctFlow
from bindoptix.particles import (
    CAPTURED,
    EXITED,
    FLOOR,
    FREE,
    STUCK,
    ChannelModel,
    Spot,
    inlet_positions,
    run,
)


def spot(x: float) -> Spot:
    poly = np.array([[x - 100, -100], [x + 100, -100], [x + 100, 100], [x - 100, 100]], float)
    return Spot(poly, math.inf)


def main() -> None:
    cell = CellModel()
    spots = tuple(spot(x) for x in (3400.0, 4200.0, 5000.0, 5800.0, 6600.0))
    n = 1500
    for q_nl_min in (0.5, 1.5, 5.0, 15.0):
        flow = DuctFlow(200.0, 100.0, q_nl_min * 1e6 / 60.0)
        u_wall = 0.5 * cell.equivalent_radius_um * float(flow.floor_shear_rate_s([0.0])[0])
        for mode, k_ns in (("stick", math.inf), ("slide", 0.0)):
            m = ChannelModel(flow, cell, nonspecific_rate_s=k_ns, spots=spots, x_out_um=10000.0)
            r = run(m, inlet_positions(flow, n, 1), 6 * 3600.0, 10.0, 7)
            left = float(np.mean((r.state == FREE) | (r.state == FLOOR)))
            print(
                f"Q={q_nl_min:4} nL/min  mean={flow.mean_speed_um_s:6.2f} um/s  "
                f"u_wall={u_wall:5.3f} um/s  {mode:5s}  spots={r.count(CAPTURED) / n:.3f}  "
                f"stuck={r.count(STUCK) / n:.3f}  out={r.count(EXITED) / n:.3f}  left={left:.3f}"
            )


if __name__ == "__main__":
    main()
