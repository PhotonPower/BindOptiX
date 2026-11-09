"""Z7: end-to-end run of one channel: protocol → capture → frames → readout → counts.

Explorative script (not part of the tests); results in docs/ergebnisse-z7.md.
Placeholder parameters (P) from Z3 and Z6.
"""

import math

import numpy as np

from bindoptix.cells import CellModel
from bindoptix.chip import ChipLayout
from bindoptix.chipsim import channel_model, to_sensor_frame
from bindoptix.frames import Mount
from bindoptix.imaging import FrameGeometry, fixed_pattern, render
from bindoptix.particles import CAPTURED
from bindoptix.protocol import Phase, phase_starts, run_particles
from bindoptix.readout import Tracker, counts_per_spot, detect
from bindoptix.signal import Camera


def main() -> None:
    layout = ChipLayout(grid_angle_rad=math.atan(1 / 5))
    k = 2
    ch = layout.all_channels()[k]
    nl = 1e6 / 60.0  # 1 nL/min in µm³/s
    phases = [
        Phase("load", 40 * 60.0, 5 * nl, cells_per_s=100 / (40 * 60.0)),
        Phase("incubate", 60 * 60.0, 0.0),
        Phase("wash", 10 * 60.0, 20 * nl),
    ]

    def make(q: float):  # type: ignore[no-untyped-def]
        m, _ = channel_model(layout, k, CellModel(), q, Mount.sensor_floor(), math.inf)
        return m

    res = run_particles(make, phases, 10.0, 11, id_offset=k * 1_000_000_000)
    final = res[-1]
    _, cs = channel_model(layout, k, CellModel(), 5 * nl, Mount.sensor_floor())
    spots = [to_sensor_frame(ch, s.polygon_um) for s in cs.spots]
    got = final.state == CAPTURED
    cells_xy = to_sensor_frame(ch, final.position_um[got, :2])
    t_cap = final.event_time_s[got]

    rate_on, sites, off, n_on = 12.6, 1e4, 0.01, 200  # Z6 placeholders (P)
    cam = Camera(2.44, 0.8, 1.5, 0.2, 1.0, 1.2)
    allxy = np.vstack(spots)
    lo, hi = allxy.min(axis=0) - 60.0, allxy.max(axis=0) + 60.0
    geom = FrameGeometry(
        lo[0], lo[1], int((hi[0] - lo[0]) / 2.44), int((hi[1] - lo[1]) / 2.44), 2.44
    )
    pattern = fixed_pattern(geom, 0.3, 10.0, 5)
    areas = tuple((p, sites * off * rate_on * cam.exposure_s) for p in spots)
    base = render(geom, 1.2, areas=areas, area_pattern=pattern, uniform_photons_per_um2=0.2)
    bg = cam.quantum_efficiency * base + cam.dark_e_s * cam.exposure_s
    g = np.random.default_rng(3)
    ref = cam.sample(base, g)
    t_end = phase_starts(phases)[-1]
    # every 10 min, plus a second frame 1 min after the end to confirm late captures
    times = np.append(np.arange(10, int(t_end / 60) + 1, 10) * 60.0, t_end + 60.0)
    tracker = Tracker(4.0, 2)
    for f, t in enumerate(times, start=1):
        on = t_cap <= t
        img = base + render(geom, 1.2, cells_xy[on], np.full(int(on.sum()), n_on * rate_on))
        tracker.update(f, detect(cam.sample(img, g), ref, bg, geom, cam, 0.05))
    counts = counts_per_spot(tracker.confirmed, spots, len(times) + 1, margin_um=5.0)
    print("ROI", geom.nx, "x", geom.ny, "px; captured cells:", int(got.sum()))
    print("t[min]  true  read")
    for f, t in enumerate(times, start=1):
        true = int(np.sum(t_cap <= t))
        print(f"{t / 60:6.0f}  {true:4d}  {int(counts[f].sum()):4d}")
    print(
        "per spot true:",
        np.bincount(final.spot[got], minlength=5).tolist(),
        "read:",
        counts[-1].tolist(),
    )


if __name__ == "__main__":
    main()
