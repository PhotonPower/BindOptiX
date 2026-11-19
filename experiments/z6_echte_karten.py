"""Z6: photon budget with the field map and the response table exported by hp-FEM.

Explorative script (not part of the tests); results in docs/ergebnisse-z6.md.
Arrangement (D4): excitation from the solution side at normal incidence, collection through the
substrate (SIL below the chip). The maps replace the placeholders of z6_photonenbudget.py; the
dye, the power, the camera and the heights of the dyes stay placeholders (P).

Usage: python experiments/z6_echte_karten.py FIELD_MAP.json RESPONSE_TABLE.json [REGION ...]
REGION: detector regions taken in by the SIL (default: every region whose name starts with
"substrate").
"""

import math
import sys
from pathlib import Path

import numpy as np

from bindoptix.ports import EmitterResponseTable, FieldMap
from bindoptix.signal import (
    Camera,
    EmitterEnsemble,
    Fluorophore,
    detected_rate_s,
    excitation_rate_s,
    spot_detection,
    uniform_field_map,
)

HEIGHTS_NM = (2.0, 3.0, 5.0, 10.0, 20.0)  # dye above the gold (P, until BioOptoMD)


def main(argv: list[str]) -> None:
    if len(argv) < 2:
        raise SystemExit(__doc__)
    dye = Fluorophore(250_000.0, 0.27, 0.640, 0.670, bleach_yield=1e-6)  # Cy5-like (P)
    power_mw, field_mm = 10.0, 5.0
    irr_per_mw = 1e-3 / (math.pi * (0.5 * field_mm * 1e-3) ** 2)  # W/m² per mW (P)
    field = FieldMap.from_json(Path(argv[0]).read_text(encoding="utf-8"))
    if field.irradiance_w_m2 is not None:
        field = field.per_port_power(irr_per_mw, "P")
    table = EmitterResponseTable.from_json(Path(argv[1]).read_text(encoding="utf-8"))
    regions = argv[2:] or [d for d in table.detectors if d.startswith("substrate")]
    sel = [table.detectors.index(r) for r in regions]
    print(f"field: {field.source.get('stack', '?')}")
    print(f"table: {table.source.get('stack', '?')}; SIL takes {regions}")

    free = uniform_field_map(irr_per_mw, dye.excitation_um)  # same irradiance, no gold
    print("h/nm  excitation/free  q      eta per region                  detected 1/s per dye")
    rates = {}
    for h in HEIGHTS_NM:
        on = EmitterEnsemble.isotropic(h * 1e-3, 300, 1)
        r_exc = excitation_rate_s(field.e_at(on.offsets_um, power_mw), on.mu_abs, dye, 1.333)
        r_free = excitation_rate_s(free.e_at(on.offsets_um, power_mw), on.mu_abs, dye, 1.333)
        q = table.quantum_yield(on.offsets_um, on.mu_em, dye.quantum_yield)
        eta = table.collection(on.offsets_um, on.mu_em).mean(axis=0)
        det = detected_rate_s([0, 0, 0], on, field, table, dye, power_mw)
        rates[h] = float(det[sel].sum())
        print(
            f"{h:4.0f}  {r_exc.mean() / r_free.mean():15.3f}  {q.mean():.3f}  "
            f"{np.array2string(eta, precision=3):30s}  {rates[h]:8.2f}"
        )

    sites, other_bg = 1e4, 0.2  # aptamers/µm²; gold PL, lacquer, glass photons/(µm² s) (P)
    for h in (HEIGHTS_NM[0], HEIGHTS_NM[2]):
        rate_on = rates[h]
        print(f"\nh = {h:.0f} nm: N_on, off fraction, t_min for P_det >= 0.95 at 0.05/mm²")
        for n_on in (20, 200, 2500):
            for off in (0.001, 0.01, 0.05):
                t_min = math.nan
                for t in (0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300):
                    cam = Camera(2.44, 0.8, 1.5, 0.2, t, 1.2)  # 5 mm on 2048 px (P)
                    d = spot_detection(
                        cam, n_on * rate_on * t, (sites * off * rate_on + other_bg) * t, 0.05
                    )
                    if d.detection_probability >= 0.95:
                        t_min = t
                        break
                print(f"{n_on:5d}  {off:5.3f}  {t_min:7.2f} s")


if __name__ == "__main__":
    main(sys.argv[1:])
