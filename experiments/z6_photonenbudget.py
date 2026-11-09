"""Z6: photon budget of one captured cell against the background of the spot.

Explorative script (not part of the tests); results in docs/ergebnisse-z6.md.
All parameters are placeholders (P).
"""

import math

from bindoptix.signal import (
    Camera,
    EmitterEnsemble,
    Fluorophore,
    detected_rate_s,
    excitation_rate_s,
    placeholder_response_table,
    spot_detection,
    uniform_field_map,
)


def main() -> None:
    dye = Fluorophore(250_000.0, 0.27, 0.640, 0.670, bleach_yield=1e-6)  # Cy5-like (P)
    power_mw, field_mm = 10.0, 5.0
    irr_per_mw = 1e-3 / (math.pi * (0.5 * field_mm * 1e-3) ** 2)  # W/m² per mW
    field = uniform_field_map(irr_per_mw, dye.excitation_um)
    table = placeholder_response_table(0.3, dye.emission_um)  # collection 30 % (P)
    on = EmitterEnsemble.isotropic(0.01, 300, 1)
    rate_on = float(detected_rate_s([0, 0, 0], on, field, table, dye, power_mw)[0])
    r_exc = float(
        excitation_rate_s(field.e_at(on.offsets_um, power_mw), on.mu_abs, dye, 1.333).mean()
    )
    sites, other_bg = 1e4, 0.2  # aptamers/µm²; gold PL, lacquer, glass photons/(µm² s) (P)
    t_bleach_h = 1 / (dye.bleach_yield * r_exc) / 3600
    print(
        f"irradiance {irr_per_mw * power_mw:.0f} W/m², excitation {r_exc:.0f} 1/s per dye, "
        f"detected {rate_on:.1f} 1/s per on-dye, bleaching time {t_bleach_h:.1f} h"
    )
    print("N_on   off    t_min for P_det >= 0.95 at 0.05 false alarms/mm² (shot and read noise)")
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
    main()
