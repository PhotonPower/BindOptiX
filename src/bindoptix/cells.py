"""Cell model: rod-shaped bacteria as prolate spheroids (stage Z3, docs/abschaetzung-einfang.md).

Translational friction of a prolate spheroid with semi-axes a ≥ b for motion along and across the
symmetry axis (happel1983), orientation-averaged mobility 1/f = (1/f∥ + 2/f⊥)/3; sedimentation
velocity v_s = Δρ V g / f, diffusion D = k_B T / f. Run-and-tumble motility with run speed v,
mean run time τ and mean cosine α of the turning angle gives D_eff = v²τ/(3(1 − α)) in free space
(lovely1975). Values marked P are placeholders (AGENTS.md rule 9).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from bindoptix.units import BOLTZMANN, STANDARD_GRAVITY_M_S2, UM_PER_M


@dataclass(frozen=True)
class CellModel:
    """Rod-shaped cell as prolate spheroid with orientation-averaged transport properties.

    Attributes:
        width_um: diameter (P: 0.5 µm for Legionella).
        length_um: length, ≥ width (P: 2.5 µm).
        delta_rho_kg_m3: density excess over the medium (P: 80 kg/m³).
        viscosity_pa_s: viscosity of the medium (L: water 20 °C, 1.002 mPa·s).
        temperature_k: temperature (293.15 K).
        swim_speed_um_s: run speed v, 0 for non-motile cells (P).
        run_time_s: mean run time τ (P).
        tumble_cos: mean cosine α of the turning angle (P: 0.33 as for E. coli).
    """

    width_um: float = 0.5
    length_um: float = 2.5
    delta_rho_kg_m3: float = 80.0
    viscosity_pa_s: float = 1.002e-3
    temperature_k: float = 293.15
    swim_speed_um_s: float = 0.0
    run_time_s: float = 1.0
    tumble_cos: float = 0.33

    def __post_init__(self) -> None:
        for name in ("width_um", "length_um", "viscosity_pa_s", "temperature_k", "run_time_s"):
            v = getattr(self, name)
            if not (math.isfinite(v) and v > 0.0):
                raise ValueError(f"{name} must be finite and positive")
        if self.length_um < self.width_um:
            raise ValueError("length_um must not be smaller than width_um (prolate)")
        if not math.isfinite(self.delta_rho_kg_m3):
            raise ValueError("delta_rho_kg_m3 must be finite")
        if not (math.isfinite(self.swim_speed_um_s) and self.swim_speed_um_s >= 0.0):
            raise ValueError("swim_speed_um_s must be finite and non-negative")
        if not -1.0 <= self.tumble_cos < 1.0:
            raise ValueError("tumble_cos must lie in [-1, 1)")

    @property
    def volume_um3(self) -> float:
        a, b = 0.5 * self.length_um, 0.5 * self.width_um
        return 4.0 / 3.0 * math.pi * a * b * b

    @property
    def equivalent_radius_um(self) -> float:
        """Radius of the sphere of equal volume [µm]."""
        return float((3.0 * self.volume_um3 / (4.0 * math.pi)) ** (1.0 / 3.0))

    def friction_axes_kg_s(self) -> tuple[float, float]:
        """Friction coefficients (f∥, f⊥) [kg/s] (happel1983); sphere limit 6πηa."""
        mu = self.viscosity_pa_s
        a = 0.5 * self.length_um / UM_PER_M
        b = 0.5 * self.width_um / UM_PER_M
        if a - b <= 1e-9 * a:
            f = 6.0 * math.pi * mu * a
            return f, f
        e = math.sqrt(1.0 - (b / a) ** 2)
        lg = math.log((1.0 + e) / (1.0 - e))
        f_par = 16.0 * math.pi * mu * a * e**3 / ((1.0 + e * e) * lg - 2.0 * e)
        f_perp = 32.0 * math.pi * mu * a * e**3 / (2.0 * e + (3.0 * e * e - 1.0) * lg)
        return f_par, f_perp

    @property
    def mobility_s_kg(self) -> float:
        """Orientation-averaged mobility (1/f∥ + 2/f⊥)/3 [s/kg]."""
        f_par, f_perp = self.friction_axes_kg_s()
        return (1.0 / f_par + 2.0 / f_perp) / 3.0

    @property
    def diffusion_um2_s(self) -> float:
        """Brownian translational diffusion coefficient [µm²/s]."""
        return BOLTZMANN * self.temperature_k * self.mobility_s_kg * UM_PER_M**2

    def sedimentation_speed_um_s(self, g_m_s2: float = STANDARD_GRAVITY_M_S2) -> float:
        """Sedimentation speed v_s = Δρ V g · mobility [µm/s] (negative: the cell rises)."""
        v_m3 = self.volume_um3 / UM_PER_M**3
        return self.delta_rho_kg_m3 * v_m3 * g_m_s2 * self.mobility_s_kg * UM_PER_M

    @property
    def effective_diffusion_um2_s(self) -> float:
        """D + v²τ/(3(1 − α)): Brownian plus run-and-tumble contribution [µm²/s]."""
        v, tau, alpha = self.swim_speed_um_s, self.run_time_s, self.tumble_cos
        return self.diffusion_um2_s + v * v * tau / (3.0 * (1.0 - alpha))
