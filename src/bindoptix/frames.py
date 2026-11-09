"""Sensor frame, pose to the RTT frame, mirror rule and gravity (docs/konventionen.md).

Sensor frame: right-handed, sensor surface at z = 0, solution at z > 0, normal +z, lengths in µm.
RTT frame (Raytatouille, fusion): lengths in mm, light along +z of the path in question.

A pose maps sensor coordinates to RTT coordinates, x_rtt = M · (x_sensor · 10⁻³) + t. M is
orthogonal with det M = +1 (rotation) or −1 (rotation combined with a mirror). Illuminating or
observing from the solution side instead of the substrate side is the mirror z → −z, not a
rotation (finding of pilot K4 in RTT-Fusion). Polar vectors (electric dipoles, E) transform with
M, axial vectors (magnetic dipoles, H, B) with det(M)·M, rank-2 tensors of polar vectors with
M·A·Mᵀ.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import ArrayLike, NDArray

from bindoptix.units import MM_PER_UM, STANDARD_GRAVITY_M_S2, UM_PER_M, UM_PER_MM

FloatArray = NDArray[np.float64]

ORTHOGONALITY_TOL = 1e-12


def _as_vec3(v: ArrayLike, name: str) -> FloatArray:
    a = np.asarray(v, dtype=np.float64)
    if a.shape[-1] != 3:
        raise ValueError(f"{name} must have a last axis of length 3, got shape {a.shape}")
    if not np.isfinite(a).all():
        raise ValueError(f"{name} contains non-finite values")
    return a


def _as_orthogonal(m: ArrayLike) -> FloatArray:
    a = np.asarray(m, dtype=np.float64)
    if a.shape != (3, 3):
        raise ValueError(f"matrix must be 3×3, got shape {a.shape}")
    if np.max(np.abs(a @ a.T - np.eye(3))) > ORTHOGONALITY_TOL:
        raise ValueError("matrix is not orthogonal within 1e-12")
    return a


@dataclass(frozen=True)
class Pose:
    """Sensor frame (µm) → RTT frame (mm).

    Attributes:
        matrix: orthogonal 3×3 matrix M, det M = ±1.
        translation_mm: position of the sensor origin in the RTT frame [mm].
    """

    matrix: FloatArray
    translation_mm: FloatArray = field(default_factory=lambda: np.zeros(3))

    def __post_init__(self) -> None:
        object.__setattr__(self, "matrix", _as_orthogonal(self.matrix))
        t = _as_vec3(self.translation_mm, "translation_mm")
        if t.shape != (3,):
            raise ValueError("translation_mm must have shape (3,)")
        object.__setattr__(self, "translation_mm", t)

    # --- construction -------------------------------------------------------------------

    @classmethod
    def identity(cls) -> Pose:
        """Sensor frame and RTT frame coincide (light along +z enters from the substrate)."""
        return cls(np.eye(3))

    @classmethod
    def from_rotation(
        cls, rotation: ArrayLike, translation_mm: ArrayLike = (0.0, 0.0, 0.0)
    ) -> Pose:
        """Proper rotation (det = +1) and translation; a mirror is rejected here on purpose."""
        r = _as_orthogonal(rotation)
        if np.linalg.det(r) < 0.0:
            raise ValueError("det = −1: a mirror is not a rotation; use Pose.mirror_z()")
        return cls(r, np.asarray(translation_mm, dtype=np.float64))

    @classmethod
    def mirror_z(cls) -> Pose:
        """Mirror z → −z: light meets the sensor from the solution side."""
        return cls(np.diag([1.0, 1.0, -1.0]))

    # --- properties -----------------------------------------------------------------------

    @property
    def determinant(self) -> float:
        """det M, +1 for a proper pose, −1 for a pose with mirror."""
        return float(np.sign(np.linalg.det(self.matrix)))

    @property
    def is_proper(self) -> bool:
        return self.determinant > 0.0

    # --- algebra --------------------------------------------------------------------------

    def inverse(self) -> InversePose:
        """RTT frame (mm) → sensor frame (µm)."""
        return InversePose(self)

    def then(self, rtt_motion: Pose) -> Pose:
        """This pose followed by a rigid motion of the RTT frame (both in mm after this pose).

        The second pose's own length scaling is not applied: ``rtt_motion`` is read as a motion
        M₂ x + t₂ in mm.
        """
        m = rtt_motion.matrix @ self.matrix
        t = rtt_motion.matrix @ self.translation_mm + rtt_motion.translation_mm
        return Pose(m, t)

    # --- transformations ------------------------------------------------------------------

    def point_to_rtt(self, p_um: ArrayLike) -> FloatArray:
        """Points in the sensor frame [µm], shape (..., 3) → RTT frame [mm]."""
        p = _as_vec3(p_um, "p_um")
        return p * MM_PER_UM @ self.matrix.T + self.translation_mm

    def polar_to_rtt(self, v: ArrayLike) -> FloatArray:
        """Polar vectors (electric dipole, E), shape (..., 3); unit unchanged."""
        return _as_vec3(v, "v") @ self.matrix.T

    def axial_to_rtt(self, v: ArrayLike) -> FloatArray:
        """Axial vectors (magnetic dipole, H, B), shape (..., 3); unit unchanged."""
        return self.determinant * (_as_vec3(v, "v") @ self.matrix.T)

    def tensor_to_rtt(self, a: ArrayLike) -> FloatArray:
        """Rank-2 tensors of polar vectors (e.g. emitter response tensors), shape (..., 3, 3)."""
        t = np.asarray(a, dtype=np.float64)
        if t.shape[-2:] != (3, 3):
            raise ValueError(f"tensor must have trailing shape (3, 3), got {t.shape}")
        return self.matrix @ t @ self.matrix.T


@dataclass(frozen=True)
class InversePose:
    """Inverse of a :class:`Pose`: RTT frame (mm) → sensor frame (µm)."""

    pose: Pose

    def point_to_sensor(self, p_mm: ArrayLike) -> FloatArray:
        """Points in the RTT frame [mm], shape (..., 3) → sensor frame [µm]."""
        p = _as_vec3(p_mm, "p_mm")
        return (p - self.pose.translation_mm) @ self.pose.matrix * UM_PER_MM

    def polar_to_sensor(self, v: ArrayLike) -> FloatArray:
        return _as_vec3(v, "v") @ self.pose.matrix

    def axial_to_sensor(self, v: ArrayLike) -> FloatArray:
        return self.pose.determinant * (_as_vec3(v, "v") @ self.pose.matrix)

    def tensor_to_sensor(self, a: ArrayLike) -> FloatArray:
        t = np.asarray(a, dtype=np.float64)
        if t.shape[-2:] != (3, 3):
            raise ValueError(f"tensor must have trailing shape (3, 3), got {t.shape}")
        return self.pose.matrix.T @ t @ self.pose.matrix


@dataclass(frozen=True)
class Mount:
    """Mounting of the sensor: direction of gravity in the sensor frame (never implicit).

    Attributes:
        gravity_direction: unit vector of gravity in the sensor frame.
        g_m_s2: magnitude of the acceleration of gravity [m/s²].
    """

    gravity_direction: FloatArray
    g_m_s2: float = STANDARD_GRAVITY_M_S2

    def __post_init__(self) -> None:
        d = _as_vec3(self.gravity_direction, "gravity_direction")
        if d.shape != (3,):
            raise ValueError("gravity_direction must have shape (3,)")
        n = float(np.linalg.norm(d))
        if n == 0.0:
            raise ValueError("gravity_direction must not be zero")
        if not np.isfinite(self.g_m_s2) or self.g_m_s2 < 0.0:
            raise ValueError("g_m_s2 must be finite and non-negative")
        object.__setattr__(self, "gravity_direction", d / n)

    @classmethod
    def sensor_floor(cls) -> Mount:
        """Horizontal channel, sensor surface is the floor: cells sediment onto the sensor."""
        return cls(np.array([0.0, 0.0, -1.0]))

    @classmethod
    def sensor_ceiling(cls) -> Mount:
        """Horizontal channel, sensor surface is the ceiling: cells sediment away from it."""
        return cls(np.array([0.0, 0.0, 1.0]))

    @property
    def gravity_um_s2(self) -> FloatArray:
        """Acceleration of gravity in the sensor frame [µm/s²]."""
        return self.gravity_direction * self.g_m_s2 * UM_PER_M

    @property
    def sedimentation_towards_sensor(self) -> float:
        """Component of the gravity direction towards the sensor surface (−z), in [−1, 1]."""
        return float(-self.gravity_direction[2])
