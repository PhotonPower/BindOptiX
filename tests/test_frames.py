"""T1 (pose round trip) and T2 (mirror rule), docs/testfaelle.md; gravity of the mount."""

import numpy as np
import pytest
from numpy.typing import NDArray

from bindoptix.frames import Mount, Pose

RNG_SEED = 20261109


def rodrigues(axis: NDArray[np.float64], angle: float) -> NDArray[np.float64]:
    """Independent rotation matrix (Rodrigues' formula), not from the code under test."""
    k = axis / np.linalg.norm(axis)
    kx = np.array([[0.0, -k[2], k[1]], [k[2], 0.0, -k[0]], [-k[1], k[0], 0.0]])
    r: NDArray[np.float64] = np.eye(3) + np.sin(angle) * kx + (1.0 - np.cos(angle)) * kx @ kx
    return r


def random_poses(n: int) -> list[Pose]:
    rng = np.random.default_rng(RNG_SEED)
    poses = []
    for i in range(n):
        r = rodrigues(rng.normal(size=3), float(rng.uniform(0.0, 2.0 * np.pi)))
        p = Pose.from_rotation(r, rng.normal(scale=10.0, size=3))
        poses.append(p if i % 2 == 0 else p.then(Pose.mirror_z()))
    return poses


# --- T1: round trip ------------------------------------------------------------------------


@pytest.mark.parametrize("pose", random_poses(8))
def test_t1_point_round_trip(pose: Pose) -> None:
    rng = np.random.default_rng(RNG_SEED + 1)
    p_um = rng.normal(scale=500.0, size=(100, 3))
    back = pose.inverse().point_to_sensor(pose.point_to_rtt(p_um))
    assert np.max(np.abs(back - p_um)) <= 1e-12 * 500.0 * 10


@pytest.mark.parametrize("pose", random_poses(8))
def test_t1_vector_and_tensor_round_trip(pose: Pose) -> None:
    rng = np.random.default_rng(RNG_SEED + 2)
    v = rng.normal(size=(50, 3))
    a = rng.normal(size=(50, 3, 3))
    inv = pose.inverse()
    assert np.max(np.abs(inv.polar_to_sensor(pose.polar_to_rtt(v)) - v)) < 1e-12
    assert np.max(np.abs(inv.axial_to_sensor(pose.axial_to_rtt(v)) - v)) < 1e-12
    assert np.max(np.abs(inv.tensor_to_sensor(pose.tensor_to_rtt(a)) - a)) < 1e-12


def test_t1_units_mm_and_um() -> None:
    p = Pose.from_rotation(np.eye(3), (1.0, 2.0, 3.0))
    np.testing.assert_allclose(p.point_to_rtt([1000.0, 0.0, -500.0]), [2.0, 2.0, 2.5], atol=1e-15)


# --- T2: mirror rule -----------------------------------------------------------------------


def test_t2_mirror_is_not_a_rotation() -> None:
    m = Pose.mirror_z()
    assert m.determinant == -1.0
    assert not m.is_proper
    with pytest.raises(ValueError, match="mirror"):
        Pose.from_rotation(np.diag([1.0, 1.0, -1.0]))


def test_t2_mirror_on_points_and_dipoles() -> None:
    m = Pose.mirror_z()
    np.testing.assert_array_equal(m.point_to_rtt([0.0, 0.0, 1000.0]), [0.0, 0.0, -1.0])
    # electric dipole (polar): normal component flips, tangential components stay
    np.testing.assert_array_equal(m.polar_to_rtt([1.0, 2.0, 3.0]), [1.0, 2.0, -3.0])
    # magnetic dipole (axial): reflection in z keeps the z component, flips x and y
    np.testing.assert_array_equal(m.axial_to_rtt([1.0, 2.0, 3.0]), [-1.0, -2.0, 3.0])


@pytest.mark.parametrize("pose", random_poses(6))
def test_t2_cross_product_is_axial(pose: Pose) -> None:
    # Independent reference: a × b is axial, so T(a × b) = T(a) × T(b) for any orthogonal T.
    rng = np.random.default_rng(RNG_SEED + 3)
    a, b = rng.normal(size=(2, 20, 3))
    lhs = pose.axial_to_rtt(np.cross(a, b))
    rhs = np.cross(pose.polar_to_rtt(a), pose.polar_to_rtt(b))
    assert np.max(np.abs(lhs - rhs)) < 1e-12


@pytest.mark.parametrize("pose", random_poses(6))
def test_t2_quadratic_form_invariant(pose: Pose) -> None:
    # Response tensors (ADR draft, data ports): μᵀ A μ must not depend on the frame.
    rng = np.random.default_rng(RNG_SEED + 4)
    a = rng.normal(size=(3, 3))
    a = a + a.T
    mu = rng.normal(size=(10, 3))
    q_sensor = np.einsum("ni,ij,nj->n", mu, a, mu)
    mu_r, a_r = pose.polar_to_rtt(mu), pose.tensor_to_rtt(a)
    q_rtt = np.einsum("ni,ij,nj->n", mu_r, a_r, mu_r)
    assert np.max(np.abs(q_sensor - q_rtt)) < 1e-11


def test_t2_mirror_twice_is_identity() -> None:
    twice = Pose.mirror_z().then(Pose.mirror_z())
    assert twice.is_proper
    np.testing.assert_array_equal(twice.matrix, np.eye(3))


# --- mount ---------------------------------------------------------------------------------


def test_mount_floor_and_ceiling() -> None:
    floor, ceiling = Mount.sensor_floor(), Mount.sensor_ceiling()
    assert floor.sedimentation_towards_sensor == 1.0
    assert ceiling.sedimentation_towards_sensor == -1.0
    np.testing.assert_allclose(floor.gravity_um_s2, [0.0, 0.0, -9.80665e6], rtol=1e-15)


def test_mount_tilted_is_normalised() -> None:
    m = Mount(np.array([1.0, 0.0, -1.0]))
    assert np.linalg.norm(m.gravity_direction) == pytest.approx(1.0, rel=1e-15)
    assert m.sedimentation_towards_sensor == pytest.approx(np.sqrt(0.5), rel=1e-15)


def test_mount_rejects_zero_direction() -> None:
    with pytest.raises(ValueError):
        Mount(np.zeros(3))
