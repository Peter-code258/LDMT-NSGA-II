"""Verify environment-change detection and that objectives differ across envs."""
import numpy as np

from problems import FDA1, FDA3, FDA4
from algorithms import EnvironmentManager


def test_environment_manager_is_change():
    em = EnvironmentManager(tau_t=10, n_t=20)
    assert em.is_change(0) is False
    assert em.is_change(5) is False
    assert em.is_change(10) is True
    assert em.is_change(20) is True
    assert em.is_change(11) is False
    assert em.env_index_at(0) == 0
    assert em.env_index_at(9) == 0
    assert em.env_index_at(10) == 1
    assert em.env_index_at(25) == 2
    assert np.isclose(em.time_at(0), 0.0)
    assert np.isclose(em.time_at(10), 1.0 / 20.0)
    assert np.isclose(em.time_at(20), 2.0 / 20.0)


def test_objectives_differ_across_envs():
    rng = np.random.default_rng(0)
    X = rng.uniform(-1, 1, size=(50, 10))
    X[:, 0] = rng.uniform(0, 1, size=50)
    for cls in (FDA1, FDA3):
        p = cls(tau_t=10, n_t=20)
        p.set_generation(0)
        F0 = p._calc_objectives(X)
        p.set_generation(10)  # one env step
        F1 = p._calc_objectives(X)
        assert not np.allclose(F0, F1), f"{cls.__name__} objectives must change"


def test_fda4_objectives_differ_across_envs():
    rng = np.random.default_rng(0)
    X = rng.uniform(0, 1, size=(50, 12))
    p = FDA4(tau_t=10, n_t=20)
    p.set_generation(0)
    F0 = p._calc_objectives(X)
    p.set_generation(10)
    F1 = p._calc_objectives(X)
    assert not np.allclose(F0, F1)


def test_fda4_unit_sphere_at_optimum():
    """When x[:,M:]=G (g=0) the objectives lie on the unit hypersphere."""
    rng = np.random.default_rng(42)
    M = 3
    p = FDA4(tau_t=10, n_t=20)
    p.set_generation(5)  # some interior time
    G = abs(np.sin(0.5 * np.pi * p.t))
    N = 300
    X = rng.uniform(0, 1, size=(N, 12))
    X[:, M:] = G  # force g = 0
    F = p._calc_objectives(X)
    norms = np.linalg.norm(F, axis=1)
    assert np.allclose(norms, 1.0, atol=1e-6), \
        f"FDA4 objectives must be unit-norm at PS, got {norms.min()}-{norms.max()}"
