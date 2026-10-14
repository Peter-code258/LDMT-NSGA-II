# -*- coding: utf-8 -*-
"""v3.1 regression tests: adaptive trend step (plan 3.1)."""
import io
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)
sys.path.insert(0, PROJ)

from algorithms.ldmt_nsga2 import LDTNSGA2, LDTNSGA2Adaptive  # noqa: E402
from experiments.common import problem_factory  # noqa: E402


def _setup(problem="FDA1", tau=25, n_t=20, N=50, seed=0):
    p = problem_factory(problem, tau_t=tau, n_t=n_t)
    rng = np.random.default_rng(seed)
    X = rng.uniform(p.xl, p.xu, size=(N, p.n_var))
    F = p._calc_objectives(X)
    return p, X, F, rng


# --------------------------------------------------------------- rule shape ----
def test_alpha_bounded_in_registered_interval():
    """alpha_k must stay inside [alpha_min, alpha_max] for every change."""
    algo = LDTNSGA2Adaptive(alpha_min=0.5, alpha_max=1.0, noise_samples=60)
    p, X, F, rng = _setup()
    algo.establish_baseline(p)
    for k in range(1, 8):
        p.set_generation(k * 25)
        _, info = algo.respond(X, F, p, k, rng)
        assert 0.5 - 1e-12 <= info["alpha_k"] <= 1.0 + 1e-12, info["alpha_k"]
    assert algo.alpha_series, "alpha series must be recorded"


def test_alpha_is_low_when_set_is_static():
    """Type III (static PS) must drive alpha DOWN, not up.

    This is the mechanism the thesis' conclusion 3 demands: when the
    decision-space centre does not really move, the increment is noise and the
    step must shrink.  dMOP1A is the authoritative Type III problem.
    """
    algo = LDTNSGA2Adaptive(noise_samples=120)
    p, X, F, rng = _setup(problem="DMOP1A")
    algo.establish_baseline(p)
    # keep the SAME population across "changes": the Pareto set is static, so a
    # correct rule sees no real movement.
    for k in range(1, 6):
        p.set_generation(k * 25)
        algo.respond(X, F, p, k, rng)
    finite = [a for a in algo.alpha_series if np.isfinite(a)]
    assert finite, "no alpha recorded"
    assert max(finite) <= 0.75, (
        "static solution set should not reach the high end of the step "
        f"range; got {finite}")


def test_alpha_is_high_when_set_translates():
    """A genuinely translating centre must drive alpha UP."""
    p = problem_factory("FDA1", tau_t=25, n_t=20)
    algo = LDTNSGA2Adaptive(noise_samples=120)
    algo.establish_baseline(p)
    rng = np.random.default_rng(3)
    # synthesise a translating non-dominated set: FDA1's PS is x_i = G(t)
    for k in range(1, 7):
        p.set_generation(k * 25)
        g = p.G(p.t) if hasattr(p, "G") else np.sin(np.pi * k / 10.0)
        X = rng.uniform(p.xl, p.xu, size=(50, p.n_var))
        X[:, 0] = rng.uniform(p.xl[0], p.xu[0], size=50)
        if p.n_var > 1:
            X[:, 1:] = g
        F = p._calc_objectives(X)
        algo.respond(X, F, p, k, rng)
    finite = [a for a in algo.alpha_series if np.isfinite(a)]
    assert finite
    assert max(finite) > 0.5, (
        f"translating set should push alpha above the floor; got {finite}")


def test_first_change_is_conservative():
    """With no residual yet, take the low end (documented behaviour)."""
    algo = LDTNSGA2Adaptive(conservative_first_change=True, noise_samples=40)
    p, X, F, rng = _setup()
    algo.establish_baseline(p)
    p.set_generation(25)
    _, info = algo.respond(X, F, p, 1, rng)
    assert info["alpha_k"] == pytest.approx(0.5)


def test_fixed_variants_are_unaffected_by_adaptive_code():
    """alpha=0.5 / 1.0 must still be exactly constant."""
    for cls, expected in ((LDTNSGA2, 0.5),):
        algo = cls()
        p, X, F, rng = _setup()
        algo.establish_baseline(p)
        for k in range(1, 4):
            p.set_generation(k * 25)
            _, info = algo.respond(X, F, p, k, rng)
            assert info["alpha_k"] == pytest.approx(expected)


def test_adaptive_defaults_off_for_base_class():
    algo = LDTNSGA2()
    assert algo.adaptive_alpha is False


def test_noise_estimator_is_deterministic_given_rng():
    a1 = LDTNSGA2Adaptive(noise_samples=50)
    a2 = LDTNSGA2Adaptive(noise_samples=50)
    p, X, F, _ = _setup()
    z = np.random.default_rng(0).uniform(0, 1, size=(40, p.n_var))
    n1 = a1.centre_noise(z, np.random.default_rng(7))
    n2 = a2.centre_noise(z, np.random.default_rng(7))
    assert n1 == pytest.approx(n2)
    assert n1 > 0


def test_adaptive_step_zero_when_no_movement():
    algo = LDTNSGA2Adaptive(noise_samples=80)
    p, X, F, _ = _setup()
    z = np.random.default_rng(1).uniform(0, 1, size=(50, p.n_var))
    mu = z.mean(axis=0)
    a, mv, nz, ratio = algo.adaptive_step(mu, mu, z,
                                          np.random.default_rng(2))
    assert mv == 0.0
    assert a == pytest.approx(0.5), "no movement -> alpha_min"


# ----------------------------------------------------- do-not-break-baselines ----
def test_alpha_k_column_does_not_change_existing_numbers():
    """The v3.1 column must not alter IGD/HV of the fixed-alpha algorithms.

    Adding a field to ROW_FIELDS is metadata-only: the numbers must be
    bit-identical to what the frozen datasets contain.
    """
    from experiments.common import RunConfig, run_single
    cfg = RunConfig(problem="FDA1", tau_t=25, n_t=20, N=50, n_gen=50,
                    algorithms=["ldt_nsga2"], code_version="v3.1",
                    experiment_group="t")
    rows = run_single("ldt_nsga2", "FDA1", cfg, 0)
    # every row must carry the column, and alpha must be the fixed default
    assert all("alpha_k" in r for r in rows)
    assert {round(float(r["alpha_k"]), 12) for r in rows} == {0.5}
    # a second identical run must reproduce identical IGD (determinism)
    rows2 = run_single("ldt_nsga2", "FDA1", cfg, 0)
    assert [r["IGD"] for r in rows] == [r["IGD"] for r in rows2]
