"""Smoke tests for v2.1 additions: PPS baseline, dMOP3 problem,
Friedman/Nemenyi statistics, transfer strategy interface, adaptive transfer,
multi-source knowledge fusion."""
import numpy as np

from algorithms import PPSNSGA2, LDTNSGA2, MomentTransferStrategy
from problems import DMOP3
from problems.true_fronts import get_true_pf
from evaluation.friedman_test import friedman_nemenyi


def test_pps_smoke_fda1():
    """PPS baseline runs on FDA1 without errors and returns finite metrics."""
    from algorithms.nsga2_core import DynamicDriver
    from problems import FDA1
    prob = FDA1(tau_t=10, n_t=5)
    algo = PPSNSGA2()
    driver = DynamicDriver(algo, prob, "FDA1", n_envs=2, pop_size=20,
                           run_seed=42)
    rows = driver.run(n_points=100)
    assert len(rows) > 0
    for r in rows:
        assert np.isfinite(r["IGD"])
        assert np.isfinite(r["HV"])


def test_dmop3_dimensions():
    """dMOP3 returns 2 objectives and n_var=10."""
    prob = DMOP3(tau_t=10, n_t=5)
    X = np.random.rand(5, 10)
    F = prob._calc_objectives(X)
    assert F.shape == (5, 2)
    assert np.all(np.isfinite(F))


def test_dmop3_true_front():
    """dMOP3 true PF has correct shape and is time-invariant."""
    pf = get_true_pf("DMOP3", n_points=200)
    assert pf.shape == (200, 2)
    assert np.all(np.isfinite(pf))
    # time-invariant: same result regardless of t
    pf2 = get_true_pf("DMOP3", n_points=200, t=0.5)
    assert np.allclose(pf, pf2)


def test_friedman_nemenyi():
    """Friedman test detects a clear ranking difference."""
    rng = np.random.default_rng(0)
    # 3 algorithms, 20 seeds; algo 0 clearly best (low values)
    data = np.column_stack([
        rng.normal(1.0, 0.1, 20),
        rng.normal(2.0, 0.1, 20),
        rng.normal(3.0, 0.1, 20),
    ])
    res = friedman_nemenyi(data, higher_is_better=False)
    assert res["friedman"]["p_value"] < 0.05
    assert res["nemenyi"]["critical_difference"] > 0
    assert len(res["nemenyi"]["significant_pairs"]) >= 1


def test_friedman_no_difference():
    """Friedman test does not reject when algorithms are equivalent."""
    rng = np.random.default_rng(1)
    data = rng.normal(1.0, 1.0, size=(30, 3))
    res = friedman_nemenyi(data, higher_is_better=False)
    assert res["friedman"]["p_value"] > 0.05


def test_transfer_strategy_interface():
    """MomentTransferStrategy produces a valid population."""
    from problems import FDA1
    prob = FDA1(tau_t=10, n_t=5)
    rng = np.random.default_rng(7)
    N, D = 30, 10
    prev_X = rng.uniform(prob.xl, prob.xu, size=(N, D))
    prev_F = prob._calc_objectives(prev_X)
    strat = MomentTransferStrategy()
    new_X, info = strat.transfer(prev_X, prev_F, prob, 1, rng)
    assert new_X.shape == (N, D)
    assert np.all(np.isfinite(new_X))
    assert "boundary_clip_ratio" in info


def test_ldt_adaptive_transfer():
    """LDT with adaptive_transfer=True runs and reports sigma_scale."""
    from algorithms.nsga2_core import DynamicDriver
    from problems import FDA1
    prob = FDA1(tau_t=10, n_t=5)
    algo = LDTNSGA2(adaptive_transfer=True, adaptive_beta=1.0)
    driver = DynamicDriver(algo, prob, "FDA1", n_envs=2, pop_size=20,
                           run_seed=42)
    rows = driver.run(n_points=100)
    assert len(rows) > 0
    for r in rows:
        assert np.isfinite(r["IGD"])


def test_ldt_multi_source():
    """LDT with multi_source=True runs and fuses history."""
    from algorithms.nsga2_core import DynamicDriver
    from problems import FDA1
    prob = FDA1(tau_t=10, n_t=5)
    algo = LDTNSGA2(multi_source=True, history_size=3)
    driver = DynamicDriver(algo, prob, "FDA1", n_envs=3, pop_size=20,
                           run_seed=42)
    rows = driver.run(n_points=100)
    assert len(rows) > 0
    for r in rows:
        assert np.isfinite(r["IGD"])


def test_ldt_pluggable_strategy():
    """LDT delegates to a TransferStrategy when provided."""
    from problems import FDA1
    prob = FDA1(tau_t=10, n_t=5)
    rng = np.random.default_rng(11)
    N, D = 30, 10
    prev_X = rng.uniform(prob.xl, prob.xu, size=(N, D))
    prev_F = prob._calc_objectives(prev_X)
    algo = LDTNSGA2(transfer_strategy=MomentTransferStrategy())
    algo.establish_baseline(prob)     # v3.0 (plan 2.6): need a baseline first
    prob.set_generation(10)           # environment changes
    new_X, info = algo.respond(prev_X, prev_F, prob, 1, rng)
    assert new_X.shape == (N, D)
    assert np.all(np.isfinite(new_X))
