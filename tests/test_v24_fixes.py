"""Regression tests for the v2.4 known-issue (K-list) fixes.

Covers:
  K2  dmop3 dead code / static behaviour (objectives independent of t)
  K3  dmop2 documented variant (still dynamic — PF must move with t)
  K4  MomentTransferStrategy normalized-space transfer
  K5  _hv_exact inclusion-exclusion correctness vs pymoo HV
  K7  Friedman tied average ranks + Nemenyi alpha table selection
  K9  boundary_clip_ratio is a REAL diagnostic (not hard-coded 0.0)
  K10 FDA4 reference front deterministic point count (dedup of 23x23 grid)
"""
import numpy as np
import pytest

from evaluation.hv import _hv_exact
from evaluation.friedman_test import _rank_per_seed, nemenyi_posthoc
from algorithms.transfer_strategy import MomentTransferStrategy
from algorithms.dnsga2_response import NSGA2Static, DNSGA2Response
from algorithms.pps_nsga2 import PPSNSGA2
from algorithms.ar_nsga2 import ARNSGA2
from problems.dmop1 import DMOP1
from problems.dmop2 import DMOP2
from problems.dmop3 import DMOP3
from problems.true_fronts import get_true_pf


# ---------------------------------------------------------------- K5: HV ----
def test_hv_exact_matches_pymoo_2d():
    """_hv_exact must equal pymoo's HV on random small 2-D sets."""
    pymoo_hv = pytest.importorskip("pymoo.indicators.hv")
    rng = np.random.default_rng(42)
    ref = np.array([1.1, 1.1])
    for _ in range(5):
        F = rng.uniform(0.0, 1.0, size=(6, 2))
        expected = float(pymoo_hv.HV(ref_point=ref)(F))
        got = _hv_exact(F, ref)
        assert np.isclose(got, expected, atol=1e-9), \
            f"_hv_exact={got} != pymoo={expected}"


def test_hv_exact_matches_pymoo_3d():
    """_hv_exact must equal pymoo's HV on random small 3-D sets."""
    pymoo_hv = pytest.importorskip("pymoo.indicators.hv")
    rng = np.random.default_rng(7)
    ref = np.array([1.1, 1.1, 1.1])
    for _ in range(3):
        F = rng.uniform(0.0, 1.0, size=(5, 3))
        expected = float(pymoo_hv.HV(ref_point=ref)(F))
        got = _hv_exact(F, ref)
        assert np.isclose(got, expected, atol=1e-9), \
            f"_hv_exact={got} != pymoo={expected}"


def test_hv_exact_single_box():
    # single point (0.5, 0.5), ref (1,1): HV = 0.25
    assert np.isclose(_hv_exact(np.array([[0.5, 0.5]]), np.array([1.0, 1.0])),
                      0.25, atol=1e-12)


# ------------------------------------------------------------- K7: Friedman ----
def test_friedman_tied_average_ranks():
    """Ties must get the AVERAGE rank (1.5, 1.5, 3), not arbitrary distinct."""
    m = np.array([[1.0, 1.0, 2.0]])  # one seed, three algorithms, tie for best
    ranks = _rank_per_seed(m, higher_is_better=False)
    assert np.allclose(ranks[0], [1.5, 1.5, 3.0]), f"got {ranks[0]}"


def test_friedman_tied_average_ranks_higher_better():
    m = np.array([[2.0, 2.0, 1.0]])  # higher is better: tie for best
    ranks = _rank_per_seed(m, higher_is_better=True)
    assert np.allclose(ranks[0], [1.5, 1.5, 3.0]), f"got {ranks[0]}"


def test_nemenyi_alpha_tables():
    rng = np.random.default_rng(0)
    rank_matrix = np.tile(np.arange(1, 5, dtype=float), (30, 1))
    cd_10 = nemenyi_posthoc(rank_matrix, alpha=0.10)["critical_difference"]
    cd_05 = nemenyi_posthoc(rank_matrix, alpha=0.05)["critical_difference"]
    # q_0.10 < q_0.05 -> tighter CD at alpha=0.10 (Demšar 2006 Table 5)
    assert cd_10 < cd_05
    # exact value check: k=4, N=30, q0.05=2.569 -> CD = 2.569*sqrt(4*5/180)
    assert np.isclose(cd_05, 2.569 * np.sqrt(4 * 5 / (6 * 30)), atol=1e-9)
    with pytest.raises(ValueError):
        nemenyi_posthoc(rank_matrix, alpha=0.5)


# ------------------------------------------- K4: MomentTransferStrategy ----
class _FakeProblem:
    def __init__(self, D=5):
        self.xl = -2.0 * np.ones(D)
        self.xu = 3.0 * np.ones(D)
        self.n_obj = 2

    def _calc_objectives(self, X):
        X = np.atleast_2d(X)
        return np.column_stack([X[:, 0], 1.0 - X[:, 0] ** 2])


def test_moment_transfer_within_bounds_real_scale():
    """K4: transfer must operate in normalized space and return REAL-scale
    points within [xl, xu] (old code mixed coordinate systems)."""
    rng = np.random.default_rng(1)
    prob = _FakeProblem(D=5)
    N, D = 40, 5
    prev_X = rng.uniform(prob.xl, prob.xu, size=(N, D))
    prev_F = np.column_stack([prev_X[:, 0], 1.0 - prev_X[:, 0] ** 2])
    strat = MomentTransferStrategy()
    for env in range(3):  # repeat so the trend branch (mu_prev) engages
        new_X, info = strat.transfer(prev_X, prev_F, prob, env, rng)
        assert new_X.shape == (N, D)
        assert np.isfinite(new_X).all()
        assert new_X.min() >= prob.xl.min() - 1e-9
        assert new_X.max() <= prob.xu.max() + 1e-9
        assert 0.0 <= info["boundary_clip_ratio"] <= 1.0


def test_moment_transfer_normalized_consistency():
    """The transfer part must centre on the NORMALIZED nd-mean: with xl=0,
    xu=1 the output mean of a large population should sit inside [0,1] and
    track the previous mean, not blow up through double normalization."""
    rng = np.random.default_rng(3)
    prob = _FakeProblem(D=4)
    prob.xl = np.zeros(4)
    prob.xu = np.ones(4)
    N = 200
    prev_X = rng.uniform(0.2, 0.4, size=(N, 4))
    prev_F = np.column_stack([prev_X[:, 0], 1.0 - prev_X[:, 0]])
    strat = MomentTransferStrategy(transfer_ratio=1.0, retain_ratio=0.0,
                                   random_ratio=0.0)
    new_X, _ = strat.transfer(prev_X, prev_F, prob, 1, rng)
    centre = new_X.mean(axis=0)
    assert np.all(centre > 0.05) and np.all(centre < 0.75), \
        f"transfer centre {centre} left the plausible region"


# ------------------------------------------------- K2/K3: DMOP behaviours ----
def test_dmop3_is_static():
    """K2: this DMOP3 is a project-custom STATIC control problem; objectives
    must NOT depend on t (the dead G variable was removed)."""
    rng = np.random.default_rng(5)
    X = rng.uniform(0.0, 1.0, size=(20, 10))
    p = DMOP3(tau_t=10, n_t=10)
    p.set_generation(0)
    F0 = p._calc_objectives(X)
    p.set_generation(57)   # t = floor(57/10)/10 = 0.5
    F1 = p._calc_objectives(X)
    assert np.allclose(F0, F1), "DMOP3 objectives changed with t"


def test_dmop3_front_shape_static():
    pf0 = get_true_pf("DMOP3", n_points=200, t=0.0)
    pf1 = get_true_pf("DMOP3", n_points=200, t=0.7)
    assert np.allclose(pf0, pf1)


def test_dmop1_pf_static_ps_moves():
    """DMOP1 (project variant): PF time-invariant, PS = G(t)."""
    rng = np.random.default_rng(9)
    p = DMOP1(tau_t=10, n_t=10)
    x0 = rng.uniform(0.0, 1.0, size=(15,))
    for g in (0, 50):   # t = 0.0 and 0.5
        p.set_generation(g)
        G = np.sin(0.5 * np.pi * p.t)
        X = np.column_stack([x0, np.full((15, 9), G)])  # on the PS
        F = p._calc_objectives(X)
        # on the PS g=1, so f2 = 1 - f1^2 regardless of t
        assert np.allclose(F[:, 1], 1.0 - F[:, 0] ** 2, atol=1e-9)


def test_dmop2_pf_moves_with_t():
    """DMOP2 (project variant): PF scales with g_pf = 1 + G(t)."""
    p = DMOP2(tau_t=10, n_t=10)
    x0 = np.linspace(0.0, 1.0, 11)
    fronts = []
    for g in (0, 50):
        p.set_generation(g)
        G = np.sin(0.5 * np.pi * p.t)
        X = np.column_stack([x0, np.full((11, 9), G)])
        F = p._calc_objectives(X)
        fronts.append(F[:, 1].copy())
    assert not np.allclose(fronts[0], fronts[1]), \
        "DMOP2 PF should move with t (g_pf = 1+G)"


# --------------------------------------------- K9: boundary_clip_ratio ----
def test_static_clip_ratio_real_value():
    prob = _FakeProblem(D=4)
    prob.xl = np.zeros(4)
    prob.xu = np.ones(4)
    rng = np.random.default_rng(11)
    # interior points -> ratio 0.0
    prev_X = rng.uniform(0.1, 0.9, size=(30, 4))
    prev_F = np.column_stack([prev_X[:, 0], 1 - prev_X[:, 0]])
    _, info = NSGA2Static().respond(prev_X, prev_F, prob, 1, rng)
    assert info["boundary_clip_ratio"] == 0.0
    # half the points ON the boundary -> ratio 0.5
    prev_X2 = prev_X.copy()
    prev_X2[:15, 0] = 0.0
    _, info2 = NSGA2Static().respond(prev_X2, prev_F, prob, 1, rng)
    assert np.isclose(info2["boundary_clip_ratio"], 0.5)


def test_response_baselines_report_real_clip_ratio():
    prob = _FakeProblem(D=4)
    prob.xl = np.zeros(4)
    prob.xu = np.ones(4)
    rng = np.random.default_rng(13)
    prev_X = rng.uniform(0.1, 0.9, size=(40, 4))
    prev_F = np.column_stack([prev_X[:, 0], 1 - prev_X[:, 0]])
    for cls in (DNSGA2Response, PPSNSGA2, ARNSGA2):
        algo = cls()
        new_X, info = algo.respond(prev_X, prev_F, prob, 1, rng)
        b = info["boundary_clip_ratio"]
        assert isinstance(b, float) and 0.0 <= b <= 1.0, \
            f"{cls.__name__} returned boundary_clip_ratio={b!r}"
        # recompute independently and compare
        expected = float(np.mean(np.any(
            (new_X <= prob.xl + 1e-9) | (new_X >= prob.xu - 1e-9), axis=1)))
        assert np.isclose(b, expected), \
            f"{cls.__name__}: reported {b} != recomputed {expected}"


# ------------------------------------------------- K10: FDA4 front size ----
def test_fda4_reference_front_dedup_count():
    """K10: FDA4 reference front = dedup of the 23x23 angular grid -> 507."""
    pf = get_true_pf("FDA4", M=3, n_points=500)
    assert pf.shape == (507, 3), f"expected (507, 3), got {pf.shape}"
    # unit sphere: every point has norm 1
    assert np.allclose(np.linalg.norm(pf, axis=1), 1.0, atol=1e-9)
