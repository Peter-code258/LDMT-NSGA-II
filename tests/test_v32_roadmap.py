"""v3.2 regression tests: plan §7 roadmap features.

Covers:
  plan 3.2  full-covariance and low-rank sampling (and diag invariance)
  plan 3.3  negative-transfer gating, budget neutrality
  plan 3.4  raw/fused centre state separation
  plan 4.2  extended change patterns (irregular / reversed / gradual / noisy)
  plan 5.3  YAML as the single configuration entry point
"""
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)
sys.path.insert(0, PROJ)

from algorithms.ldmt_nsga2 import (LDTNSGA2, LDTNSGA2Cov, LDTNSGA2LowRank,  # noqa: E402
                                   LDTNSGA2Gated, LDTNSGA2Full)
from experiments.common import RunConfig, run_single, build_algorithm  # noqa: E402
from experiments.config import load_config, available_configs, ConfigError  # noqa: E402
from problems import EXTENDED_PROBLEMS, EXTENDED_LABELS, get_true_pf  # noqa: E402


# ------------------------------------------------- plan 3.2: covariance modes ----
def test_diag_sampling_matches_v2x_formula():
    """cov_mode='diag' must be exactly mu + z*sigma (the frozen behaviour)."""
    algo = LDTNSGA2()
    mu = np.array([0.5] * 4)
    sigma = np.array([0.1] * 4)
    rng = np.random.default_rng(0)
    got = algo.sample_transfer(mu, sigma, np.diag(sigma ** 2), 20, rng)
    ref_rng = np.random.default_rng(0)
    noise = ref_rng.normal(0.0, 1.0, size=(20, 4))
    ref = np.clip(mu[None, :] + noise * sigma[None, :], 0.0, 1.0)
    assert np.allclose(got, ref)


def test_cov_estimate_is_psd_and_records_condition():
    algo = LDTNSGA2Cov()
    rng = np.random.default_rng(1)
    # correlated data on a 1-D manifold embedded in 4-D
    base = rng.normal(size=(40, 1))
    Z = np.clip(0.5 + base @ np.array([[1.0, 1.0, 1.0, 1.0]]) * 0.1, 0, 1)
    cov, cond = algo.estimate_covariance(Z)
    assert cov.shape == (4, 4)
    assert np.allclose(cov, cov.T)
    assert np.all(np.linalg.eigvalsh(cov) > 0), "covariance must be PD"
    assert cond >= 1.0


def test_full_and_lowrank_stay_in_box_and_are_finite():
    Z = np.random.default_rng(2).uniform(0, 1, size=(40, 6))
    mu = Z.mean(axis=0)
    sigma = Z.std(axis=0)
    for algo in (LDTNSGA2Cov(), LDTNSGA2LowRank(low_rank=2)):
        cov, _ = algo.estimate_covariance(Z)
        out = algo.sample_transfer(mu, sigma, cov, 50,
                                   np.random.default_rng(3))
        assert out.shape == (50, 6)
        assert np.isfinite(out).all()
        assert (out >= -1e-12).all() and (out <= 1 + 1e-12).all()


def test_lowrank_captures_structure_that_diag_cannot():
    """On genuinely low-rank data, k must concentrate the sample variance.

    Constructs a set lying on a 2-D subspace of R^5 with strongly unequal
    singular values, so a rank-2 sampler must explain almost all the variance --
    whereas a diagonal sampler would smear it across all 5 axes.
    """
    rng = np.random.default_rng(4)
    basis = rng.normal(size=(2, 5))
    basis /= np.linalg.norm(basis, axis=1, keepdims=True)
    latent = rng.normal(size=(200, 2)) * np.array([1.0, 0.15])
    Z = np.clip(0.5 + latent @ basis, 0.0, 1.0)

    algo = LDTNSGA2LowRank(low_rank=2)
    cov, _ = algo.estimate_covariance(Z)
    out = algo.sample_transfer(Z.mean(axis=0), Z.std(axis=0), cov, 2000,
                               np.random.default_rng(5))
    w = np.maximum(np.linalg.eigvalsh(np.cov(out, rowvar=False))[::-1], 1e-18)
    top2_lowrank = float(np.sum(w[:2]) / np.sum(w))

    diag = LDTNSGA2()
    out_d = diag.sample_transfer(Z.mean(axis=0), Z.std(axis=0), None, 2000,
                                 np.random.default_rng(5))
    wd = np.maximum(np.linalg.eigvalsh(np.cov(out_d, rowvar=False))[::-1], 1e-18)
    top2_diag = float(np.sum(wd[:2]) / np.sum(wd))

    assert top2_lowrank > 0.90, f"low-rank only explained {top2_lowrank:.3f}"
    assert top2_lowrank > top2_diag, (top2_lowrank, top2_diag)


def test_lowrank_residual_is_nonzero_when_k_lt_d():
    """The plan requires a residual so the sampler cannot collapse."""
    algo = LDTNSGA2LowRank(low_rank=1)
    Z = np.random.default_rng(9).uniform(0, 1, size=(60, 6))
    cov, _ = algo.estimate_covariance(Z)
    out = algo.sample_transfer(Z.mean(axis=0), Z.std(axis=0), cov, 300,
                               np.random.default_rng(10))
    w = np.linalg.eigvalsh(np.cov(out, rowvar=False))
    assert w.min() > 0, "a degenerate subspace would make the smallest eigenvalue 0"


# --------------------------------------------------- plan 3.3: gating ----------
def test_gating_is_budget_neutral():
    """rho + kappa + r must sum to 1 for any residual (plan 3.3)."""
    algo = LDTNSGA2Gated()
    for e in (None, 0.0, 0.01, 0.5, 5.0, np.inf):
        rho, kappa, r = algo.gated_ratios(e)
        assert rho + kappa + r == pytest.approx(1.0), (e, rho, kappa, r)
        assert algo.gate_rho_min - 1e-12 <= rho <= algo.transfer_ratio + 1e-12
        assert r >= -1e-12


def test_gating_shrinks_rho_when_residual_is_large():
    algo = LDTNSGA2Gated()
    rho_small, _, _ = algo.gated_ratios(1e-6)
    rho_big, _, _ = algo.gated_ratios(1e3)
    assert rho_big < rho_small, (rho_small, rho_big)
    # approaches (but need not equal) the floor
    assert rho_big == pytest.approx(algo.gate_rho_min, abs=1e-3)
    assert rho_small == pytest.approx(algo.transfer_ratio, abs=1e-3)


def test_gating_off_reproduces_nominal_ratios():
    algo = LDTNSGA2Gated(gating=False)
    rho, kappa, r = algo.gated_ratios(100.0)
    assert (rho, kappa, r) == (algo.transfer_ratio, algo.retain_ratio,
                               algo.random_ratio)


def test_gated_run_keeps_total_population_size():
    cfg = RunConfig(problem="FDA1", tau_t=25, n_t=20, N=40, n_gen=50,
                    algorithms=["ldt_gated"], code_version="v3.2",
                    experiment_group="t")
    rows = run_single("ldt_gated", "FDA1", cfg, 0)
    assert rows, "no rows"
    assert all(np.isfinite(float(r["IGD"])) for r in rows)


# --------------------------------------------- plan 3.4: fusion state ----------
def test_raw_and_fused_centres_are_separate_state():
    algo = LDTNSGA2(multi_source=True)
    assert hasattr(algo, "_mu_prev_fused")
    assert algo._mu_prev is None
    p = EXTENDED_PROBLEMS["IRREG_FDA1"](tau_t=25, n_t=20)
    rng = np.random.default_rng(0)
    X = rng.uniform(p.xl, p.xu, size=(40, p.n_var))
    F = p._calc_objectives(X)
    algo.establish_baseline(p)
    p.set_generation(25)
    algo.respond(X, F, p, 1, rng)
    assert algo._mu_prev is not None
    assert algo._mu_prev_fused is not None
    assert not np.shares_memory(algo._mu_prev, algo._mu_prev_fused)


# ------------------------------------------- plan 4.2: change patterns --------
def test_extended_problems_registered_and_labelled():
    assert set(EXTENDED_PROBLEMS) == {"IRREG_FDA1", "REV_FDA1", "GRAD_FDA1",
                                      "NOISY_FDA1"}
    assert set(EXTENDED_LABELS) == set(EXTENDED_PROBLEMS)


def test_irregular_breakpoints_are_non_uniform_and_increasing():
    p = EXTENDED_PROBLEMS["IRREG_FDA1"](tau_t=25, n_t=20)
    bp = p._build_breakpoints()
    gaps = np.diff(bp)
    assert (gaps > 0).all(), "breakpoints must strictly increase"
    assert len(set(gaps.tolist())) > 1, "schedule must be irregular"
    assert gaps.min() >= 1


def test_reversed_problem_reverses_direction():
    p = EXTENDED_PROBLEMS["REV_FDA1"](tau_t=25, n_t=20)
    ts = [p.set_generation(g) for g in range(0, 300, 25)]
    assert ts[0] <= ts[1], "starts by increasing"
    assert max(ts) > ts[-1], "must come back down after the pivot"


def test_gradual_problem_moves_within_an_environment():
    p = EXTENDED_PROBLEMS["GRAD_FDA1"](tau_t=25, n_t=20)
    a = p.set_generation(0)
    b = p.set_generation(10)
    assert b > a, "t must drift inside an environment"


def test_noisy_problem_perturbs_objectives_but_not_the_front():
    p = EXTENDED_PROBLEMS["NOISY_FDA1"](tau_t=25, n_t=20)
    X = np.atleast_2d(p.probe)
    # noise lives in _evaluate (what the optimizer observes), not in the
    # analytical _calc_objectives used to build the true front
    out_a, out_b = {}, {}
    p._evaluate(X, out_a)
    p._evaluate(X, out_b)
    assert not np.allclose(out_a["F"], out_b["F"]), \
        "noise must perturb repeated evaluations"
    clean = p._calc_objectives(X)
    assert not np.allclose(clean, out_a["F"]), "noise must actually be applied"
    # the analytical front is untouched
    pf = get_true_pf("NOISY_FDA1", M=2, n_points=50)
    assert pf.shape[1] == 2 and np.isfinite(pf).all()


def test_extended_problems_run_end_to_end():
    for name in EXTENDED_PROBLEMS:
        cfg = RunConfig(problem=name, tau_t=25, n_t=20, N=30, n_gen=25,
                        algorithms=["ldt_nsga2"], code_version="v3.2",
                        experiment_group="t")
        rows = run_single("ldt_nsga2", name, cfg, 0)
        assert rows, f"{name} produced no rows"
        assert all(np.isfinite(float(r["IGD"])) for r in rows)


# ------------------------------------------------ plan 5.3: configs -----------
def test_every_shipped_config_loads_and_resolves():
    names = available_configs()
    assert names, "no configs found"
    for n in names:
        cfg, applied = load_config(os.path.join("configs", n))
        assert cfg["problems"], f"{n}: no problems"
        assert applied == {}
        assert cfg.get("code_version"), f"{n}: missing code_version"


def test_unknown_key_is_rejected(tmp_path):
    bad = tmp_path / "bad.yaml"
    bad.write_text("problems: [FDA1]\nnot_a_real_key: 3\n", encoding="utf-8")
    with pytest.raises(ConfigError):
        load_config(str(bad))


def test_overrides_are_reported():
    cfg, applied = load_config("configs/main_experiment.yaml",
                               {"taus": "25", "N": 50})
    assert cfg["taus"] == [25]
    assert cfg["N"] == 50
    assert "taus" in applied and "N" in applied


def test_stale_key_names_are_no_longer_accepted(tmp_path):
    """The old file used tau_t_values/population_size/out_dir -- now an error."""
    stale = tmp_path / "stale.yaml"
    stale.write_text("problems: [FDA1]\ntau_t_values: [10, 25]\n", encoding="utf-8")
    with pytest.raises(ConfigError):
        load_config(str(stale))


# --------------------------------------------------------- roadmap bundling ----
def test_v32_bundle_enables_the_documented_mechanisms():
    algo = LDTNSGA2Full()
    assert algo.adaptive_alpha is True
    assert algo.gating is True
    assert algo.multi_source is True


def test_new_variants_are_not_in_canonical_list():
    from experiments.common import CANONICAL_ALGOS
    for v in ("ldt_cov", "ldt_lowrank", "ldt_gated", "ldt_v32"):
        assert v not in CANONICAL_ALGOS, f"{v} must stay out of the default list"
        assert build_algorithm(v) is not None
