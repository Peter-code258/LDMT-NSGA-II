"""v2.2 regression tests (fix plan P0-0B/0C/0D, P1-1E/3.2, P2-4.3).

- FDA4: point-wise objective values against a direct transcription of the
  PlatEMO FDA4.m CalObj formula (non-squared trig, distance slice
  ``PopDec(:, M:end)`` -> ``X[:, M-1:]``), plus deterministic unit-sphere
  reference front (P0-0C / P1-1C).
- Holm: known-values table exercising the running-maximum property that the
  pre-v2.2 local implementation lacked (P1-3.2).
- Baselines: NSGA2Static carries the population unchanged while
  DNSGA2Response injects immigrants (P0-0D behaviour separation).
- LDT kernel: normalized-space roundtrip, in-bounds transfer candidates,
  adaptive sigma-scale cap [0.5, 2.0] (P0-0B / P2-4.3).
- Shared fixed probe attached to every problem instance (P1-1E).
"""
import numpy as np
import pytest

from problems import FDA1, FDA4
from problems.true_fronts import get_true_pf, _grid_unit_simplex
from algorithms import NSGA2Static, DNSGA2Response, LDTNSGA2
from algorithms.distribution_features import (
    normalize, denormalize, extract_features)
from evaluation.statistics import holm_adjust


# ---------- P0-0C: FDA4 point-wise vs PlatEMO FDA4.m ----------

def _fda4_reference_row(x, t, M=3):
    """One-row transcription of PlatEMO FDA4.m CalObj (no trig squares):
    F = (1+g) * fliplr(cumprod([1, cos(theta)])) .* [1, sin(theta_reversed)].
    """
    G = abs(np.sin(0.5 * np.pi * t))
    g = float(np.sum((x[M - 1:] - G) ** 2))
    theta = x[:M - 1] * (np.pi / 2.0)
    F = np.empty(M)
    F[0] = (1.0 + g) * np.prod(np.cos(theta))
    for i in range(2, M + 1):
        val = (1.0 + g) * np.sin(theta[M - i])
        if M - i > 0:
            val *= np.prod(np.cos(theta[:M - i]))
        F[i - 1] = val
    return F


def test_fda4_pointwise_matches_platemo():
    rng = np.random.default_rng(7)
    p = FDA4(tau_t=10, n_t=20)
    for g_idx in (0, 3, 10):
        p.set_generation(g_idx)
        X = rng.uniform(0.0, 1.0, size=(40, p.D))
        F = p._calc_objectives(X)
        for k in range(40):
            ref = _fda4_reference_row(X[k], p.t, M=p.M)
            assert np.allclose(F[k], ref, atol=1e-12)


def test_fda4_reference_front_deterministic_unit_sphere():
    pf1 = get_true_pf("FDA4", n_points=100)
    pf2 = get_true_pf("FDA4", n_points=100)
    assert np.array_equal(pf1, pf2)
    # PlatEMO GetOptimum: points normalized onto the unit sphere
    assert np.allclose(np.linalg.norm(pf1, axis=1), 1.0, atol=1e-9)
    a, b = _grid_unit_simplex(3, 100), _grid_unit_simplex(3, 100)
    assert np.array_equal(a, b)              # deterministic, no RNG
    assert np.all(np.diff(a[:, 0]) >= -1e-12)  # sorted by f1


# ---------- P1-3.2: Holm step-down with running maximum ----------

def test_holm_adjust_known_values():
    # sorted 0.01,0.03,0.04 -> 3*0.01=0.03, 2*0.03=0.06; 0.04 must INHERIT
    # the running max 0.06 (the pre-v2.2 code wrongly returned 0.04)
    assert np.allclose(holm_adjust([0.01, 0.04, 0.03]), [0.03, 0.06, 0.06])
    assert np.allclose(holm_adjust([0.001, 0.2, 0.5]), [0.003, 0.4, 0.5])
    assert np.allclose(holm_adjust([0.9, 0.95]), [1.0, 1.0])
    assert np.allclose(holm_adjust([0.004, 0.02, 0.3, 0.04, 0.7]),
                       [0.02, 0.08, 0.6, 0.12, 0.7])
    assert holm_adjust([]) == []


# ---------- P0-0D: static vs responsive baseline separation ----------

def test_static_baseline_carries_population_unchanged():
    rng = np.random.default_rng(0)
    p = FDA1(tau_t=5, n_t=10)
    p.set_generation(5)
    X = rng.uniform(p.xl, p.xu, size=(30, p.n_var))
    F = p._calc_objectives(X)
    new_X, info = NSGA2Static().respond(X, F, p, 1, rng)
    assert np.array_equal(new_X, X), \
        "static baseline must not alter the decision variables"
    assert info["boundary_clip_ratio"] == 0.0


def test_dnsga2_inspired_injects_immigrants():
    rng = np.random.default_rng(0)
    p = FDA1(tau_t=5, n_t=10)
    p.set_generation(5)
    X = rng.uniform(p.xl, p.xu, size=(30, p.n_var))
    F = p._calc_objectives(X)
    new_X, _ = DNSGA2Response().respond(X, F, p, 1, rng)
    assert new_X.shape == X.shape
    assert not np.allclose(new_X, X), \
        "DNSGA-II-inspired must replace ~20% of the population"


# ---------- P0-0B: normalized-space LDT kernel ----------

def test_normalize_denormalize_roundtrip():
    rng = np.random.default_rng(3)
    p = FDA1(tau_t=5, n_t=10)
    X = rng.uniform(p.xl, p.xu, size=(20, p.n_var))
    assert np.allclose(denormalize(normalize(X, p.xl, p.xu), p.xl, p.xu), X)


def test_extract_features_in_unit_space():
    rng = np.random.default_rng(4)
    p = FDA1(tau_t=5, n_t=10)
    X = rng.uniform(p.xl, p.xu, size=(50, p.n_var))
    mu, sigma = extract_features(X, p.xl, p.xu)
    assert np.all(mu >= -1e-12) and np.all(mu <= 1.0 + 1e-12)
    assert np.all(sigma >= 0.0)


def test_ldt_respond_in_bounds_and_shape():
    rng = np.random.default_rng(11)
    p = FDA1(tau_t=5, n_t=10)
    X = rng.uniform(p.xl, p.xu, size=(50, p.n_var))
    F = p._calc_objectives(X)
    new_X, info = LDTNSGA2().respond(X, F, p, 1, rng)
    assert new_X.shape == X.shape
    xl, xu = np.asarray(p.xl), np.asarray(p.xu)
    assert np.all(new_X >= xl - 1e-12) and np.all(new_X <= xu + 1e-12)
    assert "Dk" in info and "sigma_scale" in info
    assert np.isfinite(info["Dk"]) and np.isfinite(info["sigma_scale"])


def test_ldt_adaptive_sigma_scale_capped():
    rng = np.random.default_rng(5)
    p = FDA1(tau_t=5, n_t=10)
    algo = LDTNSGA2(adaptive_transfer=True, adaptive_beta=1e6)
    X = rng.uniform(p.xl, p.xu, size=(50, p.n_var))
    F = p._calc_objectives(X)
    _, info0 = algo.respond(X, F, p, 0, rng)
    p.set_generation(5)  # change environment -> nonzero D_k on 2nd respond
    _, info1 = algo.respond(X, F, p, 1, rng)
    lo, hi = algo.SIGMA_SCALE_MIN, algo.SIGMA_SCALE_MAX
    assert lo <= info0["sigma_scale"] <= hi
    assert lo <= info1["sigma_scale"] <= hi
    assert info1["sigma_scale"] == pytest.approx(hi), \
        "extreme adaptive_beta must saturate at SIGMA_SCALE_MAX, not exceed it"


# ---------- P1-1E: shared fixed probe ----------

def test_shared_fixed_probe():
    a = FDA1(tau_t=5, n_t=10)
    b = FDA1(tau_t=25, n_t=10)
    assert np.array_equal(a.probe, b.probe)  # identical across configs
    xl, xu = np.asarray(a.xl), np.asarray(a.xu)
    assert np.all(a.probe >= xl) and np.all(a.probe <= xu)
    algo = LDTNSGA2()
    algo._ensure_probe(None, a)
    assert np.array_equal(algo._probe, a.probe)  # adopted by algorithms
