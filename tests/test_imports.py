"""P0-0A: import-completeness test.

Every module under algorithms/, problems/, evaluation/ and experiments/
must import without error, and the documented public API must exist with
stable canonical ids. Intended to run FIRST so a broken import anywhere
fails immediately.
"""
import importlib
import pkgutil

import algorithms
import problems
import evaluation
import experiments


def test_every_module_imports():
    for pkg in (algorithms, problems, evaluation, experiments):
        names = [m.name for m in pkgutil.iter_modules(pkg.__path__)]
        assert names, f"{pkg.__name__} exposes no modules"
        for name in names:
            mod = importlib.import_module(f"{pkg.__name__}.{name}")
            assert mod is not None


def test_algorithms_public_api_and_canonical_ids():
    from algorithms import (DynamicDriver, make_nsga2, EnvironmentManager,
                            ChangeDetector, normalize, denormalize,
                            extract_features, non_dominated_set,
                            MigrationFunction, NSGA2Static, DNSGA2Response,
                            ARNSGA2, GMMTransfer, LDTNSGA2, LDMTNSGA2,
                            PPSNSGA2, TransferStrategy, TRANSFER_REGISTRY)
    # canonical ids are unique and stable (v2.2, P1-3.1)
    ids = [NSGA2Static.name, DNSGA2Response.name, ARNSGA2.name,
           GMMTransfer.name, PPSNSGA2.name, LDTNSGA2.name]
    assert ids == ["nsga2_static", "dnsga2_inspired", "ar1_predictive",
                   "gmm_inspired", "pps", "ldt_nsga2"]
    assert len(set(ids)) == 6
    assert LDMTNSGA2 is LDTNSGA2  # backward-compat alias, same class


def test_problems_public_api_and_shared_probe():
    import numpy as np
    from problems import FDA1, FDA3, FDA4, DMOP1, DMOP2, DMOP3
    from problems.base import DynamicProblem
    for cls in (FDA1, FDA3, FDA4, DMOP1, DMOP2, DMOP3):
        p = cls(tau_t=5, n_t=10)
        assert isinstance(p, DynamicProblem)
        # v2.2 (P1-1E): every problem instance owns a FIXED shared probe
        assert p.probe.shape == (p.PROBE_SIZE, p.n_var)
        assert np.all(p.probe >= p.xl - 1e-12)
        assert np.all(p.probe <= p.xu + 1e-12)
    # identical class -> identical probe across instances/configs
    a, b = FDA1(tau_t=5, n_t=10), FDA1(tau_t=25, n_t=10)
    assert np.array_equal(a.probe, b.probe)


def test_evaluation_public_api():
    from evaluation import statistics, hv, igd, reference_front, friedman_test
    from evaluation.statistics import holm_adjust, wilcoxon_holm, summarize
    from evaluation.friedman_test import friedman_nemenyi
    assert callable(holm_adjust) and callable(wilcoxon_holm)
    assert callable(summarize) and callable(friedman_nemenyi)


def test_experiments_common_api():
    from experiments.common import (RunConfig, ExperimentRunner, run_single,
                                    problem_factory, ROW_FIELDS,
                                    CANONICAL_ALGOS, ALGORITHM_REGISTRY,
                                    LEGACY_ALGO_MAP, ALGO_DISPLAY,
                                    build_algorithm)
    assert CANONICAL_ALGOS == ["nsga2_static", "dnsga2_inspired",
                               "ar1_predictive", "gmm_inspired", "pps",
                               "ldt_nsga2"]
    # every canonical id resolves through the registry and maps back to
    # its own class id (no duplicate counting)
    for cid in CANONICAL_ALGOS:
        algo = build_algorithm(cid)
        assert getattr(algo, "name", None) == cid
    # legacy aliases all map INTO the canonical set
    for legacy, canon in LEGACY_ALGO_MAP.items():
        assert canon in CANONICAL_ALGOS
    # every canonical id has a display name.  Experimental variants
    # (ldt_alpha1, ldt_adaptive) are also displayable but are deliberately NOT
    # in CANONICAL_ALGOS, so the default experiment list is never widened.
    assert set(CANONICAL_ALGOS) <= set(ALGO_DISPLAY)
    assert "ldt_adaptive" in ALGO_DISPLAY
    assert "ldt_alpha1" in ALGO_DISPLAY
    for variant in ("ldt_adaptive", "ldt_alpha1"):
        assert variant not in CANONICAL_ALGOS, \
            f"{variant} must stay out of the default (canonical) list"
