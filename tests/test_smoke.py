"""Smoke tests: each algorithm runs FDA1 for 2 environments (report 6.5).

Report requirement: each algorithm must complete at least FDA1, N=20,
2 environments, 100 evaluations per environment end-to-end smoke test.
"""
import numpy as np
import pytest

from problems import FDA3, FDA1, FDA4
from algorithms import (LDTNSGA2, ARNSGA2, GMMTransfer, DNSGA2Response,
                         DynamicDriver)

pytestmark = pytest.mark.smoke


def _run_smoke(algo, problem_name, problem_cls, tau_t=5, n_envs=2, N=20):
    p = problem_cls(tau_t=tau_t, n_t=20)
    drv = DynamicDriver(algo, p, problem_name, n_envs=n_envs,
                         pop_size=N, run_seed=0)
    rows = drv.run(n_points=100)
    assert len(rows) == n_envs * tau_t
    igds = [r["IGD"] for r in rows]
    hvs = [r["HV"] for r in rows]
    assert np.all(np.isfinite(igds)), "IGD must be finite"
    assert np.all(np.isfinite(hvs)), "HV must be finite"
    assert all(v >= 0 for v in igds)
    return rows


def test_smoke_nsga2_fda1():
    _run_smoke(DNSGA2Response(), "FDA1", FDA1)


def test_smoke_dnsga2_fda1():
    _run_smoke(DNSGA2Response(), "FDA1", FDA1)


def test_smoke_ar1_fda1():
    _run_smoke(ARNSGA2(), "FDA1", FDA1)


def test_smoke_gmm_fda1():
    _run_smoke(GMMTransfer(), "FDA1", FDA1)


def test_smoke_ldt_fda1():
    _run_smoke(LDTNSGA2(), "FDA1", FDA1)


def test_smoke_ldt_fda3():
    _run_smoke(LDTNSGA2(), "FDA3", FDA3)


def test_smoke_ldt_fda4():
    _run_smoke(LDTNSGA2(), "FDA4", FDA4)
