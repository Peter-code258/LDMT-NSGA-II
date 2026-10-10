"""Smoke tests: each of the 4 algorithms runs FDA3 for 2 environments."""
import numpy as np
import pytest

from problems import FDA3, FDA1, FDA4
from algorithms import (LDMTNSGA2, ARNSGA2, GMMTransfer, DNSGA2Response,
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


def test_smoke_dnsga2_fda3():
    _run_smoke(DNSGA2Response(), "FDA3", FDA3)


def test_smoke_ar_fda3():
    _run_smoke(ARNSGA2(), "FDA3", FDA3)


def test_smoke_gmm_fda3():
    _run_smoke(GMMTransfer(), "FDA3", FDA3)


def test_smoke_ldmt_fda3():
    _run_smoke(LDMTNSGA2(), "FDA3", FDA3)


def test_smoke_ldmt_fda1():
    _run_smoke(LDMTNSGA2(), "FDA1", FDA1)


def test_smoke_ldmt_fda4():
    _run_smoke(LDMTNSGA2(), "FDA4", FDA4)
