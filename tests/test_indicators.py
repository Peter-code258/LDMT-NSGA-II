"""Verify IGD and HV indicators on known sets."""
import numpy as np

from evaluation.igd import igd
from evaluation.hv import hv, get_reference_point
from problems.true_fronts import get_true_pf


def test_igd_self_is_zero():
    pf = get_true_pf("FDA1", n_points=200)
    val = igd(pf, pf)
    assert val == 0.0 or val < 1e-12, f"IGD(PF,PF) should be 0, got {val}"


def test_igd_fda3_current_t():
    pf = get_true_pf("FDA3", n_points=200, t=0.3)
    val = igd(pf, pf)
    assert val < 1e-12


def test_igd_known_distance():
    F_true = np.array([[0.0, 1.0], [1.0, 0.0]])
    F_obt = np.array([[0.0, 0.0]])  # equidistant sqrt(0.5^2+0.5^2)? no: dist from
    # (0,1)->(0,0)=1 ; (1,0)->(0,0)=1 ; mean=1
    val = igd(F_obt, F_true)
    assert np.isclose(val, 1.0, atol=1e-9), f"expected 1.0 got {val}"


def test_hv_known_set():
    # 2D, ref (1.1, 1.1); a single point (0.5, 0.5) -> HV = 1.1*1.1 - 0.6*0.6
    # Actually HV of {(0.5,0.5)} wrt (1.1,1.1) = (1.1-0.5)*(1.1-0.5) = 0.36
    ref = np.array([1.1, 1.1])
    F = np.array([[0.5, 0.5]])
    val = hv(F, ref)
    assert np.isclose(val, 0.36, atol=1e-6), f"expected 0.36 got {val}"


def test_hv_reference_points():
    assert np.allclose(get_reference_point("FDA1"), [1.1, 1.1])
    assert np.allclose(get_reference_point("FDA3"), [2.2, 1.1])
    assert np.allclose(get_reference_point("FDA4"), [1.1, 1.1, 1.1])
