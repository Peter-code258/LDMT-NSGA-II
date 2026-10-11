"""Verify problem dimensions and bounds (FDA1/FDA3/FDA4)."""
import numpy as np

from problems import FDA1, FDA3, FDA4, DMOP1, DMOP2


def test_fda1_dimensions():
    p = FDA1()
    assert p.n_var == 10
    assert p.n_obj == 2
    assert np.allclose(p.xl, [0.0] + [-1.0] * 9)
    assert np.allclose(p.xu, [1.0] * 10)


def test_fda3_dimensions():
    p = FDA3()
    assert p.n_var == 10
    assert p.n_obj == 2
    assert np.allclose(p.xl, [0.0] + [-1.0] * 9)
    assert np.allclose(p.xu, [1.0] * 10)


def test_fda4_dimensions():
    p = FDA4()
    assert p.n_var == 12  # M + 9 = 3 + 9
    assert p.n_obj == 3
    assert np.allclose(p.xl, [0.0] * 12)
    assert np.allclose(p.xu, [1.0] * 12)


def test_dmop_dimensions():
    assert DMOP1().n_var == 10 and DMOP1().n_obj == 2
    assert DMOP2().n_var == 10 and DMOP2().n_obj == 2
