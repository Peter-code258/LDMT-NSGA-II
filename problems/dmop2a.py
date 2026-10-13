"""Authoritative dMOP2 (Goh & Tan 2009), Type II — v2.5 experiment.

See ``dmop1a.py`` for the provenance note and the audit trail
(文档/论文材料/_audit/dmop_definitions.md).

Definition (M = 2, D = 10)
-------------------------
    t      = floor(g / tau_t) / n_t
    G(t)   = sin(0.5 * pi * t)
    H(t)   = 0.75 * sin(0.5 * pi * t) + 1.25          in [0.5, 2.0]
    g      = 1 + 9 * sum_{i=2..D} (x_i - G(t))^2
    f1     = x_1
    f2     = g * (1 - (f1 / g) ** H(t))
    bounds x_1 in [0,1]; x_i in [-1,1] for i >= 2

Derived Pareto set / front
--------------------------
    PS(t) = { (x_1, G(t), ..., G(t)) : x_1 in [0,1] }        (moves with t)
    PF(t) = { (u, 1 - u ** H(t)) : u in [0,1] }               (moves with t)

Both the Pareto set AND the Pareto front vary with time, which makes dMOP2 a
**Type II** problem: at H = 2 the front is the convex parabola 1 - u^2, while
at H -> 0.5 it bends towards 1 - sqrt(u). The time-varying exponent is the
essential difference from the project variant ``problems/dmop2.py``, whose
exponent is frozen at 2 and whose g carries an extra +G term.
"""
import numpy as np
from .base import DynamicProblem


class DMOP2A(DynamicProblem):
    M = 2
    D = 10

    def __init__(self, tau_t=10, n_t=10):
        D = self.D
        xl = np.array([0.0] + [-1.0] * (D - 1))
        xu = np.ones(D)
        super().__init__(n_var=D, n_obj=2, xl=xl, xu=xu,
                         tau_t=tau_t, n_t=n_t)

    def _calc_objectives(self, X):
        t = self.t
        G = np.sin(0.5 * np.pi * t)
        H = 0.75 * np.sin(0.5 * np.pi * t) + 1.25
        x1 = X[:, 0]
        tail = X[:, 1:]
        with np.errstate(invalid="ignore"):
            g = 1.0 + 9.0 * np.sum((tail - G) ** 2, axis=1)
            f1 = x1
            f2 = g * (1.0 - (np.maximum(f1, 0.0) / g) ** H)
        return np.column_stack([f1, f2])
