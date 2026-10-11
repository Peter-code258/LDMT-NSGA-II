"""FDA3 dynamic multi-objective problem (Farina & Deb, 2003).

M=2, D=10. Time-VARYING Pareto front depending on ``G = abs(sin(0.5*pi*t))``.
PS: x[:,1:] = G. Formulas match PlatEMO master FDA3.m.
"""
import numpy as np
from .base import DynamicProblem


class FDA3(DynamicProblem):
    M = 2
    D = 10

    def __init__(self, tau_t=10, n_t=10):
        D = self.D
        xl = np.array([0.0] + [-1.0] * (D - 1))
        xu = np.array([1.0] * D)
        super().__init__(n_var=D, n_obj=2, xl=xl, xu=xu,
                         tau_t=tau_t, n_t=n_t)

    def _calc_objectives(self, X):
        t = self.t
        G = np.abs(np.sin(0.5 * np.pi * t))
        # x0 in [0,1]; exponent 10**(2*sin(0.5*pi*t)) keeps f1 in [0,1]
        exp_k = 10.0 ** (2.0 * np.sin(0.5 * np.pi * t))
        x0 = np.clip(X[:, 0], 0.0, 1.0)
        tail = X[:, 1:]
        with np.errstate(invalid="ignore"):
            f1 = x0 ** exp_k
            g = 1.0 + G + np.sum((tail - G) ** 2, axis=1)
            f2 = g * (1.0 - np.sqrt(np.maximum(f1, 0.0) / g))
        return np.column_stack([f1, f2])
