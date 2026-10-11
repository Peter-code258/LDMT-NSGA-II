"""FDA1 dynamic multi-objective problem (Farina & Deb, 2003).

M=2, D=10. Time-invariant Pareto front: {(u, 1 - sqrt(u)) : u in [0,1]}.
PS: x[:,1:] = sin(0.5*pi*t). Formulas match PlatEMO master FDA1.m.
"""
import numpy as np
from .base import DynamicProblem


class FDA1(DynamicProblem):
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
        x0 = X[:, 0]
        tail = X[:, 1:]
        with np.errstate(invalid="ignore"):
            g = 1.0 + np.sum((tail - np.sin(0.5 * np.pi * t)) ** 2, axis=1)
            f1 = x0
            f2 = g * (1.0 - np.sqrt(np.maximum(f1, 0.0) / g))
        return np.column_stack([f1, f2])
