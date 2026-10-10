"""dMOP1 supplementary dynamic multi-objective problem (Farina, Deb, Amato 2004).

M=2, D=10. Uses the convex-form power-2 front (differs from FDA1 which uses sqrt).
Standard formulation: G(t)=sin(0.5*pi*t); f1=x0; g=1+sum((x[1:]-G)^2);
f2 = g*(1-(f1/g)^2). Pareto front (time-invariant): {(u, 1-u^2): u in [0,1]}.
PS: x[:,1:] = G(t).

NOTE: PlatEMO master does not ship a DMOP1.m (verified via the BIMK/PlatEMO
repo). This is the widely-published standard dMOP1; flagged as the
authoritative Farina-Deb-Amato (2004) formulation.
"""
import numpy as np
from .base import DynamicProblem


class DMOP1(DynamicProblem):
    M = 2
    D = 10

    def __init__(self, tau_t=10, n_t=10):
        D = self.D
        xl = np.zeros(D)
        xu = np.ones(D)
        super().__init__(n_var=D, n_obj=2, xl=xl, xu=xu,
                         tau_t=tau_t, n_t=n_t)

    def _calc_objectives(self, X):
        t = self.t
        G = np.sin(0.5 * np.pi * t)
        x0 = np.clip(X[:, 0], 0.0, 1.0)
        tail = X[:, 1:]
        with np.errstate(invalid="ignore"):
            g = 1.0 + np.sum((tail - G) ** 2, axis=1)
            f1 = x0
            f2 = g * (1.0 - (np.maximum(f1, 0.0) / g) ** 2)
        return np.column_stack([f1, f2])
