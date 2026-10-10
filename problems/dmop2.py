"""dMOP2 supplementary dynamic multi-objective problem.

M=2, D=10. PROJECT IMPLEMENTATION (not authoritative): PlatEMO master does
not ship a DMOP2.m (verified via the BIMK/PlatEMO repo). To provide a
distinct, time-VARYING counterpart to dMOP1 we use a g that depends on G at
the optimum (mirroring FDA3's strategy): g = 1 + G + sum((x[1:]-G)^2);
f1 = x0; f2 = g*(1-(f1/g)^2). At the PS (x[1:]=G): g_pf=1+G, so the PF is
{(u, (1+G)*(1-(u/(1+G))^2)) : u in [0,1]} which shifts with t.
G(t) = sin(0.5*pi*t). PS: x[:,1:] = G.
"""
import numpy as np
from .base import DynamicProblem


class DMOP2(DynamicProblem):
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
            g = 1.0 + G + np.sum((tail - G) ** 2, axis=1)
            f1 = x0
            f2 = g * (1.0 - (np.maximum(f1, 0.0) / g) ** 2)
        return np.column_stack([f1, f2])
