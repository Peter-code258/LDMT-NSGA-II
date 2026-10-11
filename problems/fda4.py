"""FDA4 dynamic multi-objective problem (Farina & Deb, 2003).

M=3, D=M+9=12. Time-INVARIANT Pareto front: first octant of the unit
hypersphere. PS: x[:, M:] = G where G = abs(sin(0.5*pi*t)). The variable
at index M-1 (0-based) is unused, matching PlatEMO FDA4.m: position vars
are x[:, :M-1], distance vars are x[:, M:].
"""
import numpy as np
from .base import DynamicProblem


class FDA4(DynamicProblem):
    M = 3
    D = 12  # M + 9

    def __init__(self, tau_t=10, n_t=10, M=None):
        if M is not None:
            self.M = int(M)
            self.D = self.M + 9
        D, M = self.D, self.M
        xl = np.zeros(D)
        xu = np.ones(D)
        super().__init__(n_var=D, n_obj=M, xl=xl, xu=xu,
                         tau_t=tau_t, n_t=n_t)

    def _calc_objectives(self, X):
        t = self.t
        M = self.M
        N = X.shape[0]
        G = np.abs(np.sin(0.5 * np.pi * t))
        with np.errstate(invalid="ignore"):
            g = np.sum((X[:, M:] - G) ** 2, axis=1)
            theta = X[:, :M - 1] * (np.pi / 2.0)
            cos_th = np.cos(theta)                      # (N, M-1)
            base = np.ones((N, 1))
            cumprod = np.cumprod(np.hstack([base, cos_th]), axis=1)  # (N, M)
            cumprod_flip = np.fliplr(cumprod)            # (N, M)
            sin_th_rev = np.sin(theta[:, ::-1]) if M > 1 else np.zeros((N, 0))
            sin_part = np.hstack([np.ones((N, 1)), sin_th_rev])  # (N, M)
            g_col = g.reshape(-1, 1)
            F = (1.0 + g_col) * cumprod_flip * sin_part
        return F
