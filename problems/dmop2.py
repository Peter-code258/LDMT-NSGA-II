"""dMOP2 supplementary dynamic multi-objective problem (PROJECT VARIANT).

M=2, D=10. PROJECT IMPLEMENTATION (not authoritative): PlatEMO master
(commit d25e65d) ships no DMOP2.m (verified via the GitHub git-tree API).

v2.4 (K3 docstring fix): the authoritative dMOP2 of Goh & Tan (2009), IEEE
Trans. Evol. Comput. 13(1):103-127 — also Helbig 2012 PhD thesis, p.52-54,
eq. 3.39 (https://www.tech.dmu.ac.uk/~syang/PhDTheses/MardeHelbig-PhDThesis2012.pdf)
— is a Type-II problem with
    g = 1 + 9*sum((x[1:] - G(t))^2),  f2 = g*(1 - (f1/g)^{H(t)}),
    H(t) = 0.75*sin(0.5*pi*t) + 1.25,  G(t) = sin(0.5*pi*t),
i.e. BOTH the Pareto set and the Pareto front vary with time through H(t).

THIS implementation was verified to DIFFER from the authoritative dMOP2:
  - g coefficient is 1 (not 9) and g carries an extra +G term:
    g = 1 + G + sum((x[1:]-G)^2)  (mirroring FDA3's strategy)
  - front exponent is fixed at 2 (no time-varying H(t)):
    f2 = g*(1-(f1/g)^2)
  - at the PS (x[1:]=G): g_pf = 1+G, so the PF
    {(u, (1+G)*(1-(u/(1+G))^2)) : u in [0,1]} shifts with t.
G(t) = sin(0.5*pi*t). PS: x[:,1:] = G.
Results on this problem are labelled "supplementary" and must not be cited
as the standard dMOP2.
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
