"""DMOP3 dynamic multi-objective problem (PROJECT CUSTOM, static).

v2.4 (K2 docstring fix + dead-code removal): an earlier docstring attributed
this problem to "Farina, Deb, Amato 2004" and claimed the POF shifts with
time. Both statements were WRONG:

  - The dMOP suite is due to Goh & Tan (2009), IEEE Trans. Evol. Comput.
    13(1):103-127 — NOT Farina/Deb/Amato (2004) (the FDA suite).
  - The authoritative dMOP3 (Goh & Tan 2009; Helbig 2012 PhD thesis, p.52-54,
    eq. 3.40, https://www.tech.dmu.ac.uk/~syang/PhDTheses/MardeHelbig-PhDThesis2012.pdf)
    is a Type-I problem: the PARETO SET is time-varying while the PARETO
    FRONT stays static.

THIS implementation is a project-custom STATIC deceptive problem: neither
the PS nor the PF changes with time (the old code computed
G = |sin(0.5*pi*t)| and never used it — dead code, removed in v2.4).
    g = 1 + sum((x[1:])^2);  f1 = x0;
    f2 = g * (1 - sqrt(f1/g) - (f1/g)*sin(10*pi*f1)).
Pareto front (time-invariant, multi-modal/deceptive):
    {(u, 1 - sqrt(u) - u*sin(10*pi*u)) : u in [0,1]}.
PS: x[:,1:] = 0  (static).

It is kept as a static-control problem ("does the transfer mechanism hurt
when nothing moves?") and must not be cited as the standard dMOP3.
"""
import numpy as np
from .base import DynamicProblem


class DMOP3(DynamicProblem):
    M = 2
    D = 10

    def __init__(self, tau_t=10, n_t=10):
        D = self.D
        xl = np.zeros(D)
        xu = np.ones(D)
        super().__init__(n_var=D, n_obj=2, xl=xl, xu=xu,
                         tau_t=tau_t, n_t=n_t)

    def _calc_objectives(self, X):
        # v2.4 (K2): dead code removed — G = |sin(0.5*pi*t)| was computed but
        # never used (the problem is static by design; see module docstring).
        x0 = np.clip(X[:, 0], 0.0, 1.0)
        tail = X[:, 1:]
        with np.errstate(invalid="ignore"):
            g = 1.0 + np.sum(tail ** 2, axis=1)
            f1 = x0
            f2 = g * (1.0 - np.sqrt(np.maximum(f1, 0.0) / g)
                      - (f1 / g) * np.sin(10.0 * np.pi * f1))
        return np.column_stack([f1, f2])
