"""dMOP1 supplementary dynamic multi-objective problem (PROJECT VARIANT).

v2.4 (attribution fix): the dMOP benchmark suite is due to

    Goh C.-K., Tan K. C. (2009). "A Competitive-Cooperative Coevolutionary
    Paradigm for Dynamic Multiobjective Optimization." IEEE Transactions on
    Evolutionary Computation, 13(1): 103-127.

NOT Farina/Deb/Amato (2004) (the FDA suite) as an earlier version of this
docstring wrongly claimed. The authoritative dMOP1 (Goh & Tan 2009; also
Helbig 2012 PhD thesis, p.52-54, eq. 3.38,
https://www.tech.dmu.ac.uk/~syang/PhDTheses/MardeHelbig-PhDThesis2012.pdf)
uses g = 1 + 9*sum((x[1:] - G(t))^2) and a sqrt-shaped front
f2 = g*(1 - sqrt(f1/g)).

THIS implementation is a project variant and is NOT the authoritative dMOP1:
  - g coefficient is 1 (not 9): g = 1 + sum((x[1:]-G)^2)
  - front exponent is 2 (not sqrt): f2 = g*(1-(f1/g)^2), giving the
    time-invariant PF {(u, 1-u^2): u in [0,1]}
  - PS: x[:,1:] = G(t), G(t)=sin(0.5*pi*t)  (matches the original)
PlatEMO master (commit d25e65d) ships no DMOP1.m (verified via the GitHub
git-tree API), so results on this problem are labelled "supplementary" and
must not be cited as the standard dMOP1.
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
