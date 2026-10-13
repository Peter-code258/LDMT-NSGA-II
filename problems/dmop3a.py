"""Authoritative dMOP3 (Goh & Tan 2009), **Type I** — v2.5 experiment.

Provenance
----------
Goh C.-K., Tan K.-C. IEEE Transactions on Evolutionary Computation, 2009,
13(1): 103-127. DOI 10.1109/TEVC.2008.920671.  Paywalled; this implementation
follows the three independent reproductions listed in
``文档/论文材料/_audit/dmop_definitions.md`` (S1 eq. 18, S2 eqs. 2.25-2.27).

Definition (M = 2, D = 10, all x_i in [0, 1])
---------------------------------------------
    t      = floor(g / tau_t) / n_t
    G(t)   = |sin(0.5 * pi * t)|      (variant A; see the ambiguity note below)
    r      = index of the "position-related" variable, switched over time
    g      = 1 + 9 * sum_{i != r} (x_i - G(t))^2
    f1     = x_r
    f2     = g * (1 - sqrt(f1 / g))

Derived Pareto set / front — **Type I**
---------------------------------------
``g >= 1`` with equality iff ``x_i = G(t)`` for every ``i != r``, so

    PS(t) = { x : x_r in [0,1], x_i = G(t) for all i != r }   ** TIME-VARYING **

while on the PS ``g = 1`` and ``f1 = x_r``, so

    PF = { (u, 1 - sqrt(u)) : u in [0,1] }                    ** TIME-INVARIANT **

So the Pareto SET moves while the FRONT stays fixed — a Type I problem.
The difficulty specific to dMOP3 is the SWITCHING of ``r``: when ``r`` changes,
variables that were position-related become distance-related and must jump to
``G(t)`` simultaneously, which causes severe diversity loss.  That is why this
problem is a useful stress test for warm-start mechanisms.

Ambiguity note (recorded honestly)
----------------------------------
Two reproductions disagree on ``G(t)`` for dMOP3: S1 prints
``G(t) = sin(0.5*pi*t)`` while S2 prints ``G(t) = |sin(0.5*pi*t)|``.  The
original paper could not be read, so no tie-break is possible.  This
implementation uses the ABSOLUTE-VALUE variant (S2), which matches the
``abs(sin(...))`` convention of the FDA suite (FDA3/FDA4/FDA5) and comes with a
complete, self-consistent PS/PF derivation in S2.  The problem TYPE and the
Pareto FRONT are identical under either choice; only the trajectory of the
Pareto SET differs.

Switching policy for ``r``
--------------------------
``r`` is a global "position-related" variable index shared by the whole
population (as in the reference formulations), and it is switched at each
environment change in round-robin fashion.  ``set_generation`` therefore
recomputes ``r`` from the *environment index* ``floor(g / tau_t)``, not from
the raw generation count, so that ``r`` is constant within an environment and
changes exactly when the environment changes.  This keeps the problem
well-defined and reproducible.
"""
import numpy as np
from .base import DynamicProblem


class DMOP3A(DynamicProblem):
    M = 2
    D = 10

    def __init__(self, tau_t=10, n_t=10):
        D = self.D
        xl = np.zeros(D)
        xu = np.ones(D)
        self.r = 0                      # 0-based index of x_r (position var)
        super().__init__(n_var=D, n_obj=2, xl=xl, xu=xu,
                         tau_t=tau_t, n_t=n_t)

    def set_generation(self, g):
        """Advance time AND switch the position-related variable index.

        ``r`` is switched once per environment (round-robin over all D
        variables), mirroring the reference description of dMOP3.
        """
        t = super().set_generation(g)
        self.r = int((g // self.tau_t) % self.D)
        return t

    def _calc_objectives(self, X):
        G = np.abs(np.sin(0.5 * np.pi * self.t))
        r = self.r
        N = X.shape[0]
        f1 = X[:, r]
        # distance variables: all but x_r
        mask = np.ones(X.shape, dtype=bool)
        mask[:, r] = False
        dist = np.where(mask, (X - G) ** 2, 0.0)
        with np.errstate(invalid="ignore"):
            g = 1.0 + 9.0 * np.sum(dist, axis=1)
            f2 = g * (1.0 - np.sqrt(np.maximum(f1, 0.0) / g))
        return np.column_stack([f1, f2])
