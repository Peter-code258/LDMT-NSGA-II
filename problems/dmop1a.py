"""Authoritative dMOP1 (Goh & Tan 2009), **Type III** — v2.5 experiment.

Provenance
----------
Goh C.-K., Tan K.-C. "A Competitive-Cooperative Coevolutionary Paradigm for
Dynamic Multiobjective Optimization." IEEE Transactions on Evolutionary
Computation, 2009, 13(1): 103-127. DOI 10.1109/TEVC.2008.920671.

The original paper is paywalled with no open-access copy, so this
implementation follows three independent peer-reviewed reproductions that
quote its equations verbatim (all recorded in
``文档/论文材料/_audit/dmop_definitions.md``):

  S1 Helbig & Engelbrecht, "Dynamic multi-objective optimisation using PSO",
     CSIR ResearchSpace, eq. (16)
  S2 Joćko, MSc thesis, Brock University 2021, eqs. (2.21)-(2.22)
  S4 MDPI Mathematics 2022, 10(15):1587, Tables 2-3

Definition (M = 2, D = 10, **all** x_i in [0, 1])
------------------------------------------------
    t      = floor(g / tau_t) / n_t
    H(t)   = 0.75 * sin(0.5 * pi * t) + 1.25          in [0.5, 2.0]
    g      = 1 + 9 * sum_{i=2..D} x_i^2               (NOTE: no G(t)!)
    f1     = x_1
    f2     = g * (1 - (f1 / g) ** H(t))

Derived Pareto set / front — **Type III**, NOT Type I
-----------------------------------------------------
``g >= 1`` with equality iff ``x_i = 0`` for all ``i >= 2``, so

    PS = { (x_1, 0, ..., 0) : x_1 in [0,1] }        ** TIME-INVARIANT **

while on the PS ``g = 1``, ``f1 = x_1`` and therefore

    PF(t) = { (u, 1 - u ** H(t)) : u in [0,1] }    ** TIME-VARYING **

i.e. the Pareto-optimal SET is frozen while the Pareto FRONT rotates between
convex (H > 1), linear (H = 1) and concave (H < 1) as t advances.  That is the
definition of a Type III problem (POS static, POF changing).

This is the *opposite* classification from an earlier draft of this project,
which had labelled dMOP1 as Type I with a sqrt front.  The two must not be
confused: dMOP1 is the one dMOP case in which the algorithm must re-shape an
already-correct solution set rather than move it.

Note dMOP1 and dMOP2 share the same PF *family* ``1 - u^H(t)``; the difference
is that dMOP1's PS is static whereas dMOP2's PS slides to ``x_i = G(t)``.
"""
import numpy as np
from .base import DynamicProblem


class DMOP1A(DynamicProblem):
    M = 2
    D = 10

    def __init__(self, tau_t=10, n_t=10):
        D = self.D
        xl = np.zeros(D)
        xu = np.ones(D)
        super().__init__(n_var=D, n_obj=2, xl=xl, xu=xu,
                         tau_t=tau_t, n_t=n_t)

    def _calc_objectives(self, X):
        H = 0.75 * np.sin(0.5 * np.pi * self.t) + 1.25
        x1 = X[:, 0]
        tail = X[:, 1:]
        with np.errstate(invalid="ignore"):
            g = 1.0 + 9.0 * np.sum(tail ** 2, axis=1)
            f1 = x1
            f2 = g * (1.0 - (np.maximum(f1, 0.0) / g) ** H)
        return np.column_stack([f1, f2])
