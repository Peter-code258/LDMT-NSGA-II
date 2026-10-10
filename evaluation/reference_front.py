"""Dispatch helper that returns the reference front for a problem."""
import numpy as np
from problems.true_fronts import get_true_pf


def reference_front(problem_name, M=None, n_points=500, t=None):
    """Return the true Pareto front (P, M) for ``problem_name``.

    Thin wrapper around :func:`problems.true_fronts.get_true_pf`.
    """
    return np.asarray(get_true_pf(problem_name, M=M, n_points=n_points, t=t),
                       dtype=float)
