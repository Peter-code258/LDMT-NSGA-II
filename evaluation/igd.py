"""Inverted Generational Distance (IGD) indicator."""
import numpy as np
from scipy.spatial.distance import cdist


def igd(F_obtained, F_true):
    """Mean distance from each true point to the nearest obtained point.

    Both arrays are (P, M). Returns a float (lower is better; 0 = perfect).
    """
    F_obtained = np.atleast_2d(np.asarray(F_obtained, dtype=float))
    F_true = np.atleast_2d(np.asarray(F_true, dtype=float))
    if F_obtained.shape[0] == 0:
        return float("inf")
    d = cdist(F_true, F_obtained)  # (P_true, P_obtained)
    return float(np.mean(d.min(axis=1)))
