"""Reference Pareto fronts for FDA / dMOP problems.

FDA1 / dMOP1 / FDA4 have time-invariant fronts; FDA3 / dMOP2 have fronts that
shift with the current time ``t``. Results are cached per (problem, M, n_points).
"""
import numpy as np

_CACHE = {}


def _grid_unit_simplex(M, P):
    """Deterministic unit-simplex grid of P points (0-based inclusive)."""
    if M == 2:
        u = np.linspace(0.0, 1.0, P)
        return np.column_stack([u, 1.0 - u])
    # M >= 3: sample on the simplex then normalize to the first-octant sphere
    rng = np.random.default_rng(20240520)
    X = rng.dirichlet(np.ones(M), size=P)
    return X / np.sqrt(np.sum(X ** 2, axis=1, keepdims=True))


def _pf_fda1(P):
    u = np.linspace(0.0, 1.0, P)
    return np.column_stack([u, 1.0 - np.sqrt(np.maximum(u, 0.0))])


def _pf_fda3(t, P):
    G = np.abs(np.sin(0.5 * np.pi * t))
    g_pf = 1.0 + G
    u = np.linspace(0.0, 1.0, P)
    f2 = g_pf * (1.0 - np.sqrt(np.maximum(u, 0.0) / g_pf))
    return np.column_stack([u, f2])


def _pf_fda4(M, P):
    return _grid_unit_simplex(M, P)


def _pf_dmop1(P):
    u = np.linspace(0.0, 1.0, P)
    return np.column_stack([u, 1.0 - u ** 2])


def _pf_dmop2(t, P):
    G = np.sin(0.5 * np.pi * t)
    g_pf = 1.0 + G
    u = np.linspace(0.0, 1.0, P)
    f2 = g_pf * (1.0 - (np.maximum(u, 0.0) / g_pf) ** 2)
    return np.column_stack([u, f2])


def _pf_dmop3(P):
    u = np.linspace(0.0, 1.0, P)
    f2 = 1.0 - np.sqrt(u) - u * np.sin(10.0 * np.pi * u)
    return np.column_stack([u, f2])


def get_true_pf(problem_name, M=None, n_points=500, t=None):
    """Return the (P, M) reference Pareto front for a problem.

    ``t`` is required for FDA3 / dMOP2 (time-varying fronts); ignored for the
    time-invariant fronts. Cached by (name, M, n_points, t_bucket).
    """
    name = problem_name.upper()
    if M is None:
        M = 2 if name in ("FDA1", "FDA3", "DMOP1", "DMOP2", "DMOP3") else 3
    t_bucket = None if t is None else round(float(t), 6)
    key = (name, M, int(n_points), t_bucket)
    if key in _CACHE:
        return _CACHE[key]

    if name == "FDA1":
        pf = _pf_fda1(n_points)
    elif name == "FDA3":
        tt = 0.0 if t is None else float(t)
        pf = _pf_fda3(tt, n_points)
    elif name == "FDA4":
        pf = _pf_fda4(M, n_points)
    elif name == "DMOP1":
        pf = _pf_dmop1(n_points)
    elif name == "DMOP2":
        tt = 0.0 if t is None else float(t)
        pf = _pf_dmop2(tt, n_points)
    elif name == "DMOP3":
        pf = _pf_dmop3(n_points)
    else:
        raise ValueError(f"Unknown problem name: {problem_name}")

    _CACHE[key] = np.asarray(pf, dtype=float)
    return _CACHE[key]
