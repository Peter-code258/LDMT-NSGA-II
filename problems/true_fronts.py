"""Reference Pareto fronts for FDA / dMOP problems.

FDA1 / dMOP1 / FDA4 have time-invariant fronts; FDA3 / dMOP2 have fronts that
shift with the current time ``t``. Results are cached per (problem, M, n_points).
"""
import numpy as np

_CACHE = {}


def _grid_unit_simplex(M, P):
    """Deterministic coverage of the first-octant unit sphere (v2.2, P1-1C).

    Replaces the previous Dirichlet random sampling: points are now generated
    by a fixed angular grid (theta, phi over [0, pi/2]) mapped to the sphere,
    deduplicated and sorted. Deterministic across calls and machines; the
    method and point count are recorded in the front metadata. Note this is a
    project-defined deterministic sampling (equal spacing in angle, not
    equal-area); PlatEMO's GetOptimum uses its own UniformPoint rules.
    """
    if M == 2:
        u = np.linspace(0.0, 1.0, P)
        return np.column_stack([u, 1.0 - u])
    if M == 3:
        n_side = max(2, int(np.ceil(np.sqrt(P))))
        a = np.linspace(0.0, np.pi / 2.0, n_side)
        b = np.linspace(0.0, np.pi / 2.0, n_side)
        A, B = np.meshgrid(a, b)
        f1 = (np.cos(A) * np.cos(B)).ravel()
        f2 = (np.cos(A) * np.sin(B)).ravel()
        f3 = np.sin(A).ravel()
        pts = np.column_stack([f1, f2, f3])
        pts = np.unique(np.round(pts, 12), axis=0)
        order = np.lexsort((pts[:, 2], pts[:, 1], pts[:, 0]))
        return pts[order]
    # M >= 4: deterministic scalar grid is impractical; fall back to a FIXED
    # seeded Dirichlet draw (documented as project-defined, not PlatEMO).
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
