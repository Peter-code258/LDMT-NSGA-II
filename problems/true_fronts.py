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


# ---------------------------------------------------------------------------
# v2.5: authoritative dMOP1-3 (Goh & Tan 2009) reference fronts.
#
# Types and formulas were re-verified online against three independent
# peer-reviewed reproductions (see 文档/论文材料/_audit/dmop_definitions.md).
# An EARLIER draft of this project had dMOP1 mis-classified as Type I with a
# sqrt front; that was wrong and is corrected here.
#
#   dMOP1 (auth.): **Type III** -- POS static, POF time-varying.
#                  g = 1 + 9*sum_{i>=2} x_i^2  (NO G(t)),  f1 = x_1,
#                  f2 = g*(1 - (f1/g)^H(t)).
#                  PS = {(x_1, 0, ..., 0)}  (TIME-INVARIANT)
#                  PF(t) = {(u, 1 - u^H(t))}  (TIME-VARYING)
#
#   dMOP2 (auth.): **Type II** -- both POS and POF move.
#                  g = 1 + 9*sum_{i>=2} (x_i - G(t))^2,  G(t) = sin(0.5*pi*t),
#                  f2 = g*(1 - (f1/g)^H(t)).
#                  PS(t) = {(x_1, G(t), ..., G(t))}   PF(t) = {(u, 1 - u^H(t))}
#
#   dMOP3 (auth.): **Type I** -- POS moves, POF static.
#                  g = 1 + 9*sum_{i != r} (x_i - G(t))^2,  f1 = x_r,
#                  f2 = g*(1 - sqrt(f1/g)),  G(t) = |sin(0.5*pi*t)| (variant A).
#                  PS(t) = {x_r in [0,1], x_i = G(t) for i != r}
#                  PF = {(u, 1 - sqrt(u))}  (TIME-INVARIANT)
#
# H(t) = 0.75*sin(0.5*pi*t) + 1.25 for dMOP1 and dMOP2.
# ---------------------------------------------------------------------------


def _H_authoritative(t):
    return 0.75 * np.sin(0.5 * np.pi * t) + 1.25


def _pf_dmop1_authoritative(t, P):
    """Type III front of the authoritative dMOP1: {(u, 1 - u^H(t))}."""
    H = _H_authoritative(t)
    u = np.clip(np.linspace(0.0, 1.0, P), 0.0, 1.0)
    return np.column_stack([u, 1.0 - u ** H])


def _pf_dmop2_authoritative(t, P):
    """Type II front of the authoritative dMOP2: {(u, 1 - u^H(t))}."""
    H = _H_authoritative(t)
    u = np.clip(np.linspace(0.0, 1.0, P), 0.0, 1.0)
    return np.column_stack([u, 1.0 - u ** H])


def _pf_dmop3_authoritative(P):
    """Type I front of the authoritative dMOP3: {(u, 1 - sqrt(u))}, static."""
    u = np.linspace(0.0, 1.0, P)
    return np.column_stack([u, 1.0 - np.sqrt(np.maximum(u, 0.0))])



def get_true_pf(problem_name, M=None, n_points=500, t=None):
    """Return the (P, M) reference Pareto front for a problem.

    ``t`` is required for FDA3 / dMOP2 (time-varying fronts); ignored for the
    time-invariant fronts. Cached by (name, M, n_points, t_bucket).
    """
    name = problem_name.upper()
    # v3.2 (plan 4.2): the extended stress variants are built on FDA1 and share
    # its TIME-INVARIANT front, so they resolve to the same analytical front.
    if name in ("IRREG_FDA1", "REV_FDA1", "GRAD_FDA1", "NOISY_FDA1"):
        name = "FDA1"
    if M is None:
        M = 3 if name in ("FDA4", "FDA5") else 2
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
    elif name == "DMOP1A":
        tt = 0.0 if t is None else float(t)
        pf = _pf_dmop1_authoritative(tt, n_points)
    elif name == "DMOP2A":
        tt = 0.0 if t is None else float(t)
        pf = _pf_dmop2_authoritative(tt, n_points)
    elif name == "DMOP3A":
        pf = _pf_dmop3_authoritative(n_points)
    else:
        raise ValueError(f"Unknown problem name: {problem_name}")

    _CACHE[key] = np.asarray(pf, dtype=float)
    return _CACHE[key]
