"""Hypervolume (HV) indicator using pymoo's HV implementation.

Fixed reference points per problem (report 5.4.2) — raw objectives are used
(no per-run normalization) so that HV stays comparable across algorithms on
the same problem:

    FDA1  -> (1.1, 1.1)
    FDA3  -> (2.2, 1.1)   # f2 max = (1+G_max)=2 when G=1; f1 max = 1
    FDA4  -> (1.1, 1.1, 1.1)
    dMOP1 -> (1.1, 1.1)
    dMOP2 -> (2.2, 1.1)   # g_pf max = 1+1 = 2
"""
import numpy as np

try:
    from pymoo.indicators.hv import HV
    _HAVE_PYMOO_HV = True
except Exception:  # pragma: no cover
    _HAVE_PYMOO_HV = False


# Fixed reference points per problem (raw objective space).
REF_POINTS = {
    "FDA1": np.array([1.1, 1.1]),
    "FDA3": np.array([2.2, 1.1]),
    "FDA4": np.array([1.1, 1.1, 1.1]),
    "DMOP1": np.array([1.1, 1.1]),
    "DMOP2": np.array([2.2, 1.1]),
}


def get_reference_point(problem_name, M=None):
    name = problem_name.upper()
    ref = REF_POINTS.get(name)
    if ref is None:
        # Fallback: 1.1 per objective.
        ref = np.full(M if M else 2, 1.1)
    return np.array(ref, dtype=float)


def hv(F_obtained, ref_point):
    """Hypervolume of ``F_obtained`` w.r.t. ``ref_point`` (higher is better).

    Points that do not dominate the reference point are discarded. Returns 0.0
    for empty / all-dominated sets.
    """
    F = np.atleast_2d(np.asarray(F_obtained, dtype=float))
    ref = np.asarray(ref_point, dtype=float)
    if F.shape[0] == 0:
        return 0.0
    # keep only points strictly better than the ref in all objectives
    mask = np.all(F < ref[None, :], axis=1)
    F = F[mask]
    if F.shape[0] == 0:
        return 0.0
    if _HAVE_PYMOO_HV:
        try:
            indicator = HV(ref_point=ref)
            return float(indicator(F))
        except Exception:
            pass
    return _hv_exact(F, ref)


def _hv_exact(F, ref):  # pragma: no cover - fallback only
    """Simple inclusion-exclusion HV (works for low M; used as last resort)."""
    n, m = F.shape
    vol = 0.0
    for k in range(1, 2 ** n + 1):
        sign = -1 if (bin(k).count("1") % 2) else 1
        sub = np.ones(m)
        for i in range(n):
            if (k >> i) & 1:
                sub = np.minimum(sub, F[i])
        hv = np.prod(np.maximum(ref - sub, 0.0))
        vol += sign * hv
    return float(vol)
