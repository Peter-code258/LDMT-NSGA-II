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
    for empty / all-dominated sets. If pymoo's HV is unavailable, the exact
    fallback is used ONLY for tiny inputs (<= 20 points, <= 3 objectives);
    larger inputs raise RuntimeError instead of silently attempting an
    exponential-time computation (v2.2, P1-3.4).
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
    n, m = F.shape
    if n <= 20 and m <= 3:
        return _hv_exact(F, ref)
    raise RuntimeError(
        f"pymoo HV unavailable and exact fallback refused for n={n}, m={m} "
        "(exponential cost); install pymoo to compute this HV.")


_HV_FALLBACK_MAX_POINTS = 20
_HV_FALLBACK_MAX_OBJ = 3


def _hv_exact(F, ref):  # pragma: no cover - fallback only
    """Exact HV by inclusion-exclusion over boxes [F_i, ref] (low M only).

    v2.4 (K5 fix): three bugs corrected —
    1. sign was inverted: odd-size subsets must contribute +vol;
    2. box intersection uses the component-wise MAX of corner points
       (not min): ∩[F_i, ref] = [max F_i, ref];
    3. loop ran to 2**n inclusive, adding a spurious empty-subset term.
    Verified against pymoo's HV on random small sets (test_v24_fixes).
    """
    n, m = F.shape
    vol = 0.0
    for k in range(1, 2 ** n):
        corner = np.ones(m) * -np.inf
        size = 0
        for i in range(n):
            if (k >> i) & 1:
                corner = np.maximum(corner, F[i])
                size += 1
        v = float(np.prod(np.maximum(ref - corner, 0.0)))
        vol += v if (size % 2 == 1) else -v
    return float(vol)
