"""Non-parametric statistics: Wilcoxon rank-sum + Holm correction, summary."""
import numpy as np

try:
    import pandas as pd
    _HAVE_PANDAS = True
except Exception:  # pragma: no cover
    _HAVE_PANDAS = False

from scipy import stats


def summarize(df, metrics):
    """Return per-group mean/std/median/IQR for each metric."""
    if not _HAVE_PANDAS:  # pragma: no cover
        return _summarize_plain(df, metrics)
    rows = []
    group_cols = [c for c in ("algorithm", "problem", "tau_t") if c in df.columns]
    for _, g in df.groupby(group_cols):
        row = {c: g[c].iloc[0] for c in group_cols}
        for m in metrics:
            if m in g.columns:
                vals = g[m].dropna().astype(float).values
                if vals.size:
                    row[f"{m}_mean"] = float(np.mean(vals))
                    row[f"{m}_std"] = float(np.std(vals, ddof=1)) if vals.size > 1 else 0.0
                    row[f"{m}_median"] = float(np.median(vals))
                    if vals.size > 1:
                        q1, q3 = np.percentile(vals, [25, 75])
                        row[f"{m}_iqr"] = float(q3 - q1)
                    else:
                        row[f"{m}_iqr"] = 0.0
        rows.append(row)
    return pd.DataFrame(rows)


def _summarize_plain(df, metrics):  # pragma: no cover
    return [{"_note": "pandas not available", "n": len(df)}]


def holm_adjust(p_values):
    """Holm step-down adjustment (v2.2, P1-3.2).

    Sorts p-values ascending, multiplies by (m - rank), applies a running
    maximum, truncates at 1.0, and maps back to the ORIGINAL order.

    Parameters
    ----------
    p_values : sequence of floats (finite, in [0, 1])

    Returns
    -------
    list of adjusted p-values in the original order.
    """
    p = np.asarray(list(p_values), dtype=float)
    m = len(p)
    if m == 0:
        return []
    order = np.argsort(p)               # ascending
    adjusted_sorted = np.empty(m, dtype=float)
    running_max = 0.0
    for rank, idx in enumerate(order):
        candidate = (m - rank) * p[idx]
        running_max = max(running_max, candidate)
        adjusted_sorted[rank] = min(1.0, running_max)
    # map back to original order
    adjusted = np.empty(m, dtype=float)
    adjusted[order] = adjusted_sorted
    return adjusted.tolist()


def wilcoxon_holm(df, metric_col, group_col="algorithm", reference=None,
                  paired=True, seed_col="seed"):
    """Pairwise Wilcoxon tests with Holm correction (v2.2, P1-3.2).

    ``paired=True`` (default): seeds are treated as matched pairs — values
    are aligned by ``seed_col`` and scipy.stats.wilcoxon (paired) is used.
    Groups whose seed sets do not fully intersect are compared with the
    independent-samples rank-sum test instead (never silently misaligned).

    Returns a list of dicts: {A, B, p_value, Holm_adjusted_p, significant,
    test, n}.
    """
    groups = {}
    seeds = {}
    if not _HAVE_PANDAS:  # pragma: no cover
        return []
    for name, g in df.groupby(group_col):
        sub = g[[metric_col] + ([seed_col] if seed_col in g.columns else [])]
        sub = sub.dropna(subset=[metric_col])
        groups[name] = sub[metric_col].astype(float).values
        if seed_col in sub.columns:
            seeds[name] = sub[seed_col].astype(int).values
        else:
            seeds[name] = np.arange(len(sub))
    names = list(groups.keys())
    if reference is not None and reference in names:
        pairs = [(reference, b) for b in names if b != reference]
    else:
        pairs = [(names[i], names[j])
                 for i in range(len(names))
                 for j in range(i + 1, len(names))]
    raw = []
    valid_pairs = []
    test_names = []
    ns = []
    for a, b in pairs:
        x, y = groups[a], groups[b]
        if len(x) < 2 or len(y) < 2:
            continue
        # decide paired vs independent by seed alignment
        can_pair = (paired and set(seeds[a]) == set(seeds[b])
                    and len(seeds[a]) == len(set(seeds[a])))
        if can_pair:
            sa = pd.Series(groups[a], index=seeds[a]).sort_index()
            sb = pd.Series(groups[b], index=seeds[b]).sort_index()
            xa, xb = sa.values, sb.values
            try:
                _, p = stats.wilcoxon(xa, xb)
                test_names.append("wilcoxon_paired")
                ns.append(len(xa))
            except Exception:
                continue
        else:
            try:
                _, p = stats.ranksums(x, y)
                test_names.append("ranksums_independent")
                ns.append(min(len(x), len(y)))
            except Exception:
                continue
        raw.append(p)
        valid_pairs.append((a, b))
    adjusted = holm_adjust(raw)
    out = []
    for k, (a, b) in enumerate(valid_pairs):
        out.append({
            "A": a, "B": b, "p_value": float(raw[k]),
            "Holm_adjusted_p": float(adjusted[k]),
            "significant": bool(adjusted[k] < 0.05),
            "test": test_names[k], "n": ns[k],
        })
    return out
