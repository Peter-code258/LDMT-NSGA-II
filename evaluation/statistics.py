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


def wilcoxon_holm(df, metric_col, group_col="algorithm", reference=None):
    """Pairwise Wilcoxon rank-sum tests with Holm correction.

    Returns a list of dicts: {A, B, p_value, Holm_adjusted_p, significant}.
    ``reference`` (optional) limits comparisons to (algo vs reference).
    """
    groups = {}
    if not _HAVE_PANDAS:  # pragma: no cover
        return []
    for name, g in df.groupby(group_col):
        groups[name] = g[metric_col].dropna().astype(float).values
    names = list(groups.keys())
    if reference is not None and reference in names:
        pairs = [(reference, b) for b in names if b != reference]
    else:
        pairs = [(names[i], names[j])
                 for i in range(len(names))
                 for j in range(i + 1, len(names))]
    raw = []
    valid_pairs = []
    for a, b in pairs:
        x, y = groups[a], groups[b]
        if len(x) < 2 or len(y) < 2:
            continue
        try:
            _, p = stats.ranksums(x, y)
        except Exception:
            continue
        raw.append(p)
        valid_pairs.append((a, b))
    m = len(raw)
    # Holm correction
    order = np.argsort(raw)
    adj = [0.0] * m
    running = 0.0
    for rank, idx in enumerate(order):
        holm_p = min(1.0, raw[idx] * (m - rank))
        running = max(running, holm_p)
        adj[idx] = running
    out = []
    for k, (a, b) in enumerate(valid_pairs):
        out.append({
            "A": a, "B": b, "p_value": float(raw[k]),
            "Holm_adjusted_p": float(adj[k]),
            "significant": bool(adj[k] < 0.05),
        })
    return out
