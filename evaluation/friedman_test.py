"""Friedman + Nemenyi post-hoc statistical tests for multi-algorithm comparison.

Standard procedure in evolutionary multi-objective optimisation:
    1.  For each problem, rank algorithms per independent run (seed).
    2.  Friedman test checks whether the ranking distributions differ.
    3.  If significant, Nemenyi post-hoc identifies pairwise differences.

Returns a results dict that can be serialised to CSV/JSON.
"""
import numpy as np
from scipy import stats


def _rank_per_seed(metric_matrix, higher_is_better=False):
    """Rank algorithms per seed (column-wise). Lower rank = better.

    Parameters
    ----------
    metric_matrix : (n_seeds, n_algorithms) array
    higher_is_better : bool

    Returns
    -------
    ranks : (n_seeds, n_algorithms) array, 1 = best
    """
    if higher_is_better:
        # negate so that ascending rank = best
        rankable = -metric_matrix
    else:
        rankable = metric_matrix
    # v2.4 (K7 fix): average ranks for ties (scipy rankdata, method='average')
    # instead of argsort-argsort which assigned arbitrary distinct ranks.
    ranks = np.apply_along_axis(stats.rankdata, 1, rankable)
    return ranks.astype(float)


def friedman_test(metric_matrix, higher_is_better=False):
    """Friedman chi-square test on ranked data.

    Returns
    -------
    dict with keys:
        chi2, p_value, df, ranks (avg), rank_matrix
    """
    data = np.asarray(metric_matrix, dtype=float)
    if data.ndim != 2:
        raise ValueError("metric_matrix must be 2-D (seeds x algorithms)")
    n_seeds, n_algorithms = data.shape
    rank_matrix = _rank_per_seed(data, higher_is_better)
    avg_ranks = rank_matrix.mean(axis=0)
    # Friedman statistic
    sum_ranks_sq = np.sum(avg_ranks ** 2)
    chi2 = (12.0 * n_seeds / (n_algorithms * (n_algorithms + 1.0))) * (
        sum_ranks_sq - n_algorithms * (n_algorithms + 1.0) ** 2 / 4.0)
    df = n_algorithms - 1
    p_value = 1.0 - stats.chi2.cdf(chi2, df)
    return {
        "chi2": float(chi2),
        "p_value": float(p_value),
        "df": int(df),
        "avg_ranks": avg_ranks.tolist(),
        "rank_matrix": rank_matrix.tolist(),
    }


def nemenyi_posthoc(rank_matrix, alpha=0.05):
    """Nemenyi post-hoc test after Friedman.

    Returns
    -------
    dict with keys:
        critical_difference, significant_pairs (list of (i,j) 0-based)
    """
    rank_matrix = np.asarray(rank_matrix, dtype=float)
    n_seeds, n_algorithms = rank_matrix.shape
    avg_ranks = rank_matrix.mean(axis=0)
    # Nemenyi critical values q_alpha (two-tailed) for 2..10 algorithms
    # (Demšar 2006, Table 5). v2.4 (K7 fix): the alpha parameter now actually
    # selects the table; unsupported alpha raises instead of silently using
    # the 0.05 table.
    q_tables = {
        0.10: {2: 1.645, 3: 2.052, 4: 2.291, 5: 2.459, 6: 2.589,
               7: 2.693, 8: 2.780, 9: 2.855, 10: 2.920},
        0.05: {2: 1.960, 3: 2.343, 4: 2.569, 5: 2.728, 6: 2.850,
               7: 2.949, 8: 3.031, 9: 3.102, 10: 3.164},
    }
    if alpha not in q_tables:
        raise ValueError(
            f"Nemenyi q-table only available for alpha in "
            f"{sorted(q_tables)}; got alpha={alpha}")
    q_alpha_table = q_tables[alpha]
    if n_algorithms not in q_alpha_table:
        # fallback: use q_alpha for 10 (conservative)
        q_alpha = q_alpha_table[10]
    else:
        q_alpha = q_alpha_table[n_algorithms]
    critical_difference = q_alpha * np.sqrt(
        n_algorithms * (n_algorithms + 1.0) / (6.0 * n_seeds))
    significant_pairs = []
    for i in range(n_algorithms):
        for j in range(i + 1, n_algorithms):
            if abs(avg_ranks[i] - avg_ranks[j]) > critical_difference:
                significant_pairs.append((i, j))
    return {
        "critical_difference": float(critical_difference),
        "significant_pairs": significant_pairs,
        "alpha": float(alpha),
    }


def friedman_nemenyi(metric_matrix, higher_is_better=False, alpha=0.05):
    """Convenience wrapper: Friedman + Nemenyi in one call."""
    fried = friedman_test(metric_matrix, higher_is_better)
    nem = nemenyi_posthoc(fried["rank_matrix"], alpha=alpha)
    return {"friedman": fried, "nemenyi": nem}
