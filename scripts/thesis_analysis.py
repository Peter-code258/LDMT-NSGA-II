"""Thesis-grade analysis for the v2.4 + v2.5 data sets.

This script NEVER writes into ``results/summary`` (which holds the frozen
v2.4 aggregates produced by ``analyze.py``).  It reads raw per-generation CSVs
from one or more directories and writes its own outputs to
``results/summary_v25``.

METRIC DEFINITIONS (exactly as implemented here, and as documented in the
thesis section 5.1.3)
---------------------
per run = one (algorithm, problem, tau_t, N, n_t, seed) combination:

  steady-state MIGD ("MIGD")   mean over environment segments of the IGD of
                               the LAST generation of that segment, then mean
                               over seeds.  (This is the thesis main metric.)
  steady-state MHV ("MHV")     same aggregation applied to HV.
  post-change MIGD             mean over CHANGED segments (env_index >= 1) of
                               the IGD of generation 0 of that segment, then
                               mean over seeds.
  all-generation MIGD          mean of ALL per-generation IGD values of the
                               run, then mean over seeds.  (Auxiliary only.)
  recovery generations         for each changed environment e: the number of
                               generations after the change until IGD first
                               becomes <= the last-generation IGD of
                               environment e-1; tau_t if it never does.
  migration time               sum of ``migration_runtime_s`` over the run,
                               then mean over seeds (per-environment quantity
                               recorded only in the first generation of a
                               changed environment).

USAGE
-----
    python -m scripts.thesis_analysis
"""
import json
import os
import sys
from itertools import combinations

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from evaluation.statistics import holm_adjust                # noqa: E402
from scipy import stats                                      # noqa: E402

RAW_DIRS = ["results/raw", "results/raw_v25"]
OUT = "results/summary_v25"

# canonical algorithm ids -> display order for the v2.4 comparison
V24_ALGOS = ["nsga2_static", "dnsga2_inspired", "ar1_predictive",
             "gmm_inspired", "pps", "ldt_nsga2"]
V25_ALGOS = ["nsga2_static", "dnsga2_inspired", "ar1_predictive",
             "gmm_inspired", "pps", "ldt_alpha1"]
ALGO_LABEL = {
    "nsga2_static": "nsga2_static",
    "dnsga2_inspired": "dnsga2_inspired",
    "ar1_predictive": "ar1_predictive",
    "gmm_inspired": "gmm_inspired",
    "pps": "pps",
    "ldt_nsga2": "ldt_nsga2 (alpha=0.5)",
    "ldt_alpha1": "ldt_alpha1 (alpha=1.0)",
}

RUN_KEYS = ["algorithm", "problem", "tau_t", "population_size", "n_t", "seed"]
GROUP_KEYS = ["_source", "_kind"] + RUN_KEYS


# --------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------
def _kind_from_filename(f):
    """Experiment group of a raw CSV, from its filename prefix."""
    for k in ("comparison", "supplementary", "ablation", "sensitivity"):
        if f.startswith(k + "_"):
            return k
    return "other"


def load_all(raw_dirs):
    frames = []
    for d in raw_dirs:
        if not os.path.isdir(d):
            continue
        for f in sorted(os.listdir(d)):
            if not f.endswith(".csv"):
                continue
            p = os.path.join(d, f)
            try:
                df = pd.read_csv(p)
            except Exception as exc:            # pragma: no cover
                print(f"  [SKIP] {p}: {exc}")
                continue
            if df.empty or "IGD" not in df.columns:
                continue
            df["_source"] = d
            df["_file"] = f
            df["_kind"] = _kind_from_filename(f)
            frames.append(df)
    if not frames:
        raise SystemExit("no raw CSV data found")
    return pd.concat(frames, ignore_index=True)


# --------------------------------------------------------------------------
# per-run metrics
# --------------------------------------------------------------------------
def per_run_metrics(df):
    """One row per run with all aggregated metrics (see module docstring)."""
    g = df.copy()
    g["gen_in_env"] = g["generation"] % g["tau_t"]

    idx_last = g.groupby(GROUP_KEYS + ["environment_index"])["generation"].idxmax()
    last = g.loc[idx_last]
    idx_first = g.groupby(GROUP_KEYS + ["environment_index"])["generation"].idxmin()
    first = g.loc[idx_first]

    steady = last.groupby(GROUP_KEYS).agg(
        MIGD=("IGD", "mean"), MHV=("HV", "mean"),
        n_envs=("environment_index", "nunique")).reset_index()

    fc = first[first["environment_index"] > 0]
    post = fc.groupby(GROUP_KEYS).agg(
        post_change_MIGD=("IGD", "mean")).reset_index()

    allgen = g.groupby(GROUP_KEYS).agg(
        allgen_MIGD=("IGD", "mean")).reset_index()

    mig = g.groupby(GROUP_KEYS).agg(
        total_mig_time_s=("migration_runtime_s", "sum"),
        mean_gen_time_s=("generation_runtime_s", "mean"),
        n_eval_total=("n_evaluations", "max"),
        code_version=("code_version", "first")).reset_index()

    out = steady.merge(post, on=GROUP_KEYS, how="left") \
                .merge(allgen, on=GROUP_KEYS, how="left") \
                .merge(mig, on=GROUP_KEYS, how="left")
    return out


def recovery_generations(df):
    """Recovery generations per (run, environment), per analyze.py's rule."""
    rows = []
    for keys, grp in df.groupby(RUN_KEYS, sort=False):
        envs = sorted(grp["environment_index"].unique())
        if len(envs) < 2:
            continue
        for e in envs[1:]:
            pre = grp[grp["environment_index"] == e - 1]
            post = grp[grp["environment_index"] == e]
            if pre.empty or post.empty:
                continue
            pre_last = pre.loc[pre["generation"].idxmax(), "IGD"]
            post_sorted = post.sort_values("generation")
            hit = post_sorted[post_sorted["IGD"] <= pre_last]
            rec = (int(hit["generation"].min() - post_sorted["generation"].min())
                   if len(hit) else int(post_sorted["tau_t"].iloc[0]))
            rows.append({"algorithm": keys[0], "problem": keys[1],
                         "tau_t": keys[2], "seed": keys[5],
                         "env": int(e), "recovery_gens": rec})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# statistics
# --------------------------------------------------------------------------
def summarize(runs, algos):
    r = runs[runs["algorithm"].isin(algos)]
    out = r.groupby(["problem", "tau_t", "algorithm"]).agg(
        MIGD_mean=("MIGD", "mean"), MIGD_std=("MIGD", "std"),
        MIGD_med=("MIGD", "median"),
        MHV_mean=("MHV", "mean"), MHV_std=("MHV", "std"),
        post_change_MIGD=("post_change_MIGD", "mean"),
        allgen_MIGD=("allgen_MIGD", "mean"),
        mig_time_mean=("total_mig_time_s", "mean"),
        n_runs=("seed", "nunique")).reset_index()
    return out


def paired_wilcoxon_holm(runs, ref, others):
    """Paired Wilcoxon (ref vs each other) per (problem, tau), Holm corrected.

    Returns a long DataFrame: problem, tau_t, contrast, n_pairs, p_raw,
    p_holm, mean_ref, mean_other, better.
    """
    rows = []
    for (prob, tau), sub in runs.groupby(["problem", "tau_t"]):
        piv = sub.pivot_table(index="seed", columns="algorithm",
                              values="MIGD", aggfunc="mean")
        if ref not in piv.columns:
            continue
        pvals, entries = [], []
        for other in others:
            if other not in piv.columns:
                continue
            pair = piv[[ref, other]].dropna()
            if len(pair) < 5:
                continue
            try:
                stat, p = stats.wilcoxon(pair[ref].to_numpy(),
                                         pair[other].to_numpy(),
                                         zero_method="wilcox",
                                         alternative="two-sided")
            except ValueError:
                # all differences zero -> identical
                stat, p = 0.0, 1.0
            pvals.append(float(p))
            entries.append({
                "problem": prob, "tau_t": int(tau),
                "contrast": f"{ref} vs {other}", "other": other,
                "n_pairs": int(len(pair)), "p_raw": float(p),
                "mean_ref": float(pair[ref].mean()),
                "mean_other": float(pair[other].mean()),
                "median_diff": float((pair[ref] - pair[other]).median()),
                "better": ref if pair[ref].mean() < pair[other].mean()
                          else other,
            })
        if not pvals:
            continue
        adj = holm_adjust(pvals)
        for e, a in zip(entries, adj):
            e["p_holm"] = float(a)
            rows.append(e)
    return pd.DataFrame(rows)


DEMSAR_Q = {
    0.10: {2: 1.645, 3: 2.052, 4: 2.291, 5: 2.459, 6: 2.589,
           7: 2.693, 8: 2.780, 9: 2.855, 10: 2.920},
    0.05: {2: 1.960, 3: 2.343, 4: 2.569, 5: 2.728, 6: 2.850,
           7: 2.949, 8: 3.031, 9: 3.102, 10: 3.164},
}


def friedman_nemenyi(runs, algos, alpha=0.05):
    """Friedman chi2 + Nemenyi CD per (problem, tau) and overall average ranks.

    q values are Demšar (2006) Table 5; k > 10 falls back to the k = 10 row
    (conservative), matching evaluation/friedman_test.py.
    """
    q_tables = DEMSAR_Q
    per_cell = []
    for (prob, tau), sub in runs.groupby(["problem", "tau_t"]):
        piv = sub.pivot_table(index="seed", columns="algorithm",
                              values="MIGD", aggfunc="mean")
        cols = [a for a in algos if a in piv.columns]
        mat = piv[cols].dropna()
        if mat.shape[0] < 3 or len(cols) < 3:
            continue
        ranks = mat.rank(axis=1, method="average")
        avg = ranks.mean()
        k, N = len(cols), mat.shape[0]
        chi2, p = stats.friedmanchisquare(
            *[mat[c].to_numpy() for c in cols])
        q = q_tables[alpha][min(k, 10)]
        cd = q * np.sqrt(k * (k + 1) / (6.0 * N))
        per_cell.append({
            "problem": prob, "tau_t": int(tau), "k": k, "N": N,
            "chi2": float(chi2), "p": float(p), "CD": float(cd),
            **{f"rank_{c}": float(avg[c]) for c in cols},
        })
    cells = pd.DataFrame(per_cell)
    if cells.empty:
        return cells, pd.DataFrame(), pd.DataFrame()
    rank_cols = [c for c in cells.columns if c.startswith("rank_")]
    overall = cells[rank_cols].mean().rename(
        lambda s: s.replace("rank_", "")).sort_values()
    # pairwise rank gaps vs the best-ranked algorithm, with CD comparison
    gaps = []
    for _, row in cells.iterrows():
        ranks = {c.replace("rank_", ""): row[c] for c in rank_cols}
        for a, b in combinations(sorted(ranks), 2):
            gaps.append({"problem": row["problem"], "tau_t": row["tau_t"],
                         "pair": f"{a} vs {b}", "rank_gap": abs(ranks[a] - ranks[b]),
                         "CD": row["CD"],
                         "indistinguishable": abs(ranks[a] - ranks[b]) <= row["CD"]})
    return cells, overall.to_frame("mean_rank"), pd.DataFrame(gaps)


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main():
    os.makedirs(OUT, exist_ok=True)
    print("loading raw data from:", RAW_DIRS)
    raw = load_all(RAW_DIRS)
    print(f"  {len(raw):,} per-generation rows, "
          f"{raw['_file'].nunique()} files, "
          f"code_version={sorted(raw['code_version'].astype(str).unique())}")

    # ---- integrity ----
    combos = raw.groupby(RUN_KEYS).size().reset_index(name="n_rows")
    distinct = raw.groupby(["algorithm", "problem", "tau_t",
                            "population_size", "n_t", "seed"]).ngroups
    # run inventory: one row per (directory, experiment kind, algorithm,
    # problem, tau, N, n_t) with the number of distinct seeds actually present
    by_prefix = (raw.groupby(["_source", "_kind", "algorithm", "problem",
                              "tau_t", "population_size", "n_t"],
                             dropna=False)
                    .agg(n_runs=("seed", "nunique"),
                         n_rows=("IGD", "size"))
                    .reset_index())
    by_prefix.to_csv(os.path.join(OUT, "integrity_runs_by_prefix.csv"),
                     index=False)
    integrity = {
        "per_generation_rows": int(len(raw)),
        "distinct_runs": int(distinct),
        "rows_per_run_unique": sorted(combos["n_rows"].unique().tolist()),
        "files": int(raw["_file"].nunique()),
        "code_versions": sorted(raw["code_version"].astype(str).unique()),
    }
    print("integrity:", json.dumps(integrity, ensure_ascii=False))

    runs = per_run_metrics(raw)
    runs.to_csv(os.path.join(OUT, "per_run_metrics.csv"), index=False)
    rec = recovery_generations(raw)
    rec.to_csv(os.path.join(OUT, "recovery_generations_raw.csv"), index=False)

    report = []

    def emit(s=""):
        report.append(s)
        print(s)

    emit("# v2.5 thesis analysis")
    emit("")
    emit(f"- per-generation rows: {integrity['per_generation_rows']:,}")
    emit(f"- distinct runs: {integrity['distinct_runs']:,}")
    emit(f"- files: {integrity['files']}")
    emit(f"- code versions: {', '.join(integrity['code_versions'])}")
    emit("")

    # ---- v2.4 comparison (frozen data only) ----
    # MAIN_COMPARISON_KEYS below is the set of (source, experiment-group) pairs
    # that constitute the two controlled comparison data sets.  The sensitivity
    # group is EXCLUDED on purpose: its rho/alpha files are named
    # sensitivity_rho_*_tau25_rho<r> and sensitivity_alpha_*_tau25_alpha<a>,
    # which produce the SAME (algorithm, problem, tau_t, N, n_t) key as the
    # main FDA tau=25 comparison runs and would otherwise be merged into the
    # same group (producing bogus n=40 groups).
    MAIN_KINDS = ("comparison", "supplementary")
    main = runs[runs["_kind"].isin(MAIN_KINDS)]
    v24 = main[main["_source"] == "results/raw"]
    v25 = main[main["_source"] == "results/raw_v25"]

    for tag, sub, algos in (("v24", v24, V24_ALGOS),
                            ("v25", v25, V25_ALGOS)):
        if sub.empty:
            continue
        s = summarize(sub, algos)
        s.to_csv(os.path.join(OUT, f"summary_{tag}.csv"), index=False)
        emit(f"## comparison summary [{tag}]")
        emit("")
        emit("| problem | tau | algorithm | MIGD mean | MIGD std | MHV mean | "
             "post-change MIGD | mig_time_s | n |")
        emit("|---|---|---|---|---|---|---|---|---|")
        for _, r in s.sort_values(["problem", "tau_t", "MIGD_mean"]).iterrows():
            emit(f"| {r['problem']} | {int(r['tau_t'])} | "
                 f"{ALGO_LABEL.get(r['algorithm'], r['algorithm'])} | "
                 f"{r['MIGD_mean']:.4f} | {r['MIGD_std']:.4f} | "
                 f"{r['MHV_mean']:.4f} | {r['post_change_MIGD']:.4f} | "
                 f"{r['mig_time_mean']:.5f} | {int(r['n_runs'])} |")
        emit("")

        # ---- statistical tests ----
        ref = "ldt_nsga2" if tag == "v24" else "ldt_alpha1"
        others = [a for a in algos if a != ref]
        tests = paired_wilcoxon_holm(sub, ref, others)
        tests.to_csv(os.path.join(OUT, f"wilcoxon_holm_{tag}.csv"), index=False)
        emit(f"### Wilcoxon (paired) + Holm, ref = {ref} [{tag}]")
        emit("")
        emit("| problem | tau | contrast | n | mean ref | mean other | "
             "p_raw | p_holm | better |")
        emit("|---|---|---|---|---|---|---|---|---|")
        for _, r in tests.sort_values(["problem", "tau_t", "p_holm"]).iterrows():
            emit(f"| {r['problem']} | {int(r['tau_t'])} | {r['contrast']} | "
                 f"{int(r['n_pairs'])} | {r['mean_ref']:.4f} | "
                 f"{r['mean_other']:.4f} | {r['p_raw']:.3e} | "
                 f"{r['p_holm']:.3e} | {r['better']} |")
        emit("")

        cells, overall, gaps = friedman_nemenyi(sub, algos)
        if not cells.empty:
            cells.to_csv(os.path.join(OUT, f"friedman_cells_{tag}.csv"),
                         index=False)
            overall.to_csv(os.path.join(OUT, f"friedman_overall_{tag}.csv"))
            gaps.to_csv(os.path.join(OUT, f"nemenyi_gaps_{tag}.csv"), index=False)
            emit(f"### Friedman + Nemenyi [{tag}]")
            emit("")
            emit("| problem | tau | k | N | chi2 | p | CD | avg ranks |")
            emit("|---|---|---|---|---|---|---|---|")
            rank_cols = [c for c in cells.columns if c.startswith("rank_")]
            for _, r in cells.iterrows():
                rs = ", ".join(f"{c.replace('rank_', '')}={r[c]:.2f}"
                               for c in rank_cols)
                emit(f"| {r['problem']} | {int(r['tau_t'])} | {int(r['k'])} | "
                     f"{int(r['N'])} | {r['chi2']:.2f} | {r['p']:.2e} | "
                     f"{r['CD']:.4f} | {rs} |")
            emit("")
            emit("overall mean rank (all cells averaged): " +
                 ", ".join(f"{i}={v:.3f}" for i, v in overall["mean_rank"].items()))
            emit("")

    # ---- paired alpha comparison on the SAME configuration ----
    if not v25.empty:
        emit("## paired alpha=1.0 (v2.5) vs alpha=0.5 (v2.4)")
        emit("")
        merged = v25[v25["algorithm"] == "ldt_alpha1"].merge(
            v24[v24["algorithm"] == "ldt_nsga2"],
            on=["problem", "tau_t", "seed"], suffixes=("_a1", "_a05"))
        rows = []
        for (prob, tau), sub in merged.groupby(["problem", "tau_t"]):
            a1, a05 = sub["MIGD_a1"], sub["MIGD_a05"]
            stat, p = stats.wilcoxon(a1, a05, alternative="two-sided")
            rows.append({"problem": prob, "tau_t": int(tau),
                         "n": len(sub),
                         "MIGD_alpha1": a1.mean(), "std_alpha1": a1.std(),
                         "MIGD_alpha05": a05.mean(), "std_alpha05": a05.std(),
                         "improvement_pct": 100 * (a05.mean() - a1.mean()) / a05.mean(),
                         "p": float(p),
                         "MHV_alpha1": sub["MHV_a1"].mean(),
                         "MHV_alpha05": sub["MHV_a05"].mean()})
        alpha_tbl = pd.DataFrame(rows).sort_values(["problem", "tau_t"])
        alpha_tbl.to_csv(os.path.join(OUT, "alpha1_vs_alpha05.csv"), index=False)
        emit("| problem | tau | n | MIGD a=1.0 | MIGD a=0.5 | improvement | p |")
        emit("|---|---|---|---|---|---|---|")
        for _, r in alpha_tbl.iterrows():
            emit(f"| {r['problem']} | {int(r['tau_t'])} | {int(r['n'])} | "
                 f"{r['MIGD_alpha1']:.4f} | {r['MIGD_alpha05']:.4f} | "
                 f"{r['improvement_pct']:+.1f}% | {r['p']:.3e} |")
        emit("")

        # best-algorithm gap for the v2.5 table
        s = summarize(v25, V25_ALGOS)
        emit("## v2.5 winner per (problem, tau)")
        emit("")
        emit("| problem | tau | best algo | best MIGD | ldt_alpha1 MIGD | gap |")
        emit("|---|---|---|---|---|---|")
        for (prob, tau), sub in s.groupby(["problem", "tau_t"]):
            best = sub.loc[sub["MIGD_mean"].idxmin()]
            own = sub[sub["algorithm"] == "ldt_alpha1"]
            own_mean = float(own["MIGD_mean"].iloc[0]) if len(own) else float("nan")
            gap = 100 * (own_mean - best["MIGD_mean"]) / best["MIGD_mean"]
            emit(f"| {prob} | {int(tau)} | {best['algorithm']} | "
                 f"{best['MIGD_mean']:.4f} | {own_mean:.4f} | {gap:+.1f}% |")
        emit("")

    # ---- recovery summary ----
    if not rec.empty:
        rs = rec.groupby(["problem", "tau_t", "algorithm"]).agg(
            mean_recovery=("recovery_gens", "mean"),
            std_recovery=("recovery_gens", "std")).reset_index()
        rs.to_csv(os.path.join(OUT, "recovery_summary.csv"), index=False)
        emit("## recovery generations (generations until IGD <= previous "
             "environment's final IGD; tau_t if never)")
        emit("")
        emit("| problem | tau | algorithm | mean | std |")
        emit("|---|---|---|---|---|")
        for _, r in rs.sort_values(["problem", "tau_t", "mean_recovery"]).iterrows():
            emit(f"| {r['problem']} | {int(r['tau_t'])} | "
                 f"{ALGO_LABEL.get(r['algorithm'], r['algorithm'])} | "
                 f"{r['mean_recovery']:.2f} | {r['std_recovery']:.2f} |")
        emit("")

    with open(os.path.join(OUT, "analysis_report.md"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(report))
    with open(os.path.join(OUT, "integrity.json"), "w", encoding="utf-8") as f:
        json.dump(integrity, f, indent=2, ensure_ascii=False)
    print(f"\nwrote {OUT}/analysis_report.md")


if __name__ == "__main__":
    main()
