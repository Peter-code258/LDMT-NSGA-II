"""Extract every number the v2.5 defense deck needs, into one JSON file.

The deck generator (``make_defense.py``) consumes this JSON, so no chart value
is ever typed by hand into the HTML.

Outputs
-------
  文档/论文材料/_源文件/defense_data.json
  文档/论文材料/_源文件/defense_data.md      (human-readable dump for checking)

Run:  python -m scripts.export_defense_data
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from evaluation.statistics import holm_adjust                      # noqa: E402
from evaluation.friedman_test import friedman_test, nemenyi_posthoc  # noqa: E402

RAW_V24 = "results/raw"
RAW_V25 = "results/raw_v25"
OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "..", "文档", "论文材料", "_源文件")
JSON_OUT = os.path.join(OUT_DIR, "defense_data.json")
MD_OUT = os.path.join(OUT_DIR, "defense_data.md")

KEYS = ["algorithm", "problem", "tau_t", "population_size", "n_t", "seed"]
LABEL = {
    "nsga2_static": "NSGA-II (static)",
    "dnsga2_inspired": "DNSGA2-inspired",
    "ar1_predictive": "AR1-Predictive",
    "gmm_inspired": "GMM-inspired",
    "pps": "PPS",
    "ldt_nsga2": "LDT (α=0.5)",
    "ldt_alpha1": "LDT (α=1.0)",
}
V24 = ["nsga2_static", "dnsga2_inspired", "ar1_predictive", "gmm_inspired",
       "pps", "ldt_nsga2"]
V25 = ["nsga2_static", "dnsga2_inspired", "ar1_predictive", "gmm_inspired",
       "pps", "ldt_alpha1"]


def load(d, prefix, problems=None, with_study=False):
    """Load raw CSVs with a prefix.  With ``with_study``, sensitivity files get
    ``_study`` (N / n_t / rho / alpha) and ``_value`` columns parsed from the
    filename, because those configs collide on the data columns alone."""
    import re
    frames = []
    if not os.path.isdir(d):
        return pd.DataFrame()
    for f in sorted(os.listdir(d)):
        if not (f.startswith(prefix + "_") and f.endswith(".csv")):
            continue
        df = pd.read_csv(os.path.join(d, f))
        if df.empty or "IGD" not in df.columns:
            continue
        if problems and df["problem"].iloc[0] not in problems:
            continue
        if with_study:
            m = re.search(r"_(N|n_t|rho|alpha)([0-9.]+)\.csv$", f)
            if m:
                df["_study"] = m.group(1)
                df["_value"] = float(m.group(2))
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def steady(df, extra_n_t=None):
    """Steady-state MIGD/MHV + post-change MIGD + migration time per run."""
    if df.empty:
        return pd.DataFrame()
    d = df.copy()
    if extra_n_t is not None:
        d = d[d["n_t"] == extra_n_t]
    idx = d.groupby(KEYS + ["environment_index"])["generation"].idxmax()
    last = d.loc[idx]
    idf = d.groupby(KEYS + ["environment_index"])["generation"].idxmin()
    first = d.loc[idf]
    out = last.groupby(KEYS).agg(MIGD=("IGD", "mean"), MHV=("HV", "mean"),
                                 n_envs=("environment_index", "nunique")
                                 ).reset_index()
    post = first[first["environment_index"] > 0].groupby(KEYS).agg(
        post=("IGD", "mean")).reset_index()
    mig = d.groupby(KEYS).agg(mig_t=("migration_runtime_s", "sum"),
                              code_version=("code_version", "first")).reset_index()
    return out.merge(post, on=KEYS, how="left").merge(mig, on=KEYS, how="left")


def table(runs, problems, algos, taus):
    """{problem: {tau: {algo: {migd, std, mhv, post, mig, n}}}}"""
    res = {}
    for p in problems:
        res[p] = {}
        for t in taus:
            res[p][t] = {}
            for a in algos:
                s = runs[(runs.problem == p) & (runs.tau_t == t) &
                         (runs.algorithm == a)]
                if s.empty:
                    continue
                res[p][t][a] = {
                    "migd": round(float(s.MIGD.mean()), 6),
                    "std": round(float(s.MIGD.std(ddof=1)), 6),
                    "mhv": round(float(s.MHV.mean()), 6),
                    "post": round(float(s.post.mean()), 6),
                    "mig": round(float(s.mig_t.mean()), 6),
                    "n": int(s.seed.nunique()),
                }
    return res


def wilcoxon_vs(runs, problems, taus, ref, others):
    """Paired Wilcoxon + Holm, ref vs each other, per (problem, tau)."""
    out = {}
    for p in problems:
        out[p] = {}
        for t in taus:
            sub = runs[(runs.problem == p) & (runs.tau_t == t)]
            piv = sub.pivot_table(index="seed", columns="algorithm",
                                  values="MIGD", aggfunc="mean")
            if ref not in piv.columns:
                continue
            ps, entries = [], []
            for o in others:
                if o not in piv.columns:
                    continue
                j = piv[[ref, o]].dropna()
                if len(j) < 5:
                    continue
                try:
                    _, pv = stats.wilcoxon(j[ref], j[o], alternative="two-sided")
                except ValueError:
                    pv = 1.0
                ps.append(float(pv))
                entries.append({
                    "other": o, "p_raw": float(pv),
                    "mean_ref": round(float(j[ref].mean()), 6),
                    "mean_other": round(float(j[o].mean()), 6),
                    "better": ref if j[ref].mean() < j[o].mean() else o,
                    "n": int(len(j)),
                })
            if not ps:
                continue
            adj = holm_adjust(ps)
            for e, av in zip(entries, adj):
                e["p_holm"] = float(av)
                out[p][t] = out.get(p, {}).get(t, {})
                out[p][t][e["other"]] = e
    return out


def alpha_pairs(v24, v25):
    """Paired alpha=1.0 vs alpha=0.5 over FDA(v2.4 vs v2.5) + dMOP(v2.5 both)."""
    rows = []
    for (p, t), g in v25.groupby(["problem", "tau_t"]):
        a = g[g.algorithm == "ldt_alpha1"].set_index("seed")["MIGD"]
        b = None
        if g[g.algorithm == "ldt_nsga2"].shape[0]:
            b = g[g.algorithm == "ldt_nsga2"].set_index("seed")["MIGD"]
        if b is None or b.empty:
            if not v24.empty:
                s = v24[(v24.problem == p) & (v24.tau_t == t) &
                        (v24.algorithm == "ldt_nsga2")]
                b = s.set_index("seed")["MIGD"] if len(s) else None
        if b is None or len(b) == 0:
            continue
        j = pd.concat([b.rename("m05"), a.rename("m10")], axis=1).dropna()
        if len(j) < 5:
            continue
        try:
            _, pv = stats.wilcoxon(j["m05"], j["m10"], alternative="two-sided")
        except ValueError:
            pv = 1.0
        m05, m10 = float(j["m05"].mean()), float(j["m10"].mean())
        rows.append({
            "problem": p, "tau": int(t), "n": int(len(j)),
            "m05": round(m05, 6), "m10": round(m10, 6),
            "improve_pct": round(100.0 * (m05 - m10) / m05, 1),
            "p": float(pv), "winner": "a1.0" if m10 < m05 else "a0.5",
        })
    order = {"FDA1": 0, "FDA3": 1, "FDA4": 2, "DMOP1A": 3, "DMOP2A": 4,
             "DMOP3A": 5}
    rows.sort(key=lambda r: (order.get(r["problem"], 99), r["tau"]))
    return rows


def friedman(runs, problems, taus, algos):
    out = {}
    for p in problems:
        for t in taus:
            sub = runs[(runs.problem == p) & (runs.tau_t == t)]
            piv = sub.pivot_table(index="seed", columns="algorithm",
                                  values="MIGD", aggfunc="mean")
            cols = [a for a in algos if a in piv.columns]
            mat = piv[cols].dropna()
            if len(cols) < 3 or mat.shape[0] < 3:
                continue
            f = friedman_test(mat[cols].to_numpy(), higher_is_better=False)
            nem = nemenyi_posthoc(np.asarray(f["rank_matrix"]), alpha=0.05)
            out.setdefault(p, {})[t] = {
                "k": len(cols), "N": int(mat.shape[0]),
                "chi2": round(float(f["chi2"]), 2),
                "p": float(f["p_value"]),
                "CD": round(float(nem["critical_difference"]), 5),
                "ranks": {a: round(float(f["avg_ranks"][i]), 3)
                          for i, a in enumerate(cols)},
            }
    return out


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    cmp_v24 = load(RAW_V24, "comparison")
    cmp_v25 = load(RAW_V25, "comparison")
    supp_v24 = load(RAW_V24, "supplementary")
    supp_v25 = load(RAW_V25, "supplementary")
    abl = load(RAW_V24, "ablation")
    sens = load(RAW_V24, "sensitivity", with_study=True)

    r24 = steady(cmp_v24)
    r25 = steady(cmp_v25)
    # alpha pairs need the dMOP family too, and those live in the supplementary
    # files: v2.4 has alpha=0.5 there, v2.5 has BOTH alpha=1.0 and alpha=0.5.
    r24_all = steady(pd.concat([cmp_v24, supp_v24], ignore_index=True))
    r25_all = steady(pd.concat([cmp_v25, supp_v25], ignore_index=True))
    s24 = steady(supp_v24)          # project DMOP1/2/3 (tau=25 only)
    s25 = steady(supp_v25)          # authoritative DMOP1A/2A/3A (tau=25 only)
    ra = steady(abl)
    # sensitivity: alpha scan and rho scan live in the same key space, so filter
    # by the file-level suffix before aggregating
    sens_alpha = _sensitivity_slice(sens, "alpha")
    sens_rho = _sensitivity_slice(sens, "rho")
    sens_N = _sensitivity_slice(sens, "N")
    sens_nt = _sensitivity_slice(sens, "n_t")

    data = {
        "meta": {
            "version": "v2.5",
            "date": "2026-10",
            "runs_total": 6030,
            "runs_v24": 3780,
            "runs_v25": 2250,
            "tests": 63,
            "cd": 1.37668,
        },
        "labels": LABEL,
        "fda_problems": ["FDA1", "FDA3", "FDA4"],
        "taus": [10, 25, 50],
        "v24_algos": V24,
        "v25_algos": V25,
        "migd_v24": table(r24, ["FDA1", "FDA3", "FDA4"], V24, [10, 25, 50]),
        "migd_v25": table(r25, ["FDA1", "FDA3", "FDA4"], V25, [10, 25, 50]),
        "wilcoxon_v24": wilcoxon_vs(
            r24, ["FDA1", "FDA3", "FDA4"], [10, 25, 50], "ldt_nsga2",
            [a for a in V24 if a != "ldt_nsga2"]),
        "wilcoxon_v25": wilcoxon_vs(
            r25, ["FDA1", "FDA3", "FDA4"], [10, 25, 50], "ldt_alpha1",
            [a for a in V25 if a != "ldt_alpha1"]),
        "friedman_v24": friedman(r24, ["FDA1", "FDA3", "FDA4"], [10, 25, 50], V24),
        "friedman_v25": friedman(r25, ["FDA1", "FDA3", "FDA4"], [10, 25, 50], V25),
        "alpha_pairs": alpha_pairs(r24_all, r25_all),
        "ablation": _ablation(ra),
        "sens_alpha": sens_alpha,
        "sens_rho": sens_rho,
        "sens_N": sens_N,
        "sens_nt": sens_nt,
        "dmop_project": table(s24, ["DMOP1", "DMOP2", "DMOP3"], V24, [25]),
        "dmop_auth": table(s25, ["DMOP1A", "DMOP2A", "DMOP3A"], V24 + ["ldt_alpha1"], [25]),
        "migration_time": _migration(r24),
    }

    with open(JSON_OUT, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=1)

    lines = ["# defense_data (auto-generated by scripts/export_defense_data.py)", ""]
    lines.append(f"- runs_total={data['meta']['runs_total']} "
                 f"tests={data['meta']['tests']} CD={data['meta']['cd']}")
    lines.append("")
    lines.append("## alpha pairs (every configuration)")
    lines.append("| problem | tau | a=0.5 | a=1.0 | improve | p | winner |")
    lines.append("|---|---|---|---|---|---|---|")
    for r in data["alpha_pairs"]:
        lines.append(f"| {r['problem']} | {r['tau']} | {r['m05']:.4f} | "
                     f"{r['m10']:.4f} | {r['improve_pct']:+.1f}% | {r['p']:.2e} | {r['winner']} |")
    lines.append("")
    lines.append("## v2.5 FDA MIGD means")
    for p in data["fda_problems"]:
        for t in data["taus"]:
            cells = data["migd_v25"][p][str(t)] if str(t) in data["migd_v25"][p] else data["migd_v25"][p].get(t, {})
            lines.append(f"- {p} tau={t}: " + ", ".join(
                f"{LABEL.get(a, a)}={v['migd']:.4f}" for a, v in cells.items()))
    with open(MD_OUT, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    print(f"wrote {JSON_OUT} ({os.path.getsize(JSON_OUT):,} bytes)")
    print(f"wrote {MD_OUT}")
    print(f"  alpha_pairs: {len(data['alpha_pairs'])} rows, "
          f"a1.0 wins {sum(1 for r in data['alpha_pairs'] if r['winner']=='a1.0')}")
    return 0


def _sensitivity_slice(sens, kind):
    """Aggregate one sensitivity study using the filename-derived _study tag.

    ``kind`` is one of N / n_t / rho / alpha.  For ``alpha`` the CSV files also
    contain the alpha=0.5 default run (``..._alpha0.5``) so the sweep is fully
    self-contained; for rho the same holds for rho=0.5.
    """
    if sens.empty or "_study" not in sens.columns:
        return []
    d = sens[sens["_study"] == kind]
    if d.empty:
        return []
    rows = []
    for (p, val), g in d.groupby(["problem", "_value"]):
        idx = g.groupby(KEYS + ["environment_index"])["generation"].idxmax()
        last = g.loc[idx]
        per_run = last.groupby(KEYS).agg(MIGD=("IGD", "mean")).reset_index()
        rows.append({"problem": p, "value": float(val),
                     "migd": round(float(per_run.MIGD.mean()), 6),
                     "std": round(float(per_run.MIGD.std(ddof=1)), 6),
                     "n": int(per_run.seed.nunique())})
    rows.sort(key=lambda r: (r["problem"], r["value"]))
    return rows


def _ablation(ra):
    if ra.empty:
        return {}
    out = {}
    for p in ["FDA1", "FDA3", "FDA4"]:
        out[p] = {}
        for a in ["ldt_nsga2", "LDT_AblationA", "LDT_AblationB", "LDT_AblationC"]:
            s = ra[(ra.problem == p) & (ra.algorithm == a)]
            if s.empty:
                continue
            out[p][a] = {"migd": round(float(s.MIGD.mean()), 6),
                         "std": round(float(s.MIGD.std(ddof=1)), 6),
                         "post": round(float(s.post.mean()), 6)}
    return out


def _migration(r24):
    if r24.empty:
        return {}
    out = {}
    for (p, t), g in r24.groupby(["problem", "tau_t"]):
        out.setdefault(p, {})[int(t)] = {
            a: round(float(v), 6)
            for a, v in g.groupby("algorithm")["mig_t"].mean().items()}
    return out


if __name__ == "__main__":
    sys.exit(main())
