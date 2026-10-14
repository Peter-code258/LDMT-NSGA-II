# -*- coding: utf-8 -*-
"""v3.1 analysis: does the online adaptive trend step actually help?

Reads ``results/raw_v31`` and produces, for the CALIBRATION phase, a ranking of
the registered (alpha_min, alpha_max, gamma) candidates, and for the
COMPARISON phase a paired comparison of the selected adaptive rule against the
fixed-alpha variants (0 / 0.5 / 1.0) and the reference baselines on the
held-out problems.

Metrics (definitions follow the thesis):
  steady MIGD   mean over environment segments of the LAST generation's IGD,
                then averaged over seeds
  post-change   mean over changed segments of the FIRST generation's IGD
  recovery      mean generations until IGD returns to the pre-change level
  migration_s   total time spent inside the response strategy
  alpha stats   from the per-row ``alpha_k`` column (parameter stability)

Outputs ``results/summary_v31/`` (tables + a markdown report + JSON).
"""
import csv
import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from evaluation.statistics import holm_adjust  # noqa: E402

OUT_ROOT = "results/raw_v31"
SUMMARY = "results/summary_v31"

# Known arm labels, longest first so prefixes cannot be mis-split.
KNOWN_ARMS = [
    "ldt_adaptive", "ldt_alpha0", "ldt_half", "ldt_full", "pps",
    "nsga2_static",
]


def split_arm(base):
    """Split ``<arm>_<problem>_tau<k>`` into (arm, problem, tau)."""
    for arm in KNOWN_ARMS:
        if base.startswith(arm + "_"):
            rest = base[len(arm) + 1:]
            problem, tau_s = rest.rsplit("_tau", 1)
            return arm, problem, int(tau_s)
    # calibration candidates are named adapt_<a>_<g>
    i = base.find("_DMOP")
    if i < 0:
        i = base.find("_FDA")
    arm, rest = base[:i], base[i + 1:]
    problem, tau_s = rest.rsplit("_tau", 1)
    return arm, problem, int(tau_s)

# the two dMOP problems with a verified dynamic type, for the Type III check
DMOP_TYPES = {"DMOP1A": "III", "DMOP2A": "II", "DMOP3A": "I",
              "FDA1": "I", "FDA3": "II", "FDA4": "I"}


def read_run(path):
    """Return {seed: {'rows': [...]}} for one CSV."""
    per_seed = defaultdict(list)
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r.get("seed"):
                per_seed[int(r["seed"])].append(r)
    return per_seed


def run_metrics(rows):
    """Steady MIGD / post-change MIGD / recovery / migration time / alpha."""
    tau = int(rows[0]["tau_t"])
    by_env = defaultdict(dict)
    mig = 0.0
    alphas = []
    for r in rows:
        e = int(r["environment_index"])
        g = int(r["generation"])
        by_env[e][g] = float(r["IGD"])
        mig += float(r["migration_runtime_s"] or 0.0)
        if r.get("alpha_k") not in (None, ""):
            alphas.append(float(r["alpha_k"]))
    envs = sorted(by_env)
    last = [by_env[e][max(by_env[e])] for e in envs]
    steady = float(np.mean(last))
    # post-change: first generation of every CHANGED environment
    post = [by_env[e][min(by_env[e])] for e in envs if e > 0]
    post_change = float(np.mean(post)) if post else float("nan")
    # recovery: generations until IGD <= previous environment's last generation
    rec = []
    for i in range(1, len(envs)):
        target = last[i - 1]
        gens = sorted(by_env[envs[i]])
        hit = None
        for g in gens:
            if by_env[envs[i]][g] <= target:
                hit = g - gens[0]
                break
        rec.append(tau if hit is None else hit)
    return {
        "steady": steady,
        "post_change": post_change,
        "recovery": float(np.mean(rec)) if rec else float("nan"),
        "migration_s": mig,
        "alpha_mean": float(np.mean(alphas)) if alphas else float("nan"),
        "alpha_std": float(np.std(alphas)) if alphas else float("nan"),
        "alpha_min_seen": float(np.min(alphas)) if alphas else float("nan"),
        "alpha_max_seen": float(np.max(alphas)) if alphas else float("nan"),
    }


def load_phase(phase):
    """{arm: {problem: {tau: {seed: metrics}}}.

    The ``comparison`` and ``endpoint`` phases use the same arms and the same
    seeds, so they are merged: that is what makes the adaptive arm comparable
    against every fixed-alpha arm on ALL six problems with a paired test.
    """
    phases = [phase]
    if phase == "comparison":
        phases.append("endpoint")
    data = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    for ph in phases:
        for f in sorted(glob.glob(os.path.join(OUT_ROOT, f"v31_{ph}_*.csv"))):
            base = os.path.basename(f)[len(f"v31_{ph}_"):-4]
            arm, problem, tau = split_arm(base)
            for seed, rows in read_run(f).items():
                if len(rows) < 2:
                    continue
                data[arm][problem][tau][seed] = run_metrics(rows)
    return data


def mean_metric(data, metric):
    """Mean over problems/taus of the per-(problem,tau) seed-mean."""
    vals = []
    for arm, probs in data.items():
        pass
    return vals


def table_for(data, metric, keys=None):
    """{(arm, problem, tau): mean value}"""
    out = {}
    for arm, probs in data.items():
        for p, taus in probs.items():
            for t, seeds in taus.items():
                if keys and (arm, p, t) not in keys:
                    continue
                out[(arm, p, t)] = float(np.mean(
                    [m[metric] for m in seeds.values()]))
    return out


def ranks_by_config(data, metric):
    """Per (problem, tau) rank of each arm by mean metric (1 = best)."""
    rows = table_for(data, metric)
    by_cfg = defaultdict(dict)
    for (arm, p, t), v in rows.items():
        if np.isfinite(v):
            by_cfg[(p, t)][arm] = v
    rank_sum = defaultdict(float)
    rank_cnt = defaultdict(int)
    for cfg, vals in by_cfg.items():
        order = sorted(vals, key=lambda a: vals[a])
        for i, a in enumerate(order, 1):
            rank_sum[a] += i
            rank_cnt[a] += 1
    return {a: rank_sum[a] / rank_cnt[a] for a in rank_sum}, by_cfg


def paired_wilcoxon(data, metric, arm_a, arm_b):
    """Per-config paired Wilcoxon; returns list of (problem, tau, p, delta%)."""
    rows = []
    for p in sorted(data[arm_a]):
        for t in sorted(data[arm_a][p]):
            if p not in data[arm_b] or t not in data[arm_b][p]:
                continue
            sa = data[arm_a][p][t]
            sb = data[arm_b][p][t]
            seeds = sorted(set(sa) & set(sb))
            if len(seeds) < 5:
                continue
            va = np.array([sa[s][metric] for s in seeds])
            vb = np.array([sb[s][metric] for s in seeds])
            if not (np.all(np.isfinite(va)) and np.all(np.isfinite(vb))):
                continue
            try:
                st, pv = stats.wilcoxon(va, vb)
            except ValueError:
                pv = 1.0
            delta = 100.0 * (np.mean(vb) - np.mean(va)) / np.mean(vb)
            rows.append((p, t, float(pv), float(delta),
                         float(np.mean(va)), float(np.mean(vb))))
    return rows


def main():
    os.makedirs(SUMMARY, exist_ok=True)
    report = []
    out_json = {}

    # ─────────────────────── calibration ───────────────────────
    cal = load_phase("calibration")
    print(f"[calibration] arms={sorted(cal)}")
    rank, by_cfg = ranks_by_config(cal, "steady")
    cal_tbl = table_for(cal, "steady")
    lines = ["### 校准阶段（FDA1, DMOP1A, DMOP3A；τ∈{10,25}；10 种子）", "",
             "| 候选（alpha_min, gamma） | 平均稳态 MIGD | 平均秩 | 含 DMOP1A 的配置 |",
             "|---|---|---|---|"]
    per_problem = defaultdict(dict)
    for (arm, p, t), v in cal_tbl.items():
        per_problem[arm].setdefault(p, []).append(v)
    for arm in sorted(rank, key=lambda a: rank[a]):
        allv = [v for p in per_problem[arm].values() for v in p]
        dm = np.mean(per_problem[arm].get("DMOP1A", [np.nan]))
        lines.append(f"| `{arm}` | {np.mean(allv):.5f} | **{rank[arm]:.2f}** | "
                     f"{dm:.5f} |")
    lines.append("")
    best_arm = min(rank, key=lambda a: rank[a])
    lines.append(f"校准阶段按平均秩选出：**`{best_arm}`**"
                 f"（平均秩 {rank[best_arm]:.2f}，平均稳态 MIGD "
                 f"{np.mean([v for p in per_problem[best_arm].values() for v in p]):.5f}）")
    lines.append("")
    report += lines
    out_json["calibration_rank"] = {a: rank[a] for a in rank}
    out_json["calibration_selected"] = best_arm

    # ─────────────────────── comparison ───────────────────────
    comp = load_phase("comparison")
    print(f"[comparison] arms={sorted(comp)}")
    rank_c, cfg_c = ranks_by_config(comp, "steady")
    lines = ["### 对照阶段（held-out：FDA3, FDA4, DMOP2A；τ∈{10,25,50}；30 种子）", "",
             "| 算法 | 平均稳态 MIGD | 平均秩 | 变化后 MIGD | 恢复代数 | 迁移耗时(s) |",
             "|---|---|---|---|---|---|"]
    for arm in sorted(rank_c, key=lambda a: rank_c[a]):
        t_steady = table_for(comp, "steady")
        t_post = table_for(comp, "post_change")
        t_rec = table_for(comp, "recovery")
        t_mig = table_for(comp, "migration_s")
        def mm(tbl, a):
            vs = [v for (x, _, _), v in tbl.items() if x == a and np.isfinite(v)]
            return float(np.mean(vs)) if vs else float("nan")
        lines.append(f"| `{arm}` | {mm(t_steady, arm):.5f} | **{rank_c[arm]:.2f}** "
                     f"| {mm(t_post, arm):.5f} | {mm(t_rec, arm):.2f} "
                     f"| {mm(t_mig, arm):.4f} |")
    lines.append("")
    report += lines
    out_json["comparison_rank"] = {a: rank_c[a] for a in rank_c}

    # per-problem detail
    lines = ["### 对照阶段逐问题稳态 MIGD", "",
             "| 问题 | 类型 | " + " | ".join(sorted(rank_c)) + " |",
             "|---|---|" + "---|" * len(rank_c)]
    t_steady = table_for(comp, "steady")
    for p in sorted({k[1] for k in t_steady}):
        cells = []
        for a in sorted(rank_c):
            vs = [v for (x, pp, _), v in t_steady.items() if x == a and pp == p]
            cells.append(f"{np.mean(vs):.5f}" if vs else "-")
        lines.append(f"| {p} | {DMOP_TYPES.get(p, '?')} | " + " | ".join(cells) + " |")
    lines.append("")
    lines += ["### alpha_k 参数稳定性（自适应臂）", ""]
    t_am = table_for(comp, "alpha_mean")
    t_as = table_for(comp, "alpha_std")
    t_amin = table_for(comp, "alpha_min_seen")
    t_amax = table_for(comp, "alpha_max_seen")
    lines.append("| 问题 | τ | mean | std | min | max |")
    lines.append("|---|---|---|---|---|---|")
    for (arm, p, t) in sorted(t_am):
        if arm != "ldt_adaptive":
            continue
        lines.append(f"| {p} | {t} | {t_am[(arm,p,t)]:.3f} | {t_as[(arm,p,t)]:.3f} "
                     f"| {t_amin[(arm,p,t)]:.3f} | {t_amax[(arm,p,t)]:.3f} |")
    lines.append("")
    report += lines

    # paired tests vs the adaptive arm
    if "ldt_adaptive" in comp:
        for other in sorted(a for a in comp if a != "ldt_adaptive"):
            rows = paired_wilcoxon(comp, "steady", "ldt_adaptive", other)
            if not rows:
                continue
            ps = [r[2] for r in rows]
            adj = holm_adjust(ps)
            lines = [f"#### 自适应 vs `{other}`（配对 Wilcoxon + Holm）", "",
                     "| 问题 | τ | 自适应 | 对照 | 改进 | p_raw | p_holm | 显著 |",
                     "|---|---|---|---|---|---|---|---|"]
            n_sig = 0
            for (p, t, pv, delta, va, vb), pa in zip(rows, adj):
                sig = "✅" if pa < 0.05 else ""
                if pa < 0.05:
                    n_sig += 1
                lines.append(f"| {p} | {t} | {va:.5f} | {vb:.5f} | {delta:+.1f}% "
                             f"| {pv:.3g} | {pa:.3g} | {sig} |")
            lines.append("")
            lines.append(f"共 {len(rows)} 个配置，其中 **{n_sig}** 个差异显著（Holm 校正后 p<0.05）。")
            lines.append("")
            report += lines
            out_json[f"vs_{other}_significant"] = n_sig
            out_json[f"vs_{other}_n"] = len(rows)

    md = "\n".join(report)
    with open(os.path.join(SUMMARY, "adaptive_alpha_report.md"), "w",
              encoding="utf-8") as fh:
        fh.write("# v3.1 自适应趋势步长：实验结果\n\n" + md)
    with open(os.path.join(SUMMARY, "adaptive_alpha.json"), "w",
              encoding="utf-8") as fh:
        json.dump(out_json, fh, ensure_ascii=False, indent=2)
    print(md)
    print(f"\nwrote {SUMMARY}/adaptive_alpha_report.md")


if __name__ == "__main__":
    main()
