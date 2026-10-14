# -*- coding: utf-8 -*-
"""v3.2 analysis: do the §7 roadmap mechanisms actually help?

Compares the seven arms run in the v3.2 phase (1890 runs, 30 seeds, identical
budgets) on FDA1/FDA3/FDA4 x tau in {10,25,50}:

    ldt_half      alpha=0.5 (the v2.x default)
    ldt_full      alpha=1.0
    ldt_adaptive  online adaptive step (v3.1)
    ldt_cov       shrinkage full covariance (plan 3.2)
    ldt_lowrank   low-rank / PCA subspace (plan 3.2)
    ldt_gated     residual-gated transfer budget (plan 3.3)
    ldt_v32       adaptive + gating + multi-source bundle

Reports steady MIGD, post-change MIGD, recovery generations, migration time,
rank, and paired Wilcoxon + Holm against the v2.x default (ldt_half).
"""
import csv
import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)
sys.path.insert(0, PROJ)

from evaluation.statistics import holm_adjust  # noqa: E402
from scripts.igd_area import area_for_run  # noqa: E402

RAW = os.path.join(PROJ, "results", "raw_v31")
OUT = os.path.join(PROJ, "results", "summary_v32")

LABELS = {
    "ldt_half": "LDT \u03b1=0.5（v2.x 默认）",
    "ldt_full": "LDT \u03b1=1.0",
    "ldt_adaptive": "LDT 自适应步长",
    "ldt_cov": "LDT 完整协方差",
    "ldt_lowrank": "LDT 低秩子空间",
    "ldt_gated": "LDT 迁移门控",
    "ldt_v32": "LDT v3.2 组合",
}


def run_metrics(rows):
    tau = int(rows[0]["tau_t"])
    by = defaultdict(dict)
    mig = 0.0
    for r in rows:
        by[int(r["environment_index"])][int(r["generation"])] = float(r["IGD"])
        mig += float(r["migration_runtime_s"] or 0.0)
    envs = sorted(by)
    last = [by[e][max(by[e])] for e in envs]
    steady = float(np.mean(last))
    post = [by[e][min(by[e])] for e in envs if e > 0]
    rec = []
    for i in range(1, len(envs)):
        target = last[i - 1]
        gens = sorted(by[envs[i]])
        hit = next((g - gens[0] for g in gens if by[envs[i]][g] <= target), None)
        rec.append(tau if hit is None else hit)
    return {
        "steady": steady,
        "post_change": float(np.mean(post)) if post else float("nan"),
        "recovery": float(np.mean(rec)) if rec else float("nan"),
        "migration_s": mig,
        "area": area_for_run(rows),
    }


def load():
    data = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    for f in sorted(glob.glob(os.path.join(RAW, "v31_v32_*.csv"))):
        base = os.path.basename(f)[len("v31_v32_"):-4]
        arm, rest = base.split("_", 1)
        while not rest.startswith(("FDA", "DMOP")) and "_" in rest:
            arm += "_" + rest.split("_", 1)[0]
            rest = rest.split("_", 1)[1]
        problem, tau_s = rest.rsplit("_tau", 1)
        per = defaultdict(list)
        with open(f, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                if r.get("seed"):
                    per[int(r["seed"])].append(r)
        for seed, rows in per.items():
            if len(rows) >= 4:
                data[arm][problem][int(tau_s)][seed] = run_metrics(rows)
    return data


def main():
    data = load()
    os.makedirs(OUT, exist_ok=True)
    arms = sorted(data)
    print("arms:", arms)
    metrics = ("steady", "post_change", "recovery", "migration_s", "area")

    # per-arm aggregate
    agg = {a: {m: [] for m in metrics} for a in arms}
    ranks = defaultdict(lambda: defaultdict(list))
    for (p, t) in sorted({(p, t) for a in arms for p in data[a] for t in data[a][p]}):
        per_arm = {}
        for a in arms:
            if p in data[a] and t in data[a][p]:
                seeds = data[a][p][t]
                per_arm[a] = {m: float(np.mean([s[m] for s in seeds.values()]))
                              for m in metrics}
        if len(per_arm) < 2:
            continue
        for i, a in enumerate(sorted(per_arm, key=lambda k: per_arm[k]["steady"]), 1):
            ranks[a]["steady"].append(i)
        for a, vals in per_arm.items():
            for m in metrics:
                agg[a][m].append(vals[m])

    lines = ["# v3.2：研究路线图各机制的实测结果", "",
             f"数据：`results/raw_v31/v31_v32_*.csv`，{len(arms)} 臂 × "
             "FDA1/FDA3/FDA4 × τ∈{10,25,50} × 30 种子。",
             "全部臂共用同一预算 $N\\tau_t$ 每环境，配对同种子。", "",
             "| 算法 | 稳态 MIGD | 平均秩 | 变化后 MIGD | 恢复代数 | 迁移耗时(s) | IGD 面积(w=5) |",
             "|---|---|---|---|---|---|---|"]
    summary = {}
    for a in sorted(arms, key=lambda k: float(np.mean(ranks[k]["steady"] or [99]))):
        r = float(np.mean(ranks[a]["steady"])) if ranks[a]["steady"] else float("nan")
        row = {m: float(np.mean(agg[a][m])) for m in metrics if agg[a][m]}
        summary[a] = dict(row, rank=r)
        lines.append(
            f"| `{a}` {LABELS.get(a, '')} | {row['steady']:.5f} | **{r:.2f}** "
            f"| {row['post_change']:.5f} | {row['recovery']:.2f} "
            f"| {row['migration_s']:.4f} | {row['area']:.5f} |")
    lines.append("")

    # paired tests vs the v2.x default
    base = "ldt_half"
    for other in [a for a in arms if a != base]:
        ps, rows_out = [], []
        for p in sorted(data[base]):
            for t in sorted(data[base][p]):
                if p not in data[other] or t not in data[other][p]:
                    continue
                common = sorted(set(data[base][p][t]) & set(data[other][p][t]))
                if len(common) < 5:
                    continue
                x = np.array([data[other][p][t][s]["steady"] for s in common])
                y = np.array([data[base][p][t][s]["steady"] for s in common])
                if not (np.isfinite(x).all() and np.isfinite(y).all()):
                    continue
                try:
                    _, pv = stats.wilcoxon(x, y)
                except ValueError:
                    pv = 1.0
                ps.append(pv)
                rows_out.append((p, t, x.mean(), y.mean(),
                                 100.0 * (y.mean() - x.mean()) / y.mean(), pv))
        if not ps:
            continue
        adj = holm_adjust(ps)
        n_sig = sum(1 for pa in adj if pa < 0.05)
        n_better = sum(1 for (_, _, xa, ya, d, _), pa in zip(rows_out, adj)
                       if pa < 0.05 and d > 0)
        n_worse = sum(1 for (_, _, xa, ya, d, _), pa in zip(rows_out, adj)
                      if pa < 0.05 and d < 0)
        lines.append(f"### `{other}` vs `{base}`（配对 Wilcoxon + Holm）")
        lines.append("")
        lines.append(f"共 {len(rows_out)} 个配置：**{n_better}** 个显著更好，"
                     f"**{n_worse}** 个显著更差，其余不显著。")
        lines.append("")
        lines.append("| 问题 | τ | 本臂 | 默认 | 改进 | p_holm | 显著 |")
        lines.append("|---|---|---|---|---|---|---|")
        for (p, t, xa, ya, d, pv), pa in zip(rows_out, adj):
            lines.append(f"| {p} | {t} | {xa:.5f} | {ya:.5f} | {d:+.1f}% "
                         f"| {pa:.3g} | {'✅' if pa < 0.05 else ''} |")
        lines.append("")
        summary.setdefault(other, {})["vs_default_better"] = n_better
        summary.setdefault(other, {})["vs_default_worse"] = n_worse

    md = "\n".join(lines)
    with open(os.path.join(OUT, "v32_mechanisms.md"), "w", encoding="utf-8") as fh:
        fh.write(md)
    with open(os.path.join(OUT, "v32_mechanisms.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2, default=str)
    print(md)
    print(f"\nwrote {OUT}/v32_mechanisms.md")


if __name__ == "__main__":
    main()
