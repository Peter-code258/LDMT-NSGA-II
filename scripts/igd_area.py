# -*- coding: utf-8 -*-
"""Post-change IGD-curve AREA — the new metric required by plan 4.3.

"新增指标：变化后前 w 代 IGD 曲线面积，用于量化恢复期间的累计损失"

Steady-state MIGD only looks at the last generation of each environment, so it
cannot see HOW PAINFUL the recovery was: two algorithms can end an environment
equally well while one spent five generations lost.  This metric integrates the
IGD excess over the first ``w`` generations after each change:

    A_w = (1 / |changes|) * sum_over_changes  (1/w) * sum_{i=1..w} max(0, IGD_i - IGD_ref)

where ``IGD_ref`` is the previous environment's last-generation IGD (the level
the algorithm had already achieved).  Only the EXCESS above that level counts,
so the metric is "cumulative loss relative to not having been disturbed", not a
restatement of convergence quality.  Lower is better.

Reported per (algorithm, problem, tau) with the usual paired statistics.

Usage:
    python -m scripts.igd_area --csv-dir results/raw --out results/summary_v32
    python -m scripts.igd_area --csv-dir results/raw_v31 --pattern 'v31_v32_*'
"""
import argparse
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

DEFAULT_W = 5


def area_for_run(rows, w=DEFAULT_W):
    """Area under the post-change IGD excess, averaged over changes."""
    tau = int(rows[0]["tau_t"])
    w = max(1, min(int(w), tau))
    by_env = defaultdict(dict)
    for r in rows:
        by_env[int(r["environment_index"])][int(r["generation"])] = float(r["IGD"])
    envs = sorted(by_env)
    if len(envs) < 2:
        return float("nan")
    losses = []
    for i in range(1, len(envs)):
        ref = by_env[envs[i - 1]][max(by_env[envs[i - 1]])]
        gens = sorted(by_env[envs[i]])[:w]
        excess = [max(0.0, by_env[envs[i]][g] - ref) for g in gens]
        losses.append(float(np.mean(excess)))
    return float(np.mean(losses))


def collect(csv_dir, pattern, w=DEFAULT_W):
    """{(arm, problem, tau): {seed: area}} plus rows cache for failure checks."""
    out = defaultdict(dict)
    for f in sorted(glob.glob(os.path.join(csv_dir, pattern))):
        base = os.path.basename(f)[:-4]
        if "_tau" not in base:
            continue
        head, tau_s = base.rsplit("_tau", 1)
        try:
            tau = int(tau_s)
        except ValueError:
            continue
        per_seed = defaultdict(list)
        with open(f, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                if r.get("seed"):
                    per_seed[int(r["seed"])].append(r)
        for seed, rows in per_seed.items():
            if len(rows) < 4:
                continue
            arm = rows[0]["algorithm"]
            problem = rows[0]["problem"]
            out[(arm, problem, tau)][seed] = area_for_run(rows, w)
    return out


def main():
    ap = argparse.ArgumentParser(description="post-change IGD-area metric")
    ap.add_argument("--csv-dir", default="results/raw")
    ap.add_argument("--pattern", default="*.csv")
    ap.add_argument("--w", type=int, default=DEFAULT_W)
    ap.add_argument("--out", default="results/summary_v32")
    ap.add_argument("--reference", default=None,
                    help="arm to compare against (default: ldt_nsga2)")
    args = ap.parse_args()

    data = collect(args.csv_dir, args.pattern, args.w)
    if not data:
        print(f"no CSVs matched {os.path.join(args.csv_dir, args.pattern)}")
        return 2
    os.makedirs(args.out, exist_ok=True)
    ref_arm = args.reference or "ldt_nsga2"

    table = {k: float(np.mean(list(v.values()))) for k, v in data.items()}
    lines = [f"# 变化后前 {args.w} 代 IGD 曲线面积（指标定义见脚本 docstring）", "",
             f"数据源 `{args.csv_dir}/{args.pattern}`；越小越好。", "",
             "| 算法 | 问题 | τ | 面积 | n |", "|---|---|---|---|---|"]
    for (arm, prob, tau) in sorted(table):
        n = len(data[(arm, prob, tau)])
        lines.append(f"| {arm} | {prob} | {tau} | {table[(arm, prob, tau)]:.5f} | {n} |")
    lines.append("")

    # paired comparison against the reference arm
    lines += [f"## 与 `{ref_arm}` 的配对检验（Wilcoxon + Holm）", "",
              "| 问题 | τ | 本臂 | 参照 | 改进 | p_raw | 显著 |",
              "|---|---|---|---|---|---|---|"]
    ps, rows_out = [], []
    for key in sorted(table):
        arm, prob, tau = key
        if arm == ref_arm:
            continue
        rkey = (ref_arm, prob, tau)
        if rkey not in data:
            continue
        common = sorted(set(data[key]) & set(data[rkey]))
        if len(common) < 5:
            continue
        a = np.array([data[key][s] for s in common])
        b = np.array([data[rkey][s] for s in common])
        if not (np.isfinite(a).all() and np.isfinite(b).all()):
            continue
        try:
            _, pv = stats.wilcoxon(a, b)
        except ValueError:
            pv = 1.0
        delta = 100.0 * (b.mean() - a.mean()) / b.mean() if b.mean() else 0.0
        ps.append(pv)
        rows_out.append((prob, tau, a.mean(), b.mean(), delta, pv))
    if ps:
        from evaluation.statistics import holm_adjust
        adj = holm_adjust(ps)
        n_sig = 0
        for (prob, tau, va, vb, delta, pv), pa in zip(rows_out, adj):
            sig = pa < 0.05
            n_sig += int(sig)
            lines.append(f"| {prob} | {tau} | {va:.5f} | {vb:.5f} | {delta:+.1f}% "
                         f"| {pv:.3g} | {'✅' if sig else ''} |")
        lines.append("")
        lines.append(f"共 {len(rows_out)} 个配置，**{n_sig}** 个显著（Holm 校正后 p<0.05）。")
    else:
        lines.append("（无可比较的配对配置）")
    lines.append("")

    path = os.path.join(args.out, f"igd_area_w{args.w}.md")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    with open(os.path.join(args.out, f"igd_area_w{args.w}.json"), "w",
              encoding="utf-8") as fh:
        json.dump({f"{a}|{p}|{t}": v for (a, p, t), v in table.items()}, fh,
                  ensure_ascii=False, indent=2)
    print("\n".join(lines))
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
