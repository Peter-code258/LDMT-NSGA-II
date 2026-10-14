# -*- coding: utf-8 -*-
"""Multi-resolution convergence check (plan 4.3).

Plan 4.3 asks for "at least a multi-resolution convergence check on the key
problems, confirming that densifying the reference front does not change the
relative ordering of the algorithms".

This script answers exactly that.  It recomputes IGD for a stored per-generation
population ... which we do NOT have (only IGD/HV were written to CSV).  So the
check is done the honest way: it re-RUNS a small, fixed cohort at several
reference-front resolutions and reports whether the algorithm ORDER is stable.

Ordering stability is the property that matters: absolute IGD values are
expected to change with resolution, but the ranking must not.

Usage:
    python -m scripts.multiresolution_check --problems FDA1,FDA3 --seeds 0-4
"""
import argparse
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)
sys.path.insert(0, PROJ)

from algorithms.ldmt_nsga2 import LDTNSGA2, LDTNSGA2Alpha1  # noqa: E402
from algorithms.dnsga2_response import NSGA2Static  # noqa: E402
from algorithms.nsga2_core import DynamicDriver  # noqa: E402
from experiments.common import problem_factory  # noqa: E402

RESOLUTIONS = [250, 500, 1000, 2000]
ARMS = [("static", NSGA2Static), ("ldt_half", LDTNSGA2),
        ("ldt_full", LDTNSGA2Alpha1)]


def steady_migd(problem_name, tau, seed, n_points, algo_factory, n_envs=4,
                N=60):
    p = problem_factory(problem_name, tau_t=tau, n_t=20)
    algo = algo_factory()
    rows = DynamicDriver(algo, p, problem_name, n_envs=n_envs, pop_size=N,
                         run_seed=seed, code_version="v32-mr").run(
        n_points=n_points)
    by_env = defaultdict(dict)
    for r in rows:
        by_env[int(r["environment_index"])][int(r["generation"])] = float(r["IGD"])
    return float(np.mean([by_env[e][max(by_env[e])] for e in by_env]))


def main():
    ap = argparse.ArgumentParser(description="multi-resolution IGD check")
    ap.add_argument("--problems", default="FDA1,FDA3")
    ap.add_argument("--tau", type=int, default=25)
    ap.add_argument("--seeds", default="0-4")
    ap.add_argument("--out", default="results/summary_v32")
    args = ap.parse_args()

    problems = [p.strip() for p in args.problems.split(",") if p.strip()]
    if "-" in args.seeds:
        lo, hi = args.seeds.split("-")
        seeds = list(range(int(lo), int(hi) + 1))
    else:
        seeds = [int(x) for x in args.seeds.split(",")]

    os.makedirs(args.out, exist_ok=True)
    lines = ["# 多分辨率收敛检查（方案 4.3）", "",
             "参考前沿点数加密后，**算法之间的相对排序是否改变**。", ""]
    stable_all = True
    for prob in problems:
        table = {}
        for name, fac in ARMS:
            table[name] = []
            for npts in RESOLUTIONS:
                vals = [steady_migd(prob, args.tau, s, npts, fac)
                        for s in seeds]
                table[name].append(float(np.mean(vals)))
        lines.append(f"## {prob}（τ={args.tau}，{len(seeds)} 种子，"
                     f"n_envs=4，N=60）")
        lines.append("")
        lines.append("| 算法 | " + " | ".join(f"{n} 点" for n in RESOLUTIONS)
                     + " | 排名是否恒定 |")
        lines.append("|---|" + "---|" * (len(RESOLUTIONS) + 1))
        ranks = {n: [] for n in table}
        for i, npts in enumerate(RESOLUTIONS):
            order = sorted(table, key=lambda k: table[k][i])
            for pos, k in enumerate(order, 1):
                ranks[k].append(pos)
        for name in table:
            r = ranks[name]
            stable = len(set(r)) == 1
            stable_all = stable_all and stable
            cells = " | ".join(f"{v:.5f}" for v in table[name])
            lines.append(f"| {name} | {cells} | "
                         f"{'✅ 恒定' if stable else '❌ 改变 ' + str(r)} |")
        lines.append("")
        print(f"[{prob}] ranks per resolution: {ranks}")

    lines.append("## 结论")
    lines.append("")
    if stable_all:
        lines.append("在全部问题上，参考前沿从 250 点加密到 2000 点时"
                     "**算法排序保持恒定**，因此本文的相对结论不是参考前沿分辨率的产物。")
    else:
        lines.append("⚠️ 至少一个问题的排序随分辨率改变——该问题上的相对结论"
                     "**必须**标注为依赖参考前沿分辨率。")
    lines.append("")

    path = os.path.join(args.out, "multiresolution_check.md")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    print()
    print("\n".join(lines))
    print(f"wrote {path}")
    return 0 if stable_all else 1


if __name__ == "__main__":
    sys.exit(main())
