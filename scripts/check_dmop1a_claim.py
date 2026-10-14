# -*- coding: utf-8 -*-
"""Cross-check the v2.5 DMOP1A claim that the adaptive result contradicts.

The v2.5 thesis states (conclusion 3):
    "on the authoritative Type III problem dMOP1A, alpha=0.5 beats alpha=1.0 by
     127%  (0.0251 vs 0.0570, p=1.86e-9)"

The v3.1 paired comparison finds the opposite magnitude: adaptive (which takes
alpha ~= 1.0 on this problem) vs a fixed alpha=1.0 differ by 0.1%, and vs a
fixed alpha=0.5 by 0.4%.  This script re-derives the *v2.5* numbers from the
frozen ``results/raw_v25`` CSVs using the SAME metric code as the v3.1 analysis,
so the comparison is apples-to-apples.
"""
import csv
import glob
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)
sys.path.insert(0, PROJ)
from scripts.adaptive_alpha_analysis import run_metrics, read_run  # noqa: E402


def collect(paths):
    out = defaultdict(dict)          # arm -> seed -> metrics
    for f in paths:
        base = os.path.basename(f)
        arm = None
        for a in ("ldt_alpha1", "ldt_nsga2", "ldt_adaptive"):
            if f"_{a}_" in base:
                arm = a
                break
        if arm is None:
            continue
        for seed, rows in read_run(f).items():
            if len(rows) >= 2:
                out[arm][seed] = run_metrics(rows)
    return out


print("=== DMOP1A, tau=25, metric = steady-state MIGD (mean over seeds) ===")
for label, pat in (
    ("frozen v2.5 (results/raw_v25)", "results/raw_v25/*DMOP1A*tau25*.csv"),
    ("frozen v2.4 (results/raw)",     "results/raw/*DMOP1A*"),
    ("new v3.1  (results/raw_v31)",   "results/raw_v31/*DMOP1A*tau25*.csv"),
):
    files = sorted(glob.glob(os.path.join(PROJ, pat)))
    if not files:
        print(f"  {label}: no files")
        continue
    data = collect(files)
    print(f"  {label}: {len(files)} file(s)")
    for arm in sorted(data):
        seeds = sorted(data[arm])
        vals = [data[arm][s]["steady"] for s in seeds]
        print(f"     {arm:<12} n={len(seeds):<3} mean={np.mean(vals):.5f} "
              f"sd={np.std(vals):.5f}")
    # paired test if both present
    if "ldt_alpha1" in data and "ldt_nsga2" in data:
        common = sorted(set(data["ldt_alpha1"]) & set(data["ldt_nsga2"]))
        a = np.array([data["ldt_alpha1"][s]["steady"] for s in common])
        b = np.array([data["ldt_nsga2"][s]["steady"] for s in common])
        from scipy import stats
        try:
            _, p = stats.wilcoxon(a, b)
        except ValueError:
            p = 1.0
        d = 100.0 * (b.mean() - a.mean()) / b.mean()
        print(f"     paired n={len(common)}: alpha1={a.mean():.5f} "
              f"half={b.mean():.5f}  alpha1-vs-half = {d:+.1f}%  p={p:.3g}")
    print()
