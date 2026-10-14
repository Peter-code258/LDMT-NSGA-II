# -*- coding: utf-8 -*-
"""Audit result staleness: frozen vs freshly recomputed, per (problem, arm, tau).

Sources
-------
frozen : results/raw_v25   (the v2.4/v2.5 artefacts quoted by the thesis)
fresh  : results/raw_v31   (recomputed under identical seeds/budget)

Arm mapping: the thesis' ldt_nsga2 (alpha=0.5) == v31 arm `ldt_half`;
the thesis' ldt_alpha1 (alpha=1.0) == v31 arm `ldt_full`.
"""
import csv
import glob
import os
import sys
from collections import defaultdict

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)


def steady_map(path):
    """{seed: steady-state MIGD}"""
    by = defaultdict(dict)
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if not r.get("seed"):
                continue
            by[int(r["seed"])][(int(r["environment_index"]),
                                int(r["generation"]))] = float(r["IGD"])
    out = {}
    for s, d in by.items():
        e = defaultdict(dict)
        for (env, g), v in d.items():
            e[env][g] = v
        out[s] = float(np.mean([e[k][max(e[k])] for k in e]))
    return out


ARM = {"ldt_nsga2": ("ldt_half", "alpha=0.5"),
       "ldt_alpha1": ("ldt_full", "alpha=1.0"),
       "ldt_adaptive": ("ldt_adaptive", "adaptive")}

KNOWN_PROBLEMS = ["DMOP1A", "DMOP2A", "DMOP3A", "DMOP1", "DMOP2", "DMOP3",
                  "FDA1", "FDA3", "FDA4", "FDA5"]


def parse_frozen(base):
    """'<kind>_<alg>_<problem>_tau<k>' -> (alg, problem, tau).

    Algorithm ids contain underscores, so match against the known id list
    instead of splitting blindly.
    """
    if "_tau" not in base:
        return None
    head, tau_s = base.rsplit("_tau", 1)
    try:
        tau = int(tau_s)
    except ValueError:
        return None
    for alg in ARM:
        if f"_{alg}_" in head:
            return alg, head.split(f"_{alg}_", 1)[1], tau
    for prob in KNOWN_PROBLEMS:
        if head.endswith("_" + prob):
            return None, prob, tau
    return None


print(f"{'problem':<9} {'tau':>4} {'arm':<10} {'frozen':>9} {'fresh':>9} "
      f"{'diff%':>8} {'p':>10}  verdict")
print("-" * 78)

rows = []
frozen_files = sorted(
    glob.glob(os.path.join(PROJ, "results/raw_v25", "*_DMOP*_tau*.csv")) +
    glob.glob(os.path.join(PROJ, "results/raw_v25", "*_FDA*_tau*.csv")))
print(f"(scanning {len(frozen_files)} frozen files)")
for frozen_path in frozen_files:
    base = os.path.basename(frozen_path)[:-4]
    parsed = parse_frozen(base)
    if not parsed or parsed[0] is None:
        continue
    alg, problem, tau = parsed
    v31_arm, label = ARM[alg]

    phase = "endpoint" if problem in ("DMOP1A", "FDA1") else "comparison"
    fresh_path = os.path.join(
        PROJ, "results/raw_v31", f"v31_{phase}_{v31_arm}_{problem}_tau{tau}.csv")
    if not os.path.exists(fresh_path):
        continue

    fr = steady_map(frozen_path)
    nw = steady_map(fresh_path)
    common = sorted(set(fr) & set(nw))
    if len(common) < 5:
        continue
    a = np.array([fr[s] for s in common])
    b = np.array([nw[s] for s in common])
    diff = 100.0 * (a.mean() - b.mean()) / b.mean()
    try:
        _, p = stats.wilcoxon(a, b)
    except ValueError:
        p = 1.0
    verdict = "STALE" if p < 0.01 and abs(diff) > 5 else "ok"
    rows.append((problem, tau, alg, label, a.mean(), b.mean(), diff, p, verdict))
    print(f"{problem:<9} {tau:>4} {alg:<10} {a.mean():>9.5f} {b.mean():>9.5f} "
          f"{diff:>+8.1f} {p:>10.3g}  {verdict}")

print()
stale = [r for r in rows if r[-1] == "STALE"]
print(f"checked {len(rows)} (problem, arm) combinations; "
      f"{len(stale)} look stale:")
for r in stale:
    print(f"  {r[0]} tau={r[1]} {r[2]} ({r[3]}): frozen {r[4]:.5f} vs "
          f"fresh {r[5]:.5f}  ({r[6]:+.1f}%)")
