# -*- coding: utf-8 -*-
"""Calibrate the movement/noise ratio of the adaptive rule.

Plan 3.1 requires the threshold gamma to be chosen on a calibration run
(under-defaulted for the registered interval), NOT on the final test set.
This script reports the ratio distribution on:

  * a TRANSLATING Pareto-set problem (FDA1)  -> expect ratio >> 1
  * a STATIC   Pareto-set problem (DMOP1A)   -> expect ratio ~ O(1)

so gamma = 1.0 can be judged rather than guessed.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)
sys.path.insert(0, PROJ)

from algorithms.ldmt_nsga2 import LDTNSGA2Adaptive  # noqa: E402
from algorithms.nsga2_core import DynamicDriver, make_nsga2  # noqa: E402
from experiments.common import problem_factory  # noqa: E402


def run_and_collect(problem_name, tau=25, n_envs=6, N=100, seed=0,
                    noise_samples=200):
    """Drive real NSGA-II environments and log the adaptation signals."""
    p = problem_factory(problem_name, tau_t=tau, n_t=20)
    algo = LDTNSGA2Adaptive(noise_samples=noise_samples)
    rows = DynamicDriver(algo, p, problem_name, n_envs=n_envs, pop_size=N,
                         run_seed=seed, code_version="calib").run(n_points=100)
    return algo, rows


print("problem   change  movement   noise     ratio   alpha_k")
print("-" * 60)
summary = {}
for name in ("FDA1", "FDA4", "DMOP1A", "DMOP2A", "DMOP3A"):
    algo, rows = run_and_collect(name)
    for i, (a, mv, nz, r) in enumerate(zip(algo.alpha_series,
                                           algo.movement_series,
                                           algo.noise_series,
                                           algo.ratio_series), 1):
        print(f"{name:<9} {i:<7} {mv:8.5f} {nz:8.5f} {r:8.3f}   {a:.3f}")
    fin = [a for a in algo.alpha_series if np.isfinite(a)]
    rr = [r for r in algo.ratio_series if np.isfinite(r)]
    summary[name] = (np.mean(fin), np.mean(rr), np.max(rr))
    print()

print("=== summary (mean alpha_k, mean ratio, max ratio) ===")
for k, v in summary.items():
    print(f"  {k:<9} alpha={v[0]:.3f}  ratio_mean={v[1]:8.3f}  ratio_max={v[2]:8.3f}")
