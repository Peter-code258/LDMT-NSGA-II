"""Ablation study per opening report section 7.7.

Variants:
  - Full LDT        : all modules enabled (rho=0.5, kappa=0.2, random=0.3, alpha=0.5)
  - Ablation-A      : disable trend extrapolation (alpha=0, mu_hat = mu_e)
  - Ablation-B      : disable distribution transfer (use old pop + random candidates)
  - Ablation-C      : kappa=0 (no old population retention)
"""
import argparse
import os
import sys
import csv

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from algorithms.ldmt_nsga2 import LDTNSGA2
from algorithms.distribution_features import (
    normalize, denormalize, extract_features, non_dominated_set)
from algorithms.environment_manager import ChangeDetector
from experiments.common import (RunConfig, ExperimentRunner, problem_factory,
                                 ROW_FIELDS)


class LDTAblationA(LDTNSGA2):
    """Ablation-A: disable trend extrapolation (alpha=0, mu_hat = mu_e)."""
    name = "LDT_AblationA"

    def __init__(self, **kw):
        super().__init__(trend_alpha=0.0, **kw)


class LDTAblationB(LDTNSGA2):
    """Ablation-B: disable distribution transfer, use old pop + random mix."""
    name = "LDT_AblationB"

    def respond(self, prev_X, prev_F, problem, env_index, rng):
        N, D = prev_X.shape
        xl = np.asarray(problem.xl, dtype=float)
        xu = np.asarray(problem.xu, dtype=float)
        nd_mask = non_dominated_set(prev_X, prev_F)
        nd_X = prev_X[nd_mask] if nd_mask.sum() >= max(10, N // 4) else prev_X
        n_nd = int(nd_X.shape[0])
        mu_e, sigma_e = extract_features(nd_X, xl, xu)
        self._ensure_probe(prev_X, problem)
        Dk = self._detect_change(problem)
        # no distribution transfer: reuse old population directly
        n_transfer = int(round(self.transfer_ratio * N))
        n_retain = int(round(self.retain_ratio * N))
        n_random = int(round(self.random_ratio * N))
        # transfer part = random samples from uniform (no distribution info)
        transfer_part = rng.uniform(xl, xu, size=(n_transfer, D))
        idx = rng.permutation(N)
        retain_part = prev_X[idx[:n_retain]].copy()
        random_part = rng.uniform(xl, xu, size=(n_random, D))
        new_X = np.concatenate([transfer_part, retain_part, random_part], axis=0)
        new_X = self._resize(new_X, N, xl, xu, rng)
        bclip = float(np.mean(np.any(
            (new_X <= xl + 1e-9) | (new_X >= xu - 1e-9), axis=1)))
        return new_X, {
            "lambda_k": 1.0,
            "boundary_clip_ratio": bclip,
            "n_nd": n_nd,
            "Dk": float(Dk),
        }


class LDTAblationC(LDTNSGA2):
    """Ablation-C: kappa=0 (no old population retention)."""
    name = "LDT_AblationC"

    def __init__(self, **kw):
        super().__init__(retain_ratio=0.0, random_ratio=0.5, **kw)


ABLATION_VARIANTS = {
    "full": LDTNSGA2,
    "ablation_a": LDTAblationA,
    "ablation_b": LDTAblationB,
    "ablation_c": LDTAblationC,
}


def run_ablation_variant(variant_name, problem_name, cfg, seed):
    from algorithms import DynamicDriver
    problem = problem_factory(problem_name, tau_t=cfg.tau_t, n_t=cfg.n_t)
    algo = ABLATION_VARIANTS[variant_name]()
    n_envs = max(1, cfg.n_gen // cfg.tau_t)
    driver = DynamicDriver(algo, problem, problem_name, n_envs=n_envs,
                           pop_size=cfg.N, run_seed=seed)
    return driver.run(n_points=500)


def main():
    ap = argparse.ArgumentParser(description="LDT-NSGA-II ablation")
    ap.add_argument("--problems", default="FDA1,FDA3,FDA4")
    ap.add_argument("--tau", type=int, default=25)
    ap.add_argument("--seeds", default="0-29")
    ap.add_argument("--out", default="results/raw")
    ap.add_argument("--n-gen", type=int, default=300)
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--n-t", type=int, default=20)
    ap.add_argument("--limit-seeds", type=int, default=None)
    args = ap.parse_args()

    problems = [p.strip() for p in args.problems.split(",")]
    seeds = ([int(s) for s in args.seeds.split(",")]
             if "-" not in args.seeds
             else list(range(int(args.seeds.split("-")[0]),
                              int(args.seeds.split("-")[1]) + 1)))
    if args.limit_seeds is not None and args.limit_seeds > 0:
        seeds = seeds[:args.limit_seeds]

    cfg = RunConfig(problem=problems[0], tau_t=args.tau, n_t=args.n_t,
                    N=args.n, n_gen=args.n_gen, seeds=range(30))
    runner = ExperimentRunner(out_dir=args.out)
    print(f"Ablation: variants={list(ABLATION_VARIANTS)} problems={problems} "
          f"tau={args.tau} seeds={seeds}")
    for v in ABLATION_VARIANTS:
        for p in problems:
            for s in seeds:
                path = os.path.join(args.out, f"ablation_{v}_{p}_tau{args.tau}.csv")
                if os.path.exists(path):
                    done = set()
                    with open(path) as f:
                        for row in csv.DictReader(f):
                            if row.get("seed"):
                                done.add(int(row["seed"]))
                    if s in done:
                        continue
                try:
                    rows = run_ablation_variant(v, p, cfg, s)
                    write_header = not os.path.exists(path)
                    with open(path, "a", newline="") as f:
                        w = csv.DictWriter(f, fieldnames=ROW_FIELDS)
                        if write_header:
                            w.writeheader()
                        for r in rows:
                            w.writerow({k: r.get(k, "") for k in ROW_FIELDS})
                except Exception as e:
                    print(f"  [FAIL] {v}/{p}/seed{s}: {e}")
                    continue
    print("Ablation done. CSVs in", args.out)


if __name__ == "__main__":
    main()
