"""Ablation study: FDA1/3/4, tau=25, 30 seeds.

Variants (report 5.5):
  - full                          : standard LDMT-NSGA-II.
  - no_feature_simplification     : use a GMM / covariance distribution
                                    description instead of mean-std, but keep
                                    the LDMT transfer + lambda adaptation.
  - no_transfer                   : keep feature extraction but drop the
                                    center/scale mapping (reuse prev pop) +
                                    standard random response.
"""
import argparse
import os
import sys
import copy

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from algorithms.ldmt_nsga2 import LDMTNSGA2
from algorithms.distribution_features import (
    normalize, denormalize, extract_features, non_dominated_set,
    transfer_population)
from experiments.common import (RunConfig, ExperimentRunner, run_single,
                                 ALGORITHM_REGISTRY, ROW_FIELDS)
import csv


class LDMTNoFeatureSimplification(LDMTNSGA2):
    """GMM/covariance distribution description instead of mean-std, but the
    LDMT transfer + lambda adaptation is preserved."""
    name = "LDMT_no_feature"

    def __init__(self, **kw):
        super().__init__(**kw)
        self._gmm = None

    def respond(self, prev_X, prev_F, problem, env_index, rng):
        from sklearn.mixture import GaussianMixture
        N, D = prev_X.shape
        xl = np.asarray(problem.xl, dtype=float)
        xu = np.asarray(problem.xu, dtype=float)
        nd_mask = non_dominated_set(prev_X, prev_F)
        nd_X = prev_X[nd_mask] if nd_mask.sum() >= max(10, N // 4) else prev_X
        mu_prev, sigma_prev = extract_features(nd_X, xl, xu)
        self._ensure_probe(prev_X, D)
        Dk = self._detect_change(problem)
        mu_prev2 = self._mu_prev2 if self._mu_prev2 is not None else mu_prev
        sigma_prev2 = (self._sigma_prev2 if self._sigma_prev2 is not None
                        else sigma_prev)
        D_prev = self._D_prev if self._D_prev is not None else 0.0
        mu_hat, sigma_hat, lam_k = self.migration.predict(
            mu_prev, mu_prev2, sigma_prev, sigma_prev2, Dk, D_prev)
        self.migration.record(mu_prev, sigma_prev, Dk)
        self._mu_prev2 = mu_prev
        self._sigma_prev2 = sigma_prev
        self._D_prev = Dk
        # GMM on normalized ND_X; shift samples by (mu_hat - mu_prev)
        Z = np.clip(normalize(nd_X, xl, xu), 0, 1)
        k = max(1, min(3, len(nd_X) - 1))
        try:
            gmm = GaussianMixture(n_components=k, covariance_type="diag",
                                   reg_covar=1e-4, random_state=42)
            gmm.fit(Z)
            samples = gmm.sample(N)[0]
        except Exception:
            samples = rng.uniform(0, 1, size=(N, D))
        samples = np.clip(samples + (mu_hat - mu_prev), 0, 1)
        transferred = denormalize(samples, xl, xu)
        # 70/20/10 composition (same as full LDMT)
        new_X = self._compose(prev_X, transferred, xl, xu, rng)
        bclip = float(np.mean(np.any(
            (new_X <= xl + 1e-9) | (new_X >= xu - 1e-9), axis=1)))
        return new_X, {"lambda_k": float(lam_k),
                        "boundary_clip_ratio": bclip,
                        "n_nd": int(nd_X.shape[0]), "Dk": float(Dk)}

    def _compose(self, prev_X, transferred, xl, xu, rng):
        N, D = prev_X.shape
        n_transfer = int(round(self.transfer_ratio * N))
        n_mut = int(round(self.local_mut_ratio * N))
        n_imm = int(round(self.immigrant_ratio * N))
        idx = rng.permutation(N)
        transfer_part = transferred[idx[:n_transfer]]
        mut_part = transferred[idx[n_transfer:n_transfer + n_mut]].copy()
        if mut_part.shape[0] > 0:
            n_dims = max(1, D // 3)
            dim_sel = rng.choice(D, size=n_dims, replace=False)
            mut_part[:, dim_sel] += rng.normal(0, 0.05, size=mut_part.shape)[:, dim_sel]
            mut_part = np.clip(mut_part, xl, xu)
        imm_part = rng.uniform(xl, xu, size=(n_imm, D))
        new_X = np.concatenate([transfer_part, mut_part, imm_part], axis=0)
        return self._resize(new_X, N, xl, xu, rng)


class LDMTNoTransfer(LDMTNSGA2):
    """Keep feature extraction but drop the center/scale mapping — reuse the
    previous population + standard random response."""
    name = "LDMT_no_transfer"

    def respond(self, prev_X, prev_F, problem, env_index, rng):
        N, D = prev_X.shape
        xl = np.asarray(problem.xl, dtype=float)
        xu = np.asarray(problem.xu, dtype=float)
        nd_mask = non_dominated_set(prev_X, prev_F)
        nd_X = prev_X[nd_mask] if nd_mask.sum() >= max(10, N // 4) else prev_X
        mu_prev, sigma_prev = extract_features(nd_X, xl, xu)  # features kept
        self._ensure_probe(prev_X, D)
        Dk = self._detect_change(problem)
        # no transfer: reuse prev pop directly (no center/scale mapping)
        transferred = prev_X.copy()
        new_X = self._compose(prev_X, transferred, xl, xu, rng)
        bclip = 0.0
        return new_X, {"lambda_k": 1.0, "boundary_clip_ratio": bclip,
                        "n_nd": int(nd_X.shape[0]), "Dk": float(Dk)}

    def _compose(self, prev_X, transferred, xl, xu, rng):
        N, D = prev_X.shape
        n_transfer = int(round(self.transfer_ratio * N))
        n_mut = int(round(self.local_mut_ratio * N))
        n_imm = int(round(self.immigrant_ratio * N))
        idx = rng.permutation(N)
        transfer_part = transferred[idx[:n_transfer]]
        mut_part = transferred[idx[n_transfer:n_transfer + n_mut]].copy()
        if mut_part.shape[0] > 0:
            n_dims = max(1, D // 3)
            dim_sel = rng.choice(D, size=n_dims, replace=False)
            mut_part[:, dim_sel] += rng.normal(0, 0.05, size=mut_part.shape)[:, dim_sel]
            mut_part = np.clip(mut_part, xl, xu)
        imm_part = rng.uniform(xl, xu, size=(n_imm, D))
        new_X = np.concatenate([transfer_part, mut_part, imm_part], axis=0)
        return self._resize(new_X, N, xl, xu, rng)


ABLATION_VARIANTS = {
    "full": LDMTNSGA2,
    "no_feature_simplification": LDMTNoFeatureSimplification,
    "no_transfer": LDMTNoTransfer,
}


def run_ablation_variant(variant_name, problem_name, cfg, seed):
    from algorithms import DynamicDriver
    from experiments.common import problem_factory
    problem = problem_factory(problem_name, tau_t=cfg.tau_t, n_t=cfg.n_t)
    algo = ABLATION_VARIANTS[variant_name]()
    n_envs = max(1, cfg.n_gen // cfg.tau_t)
    driver = DynamicDriver(algo, problem, problem_name, n_envs=n_envs,
                           pop_size=cfg.N, run_seed=seed)
    return driver.run(n_points=500)


def main():
    ap = argparse.ArgumentParser(description="LDMT-NSGA-II ablation")
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
