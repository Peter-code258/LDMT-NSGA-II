"""Common experiment utilities: RunConfig, problem factory, single-run driver,
and ExperimentRunner with resumable CSV output + validation."""
import os
import csv
import dataclasses
from typing import List, Optional

import numpy as np

from problems import FDA1, FDA3, FDA4, DMOP1, DMOP2, DMOP3
from algorithms import (DNSGA2Response, ARNSGA2, GMMTransfer, LDTNSGA2,
                         PPSNSGA2, DynamicDriver)

PROBLEM_REGISTRY = {
    "FDA1": FDA1,
    "FDA3": FDA3,
    "FDA4": FDA4,
    "DMOP1": DMOP1,
    "DMOP2": DMOP2,
    "DMOP3": DMOP3,
}

ALGORITHM_REGISTRY = {
    "NSGA-II": DNSGA2Response,  # static baseline (no migration)
    "DNSGA-II-inspired": DNSGA2Response,
    "AR1-Predictive": ARNSGA2,
    "GMM-inspired": GMMTransfer,
    "LDT-NSGA-II": LDTNSGA2,
    "PPS": PPSNSGA2,  # population prediction strategy baseline
    # backward-compatible aliases
    "DNSGA2": DNSGA2Response,
    "AR_NSGA2": ARNSGA2,
    "GMM_Transfer": GMMTransfer,
    "LDMT_NSGA2": LDTNSGA2,
}

ROW_FIELDS = [
    "algorithm", "problem", "seed", "population_size", "total_generations",
    "n_t", "tau_t", "environment_index", "generation", "IGD", "HV",
    "generation_runtime_s", "migration_runtime_s", "n_evaluations",
    "code_version", "lambda_k", "boundary_clip_ratio", "n_nd",
]


@dataclasses.dataclass
class RunConfig:
    problem: str
    tau_t: int = 25
    n_t: int = 20
    N: int = 100
    n_gen: int = 300
    algorithms: List[str] = dataclasses.field(
        default_factory=lambda: list(ALGORITHM_REGISTRY.keys()))
    seeds: range = range(30)
    code_version: str = "v1.0"
    # report 7.4: evaluations per environment = 5000
    evals_per_env: int = 5000


def problem_factory(problem_name, tau_t=25, n_t=20):
    cls = PROBLEM_REGISTRY[problem_name.upper()]
    return cls(tau_t=tau_t, n_t=n_t)


def build_algorithm(name, **kw):
    key = _ci_lookup(ALGORITHM_REGISTRY, name)
    return ALGORITHM_REGISTRY[key](**kw)


def _ci_lookup(registry, name):
    n = name.lower()
    for k in registry:
        if k.lower() == n:
            return k
    raise KeyError(f"Unknown algorithm: {name}")


def run_single(algorithm_name, problem_name, cfg, seed):
    """Run one (algorithm, problem, seed) and return a list of row dicts."""
    problem = problem_factory(problem_name, tau_t=cfg.tau_t, n_t=cfg.n_t)
    algo = build_algorithm(algorithm_name)
    n_envs = max(1, cfg.n_gen // cfg.tau_t)
    driver = DynamicDriver(algo, problem, problem_name, n_envs=n_envs,
                           pop_size=cfg.N, run_seed=seed)
    rows = driver.run(n_points=500)
    # validate
    for r in rows:
        for f in ("IGD", "HV"):
            v = r[f]
            if not np.isfinite(v):
                raise ValueError(
                    f"Non-finite {f} in {algorithm_name}/{problem_name}"
                    f"/seed{seed}: {v}")
    return rows


class ExperimentRunner:
    """Loops over (problem, tau, algo, seed) combos and writes resumable CSV."""

    def __init__(self, out_dir="results/raw"):
        self.out_dir = out_dir
        os.makedirs(out_dir, exist_ok=True)

    def _path(self, exp_name, algo, problem, tau):
        return os.path.join(
            self.out_dir,
            f"{exp_name}_{algo}_{problem}_tau{tau}.csv")

    def _existing_seeds(self, path):
        if not os.path.exists(path):
            return set()
        seeds = set()
        with open(path, "r", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                if "seed" in row and row["seed"]:
                    seeds.add(int(row["seed"]))
        return seeds

    def _write_rows(self, path, rows):
        write_header = not os.path.exists(path)
        with open(path, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=ROW_FIELDS)
            if write_header:
                writer.writeheader()
            for r in rows:
                writer.writerow({k: r.get(k, "") for k in ROW_FIELDS})

    def run(self, exp_name, cfg, algorithms=None, problems=None,
            taus=None, seeds=None, limit_seeds=None,
            progress_label=""):
        algos = algorithms or cfg.algorithms
        probs = problems or [cfg.problem]
        ts = taus or [cfg.tau_t]
        if seeds is None:
            seeds = list(cfg.seeds)
        if limit_seeds is not None and limit_seeds > 0:
            seeds = list(seeds)[:limit_seeds]
        combos = [(p, t, a, s) for p in probs for t in ts
                  for a in algos for s in seeds]
        for i, (p, t, a, s) in enumerate(combos):
            path = self._path(exp_name, a, p, t)
            done = self._existing_seeds(path)
            if s in done:
                continue
            sub_cfg = RunConfig(problem=p, tau_t=t, n_t=cfg.n_t, N=cfg.N,
                                n_gen=cfg.n_gen,
                                algorithms=cfg.algorithms,
                                code_version=cfg.code_version)
            try:
                rows = run_single(a, p, sub_cfg, s)
            except Exception as e:
                print(f"  [FAIL] {a}/{p}/tau{t}/seed{s}: {e}")
                continue
            self._write_rows(path, rows)
            if (i + 1) % 5 == 0 or i + 1 == len(combos):
                print(f"  {progress_label} {i+1}/{len(combos)} combos done")
        return self
