"""Common experiment utilities: RunConfig, problem factory, single-run driver,
and ExperimentRunner with resumable CSV output + validation.

v2.2 (P1-3.1 / 3.6):
- Registry maps user-facing CLI names to algorithm CLASSES; each class
  carries a unique canonical id (``.name``). CANONICAL_ALGOS is the default
  experiment list (no aliases -> no duplicate counting).
- run_single returns (rows, meta); ExperimentRunner appends per-run metadata
  (seed, git commit, python/dependency versions, timings, n_eval_total,
  status) to ``run_meta.jsonl`` and failures to ``failed_runs.jsonl``.
"""
import os
import csv
import json
import time
import platform
import subprocess
import dataclasses
from typing import List, Optional

import numpy as np

from problems import (FDA1, FDA3, FDA4, DMOP1, DMOP2, DMOP3,
                      DMOP1A, DMOP2A, DMOP3A)
from algorithms import (NSGA2Static, DNSGA2Response, ARNSGA2, GMMTransfer,
                         LDTNSGA2, LDTNSGA2Alpha1, PPSNSGA2, DynamicDriver)

PROBLEM_REGISTRY = {
    "FDA1": FDA1,
    "FDA3": FDA3,
    "FDA4": FDA4,
    "DMOP1": DMOP1,
    "DMOP2": DMOP2,
    "DMOP3": DMOP3,
    # v2.5: authoritative dMOP1/2/3 of Goh & Tan (2009)
    "DMOP1A": DMOP1A,
    "DMOP2A": DMOP2A,
    "DMOP3A": DMOP3A,
}

# CLI name -> class. Each class has ONE canonical id (.name attribute).
ALGORITHM_REGISTRY = {
    # canonical ids are valid CLI names (v2.2 default experiment list)
    "nsga2_static": NSGA2Static,
    "dnsga2_inspired": DNSGA2Response,
    "ar1_predictive": ARNSGA2,
    "gmm_inspired": GMMTransfer,
    "pps": PPSNSGA2,
    "ldt_nsga2": LDTNSGA2,
    # v2.5 experimental variant: LDT with FULL-STEP trend extrapolation
    # (alpha=1.0). Deliberately NOT part of CANONICAL_ALGOS so that every
    # existing v2.4 artefact and default experiment list stays unchanged.
    "ldt_alpha1": LDTNSGA2Alpha1,
    # human-readable CLI names
    "NSGA-II": NSGA2Static,               # static baseline (no migrants)
    "DNSGA-II-inspired": DNSGA2Response,  # random-immigrant response
    "AR1-Predictive": ARNSGA2,
    "GMM-inspired": GMMTransfer,
    "LDT-NSGA-II": LDTNSGA2,
    "PPS": PPSNSGA2,
    # backward-compatible aliases (CLI input only, NOT in default lists)
    "DNSGA2": DNSGA2Response,
    "AR_NSGA2": ARNSGA2,
    "GMM_Transfer": GMMTransfer,
    "LDMT_NSGA2": LDTNSGA2,
    "LDT_NSGA2": LDTNSGA2,
}

# Canonical experiment list (one entry per distinct behaviour).
CANONICAL_ALGOS = ["nsga2_static", "dnsga2_inspired", "ar1_predictive",
                   "gmm_inspired", "pps", "ldt_nsga2"]

# canonical id -> display name (for plots / tables)
ALGO_DISPLAY = {
    "nsga2_static": "NSGA-II (static)",
    "dnsga2_inspired": "DNSGA-II-inspired",
    "ar1_predictive": "AR1-Predictive",
    "gmm_inspired": "GMM-inspired",
    "pps": "PPS",
    "ldt_nsga2": "LDT-NSGA-II",
}

# legacy CSV labels -> canonical ids (analyze-time normalization only)
LEGACY_ALGO_MAP = {
    "DNSGA2": "dnsga2_inspired",      # pre-v2.2 runs used the immigrant variant
    "NSGA-II": "nsga2_static",
    "DNSGA-II-inspired": "dnsga2_inspired",
    "AR_NSGA2": "ar1_predictive",
    "AR1-Predictive": "ar1_predictive",
    "GMM_Transfer": "gmm_inspired",
    "GMM-inspired": "gmm_inspired",
    "LDMT_NSGA2": "ldt_nsga2",
    "LDT_NSGA2": "ldt_nsga2",
    "LDT-NSGA-II": "ldt_nsga2",
    "PPS": "pps",
    "PPS_NSGA2": "pps",
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
        default_factory=lambda: list(CANONICAL_ALGOS))
    seeds: range = range(30)
    code_version: str = "v2.4"
    # v2.4 (K6 fix): the unused ``evals_per_env`` field was removed; the real
    # evaluation budget per environment is N * tau_t (driver-controlled).


def problem_factory(problem_name, tau_t=25, n_t=20):
    cls = PROBLEM_REGISTRY[problem_name.upper()]
    return cls(tau_t=tau_t, n_t=n_t)


def canonical_id(name):
    """Resolve any CLI/alias name to the class' canonical id."""
    cls = ALGORITHM_REGISTRY[_ci_lookup(ALGORITHM_REGISTRY, name)]
    return getattr(cls, "name", cls.__name__)


def build_algorithm(name, **kw):
    key = _ci_lookup(ALGORITHM_REGISTRY, name)
    return ALGORITHM_REGISTRY[key](**kw)


def _ci_lookup(registry, name):
    n = name.lower()
    for k in registry:
        if k.lower() == n:
            return k
    raise KeyError(f"Unknown algorithm: {name}")


def git_commit():
    """Best-effort git commit hash (empty string when unavailable)."""
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
            timeout=5, cwd=os.path.dirname(os.path.abspath(__file__)))
        if out.returncode == 0:
            return out.stdout.strip()[:12]
    except Exception:
        pass
    return "unknown"


def environment_meta():
    """Runtime environment metadata for reproducibility (P1-3.6)."""
    import sys
    try:
        import pymoo
        pymoo_v = getattr(pymoo, "__version__", "unknown")
    except Exception:
        pymoo_v = "missing"
    try:
        import sklearn
        sklearn_v = getattr(sklearn, "__version__", "unknown")
    except Exception:
        sklearn_v = "missing"
    return {
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "pymoo": pymoo_v,
        "sklearn": sklearn_v,
        "platform": platform.platform(),
        "processor": platform.processor(),
    }


def run_single(algorithm_name, problem_name, cfg, seed, return_meta=False):
    """Run one (algorithm, problem, seed); return rows (+ meta dict)."""
    t0 = time.perf_counter()
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
    meta = {
        "algorithm": getattr(algo, "name", algo.__class__.__name__),
        "algorithm_cli_name": algorithm_name,
        "problem": problem_name,
        "tau_t": cfg.tau_t,
        "n_t": cfg.n_t,
        "N": cfg.N,
        "n_gen": cfg.n_gen,
        "seed": int(seed),
        "n_eval_total": int(driver.n_eval_total),
        "n_eval_per_env": [int(v) for v in driver.n_eval_per_env],
        "runtime_s": round(time.perf_counter() - t0, 3),
        "git_commit": git_commit(),
        "code_version": cfg.code_version,
        "status": "ok",
        "environment": environment_meta(),
    }
    return (rows, meta) if return_meta else rows


class ExperimentRunner:
    """Loops over (problem, tau, algo, seed) combos and writes resumable CSV.

    v2.2: records per-run metadata to ``run_meta.jsonl`` and failures to
    ``failed_runs.jsonl``; ``run()`` returns a summary dict with expected vs
    completed combo counts so callers can exit non-zero on failures.
    """

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

    def _append_jsonl(self, name, record):
        with open(os.path.join(self.out_dir, name), "a",
                  newline="", encoding="utf-8") as f:
            f.write(json.dumps(record, ensure_ascii=False, default=str) + "\n")

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
        n_done = 0
        n_skipped = 0
        failures = []
        for i, (p, t, a, s) in enumerate(combos):
            path = self._path(exp_name, a, p, t)
            done = self._existing_seeds(path)
            if s in done:
                n_skipped += 1
                continue
            sub_cfg = RunConfig(problem=p, tau_t=t, n_t=cfg.n_t, N=cfg.N,
                                n_gen=cfg.n_gen,
                                algorithms=cfg.algorithms,
                                code_version=cfg.code_version)
            try:
                rows, meta = run_single(a, p, sub_cfg, s, return_meta=True)
            except Exception as e:
                print(f"  [FAIL] {a}/{p}/tau{t}/seed{s}: {e}")
                failures.append({"algorithm": a, "problem": p, "tau_t": t,
                                 "seed": s, "error": str(e)})
                self._append_jsonl("failed_runs.jsonl",
                                   {"exp": exp_name, "algorithm": a,
                                    "problem": p, "tau_t": t, "seed": s,
                                    "error": str(e), "git_commit":
                                    git_commit(), "code_version":
                                    cfg.code_version})
                continue
            self._write_rows(path, rows)
            meta["exp"] = exp_name
            self._append_jsonl("run_meta.jsonl", meta)
            n_done += 1
            if (i + 1) % 5 == 0 or i + 1 == len(combos):
                print(f"  {progress_label} {i+1}/{len(combos)} combos done")
        summary = {
            "expected": len(combos),
            "completed": n_done,
            "skipped_resume": n_skipped,
            "failed": len(failures),
            "failures": failures,
        }
        print(f"  [SUMMARY] expected={summary['expected']} "
              f"completed={summary['completed']} "
              f"skipped={summary['skipped_resume']} "
              f"failed={summary['failed']}")
        return summary
