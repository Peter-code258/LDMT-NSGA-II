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
import hashlib
import platform
import subprocess
import dataclasses
from typing import List, Optional

import numpy as np

from problems import (FDA1, FDA3, FDA4, DMOP1, DMOP2, DMOP3,
                      DMOP1A, DMOP2A, DMOP3A,
                      IrregularFDA1, ReversedFDA1, GradualFDA1, NoisyFDA1)
from algorithms import (NSGA2Static, DNSGA2Response, ARNSGA2, GMMTransfer,
                         LDTNSGA2, LDTNSGA2Alpha1, LDTNSGA2Adaptive,
                         LDTNSGA2Cov, LDTNSGA2LowRank, LDTNSGA2Gated,
                         LDTNSGA2Full, PPSNSGA2, DynamicDriver)

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
    # v3.2 (plan 4.2): project-defined stress variants of FDA1.  Each relaxes
    # exactly one temporal or measurement assumption; they are NOT literature
    # benchmarks and must be labelled as project variants in any write-up.
    "IRREG_FDA1": IrregularFDA1,
    "REV_FDA1": ReversedFDA1,
    "GRAD_FDA1": GradualFDA1,
    "NOISY_FDA1": NoisyFDA1,
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
    # v3.1 (plan 3.1): LDT with the ONLINE adaptive trend step.  Also NOT part
    # of CANONICAL_ALGOS, so no existing artefact or default list changes.
    "ldt_adaptive": LDTNSGA2Adaptive,
    # v3.2 (plan 3.2/3.3): correlated / low-rank modelling and transfer gating.
    # Also OUTSIDE CANONICAL_ALGOS so no default list changes.
    "ldt_cov": LDTNSGA2Cov,
    "ldt_lowrank": LDTNSGA2LowRank,
    "ldt_gated": LDTNSGA2Gated,
    "ldt_v32": LDTNSGA2Full,
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
    "ldt_alpha1": "LDT-NSGA-II (alpha=1.0)",
    "ldt_adaptive": "LDT-NSGA-II (adaptive alpha)",
    "ldt_cov": "LDT-NSGA-II (full covariance)",
    "ldt_lowrank": "LDT-NSGA-II (low-rank)",
    "ldt_gated": "LDT-NSGA-II (gated transfer)",
    "ldt_v32": "LDT-NSGA-II (v3.2 bundle)",
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
    "code_version", "run_id", "lambda_k", "boundary_clip_ratio", "n_nd",
    # v3.1 (plan 3.1): the trend step actually applied in this environment,
    # so online adaptation is auditable per environment (parameter stability).
    "alpha_k",
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
    # v3.0: default code version for new experiments.  Must be overridden by
    # each driver so the stamp on disk matches the code that really ran.
    code_version: str = "v3.0"
    # v3.0 (plan 2.2/5.1): identifies WHICH experiment a run belongs to, so
    # different experiment groups and parameter sweeps can never be merged by
    # accident (the earlier "sensitivity tau25 merged with main FDA tau25"
    # bug).  A run is uniquely keyed by (experiment_group, algorithm, problem,
    # tau_t, N, n_t, seed) plus any non-default sweep parameters.
    experiment_group: str = "unspecified"
    # Parameters that only some experiments set.  ``None`` means "not
    # applicable to this experiment" and is recorded as an explicit null, so a
    # sweep value can never be confused with "left at default".
    alpha: float = None
    rho: float = None
    extra: dict = dataclasses.field(default_factory=dict)
    # v3.1: constructor kwargs for the algorithm under test, so variant
    # parameters (e.g. the adaptive rule's alpha_min/gamma) travel WITH the
    # configuration instead of needing ad-hoc subclasses in the registry.
    algo_kwargs: dict = dataclasses.field(default_factory=dict)
    # v2.4 (K6 fix): the unused ``evals_per_env`` field was removed; the real
    # evaluation budget per environment is N * tau_t (driver-controlled).

    def config_dict(self):
        """Every parameter that can change a run's behaviour."""
        d = {
            "problem": self.problem,
            "tau_t": int(self.tau_t),
            "n_t": int(self.n_t),
            "N": int(self.N),
            "n_gen": int(self.n_gen),
            "algorithms": sorted(self.algorithms),
            "seed_range": [min(self.seeds), max(self.seeds),
                           len(list(self.seeds))],
            "code_version": self.code_version,
            "experiment_group": self.experiment_group,
            "alpha": self.alpha,
            "rho": self.rho,
        }
        d.update(self.extra or {})
        if self.algo_kwargs:
            d["algo_kwargs"] = {k: self.algo_kwargs[k]
                                for k in sorted(self.algo_kwargs)}
        return d

    def config_hash(self):
        """Stable hash of the behaviour-affecting configuration."""
        blob = json.dumps(self.config_dict(), sort_keys=True,
                          ensure_ascii=True, default=str)
        return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]


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


def git_dirty():
    """True when the working tree has uncommitted changes (or git is absent).

    plan 2.1/5.1: a run that was produced from a dirty tree cannot be traced
    to a single commit, so the fact must be recorded rather than implied.
    """
    try:
        out = subprocess.run(
            ["git", "status", "--porcelain"], capture_output=True, text=True,
            timeout=5, cwd=os.path.dirname(os.path.abspath(__file__)))
        if out.returncode == 0:
            return bool(out.stdout.strip())
    except Exception:
        pass
    return None


def make_run_id(cfg, algorithm_name, problem_name, seed):
    """Immutable identity for one run (plan 5.1).

    Stable across processes and machines: it is derived from the
    behaviour-affecting configuration, not from a clock or a random source, so
    re-running the same command reproduces the same id and duplicates are
    detectable.
    """
    blob = json.dumps({
        "config_hash": cfg.config_hash(),
        "algorithm": str(algorithm_name),
        "problem": str(problem_name),
        "seed": int(seed),
    }, sort_keys=True)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:20]


def run_single(algorithm_name, problem_name, cfg, seed, return_meta=False):
    """Run one (algorithm, problem, seed); return rows (+ meta dict)."""
    t0 = time.perf_counter()
    problem = problem_factory(problem_name, tau_t=cfg.tau_t, n_t=cfg.n_t)
    algo = build_algorithm(algorithm_name, **dict(cfg.algo_kwargs or {}))
    n_envs = max(1, cfg.n_gen // cfg.tau_t)
    run_id = make_run_id(cfg, algorithm_name, problem_name, seed)
    driver = DynamicDriver(algo, problem, problem_name, n_envs=n_envs,
                           pop_size=cfg.N, run_seed=seed,
                           code_version=cfg.code_version)
    driver.run_id = run_id
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
        "run_id": run_id,
        "config_hash": cfg.config_hash(),
        "config": cfg.config_dict(),
        "experiment_group": cfg.experiment_group,
        "algorithm": getattr(algo, "name", algo.__class__.__name__),
        "algorithm_cli_name": algorithm_name,
        "problem": problem_name,
        "tau_t": cfg.tau_t,
        "n_t": cfg.n_t,
        "N": cfg.N,
        "n_gen": cfg.n_gen,
        "seed": int(seed),
        "n_gen_rows": len(rows),
        "n_eval_total": int(driver.n_eval_total),
        "n_eval_per_env": [int(v) for v in driver.n_eval_per_env],
        "runtime_s": round(time.perf_counter() - t0, 3),
        "git_commit": git_commit(),
        "git_dirty": git_dirty(),
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
        # Set by run() so completeness checks know how many rows a full run has.
        self._expected_rows = None

    def _path(self, exp_name, algo, problem, tau):
        return os.path.join(
            self.out_dir,
            f"{exp_name}_{algo}_{problem}_tau{tau}.csv")

    def _seed_row_counts(self, path):
        """Rows-per-seed for one CSV (empty dict when the file is absent)."""
        if not os.path.exists(path):
            return {}
        counts = {}
        with open(path, "r", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                if "seed" in row and row["seed"]:
                    s = int(row["seed"])
                    counts[s] = counts.get(s, 0) + 1
        return counts

    def _existing_seeds(self, path):
        """Seeds present in the CSV *with a complete row set*.

        plan 2.2/5.2: a seed whose rows are truncated (interrupted run, disk
        full, killed process) must NOT count as done -- otherwise a resumed
        experiment silently reports success on partial data.
        """
        counts = self._seed_row_counts(path)
        if self._expected_rows is None:
            return set(counts)
        return {s for s, n in counts.items() if n >= self._expected_rows}

    def _partial_seeds(self, path):
        """Seeds present but with fewer rows than a complete run."""
        if self._expected_rows is None:
            return {}
        counts = self._seed_row_counts(path)
        return {s: n for s, n in counts.items() if n < self._expected_rows}

    def _write_rows(self, path, rows):
        write_header = not os.path.exists(path)
        with open(path, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=ROW_FIELDS)
            if write_header:
                writer.writeheader()
            for r in rows:
                writer.writerow({k: r.get(k, "") for k in ROW_FIELDS})

    def _drop_partial_seeds(self, path):
        """Rewrite ``path`` keeping only fully-complete seeds.

        v3.0 (plan 5.2): a truncated seed must be re-run.  Because rows are
        appended, re-running would otherwise leave the old partial rows in
        place and the file would hold a duplicated, inconsistent run for that
        seed.  Dropping the partial rows first makes the re-run clean.
        """
        if not os.path.exists(path) or self._expected_rows is None:
            return 0
        partial = self._partial_seeds(path)
        if not partial:
            return 0
        keep = []
        with open(path, "r", newline="") as f:
            reader = csv.DictReader(f)
            fields = reader.fieldnames
            for row in reader:
                sd = row.get("seed")
                if not sd:
                    continue
                if int(sd) not in partial:
                    keep.append(row)
        with open(path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields or ROW_FIELDS)
            writer.writeheader()
            for row in keep:
                writer.writerow(row)
        return len(partial)

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
        # v3.0 (plan 5.2): a complete run writes n_gen rows (one per
        # generation); anything short of that is a truncated run.
        self._expected_rows = int(cfg.n_gen)
        combos = [(p, t, a, s) for p in probs for t in ts
                  for a in algos for s in seeds]
        n_done = 0
        n_skipped = 0
        failures = []
        partial_before = {}
        cleaned = set()
        for i, (p, t, a, s) in enumerate(combos):
            path = self._path(exp_name, a, p, t)
            # v3.0 (plan 5.2): clear truncated rows for this file once, so a
            # re-run of a partial seed appends cleanly instead of duplicating.
            if path not in cleaned:
                n_bad = self._drop_partial_seeds(path)
                if n_bad:
                    print(f"  [RESUME] dropped {n_bad} truncated seed(s) "
                          f"from {os.path.basename(path)}")
                cleaned.add(path)
            done = self._existing_seeds(path)
            if s in done:
                n_skipped += 1
                continue
            sub_cfg = RunConfig(problem=p, tau_t=t, n_t=cfg.n_t, N=cfg.N,
                                n_gen=cfg.n_gen,
                                algorithms=cfg.algorithms,
                                code_version=cfg.code_version,
                                experiment_group=cfg.experiment_group,
                                alpha=cfg.alpha, rho=cfg.rho,
                                extra=dict(cfg.extra or {}),
                                algo_kwargs=dict(cfg.algo_kwargs or {}))
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
                                    cfg.code_version,
                                    "experiment_group": cfg.experiment_group})
                continue
            self._write_rows(path, rows)
            meta["exp"] = exp_name
            self._append_jsonl("run_meta.jsonl", meta)
            n_done += 1
            if (i + 1) % 5 == 0 or i + 1 == len(combos):
                print(f"  {progress_label} {i+1}/{len(combos)} combos done")
        # v3.0 (plan 5.2): report every category explicitly.  ``missing`` is
        # the count of combos with no data at all; ``partial`` counts seeds
        # that exist but are truncated.  Either being non-zero means the
        # experiment as a whole is NOT complete.
        n_missing = 0
        partial_before = {}
        for (p, t, a, s) in combos:
            path = self._path(exp_name, a, p, t)
            if not os.path.exists(path):
                n_missing += 1
                continue
            n_rows = self._seed_row_counts(path).get(s, 0)
            if n_rows == 0:
                n_missing += 1
            elif self._expected_rows is not None and n_rows < self._expected_rows:
                partial_before[f"{a}/{p}/tau{t}/seed{s}"] = n_rows
        summary = {
            "experiment_group": cfg.experiment_group,
            "config_hash": cfg.config_hash(),
            "code_version": cfg.code_version,
            "expected": len(combos),
            "completed": n_done,
            "skipped_resume": n_skipped,
            "failed": len(failures),
            "missing": n_missing,
            "partial": len(partial_before),
            "partial_detail": dict(sorted(partial_before.items())),
            "complete": (len(failures) == 0 and n_missing == 0
                         and not partial_before),
            "failures": failures,
        }
        print(f"  [SUMMARY] group={summary['experiment_group']} "
              f"expected={summary['expected']} "
              f"completed={summary['completed']} "
              f"skipped={summary['skipped_resume']} "
              f"failed={summary['failed']} "
              f"missing={summary['missing']} "
              f"partial={summary['partial']} "
              f"complete={summary['complete']}")
        return summary
