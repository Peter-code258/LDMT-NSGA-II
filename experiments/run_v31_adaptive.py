# -*- coding: utf-8 -*-
"""v3.1 experiment driver: online adaptive trend step (plan 3.1).

Two phases, deliberately separated so the rule's threshold is chosen on a
CALIBRATION set and only validated afterwards -- plan 3.1 requires that
"threshold choice and the final test set stay separate".

    calibration : grid of (alpha_min, alpha_max, gamma) candidates evaluated on
                  three problems bracketing the two regimes (one translating
                  set, one static set).  No held-out problem is touched.
    comparison  : the selected candidate plus the fixed-alpha variants
                  (alpha = 0 / 0.5 / 1.0) and the reference baselines, run on
                  the DISJOINT held-out problems, paired over identical seeds.

Design notes
------------
* Output goes to ``results/raw_v31`` only.  Writing into ``results/raw`` or
  ``results/raw_v25`` (the frozen v2.4/v2.5 evidence) is refused.
* Every run uses the same budget (N=100, n_gen=300, N*tau_t evaluations per
  environment), so all arms are budget-fair.
* Runs are parallelised across the whole (candidate, problem, tau, seed)
  space; each worker is a fresh process.

Usage:
    python -m experiments.run_v31_adaptive --phase calibration
    python -m experiments.run_v31_adaptive --phase comparison --candidate 0.5,1.0
    python -m experiments.run_v31_adaptive --phase both
"""
import argparse
import csv
import itertools
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from experiments.common import RunConfig, ROW_FIELDS, run_single  # noqa: E402

CODE_VERSION = "v3.1"
DEFAULT_OUT = "results/raw_v31"
FROZEN_DIRS = ("results/raw", "results/raw_v25")

# ── registered candidate grid for the adaptive rule ────────────────────────
# alpha_max is pinned at 1.0: that endpoint is already supported by the v2.5
# evidence (full step wins on 11 of 12 translating configurations).
CANDIDATES = [
    # (name, alpha_min, alpha_max, gamma)
    ("a0.50_g1.0", 0.50, 1.0, 1.0),   # the initially registered rule
    ("a0.25_g1.0", 0.25, 1.0, 1.0),
    ("a0.00_g1.0", 0.00, 1.0, 1.0),
    ("a0.50_g0.5", 0.50, 1.0, 0.5),
    ("a0.25_g0.5", 0.25, 1.0, 0.5),
    ("a0.75_g1.0", 0.75, 1.0, 1.0),   # deliberately too conservative
    ("a0.00_g0.5", 0.00, 1.0, 0.5),
]

CALIBRATION_PROBLEMS = ["FDA1", "DMOP1A", "DMOP3A"]
CALIBRATION_TAUS = [10, 25]
CALIBRATION_SEEDS = list(range(10))

HELDOUT_PROBLEMS = ["FDA3", "FDA4", "DMOP2A"]
HELDOUT_TAUS = [10, 25, 50]

# endpoint validation: the Type III problem (static PS -> step must shrink) and
# one translating problem (step must stay full).  Plan 3.1 asks explicitly for
# both endpoints to be checked, and neither was in the held-out set above.
ENDPOINT_PROBLEMS = ["DMOP1A", "FDA1"]
ENDPOINT_TAUS = [10, 25, 50]

# v3.2 (plan 3.2/3.3): mechanism comparison for the §7 roadmap features.  Every
# arm differs from the plain LDT only in the named mechanism, so the comparison
# is attributable; all arms share budget N*tau_t per environment.
V32_PROBLEMS = ["FDA1", "FDA3", "FDA4"]
V32_TAUS = [10, 25, 50]
V32_ARMS = ["ldt_half", "ldt_full", "ldt_adaptive", "ldt_cov", "ldt_lowrank",
            "ldt_gated", "ldt_v32"]

COMPARISON_ALGOS = ["ldt_adaptive", "ldt_alpha0", "ldt_half", "ldt_full",
                    "pps", "nsga2_static"]


def _guard(out):
    norm = os.path.normpath(out).replace("\\", "/")
    for frozen in FROZEN_DIRS:
        if norm == frozen or norm.startswith(frozen + "/"):
            raise SystemExit(
                f"refusing to write into the frozen directory {frozen!r} "
                f"(got {out!r}); v3.1 data must go to results/raw_v31")
    return out


def _install_variants():
    """Register the v3.1 arms in the shared algorithm registry (once/worker).

    Keys are unique and deliberately kept out of CANONICAL_ALGOS so no existing
    experiment or default list changes.  All variant parameters travel through
    ``RunConfig.algo_kwargs`` -- no ad-hoc subclasses.
    """
    from algorithms import LDTNSGA2, LDTNSGA2Alpha1, LDTNSGA2Adaptive
    from experiments import common as C

    C.ALGORITHM_REGISTRY.setdefault("ldt_alpha0", LDTNSGA2)
    C.ALGORITHM_REGISTRY.setdefault("ldt_half", LDTNSGA2)
    C.ALGORITHM_REGISTRY.setdefault("ldt_full", LDTNSGA2Alpha1)
    C.ALGORITHM_REGISTRY.setdefault("ldt_adaptive", LDTNSGA2Adaptive)
    return C


# arm -> (registry key, canonical algorithm id, extra constructor kwargs)
ARM_SPEC = {
    "ldt_adaptive": ("ldt_adaptive", "ldt_adaptive", {}),
    "ldt_half": ("ldt_half", "ldt_nsga2", {"trend_alpha": 0.5}),
    "ldt_full": ("ldt_full", "ldt_alpha1", {"trend_alpha": 1.0}),
    "ldt_alpha0": ("ldt_alpha0", "ldt_nsga2", {"trend_alpha": 0.0}),
    "pps": ("pps", "pps", {}),
    "nsga2_static": ("nsga2_static", "nsga2_static", {}),
    # v3.2 roadmap mechanisms (plan 3.2 / 3.3)
    "ldt_cov": ("ldt_cov", "ldt_cov", {}),
    "ldt_lowrank": ("ldt_lowrank", "ldt_lowrank", {"low_rank": 3}),
    "ldt_gated": ("ldt_gated", "ldt_gated", {}),
    "ldt_v32": ("ldt_v32", "ldt_v32", {}),
}


def _run_one(job):
    """Run one (arm, problem, tau, seed).  Executed inside a worker."""
    arm, problem, tau, seed, phase, amin, gam = job
    C = _install_variants()
    if arm.startswith("adapt_"):
        cname = arm[len("adapt_"):]
        cand = dict((c[0], c) for c in CANDIDATES)[cname]
        key, canon = "ldt_adaptive", "ldt_adaptive"
        kw = dict(alpha_min=cand[1], alpha_max=cand[2], gamma=cand[3])
    else:
        key, canon, base = ARM_SPEC[arm]
        kw = dict(base)
        if arm == "ldt_adaptive":
            kw = dict(alpha_min=amin, alpha_max=1.0, gamma=gam)

    cfg = RunConfig(problem=problem, tau_t=tau, n_t=20, N=100, n_gen=300,
                    algorithms=[canon], seeds=range(30),
                    code_version=CODE_VERSION,
                    experiment_group=f"v31_{phase}",
                    algo_kwargs=kw,
                    extra={"v31_arm": arm})
    rows, meta = run_single(key, problem, cfg, seed, return_meta=True)
    for r in rows:
        r["algorithm"] = canon
    meta["algorithm"] = canon
    meta["v31_arm"] = arm
    return arm, rows, meta


def _seed_path(out, phase, alias, problem, tau):
    return os.path.join(out, f"v31_{phase}_{alias}_{problem}_tau{tau}.csv")


def _done_seeds(path):
    done = set()
    if os.path.exists(path):
        with open(path, newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                if row.get("seed"):
                    done.add(int(row["seed"]))
    return done


def main():
    ap = argparse.ArgumentParser(description="v3.1 adaptive-alpha study")
    ap.add_argument("--phase", default="both",
                    choices=["calibration", "comparison", "endpoint", "v32",
                             "both", "all"])
    ap.add_argument("--out", default=DEFAULT_OUT)
    ap.add_argument("--candidate", default="0.5,1.0",
                    help="alpha_min,gamma selected from the calibration phase")
    ap.add_argument("--workers", type=int,
                    default=min(12, max(1, (os.cpu_count() or 4) - 4)),
                    help="parallel worker processes (a large pool can be killed "
                         "by the OS; 8-12 is stable on this machine)")
    ap.add_argument("--seeds", default=None, help="e.g. 0-29")
    args = ap.parse_args()

    amin, gam = (float(x) for x in args.candidate.split(","))
    out = _guard(args.out)
    os.makedirs(out, exist_ok=True)

    seeds_override = None
    if args.seeds:
        s = args.seeds
        if "-" in s:
            lo, hi = s.split("-")
            seeds_override = list(range(int(lo), int(hi) + 1))
        else:
            seeds_override = [int(x) for x in s.split(",")]

    phases = {"calibration": ["calibration"],
              "comparison": ["comparison"],
              "endpoint": ["endpoint"],
              "v32": ["v32"],
              "both": ["calibration", "comparison"],
              "all": ["calibration", "comparison", "endpoint", "v32"]}[args.phase]
    print(f"[v3.1] code_version={CODE_VERSION} out={out} "
          f"workers={args.workers} selected=(alpha_min={amin}, gamma={gam})")

    rc = 0
    for phase in phases:
        if phase == "calibration":
            probs, taus = CALIBRATION_PROBLEMS, CALIBRATION_TAUS
            these = seeds_override if seeds_override is not None \
                else CALIBRATION_SEEDS
            arms = [f"adapt_{c[0]}" for c in CANDIDATES]
        elif phase == "endpoint":
            probs, taus = ENDPOINT_PROBLEMS, ENDPOINT_TAUS
            these = seeds_override if seeds_override is not None \
                else list(range(30))
            arms = COMPARISON_ALGOS
        elif phase == "v32":
            probs, taus = V32_PROBLEMS, V32_TAUS
            these = seeds_override if seeds_override is not None \
                else list(range(30))
            arms = V32_ARMS
        else:
            probs, taus = HELDOUT_PROBLEMS, HELDOUT_TAUS
            these = seeds_override if seeds_override is not None \
                else list(range(30))
            arms = COMPARISON_ALGOS

        jobs = [(a, p, t, s, phase, amin, gam)
                for a, p, t, s in itertools.product(arms, probs, taus, these)]
        todo = [j for j in jobs
                if j[3] not in _done_seeds(_seed_path(out, phase, j[0], j[1], j[2]))]
        print(f"\n[v3.1/{phase}] problems={probs} taus={taus} arms={len(arms)} "
              f"seeds={len(these)} -> {len(jobs)} runs ({len(todo)} to do)")
        if not todo:
            print(f"[v3.1/{phase}] nothing to do")
            continue

        t0 = time.perf_counter()
        n_done = n_fail = 0
        failures = []
        # An over-large pool can be killed by the OS ("terminated abruptly"),
        # which would lose the whole batch.  Retry at a decreasing pool size,
        # re-deriving the outstanding work from disk each time.
        pool_size = max(1, args.workers)
        attempted = set()
        while True:
            remaining = []
            for j in todo:
                key = (j[0], j[1], j[2], j[3])
                if key in attempted:
                    continue
                if j[3] in _done_seeds(
                        _seed_path(out, phase, j[0], j[1], j[2])):
                    continue
                remaining.append(j)
            if not remaining:
                break
            print(f"  [pool] running {len(remaining)} job(s) with "
                  f"{pool_size} worker(s)")
            chunk_ok = 0
            try:
                with ProcessPoolExecutor(max_workers=pool_size) as ex:
                    futs = {ex.submit(_run_one, j): j for j in remaining}
                    for fut in as_completed(futs):
                        j = futs[fut]
                        key = (j[0], j[1], j[2], j[3])
                        attempted.add(key)
                        try:
                            alias, rows, meta = fut.result()
                        except Exception as e:              # noqa: BLE001
                            n_fail += 1
                            failures.append({"arm": j[0], "problem": j[1],
                                             "tau": j[2], "seed": j[3],
                                             "error": str(e)})
                            continue
                        path = _seed_path(out, phase, alias, j[1], j[2])
                        header = not os.path.exists(path)
                        with open(path, "a", newline="",
                                  encoding="utf-8") as fh:
                            w = csv.DictWriter(fh, fieldnames=ROW_FIELDS)
                            if header:
                                w.writeheader()
                            for r in rows:
                                w.writerow({k: r.get(k, "") for k in ROW_FIELDS})
                        meta["exp"] = f"v31_{phase}_{alias}"
                        with open(os.path.join(out, "run_meta.jsonl"), "a",
                                  newline="", encoding="utf-8") as fh:
                            fh.write(json.dumps(meta, ensure_ascii=False,
                                                default=str) + "\n")
                        n_done += 1
                        chunk_ok += 1
                        if n_done % 50 == 0:
                            el = time.perf_counter() - t0
                            print(f"  {n_done} done ({el:.0f}s, "
                                  f"{n_done/max(el,1e-9):.1f} runs/s)")
            except Exception as e:                          # noqa: BLE001
                print(f"  [pool] pool failed at size {pool_size}: "
                      f"{type(e).__name__}: {e}")
            if chunk_ok == 0 and pool_size == 1:
                print("  [pool] single-worker pass made no progress; aborting")
                break
            if chunk_ok == 0 or pool_size > 2:
                pool_size = max(1, pool_size // 2)
            # loop re-derives the remaining work from disk
        el = time.perf_counter() - t0
        print(f"[v3.1/{phase}] completed={n_done} failed={n_fail} "
              f"in {el:.0f}s")
        if failures:
            with open(os.path.join(out, "failed_runs.jsonl"), "a",
                      newline="", encoding="utf-8") as fh:
                for f in failures:
                    fh.write(json.dumps(f, ensure_ascii=False, default=str)
                             + "\n")
            rc = 1
    print("\n[v3.1] done. CSVs in", out)
    return rc


if __name__ == "__main__":
    sys.exit(main())
