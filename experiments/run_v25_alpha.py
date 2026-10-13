"""v2.5 experiment driver: alpha=1.0 (full-step) LDT vs the five v2.4 baselines.

WHY A SEPARATE DRIVER
---------------------
``experiments/run_comparison.py`` always writes files named
``comparison_<algo>_<problem>_tau<tau>.csv`` into ``--out``.  ``analyze.py``
loads ``results/raw/comparison_*.csv`` with a glob, so writing the new runs
into the default output directory would silently contaminate the frozen v2.4
summary.  This driver is therefore byte-compatible in output format but is
always invoked with a SEPARATE ``--out`` directory (``results/raw_v25``).
The v2.4 artefacts under ``results/raw`` and ``results/summary`` are never
read or written by this script.

WHAT IT RUNS (per problem x tau x seed)
---------------------------------------
exactly one fresh DynamicDriver run per algorithm:
    nsga2_static, dnsga2_inspired, ar1_predictive, gmm_inspired, pps,
    ldt_alpha1                       <- LDT-NSGA-II with trend_alpha = 1.0

All six share the identical NSGA-II kernel, operators, probe set, budget and
seed, so the only difference from the v2.4 comparison is that the
``ldt_nsga2`` entry (alpha=0.5) is replaced by ``ldt_alpha1`` (alpha=1.0).
The alpha=0.5 reference rows already exist in ``results/raw`` and are reused
for the paired comparison, so no run is duplicated.

EXAMPLES
--------
    python -m experiments.run_v25_alpha --problems FDA1,FDA3,FDA4 \
        --taus 10,25,50 --out results/raw_v25
    python -m experiments.run_v25_alpha --problems DMOP1,DMOP2,DMOP3 \
        --out results/raw_v25 --phase supplementary
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from experiments.common import (RunConfig, ExperimentRunner, CANONICAL_ALGOS)

# v2.4 comparison run order, with the LDT entry swapped for the alpha=1.0
# variant.  Keeping the five baselines means the new directory is a complete,
# self-contained comparison table (no cross-directory dependency for the
# baselines themselves).
V25_ALGOS = ["nsga2_static", "dnsga2_inspired", "ar1_predictive",
             "gmm_inspired", "pps", "ldt_alpha1"]


def _parse_seeds(s):
    if "-" in s:
        lo, hi = s.split("-")
        return list(range(int(lo), int(hi) + 1))
    return [int(x) for x in s.split(",")]


def main():
    ap = argparse.ArgumentParser(
        description="v2.5 LDT alpha=1.0 comparison / supplementary runs")
    ap.add_argument("--problems", default="FDA1,FDA3,FDA4")
    ap.add_argument("--taus", default="10,25,50")
    ap.add_argument("--seeds", default="0-29")
    ap.add_argument("--out", default="results/raw_v25")
    ap.add_argument("--n-gen", type=int, default=300)
    ap.add_argument("--n", type=int, default=100, help="population size")
    ap.add_argument("--n-t", type=int, default=20)
    ap.add_argument("--phase", default="comparison",
                    choices=["comparison", "supplementary"],
                    help="CSV file prefix; also tags run_meta.jsonl")
    ap.add_argument("--algos", default=",".join(V25_ALGOS))
    ap.add_argument("--limit-seeds", type=int, default=None)
    args = ap.parse_args()

    problems = [p.strip() for p in args.problems.split(",") if p.strip()]
    taus = [int(t) for t in args.taus.split(",") if t.strip()]
    seeds = _parse_seeds(args.seeds)
    algos = [a.strip() for a in args.algos.split(",") if a.strip()]

    # Sanity guard: refuse to write a v2.5 run into the frozen v2.4 directory.
    out_norm = os.path.normpath(args.out).replace("\\", "/").rstrip("/")
    if out_norm.endswith("results/raw") or out_norm.endswith("results/raw_legacy_v21") \
            or "raw_legacy" in out_norm:
        sys.exit("REFUSED: --out must not point at a frozen v2.4/v2.1 result "
                 "directory. Use results/raw_v25.")

    cfg = RunConfig(problem=problems[0] if problems else "FDA1",
                    tau_t=taus[0] if taus else 25, n_t=args.n_t, N=args.n,
                    n_gen=args.n_gen, algorithms=algos, seeds=range(30),
                    code_version="v2.5")

    shown = (seeds[:args.limit_seeds]
             if (args.limit_seeds and args.limit_seeds > 0) else seeds)
    print(f"[v2.5/{args.phase}] problems={problems} taus={taus} "
          f"algos={algos} seeds={shown} out={args.out}")

    runner = ExperimentRunner(out_dir=args.out)
    summary = runner.run(args.phase, cfg, algorithms=algos,
                         problems=problems, taus=taus, seeds=seeds,
                         limit_seeds=args.limit_seeds,
                         progress_label=f"[v2.5/{args.phase}]")
    print("Done. CSVs in", args.out)
    if summary.get("failed", 0) > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
