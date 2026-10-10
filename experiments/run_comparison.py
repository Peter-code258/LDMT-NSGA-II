"""Main comparison experiment: FDA1/FDA3/FDA4 x tau{10,25,50} x 4 algos x 30 seeds.

CLI example:
    python -m experiments.run_comparison --problems FDA1,FDA3,FDA4 \
        --taus 10,25,50 --seeds 0-29 --out results/raw

Quick smoke:
    python -m experiments.run_comparison --problems FDA3 --taus 25 \
        --seeds 0 --out results/raw --n-gen 50 --limit-seeds 0
"""
import argparse
import os
import sys

# ensure project root importable when run as a script
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from experiments.common import (RunConfig, ExperimentRunner,
                                 ALGORITHM_REGISTRY)


def _parse_seeds(s):
    if "-" in s:
        lo, hi = s.split("-")
        return list(range(int(lo), int(hi) + 1))
    return [int(x) for x in s.split(",")]


def main():
    ap = argparse.ArgumentParser(description="LDMT-NSGA-II main comparison")
    ap.add_argument("--problems", default="FDA1,FDA3,FDA4")
    ap.add_argument("--taus", default="10,25,50")
    ap.add_argument("--seeds", default="0-29")
    ap.add_argument("--out", default="results/raw")
    ap.add_argument("--n-gen", type=int, default=300)
    ap.add_argument("--n", type=int, default=100, help="population size")
    ap.add_argument("--n-t", type=int, default=20)
    ap.add_argument("--algos", default=",".join(ALGORITHM_REGISTRY.keys()))
    ap.add_argument("--limit-seeds", type=int, default=None,
                    help="take only first N seeds (quick tests)")
    args = ap.parse_args()

    problems = [p.strip() for p in args.problems.split(",") if p.strip()]
    taus = [int(t) for t in args.taus.split(",") if t.strip()]
    seeds = _parse_seeds(args.seeds)
    algos = [a.strip() for a in args.algos.split(",") if a.strip()]

    cfg = RunConfig(problem=problems[0] if problems else "FDA1",
                    tau_t=taus[0] if taus else 25, n_t=args.n_t, N=args.n,
                    n_gen=args.n_gen, algorithms=algos, seeds=range(30),
                    code_version="v1.0")

    runner = ExperimentRunner(out_dir=args.out)
    try:
        from tqdm import tqdm
        have_tqdm = True
    except Exception:
        have_tqdm = False

    shown_seeds = (seeds[:args.limit_seeds] if (args.limit_seeds and args.limit_seeds > 0)
                   else seeds)
    print(f"Comparison: problems={problems} taus={taus} "
          f"algos={algos} seeds={shown_seeds}")
    runner.run("comparison", cfg, algorithms=algos, problems=problems,
               taus=taus, seeds=seeds, limit_seeds=args.limit_seeds,
               progress_label="[comparison]")
    print("Done. CSVs in", args.out)


if __name__ == "__main__":
    main()
