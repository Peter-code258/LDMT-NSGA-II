"""Supplementary experiment: dMOP family, tau=25, n_t=20, N=100, n_gen=300,
all 6 canonical algorithms, 30 seeds.

v3.0 (plan 2.4): the default problem list is now the full PROJECT-VARIANT set
``DMOP1,DMOP2,DMOP3``.  Previously it defaulted to ``DMOP1,DMOP2`` while the
documented full command passed three problems, so a bare invocation (including
the one issued by ``main.py``) silently produced 360 rows instead of 540 while
appearing to have completed.  The project variants (``DMOP1/2/3``) stay
strictly separate from the authoritative definitions (``DMOP1A/2A/3A``), which
belong to the v2.5+/v3.0 alpha study.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from experiments.common import (RunConfig, ExperimentRunner,
                                 CANONICAL_ALGOS)

DEFAULT_PROBLEMS = "DMOP1,DMOP2,DMOP3"


def main():
    ap = argparse.ArgumentParser(description="LDT-NSGA-II supplementary dMOP")
    ap.add_argument("--problems", default=DEFAULT_PROBLEMS)
    ap.add_argument("--tau", type=int, default=25)
    ap.add_argument("--n-gen", type=int, default=300)
    ap.add_argument("--seeds", default="0-29")
    ap.add_argument("--out", default="results/raw")
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--n-t", type=int, default=20)
    # v2.2: default to the canonical list ONLY (registry keys include
    # aliases that would otherwise run the same algorithm repeatedly)
    ap.add_argument("--algos", default=",".join(CANONICAL_ALGOS))
    ap.add_argument("--limit-seeds", type=int, default=None)
    args = ap.parse_args()

    problems = [p.strip() for p in args.problems.split(",") if p.strip()]
    seeds = ([int(s) for s in args.seeds.split(",")]
             if "-" not in args.seeds
             else list(range(int(args.seeds.split("-")[0]),
                              int(args.seeds.split("-")[1]) + 1)))
    if args.limit_seeds is not None and args.limit_seeds > 0:
        seeds = seeds[:args.limit_seeds]
    algos = [a.strip() for a in args.algos.split(",")]

    cfg = RunConfig(problem=problems[0], tau_t=args.tau, n_t=args.n_t,
                    N=args.n, n_gen=args.n_gen, algorithms=algos,
                    seeds=range(30),
                    code_version="v3.0", experiment_group="supplementary")
    runner = ExperimentRunner(out_dir=args.out)
    # v3.0 (plan 2.4): print the expected combination count up front so a
    # mismatch between intent and execution is visible immediately.
    n_expected = len(problems) * len(algos) * len(seeds)
    print(f"Supplementary: problems={problems} tau={args.tau} "
          f"algos={algos} seeds={seeds}")
    print(f"  expected combos = {len(problems)} problems x {len(algos)} algos "
          f"x {len(seeds)} seeds = {n_expected}")
    summary = runner.run("supplementary", cfg, algorithms=algos,
                         problems=problems, taus=[args.tau], seeds=seeds,
                         limit_seeds=args.limit_seeds,
                         progress_label="[supplementary]")
    print(f"Supplementary done. CSVs in {args.out}")
    # v3.0 (plan 5.2): a non-zero exit signals incomplete data, so a wrapping
    # pipeline cannot mistake a partial sweep for a finished one.
    return 0 if summary.get("complete") else 1


if __name__ == "__main__":
    sys.exit(main())
