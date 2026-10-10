"""Supplementary experiment: dMOP1/dMOP2, tau=25, n_t=20, N=100, n_gen=300,
all 4 algorithms, 30 seeds."""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from experiments.common import (RunConfig, ExperimentRunner,
                                 ALGORITHM_REGISTRY)


def main():
    ap = argparse.ArgumentParser(description="LDMT-NSGA-II supplementary dMOP")
    ap.add_argument("--problems", default="DMOP1,DMOP2")
    ap.add_argument("--tau", type=int, default=25)
    ap.add_argument("--n-gen", type=int, default=300)
    ap.add_argument("--seeds", default="0-29")
    ap.add_argument("--out", default="results/raw")
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--n-t", type=int, default=20)
    ap.add_argument("--algos", default=",".join(ALGORITHM_REGISTRY.keys()))
    ap.add_argument("--limit-seeds", type=int, default=None)
    args = ap.parse_args()

    problems = [p.strip() for p in args.problems.split(",")]
    seeds = ([int(s) for s in args.seeds.split(",")]
             if "-" not in args.seeds
             else list(range(int(args.seeds.split("-")[0]),
                              int(args.seeds.split("-")[1]) + 1)))
    if args.limit_seeds is not None and args.limit_seeds > 0:
        seeds = seeds[:args.limit_seeds]
    algos = [a.strip() for a in args.algos.split(",")]

    cfg = RunConfig(problem=problems[0], tau_t=args.tau, n_t=args.n_t,
                    N=args.n, n_gen=args.n_gen, algorithms=algos,
                    seeds=range(30))
    runner = ExperimentRunner(out_dir=args.out)
    print(f"Supplementary: problems={problems} tau={args.tau} "
          f"algos={algos} seeds={seeds}")
    runner.run("supplementary", cfg, algorithms=algos, problems=problems,
               taus=[args.tau], seeds=seeds,
               limit_seeds=args.limit_seeds,
               progress_label="[supplementary]")
    print("Supplementary done. CSVs in", args.out)


if __name__ == "__main__":
    main()
