"""Run an experiment from a YAML definition (v3.2, plan 5.3).

YAML is the source of truth; CLI flags override it.  Every override is echoed
and recorded, so a run can always be traced back to (file + overrides).

    python -m experiments.run_from_config configs/main_experiment.yaml
    python -m experiments.run_from_config configs/main_experiment.yaml --taus 25
    python -m experiments.run_from_config --list
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from experiments.common import RunConfig, ExperimentRunner, CANONICAL_ALGOS  # noqa: E402
from experiments.config import (load_config, describe, available_configs,   # noqa: E402
                                ConfigError)


def main():
    ap = argparse.ArgumentParser(description="config-driven experiment runner")
    ap.add_argument("config", nargs="?", help="path to a configs/*.yaml file")
    ap.add_argument("--list", action="store_true", help="list available configs")
    # overrides -- None means "use the value from the file"
    ap.add_argument("--problems")
    ap.add_argument("--taus")
    ap.add_argument("--seeds")
    ap.add_argument("--algos")
    ap.add_argument("--N", type=int)
    ap.add_argument("--n-gen", type=int)
    ap.add_argument("--n-t", type=int)
    ap.add_argument("--out")
    ap.add_argument("--limit-seeds", type=int)
    ap.add_argument("--dry-run", action="store_true",
                    help="resolve and print the plan without running")
    args = ap.parse_args()

    if args.list or not args.config:
        print("available configs:")
        for c in available_configs():
            print("  configs/" + c)
        return 0

    overrides = {
        "problems": args.problems,
        "taus": args.taus,
        "seeds": args.seeds,
        "algorithms": args.algos,
        "N": args.N,
        "n_gen": args.n_gen,
        "n_t": args.n_t,
        "out": args.out,
        "limit_seeds": args.limit_seeds,
    }
    try:
        cfg, applied = load_config(args.config, overrides)
    except ConfigError as e:
        print(f"[config] ERROR: {e}", file=sys.stderr)
        return 2

    print(f"[config] file      : {os.path.relpath(cfg['_path'])}")
    print(f"[config] name      : {cfg.get('name')}")
    print(f"[config] effective : {describe(cfg)}")
    if applied:
        print("[config] overrides applied:")
        for k, (old, new) in sorted(applied.items()):
            print(f"           {k}: {old!r} -> {new!r}")
    else:
        print("[config] overrides applied: (none)")

    problems = cfg.get("problems") or []
    algos = cfg.get("algorithms") or list(CANONICAL_ALGOS)
    taus = cfg.get("taus") or [cfg.get("tau") or 25]
    seeds = cfg.get("seeds") or list(range(30))
    if cfg.get("limit_seeds"):
        seeds = seeds[:cfg["limit_seeds"]]
    if not problems:
        print("[config] ERROR: no problems defined", file=sys.stderr)
        return 2

    n_expected = len(problems) * len(taus) * len(algos) * len(seeds)
    print(f"[config] plan      : {len(problems)} problems x {len(taus)} taus "
          f"x {len(algos)} algos x {len(seeds)} seeds = {n_expected} runs")
    if args.dry_run:
        print("[config] dry run -- nothing executed")
        return 0

    rc = RunConfig(problem=problems[0], tau_t=taus[0],
                   n_t=cfg.get("n_t", 20), N=cfg.get("N", 100),
                   n_gen=cfg.get("n_gen", 300), algorithms=algos,
                   seeds=range(30), code_version=cfg.get("code_version", "v3.2"),
                   experiment_group=cfg.get("experiment_group", "config"),
                   algo_kwargs=cfg.get("algo_kwargs") or {},
                   extra={"config_file": os.path.basename(cfg["_path"]),
                          "overrides": {k: list(v) for k, v in applied.items()}})
    runner = ExperimentRunner(out_dir=cfg.get("out", "results/raw"))
    summary = runner.run(cfg.get("experiment_group", "config"), rc,
                         algorithms=algos, problems=problems, taus=taus,
                         seeds=seeds, progress_label=f"[{cfg.get('name')}]")
    print("[config] summary:", {k: summary[k] for k in
                                ("expected", "completed", "skipped_resume",
                                 "failed", "missing", "partial", "complete")})
    return 0 if summary.get("complete") else 1


if __name__ == "__main__":
    sys.exit(main())
