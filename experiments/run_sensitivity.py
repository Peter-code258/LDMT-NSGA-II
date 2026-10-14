"""Sensitivity study per opening report section 7.8.

Studies:
  (a) Population size N in {50, 100, 200}
  (b) Change magnitude n_t in {10, 20, 40}
  (c) Transfer ratio rho in {0.25, 0.5, 0.75}
  (d) Trend coefficient alpha in {0, 0.25, 0.5, 0.75, 1.0}

All on FDA1/FDA3/FDA4, tau=25, 30 seeds each.
"""
import argparse
import os
import sys
import csv

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from algorithms import DynamicDriver, LDTNSGA2
from experiments.common import (RunConfig, ExperimentRunner, problem_factory,
                                 ROW_FIELDS, run_single, git_commit,
                                 git_dirty, environment_meta, make_run_id)

CODE_VERSION = "v3.0"


def run_ldt(problem_name, tau_t, n_t, N, n_gen, seed,
            code_version=CODE_VERSION, **algo_kw):
    problem = problem_factory(problem_name, tau_t=tau_t, n_t=n_t)
    algo = LDTNSGA2(**algo_kw)
    n_envs = max(1, n_gen // tau_t)
    driver = DynamicDriver(algo, problem, problem_name, n_envs=n_envs,
                           pop_size=N, run_seed=seed,
                           code_version=code_version)
    return driver.run(n_points=500)


def _parse_seeds(s):
    if "-" in s:
        lo, hi = s.split("-")
        return list(range(int(lo), int(hi) + 1))
    return [int(x) for x in s.split(",")]


def _run_grid(study, params, problems, seeds, tau, n_gen, out, **algo_kw):
    """One parameter sweep.

    v3.0 (plan 2.2/5.1): the swept variable is recorded explicitly in both the
    per-run metadata and the run identity.  Previously the sensitivity driver
    wrote no metadata at all and stamped rows "unknown", which is how the
    ``sensitivity_*_tau25`` files became indistinguishable from the main
    ``comparison_*_tau25`` files during analysis.
    """
    runner = ExperimentRunner(out_dir=out)
    runner._expected_rows = int(n_gen)
    for p in problems:
        for val in params:
            for s in seeds:
                tag = f"sensitivity_{study}_{p}_tau{tau}_{study}{val}"
                path = os.path.join(out, f"{tag}.csv")
                if os.path.exists(path):
                    done = set()
                    with open(path) as f:
                        for row in csv.DictReader(f):
                            if row.get("seed"):
                                done.add(int(row["seed"]))
                    if s in done:
                        continue
                if study == "N":
                    N, n_t = val, 20
                elif study == "n_t":
                    N, n_t = 100, val
                else:
                    N, n_t = 100, 20
                # the swept variable, so the run is uniquely identifiable
                sweep = {}
                if study == "N":
                    sweep = {"sweep": "N", "N": N}
                elif study == "n_t":
                    sweep = {"sweep": "n_t", "n_t": n_t}
                elif study == "rho":
                    sweep = {"sweep": "rho", "rho": float(val)}
                elif study == "alpha":
                    sweep = {"sweep": "alpha", "alpha": float(val)}
                cfg = RunConfig(
                    problem=p, tau_t=tau, n_t=n_t, N=N, n_gen=n_gen,
                    algorithms=["ldt_nsga2"], seeds=range(30),
                    code_version=CODE_VERSION,
                    experiment_group=f"sensitivity_{study}",
                    alpha=float(val) if study == "alpha" else None,
                    rho=float(val) if study == "rho" else None,
                    extra=sweep)
                try:
                    rows = run_ldt(p, tau, n_t, N, n_gen, s, **algo_kw)
                    write_header = not os.path.exists(path)
                    with open(path, "a", newline="") as f:
                        w = csv.DictWriter(f, fieldnames=ROW_FIELDS)
                        if write_header:
                            w.writeheader()
                        for r in rows:
                            w.writerow({k: r.get(k, "") for k in ROW_FIELDS})
                    runner._append_jsonl("run_meta.jsonl", {
                        "run_id": make_run_id(cfg, "ldt_nsga2", p, s),
                        "config_hash": cfg.config_hash(),
                        "config": cfg.config_dict(),
                        "experiment_group": cfg.experiment_group,
                        "algorithm": "ldt_nsga2",
                        "problem": p, "tau_t": tau, "n_t": n_t, "N": N,
                        "n_gen": n_gen, "seed": int(s),
                        "n_gen_rows": len(rows),
                        "exp": tag,
                        "status": "ok",
                        "code_version": CODE_VERSION,
                        "git_commit": git_commit(),
                        "git_dirty": git_dirty(),
                        "environment": environment_meta(),
                    })
                except Exception as e:
                    print(f"  [FAIL] {tag}/seed{s}: {e}")
                    runner._append_jsonl("failed_runs.jsonl", {
                        "exp": tag, "algorithm": "ldt_nsga2", "problem": p,
                        "tau_t": tau, "seed": int(s), "error": str(e),
                        "experiment_group": cfg.experiment_group,
                        "code_version": CODE_VERSION,
                        "git_commit": git_commit()})


def main():
    ap = argparse.ArgumentParser(description="LDT-NSGA-II sensitivity")
    ap.add_argument("--problems", default="FDA1,FDA3,FDA4")
    ap.add_argument("--tau", type=int, default=25)
    ap.add_argument("--n-gen", type=int, default=300)
    ap.add_argument("--seeds", default="0-29")
    ap.add_argument("--out", default="results/raw")
    ap.add_argument("--limit-seeds", type=int, default=None)
    args = ap.parse_args()

    problems = [p.strip() for p in args.problems.split(",")]
    seeds = _parse_seeds(args.seeds)
    if args.limit_seeds is not None:
        seeds = seeds[:args.limit_seeds]

    print("Sensitivity (a): N in {50,100,200}")
    _run_grid("N", [50, 100, 200], problems, seeds, args.tau,
              args.n_gen, args.out)
    print("Sensitivity (b): n_t in {10,20,40}")
    _run_grid("n_t", [10, 20, 40], problems, seeds, args.tau,
              args.n_gen, args.out)
    print("Sensitivity (c): rho in {0.25,0.5,0.75}")
    for rho in [0.25, 0.5, 0.75]:
        _run_grid("rho", [rho], problems, seeds, args.tau,
                  args.n_gen, args.out,
                  transfer_ratio=rho, retain_ratio=0.2,
                  random_ratio=1.0-rho-0.2)
    print("Sensitivity (d): alpha in {0,0.25,0.5,0.75,1.0}")
    for alpha in [0, 0.25, 0.5, 0.75, 1.0]:
        _run_grid("alpha", [alpha], problems, seeds, args.tau,
                  args.n_gen, args.out, trend_alpha=alpha)
    print("Sensitivity done. CSVs in", args.out)


if __name__ == "__main__":
    main()
