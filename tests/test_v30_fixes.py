"""Regression tests for the v3.0 changes driven by the V2.0 optimisation plan.

Covers:
  plan 2.1  code_version is injected (not hard-coded) and reaches every row
  plan 2.1  run metadata carries git commit, dirty flag, config hash, run_id
  plan 2.2  different experiment groups / sweep values get different run ids
  plan 2.2  a run's unique key includes the experiment group
  plan 2.3  build_report.py resolves paths from its own location
  plan 2.4  run_supplementary defaults to the full project-variant problem set
  plan 5.2  truncated runs are not counted as complete
"""
import io
import json
import os
import subprocess
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)

from experiments.common import (RunConfig, ExperimentRunner, run_single,
                                 make_run_id, git_commit, git_dirty,
                                 ROW_FIELDS)


# --------------------------------------------------- plan 2.1: code_version ----
def test_code_version_is_injected_not_hardcoded():
    """The driver must stamp the version it is TOLD, not a literal.

    Before v3.0 ``DynamicDriver`` hard-coded ``self.code_version = "v2.4"``, so
    every CSV said v2.4 regardless of the code that produced it.
    """
    from algorithms.nsga2_core import DynamicDriver
    from algorithms.dnsga2_response import NSGA2Static
    from experiments.common import problem_factory

    problem = problem_factory("FDA1", tau_t=25, n_t=20)
    algo = NSGA2Static()
    driver = DynamicDriver(algo, problem, "FDA1", n_envs=1, pop_size=20,
                           run_seed=0, code_version="v9.9-test")
    rows = driver.run(n_points=50)
    assert rows, "driver produced no rows"
    assert {r["code_version"] for r in rows} == {"v9.9-test"}


def test_driver_default_version_is_unknown_not_v24():
    """Omitting the version must NOT silently claim to be v2.4."""
    from algorithms.nsga2_core import DynamicDriver
    from algorithms.dnsga2_response import NSGA2Static
    from experiments.common import problem_factory

    problem = problem_factory("FDA1", tau_t=25, n_t=20)
    driver = DynamicDriver(NSGA2Static(), problem, "FDA1", n_envs=1,
                           pop_size=20, run_seed=0)
    assert driver.code_version == "unknown"


def test_run_id_present_in_every_row():
    """plan 5.1: every per-generation row is traceable to one run."""
    assert "run_id" in ROW_FIELDS
    cfg = RunConfig(problem="FDA1", tau_t=25, n_t=20, N=20, n_gen=25,
                    algorithms=["nsga2_static"], experiment_group="t")
    rows = run_single("nsga2_static", "FDA1", cfg, 0)
    ids = {r["run_id"] for r in rows}
    assert len(ids) == 1, f"expected one run_id, got {ids}"
    assert ids.pop(), "run_id must not be empty"


# ---------------------------------------- plan 2.1: run metadata completeness ----
def test_run_meta_has_traceability_fields():
    cfg = RunConfig(problem="FDA1", tau_t=25, n_t=20, N=20, n_gen=25,
                    algorithms=["nsga2_static"], experiment_group="grp-A",
                    alpha=0.75)
    rows, meta = run_single("nsga2_static", "FDA1", cfg, 0, return_meta=True)
    for key in ("run_id", "config_hash", "config", "experiment_group",
                "git_commit", "git_dirty", "code_version", "n_gen_rows",
                "status", "environment"):
        assert key in meta, f"missing metadata field: {key}"
    assert meta["status"] == "ok"
    assert meta["code_version"] == cfg.code_version
    assert meta["n_gen_rows"] == len(rows)
    assert meta["experiment_group"] == "grp-A"
    # the sweep parameter is recorded explicitly, not implied
    assert meta["config"]["alpha"] == 0.75
    assert meta["git_dirty"] in (True, False, None)


def test_config_hash_changes_with_behaviour():
    """plan 2.2: two runs that differ behaviourally must not share an id."""
    base = dict(problem="FDA1", tau_t=25, n_t=20, N=20, n_gen=25,
                algorithms=["ldt_nsga2"], experiment_group="g")
    a = RunConfig(**base)
    b = RunConfig(**{**base, "alpha": 1.0})
    c = RunConfig(**{**base, "experiment_group": "other"})
    assert a.config_hash() != b.config_hash(), "alpha must affect the hash"
    assert a.config_hash() != c.config_hash(), "group must affect the hash"
    d = RunConfig(**base)
    assert a.config_hash() == d.config_hash(), "hash must be deterministic"


def test_run_id_distinguishes_experiment_groups():
    """plan 2.2: the same (algo, problem, seed) in two groups is two runs.

    This is the bug that merged the sensitivity sweep into the main
    comparison during analysis.
    """
    base = dict(problem="FDA1", tau_t=25, n_t=20, N=20, n_gen=25,
                algorithms=["ldt_nsga2"])
    g1 = RunConfig(**base, experiment_group="comparison")
    g2 = RunConfig(**base, experiment_group="sensitivity_alpha", alpha=1.0)
    id1 = make_run_id(g1, "ldt_nsga2", "FDA1", 0)
    id2 = make_run_id(g2, "ldt_nsga2", "FDA1", 0)
    assert id1 != id2
    # and it is stable across calls
    assert id1 == make_run_id(g1, "ldt_nsga2", "FDA1", 0)


def test_run_id_distinguishes_sweep_values():
    base = dict(problem="FDA1", tau_t=25, n_t=20, N=100, n_gen=25,
                algorithms=["ldt_nsga2"],
                experiment_group="sensitivity_alpha")
    ids = {make_run_id(RunConfig(**base, alpha=a), "ldt_nsga2", "FDA1", 0)
           for a in (0.0, 0.5, 1.0)}
    assert len(ids) == 3, f"sweep values collapsed to {ids}"


# ------------------------------------------------ plan 5.2: truncated runs ----
def test_truncated_run_is_not_counted_complete(tmp_path):
    """plan 5.2: a half-written seed must be re-run, not treated as done."""
    out = str(tmp_path)
    runner = ExperimentRunner(out_dir=out)
    runner._expected_rows = 30
    path = os.path.join(out, "comparison_nsga2_static_FDA1_tau25.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        f.write(",".join(ROW_FIELDS) + "\n")
        # seed 0 only got 10 of 30 rows
        for g in range(10):
            row = {k: "" for k in ROW_FIELDS}
            row.update({"seed": 0, "generation": g})
            f.write(",".join(str(row[k]) for k in ROW_FIELDS) + "\n")
    assert 0 not in runner._existing_seeds(path), \
        "a truncated seed must not count as complete"
    assert runner._partial_seeds(path) == {0: 10}


def test_complete_run_is_counted_complete(tmp_path):
    out = str(tmp_path)
    runner = ExperimentRunner(out_dir=out)
    runner._expected_rows = 5
    path = os.path.join(out, "x.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        f.write(",".join(ROW_FIELDS) + "\n")
        for g in range(5):
            row = {k: "" for k in ROW_FIELDS}
            row.update({"seed": 3, "generation": g})
            f.write(",".join(str(row[k]) for k in ROW_FIELDS) + "\n")
    assert runner._existing_seeds(path) == {3}


def test_partial_seed_rows_are_dropped_before_rerun(tmp_path):
    """Re-running a truncated seed must not leave duplicate rows behind."""
    out = str(tmp_path)
    runner = ExperimentRunner(out_dir=out)
    runner._expected_rows = 4
    path = os.path.join(out, "x.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        f.write(",".join(ROW_FIELDS) + "\n")
        for seed, n in ((0, 4), (1, 2)):        # seed 1 is truncated
            for g in range(n):
                row = {k: "" for k in ROW_FIELDS}
                row.update({"seed": seed, "generation": g})
                f.write(",".join(str(row[k]) for k in ROW_FIELDS) + "\n")
    dropped = runner._drop_partial_seeds(path)
    assert dropped == 1
    assert runner._existing_seeds(path) == {0}
    assert runner._seed_row_counts(path) == {0: 4}


# -------------------------------------------- plan 2.4: supplementary default ----
def test_supplementary_default_covers_three_project_variants():
    """plan 2.4: a bare invocation must run DMOP1, DMOP2 AND DMOP3."""
    sys.path.insert(0, PROJ)
    from experiments import run_supplementary as rs
    probs = [p.strip() for p in rs.DEFAULT_PROBLEMS.split(",")]
    assert probs == ["DMOP1", "DMOP2", "DMOP3"], probs
    # project variants must stay distinct from the authoritative definitions
    assert not any(p.endswith("A") for p in probs)


# ------------------------------------------------- plan 2.3: report paths ----
def test_build_report_has_no_hardcoded_drive_path():
    src = io.open(os.path.join(PROJ, "build_report.py"),
                  encoding="utf-8").read()
    head = src[:2500]
    assert "Path(__file__)" in head, "ROOT must derive from the script location"
    assert 'Path(r"D:' not in head, "hard-coded Windows path still present"
    assert 'Path(r"C:' not in head, "hard-coded Windows path still present"


def test_build_report_fails_loudly_on_missing_inputs(tmp_path):
    """plan 2.3: a clean clone must get an actionable error, not a traceback."""
    env = dict(os.environ)
    env["LDT_ROOT"] = str(tmp_path)          # empty tree -> inputs missing
    env["LDT_OUT"] = str(tmp_path / "out.html")
    r = subprocess.run([sys.executable, os.path.join(PROJ, "build_report.py")],
                       capture_output=True, text=True, env=env, cwd=PROJ,
                       encoding="utf-8", errors="replace")
    assert r.returncode != 0
    combined = (r.stdout or "") + (r.stderr or "")
    assert "missing required input file" in combined
