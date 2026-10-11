"""Shared NSGA-II kernel and dynamic driver.

All four algorithms share the same pymoo NSGA-II operators. Each algorithm
class implements ``respond(prev_X, prev_F, problem, env_index, rng)`` producing
the warm-start decision variables for the next environment (optionally
returning an ``(X, info)`` tuple). The :class:`DynamicDriver` runs the
per-environment segment loop generation by generation so IGD/HV can be
recorded EACH generation (needed for recovery curves).

For fair comparison a FRESH NSGA2 instance is created each environment segment
(internal state is clean), seeded with ``run_seed + env_index`` so each
segment is reproducible yet distinct.
"""
import time
import numpy as np
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.operators.sampling.rnd import FloatRandomSampling

from evaluation.igd import igd
from evaluation.hv import hv, get_reference_point
from problems.true_fronts import get_true_pf


def make_nsga2(pop_size, sampling=None):
    """Return a fresh NSGA2 with the fixed operators (SBX/PM, eta=15/20)."""
    if sampling is None:
        sampling = FloatRandomSampling()
    return NSGA2(
        pop_size=pop_size,
        sampling=sampling,
        crossover=SBX(prob=0.9, eta=15),
        mutation=PM(eta=20),
        eliminate_duplicates=True,
    )


def _row_metrics(F, problem_name, problem, n_points):
    """Compute IGD/HV for the current population objectives."""
    M = F.shape[1]
    t = problem.t
    if problem_name.upper() in ("FDA3", "DMOP2"):
        ref = get_true_pf(problem_name, M=M, n_points=n_points, t=t)
    else:
        ref = get_true_pf(problem_name, M=M, n_points=n_points)
    ref_pt = get_reference_point(problem_name, M=M)
    return igd(F, ref), hv(F, ref_pt)


class DynamicDriver:
    """Runs one algorithm over all environments of a problem.

    Parameters
    ----------
    algorithm : object with ``respond(prev_X, prev_F, problem, env_index, rng)``
    problem : DynamicProblem (already constructed with tau_t, n_t)
    problem_name : str (for IGD/HV front lookup)
    n_envs : int, number of environments (= n_gen // tau_t)
    pop_size : int (N)
    run_seed : int (base seed; segment k uses run_seed + k)
    """

    def __init__(self, algorithm, problem, problem_name, n_envs,
                 pop_size, run_seed):
        self.algorithm = algorithm
        self.problem = problem
        self.problem_name = problem_name
        self.n_envs = int(n_envs)
        self.pop_size = int(pop_size)
        self.run_seed = int(run_seed)
        self.tau_t = problem.tau_t
        self.n_t = problem.n_t
        self.n_var = problem.n_var
        self.M = problem.n_obj
        self.xl = np.asarray(problem.xl, dtype=float)
        self.xu = np.asarray(problem.xu, dtype=float)
        self.rows = []

    def run(self, n_points=500):
        rng = np.random.default_rng(self.run_seed)
        prev_X = None
        prev_F = None
        for k in range(self.n_envs):
            g_start = k * self.tau_t
            self.problem.set_generation(g_start)
            # ---- warm-start decision variables ----
            t0_mig = time.perf_counter()
            if k == 0:
                warm_X = None  # let NSGA2 use random sampling
                migration_runtime = 0.0
                lam_k = 1.0
                bclip = 0.0
                n_nd = 0
            else:
                result = self.algorithm.respond(prev_X, prev_F,
                                                 self.problem, k, rng)
                if isinstance(result, tuple):
                    warm_X, info = result
                    lam_k = float(info.get("lambda_k", 1.0))
                    bclip = float(info.get("boundary_clip_ratio", 0.0))
                    n_nd = int(info.get("n_nd", 0))
                else:
                    warm_X = result
                    lam_k = 1.0
                    bclip = 0.0
                    n_nd = 0
                migration_runtime = time.perf_counter() - t0_mig
            # ---- run tau_t generations generation-by-generation ----
            algo = make_nsga2(
                self.pop_size,
                sampling=warm_X if warm_X is not None else None)
            algo.setup(self.problem, seed=self.run_seed + k,
                       termination=("n_gen", self.tau_t))
            for gen_in_env in range(self.tau_t):
                g_global = g_start + gen_in_env
                self.problem.set_generation(g_global)
                t_seg0 = time.perf_counter()
                algo.next()
                gen_runtime = time.perf_counter() - t_seg0
                pop = algo.pop
                F = pop.get("F")
                igd_val, hv_val = _row_metrics(
                    F, self.problem_name, self.problem, n_points)
                # Migration time is a per-environment quantity; record it only
                # at the first generation of a (changed) environment so that
                # summing the column over generations yields the TRUE total
                # migration time (report 5.4.3).
                mig_recorded = migration_runtime if (
                    gen_in_env == 0 and k > 0) else 0.0
                # cumulative function evaluations reported by pymoo
                n_eval = int(getattr(algo, "n_eval", (g_global + 1)
                                    * self.pop_size))
                self.rows.append({
                    "algorithm": getattr(self.algorithm, "name",
                                         self.algorithm.__class__.__name__),
                    "problem": self.problem_name,
                    "seed": self.run_seed,
                    "population_size": self.pop_size,
                    "total_generations": self.n_envs * self.tau_t,
                    "n_t": self.n_t,
                    "tau_t": self.tau_t,
                    "environment_index": k,
                    "generation": g_global,
                    "IGD": igd_val,
                    "HV": hv_val,
                    "generation_runtime_s": gen_runtime,
                    "migration_runtime_s": mig_recorded,
                    "n_evaluations": n_eval,
                    "code_version": "v1.0",
                    "lambda_k": lam_k,
                    "boundary_clip_ratio": bclip,
                    "n_nd": n_nd,
                })
            pop = algo.pop
            prev_X = pop.get("X")
            prev_F = pop.get("F")
        return self.rows
