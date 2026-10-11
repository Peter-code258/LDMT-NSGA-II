"""Static NSGA-II and DNSGA-II-inspired (response-based) baselines.

v2.2 (P0-2.4): the two baselines are now SEPARATE classes with distinct,
tested behaviour:

- :class:`NSGA2Static`  — static-response benchmark: on environment change
  the previous decision variables are carried over unchanged (pymoo
  re-evaluates them under the new environment); NO immigrants, NO
  prediction, NO migration. Canonical id ``nsga2_static``.

- :class:`DNSGA2Response` — DNSGA-II-*inspired* project implementation:
  keeps a fixed share of the previous population and injects
  ``immigrant_ratio`` random (or heavily-mutated) immigrants. This is NOT a
  reproduction of the original DNSGA-II paper. Canonical id
  ``dnsga2_inspired``.
"""
import numpy as np


class NSGA2Static:
    """Static-response baseline: carry over the population, no immigrants."""
    name = "nsga2_static"

    def __init__(self, pop_size=None):
        self.pop_size = pop_size

    def respond(self, prev_X, prev_F, problem, env_index, rng):
        # Keep the previous decision variables EXACTLY as they were; pymoo
        # re-evaluates them under the new environment. No new individuals.
        new_X = np.asarray(prev_X, dtype=float).copy()
        n_nd = int(np.sum(np.isfinite(prev_F).all(axis=1))) \
            if prev_F is not None else 0
        # v2.4 (K9 fix): real diagnostic instead of hard-coded 0.0 — fraction
        # of individuals sitting ON the variable bounds after the response.
        xl = np.asarray(problem.xl, dtype=float)
        xu = np.asarray(problem.xu, dtype=float)
        bclip = float(np.mean(np.any(
            (new_X <= xl + 1e-9) | (new_X >= xu - 1e-9), axis=1)))
        return new_X, {"lambda_k": 1.0, "boundary_clip_ratio": bclip,
                       "n_nd": n_nd}


class DNSGA2Response:
    """DNSGA-II-inspired response baseline with random/mutant immigrants."""
    name = "dnsga2_inspired"

    def __init__(self, strategy="random", immigrant_ratio=0.2, pop_size=None):
        self.strategy = strategy
        self.immigrant_ratio = float(immigrant_ratio)
        self.pop_size = pop_size

    def _n_immigrants(self, N):
        return int(round(self.immigrant_ratio * N))

    def respond(self, prev_X, prev_F, problem, env_index, rng):
        N = prev_X.shape[0]
        D = prev_X.shape[1]
        xl = np.asarray(problem.xl, dtype=float)
        xu = np.asarray(problem.xu, dtype=float)
        n_imm = self._n_immigrants(N)
        n_keep = N - n_imm
        # keep the first n_keep previous individuals (re-evaluated by pymoo)
        kept = prev_X[:n_keep].copy()
        if self.strategy == "random":
            immigrants = rng.uniform(xl, xu, size=(n_imm, D))
        elif self.strategy == "mutation":
            base = prev_X[:n_imm].copy()
            noise = rng.normal(0.0, 0.1 * (xu - xl), size=base.shape)
            immigrants = np.clip(base + noise, xl, xu)
        else:
            immigrants = rng.uniform(xl, xu, size=(n_imm, D))
        new_X = np.concatenate([kept, immigrants], axis=0)
        new_X = self._resize(new_X, N, xl, xu, rng)
        # v2.4 (K9 fix): real diagnostic instead of hard-coded 0.0.
        bclip = float(np.mean(np.any(
            (new_X <= xl + 1e-9) | (new_X >= xu - 1e-9), axis=1)))
        return new_X, {"lambda_k": 1.0, "boundary_clip_ratio": bclip,
                        "n_nd": int(np.sum(np.isfinite(prev_F).all(axis=1)))}

    def _resize(self, X, N, xl, xu, rng):
        D = X.shape[1]
        if X.shape[0] > N:
            return X[:N]
        while X.shape[0] < N:
            X = np.concatenate(
                [X, rng.uniform(xl, xu, size=(1, D))], axis=0)
        return X
