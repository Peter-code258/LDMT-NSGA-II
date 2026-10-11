"""DNSGA-II-A (response-based) baseline.

On environment change, re-evaluate the previous population (pymoo handles
re-evaluation under the new env) and inject ``immigrant_ratio`` of random
(resp. heavily-mutated) immigrants.
"""
import numpy as np


class DNSGA2Response:
    name = "DNSGA2"

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
        return new_X, {"lambda_k": 1.0, "boundary_clip_ratio": 0.0,
                        "n_nd": int(np.sum(np.isfinite(prev_F).all(axis=1)))}

    def _resize(self, X, N, xl, xu, rng):
        D = X.shape[1]
        if X.shape[0] > N:
            return X[:N]
        while X.shape[0] < N:
            X = np.concatenate(
                [X, rng.uniform(xl, xu, size=(1, D))], axis=0)
        return X
