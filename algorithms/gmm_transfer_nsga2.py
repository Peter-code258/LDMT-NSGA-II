"""GMM full-distribution transfer NSGA-II baseline.

On environment change, fit a Gaussian Mixture on the normalized previous
population (or its non-dominated subset), sample N new individuals from it,
denormalize, clip to bounds, and add ~10% random immigrants for diversity.

NOTE: this is a PROJECT COMPARISON IMPLEMENTATION, not a reproduction of any
specific paper's KTMM. It stands in for the "full-distribution transfer"
family of methods (learn the whole distribution, re-sample).
"""
import numpy as np

from .distribution_features import (normalize, denormalize, non_dominated_set)
from .environment_manager import ChangeDetector


class GMMTransfer:
    name = "gmm_inspired"

    def __init__(self, n_components=None, reg_covar=1e-4,
                 immigrant_ratio=0.1, max_components=3):
        self.n_components = n_components
        self.reg_covar = float(reg_covar)
        self.immigrant_ratio = float(immigrant_ratio)
        self.max_components = int(max_components)
        self._probe = None
        self._prev_probe_F = None

    def _ensure_probe(self, prev_X, problem):
        shared = getattr(problem, "probe", None)
        if shared is not None:
            self._probe = np.asarray(shared, dtype=float)
            return self._probe
        if self._probe is None or self._probe.shape[1] != prev_X.shape[1]:
            n = min(10, prev_X.shape[0])
            self._probe = np.asarray(prev_X[:n], dtype=float).copy()
        return self._probe

    def _detect_change(self, problem):
        probe = self._probe
        if probe is None:
            return 0.0
        new_F = problem._calc_objectives(np.atleast_2d(probe))
        prev_F = self._prev_probe_F
        if prev_F is None:
            self._prev_probe_F = new_F
            return 0.0
        cd = ChangeDetector(probe, problem.n_obj)
        Dk = cd.detect(prev_F, new_F)
        self._prev_probe_F = new_F
        return Dk

    def respond(self, prev_X, prev_F, problem, env_index, rng):
        from sklearn.mixture import GaussianMixture
        N, D = prev_X.shape
        xl = np.asarray(problem.xl, dtype=float)
        xu = np.asarray(problem.xu, dtype=float)
        nd_mask = non_dominated_set(prev_X, prev_F)
        nd_X = prev_X[nd_mask] if nd_mask.sum() >= max(3, N // 10) else prev_X
        Z = normalize(nd_X, xl, xu)
        Z = np.clip(Z, 0.0, 1.0)
        k = self.n_components or min(self.max_components, max(1, N // 10))
        k = max(1, min(k, len(nd_X) - 1))
        self._ensure_probe(prev_X, problem)
        Dk = self._detect_change(problem)
        # P2-4.1: reproducible per-(run, env) GMM seed derived from the
        # run's rng stream instead of a hard-coded random_state=42
        gmm_seed = int(rng.integers(0, 2 ** 31 - 1))
        samples = None
        used_cov = None
        last_err = None
        fallback_used = False
        for cov_type in ("diag", "full", "spherical"):
            try:
                gmm = GaussianMixture(
                    n_components=k, covariance_type=cov_type,
                    reg_covar=self.reg_covar, random_state=gmm_seed)
                gmm.fit(Z)
                samples = gmm.sample(N)[0]
                used_cov = cov_type
                break
            except Exception as e:
                last_err = e
                continue
        if samples is None:
            # fallback: uniform random sampling (recorded, never silent)
            fallback_used = True
            used_cov = f"fallback_uniform (fit error: {last_err})"
            samples = rng.uniform(0, 1, size=(N, D))
        samples = np.clip(samples, 0.0, 1.0)
        new_X = denormalize(samples, xl, xu)
        # immigrants (K1 fix): reserve room BEFORE appending so random
        # immigrants survive _resize. Previously "sample N rows + append
        # immigrants + truncate to first N" dropped every immigrant, making
        # immigrant_ratio a dead parameter (pure-transfer behaviour).
        n_imm = int(round(self.immigrant_ratio * N))
        immigrants = rng.uniform(xl, xu, size=(n_imm, D))
        new_X = np.concatenate([new_X[:N - n_imm], immigrants], axis=0)
        new_X = self._resize(new_X, N, xl, xu, rng)
        bclip = float(np.mean(np.any(
            (new_X[:, :D] <= xl + 1e-9) | (new_X[:, :D] >= xu - 1e-9), axis=1)))
        return new_X, {"lambda_k": 1.0, "boundary_clip_ratio": bclip,
                        "n_nd": int(nd_mask.sum()), "Dk": Dk,
                        "gmm_cov_type": used_cov,
                        "gmm_fallback": bool(fallback_used),
                        "gmm_seed": gmm_seed}

    def _resize(self, X, N, xl, xu, rng):
        D = X.shape[1]
        if X.shape[0] > N:
            return X[:N]
        while X.shape[0] < N:
            X = np.concatenate(
                [X, rng.uniform(xl, xu, size=(1, D))], axis=0)
        return X
