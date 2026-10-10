"""LDMT-NSGA-II: Lightweight Distribution-Migration dynamic optimizer.

On each environment change the algorithm:
  1. selects the non-dominated subset ND_X of the previous population
     (falls back to the whole population if |ND_X| is too small);
  2. extracts mean-std features (mu_prev, sigma_prev) from ND_X;
  3. computes the probe-based change magnitude D_k (via ChangeDetector);
  4. predicts the next distribution (mu_hat, sigma_hat) with the
     :class:`MigrationFunction` (AR(1)-style translation with adaptively
     clipped lambda_k);
  5. transfers the WHOLE previous population to the predicted distribution;
  6. composes the new population with the 70/20/10 split:
        - 70% transferred individuals,
        - 20% local-mutated individuals (small Gaussian noise on a random
          subset of dims),
        - 10% uniform random immigrants.

The 70/20/10 split and alpha=1.0 are the default-to-be-frozen configuration.
"""
import numpy as np

from .distribution_features import (
    normalize, denormalize, extract_features, non_dominated_set,
    transfer_population, MigrationFunction)
from .environment_manager import ChangeDetector


class LDMTNSGA2:
    name = "LDMT_NSGA2"

    def __init__(self, transfer_ratio=0.7, local_mut_ratio=0.2,
                 immigrant_ratio=0.1, alpha=1.0, n_probes=10):
        self.transfer_ratio = float(transfer_ratio)
        self.local_mut_ratio = float(local_mut_ratio)
        self.immigrant_ratio = float(immigrant_ratio)
        self.alpha = float(alpha)
        self.n_probes = int(n_probes)
        self.migration = MigrationFunction(alpha=alpha)
        self._probe = None
        self._prev_probe_F = None
        # store previous env features for the AR(1)-style extrapolation
        self._mu_prev2 = None
        self._sigma_prev2 = None
        self._D_prev = None

    # ---- probe / change detection ----
    def _ensure_probe(self, prev_X, D):
        if self._probe is None or self._probe.shape[1] != D:
            n = min(self.n_probes, prev_X.shape[0])
            self._probe = np.clip(prev_X[:n], 0, 1).copy()
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

    # ---- main respond ----
    def respond(self, prev_X, prev_F, problem, env_index, rng):
        N, D = prev_X.shape
        xl = np.asarray(problem.xl, dtype=float)
        xu = np.asarray(problem.xu, dtype=float)

        # 1. non-dominated subset
        nd_mask = non_dominated_set(prev_X, prev_F)
        nd_X = prev_X[nd_mask]
        min_nd = max(10, N // 4)
        if nd_X.shape[0] < min_nd:
            nd_X = prev_X.copy()
        n_nd = int(nd_X.shape[0])

        # 2. features from ND_X
        mu_prev, sigma_prev = extract_features(nd_X, xl, xu)

        # 3. probe change magnitude D_k
        self._ensure_probe(prev_X, D)
        Dk = self._detect_change(problem)

        # 4. predict next distribution
        mu_prev2 = self._mu_prev2 if self._mu_prev2 is not None else mu_prev
        sigma_prev2 = (self._sigma_prev2 if self._sigma_prev2 is not None
                        else sigma_prev)
        D_prev = self._D_prev if self._D_prev is not None else 0.0
        mu_hat, sigma_hat, lam_k = self.migration.predict(
            mu_prev, mu_prev2, sigma_prev, sigma_prev2, Dk, D_prev)
        # record into the migration history (for the lambda threshold logic)
        self.migration.record(mu_prev, sigma_prev, Dk)
        # update the prev2 slots for the NEXT change
        self._mu_prev2 = mu_prev
        self._sigma_prev2 = sigma_prev
        self._D_prev = Dk

        # 5. transfer the WHOLE previous population to the predicted dist
        transferred = transfer_population(
            prev_X, xl, xu, mu_prev, sigma_prev, mu_hat, sigma_hat)
        bclip = float(np.mean(np.any(
            (transferred <= xl + 1e-9) | (transferred >= xu - 1e-9),
            axis=1)))

        # 6/7/8. compose with 70/20/10 split
        n_transfer = int(round(self.transfer_ratio * N))
        n_mut = int(round(self.local_mut_ratio * N))
        n_imm = int(round(self.immigrant_ratio * N))
        # transfer part
        idx = rng.permutation(N)
        transfer_part = transferred[idx[:n_transfer]]
        # local mutation part: take from transferred, add small noise
        mut_part = transferred[idx[n_transfer:n_transfer + n_mut]].copy()
        if mut_part.shape[0] > 0:
            n_dims = max(1, D // 3)
            dim_sel = rng.choice(D, size=n_dims, replace=False)
            noise = rng.normal(0.0, 0.05, size=mut_part.shape)
            mut_part[:, dim_sel] += noise[:, dim_sel]
            mut_part = np.clip(mut_part, xl, xu)
        # immigrants
        imm_part = rng.uniform(xl, xu, size=(n_imm, D))
        new_X = np.concatenate([transfer_part, mut_part, imm_part], axis=0)
        new_X = self._resize(new_X, N, xl, xu, rng)

        return new_X, {
            "lambda_k": float(lam_k),
            "boundary_clip_ratio": bclip,
            "n_nd": n_nd,
            "Dk": float(Dk),
        }

    def _resize(self, X, N, xl, xu, rng):
        D = X.shape[1]
        if X.shape[0] > N:
            return X[:N]
        while X.shape[0] < N:
            X = np.concatenate(
                [X, rng.uniform(xl, xu, size=(1, D))], axis=0)
        return X
