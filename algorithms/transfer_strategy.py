"""Pluggable transfer strategies for dynamic multi-objective optimisation.

The TransferStrategy interface decouples the "how to migrate" logic from the
algorithm main loop, making it easy to swap between simplified moment transfer,
GMM-based transfer, kernel subspace alignment, etc.
"""
from abc import ABC, abstractmethod
import numpy as np

from .distribution_features import (
    normalize, denormalize, extract_features, non_dominated_set)


class TransferStrategy(ABC):
    """Abstract base class for population transfer strategies."""

    @abstractmethod
    def transfer(self, prev_X, prev_F, problem, env_index, rng, **kw):
        """Return the migrated population for the new environment.

        Parameters
        ----------
        prev_X : (N, D) array  – previous population decision variables
        prev_F : (N, M) array  – previous population objectives
        problem : DynamicProblem
        env_index : int
        rng : numpy.random.Generator
        kw : dict  – strategy-specific parameters

        Returns
        -------
        new_X : (N, D) array
        info : dict  – diagnostic information
        """
        raise NotImplementedError


class MomentTransferStrategy(TransferStrategy):
    """Simplified mean+std moment transfer (LDT-NSGA-II default).

    Uses per-dimension mean mu_e and std sigma_e from the non-dominated subset,
    applies a linear trend centre, and samples Gaussian perturbations.
    """

    def __init__(self, transfer_ratio=0.5, retain_ratio=0.2,
                 random_ratio=0.3, trend_alpha=0.5):
        self.transfer_ratio = float(transfer_ratio)
        self.retain_ratio = float(retain_ratio)
        self.random_ratio = float(random_ratio)
        self.trend_alpha = float(trend_alpha)
        assert abs(self.transfer_ratio + self.retain_ratio +
                   self.random_ratio - 1.0) < 1e-9
        assert self.transfer_ratio + self.retain_ratio <= 1.0 + 1e-9
        self._mu_prev = None
        self._sigma_prev = None

    def transfer(self, prev_X, prev_F, problem, env_index, rng, **kw):
        N, D = prev_X.shape
        xl = np.asarray(problem.xl, dtype=float)
        xu = np.asarray(problem.xu, dtype=float)
        # v2.4 (K4 fix): extract_features returns features in NORMALIZED
        # [0,1]^D space; the previous version mixed coordinate systems
        # (real-scale sigma_min floor, clip to real bounds, then a second
        # normalization). Everything below now stays in normalized space
        # until the single denormalize at the end — mirroring ldmt_nsga2.
        # 1. extract features from non-dominated subset (normalized space)
        nd_mask = non_dominated_set(prev_X, prev_F)
        nd_X = prev_X[nd_mask] if nd_mask.sum() > 0 else prev_X
        mu_e, sigma_e = extract_features(nd_X, xl, xu)   # already in [0,1]
        sigma_e = np.maximum(sigma_e, 1e-3)              # normalized floor
        # 2. trend-centre prediction (normalized space)
        if self._mu_prev is not None:
            mu_hat = np.clip(mu_e + self.trend_alpha * (mu_e - self._mu_prev),
                             0.0, 1.0)
        else:
            mu_hat = mu_e.copy()
        self._mu_prev = mu_e.copy()
        self._sigma_prev = sigma_e.copy()
        # 3. sample transfer candidates in [0,1]^D, denormalize once
        n_transfer = int(round(self.transfer_ratio * N))
        if n_transfer > 0:
            z = rng.normal(0.0, 1.0, size=(n_transfer, D))
            transfer_norm = np.clip(
                mu_hat[None, :] + z * sigma_e[None, :], 0.0, 1.0)
            transfer_part = denormalize(transfer_norm, xl, xu)
            bad = ~np.isfinite(transfer_part).all(axis=1)
            if bad.any():
                transfer_part[bad] = rng.uniform(xl, xu, size=(bad.sum(), D))
        else:
            transfer_part = np.empty((0, D))
        # 4. retain old population
        n_retain = int(round(self.retain_ratio * N))
        if n_retain > 0:
            idx = rng.permutation(N)
            retain_part = prev_X[idx[:n_retain]].copy()
        else:
            retain_part = np.empty((0, D))
        # 5. uniform random immigrants
        n_random = int(round(self.random_ratio * N))
        if n_random > 0:
            random_part = rng.uniform(xl, xu, size=(n_random, D))
        else:
            random_part = np.empty((0, D))
        new_X = np.concatenate([transfer_part, retain_part, random_part],
                               axis=0)
        new_X = self._resize(new_X, N, xl, xu, rng)
        bclip = float(np.mean(np.any(
            (new_X <= xl + 1e-9) | (new_X >= xu - 1e-9), axis=1)))
        return new_X, {"lambda_k": 1.0, "boundary_clip_ratio": bclip,
                       "n_nd": int(nd_mask.sum())}

    def _resize(self, X, N, xl, xu, rng):
        D = X.shape[1]
        if X.shape[0] > N:
            return X[:N]
        while X.shape[0] < N:
            X = np.concatenate(
                [X, rng.uniform(xl, xu, size=(1, D))], axis=0)
        return X


class GMMTransferStrategy(TransferStrategy):
    """Placeholder for full GMM-based distribution transfer.

    Not yet implemented; raises NotImplementedError.
    """

    def transfer(self, prev_X, prev_F, problem, env_index, rng, **kw):
        raise NotImplementedError(
            "GMMTransferStrategy is not yet implemented. Use MomentTransferStrategy.")


# Registry for easy lookup
TRANSFER_REGISTRY = {
    "moment": MomentTransferStrategy,
    "gmm": GMMTransferStrategy,
}
