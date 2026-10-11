"""LDT-NSGA-II: Lightweight Distribution-Transfer dynamic optimizer.

Implements the algorithm described in the opening report (section 5):

5.1 Distribution features:
    - Per-dimension mean mu_e = mean(P_e), per-dimension std sigma_e = std(P_e, ddof=0)
    - Extracted from the previous population (or non-dominated subset)
    - sigma_min = 0.001 * (u - l) to avoid degenerate dimensions

5.2 Trend centre and sampling:
    - If >= 2 environments: mu_hat = clip(mu_e + alpha * (mu_e - mu_{e-1}), l, u)
    - Else: mu_hat = mu_e
    - alpha default 0.5, sensitivity tested {0, 0.25, 0.5, 0.75, 1.0}
    - Transfer candidates: x_tr = clip(mu_hat + z * sigma_e, l, u), z ~ N(0, I)

    New population composition:
        - transfer_ratio rho = 0.5 (migrated candidates)
        - retain_ratio kappa = 0.2 (old population retained)
        - random_ratio = 1 - rho - kappa = 0.3 (uniform random)
    - rho + kappa <= 1 enforced
    - Non-finite migrated values replaced with uniform random samples

v2.1 extensions (optimisation report):
    - adaptive_transfer: scale sigma_e by change magnitude D_k
        sigma_eff = sigma_e * (1 + beta * D_k_norm), beta=adaptive_beta
        large change -> wider exploration; small change -> tight exploitation
    - multi_source: fuse K recent environment features weighted by similarity
        w_i ∝ exp(-||mu_e - mu_i||_2 / (sqrt(D) + eps))
        mu_fused = sum(w_i * mu_i) / sum(w_i);  same for sigma
    - transfer_strategy: pluggable TransferStrategy instance (report 3.1)
        when provided, delegates the migration pipeline to the strategy object

The 50/20/30 split and alpha=0.5 are the default-to-be-frozen configuration
per the opening report (section 7.4).
"""
import numpy as np

from .distribution_features import (
    normalize, denormalize, extract_features, non_dominated_set)
from .environment_manager import ChangeDetector


class LDTNSGA2:
    """LDT-NSGA-II algorithm class (report section 5, 7.4 + v2.1 extensions)."""
    name = "LDT_NSGA2"

    def __init__(self, transfer_ratio=0.5, retain_ratio=0.2,
                 random_ratio=0.3, trend_alpha=0.5, n_probes=10,
                 adaptive_transfer=False, adaptive_beta=1.0,
                 multi_source=False, history_size=3,
                 transfer_strategy=None):
        # enforce rho + kappa <= 1
        self.transfer_ratio = float(transfer_ratio)
        self.retain_ratio = float(retain_ratio)
        self.random_ratio = float(random_ratio)
        assert abs(self.transfer_ratio + self.retain_ratio +
                   self.random_ratio - 1.0) < 1e-9, \
            "transfer + retain + random ratios must sum to 1.0"
        assert self.transfer_ratio + self.retain_ratio <= 1.0 + 1e-9, \
            "transfer + retain must be <= 1.0"
        self.trend_alpha = float(trend_alpha)
        self.n_probes = int(n_probes)
        # v2.1: adaptive transfer strength (report 1.1)
        self.adaptive_transfer = bool(adaptive_transfer)
        self.adaptive_beta = float(adaptive_beta)
        # v2.1: multi-source knowledge fusion (report 1.2)
        self.multi_source = bool(multi_source)
        self.history_size = int(history_size)
        self._feature_history = []  # list of (mu, sigma) in original space
        # v2.1: pluggable transfer strategy (report 3.1)
        self.transfer_strategy = transfer_strategy
        self._probe = None
        self._prev_probe_F = None
        # store previous env features for trend extrapolation
        self._mu_prev = None
        self._sigma_prev = None
        self._dk_max = 1e-12  # running max of D_k for normalisation

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

    # ---- v2.1: multi-source fusion ----
    def _fuse_features(self, mu_e, sigma_e):
        """Fuse current features with K historical features (report 1.2).

        Similarity weights: w_i = exp(-||mu_e - mu_i|| / sqrt(D)).
        Returns fused (mu, sigma).
        """
        if not self.multi_source or not self._feature_history:
            return mu_e, sigma_e
        D = mu_e.shape[0]
        mus = [h[0] for h in self._feature_history] + [mu_e]
        sigs = [h[1] for h in self._feature_history] + [sigma_e]
        weights = []
        for m in mus[:-1]:  # history only
            dist = np.linalg.norm(mu_e - m)
            weights.append(np.exp(-dist / (np.sqrt(D) + 1e-12)))
        weights.append(1.0)  # current environment gets base weight 1
        weights = np.asarray(weights, dtype=float)
        weights /= weights.sum()
        mu_fused = np.sum([w * m for w, m in zip(weights, mus)], axis=0)
        sigma_fused = np.sum([w * s for w, s in zip(weights, sigs)], axis=0)
        return mu_fused, sigma_fused

    # ---- main respond ----
    def respond(self, prev_X, prev_F, problem, env_index, rng):
        N, D = prev_X.shape
        xl = np.asarray(problem.xl, dtype=float)
        xu = np.asarray(problem.xu, dtype=float)

        # v2.1: delegate to pluggable strategy when provided (report 3.1)
        if self.transfer_strategy is not None:
            self._ensure_probe(prev_X, D)
            Dk = self._detect_change(problem)
            new_X, info = self.transfer_strategy.transfer(
                prev_X, prev_F, problem, env_index, rng)
            info["Dk"] = float(Dk)
            return new_X, info

        # 1. extract features from previous population (non-dominated subset)
        nd_mask = non_dominated_set(prev_X, prev_F)
        nd_X = prev_X[nd_mask]
        min_nd = max(10, N // 4)
        if nd_X.shape[0] < min_nd:
            nd_X = prev_X.copy()
        n_nd = int(nd_X.shape[0])

        mu_e, sigma_e = extract_features(nd_X, xl, xu)
        # sigma_min = 0.001 * (u - l) per report 5.1
        sigma_min = 0.001 * (xu - xl)
        sigma_e = np.maximum(sigma_e, sigma_min)

        # 2. probe change magnitude D_k
        self._ensure_probe(prev_X, D)
        Dk = self._detect_change(problem)

        # v2.1: multi-source knowledge fusion (report 1.2)
        mu_eff, sigma_eff = self._fuse_features(mu_e, sigma_e)

        # 3. trend-centre prediction (report 5.2)
        if self._mu_prev is not None:
            mu_hat = np.clip(
                mu_eff + self.trend_alpha * (mu_eff - self._mu_prev),
                xl, xu)
        else:
            mu_hat = mu_eff.copy()

        # v2.1: adaptive transfer strength (report 1.1)
        # large D_k -> widen sampling sigma; small D_k -> tighten
        if self.adaptive_transfer:
            self._dk_max = max(self._dk_max, abs(Dk))
            dk_norm = abs(Dk) / self._dk_max if self._dk_max > 0 else 0.0
            sigma_scale = 1.0 + self.adaptive_beta * dk_norm
            sigma_eff = sigma_eff * sigma_scale
        else:
            sigma_scale = 1.0

        # store for next change (use un-fused current features as history)
        self._mu_prev = mu_e.copy()
        self._sigma_prev = sigma_e.copy()
        # v2.1: maintain history buffer for multi-source fusion
        if self.multi_source:
            self._feature_history.append((mu_e.copy(), sigma_e.copy()))
            if len(self._feature_history) > self.history_size:
                self._feature_history.pop(0)

        # 4. sample transfer candidates: x_tr = clip(mu_hat + z * sigma_e, l, u)
        n_transfer = int(round(self.transfer_ratio * N))
        if n_transfer > 0:
            z = rng.normal(0.0, 1.0, size=(n_transfer, D))
            # broadcast sigma_e (normalized space) -> denormalize scale
            # sample in normalized [0,1] space then denormalize
            mu_hat_norm = normalize(np.atleast_2d(mu_hat), xl, xu)[0]
            sigma_e_norm = sigma_eff / (xu - xl + 1e-12)
            transfer_norm = np.clip(
                mu_hat_norm[None, :] + z * sigma_e_norm[None, :],
                0.0, 1.0)
            transfer_part = denormalize(transfer_norm, xl, xu)
            # replace non-finite with uniform random
            bad = ~np.isfinite(transfer_part).all(axis=1)
            if bad.any():
                transfer_part[bad] = rng.uniform(xl, xu, size=(bad.sum(), D))
        else:
            transfer_part = np.empty((0, D))

        # 5. retain old population
        n_retain = int(round(self.retain_ratio * N))
        if n_retain > 0:
            idx = rng.permutation(N)
            retain_part = prev_X[idx[:n_retain]].copy()
        else:
            retain_part = np.empty((0, D))

        # 6. uniform random immigrants
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

        return new_X, {
            "lambda_k": 1.0,  # not used in LDT per report
            "boundary_clip_ratio": bclip,
            "n_nd": n_nd,
            "Dk": float(Dk),
            "sigma_scale": float(sigma_scale),
        }

    def _resize(self, X, N, xl, xu, rng):
        D = X.shape[1]
        if X.shape[0] > N:
            return X[:N]
        while X.shape[0] < N:
            X = np.concatenate(
                [X, rng.uniform(xl, xu, size=(1, D))], axis=0)
        return X


# Backward-compatible alias for existing experiment scripts
LDMTNSGA2 = LDTNSGA2
