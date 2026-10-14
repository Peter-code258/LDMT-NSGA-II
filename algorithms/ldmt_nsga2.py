"""LDT-NSGA-II: Lightweight Distribution-Transfer dynamic optimizer.

v2.2 (source-review fixes):
- P0-0B: ALL distribution features, trend centre, sigma, sampling and
  clipping are computed in NORMALIZED [0,1]^D space; denormalize is applied
  exactly once when emitting decision variables. Previous versions mixed
  normalized means with raw bounds xl/xu (coordinate-system bug).
- P1-1E: change probes come from the shared, fixed probe set attached to the
  problem (real bounds [xl, xu], fixed seed) instead of clipping prev_X rows.
- P2: adaptive sigma scale is capped to [0.5, 2.0]; multi-source fusion uses
  mixture second-moment variance (not a plain sigma average) plus a
  similarity weight floor.

Core method (opening report section 5):
    mu_e, sigma_e    per-dimension mean/std of the non-dominated subset,
                     in normalized [0,1]^D space
    sigma_min        1e-3 in normalized units (approx. 0.001 * span)
    trend centre     mu_hat = clip(mu_e + alpha * (mu_e - mu_prev), 0, 1)
    sampling         x_norm = clip(mu_hat + z * sigma_eff, 0, 1), z ~ N(0, I)
    composition      rho=0.5 transfer / kappa=0.2 retained / 0.3 random
"""
import numpy as np

from .distribution_features import (
    normalize, denormalize, extract_features, non_dominated_set)
from .environment_manager import ChangeDetector


class LDTNSGA2:
    """LDT-NSGA-II algorithm class (v2.2, normalized-space kernel)."""
    name = "ldt_nsga2"

    # caps for the adaptive sigma multiplier (P2-4.3)
    SIGMA_SCALE_MIN = 0.5
    SIGMA_SCALE_MAX = 2.0
    # minimum fused weight for a history source (P2-4.2)
    FUSION_WEIGHT_FLOOR = 0.05

    def __init__(self, transfer_ratio=0.5, retain_ratio=0.2,
                 random_ratio=0.3, trend_alpha=0.5, n_probes=10,
                 adaptive_transfer=False, adaptive_beta=1.0,
                 multi_source=False, history_size=3,
                 transfer_strategy=None):
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
        self.adaptive_transfer = bool(adaptive_transfer)
        self.adaptive_beta = float(adaptive_beta)
        self.multi_source = bool(multi_source)
        self.history_size = int(history_size)
        # normalized-space feature history [(mu, sigma), ...]
        self._feature_history = []
        self.transfer_strategy = transfer_strategy
        self._probe = None
        self._prev_probe_F = None
        self._mu_prev = None    # normalized mean of the last environment
        self._sigma_prev = None  # normalized std of the last environment
        self._dk_sum = 0.0      # running sum of |D_k| for normalization
        self._dk_n = 0

    # ---- shared fixed probe (P1-1E) ----
    def _ensure_probe(self, prev_X, problem):
        """Use the problem-owned fixed probe (real bounds, fixed seed)."""
        shared = getattr(problem, "probe", None)
        if shared is not None:
            self._probe = np.asarray(shared, dtype=float)
            return self._probe
        # fallback for problems without a probe attribute
        if self._probe is None or self._probe.shape[1] != prev_X.shape[1]:
            n = min(self.n_probes, prev_X.shape[0])
            self._probe = np.asarray(prev_X[:n], dtype=float).copy()
        return self._probe

    def establish_baseline(self, problem):
        """Record the probe's objective values in the CURRENT environment.

        v3.0 (plan 2.6): the change magnitude D_k is a difference between two
        environments' probe readings, so it is only defined once a reading
        exists.  Previously the first call to :meth:`_detect_change` happened
        during the first *response* -- i.e. after the environment had already
        changed -- so it had nothing to compare against and returned 0.0,
        disguising "unknown" as "no change" and corrupting the adaptive-sigma
        diagnostic at the first transition.

        The driver now calls this at the start of every environment (including
        environment 0), so a real D_k is available from the first change
        onwards.  The probe is problem-owned and fixed, so this does not
        consume the decision variables and is safe to call before the search.
        """
        self._ensure_probe(np.empty((0, problem.n_var)), problem)
        if self._probe is None:
            return None
        self._prev_probe_F = problem._calc_objectives(
            np.atleast_2d(self._probe))
        return None

    def _detect_change(self, problem):
        """Relative probe change vs. the previous environment, or None.

        Returning ``None`` means "not measurable" (no baseline yet) and is
        deliberately distinct from ``0.0`` meaning "measured, no change".
        """
        probe = self._probe
        if probe is None:
            return None
        new_F = problem._calc_objectives(np.atleast_2d(probe))
        prev_F = self._prev_probe_F
        if prev_F is None:
            # no baseline: record one, but report "unavailable", not "zero"
            self._prev_probe_F = new_F
            return None
        cd = ChangeDetector(probe, problem.n_obj)
        Dk = cd.detect(prev_F, new_F)
        self._prev_probe_F = new_F
        return Dk

    # ---- P2-4.2: multi-source fusion with mixture second-moment variance ----
    def _fuse_features(self, mu_e, sigma_e):
        """Fuse current normalized features with K historical ones.

        Weights: w_i = exp(-||mu_e - mu_i|| / sqrt(D)); sources with
        w < FUSION_WEIGHT_FLOOR are dropped (except the current env, weight
        fixed at 1 before renormalization). Variance is fused via the mixture
        second moment: var_fused = sum_i w_i * (sigma_i^2 + (mu_i - mu_f)^2),
        sigma_fused = sqrt(var_fused). Returns (mu, sigma) in [0,1] units.
        """
        if not self.multi_source or not self._feature_history:
            return mu_e, sigma_e
        D = mu_e.shape[0]
        mus = [h[0] for h in self._feature_history]
        sigs = [h[1] for h in self._feature_history]
        weights = []
        for m in mus:
            w = float(np.exp(-np.linalg.norm(mu_e - m) / (np.sqrt(D) + 1e-12)))
            weights.append(w)
        kept = [(w, m, s) for w, m, s in zip(weights, mus, sigs)
                if w >= self.FUSION_WEIGHT_FLOOR]
        if not kept:  # all history too dissimilar -> current only
            return mu_e, sigma_e
        w_all = np.array([k[0] for k in kept] + [1.0])
        w_all = w_all / w_all.sum()
        mu_all = np.array([k[1] for k in kept] + [mu_e])
        sg_all = np.array([k[2] for k in kept] + [sigma_e])
        mu_f = np.sum(w_all[:, None] * mu_all, axis=0)
        var_f = np.sum(w_all[:, None] *
                       (sg_all ** 2 + (mu_all - mu_f[None, :]) ** 2), axis=0)
        return mu_f, np.sqrt(np.maximum(var_f, 1e-12))

    # ---- main respond: everything in normalized [0,1]^D ----
    def respond(self, prev_X, prev_F, problem, env_index, rng):
        N, D = prev_X.shape
        xl = np.asarray(problem.xl, dtype=float)
        xu = np.asarray(problem.xu, dtype=float)

        # v2.1: delegate to a pluggable strategy when provided
        if self.transfer_strategy is not None:
            self._ensure_probe(prev_X, problem)
            Dk = self._detect_change(problem)
            new_X, info = self.transfer_strategy.transfer(
                prev_X, prev_F, problem, env_index, rng)
            info["Dk"] = float(Dk)
            return new_X, info

        # 1. features of the non-dominated subset (normalized space)
        nd_mask = non_dominated_set(prev_X, prev_F)
        nd_X = prev_X[nd_mask]
        if nd_X.shape[0] < max(10, N // 4):
            nd_X = prev_X
        n_nd = int(nd_X.shape[0])
        mu_e, sigma_e = extract_features(nd_X, xl, xu)   # already in [0,1]
        sigma_e = np.maximum(sigma_e, 1e-3)              # normalized floor

        # 2. probe change magnitude D_k (shared fixed probe, real bounds)
        self._ensure_probe(prev_X, problem)
        Dk = self._detect_change(problem)

        # 3. multi-source fusion (normalized space)
        mu_eff, sigma_eff = self._fuse_features(mu_e, sigma_e)

        # 4. trend-centre prediction, clipped in [0,1]
        if self._mu_prev is not None:
            mu_hat = np.clip(
                mu_eff + self.trend_alpha * (mu_eff - self._mu_prev),
                0.0, 1.0)
        else:
            mu_hat = np.clip(mu_eff, 0.0, 1.0)

        # 5. adaptive sigma scale (capped, normalized by running mean of |D_k|)
        sigma_scale = 1.0
        if self.adaptive_transfer and Dk is not None:
            self._dk_sum += abs(Dk)
            self._dk_n += 1
            dk_ref = self._dk_sum / self._dk_n
            dk_norm = abs(Dk) / dk_ref if dk_ref > 1e-12 else 0.0
            sigma_scale = float(np.clip(
                1.0 + self.adaptive_beta * dk_norm,
                self.SIGMA_SCALE_MIN, self.SIGMA_SCALE_MAX))
            sigma_eff = sigma_eff * sigma_scale

        # 6. store current features (normalized) for the next change
        self._mu_prev = mu_e.copy()
        self._sigma_prev = sigma_e.copy()
        if self.multi_source:
            self._feature_history.append((mu_e.copy(), sigma_e.copy()))
            if len(self._feature_history) > self.history_size:
                self._feature_history.pop(0)

        # 7. sample transfer candidates in [0,1]^D, denormalize once
        n_transfer = int(round(self.transfer_ratio * N))
        if n_transfer > 0:
            noise = rng.normal(0.0, 1.0, size=(n_transfer, D))
            transfer_norm = np.clip(mu_hat[None, :] + noise * sigma_eff[None, :],
                                    0.0, 1.0)
            transfer_part = denormalize(transfer_norm, xl, xu)
            bad = ~np.isfinite(transfer_part).all(axis=1)
            if bad.any():
                transfer_part[bad] = rng.uniform(xl, xu, size=(bad.sum(), D))
        else:
            transfer_part = np.empty((0, D))

        # 8. retain a slice of the old population
        n_retain = int(round(self.retain_ratio * N))
        if n_retain > 0:
            idx = rng.permutation(N)
            retain_part = prev_X[idx[:n_retain]].copy()
        else:
            retain_part = np.empty((0, D))

        # 9. uniform random immigrants
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
            "lambda_k": 1.0,
            "boundary_clip_ratio": bclip,
            "n_nd": n_nd,
            # v3.0 (plan 2.6): None = "not measurable yet"; never a fake 0.0
            "Dk": (None if Dk is None else float(Dk)),
            "sigma_scale": sigma_scale,
        }

    def _resize(self, X, N, xl, xu, rng):
        D = X.shape[1]
        if X.shape[0] > N:
            return X[:N]
        while X.shape[0] < N:
            X = np.concatenate(
                [X, rng.uniform(xl, xu, size=(1, D))], axis=0)
        return X


class LDTNSGA2Alpha1(LDTNSGA2):
    """LDT-NSGA-II with the FULL-STEP trend extrapolation (alpha = 1.0).

    Motivation (v2.5 experiment): the v2.4 sensitivity study found that the
    default half-step coefficient alpha=0.5 is systematically outperformed by
    the full-step first-order difference alpha=1.0 on all three FDA problems
    (see ``sensitivity_alpha_*`` CSVs). This subclass pins alpha=1.0 with a
    distinct canonical id so the effect of the extrapolation STEP SIZE can be
    separated from the effect of the moment-transfer PARADIGM itself, without
    altering any existing v2.4 artefact.

    Behaviour is otherwise byte-for-byte identical to :class:`LDTNSGA2`: same
    feature extraction, same sigma floor, same 50/20/30 composition, same
    single denormalization. Only ``trend_alpha`` differs.
    """
    name = "ldt_alpha1"

    def __init__(self, trend_alpha=1.0, **kwargs):
        super().__init__(trend_alpha=trend_alpha, **kwargs)


# Backward-compatible alias for existing experiment scripts
LDMTNSGA2 = LDTNSGA2
