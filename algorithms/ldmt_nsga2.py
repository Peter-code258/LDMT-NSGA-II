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
                 transfer_strategy=None, adaptive_alpha=False,
                 alpha_min=0.5, alpha_max=1.0, gamma=1.0,
                 noise_samples=200, conservative_first_change=True,
                 cov_mode="diag", low_rank=None, rank_energy=0.95,
                 shrinkage="auto", gating=False, gate_eta=1.0,
                 gate_rho_min=0.1, gate_rho_max=0.5,
                 detect_change=False, detect_threshold=None,
                 detect_window=3):
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
        # ---- v3.1 (plan 3.1): online trend-step adaptation ----
        # alpha_k = clip(alpha_min + (alpha_max-alpha_min) * s_k, lo, hi)
        # where s_k is a dimensionless "movement dominates noise?" score built
        # from the DECISION-space centre displacement and a bootstrap estimate
        # of how uncertain that centre is.  Plan 3.1 requires this separation
        # explicitly: the probe change magnitude D_k measures OBJECTIVE-space
        # change and must not be used to infer decision-space displacement.
        self.adaptive_alpha = bool(adaptive_alpha)
        self.alpha_min = float(alpha_min)
        self.alpha_max = float(alpha_max)
        assert 0.0 <= self.alpha_min <= self.alpha_max <= 1.0, \
            "require 0 <= alpha_min <= alpha_max <= 1"
        self.gamma = float(gamma)
        self.noise_samples = int(noise_samples)
        # With no previous prediction residual (the very first change) there is
        # no evidence that movement is sustained, so use the conservative end.
        self.conservative_first_change = bool(conservative_first_change)
        # ---- v3.2 (plan 3.2): correlated / low-rank distribution modelling ----
        # "diag" (default, the v2.x behaviour) samples each dimension
        # independently from sigma_e.  "full" samples from a shrinkage-regularised
        # full covariance; "lowrank" samples inside the leading PCA subspace and
        # adds a small isotropic residual so the sampler cannot collapse.
        # The plan warns explicitly not to assume the full covariance wins: the
        # sample N_nd is small relative to D, so the estimate must be shrunk.
        self.cov_mode = str(cov_mode).lower()
        assert self.cov_mode in ("diag", "full", "lowrank"), self.cov_mode
        self.low_rank = None if low_rank in (None, 0, "0") else int(low_rank)
        self.rank_energy = float(rank_energy)
        assert 0.0 < self.rank_energy <= 1.0, self.rank_energy
        self.shrinkage = shrinkage
        # ---- v3.2 (plan 3.3): negative-transfer gating ----
        # rho_k adapts to the recent prediction residual.  NOTE the plan's
        # constraint: any candidate-validation budget must be charged to the
        # total evaluation budget.  This implementation therefore adds NO extra
        # objective evaluations -- it re-allocates the EXISTING population
        # budget between transferred / retained / random individuals, so the
        # per-environment budget stays exactly N * tau_t.
        self.gating = bool(gating)
        self.gate_eta = float(gate_eta)
        self.gate_rho_min = float(gate_rho_min)
        self.gate_rho_max = float(gate_rho_max)
        assert 0.0 <= self.gate_rho_min <= self.gate_rho_max <= 1.0
        # ---- v3.2 (plan 4.1): online change detection ----
        # When on, the response does not assume the change instant is known:
        # the algorithm decides whether a change happened from the probe series
        # it observes, and the run reports detection latency / false alarms /
        # misses against the driver's ground truth.
        self.detect_change = bool(detect_change)
        self.detect_threshold = detect_threshold
        self.detect_window = int(detect_window)
        self._residuals = []            # historical prediction residuals
        self._rho_series = []           # rho_k actually used
        self._cov_cond = []             # condition number of the covariance
        # diagnostics
        self.alpha_series = []       # alpha_k actually used, per change
        self.movement_series = []    # normalized movement signal
        self.noise_series = []       # bootstrap noise floor
        self.ratio_series = []       # movement / noise
        # normalized-space feature history [(mu, sigma), ...]
        self._feature_history = []
        self.transfer_strategy = transfer_strategy
        self._probe = None
        self._prev_probe_F = None
        self._mu_prev = None    # normalized mean of the last environment
        self._sigma_prev = None  # normalized std of the last environment
        self._dk_sum = 0.0      # running sum of |D_k| for normalization
        self._dk_n = 0
        self._n_changes = 0     # how many environment changes have occurred
        self._alpha_used = float(trend_alpha)   # last alpha actually applied
        # v3.2 state
        self._mu_prev_fused = None      # fused centre of the last environment
        self._sigma_prev_fused = None
        self._mu_hat_prev = None        # last predicted centre (for residuals)
        self._lowrank_energy = None     # variance fraction kept by low-rank

    # ---- v3.1 (plan 3.1): online trend-step adaptation ----
    def centre_noise(self, nd_X_norm, rng):
        """Bootstrap uncertainty of the centre estimate, per-dimension.

        Plan 3.1 asks for "centre-estimate uncertainty or bootstrap stability"
        as the NOISE signal.  We report the mean (over dimensions) of the
        bootstrap standard deviation of each coordinate of the centre, in
        normalized units, so it is directly comparable to a per-dimension
        displacement.
        """
        n = nd_X_norm.shape[0]
        if n < 2:
            return float("inf")
        B = max(2, self.noise_samples)
        # B x n resample indices with replacement
        idx = rng.integers(0, n, size=(B, n))
        means = nd_X_norm[idx].mean(axis=1)          # B x D
        sd = means.std(axis=0, ddof=1)               # D
        return float(np.mean(sd))

    def adaptive_step(self, mu_e, mu_prev, nd_X_norm, rng):
        """Return ``(alpha_k, movement, noise, ratio)`` for this change.

        movement  = ||mu_e - mu_prev||_2 / sqrt(D)     (per-dimension units)
        noise     = mean_j bootstrap-SD(centre_j)      (per-dimension units)
        ratio     = movement / max(noise, eps)
        s_k       = clip((ratio - gamma) / (2*gamma), 0, 1)
        alpha_k   = alpha_min + (alpha_max - alpha_min) * s_k

        Interpretation, matching the plan's response table:
          * ratio ~<= gamma  -> movement is within noise -> s_k = 0 -> alpha_min
            (avoids amplifying sampling noise into a systematic offset)
          * ratio >> gamma   -> sustained real movement -> s_k -> 1 -> alpha_max
          * gamma is the registered "one noise unit" threshold; it is swept in
            a separate calibration run and must NOT be chosen on the final
            test set (plan 3.1).
        """
        D = mu_e.shape[0]
        movement = float(np.linalg.norm(mu_e - mu_prev) / np.sqrt(D))
        noise = self.centre_noise(nd_X_norm, rng)
        if not np.isfinite(noise) or noise <= 0.0:
            ratio = float("inf") if movement > 0 else 0.0
        else:
            ratio = movement / noise
        span = max(2.0 * self.gamma, 1e-12)
        s = float(np.clip((ratio - self.gamma) / span, 0.0, 1.0))
        alpha = self.alpha_min + (self.alpha_max - self.alpha_min) * s
        return alpha, movement, noise, ratio

    # ---- v3.2 (plan 3.2): covariance estimation and correlated sampling ----
    def estimate_covariance(self, nd_X_norm):
        """Return ``(cov, cond)`` for the non-dominated set, in normalized space.

        ``diag``    -> diagonal matrix from the per-dimension variance
        ``full``    -> shrinkage-regularised full covariance (Ledoit-Wolf when
                       scikit-learn is available, otherwise a simple ridge)
        ``lowrank`` -> the full covariance is still returned; the RANK is
                       applied at sampling time via PCA truncation

        The condition number is recorded because the plan asks to watch whether
        the cost of a full covariance is repaid -- a badly conditioned estimate
        is the main reason it would not be.
        """
        Z = np.asarray(nd_X_norm, dtype=float)
        n, D = Z.shape
        if self.cov_mode == "diag" or n < 2:
            return np.diag(np.maximum(Z.var(axis=0, ddof=0), 1e-12)), 1.0
        if self.shrinkage == "auto":
            try:
                from sklearn.covariance import LedoitWolf
                cov = LedoitWolf().fit(Z).covariance_
            except Exception:                       # pragma: no cover
                cov = np.cov(Z, rowvar=False, ddof=0)
        elif isinstance(self.shrinkage, (int, float)):
            s = float(self.shrinkage)
            cov = (1.0 - s) * np.cov(Z, rowvar=False, ddof=0) + \
                s * np.trace(np.cov(Z, rowvar=False, ddof=0)) / D * np.eye(D)
        else:
            cov = np.cov(Z, rowvar=False, ddof=0)
        cov = np.atleast_2d(cov)
        # never let the sampler see a non-PSD matrix
        cov = 0.5 * (cov + cov.T)
        w = np.linalg.eigvalsh(cov)
        w = np.maximum(w, 1e-12)
        cond = float(w.max() / w.min())
        return cov + np.eye(D) * 1e-12, cond

    def sample_transfer(self, mu_hat, sigma_e, cov, n, rng):
        """Draw ``n`` transfer candidates in normalized space.

        ``diag``    reproduces the v2.x behaviour exactly (per-dimension
                    independent Gaussian), which is what keeps every frozen
                    artefact reproducible.
        ``full``    multivariate normal with the estimated covariance.
        ``lowrank`` samples in the leading PCA subspace that explains
                    ``rank_energy`` of the variance (or ``low_rank``
                    components), plus an isotropic residual.
        """
        if n <= 0:
            return np.empty((0, mu_hat.shape[0]))
        D = mu_hat.shape[0]
        if self.cov_mode == "diag":
            noise = rng.normal(0.0, 1.0, size=(n, D))
            return np.clip(mu_hat[None, :] + noise * sigma_e[None, :], 0.0, 1.0)
        if self.cov_mode == "full":
            L = self._chol_or_eigh(cov)
            z = rng.standard_normal(size=(n, D))
            return np.clip(mu_hat[None, :] + z @ L.T, 0.0, 1.0)
        # lowrank: keep the leading components, add isotropic residual
        w, V = np.linalg.eigh(cov)
        order = np.argsort(w)[::-1]
        w, V = w[order], V[:, order]
        k = self.low_rank
        if k is None:
            c = np.cumsum(w) / np.sum(w)
            k = int(np.searchsorted(c, self.rank_energy) + 1)
        k = max(1, min(int(k), D))
        retained = float(np.sum(w[:k]) / np.sum(w))
        resid = float(np.sum(w[k:])) / max(D - k, 1) if k < D else 0.0
        z_main = rng.standard_normal(size=(n, k))
        main = (z_main * np.sqrt(np.maximum(w[:k], 0.0))[None, :]) @ V[:, :k].T
        if resid > 0:
            iso = rng.standard_normal(size=(n, D)) * np.sqrt(resid)
        else:
            iso = 0.0
        self._lowrank_energy = retained
        return np.clip(mu_hat[None, :] + main + iso, 0.0, 1.0)

    @staticmethod
    def _chol_or_eigh(cov):
        """Cholesky factor, falling back to an eigen-based square root."""
        try:
            return np.linalg.cholesky(cov)
        except np.linalg.LinAlgError:
            w, V = np.linalg.eigh(cov)
            return V * np.sqrt(np.maximum(w, 1e-12))[None, :]

    # ---- v3.2 (plan 3.3): negative-transfer gating ----
    def gated_ratios(self, residual):
        """Return ``(rho_k, kappa_k, r_k)`` for this environment.

        ``residual`` is the relative error of the previous prediction,
        ``||mu_hat_{k-1} - mu_e_{k-1}|| / (||delta mu|| + eps)``.  A large
        residual means the trend prediction was wrong, i.e. transferred
        individuals are likely misplaced, so the migration budget shrinks
        towards ``gate_rho_min`` and the freed slots go to random immigrants.
        With no history the nominal ``transfer_ratio`` is used.

        Budget neutrality (plan 3.3): the three ratios always sum to 1, so the
        number of individuals -- and therefore the number of objective
        evaluations per environment -- is unchanged; only the MIX moves.
        """
        rho0 = self.transfer_ratio
        if not self.gating or residual is None or not np.isfinite(residual):
            return rho0, self.retain_ratio, self.random_ratio
        e = np.clip(residual, 0.0, None)
        # eta scales how aggressively error is punished
        shrink = e / (e + self.gate_eta)
        rho = rho0 - (rho0 - self.gate_rho_min) * shrink
        rho = float(np.clip(rho, self.gate_rho_min, min(self.gate_rho_max, rho0)))
        kappa = self.retain_ratio
        r = max(0.0, 1.0 - rho - kappa)
        return rho, kappa, r

    def prediction_residual(self, mu_hat_prev, mu_e_now):
        """Relative error of the previous environment's predicted centre."""
        if mu_hat_prev is None or mu_e_now is None:
            return None
        D = mu_e_now.shape[0]
        err = float(np.linalg.norm(mu_hat_prev - mu_e_now) / np.sqrt(D))
        return err

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
        #    v3.1 (plan 3.1): the step size is chosen ONLINE from decision-space
        #    movement evidence.  NOTE the alpha used for the differential is
        #    applied to mu_e (the raw centre), matching the recorded state
        #    history _mu_prev; mu_eff only shifts the anchor.
        if self.adaptive_alpha and self._mu_prev is not None:
            if self.conservative_first_change and self._n_changes == 0:
                # no residual history yet -> no evidence of sustained movement
                alpha_k = self.alpha_min
                movement, noise, ratio = float("nan"), float("nan"), float("nan")
            else:
                alpha_k, movement, noise, ratio = self.adaptive_step(
                    mu_e, self._mu_prev,
                    normalize(nd_X, xl, xu), rng)
            self.alpha_series.append(alpha_k)
            self.movement_series.append(movement)
            self.noise_series.append(noise)
            self.ratio_series.append(ratio)
            self._alpha_used = alpha_k
        else:
            alpha_k = self.trend_alpha
            self._alpha_used = self.trend_alpha

        if self._mu_prev is not None:
            mu_hat = np.clip(
                mu_eff + alpha_k * (mu_e - self._mu_prev),
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
        #    v3.2 (plan 3.4): the RAW and FUSED centres are now kept in separate
        #    state variables.  The trend history must use one consistent kind of
        #    centre; previously _mu_prev stored the raw centre while the trend
        #    term could be anchored on the fused centre, mixing the two.
        mu_current_raw = mu_e.copy()
        sigma_current_raw = sigma_e.copy()
        mu_current_fused = mu_eff.copy()
        sigma_current_fused = sigma_eff.copy()
        self._mu_prev = mu_current_raw
        self._sigma_prev = sigma_current_raw
        # the trend anchor is always drawn from the SAME pool as the history
        self._mu_prev_fused = mu_current_fused
        self._sigma_prev_fused = sigma_current_fused
        if self.multi_source:
            self._feature_history.append(
                (mu_current_raw.copy(), sigma_current_raw.copy()))
            if len(self._feature_history) > self.history_size:
                self._feature_history.pop(0)
        self._n_changes += 1

        # 7. transfer budget: nominal, or gated by the last prediction residual
        residual = self.prediction_residual(self._mu_hat_prev, mu_e)
        rho_k, kappa_k, r_k = self.gated_ratios(residual)
        self._residuals.append(residual)
        self._rho_series.append(rho_k)
        self._mu_hat_prev = mu_hat.copy()

        # 8. sample transfer candidates in [0,1]^D, denormalize once
        n_transfer = int(round(rho_k * N))
        cov_cond = 1.0
        if n_transfer > 0:
            nd_norm = normalize(nd_X, xl, xu)
            cov, cov_cond = self.estimate_covariance(nd_norm)
            self._cov_cond.append(cov_cond)
            transfer_norm = self.sample_transfer(
                mu_hat, sigma_eff, cov, n_transfer, rng)
            transfer_part = denormalize(transfer_norm, xl, xu)
            bad = ~np.isfinite(transfer_part).all(axis=1)
            if bad.any():
                transfer_part[bad] = rng.uniform(xl, xu, size=(bad.sum(), D))
        else:
            transfer_part = np.empty((0, D))

        # 9. retain a slice of the old population
        n_retain = int(round(kappa_k * N))
        if n_retain > 0:
            idx = rng.permutation(N)
            retain_part = prev_X[idx[:n_retain]].copy()
        else:
            retain_part = np.empty((0, D))

        # 10. uniform random immigrants
        #     With gating OFF this reproduces the v2.x arithmetic EXACTLY
        #     (n_random = round(random_ratio * N)), which is what keeps every
        #     frozen artefact reproducible.  With gating ON the three parts must
        #     partition N exactly, so the remainder is used instead -- budget
        #     neutrality is the whole point of plan 3.3.
        if self.gating:
            n_random = max(0, N - n_transfer - n_retain)
        else:
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
            # v3.1 (plan 3.1): which step was actually taken, so the online
            # adaptation is auditable per environment (parameter stability).
            "alpha_k": float(self._alpha_used),
            # v3.2: distribution / gating / conditioning diagnostics
            "cov_mode": self.cov_mode,
            "cov_cond": float(cov_cond),
            "rho_k": float(rho_k),
            "kappa_k": float(kappa_k),
            "r_k": float(r_k),
            "pred_residual": (None if residual is None
                              else float(residual)),
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


class LDTNSGA2Adaptive(LDTNSGA2):
    """LDT-NSGA-II with the ONLINE adaptive trend step (plan 3.1).

    Motivation.  The v2.4/v2.5 experiments established two facts that no single
    fixed coefficient can satisfy at once:

    * on problems whose Pareto set TRANSLATES regularly (FDA1/FDA3/FDA4,
      dMOP2A, dMOP3A), the full step ``alpha = 1.0`` wins on 11 of 12
      configurations -- the default half step systematically under-extrapolates;
    * on the authoritative Type III problem dMOP1A, whose Pareto set is STATIC
      and only the front rotates, the half step beats the full step by 127%
      (0.0251 vs 0.0570), because the centre increment consists almost entirely
      of sampling noise and a full step amplifies it into a systematic offset.

    This class therefore chooses ``alpha_k`` ONLINE, per environment change:

        movement = ||mu_e - mu_prev||_2 / sqrt(D)          (decision space)
        noise    = mean_j bootstrap-SD of centre coordinate j
        ratio    = movement / noise
        s_k      = clip((ratio - gamma) / (2*gamma), 0, 1)
        alpha_k  = alpha_min + (alpha_max - alpha_min) * s_k

    Plan 3.1 insists that the DECISION-space displacement be separated from the
    OBJECTIVE-space change ``D_k``; this implementation uses only the former for
    the step choice (``D_k`` still feeds the optional sigma scaling).

    Defaults: ``alpha_min = 0.5`` (both registered endpoints are supported by
    the existing evidence) and ``alpha_max = 1.0``; ``gamma = 1.0`` means "one
    noise unit".  ``conservative_first_change`` takes the low end on the very
    first change, where no prediction residual exists yet.
    """
    name = "ldt_adaptive"

    def __init__(self, adaptive_alpha=True, alpha_min=0.5, alpha_max=1.0,
                 gamma=1.0, **kwargs):
        super().__init__(adaptive_alpha=adaptive_alpha, alpha_min=alpha_min,
                         alpha_max=alpha_max, gamma=gamma, **kwargs)


class LDTNSGA2Cov(LDTNSGA2):
    """LDT-NSGA-II with a shrinkage-regularised FULL covariance (plan 3.2).

    The v2.x kernel models the non-dominated set with per-dimension mean and
    standard deviation, i.e. a DIAGONAL Gaussian.  It cannot represent
    correlations between decision variables, so when the Pareto set lies on a
    coupled manifold the independent sampler scatters candidates off it.

    ``full`` mode samples ``N(mu_hat, Cov)`` with ``Cov`` estimated by
    Ledoit-Wolf shrinkage (falling back to a ridge).  The plan explicitly warns
    NOT to assume this wins: the non-dominated sample is small relative to D, so
    an unshrunk estimate is ill-conditioned, costly, and can be worse than the
    diagonal model.  ``cov_cond`` is recorded per environment so the argument
    can be settled with data rather than assertion.
    """
    name = "ldt_cov"

    def __init__(self, cov_mode="full", **kwargs):
        super().__init__(cov_mode=cov_mode, **kwargs)


class LDTNSGA2LowRank(LDTNSGA2):
    """LDT-NSGA-II with a low-rank (PCA-subspace) distribution (plan 3.2).

    Samples inside the leading eigenspace that explains ``rank_energy`` of the
    variance (or a fixed ``low_rank``), plus a small isotropic residual so the
    sampler cannot collapse onto a degenerate subspace.  Intended for the
    high-dimensional end, where the full covariance is neither estimable nor
    affordable; the plan notes the rank / explained-variance threshold is itself
    a hyper-parameter that must be reported.
    """
    name = "ldt_lowrank"

    def __init__(self, cov_mode="lowrank", **kwargs):
        super().__init__(cov_mode=cov_mode, **kwargs)


class LDTNSGA2Gated(LDTNSGA2):
    """LDT-NSGA-II with negative-transfer gating (plan 3.3).

    The migration budget ``rho_k`` shrinks towards ``gate_rho_min`` as the
    previous trend prediction turns out to have been wrong, and the freed slots
    become random immigrants.  The retained fraction is held at
    ``retain_ratio``; the three always sum to 1 so the number of objective
    evaluations per environment is unchanged.

    ``gate_mode="residual"``  uses only the historical prediction residual.
    ``gate_mode="joint"``     also considers the objective-space change D_k,
    which is the plan's strategy B.
    """
    name = "ldt_gated"

    def __init__(self, gating=True, gate_mode="residual", **kwargs):
        super().__init__(gating=gating, **kwargs)
        self.gate_mode = str(gate_mode)
        assert self.gate_mode in ("residual", "joint"), self.gate_mode
        self._last_Dk = None

    def gated_ratios(self, residual):
        rho, kappa, r = super().gated_ratios(residual)
        if self.gate_mode == "joint" and self._dk_n > 0:
            # a large objective-space change with a small decision-space
            # residual is the "front-only moved" case: shrink a little further
            dk_ref = self._dk_sum / max(self._dk_n, 1)
            dk_rel = abs(self._last_Dk or 0.0) / dk_ref if dk_ref > 1e-12 else 0.0
            if dk_rel > 1.0:
                rho = max(self.gate_rho_min, rho * (1.0 / dk_rel))
                r = max(0.0, 1.0 - rho - kappa)
        return rho, kappa, r

    def respond(self, prev_X, prev_F, problem, env_index, rng):
        new_X, info = super().respond(prev_X, prev_F, problem, env_index, rng)
        self._last_Dk = info.get("Dk")
        return new_X, info


class LDTNSGA2Full(LDTNSGA2):
    """v3.2 bundle: adaptive step + gating + multi-source, fixed 50/20/30 base.

    A single switchable object so the §7 roadmap features can be ablated one at
    a time (plan 3.1/3.3/3.4) instead of being an all-or-nothing add-on.
    """
    name = "ldt_v32"

    def __init__(self, **kwargs):
        kwargs.setdefault("adaptive_alpha", True)
        kwargs.setdefault("gating", True)
        kwargs.setdefault("multi_source", True)
        super().__init__(**kwargs)


# Backward-compatible alias for existing experiment scripts
LDMTNSGA2 = LDTNSGA2
