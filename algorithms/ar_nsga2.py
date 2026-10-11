"""AR(1)-prediction NSGA-II baseline.

Maintains a short history of the non-dominated-population mean (per dim) and
fits an AR(1) model ``mu_k = a + b*mu_{k-1}`` per dimension via least squares.
On change, predicts ``mu_hat``, shifts the previous population, and injects
20% random immigrants. Falls back to first-difference extrapolation when the
AR(1) fit is unstable.
"""
import numpy as np

from .distribution_features import (normalize, denormalize, extract_features,
                                    non_dominated_set)
from .environment_manager import ChangeDetector


class ARNSGA2:
    name = "ar1_predictive"

    def __init__(self, history_len=5, immigrant_ratio=0.2, n_probes=10):
        self.history_len = int(history_len)
        self.immigrant_ratio = float(immigrant_ratio)
        # v2.4: n_probes was accepted but never stored — _ensure_probe would
        # crash with AttributeError on problems without a shared probe.
        self.n_probes = int(n_probes)
        self.mu_history = []
        self.sigma_history = []
        self.fallback_count = 0
        self._probe = None
        self._prev_probe_F = None

    def _ensure_probe(self, prev_X, problem):
        shared = getattr(problem, "probe", None)
        if shared is not None:
            self._probe = np.asarray(shared, dtype=float)
            return self._probe
        if self._probe is None or self._probe.shape[1] != prev_X.shape[1]:
            n = min(self.n_probes, prev_X.shape[0])
            self._probe = np.asarray(prev_X[:n], dtype=float).copy()
        return self._probe

    def _detect_change(self, problem):
        probe = self._probe
        if probe is None:
            return 0.0
        # evaluate probe under current (new) env
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
        N, D = prev_X.shape
        xl = np.asarray(problem.xl, dtype=float)
        xu = np.asarray(problem.xu, dtype=float)
        # non-dominated subset for feature extraction
        nd_mask = non_dominated_set(prev_X, prev_F)
        nd_X = prev_X[nd_mask] if nd_mask.sum() > 0 else prev_X
        mu_prev, sigma_prev = extract_features(nd_X, xl, xu)
        self._ensure_probe(prev_X, problem)
        Dk = self._detect_change(problem)
        # predict mu_hat
        self.mu_history.append(mu_prev)
        self.sigma_history.append(sigma_prev)
        if len(self.mu_history) > self.history_len:
            self.mu_history.pop(0)
            self.sigma_history.pop(0)
        fallback = False
        if len(self.mu_history) >= 3:
            # Per-dimension AR(1) OLS: mu_k[d] = a[d] + b[d]*mu_{k-1}[d].
            # (A joint multivariate regression would mix unrelated dimensions
            #  and produce meaningless centres — fixed.)
            hist = np.array(self.mu_history)
            x_prev = hist[:-1]
            x_cur = hist[1:]
            xp_mu = x_prev.mean(axis=0)
            xc_mu = x_cur.mean(axis=0)
            denom = np.sum((x_prev - xp_mu) ** 2, axis=0)
            # where denom is ~0 the slope is undefined -> fall back to 1.0
            safe = denom > 1e-12
            b = np.where(safe,
                         np.sum((x_prev - xp_mu) * (x_cur - xc_mu), axis=0)
                         / np.where(safe, denom, 1.0),
                         1.0)
            a = xc_mu - b * xp_mu
            if np.any(np.isnan(b)) or np.any(np.abs(b) > 3.0):
                fallback = True
                mu_hat = mu_prev + (mu_prev - self.mu_history[-2])
            else:
                mu_hat = a + b * mu_prev
        else:
            fallback = True
            if len(self.mu_history) >= 2:
                mu_hat = mu_prev + (mu_prev - self.mu_history[-2])
            else:
                mu_hat = mu_prev
        if fallback:
            self.fallback_count += 1
        sigma_hat = sigma_prev
        # shift prev population to predicted mean
        Z = normalize(prev_X, xl, xu)
        z_new = np.clip(mu_hat + (Z - mu_prev), 0.0, 1.0)
        new_X = denormalize(z_new, xl, xu)
        # immigrants (K1 fix): reserve room BEFORE appending so random
        # immigrants survive _resize. Previously "shift N rows + append
        # immigrants + truncate to first N" dropped every immigrant, making
        # immigrant_ratio a dead parameter (pure-prediction behaviour).
        n_imm = int(round(self.immigrant_ratio * N))
        immigrants = rng.uniform(xl, xu, size=(n_imm, D))
        new_X = np.concatenate([new_X[:N - n_imm], immigrants], axis=0)
        new_X = self._resize(new_X, N, xl, xu, rng)
        # v2.4 (K9 fix): real diagnostic instead of hard-coded 0.0 — fraction
        # of individuals sitting ON the variable bounds after the response.
        bclip = float(np.mean(np.any(
            (new_X <= xl + 1e-9) | (new_X >= xu - 1e-9), axis=1)))
        return new_X, {"lambda_k": 1.0, "boundary_clip_ratio": bclip,
                         "n_nd": int(nd_mask.sum()), "fallback": fallback,
                         "Dk": Dk}

    def _resize(self, X, N, xl, xu, rng):
        D = X.shape[1]
        if X.shape[0] > N:
            return X[:N]
        while X.shape[0] < N:
            X = np.concatenate(
                [X, rng.uniform(xl, xu, size=(1, D))], axis=0)
        return X
