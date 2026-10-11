"""PPS (Population Prediction Strategy) baseline.

Simplified version for undergraduate-level comparison:
- Maintains a ring buffer of the last K non-dominated population centres
- Fits a linear trend per dimension:  mu_hat = 2*mu_{e-1} - mu_{e-2}
  (i.e. first-difference extrapolation, equivalent to AR(1) with b=1)
- Shifts the previous population by the predicted delta
- Injects 20 % random immigrants

Reference:
    Zhou A, Jin Y, Zhang Q. A population prediction strategy for evolutionary
    dynamic multiobjective optimization. IEEE Trans Cybern, 2014.
    DOI: 10.1109/TCYB.2013.2276926
"""
import numpy as np

from .distribution_features import (normalize, denormalize, extract_features,
                                    non_dominated_set)
from .environment_manager import ChangeDetector


class PPSNSGA2:
    """PPS-style population-centre prediction for dynamic MOO.

    v2.2: all feature math stays in normalized [0,1]^D space (the previous
    version re-normalized already-normalized means — a coordinate bug); the
    change probe is the problem-owned fixed probe in real bounds.
    """
    name = "pps"

    def __init__(self, history_len=3, immigrant_ratio=0.2, n_probes=10):
        self.history_len = int(history_len)
        self.immigrant_ratio = float(immigrant_ratio)
        self.n_probes = int(n_probes)
        self.mu_history = []  # normalized per-dim centres
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
        # non-dominated subset for feature extraction (normalized space)
        nd_mask = non_dominated_set(prev_X, prev_F)
        nd_X = prev_X[nd_mask] if nd_mask.sum() > 0 else prev_X
        mu_prev, sigma_prev = extract_features(nd_X, xl, xu)  # in [0,1]
        self._ensure_probe(prev_X, problem)
        Dk = self._detect_change(problem)
        # store history (normalized)
        self.mu_history.append(mu_prev)
        if len(self.mu_history) > self.history_len:
            self.mu_history.pop(0)
        # predict next centre: first-difference extrapolation in [0,1]
        if len(self.mu_history) >= 2:
            mu_hat = np.clip(2.0 * self.mu_history[-1] - self.mu_history[-2],
                             0.0, 1.0)
        else:
            mu_hat = np.clip(mu_prev, 0.0, 1.0)
        # shift the whole normalized population by the predicted delta
        Z = normalize(prev_X, xl, xu)
        z_new = np.clip(Z + (mu_hat - mu_prev)[None, :], 0.0, 1.0)
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
                       "n_nd": int(nd_mask.sum()), "Dk": Dk}

    def _resize(self, X, N, xl, xu, rng):
        D = X.shape[1]
        if X.shape[0] > N:
            return X[:N]
        while X.shape[0] < N:
            X = np.concatenate(
                [X, rng.uniform(xl, xu, size=(1, D))], axis=0)
        return X
