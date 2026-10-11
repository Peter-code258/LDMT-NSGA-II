"""Environment time tracking and change detection."""
import numpy as np


class EnvironmentManager:
    """Tracks the FDA discretely-sampled gradual environment.

    ``time_at(g)`` -> t; ``is_change(g)`` True when floor(g/tau_t) increments;
    ``env_index_at(g)`` -> the 0-based environment index k.
    """

    def __init__(self, tau_t, n_t):
        self.tau_t = int(tau_t)
        self.n_t = int(n_t)

    def time_at(self, g):
        return np.floor(g / self.tau_t) / self.n_t

    def env_index_at(self, g):
        return int(np.floor(g / self.tau_t))

    def is_change(self, g):
        """True at the FIRST generation of a new environment (g % tau_t == 0
        and g > 0)."""
        return g > 0 and (g % self.tau_t == 0)


class ChangeDetector:
    """Diagnoses the normalized objective-space change magnitude D_k.

    Holds a fixed probe set in normalized [0,1]^D; ``detect`` compares the
    previous objective vector ``prev_F`` (P, n_obj) against the re-evaluated
    ``new_F`` under the new environment. Returns a scalar D_k = mean abs diff
    per objective, normalized by per-objective scale. Probes are never used as
    solutions.
    """

    def __init__(self, probe_set, n_obj):
        self.probe = np.asarray(probe_set, dtype=float)
        self.n_obj = int(n_obj)

    def detect(self, prev_F, new_F):
        prev_F = np.atleast_2d(np.asarray(prev_F, dtype=float))
        new_F = np.atleast_2d(np.asarray(new_F, dtype=float))
        if prev_F.size == 0 or new_F.size == 0:
            return 0.0
        scale = np.abs(prev_F).max(axis=0)
        scale = np.where(scale > 1e-12, scale, 1.0)
        d = np.mean(np.abs(new_F - prev_F) / scale[None, :])
        return float(d)
