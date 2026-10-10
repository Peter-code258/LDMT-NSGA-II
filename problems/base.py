"""Base class for dynamic multi-objective test problems (FDA / dMOP).

Time rule (Farina 2004 / PlatEMO master): ``t = floor(g / tau_t) / n_t``
where ``g`` is the 0-based generation index. Subclasses implement
``_calc_objectives`` using the current ``self.t``.
"""
import numpy as np
from pymoo.core.problem import Problem


class DynamicProblem(Problem):
    """Base class holding the FDA discrete-time parameters."""

    def __init__(self, n_var, n_obj, xl, xu, tau_t=10, n_t=10):
        self.tau_t = int(tau_t)
        self.n_t = int(n_t)
        self._g = 0
        self.t = self.time_at(0)
        super().__init__(n_var=n_var, n_obj=n_obj,
                         xl=np.asarray(xl, dtype=float),
                         xu=np.asarray(xu, dtype=float))

    # ---- time helpers ----
    def time_at(self, g):
        """Return the continuous time ``t`` for a 0-based generation ``g``."""
        return np.floor(g / self.tau_t) / self.n_t

    def set_generation(self, g):
        """Advance the problem to generation ``g`` (0-based) and refresh ``t``."""
        self._g = int(g)
        self.t = self.time_at(self._g)
        return self.t

    @property
    def generation(self):
        return self._g

    def _evaluate(self, X, out, *args, **kwargs):
        out["F"] = self._calc_objectives(np.atleast_2d(X).astype(float))

    def _calc_objectives(self, X):  # pragma: no cover - abstract
        raise NotImplementedError
