"""Extended change patterns for dynamic multi-objective benchmarks (plan 4.2).

Plan 4.2 asks for change patterns beyond the regular Type I/II/III sweep, so
that "full-step extrapolation is best" can be stress-tested rather than
assumed.  The existing six problems all advance time as

    t = floor(g / tau_t) / n_t

i.e. equally spaced, monotone, abrupt, noise-free environment changes.  Each
class below relaxes exactly ONE of those assumptions so the effect is
attributable:

  ``IrregularFDA1``   non-equally-spaced changes (jittered breakpoints)
  ``ReversedFDA1``    the PS reverses direction partway through the run
  ``GradualFDA1``     each environment ramps t continuously (gradual, not abrupt)
  ``NoisyFDA1``       objective evaluations carry measurement noise

All inherit ``time_at`` / ``_evaluate`` from the corresponding base problem, so
the objective and true-front definitions are IDENTICAL to their regular
counterparts.  Only the temporal or measurement model differs -- which is what
makes them usable as controlled variants rather than new problems.

Provenance note: these are project-defined stress variants.  They are NOT
literature benchmarks and must not be cited as such; the thesis reports them
under an explicit "project variant" label.
"""
import numpy as np

from .fda1 import FDA1


class IrregularFDA1(FDA1):
    """FDA1 with non-equally-spaced change points.

    The nominal environment length is ``tau_t``; actual breakpoints are a fixed
    pseudo-random jitter of at most ``jitter`` (fraction of tau_t) around the
    nominal grid, drawn once from ``IRREGULAR_SEED`` so every algorithm and every
    seed sees the SAME change schedule.  Only the schedule is irregular, not the
    dynamics.

    ``t`` is the environment INDEX divided by ``n_t`` (so the PS magnitude per
    step is unchanged); what becomes irregular is WHEN the steps happen.
    """

    IRREGULAR_SEED = 4242
    jitter = 0.4

    def _build_breakpoints(self):
        """Deterministic jittered change schedule (built once, lazily).

        Built lazily because ``DynamicProblem.__init__`` calls ``time_at(0)``
        before the subclass constructor body runs.
        """
        if getattr(self, "_breakpoints", None) is not None:
            return self._breakpoints
        rng = np.random.default_rng(self.IRREGULAR_SEED)
        n_env = max(2, int(np.ceil(300 / max(self.tau_t, 1))) + 2)
        span = self.tau_t * self.jitter
        offsets = rng.uniform(-span, span, size=n_env)
        bp = np.round(np.arange(n_env) * self.tau_t + offsets).astype(int)
        bp[0] = 0
        # strictly increasing, each environment at least 1 generation long
        bp = np.maximum(bp, np.arange(n_env))
        self._breakpoints = np.maximum.accumulate(bp)
        return self._breakpoints

    def time_at(self, g):
        bp = self._build_breakpoints()
        idx = int(np.searchsorted(bp, int(g), side="right") - 1)
        return max(0, idx) / float(self.n_t)


class ReversedFDA1(FDA1):
    """FDA1 whose PS reverses direction halfway through the run.

    For the first ``reverse_fraction`` of the run the environment index counts
    up as usual; afterwards it counts back down, so ``mu_e - mu_prev`` keeps its
    magnitude but flips SIGN.  A fixed full-step extrapolator therefore
    over-shoots exactly when the reversal happens, which is the plan's
    "trend reversal" stress case.
    """

    reverse_fraction = 0.5

    def time_at(self, g):
        total = max(1, int(300 // max(self.tau_t, 1)))
        idx = int(g) // max(self.tau_t, 1)
        pivot = int(round(total * self.reverse_fraction))
        if idx > pivot:
            idx = max(0, 2 * pivot - idx)
        return idx / float(self.n_t)


class GradualFDA1(FDA1):
    """FDA1 where t ramps CONTINUOUSLY inside each environment.

    ``t = g / (tau_t * n_t)`` instead of ``floor(g / tau_t) / n_t``.  The
    environment segments still exist (the driver still resets the optimizer
    every ``tau_t`` generations), but the objective function drifts smoothly
    instead of jumping.  This separates "response to a discrete change" from
    "tracking a moving target".
    """

    def time_at(self, g):
        return float(g) / float(max(1, self.tau_t) * max(1, self.n_t))


class NoisyFDA1(FDA1):
    """FDA1 with multiplicative-plus-additive Gaussian measurement noise.

    Noise is applied in ``_evaluate`` (i.e. to what the optimizer and the
    IGD/HV readers observe) but NOT to the analytical true front, so the
    reference front stays exact while every evaluation is perturbed.  The noise
    stream is seeded per problem instance so a run is still reproducible.

    ``noise_sigma`` is relative to the objective scale of FDA1 (f1 in [0,1],
    f2 in [0,1]); 5% is large enough to matter for D_k detection but small
    enough that the problem stays solvable.
    """

    noise_sigma = 0.05
    NOISE_SEED = 909

    def __init__(self, *args, noise_sigma=None, **kwargs):
        super().__init__(*args, **kwargs)
        if noise_sigma is not None:
            self.noise_sigma = float(noise_sigma)
        self._noise_rng = np.random.default_rng(self.NOISE_SEED)

    def _evaluate(self, X, out, *args, **kwargs):
        super()._evaluate(X, out, *args, **kwargs)
        F = np.asarray(out["F"], dtype=float)
        scale = np.maximum(np.abs(F), 1e-6)
        out["F"] = F + self._noise_rng.normal(
            0.0, self.noise_sigma, size=F.shape) * scale


#: registry of the project-defined stress variants (plan 4.2)
EXTENDED_PROBLEMS = {
    "IRREG_FDA1": IrregularFDA1,
    "REV_FDA1": ReversedFDA1,
    "GRAD_FDA1": GradualFDA1,
    "NOISY_FDA1": NoisyFDA1,
}

#: what each variant relaxes, for tables and the thesis
EXTENDED_LABELS = {
    "IRREG_FDA1": "非等间隔变化",
    "REV_FDA1": "趋势方向反转",
    "GRAD_FDA1": "渐变（连续漂移）",
    "NOISY_FDA1": "目标评估带噪声",
}
