"""Distribution features (mean/scale), non-dominated selection, and the
LDMT migration function (AR(1)-style distribution translation with an
adaptively clipped lambda_k)."""
import numpy as np

SIGMA_MIN = 1e-3
_EPS = 1e-8


def normalize(X, xl, xu):
    xl = np.asarray(xl, dtype=float)
    xu = np.asarray(xu, dtype=float)
    span = np.where(xu - xl > _EPS, xu - xl, 1.0)
    return (np.asarray(X, dtype=float) - xl) / span


def denormalize(Z, xl, xu):
    xl = np.asarray(xl, dtype=float)
    xu = np.asarray(xu, dtype=float)
    span = np.where(xu - xl > _EPS, xu - xl, 1.0)
    return np.clip(Z, 0.0, 1.0) * span + xl


def extract_features(X, xl, xu):
    """Return (mu, sigma) of the normalized population (D,)."""
    Z = normalize(X, xl, xu)
    mu = Z.mean(axis=0)
    sigma = Z.std(axis=0, ddof=0)
    return mu, sigma


def non_dominated_set(X, F):
    """Boolean mask of non-dominated individuals (minimization)."""
    F = np.atleast_2d(np.asarray(F, dtype=float))
    n = F.shape[0]
    if n == 0:
        return np.zeros(0, dtype=bool)
    keep = np.ones(n, dtype=bool)
    for i in range(n):
        if not keep[i]:
            continue
        # j is dominated by i if F[j] >= F[i] in all and > in one
        dom = np.all(F <= F[i][None, :] + _EPS, axis=1) & \
            np.any(F < F[i][None, :] - _EPS, axis=1)
        dom[i] = False
        keep[dom] = False
    return keep


class MigrationFunction:
    """Holds the history of (mu_k, sigma_k, D_k) and predicts the next
    distribution.

    lambda_k (adaptive step) formula (documented, default-to-frozen):
        if history has >= 2 probe magnitudes:
            lambda_k = clip(1.0 + alpha * (D_k - D_prev) / (D_prev + 1e-9),
                            0.5, 1.5)
        else (first change): lambda_k = 1.0

    mu_hat  = mu_prev  + lambda_k * (mu_prev  - mu_prev2)
    sigma_hat = max(sigma_min, sigma_prev + lambda_k*(sigma_prev - sigma_prev2))
    """

    def __init__(self, alpha=1.0, sigma_min=SIGMA_MIN):
        self.alpha = float(alpha)
        self.sigma_min = float(sigma_min)
        self.mu_history = []
        self.sigma_history = []
        self.D_history = []

    def record(self, mu, sigma, D_k):
        self.mu_history.append(np.asarray(mu, dtype=float))
        self.sigma_history.append(np.asarray(sigma, dtype=float))
        self.D_history.append(float(D_k))

    def predict(self, mu_prev, mu_prev2, sigma_prev, sigma_prev2,
                D_k, D_prev):
        if len(self.D_history) >= 2 and D_prev > 1e-12:
            lam = 1.0 + self.alpha * (D_k - D_prev) / (D_prev + 1e-9)
            lam = float(np.clip(lam, 0.5, 1.5))
        else:
            lam = 1.0
        mu_hat = mu_prev + lam * (mu_prev - mu_prev2)
        sigma_hat = np.maximum(
            self.sigma_min,
            sigma_prev + lam * (sigma_prev - sigma_prev2))
        return mu_hat, sigma_hat, lam


def transfer_population(prev_X, xl, xu, mu_prev, sigma_prev,
                        mu_hat, sigma_hat):
    """Re-center & re-scale prev_X to the predicted distribution.

    z_new = clip(mu_hat + (sigma_hat/(sigma_prev+eps))*(Z - mu_prev), 0, 1).
    """
    Z = normalize(prev_X, xl, xu)
    ratio = sigma_hat / (sigma_prev + _EPS)
    z_new = mu_hat + ratio * (Z - mu_prev)
    z_new = np.clip(z_new, 0.0, 1.0)
    return denormalize(z_new, xl, xu)
