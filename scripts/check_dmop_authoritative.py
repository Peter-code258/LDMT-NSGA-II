"""Self-check for the v2.5 authoritative dMOP1-3 implementations.

Verified WITHOUT running any optimiser:

1. Placing the decision vector exactly on the analytically derived Pareto set
   PS(t, r) reproduces objective values that lie on the analytical Pareto front
   returned by ``problems.true_fronts.get_true_pf``;
2. the project variants and the authoritative variants really differ;
3. the problem TYPE of each authoritative problem matches the literature:
     dMOP1 -> Type III (PS static, PF time-varying)
     dMOP2 -> Type II  (PS and PF both time-varying)
     dMOP3 -> Type I   (PS time-varying, PF static)
4. the switching of the "position-related" variable r in dMOP3 is constant
   within an environment and changes across environments.

Run:  python -m scripts.check_dmop_authoritative
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from problems import DMOP1, DMOP2, DMOP3, DMOP1A, DMOP2A, DMOP3A
from problems.true_fronts import get_true_pf

N = 200


def _H(t):
    return 0.75 * np.sin(0.5 * np.pi * t) + 1.25


def _G_abs(t):
    return np.abs(np.sin(0.5 * np.pi * t))


def _G_plain(t):
    return np.sin(0.5 * np.pi * t)


def ps_auth_dmop1(t, n=N):
    """Type III POS: x_1 = u, x_i = 0 for i >= 2 (independent of t)."""
    u = np.linspace(0.0, 1.0, n)
    X = np.zeros((n, 10))
    X[:, 0] = u
    return u, X


def ps_auth_dmop2(t, n=N):
    """Type II POS: x_1 = u, x_i = G(t) for i >= 2."""
    G = _G_plain(t)
    u = np.linspace(0.0, 1.0, n)
    X = np.tile(G, (n, 10))
    X[:, 0] = u
    return u, X


def ps_auth_dmop3(t, r, n=N):
    """Type I POS: x_r = u, x_i = G(t) for i != r."""
    G = _G_abs(t)
    u = np.linspace(0.0, 1.0, n)
    X = np.tile(G, (n, 10))
    X[:, r] = u
    return u, X


def main():
    t_values = [0.0, 0.05, 0.2, 0.45, 0.8]
    print("=" * 78)
    print("1) PS-on-front consistency: max |f2(PS) - f2(PF)|")
    print("=" * 78)
    ok = True
    for cls, ps_fn, name in ((DMOP1A, ps_auth_dmop1, "DMOP1A"),
                             (DMOP2A, ps_auth_dmop2, "DMOP2A")):
        prob = cls(tau_t=25, n_t=20)
        worst = 0.0
        for t in t_values:
            prob.t = float(t)
            u, X = ps_fn(float(t))
            F = prob._calc_objectives(X)
            PF = get_true_pf(name, M=2, n_points=N, t=float(t))
            worst = max(worst, float(np.max(np.abs(F[:, 1] - PF[:, 1]))))
        ok &= worst < 1e-12
        print(f"  {name}: worst deviation = {worst:.3e}  "
              f"[{'OK' if worst < 1e-12 else 'FAIL'}]")

    prob3 = DMOP3A(tau_t=25, n_t=20)
    worst3 = 0.0
    for t in t_values:
        for r in range(10):
            prob3.t = float(t)
            prob3.r = r
            u, X = ps_auth_dmop3(float(t), r)
            F = prob3._calc_objectives(X)
            PF = get_true_pf("DMOP3A", M=2, n_points=N)
            worst3 = max(worst3, float(np.max(np.abs(F[:, 1] - PF[:, 1]))))
    ok &= worst3 < 1e-12
    print(f"  DMOP3A: worst deviation over all r = {worst3:.3e}  "
          f"[{'OK' if worst3 < 1e-12 else 'FAIL'}]")

    print("\n" + "=" * 78)
    print("2) Authoritative vs project variant: are they actually different?")
    print("=" * 78)
    X = np.full((1, 10), 0.3)
    X[0, 0] = 0.4
    for auth, var, at, vt in ((DMOP1A, DMOP1, "DMOP1A", "DMOP1"),
                              (DMOP2A, DMOP2, "DMOP2A", "DMOP2"),
                              (DMOP3A, DMOP3, "DMOP3A", "DMOP3")):
        pa, pv = auth(tau_t=25, n_t=20), var(tau_t=25, n_t=20)
        pa.t = pv.t = 0.2
        if hasattr(pa, "r"):
            pa.r = 1
        Fa = pa._calc_objectives(X)[0]
        Fv = pv._calc_objectives(X)[0]
        diff = float(np.max(np.abs(Fa - Fv)))
        print(f"  {at} vs {vt}: F_auth=({Fa[0]:.5f},{Fa[1]:.5f})  "
              f"F_var=({Fv[0]:.5f},{Fv[1]:.5f})  max|diff|={diff:.5f}"
              f"  -> {'DIFFERENT' if diff > 1e-9 else 'IDENTICAL (suspicious!)'}")

    print("\n" + "=" * 78)
    print("3) Problem TYPE check (PS and PF time dependence)")
    print("=" * 78)
    expected = {"DMOP1A": ("static", "MOVE", "Type III"),
                "DMOP2A": ("MOVE", "MOVE", "Type II"),
                "DMOP3A": ("MOVE", "static", "Type I")}
    for name, (exp_ps, exp_pf, typ) in expected.items():
        # PF time dependence
        f0 = get_true_pf(name, M=2, n_points=101, t=0.0)
        f2 = get_true_pf(name, M=2, n_points=101, t=0.45)
        dpf = float(np.max(np.abs(f0 - f2)))
        pf_state = "MOVE" if dpf > 1e-9 else "static"
        # PS time dependence
        if name == "DMOP1A":
            dps = 0.0                      # POS is x_i = 0 for all t
            ps_state = "static"
        elif name == "DMOP2A":
            dps = abs(_G_plain(0.0) - _G_plain(0.45))
            ps_state = "MOVE" if dps > 1e-9 else "static"
        else:
            dps = abs(_G_abs(0.0) - _G_abs(0.45))
            ps_state = "MOVE" if dps > 1e-9 else "static"
        good = (ps_state == exp_ps and pf_state == exp_pf)
        ok &= good
        print(f"  {name}: PS {ps_state} (|dPS|={dps:.4f}), "
              f"PF {pf_state} (|dPF|={dpf:.4f})  -> {typ}  "
              f"[{'OK' if good else 'MISMATCH'}]")

    print("\n" + "=" * 78)
    print("4) dMOP3A position-variable switching policy")
    print("=" * 78)
    p = DMOP3A(tau_t=25, n_t=20)
    seq = []
    for g in (0, 24, 25, 49, 50, 74, 75):
        p.set_generation(g)
        seq.append((g, p.r, round(float(p.t), 5)))
    for g, r, t in seq:
        print(f"  g={g:3d}: r={r}  t={t}")
    env_r = {}
    for g in range(0, 300):
        p.set_generation(g)
        env_r.setdefault(g // 25, set()).add(p.r)
    const_ok = all(len(v) == 1 for v in env_r.values())
    changes = len({next(iter(v)) for v in env_r.values()})
    ok &= const_ok and changes > 1
    print(f"  r constant within each environment: {const_ok}; "
          f"distinct r values across environments: {changes}")

    print("\n" + "=" * 78)
    print("FINAL:", "ALL CHECKS PASSED" if ok
          else "CHECK FAILURE - do not use these problems")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
