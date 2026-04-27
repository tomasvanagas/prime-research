"""B2 - Automorphic L-function basis for pi(x) reconstruction.

Frontier attack from ATTACK_VECTORS.md §B2. Tests whether linear functionals
of Hecke partial sums of the level-1 weight-12 cusp form Delta can
reconstruct y(x) = pi(x) - Li(x) with sub-sqrt(x) residual.

a(n) = tau(n)/n^{11/2}, |a(p)| <= 2 (Deligne).

PRE-STATED FALSIFICATION (A/B/C, see results.md for full):
    A-grade: F_tau test residual scales sub-sqrt(x) AND beats F_random
             matched-multiplicative null by > 3 sigma.
    B-grade: F_tau test residual matches F_random matched-multiplicative
             null within sample noise (Sato-Tate killing made quantitative);
             OR explicit numerical bound F_tau / F_random ratio at finite N.
    C/E: collapses to E7.1 / E1.10 / E3.13 (GUE-killing arguments) in a
             new clothing.

Cross-domain: Hecke eigenvalue partial sums; Selberg trace formula intuition.

Mathematician channel: Sarnak.
"""

from __future__ import annotations
import json
import time
from pathlib import Path

import numpy as np
import mpmath
from mpmath import mpc
from sympy import primepi, factorint, isprime


# -------- Ramanujan tau via eta^24 expansion --------

def conv_obj(a: np.ndarray, b: np.ndarray, max_len: int) -> np.ndarray:
    a = np.asarray(a, dtype=object); b = np.asarray(b, dtype=object)
    res = np.zeros(max_len, dtype=object)
    if len(a) > len(b):
        a, b = b, a
    for i in range(len(a)):
        ai = a[i]
        if ai == 0: continue
        kmax = min(len(b), max_len - i)
        if kmax > 0:
            res[i:i+kmax] = res[i:i+kmax] + ai * b[:kmax]
    return res


def tau_via_eta24(N: int) -> np.ndarray:
    eta = np.zeros(N+1, dtype=object); eta[0] = 1
    k = 1
    while True:
        p1 = k*(3*k - 1)//2; p2 = k*(3*k + 1)//2
        if p1 > N: break
        sign = (-1)**k
        eta[p1] = sign
        if p2 <= N: eta[p2] = sign
        k += 1
    eta2 = conv_obj(eta, eta, N+1)
    eta4 = conv_obj(eta2, eta2, N+1)
    eta8 = conv_obj(eta4, eta4, N+1)
    eta16 = conv_obj(eta8, eta8, N+1)
    eta24 = conv_obj(eta16, eta8, N+1)
    tau = np.zeros(N+1, dtype=object)
    tau[1:N+1] = eta24[0:N]
    return tau


# -------- a(n) for various ensembles --------

def normalised_a(tau: np.ndarray) -> np.ndarray:
    """a(n) = tau(n) / n^{11/2}."""
    N = len(tau) - 1
    a = np.zeros(N+1)
    for n in range(1, N+1):
        a[n] = float(tau[n]) / (n**5.5)
    return a


def make_random_iid_satotate(N: int, seed: int) -> np.ndarray:
    """Independent Sato-Tate-distributed a_rand(n) for ALL n in [1, N].
    Density (2/pi) sin^2(theta) -> a = 2 cos(theta) on [-2, 2]."""
    rng = np.random.default_rng(seed)
    # Sample theta with Sato-Tate density: rejection on uniform with bound 2/pi
    a = np.zeros(N+1)
    samples_needed = N
    samples_collected = 0
    while samples_collected < N:
        batch = max(N - samples_collected, 1024)
        theta = rng.uniform(0, np.pi, size=batch)
        u = rng.uniform(0, 1, size=batch)
        # accept if u < sin^2(theta) (since envelope max is 1; pdf is (2/pi)*sin^2)
        accept = u < np.sin(theta)**2
        good = theta[accept]
        nshould = min(N - samples_collected, len(good))
        a[1+samples_collected:1+samples_collected+nshould] = 2*np.cos(good[:nshould])
        samples_collected += nshould
    return a


def make_random_multiplicative(N: int, seed: int) -> np.ndarray:
    """Random multiplicative function: pick a_rand(p) ~ Sato-Tate at primes,
    extend via Hecke recursion at prime powers, then multiplicativity at
    composites. Matches structural form of true a(n)."""
    rng = np.random.default_rng(seed)
    # 1) sample a(p) ~ Sato-Tate at primes
    a_p = {}
    primes = []
    for p in range(2, N+1):
        if isprime(p):
            primes.append(p)
            # sample theta with Sato-Tate density via rejection
            while True:
                theta = rng.uniform(0, np.pi)
                u = rng.uniform(0, 1)
                if u < np.sin(theta)**2:
                    a_p[p] = 2.0 * np.cos(theta)
                    break
    # 2) Hecke recursion at prime powers
    a_pk = {p: {1: a_p[p]} for p in primes}
    for p in primes:
        # for prime powers: a(p^{k+1}) = a(p)*a(p^k) - a(p^{k-1})
        # NB: this is the NORMALISED Hecke recursion (since we use a not tau).
        # Using unnormalised tau: tau(p^{k+1}) = tau(p) tau(p^k) - p^11 tau(p^{k-1}).
        # After dividing by n^{11/2} = (p^{k+1})^{11/2}:
        # a(p^{k+1}) = a(p) * a(p^k) - a(p^{k-1})  (the p^11 factor cancels)
        ap1 = a_p[p]
        a_pk[p][0] = 1.0  # a(p^0) = a(1) = 1
        k = 1
        while p**(k+1) <= N:
            a_pk[p][k+1] = ap1 * a_pk[p][k] - a_pk[p][k-1]
            k += 1
    # 3) multiplicativity at composites
    a = np.zeros(N+1)
    a[1] = 1.0
    for n in range(2, N+1):
        f = factorint(n)
        v = 1.0
        for p, k in f.items():
            v *= a_pk[p][k]
        a[n] = v
    return a


# -------- Li and zeta zeros --------

def li_real(x: float) -> float:
    return float(mpmath.li(x))


def li_complex(x: float, rho: complex) -> complex:
    z = mpc(rho) * mpmath.log(x)
    return complex(mpmath.ei(z))


# -------- Feature matrices --------

def build_F_twisted(x_grid: np.ndarray, weights: np.ndarray,
                    t_grid: np.ndarray) -> np.ndarray:
    assert np.all(np.diff(x_grid) > 0), "x_grid must be strictly increasing"
    N = len(weights) - 1
    log_n = np.zeros(N+1)
    log_n[1:] = np.log(np.arange(1, N+1))
    K = len(t_grid)
    F = np.zeros((len(x_grid), 2*K))
    n_idx = 1
    cum_cos = np.zeros(K); cum_sin = np.zeros(K)
    for i, x in enumerate(x_grid):
        x_int = int(x)
        while n_idx <= min(x_int, N):
            phase = t_grid * log_n[n_idx]
            cum_cos += weights[n_idx] * np.cos(phase)
            cum_sin += weights[n_idx] * np.sin(phase)
            n_idx += 1
        F[i, 0::2] = cum_cos
        F[i, 1::2] = cum_sin
    return F


def build_F_zeta(x_grid: np.ndarray, gamma_grid: np.ndarray) -> np.ndarray:
    K = len(gamma_grid)
    F = np.zeros((len(x_grid), 2*K))
    for i, x in enumerate(x_grid):
        for k, g in enumerate(gamma_grid):
            v = li_complex(float(x), complex(0.5, g))
            F[i, 2*k] = v.real
            F[i, 2*k+1] = v.imag
    return F


# -------- Regression --------

def kfold_cv_residual(F, y, k_folds=5, seed=42):
    M = len(y)
    rng = np.random.default_rng(seed)
    perm = rng.permutation(M)
    fs = M // k_folds
    oos = np.zeros(M)
    one = np.ones((M,1))
    F1 = np.hstack([one, F])
    rms_train = []
    rms_test = []
    for k in range(k_folds):
        ti = perm[k*fs : (k+1)*fs if k < k_folds-1 else M]
        tr = np.setdiff1d(np.arange(M), ti)
        c, *_ = np.linalg.lstsq(F1[tr], y[tr], rcond=None)
        oos[ti] = y[ti] - F1[ti] @ c
        rms_train.append(np.sqrt(np.mean((y[tr] - F1[tr] @ c)**2)))
        rms_test.append(np.sqrt(np.mean(oos[ti]**2)))
    return oos, {"rms_oos": float(np.sqrt(np.mean(oos**2))),
                 "rms_train_per_fold": rms_train,
                 "rms_test_per_fold": rms_test}


def in_sample_residual(F, y):
    M = len(y); one = np.ones((M,1))
    F1 = np.hstack([one, F])
    c, *_ = np.linalg.lstsq(F1, y, rcond=None)
    r = y - F1 @ c
    return r, float(np.sqrt(np.mean(r**2)))


def fit_residual_scaling(x_grid, residual):
    abs_r = np.abs(residual)
    mask = abs_r > 1e-12
    if np.sum(mask) < 5:
        return float('nan'), float('nan')
    logx = np.log(x_grid[mask]); logr = np.log(abs_r[mask])
    A = np.vstack([logx, np.ones_like(logx)]).T
    sol, *_ = np.linalg.lstsq(A, logr, rcond=None)
    beta, c0 = sol
    yhat = A @ sol
    ss_res = np.sum((logr - yhat)**2); ss_tot = np.sum((logr - logr.mean())**2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float('nan')
    return float(beta), float(r2)


# -------- SVD diagnostics --------

def svd_diagnostics(F):
    """Return singular values and effective rank (# above 1% of max)."""
    s = np.linalg.svd(F, compute_uv=False)
    eff_rank = int(np.sum(s > 0.01 * s.max()))
    return {"singular_values": s.tolist(), "max_sv": float(s.max()),
            "min_sv": float(s.min()), "condition_number": float(s.max()/max(s.min(), 1e-30)),
            "effective_rank_1pct": eff_rank, "n_features": int(F.shape[1])}


# -------- Main --------

def run(N_tau=10000, K_tau=8, K_zeta=8, M_anchors=200,
        n_random_iid_seeds=10, n_random_mult_seeds=10,
        cv_folds=5, out_path='b2_main.json'):
    print(f"== B2 Hecke vs zeta-zero CV regression ==")
    print(f"  N_tau={N_tau}, K_tau={K_tau}, K_zeta={K_zeta}, M={M_anchors}, "
          f"CV={cv_folds}, iid_seeds={n_random_iid_seeds}, "
          f"mult_seeds={n_random_mult_seeds}")

    t0 = time.time()
    print("(1) tau via eta^24...")
    tau = tau_via_eta24(N_tau)
    print(f"    {time.time()-t0:.1f}s")
    a_tau = normalised_a(tau)
    # sanity: a(p) std should be near 1 (Sato-Tate semicircle)
    primes_le_N = [p for p in range(2, N_tau+1) if isprime(p)]
    a_p_arr = np.array([a_tau[p] for p in primes_le_N])
    print(f"    a(p) for {len(primes_le_N)} primes: max|a|={np.abs(a_p_arr).max():.3f}, "
          f"std={a_p_arr.std():.3f} (Sato-Tate predicts ~1)")

    # 2) anchors
    x_lo = max(500, N_tau // 50)
    x_grid = np.unique(np.round(np.exp(np.linspace(
        np.log(x_lo), np.log(N_tau-1), M_anchors))).astype(int))
    M = len(x_grid)
    print(f"(2) M = {M} anchors in [{x_grid[0]}, {x_grid[-1]}]")

    # 3) y
    print("(3) Computing pi(x) - Li(x)...")
    t0 = time.time()
    y = np.zeros(M); pi_arr = np.zeros(M, dtype=int); li_arr = np.zeros(M)
    for i, x in enumerate(x_grid):
        pi_arr[i] = int(primepi(int(x)))
        li_arr[i] = li_real(float(x))
        y[i] = pi_arr[i] - li_arr[i]
    print(f"    {time.time()-t0:.1f}s; y rms={np.sqrt(np.mean(y**2)):.3f}")

    # 4) frequency grids
    t_grid = np.exp(np.linspace(np.log(1.0), np.log(50.0), K_tau))
    print(f"(4) t_grid = {t_grid}")

    print("(5) Riemann zeta zeros...")
    t0 = time.time()
    mpmath.mp.dps = 25
    gamma_zeta = np.array([float(mpmath.zetazero(k).imag) for k in range(1, K_zeta+1)])
    print(f"    {time.time()-t0:.1f}s; gamma_max={gamma_zeta[-1]:.2f}")

    # 5) features
    print("(6) Build F_tau, F_zeta...")
    t0 = time.time()
    F_tau = build_F_twisted(x_grid, a_tau, t_grid)
    F_zeta = build_F_zeta(x_grid, gamma_zeta)
    print(f"    {time.time()-t0:.1f}s; F_tau {F_tau.shape}, F_zeta {F_zeta.shape}")

    print("(7) Build F_random_iid (matched-cardinality, no multiplicative)...")
    F_iid = []
    for s in range(n_random_iid_seeds):
        a_iid = make_random_iid_satotate(N_tau, seed=s)
        F_iid.append(build_F_twisted(x_grid, a_iid, t_grid))

    print(f"(8) Build F_random_mult ({n_random_mult_seeds} matched-multiplicative)...")
    F_mult = []
    t0 = time.time()
    for s in range(n_random_mult_seeds):
        a_m = make_random_multiplicative(N_tau, seed=1000+s)
        F_mult.append(build_F_twisted(x_grid, a_m, t_grid))
    print(f"    {time.time()-t0:.1f}s")

    # 6) Regressions: in-sample + CV
    print("(9) Running OLS+CV...")
    r_tau_is, rms_tau_is = in_sample_residual(F_tau, y)
    r_zeta_is, rms_zeta_is = in_sample_residual(F_zeta, y)
    r_both_is, rms_both_is = in_sample_residual(np.hstack([F_tau, F_zeta]), y)

    r_tau_oos, m_tau = kfold_cv_residual(F_tau, y, cv_folds)
    r_zeta_oos, m_zeta = kfold_cv_residual(F_zeta, y, cv_folds)
    r_both_oos, m_both = kfold_cv_residual(np.hstack([F_tau, F_zeta]), y, cv_folds)

    rms_iid_oos = []; r_iid_oos = []
    for F_r in F_iid:
        r, m = kfold_cv_residual(F_r, y, cv_folds)
        rms_iid_oos.append(m["rms_oos"]); r_iid_oos.append(r)
    rms_mult_oos = []; r_mult_oos = []
    for F_r in F_mult:
        r, m = kfold_cv_residual(F_r, y, cv_folds)
        rms_mult_oos.append(m["rms_oos"]); r_mult_oos.append(r)

    # 7) scaling
    beta_y, R2_y = fit_residual_scaling(x_grid, y)
    beta_tau_oos, R2_tau_oos = fit_residual_scaling(x_grid, r_tau_oos)
    beta_zeta_oos, R2_zeta_oos = fit_residual_scaling(x_grid, r_zeta_oos)
    beta_both_oos, R2_both_oos = fit_residual_scaling(x_grid, r_both_oos)
    betas_iid = [fit_residual_scaling(x_grid, r)[0] for r in r_iid_oos]
    betas_mult = [fit_residual_scaling(x_grid, r)[0] for r in r_mult_oos]

    # 8) SVD
    svd_tau = svd_diagnostics(F_tau)
    svd_zeta = svd_diagnostics(F_zeta)
    svd_iid = [svd_diagnostics(F_r) for F_r in F_iid]
    svd_mult = [svd_diagnostics(F_r) for F_r in F_mult]

    # 9) Z-scores: rms and beta of F_tau vs each null
    iid_arr = np.array(rms_iid_oos)
    mult_arr = np.array(rms_mult_oos)
    z_iid = (m_tau["rms_oos"] - iid_arr.mean()) / max(iid_arr.std(), 1e-9)
    z_mult = (m_tau["rms_oos"] - mult_arr.mean()) / max(mult_arr.std(), 1e-9)
    z_iid_beta = (beta_tau_oos - np.mean(betas_iid)) / max(np.std(betas_iid), 1e-9)
    z_mult_beta = (beta_tau_oos - np.mean(betas_mult)) / max(np.std(betas_mult), 1e-9)

    # Print
    print()
    print("=== RESULTS ===")
    print(f"  Baseline:     rms|y|         = {np.sqrt(np.mean(y**2)):.3f}, beta_y={beta_y:.3f}")
    print(f"  F_tau:        rms_is = {rms_tau_is:.3f}, rms_oos = {m_tau['rms_oos']:.3f}, "
          f"beta_oos={beta_tau_oos:.3f}, eff_rank={svd_tau['effective_rank_1pct']}/{svd_tau['n_features']}")
    print(f"  F_zeta:       rms_is = {rms_zeta_is:.3f}, rms_oos = {m_zeta['rms_oos']:.3f}, "
          f"beta_oos={beta_zeta_oos:.3f}, eff_rank={svd_zeta['effective_rank_1pct']}/{svd_zeta['n_features']}")
    print(f"  F_tau+F_zeta: rms_is = {rms_both_is:.3f}, rms_oos = {m_both['rms_oos']:.3f}, "
          f"beta_oos={beta_both_oos:.3f}")
    print(f"  F_iid:        rms_oos = {iid_arr.mean():.3f} +/- {iid_arr.std():.3f}, "
          f"beta_oos = {np.mean(betas_iid):.3f} +/- {np.std(betas_iid):.3f}, "
          f"mean eff_rank={np.mean([s['effective_rank_1pct'] for s in svd_iid]):.1f}")
    print(f"  F_mult:       rms_oos = {mult_arr.mean():.3f} +/- {mult_arr.std():.3f}, "
          f"beta_oos = {np.mean(betas_mult):.3f} +/- {np.std(betas_mult):.3f}, "
          f"mean eff_rank={np.mean([s['effective_rank_1pct'] for s in svd_mult]):.1f}")
    print(f"  Z(F_tau vs F_iid)  rms = {z_iid:.2f}, beta = {z_iid_beta:.2f}")
    print(f"  Z(F_tau vs F_mult) rms = {z_mult:.2f}, beta = {z_mult_beta:.2f}")

    out = {
        "config": {"N_tau": N_tau, "K_tau": K_tau, "K_zeta": K_zeta,
                   "M_anchors": M, "x_lo": int(x_grid[0]), "x_hi": int(x_grid[-1]),
                   "cv_folds": cv_folds, "n_iid_seeds": n_random_iid_seeds,
                   "n_mult_seeds": n_random_mult_seeds,
                   "n_features_tau": int(F_tau.shape[1]),
                   "n_features_zeta": int(F_zeta.shape[1])},
        "x_grid": x_grid.tolist(),
        "pi_x": pi_arr.tolist(),
        "li_x": li_arr.tolist(),
        "y": y.tolist(),
        "t_grid_tau": t_grid.tolist(),
        "gamma_zeta": gamma_zeta.tolist(),
        "rms_in_sample": {"y": float(np.sqrt(np.mean(y**2))), "tau": rms_tau_is,
                          "zeta": rms_zeta_is, "both": rms_both_is},
        "rms_oos": {"tau": m_tau["rms_oos"], "zeta": m_zeta["rms_oos"],
                    "both": m_both["rms_oos"],
                    "iid_seeds": rms_iid_oos, "iid_mean": float(iid_arr.mean()),
                    "iid_std": float(iid_arr.std()),
                    "mult_seeds": rms_mult_oos, "mult_mean": float(mult_arr.mean()),
                    "mult_std": float(mult_arr.std())},
        "beta_oos": {"y": beta_y, "tau": beta_tau_oos, "zeta": beta_zeta_oos,
                     "both": beta_both_oos,
                     "iid_seeds": betas_iid, "iid_mean": float(np.mean(betas_iid)),
                     "iid_std": float(np.std(betas_iid)),
                     "mult_seeds": betas_mult, "mult_mean": float(np.mean(betas_mult)),
                     "mult_std": float(np.std(betas_mult))},
        "R2_oos": {"y": R2_y, "tau": R2_tau_oos, "zeta": R2_zeta_oos,
                   "both": R2_both_oos},
        "z_scores": {"rms_tau_vs_iid": float(z_iid),
                     "rms_tau_vs_mult": float(z_mult),
                     "beta_tau_vs_iid": float(z_iid_beta),
                     "beta_tau_vs_mult": float(z_mult_beta)},
        "svd": {"tau": svd_tau, "zeta": svd_zeta,
                "iid_seeds": svd_iid, "mult_seeds": svd_mult},
        "residuals_oos": {"tau": r_tau_oos.tolist(),
                          "zeta": r_zeta_oos.tolist(),
                          "both": r_both_oos.tolist()},
    }
    Path(out_path).write_text(json.dumps(out, indent=1))
    print(f"\nSaved {out_path}")
    return out


if __name__ == "__main__":
    import sys
    N_tau = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
    K_tau = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    K_zeta = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    M = int(sys.argv[4]) if len(sys.argv) > 4 else 200
    iid_s = int(sys.argv[5]) if len(sys.argv) > 5 else 10
    mult_s = int(sys.argv[6]) if len(sys.argv) > 6 else 10
    out = sys.argv[7] if len(sys.argv) > 7 else "b2_main.json"
    run(N_tau=N_tau, K_tau=K_tau, K_zeta=K_zeta, M_anchors=M,
        n_random_iid_seeds=iid_s, n_random_mult_seeds=mult_s,
        cv_folds=5, out_path=out)
