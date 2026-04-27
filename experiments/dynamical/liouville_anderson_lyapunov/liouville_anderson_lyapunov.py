"""G1 — Liouville Anderson Lyapunov: spectral signature of Mobius pseudorandomness.

ATTACK_VECTORS.md §G1 (current highest-leverage frontier; multiplicative
regime beyond the W-trick wall). Cross-domain technique: Anderson
localisation theory (already imported S88, E2.14) + Mobius/nilsequence
orthogonality (Green-Tao arXiv:0807.1736).

Setup. Discrete 1D Schrodinger operator on Z:
    H psi(n) = -psi(n+1) - psi(n-1) + V(n) psi(n) = E psi(n)
Recurrence: psi(n+1) = (V(n) - E) psi(n) - psi(n-1).
Transfer matrix T_n(E) = [[V(n) - E, -1], [1, 0]],  det T_n = +1.
Lyapunov:  gamma(E) = lim_{N->infty} (1/N) log ||T_N ... T_1||.

Crucial difference from S88. S88 used V = chi_P in {0, 1} (density 1/log N)
or Liouville in (1-lambda)/2 in {0, 1} (density ~1/2). Both encodings
shift the energy origin and inflate the variance by a density factor.

THIS experiment uses the *centered* multiplicative encoding
    V(n) = lambda(n) in {-1, +1}     (mean -> 0 by PNT for lambda)
which has variance = 1 exactly and zero mean asymptotically. The
Pastur-Figotin perturbation prediction inside the band is then
    gamma_PF(E) = 1 / (8 sin^2 k),  E = -2 cos k.
Exactly matched by Rademacher (i.i.d. {-1, +1}). Any deviation
gamma_lambda(E) - gamma_Rademacher(E) at any energy that is sustained
across seed averaging is the spectral signature of multiplicative
structure in lambda that Rademacher does not have.

A-grade. |z| > 5 sustained at some energy E_0, NOT removed by the
multiplicative W-trick lambda_{W,b}(n) := lambda(W*n + b).

B-grade. |z| <= 3 across all energies => 38th pseudorandomness
measure at the spectral level confirming Mobius/nilsequence
orthogonality.

Usage:
  python liouville_anderson_lyapunov.py --N 100000 --seeds 50 --energies 51
  python liouville_anderson_lyapunov.py --quick
"""
import argparse
import json
import os
import time

import numpy as np


def liouville_pm1(N):
    """Return numpy array of lambda(n) in {-1, +1} for n = 1..N.

    Computed via smallest-prime-factor sieve in O(N log log N), then
    Omega(n) parity gives lambda(n) = (-1)^Omega(n).
    """
    spf = np.arange(N + 1, dtype=np.int64)
    for p in range(2, int(N ** 0.5) + 1):
        if spf[p] == p:
            for k in range(p * p, N + 1, p):
                if spf[k] == k:
                    spf[k] = p
    omega = np.zeros(N + 1, dtype=np.int64)
    for n in range(2, N + 1):
        m, o = n, 0
        while m > 1:
            p = spf[m]
            while m % p == 0:
                m //= p
                o += 1
        omega[n] = o
    omega[1] = 0  # lambda(1) = +1
    lam = np.where(omega % 2 == 0, 1.0, -1.0)
    return lam[1:].astype(np.float64)


def lyapunov_vec(V, energies, renorm_period=32):
    """Vectorised Lyapunov over an energy grid.

    Returns gamma(E) for each E with periodic L^2 renormalisation.
    """
    N = len(V)
    K = len(energies)
    psi_curr = np.ones(K)
    psi_prev = np.zeros(K)
    log_norm = np.zeros(K)
    Earr = np.asarray(energies, dtype=np.float64)
    for n in range(N):
        Vn = V[n]
        psi_next = (Vn - Earr) * psi_curr - psi_prev
        psi_prev = psi_curr
        psi_curr = psi_next
        if (n + 1) % renorm_period == 0:
            norm = np.sqrt(psi_curr * psi_curr + psi_prev * psi_prev)
            mask = norm > 0
            log_norm[mask] += np.log(norm[mask])
            psi_curr[mask] /= norm[mask]
            psi_prev[mask] /= norm[mask]
    norm = np.sqrt(psi_curr * psi_curr + psi_prev * psi_prev)
    mask = norm > 0
    log_norm[mask] += np.log(norm[mask])
    return log_norm / N


def run(N, n_seeds, n_energies, E_lo=-2.95, E_hi=2.95, out_path=None):
    """Main G1 experiment.

    V_lambda has variance 1 exactly and asymptotic mean 0, so the
    matched-variance baseline is Rademacher i.i.d. uniform on {-1, +1}.
    """
    energies = np.linspace(E_lo, E_hi, n_energies)
    print(f"[g1] N={N}, n_seeds={n_seeds}, n_energies={n_energies}")
    print(f"[g1] energy grid: [{E_lo}, {E_hi}]")

    t0 = time.time()
    V_lam = liouville_pm1(N)
    mean_lam = float(V_lam.mean())
    var_lam = float(V_lam.var())
    print(f"[g1] lambda mean = {mean_lam:.6e}, var = {var_lam:.6f}")
    g_lam = lyapunov_vec(V_lam, energies)
    print(f"[g1] gamma_lambda computed in {time.time() - t0:.2f}s")

    # Rademacher baseline: i.i.d. uniform {-1, +1}.
    bern_curves = np.zeros((n_seeds, n_energies))
    for seed in range(n_seeds):
        rng = np.random.default_rng(seed * 23 + 5)
        V_r = rng.choice([-1.0, 1.0], size=N).astype(np.float64)
        bern_curves[seed] = lyapunov_vec(V_r, energies)
        if (seed + 1) % 10 == 0:
            print(f"[g1] Rademacher seed {seed+1}/{n_seeds} done at {time.time()-t0:.1f}s")
    g_bern_mean = bern_curves.mean(axis=0)
    g_bern_std = bern_curves.std(axis=0, ddof=1)

    # Per-energy z-scores: signed deviation of lambda Lyapunov vs Rademacher.
    z_scores = (g_lam - g_bern_mean) / np.maximum(g_bern_std, 1e-12)
    max_abs_z = float(np.abs(z_scores).max())
    argmax_z = int(np.abs(z_scores).argmax())
    print(f"[g1] max |z| = {max_abs_z:.3f} at E = {energies[argmax_z]:.4f}")

    # Pastur-Figotin prediction for comparison.
    pf = []
    for E in energies:
        if -2.0 < E < 2.0:
            k = np.arccos(-E / 2.0)
            pf.append(1.0 / (8.0 * (np.sin(k) ** 2)))
        else:
            pf.append(float("nan"))

    out = {
        "N": N,
        "n_seeds": n_seeds,
        "n_energies": n_energies,
        "energies": energies.tolist(),
        "lambda_mean": mean_lam,
        "lambda_var": var_lam,
        "gamma_lambda": g_lam.tolist(),
        "gamma_rademacher_mean": g_bern_mean.tolist(),
        "gamma_rademacher_std": g_bern_std.tolist(),
        "gamma_rademacher_curves": bern_curves.tolist(),
        "z_scores": z_scores.tolist(),
        "max_abs_z": max_abs_z,
        "argmax_z_energy": float(energies[argmax_z]),
        "pastur_figotin_pred": pf,
    }

    elapsed = time.time() - t0
    out["elapsed_seconds"] = elapsed
    print(f"[g1] total elapsed: {elapsed:.1f}s")

    if out_path:
        with open(out_path, "w") as f:
            json.dump(out, f, indent=2)
        print(f"[g1] wrote {out_path}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=100000)
    ap.add_argument("--seeds", type=int, default=50)
    ap.add_argument("--energies", type=int, default=51)
    ap.add_argument("--E_lo", type=float, default=-2.95)
    ap.add_argument("--E_hi", type=float, default=2.95)
    ap.add_argument("--quick", action="store_true",
                    help="Smoke test: N=10000, 10 seeds, 21 energies.")
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    if args.quick:
        args.N = 10000
        args.seeds = 10
        args.energies = 21

    out_path = args.out
    if out_path is None:
        out_path = os.path.join(
            os.path.dirname(__file__),
            f"results_N{args.N}_s{args.seeds}_e{args.energies}.json",
        )

    run(N=args.N, n_seeds=args.seeds, n_energies=args.energies,
        E_lo=args.E_lo, E_hi=args.E_hi, out_path=out_path)


if __name__ == "__main__":
    main()
