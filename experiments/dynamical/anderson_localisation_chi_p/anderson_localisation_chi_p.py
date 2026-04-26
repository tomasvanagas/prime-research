"""Anderson localisation: Lyapunov exponent of the discrete 1D Schrödinger
operator with prime-indicator potential V(n) = chi_P(n).

ATTACK_VECTORS.md §C4. Cross-domain technique: Anderson localization
(Aizenman-Warzel 2015; Furstenberg-Kifer for products of SL(2,R)
matrices). Spectral statistic NOT in the project's 35+ pseudorandomness
battery (which are all local: correlations, entropy, LFSR, tensor rank).
The Lyapunov exponent gamma(E) is global and spectral.

Setup:
  H psi(n) = -psi(n+1) - psi(n-1) + V(n) psi(n) = E psi(n)
  Recurrence: psi(n+1) = (V(n) - E) psi(n) - psi(n-1)
  Transfer matrix: T_n = [[V(n) - E, -1], [1, 0]],  det T_n = +1.
  Lyapunov:   gamma(E) = lim_{N -> infty} (1/N) log ||T_N ... T_1||.

Free spectrum (V == 0): plane waves psi(n) = e^{ikn} with E = -2 cos(k),
so the band is E in [-2, 2]. For sparse V in {0, 1} at density rho, perturbation
theory (Pastur-Figotin) gives gamma(E) ~ (rho * sigma^2) / (8 sin^2(k))
inside the band.

Usage:
  python anderson_localisation_chi_p.py --N 100000 --seeds 50 --energies 51
  python anderson_localisation_chi_p.py --quick       (fast smoke test)

Output: JSON with per-energy gamma curves for chi_P, Bernoulli ensemble,
and Liouville (1 - lambda)/2.

Falsification: if gamma_prime(E) lies inside the Bernoulli noise band
([mean - 3*std, mean + 3*std] over 50 seeds) at every E, chi_P is
spectrally indistinguishable from Bernoulli at this scale: 36th
pseudorandomness measure with no deviation.
"""
import argparse
import json
import os
import sys
import time

import numpy as np


def chi_p_array(N):
    """V(n) = 1 if n is prime else 0, for n = 1, 2, ..., N."""
    sieve = np.ones(N + 1, dtype=np.uint8)
    sieve[0] = 0
    if N >= 1:
        sieve[1] = 0
    for p in range(2, int(N ** 0.5) + 1):
        if sieve[p]:
            sieve[p * p :: p] = 0
    return sieve[1 : N + 1].astype(np.float64)


def liouville_array(N):
    """V(n) = (1 - lambda(n))/2 in {0, 1}, density 1/2 conjecturally.

    lambda(n) = (-1)^Omega(n), Omega = total prime-factor count w/ multiplicity.
    """
    # Compute Omega(n) via sieve: smallest-prime-factor sieve.
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
    lam = np.where(omega == 0, 1, np.where(omega % 2 == 0, 1, -1))
    return ((1 - lam[1:]) // 2).astype(np.float64)


def lyapunov_vec(V, energies, renorm_period=32):
    """Vectorized Lyapunov over an energy grid.

    Returns gamma(E) for each E. Uses periodic L2 renormalization of the
    state pair (psi_curr, psi_prev) to avoid overflow.
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


def run(N, n_seeds, n_energies, E_lo=-1.95, E_hi=2.95, out_path=None,
        amp=1.0, do_liouville=True):
    """Main experiment.

    amp scales the potential V -> amp * V to test perturbation-theory
    dependence (gamma scales as amp^2 within the band).
    """
    energies = np.linspace(E_lo, E_hi, n_energies)
    print(f"[run] N={N}, n_seeds={n_seeds}, n_energies={n_energies}, amp={amp}")
    print(f"[run] energy grid: [{E_lo}, {E_hi}]")

    t0 = time.time()
    V_prime = amp * chi_p_array(N)
    rho_prime = float((V_prime > 0).mean())
    print(f"[run] chi_P density = {rho_prime:.6f} (pi(N)/N)")
    g_prime = lyapunov_vec(V_prime, energies)
    print(f"[run] gamma_prime computed in {time.time() - t0:.2f}s")

    # Bernoulli ensemble at matched density.
    bern_curves = np.zeros((n_seeds, n_energies))
    rho = rho_prime / amp  # un-amplitude the density
    for seed in range(n_seeds):
        rng = np.random.default_rng(seed * 17 + 1)
        V_b = (rng.random(N) < rho).astype(np.float64) * amp
        bern_curves[seed] = lyapunov_vec(V_b, energies)
        if (seed + 1) % 5 == 0:
            print(f"[run] Bernoulli seed {seed+1}/{n_seeds} done at {time.time()-t0:.1f}s")
    g_bern_mean = bern_curves.mean(axis=0)
    g_bern_std = bern_curves.std(axis=0, ddof=1)

    # Per-energy z-scores: how many sigma is chi_P from Bernoulli mean?
    z_scores = (g_prime - g_bern_mean) / np.maximum(g_bern_std, 1e-12)
    max_abs_z = float(np.abs(z_scores).max())
    argmax_z = int(np.abs(z_scores).argmax())
    print(f"[run] max |z| = {max_abs_z:.3f} at E = {energies[argmax_z]:.4f}")

    out = {
        "N": N,
        "n_seeds": n_seeds,
        "n_energies": n_energies,
        "amp": amp,
        "energies": energies.tolist(),
        "rho_prime": rho_prime,
        "gamma_prime": g_prime.tolist(),
        "gamma_bern_mean": g_bern_mean.tolist(),
        "gamma_bern_std": g_bern_std.tolist(),
        "gamma_bern_curves": bern_curves.tolist(),
        "z_scores": z_scores.tolist(),
        "max_abs_z": max_abs_z,
        "argmax_z_energy": float(energies[argmax_z]),
    }

    # Pastur-Figotin perturbation prediction inside band
    # gamma(E) ~ (rho * (1-rho) * amp^2) / (8 sin^2(k)),  E = -2 cos(k)
    # rho * (1 - rho) is the variance of Bernoulli(rho) re V in {0, amp};
    # actually V variance = amp^2 * rho * (1 - rho).
    pf_pred = []
    for E in energies:
        if -2.0 < E < 2.0:
            k = np.arccos(-E / 2.0)
            pf = (rho * (1 - rho) * amp ** 2) / (8 * (np.sin(k) ** 2))
        else:
            pf = float("nan")
        pf_pred.append(pf)
    out["gamma_pastur_figotin_pred"] = pf_pred

    if do_liouville and N <= 200000:
        # Liouville is also density ~ 1/2 in {0, 1}; compare to its own
        # Bernoulli(1/2) baseline.
        print(f"[run] computing Liouville at N={N}...")
        V_liou = amp * liouville_array(N)
        rho_liou = float((V_liou > 0).mean())
        g_liou = lyapunov_vec(V_liou, energies)
        # Bernoulli(1/2) reference, fewer seeds
        n_seeds_liou = max(10, n_seeds // 2)
        bern12 = np.zeros((n_seeds_liou, n_energies))
        for seed in range(n_seeds_liou):
            rng = np.random.default_rng(seed * 19 + 7)
            V_b = (rng.random(N) < rho_liou).astype(np.float64) * amp
            bern12[seed] = lyapunov_vec(V_b, energies)
        out["liouville"] = {
            "rho": rho_liou,
            "gamma": g_liou.tolist(),
            "gamma_bern_mean": bern12.mean(axis=0).tolist(),
            "gamma_bern_std": bern12.std(axis=0, ddof=1).tolist(),
            "z_scores": ((g_liou - bern12.mean(axis=0))
                         / np.maximum(bern12.std(axis=0, ddof=1), 1e-12)).tolist(),
        }
        liou_z = np.abs(out["liouville"]["z_scores"]).max()
        print(f"[run] Liouville max |z| = {liou_z:.3f}")

    elapsed = time.time() - t0
    out["elapsed_seconds"] = elapsed
    print(f"[run] total elapsed: {elapsed:.1f}s")

    if out_path:
        with open(out_path, "w") as f:
            json.dump(out, f, indent=2)
        print(f"[run] wrote {out_path}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=100000)
    ap.add_argument("--seeds", type=int, default=50)
    ap.add_argument("--energies", type=int, default=51)
    ap.add_argument("--E_lo", type=float, default=-1.95)
    ap.add_argument("--E_hi", type=float, default=2.95)
    ap.add_argument("--amp", type=float, default=1.0)
    ap.add_argument("--no-liouville", action="store_true")
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
            f"results_N{args.N}_s{args.seeds}_e{args.energies}_a{args.amp}.json",
        )

    run(N=args.N, n_seeds=args.seeds, n_energies=args.energies,
        E_lo=args.E_lo, E_hi=args.E_hi, out_path=out_path, amp=args.amp,
        do_liouville=not args.no_liouville)


if __name__ == "__main__":
    main()
