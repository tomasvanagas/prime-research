"""Parity-matched control for the Lyapunov-exponent comparison.

The base run found gamma_prime(E) deviating from gamma_bern_mean(E) by up
to 88 sigma at N = 10^5. PRIMES are concentrated on ODD indices (every
prime > 2 is odd). A naive Bernoulli baseline at matched density rho =
pi(N)/N puts mass on both parities. The deviation may be entirely due
to this trivial parity confound.

Controls (all at matched count = pi(N)):
  C0  Bernoulli(rho) on all indices.                    [base run, deviates 88 sigma]
  C1  Random subset of ODD indices of size pi(N) - 1, plus V(2) = 1.
      Same parity profile as chi_P. Random multiplicative structure.
  C2  Random subset of ODD indices of size pi(N) - 1, plus V(2) = 1,
      with adjacent-pair correlation matched to twin-prime density.
      (Stretch: not run by default.)
  C3  V(n) = chi_P(n) but shuffled within odd indices (preserves parity
      profile and exact count, destroys ordering).

Tests:
  T1  gamma_prime vs C1 across n_seeds runs of C1: residual z-score?
      If |z| < 3 everywhere -> the 88 sigma is a parity artefact.
      If |z| >> 3 -> deeper arithmetic structure.
  T2  Compare gamma_C1_mean to gamma_C0_mean: by HOW MUCH does parity
      restriction shift the Lyapunov curve? Quantifies the parity
      confound size.

Usage:
  python parity_control.py --N 100000 --seeds 50
"""
import argparse
import json
import os
import time

import numpy as np

from anderson_localisation_chi_p import (
    chi_p_array,
    lyapunov_vec,
)


def random_odd_subset_potential(N, n_primes_minus_1, rng, amp=1.0):
    """Random subset of odd indices of given size; plus V(2) = 1."""
    odd_idx = np.arange(3, N + 1, 2)  # odd indices >= 3 (1-indexed)
    chosen = rng.choice(len(odd_idx), size=n_primes_minus_1, replace=False)
    V = np.zeros(N, dtype=np.float64)
    # Convert 1-indexed positions to 0-indexed array positions.
    V[odd_idx[chosen] - 1] = amp
    V[2 - 1] = amp  # the prime 2
    return V


def shuffled_chi_p_within_odds(N, rng, amp=1.0):
    """chi_P shuffled within odd indices: same parity profile, exact count,
    destroys ordering. V(2)=1 fixed."""
    chi = chi_p_array(N) * amp
    odd_pos = np.arange(2, N, 2)  # 0-indexed positions of odd integers >= 3
    odd_vals = chi[odd_pos].copy()
    rng.shuffle(odd_vals)
    chi[odd_pos] = odd_vals
    return chi


def run(N, n_seeds, n_energies, E_lo=-1.95, E_hi=2.95, amp=1.0, out_path=None):
    energies = np.linspace(E_lo, E_hi, n_energies)
    print(f"[parity] N={N}, n_seeds={n_seeds}, n_energies={n_energies}, amp={amp}")

    t0 = time.time()
    V_prime = amp * chi_p_array(N)
    n_primes = int((V_prime > 0).sum())
    print(f"[parity] pi(N) = {n_primes}, density = {n_primes/N:.6f}")
    g_prime = lyapunov_vec(V_prime, energies)
    print(f"[parity] gamma_prime in {time.time()-t0:.2f}s")

    # C1: random odd-subset potential (matched count, matched parity).
    c1_curves = np.zeros((n_seeds, n_energies))
    for s in range(n_seeds):
        rng = np.random.default_rng(s * 23 + 5)
        V = random_odd_subset_potential(N, n_primes - 1, rng, amp=amp)
        c1_curves[s] = lyapunov_vec(V, energies)
        if (s + 1) % 10 == 0:
            print(f"[parity] C1 seed {s+1}/{n_seeds} at {time.time()-t0:.1f}s")
    c1_mean = c1_curves.mean(axis=0)
    c1_std = c1_curves.std(axis=0, ddof=1)
    z_c1 = (g_prime - c1_mean) / np.maximum(c1_std, 1e-12)
    max_z_c1 = float(np.abs(z_c1).max())
    arg_z_c1 = int(np.abs(z_c1).argmax())
    print(f"[parity] vs C1 (parity-matched): max |z| = {max_z_c1:.3f}"
          f" at E = {energies[arg_z_c1]:.4f}")

    # C3: chi_P shuffled within odd indices (exact parity profile, exact
    # multiplicities, destroys ordering).
    c3_curves = np.zeros((n_seeds, n_energies))
    for s in range(n_seeds):
        rng = np.random.default_rng(s * 29 + 11)
        V = shuffled_chi_p_within_odds(N, rng, amp=amp)
        c3_curves[s] = lyapunov_vec(V, energies)
    c3_mean = c3_curves.mean(axis=0)
    c3_std = c3_curves.std(axis=0, ddof=1)
    z_c3 = (g_prime - c3_mean) / np.maximum(c3_std, 1e-12)
    max_z_c3 = float(np.abs(z_c3).max())
    arg_z_c3 = int(np.abs(z_c3).argmax())
    print(f"[parity] vs C3 (shuffled within odds): max |z| = {max_z_c3:.3f}"
          f" at E = {energies[arg_z_c3]:.4f}")

    # Compare C1 to a fully-Bernoulli baseline (C0).
    c0_curves = np.zeros((min(n_seeds, 20), n_energies))
    rho = n_primes / N
    for s in range(c0_curves.shape[0]):
        rng = np.random.default_rng(s * 31 + 17)
        V = (rng.random(N) < rho).astype(np.float64) * amp
        c0_curves[s] = lyapunov_vec(V, energies)
    c0_mean = c0_curves.mean(axis=0)

    out = {
        "N": N,
        "n_seeds": n_seeds,
        "amp": amp,
        "n_primes": n_primes,
        "energies": energies.tolist(),
        "gamma_prime": g_prime.tolist(),
        "C1_mean": c1_mean.tolist(),
        "C1_std": c1_std.tolist(),
        "C1_z": z_c1.tolist(),
        "C1_max_abs_z": max_z_c1,
        "C1_argmax_E": float(energies[arg_z_c1]),
        "C3_mean": c3_mean.tolist(),
        "C3_std": c3_std.tolist(),
        "C3_z": z_c3.tolist(),
        "C3_max_abs_z": max_z_c3,
        "C3_argmax_E": float(energies[arg_z_c3]),
        "C0_mean": c0_mean.tolist(),  # for reference: full-Bernoulli baseline
        "elapsed": time.time() - t0,
    }
    if out_path:
        with open(out_path, "w") as f:
            json.dump(out, f, indent=2)
        print(f"[parity] wrote {out_path}")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=100000)
    ap.add_argument("--seeds", type=int, default=50)
    ap.add_argument("--energies", type=int, default=51)
    ap.add_argument("--amp", type=float, default=1.0)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()
    out_path = args.out
    if out_path is None:
        out_path = os.path.join(
            os.path.dirname(__file__),
            f"parity_N{args.N}_s{args.seeds}_e{args.energies}_a{args.amp}.json",
        )
    run(N=args.N, n_seeds=args.seeds, n_energies=args.energies,
        amp=args.amp, out_path=out_path)
