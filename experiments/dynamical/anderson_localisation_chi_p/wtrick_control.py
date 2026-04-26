"""W-tricked control. Reduces chi_P vs random comparison after removing
small-modulus residue-class structure.

The primary run found chi_P gamma deviates from parity-matched (C1)
random by 33 sigma at E = 1.088 ~ -2 cos(2 pi / 3) = +1, the MOD-3
RESONANCE. Hypothesis: the deviation is just the mod-6 sieving
(primes >3 are coprime to 6).

W-trick (Green-Tao):
  W = product of small primes (W = 2*3 = 6, or W = 2*3*5 = 30, etc.)
  Restrict V to support {n : gcd(n, W) = 1}.
  Match the count to chi_P's count.

Controls:
  CW6  random subset of {n : gcd(n, 6) = 1} ∩ [1, N], size pi(N) - 2.
  CW30 random subset of {n : gcd(n, 30) = 1} ∩ [1, N], size pi(N) - 3.
  CW210 random subset of {n : gcd(n, 210) = 1} ∩ [1, N], size pi(N) - 4.

If gamma_prime matches CW210 within sample noise -> the localization
deviation is just the small-mod sieving of primes.

If gamma_prime STILL deviates from CW210 -> deeper structure (twin
primes? prime k-tuples? Hardy-Littlewood density?).
"""
import argparse
import json
import os
import time

import numpy as np

from anderson_localisation_chi_p import chi_p_array, lyapunov_vec


def coprime_mask(N, W):
    """Boolean array: mask[i-1] = True iff gcd(i, W) = 1, for i in [1, N]."""
    n = np.arange(1, N + 1, dtype=np.int64)
    mask = np.ones(N, dtype=bool)
    for p in [2, 3, 5, 7, 11, 13]:
        if W % p == 0:
            mask &= (n % p != 0)
    return mask


def w_tricked_random(N, W, n_target, rng, amp=1.0):
    """Random subset of {n : gcd(n, W) = 1} of size n_target. Plus the
    primes p <= W' that are removed (they sit on positions divisible by
    themselves; we add them back only if we want to match chi_P's full
    count). For simplicity here we MATCH the coprime portion of chi_P
    (i.e., primes > max prime factor of W).
    """
    mask = coprime_mask(N, W)
    eligible = np.where(mask)[0]  # 0-indexed positions of integers coprime to W
    if n_target > len(eligible):
        raise ValueError(f"n_target={n_target} > eligible={len(eligible)}")
    chosen = rng.choice(len(eligible), size=n_target, replace=False)
    V = np.zeros(N, dtype=np.float64)
    V[eligible[chosen]] = amp
    return V


def run(N, n_seeds, n_energies, E_lo=-1.95, E_hi=2.95, amp=1.0, out_path=None):
    energies = np.linspace(E_lo, E_hi, n_energies)
    print(f"[wtrick] N={N}, n_seeds={n_seeds}, n_energies={n_energies}")

    t0 = time.time()
    V_prime = amp * chi_p_array(N)
    n_primes = int((V_prime > 0).sum())
    # chi_P contribution at coprime-to-W positions (= primes > max(prime
    # factors of W)):
    coprime_to_6 = ((V_prime > 0) & coprime_mask(N, 6)).sum()
    coprime_to_30 = ((V_prime > 0) & coprime_mask(N, 30)).sum()
    coprime_to_210 = ((V_prime > 0) & coprime_mask(N, 210)).sum()
    print(f"[wtrick] pi(N) = {n_primes}; coprime6={coprime_to_6} (=pi-2),"
          f" coprime30={coprime_to_30}, coprime210={coprime_to_210}")
    g_prime = lyapunov_vec(V_prime, energies)
    print(f"[wtrick] gamma_prime in {time.time()-t0:.2f}s")

    results = {
        "N": N, "n_seeds": n_seeds, "amp": amp,
        "n_primes": n_primes,
        "energies": energies.tolist(),
        "gamma_prime": g_prime.tolist(),
    }

    coprime_to_2310 = ((V_prime > 0) & coprime_mask(N, 2310)).sum()
    print(f"[wtrick] coprime2310={coprime_to_2310}")
    for W, n_match in [(6, coprime_to_6), (30, coprime_to_30),
                       (210, coprime_to_210), (2310, coprime_to_2310)]:
        curves = np.zeros((n_seeds, n_energies))
        for s in range(n_seeds):
            rng = np.random.default_rng(s * 37 + W)
            V = w_tricked_random(N, W, n_match, rng, amp=amp)
            curves[s] = lyapunov_vec(V, energies)
        mean = curves.mean(axis=0)
        std = curves.std(axis=0, ddof=1)
        # Compare to gamma_prime restricted to its coprime-to-W subset
        # (since chi_P contains primes 2, 3, 5, 7 etc that are NOT coprime
        # to W; for fair comparison we should add them back to the random
        # control). Here we add them as fixed delta.
        # Add the small primes p | W to V to make a "chi_P*" matching:
        fixed = np.zeros(N, dtype=np.float64)
        small = [p for p in [2, 3, 5, 7, 11, 13] if W % p == 0]
        for p in small:
            if p <= N:
                fixed[p - 1] = amp
        curves_fp = np.zeros_like(curves)
        for s in range(n_seeds):
            rng = np.random.default_rng(s * 37 + W)
            V = w_tricked_random(N, W, n_match, rng, amp=amp) + fixed
            curves_fp[s] = lyapunov_vec(V, energies)
        mean_fp = curves_fp.mean(axis=0)
        std_fp = curves_fp.std(axis=0, ddof=1)
        z_fp = (g_prime - mean_fp) / np.maximum(std_fp, 1e-12)
        max_z = float(np.abs(z_fp).max())
        arg_z = int(np.abs(z_fp).argmax())
        print(f"[wtrick] W={W}: max |z|(prime vs CW{W}+small) = {max_z:.3f}"
              f" at E = {energies[arg_z]:.4f}")
        results[f"CW{W}"] = {
            "n_match": int(n_match),
            "small_primes_added": small,
            "mean": mean_fp.tolist(),
            "std": std_fp.tolist(),
            "z": z_fp.tolist(),
            "max_abs_z": max_z,
            "argmax_E": float(energies[arg_z]),
        }

    results["elapsed"] = time.time() - t0
    if out_path:
        with open(out_path, "w") as f:
            json.dump(results, f, indent=2)
        print(f"[wtrick] wrote {out_path}")
    return results


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=100000)
    ap.add_argument("--seeds", type=int, default=30)
    ap.add_argument("--energies", type=int, default=51)
    ap.add_argument("--amp", type=float, default=1.0)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()
    out_path = args.out or os.path.join(
        os.path.dirname(__file__),
        f"wtrick_N{args.N}_s{args.seeds}_e{args.energies}.json",
    )
    run(N=args.N, n_seeds=args.seeds, n_energies=args.energies,
        amp=args.amp, out_path=out_path)
