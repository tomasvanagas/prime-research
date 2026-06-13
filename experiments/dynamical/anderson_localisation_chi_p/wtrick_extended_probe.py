"""S416 re-verify-closure probe for E2.14 (Anderson Lyapunov of chi_P).

Adversarial extension of `wtrick_control.py` testing **Falsifier #1** from
the S88 closure (`anderson_localisation_chi_p_results.md` §"What Would
Falsify This"):

    W-trick saturation. A residual deviation gamma_prime - gamma_CW of
    >> 5 sigma persists at any W = primorial(k) no matter how large k,
    and cannot be explained by a k-tuple Hardy-Littlewood constant.

The S88 cascade reached W=2310 (= 2*3*5*7*11) with max |z|=4.0 at
N=2*10^5, 30 seeds. HL prediction: continuing the cascade to
W=30030 (= 2*3*5*7*11*13) and W=510510 (= 2*3*5*7*11*13*17) the
deviation should keep decaying as k-tuple HL singular series corrections
shrink. If max |z| at W=30030 saturates at >= 4 sigma instead of
decaying, that is a non-HL signal — **the missed angle the closure
would have failed to anticipate**.

Run with:
    python wtrick_extended_probe.py --N 200000 --seeds 30 --energies 31
"""
import argparse
import json
import os
import time

import numpy as np

from anderson_localisation_chi_p import chi_p_array, lyapunov_vec


PRIMORIAL_PRIMES = {
    2: [2],
    6: [2, 3],
    30: [2, 3, 5],
    210: [2, 3, 5, 7],
    2310: [2, 3, 5, 7, 11],
    30030: [2, 3, 5, 7, 11, 13],
    510510: [2, 3, 5, 7, 11, 13, 17],
}


def coprime_mask(N, W):
    """Boolean array: mask[i-1] = True iff gcd(i, W) = 1, for i in [1, N]."""
    n = np.arange(1, N + 1, dtype=np.int64)
    mask = np.ones(N, dtype=bool)
    for p in PRIMORIAL_PRIMES[W]:
        mask &= (n % p != 0)
    return mask


def w_tricked_random(N, W, n_target, rng, amp=1.0):
    mask = coprime_mask(N, W)
    eligible = np.where(mask)[0]
    if n_target > len(eligible):
        raise ValueError(f"n_target={n_target} > eligible={len(eligible)}")
    chosen = rng.choice(len(eligible), size=n_target, replace=False)
    V = np.zeros(N, dtype=np.float64)
    V[eligible[chosen]] = amp
    return V


def run(N, n_seeds, n_energies, W_list, E_lo=-1.95, E_hi=2.95, amp=1.0,
        out_path=None):
    energies = np.linspace(E_lo, E_hi, n_energies)
    print(f"[probe] N={N}, n_seeds={n_seeds}, n_energies={n_energies}, "
          f"W_list={W_list}")

    t0 = time.time()
    V_prime = amp * chi_p_array(N)
    n_primes = int((V_prime > 0).sum())
    print(f"[probe] pi(N) = {n_primes}, density = {n_primes/N:.5f}")
    g_prime = lyapunov_vec(V_prime, energies)
    print(f"[probe] gamma_prime computed in {time.time()-t0:.2f}s")

    results = {
        "N": N,
        "n_seeds": n_seeds,
        "amp": amp,
        "n_primes": n_primes,
        "energies": energies.tolist(),
        "gamma_prime": g_prime.tolist(),
        "W_list": W_list,
    }

    for W in W_list:
        small = PRIMORIAL_PRIMES[W]
        # Match coprime-to-W portion of chi_P (= primes > max(small))
        coprime_count = int(((V_prime > 0) & coprime_mask(N, W)).sum())
        max_small = max(small)
        n_small_in_chi = sum(1 for p in small if p <= N)
        print(f"[probe] W={W}: coprime={coprime_count}, "
              f"small_primes_added={small} ({n_small_in_chi} fit in [1, N])")

        # Build random control: random subset of {n : gcd(n,W)=1} matching
        # coprime_count, plus fixed deltas at the small primes themselves.
        fixed = np.zeros(N, dtype=np.float64)
        for p in small:
            if p <= N:
                fixed[p - 1] = amp

        curves = np.zeros((n_seeds, n_energies))
        t1 = time.time()
        for s in range(n_seeds):
            rng = np.random.default_rng(s * 41 + W)
            V_rand = w_tricked_random(N, W, coprime_count, rng, amp=amp)
            V = V_rand + fixed
            curves[s] = lyapunov_vec(V, energies)
            if (s + 1) % 5 == 0:
                print(f"[probe] W={W} seed {s+1}/{n_seeds} done at "
                      f"{time.time()-t1:.1f}s")
        mean = curves.mean(axis=0)
        std = curves.std(axis=0, ddof=1)
        z = (g_prime - mean) / np.maximum(std, 1e-12)
        max_z = float(np.abs(z).max())
        arg_z = int(np.abs(z).argmax())
        print(f"[probe] W={W}: max |z| = {max_z:.3f} at "
              f"E = {energies[arg_z]:.4f}")
        results[f"CW{W}"] = {
            "n_match": coprime_count,
            "small_primes_added": small,
            "mean": mean.tolist(),
            "std": std.tolist(),
            "z": z.tolist(),
            "max_abs_z": max_z,
            "argmax_E": float(energies[arg_z]),
        }

    results["elapsed"] = time.time() - t0
    if out_path:
        with open(out_path, "w") as f:
            json.dump(results, f, indent=2)
        print(f"[probe] wrote {out_path}")
    return results


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=200000)
    ap.add_argument("--seeds", type=int, default=30)
    ap.add_argument("--energies", type=int, default=31)
    ap.add_argument("--amp", type=float, default=1.0)
    ap.add_argument("--W", type=int, nargs="+",
                    default=[2310, 30030, 510510])
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()
    out_path = args.out or os.path.join(
        os.path.dirname(__file__),
        f"wtrick_extended_N{args.N}_s{args.seeds}_e{args.energies}.json",
    )
    run(N=args.N, n_seeds=args.seeds, n_energies=args.energies,
        W_list=args.W, amp=args.amp, out_path=out_path)
