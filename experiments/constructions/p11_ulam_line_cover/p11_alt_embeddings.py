"""
Thread 11 / Slot 2 — line cover under alternative 2D embeddings.

Slot 1 (S484) measured L_Ulam(primes, N) on the Ulam spiral and found
both primes and matched-baseline random points scale as ~sqrt(N), with
prime advantage L_p/L_r shrinking as N^{-0.24} (Ulam embedding gives
no constant-factor compression at large N).

Slot 2 asks: do other 2D embeddings give qualitatively different
scaling? Two families:

  (R) Residue-class grid:        Phi_R(n; p) = (n mod p, n // p)
       Grid is p columns wide x N/p rows tall.
       Vertical lines = arithmetic progressions a (mod p).
  (Q) Polynomial-image grid:     Phi_Q(n; p) = (n, n^2 mod p)
       Strip N wide x p tall.
       Horizontal lines = level sets of n^2 mod p (QR-classes).

For each embedding we compute the bounded-direction greedy line cover
L_Phi(primes, N) and a matched-baseline L_Phi(random, N), at
N in {10^4, 10^5, 10^6}, p in {2, 3, 5, 7, 11}.

Imports the bounded_greedy_line_cover from p11_ulam_bounded.

CLI:
  python p11_alt_embeddings.py --N 1000000 --K 5 --baseline-trials 3
  python p11_alt_embeddings.py --N 100000  --K 10 --baseline-trials 5
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import random
import sys
import time

import numpy as np

# Import line-cover machinery from slot-1 evaluator.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p11_ulam_bounded import bounded_greedy_line_cover, sieve_primes


def residue_class_coords(N: int, p: int) -> tuple[np.ndarray, np.ndarray]:
    """Phi_R(n; p) = (n mod p, n // p). Returns X, Y arrays of size N+1.
    X[0], Y[0] unused."""
    n = np.arange(N + 1, dtype=np.int64)
    X = (n % p).astype(np.int64)
    Y = (n // p).astype(np.int64)
    return X, Y


def polynomial_image_coords(N: int, p: int) -> tuple[np.ndarray, np.ndarray]:
    """Phi_Q(n; p) = (n, n^2 mod p). Returns X, Y arrays of size N+1.
    X[0], Y[0] unused."""
    n = np.arange(N + 1, dtype=np.int64)
    X = n
    Y = (n * n) % p
    return X, Y


def run_one(name: str, X_full: np.ndarray, Y_full: np.ndarray,
            primes: np.ndarray, N: int, K: int,
            baseline_trials: int, rng: np.random.Generator,
            top_k: int = 10) -> dict:
    """Compute line cover for primes and matched baselines under given embedding."""
    pi_N = len(primes)
    prime_X = X_full[primes]
    prime_Y = Y_full[primes]
    t0 = time.time()
    L_prime, chosen_prime = bounded_greedy_line_cover(prime_X, prime_Y, K, verbose=False)
    t_prime = time.time() - t0

    L_random_list = []
    for trial in range(baseline_trials):
        sample_idxs = rng.choice(np.arange(1, N + 1), size=pi_N, replace=False)
        rx = X_full[sample_idxs]
        ry = Y_full[sample_idxs]
        L_r, _ = bounded_greedy_line_cover(rx, ry, K, verbose=False)
        L_random_list.append(L_r)
    L_random_mean = sum(L_random_list) / len(L_random_list)
    if len(L_random_list) > 1:
        var = sum((x - L_random_mean) ** 2 for x in L_random_list) / (len(L_random_list) - 1)
        L_random_std = math.sqrt(var)
    else:
        L_random_std = 0.0
    z = (L_random_mean - L_prime) / L_random_std if L_random_std > 0 else float("nan")

    chosen_sorted = sorted(chosen_prime, key=lambda x: -len(x[1]))[:top_k]

    print(f"  [{name}] N={N} pi={pi_N} K={K}: "
          f"L_p={L_prime}  L_r={L_random_mean:.1f}+/-{L_random_std:.1f}  "
          f"L_p/L_r={L_prime/max(L_random_mean,1e-9):.3f}  z={z:+.2f}  "
          f"(t_prime={t_prime:.1f}s)")
    sys.stdout.flush()

    return {
        "name": name,
        "N": N,
        "pi_N": pi_N,
        "K": K,
        "L_primes": L_prime,
        "L_random_mean": L_random_mean,
        "L_random_std": L_random_std,
        "L_random_trials": L_random_list,
        "z": z,
        "n_baseline_trials": baseline_trials,
        "top_lines": [((a, b, c), len(members)) for (a, b, c), members in chosen_sorted],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=100_000)
    parser.add_argument("--K", type=int, default=5,
                        help="bounded-direction K_max (small for narrow embeddings)")
    parser.add_argument("--baseline-trials", type=int, default=3)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--p-list", type=str, default="2,3,5,7,11",
                        help="comma-separated list of small primes p for embeddings")
    parser.add_argument("--top-k-lines", type=int, default=10)
    parser.add_argument("--out-dir", type=str, default=".")
    parser.add_argument("--skip-residue", action="store_true")
    parser.add_argument("--skip-poly", action="store_true")
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    p_list = [int(p) for p in args.p_list.split(",") if p.strip()]
    print(f"# Thread 11 / Slot 2 — alternative 2D embeddings")
    print(f"# N={args.N} K={args.K} p_list={p_list} baseline_trials={args.baseline_trials}")
    sys.stdout.flush()

    t0 = time.time()
    primes = sieve_primes(args.N)
    pi_N = len(primes)
    print(f"# pi(N)={pi_N}, sieve in {time.time()-t0:.1f}s")
    sys.stdout.flush()

    rows = []
    top_line_rows = []

    for p in p_list:
        if not args.skip_residue:
            print(f"\n## Residue-class grid mod p={p}: Phi(n) = (n mod {p}, n//{p})")
            X, Y = residue_class_coords(args.N, p)
            name = f"residue_p{p}"
            r = run_one(name, X, Y, primes, args.N, args.K,
                        args.baseline_trials, rng, args.top_k_lines)
            r["embedding"] = "residue_class"
            r["p"] = p
            rows.append(r)
            for rank, (key, count) in enumerate(r["top_lines"], 1):
                top_line_rows.append((name, args.N, rank, key[0], key[1], key[2], count))

        if not args.skip_poly:
            print(f"\n## Polynomial-image grid mod p={p}: Phi(n) = (n, n^2 mod {p})")
            X, Y = polynomial_image_coords(args.N, p)
            name = f"poly_p{p}"
            r = run_one(name, X, Y, primes, args.N, args.K,
                        args.baseline_trials, rng, args.top_k_lines)
            r["embedding"] = "polynomial_image"
            r["p"] = p
            rows.append(r)
            for rank, (key, count) in enumerate(r["top_lines"], 1):
                top_line_rows.append((name, args.N, rank, key[0], key[1], key[2], count))

    # Save summary CSV.
    summary_path = os.path.join(out_dir, f"summary_alt_N{args.N}_K{args.K}.csv")
    with open(summary_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["name", "embedding", "p", "N", "pi_N", "K",
                    "L_primes", "L_random_mean", "L_random_std",
                    "n_baseline_trials", "L_p_over_L_r", "z_score"])
        for r in rows:
            ratio = r["L_primes"] / r["L_random_mean"] if r["L_random_mean"] > 0 else float("nan")
            w.writerow([r["name"], r["embedding"], r["p"], r["N"], r["pi_N"], r["K"],
                        r["L_primes"], f"{r['L_random_mean']:.2f}",
                        f"{r['L_random_std']:.2f}", r["n_baseline_trials"],
                        f"{ratio:.4f}",
                        f"{r['z']:.3f}" if not math.isnan(r["z"]) else "NA"])

    toplines_path = os.path.join(out_dir, f"top_lines_alt_N{args.N}_K{args.K}.csv")
    with open(toplines_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["name", "N", "rank", "a", "b", "intercept", "prime_count"])
        for row in top_line_rows:
            w.writerow(row)

    print(f"\n# wrote {summary_path}")
    print(f"# wrote {toplines_path}")
    print(f"# DONE alt embeddings N={args.N}")


if __name__ == "__main__":
    main()
