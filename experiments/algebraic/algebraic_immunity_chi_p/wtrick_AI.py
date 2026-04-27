#!/usr/bin/env python3
"""
W-trick correction to AI(chi_P).

The deg-2 annihilator g(x) = (1 + x_0)(1 + x_1) of chi_P over F_2^N takes
value 1 iff n ≡ 0 mod 4. No prime > 2 has n ≡ 0 mod 4, and chi_P(2) gets
annihilated by x_1 = 1.  So AI(chi_P) = 2 is EXACTLY the mod-4 sieve fact.

If we W-trick by restricting to a residue class b mod W (gcd(b, W) = 1),
the mod-2/mod-4 annihilator vanishes because the W-tricked function
  chi_P_{W,b}(n) := chi_P(W*n + b)
is no longer supported only on "n ≢ 0 mod 4".

Question: does AI(chi_P_{W,b}) match AI(random matched-density Bernoulli)?
If yes -> AI deviation is fully explained by the W-trick (consistent with
E2.13, Gowers U^k singular series picture).
If no  -> a new polynomial-method anomaly persists beyond mod-W structure.

We test W ∈ {2, 6, 30} (cumulative small primes), b = 1 (canonical coprime
class), N up to 12 for tractability.
"""
import time
import json
from pathlib import Path
import numpy as np
from sympy import isprime, mobius
from sympy import factorint

from algebraic_immunity_chi_p import (
    algebraic_immunity_F2,
)


def chi_P_W(n, W, b):
    """chi_P(W*n + b)."""
    return 1 if isprime(W * n + b) else 0


def liouville_pos_W(n, W, b):
    m = W * n + b
    if m <= 0:
        return 0
    if m == 1:
        return 1
    f = factorint(m)
    omega = sum(f.values())
    return 1 if omega % 2 == 0 else 0


def mobius_nonzero_W(n, W, b):
    m = W * n + b
    if m <= 0:
        return 0
    return 1 if mobius(m) != 0 else 0


def truth_table_W(func, W, b, N):
    return np.fromiter((func(n, W, b) for n in range(2 ** N)),
                        dtype=np.int8, count=2 ** N)


def main():
    out_dir = Path(__file__).parent
    results = []
    print(f"{'W':>4} {'b':>3} {'N':>3} | {'rho_chi':>8} | "
          f"{'AI_chi':>6} | {'AI_lam+':>7} | {'AI_mu!=0':>8} | "
          f"{'AI_random_mean':>14} | {'AI_random_std':>13}")
    print("-" * 100)
    n_seeds = 8
    for W, b in [(1, 0), (2, 1), (6, 1), (6, 5), (30, 1), (30, 7), (30, 11)]:
        for N in range(6, 12):
            t0 = time.time()
            tt_chi = truth_table_W(chi_P_W, W, b, N)
            tt_lam = truth_table_W(liouville_pos_W, W, b, N)
            tt_mu = truth_table_W(mobius_nonzero_W, W, b, N)
            rho = float(np.mean(tt_chi))
            rho_lam = float(np.mean(tt_lam))
            rho_mu = float(np.mean(tt_mu))

            ai_chi, _ = algebraic_immunity_F2(tt_chi, N)
            ai_lam, _ = algebraic_immunity_F2(tt_lam, N)
            ai_mu,  _ = algebraic_immunity_F2(tt_mu, N)

            # Random matched-density: at chi_P density rho.
            ais_rand = []
            for s in range(n_seeds):
                rng = np.random.default_rng(s + N * 31 + W * 1009 + b * 13)
                tt_rand = (rng.random(2 ** N) < rho).astype(np.int8)
                ai_r, _ = algebraic_immunity_F2(tt_rand, N)
                ais_rand.append(ai_r)
            t = time.time() - t0
            ais_rand = np.array(ais_rand)
            print(f"{W:>4} {b:>3} {N:>3} | {rho:>8.4f} | "
                  f"{ai_chi:>6} | {ai_lam:>7} | {ai_mu:>8} | "
                  f"{ais_rand.mean():>14.2f} | {ais_rand.std():>13.2f}  ({t:.1f}s)")
            results.append({
                "W": W, "b": b, "N": N,
                "density_chi": rho,
                "density_lam": rho_lam,
                "density_mu": rho_mu,
                "AI_chi_P": int(ai_chi),
                "AI_liouville_pos": int(ai_lam),
                "AI_mobius_nonzero": int(ai_mu),
                "AI_random_mean": float(ais_rand.mean()),
                "AI_random_std": float(ais_rand.std()),
                "AI_random_min": int(ais_rand.min()),
                "AI_random_max": int(ais_rand.max()),
                "n_seeds": n_seeds,
                "time_s": t,
            })
    out = {"AI_W_tricked": results}
    with open(out_dir / "wtrick_AI_data.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {out_dir / 'wtrick_AI_data.json'}")


if __name__ == "__main__":
    main()
