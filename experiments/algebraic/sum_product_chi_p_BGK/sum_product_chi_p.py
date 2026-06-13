"""
Bourgain-Glibichuk-Konyagin sum-product gain on the prime set in F_p.

For a prime p and the set A_prime = {primes < p} (or sparser cuts at K < p),
compute |A+A| and |A.A| in F_p and compare to random subsets of matched
cardinality.

Cross-domain reference: BGK 2006 J. London Math. Soc.; Garaev 2008
(arXiv:math/0702780); Iosevich 2009 (arXiv:0904.2075); Tao-Vu 2006
"Additive Combinatorics".

Falsifiability:
  F_A: g(prime) > g(random) by a factor that GROWS with p AND has a
       closed-form HL singular-series interpretation -> A-grade.
  F_B: g(prime) - g(random) = o(1) within sample noise -> B-grade
       (38th pseudorandomness measure, joint additive-multiplicative
       BGK-saturation of primes).
  F_I: g(prime) < g(random) consistently -> B-grade
       (primes are sub-BGK; HL constraints reduce sum-product gain).

Edges cited: E2.13 (Gowers; HL singular series), E1.10 (pseudorandomness
battery), E2.14 (Anderson localisation HL detection).
"""

import numpy as np
from sympy import primerange, isprime, primitive_root, n_order
import json
import time
import sys
import argparse

def sieve_primes_up_to(N):
    """Sieve of Eratosthenes, returns sorted list of primes <= N."""
    if N < 2:
        return []
    sieve = np.ones(N + 1, dtype=bool)
    sieve[0] = sieve[1] = False
    for i in range(2, int(N**0.5) + 1):
        if sieve[i]:
            sieve[i*i::i] = False
    return np.where(sieve)[0].tolist()

def sumset_size_in_Zp(A, p):
    """Compute |A + A mod p| via FFT convolution.

    A: 1D array of residues in [0, p).
    Returns the cardinality of {a + b mod p : a, b in A}.
    """
    I = np.zeros(p, dtype=np.float64)
    I[A] = 1.0
    F = np.fft.rfft(I)
    conv = np.fft.irfft(F * F, n=p)
    # rounding tolerance
    return int(np.sum(conv > 0.5))

def productset_size_in_Fpstar(A, p, dlog_table):
    """Compute |A . A mod p| in F_p^*.

    Maps A -> log_g(A) in Z/(p-1)Z then computes |log_A + log_A| there.

    A: residues in [1, p) (assumed all nonzero).
    dlog_table: dlog_table[r] = k for r = g^k mod p, k in [0, p-2].
    Returns |A.A|.
    """
    A = np.asarray(A, dtype=np.int64)
    A = A[A != 0]
    if len(A) == 0:
        return 0
    log_A = dlog_table[A]
    n = p - 1
    I = np.zeros(n, dtype=np.float64)
    I[log_A] = 1.0
    F = np.fft.rfft(I)
    conv = np.fft.irfft(F * F, n=n)
    return int(np.sum(conv > 0.5))

def build_dlog_table(p, g):
    """Return array dlog of size p where dlog[g^k mod p] = k."""
    dlog = np.full(p, -1, dtype=np.int64)
    val = 1
    for k in range(p - 1):
        dlog[val] = k
        val = (val * g) % p
    return dlog

def gain(card_AsumA, card_AmulA, A_size):
    """BGK gain: max(|A+A|, |A.A|) / |A|^{15/14}."""
    if A_size <= 0:
        return 0.0
    return max(card_AsumA, card_AmulA) / (A_size ** (15.0/14.0))

def integer_sumset_size(A):
    """|A + A| as integer multiset turned into set (numpy-vectorised)."""
    A = np.asarray(A, dtype=np.int64)
    if len(A) == 0:
        return 0
    sums = (A[:, None] + A[None, :]).ravel()
    return int(np.unique(sums).size)

def integer_productset_size(A):
    """|A * A| as integer set (numpy-vectorised)."""
    A = np.asarray(A, dtype=np.int64)
    if len(A) == 0:
        return 0
    prods = (A[:, None] * A[None, :]).ravel()
    return int(np.unique(prods).size)

def _stats(arr):
    arr = np.asarray(arr, dtype=np.float64)
    return float(arr.mean()), float(arr.std(ddof=0))

def run_one_p(p, n_random=100, K_factor=1.0, seed=0, verbose=True):
    """Run the BGK probe for one prime p.

    Three control distributions to disentangle support / parity / multiplicative
    structure:
      B1 = uniform random in F_p* (Garaev / textbook BGK control;
            tests primes vs maximally-spread set in F_p).
      B2 = uniform random subset of integers in [2, K] (matched support;
            tests primes vs random integers from the same range).
      B3 = uniform random subset of {odd integers in [3, K] union {2}}
            (matched support + matched parity).

    K_factor: take primes <= p**(7/13) * K_factor as A_prime.
              K_factor = 1.0 puts us at the BGK boundary; larger K_factor
              broadens the support (and saturates the sumset).
    Returns dict with all measurements.
    """
    rng = np.random.default_rng(seed)
    g = int(primitive_root(p))
    if verbose:
        print(f"p={p}, primitive root g={g}", flush=True)

    t0 = time.time()
    dlog = build_dlog_table(p, g)
    t_dlog = time.time() - t0
    if verbose:
        print(f"  dlog table built in {t_dlog:.2f}s", flush=True)

    K_bgk = p ** (7.0 / 13.0)
    K = int(K_bgk * K_factor)
    K = min(K, p - 1)
    primes_below_K = sieve_primes_up_to(K)
    A_prime = np.array([q for q in primes_below_K if q != p and q < p], dtype=np.int64)
    A_size = len(A_prime)
    if verbose:
        print(f"  K = {K} (BGK boundary {K_bgk:.1f} * {K_factor}); "
              f"|A_prime| = {A_size}", flush=True)

    if A_size < 5:
        return {"p": p, "K": K, "A_size": A_size, "skipped": True}

    t0 = time.time()
    sumset_prime = sumset_size_in_Zp(A_prime, p)
    productset_prime = productset_size_in_Fpstar(A_prime, p, dlog)
    sumset_prime_int = integer_sumset_size(A_prime)
    productset_prime_int = integer_productset_size(A_prime)
    expected_prime_prod_int = A_size * (A_size + 1) // 2  # unique factorisation
    g_prime = gain(sumset_prime, productset_prime, A_size)
    t_prime = time.time() - t0
    if verbose:
        print(f"  PRIME: |A+A|_p={sumset_prime}, |A+A|_Z={sumset_prime_int}; "
              f"|A.A|_p={productset_prime}, |A.A|_Z={productset_prime_int} "
              f"(expected N(N+1)/2={expected_prime_prod_int}); "
              f"gain={g_prime:.4f} ({t_prime:.2f}s)", flush=True)

    # Pre-compute candidate pools for each control type
    pool_B1 = np.arange(1, p)  # F_p^*
    pool_B2 = np.arange(2, K + 1)  # integers in [2, K]
    odd_in_range = np.arange(3, K + 1, 2, dtype=np.int64)
    if A_size >= 1:
        pool_B3 = np.concatenate([np.array([2], dtype=np.int64), odd_in_range])
    else:
        pool_B3 = odd_in_range
    # B4: W-tricked control — match coprimality to 2, 3 (i.e., n mod 6 in {1, 5})
    # plus include the small primes {2, 3} themselves for parity match.
    coprime_6 = np.array([n for n in range(5, K + 1) if n % 2 == 1 and n % 3 != 0],
                         dtype=np.int64)
    extras = np.array([n for n in [2, 3] if n <= K], dtype=np.int64)
    pool_B4 = np.concatenate([extras, coprime_6])

    controls = {"B1_Fpstar": pool_B1, "B2_match_support": pool_B2,
                "B3_match_parity": pool_B3, "B4_W6_tricked": pool_B4}
    control_results = {}

    for name, pool in controls.items():
        if A_size > len(pool):
            control_results[name] = {"skipped": True,
                                     "pool_size": int(len(pool)),
                                     "A_size": A_size}
            continue
        ssums, sprods, gns = [], [], []
        ssums_int, sprods_int = [], []
        t0 = time.time()
        for trial in range(n_random):
            B = rng.choice(pool, size=A_size, replace=False)
            B_int = B.astype(np.int64)
            B = B_int % p
            B = B[B != 0]
            if len(B) == 0:
                continue
            ssum = sumset_size_in_Zp(B, p)
            sprod = productset_size_in_Fpstar(B, p, dlog)
            ssums.append(ssum)
            sprods.append(sprod)
            gns.append(gain(ssum, sprod, len(B)))
            # Integer sets only for matched-support controls (where pool subset of [2,K])
            if name in ("B2_match_support", "B3_match_parity", "B4_W6_tricked"):
                ssums_int.append(integer_sumset_size(B_int))
                sprods_int.append(integer_productset_size(B_int))
        t = time.time() - t0
        ssum_mean, ssum_std = _stats(ssums)
        sprod_mean, sprod_std = _stats(sprods)
        gn_mean, gn_std = _stats(gns)
        z_s = (sumset_prime - ssum_mean) / max(ssum_std, 1e-12)
        z_p = (productset_prime - sprod_mean) / max(sprod_std, 1e-12)
        z_g = (g_prime - gn_mean) / max(gn_std, 1e-12)
        rec = {
            "sumset_mean": ssum_mean, "sumset_std": ssum_std,
            "productset_mean": sprod_mean, "productset_std": sprod_std,
            "gain_mean": gn_mean, "gain_std": gn_std,
            "z_sumset": z_s, "z_productset": z_p, "z_gain": z_g,
            "n_random": n_random, "time_seconds": t,
            "pool_size": int(len(pool)),
        }
        if ssums_int:
            ssum_int_mean, ssum_int_std = _stats(ssums_int)
            sprod_int_mean, sprod_int_std = _stats(sprods_int)
            z_s_int = (sumset_prime_int - ssum_int_mean) / max(ssum_int_std, 1e-12)
            z_p_int = (productset_prime_int - sprod_int_mean) / max(sprod_int_std, 1e-12)
            rec.update({
                "sumset_int_mean": ssum_int_mean, "sumset_int_std": ssum_int_std,
                "productset_int_mean": sprod_int_mean, "productset_int_std": sprod_int_std,
                "z_sumset_int": z_s_int, "z_productset_int": z_p_int,
            })
        control_results[name] = rec
        if verbose:
            line = (f"  {name:>20s}: "
                    f"|A+A|_p={ssum_mean:8.1f}+-{ssum_std:5.2f} (z={z_s:7.2f}); "
                    f"|A.A|_p={sprod_mean:8.1f}+-{sprod_std:5.2f} (z={z_p:7.2f})")
            if ssums_int:
                line += (f"; |A.A|_Z={sprod_int_mean:8.1f}+-{sprod_int_std:5.2f} "
                         f"(z_int={z_p_int:6.2f})")
            print(line + f" ({t:.1f}s)", flush=True)

    BGK_lower_bound = (A_size ** (15.0/14.0)) / max(np.log(max(A_size, 2)), 1.0) ** (2.0/7.0)

    return {
        "p": p, "K": K, "K_factor": K_factor, "A_size": A_size,
        "g_primitive_root": g,
        "sumset_prime": sumset_prime,
        "productset_prime": productset_prime,
        "sumset_prime_int": int(sumset_prime_int),
        "productset_prime_int": int(productset_prime_int),
        "expected_prime_prod_int": int(expected_prime_prod_int),
        "gain_prime": g_prime,
        "controls": control_results,
        "BGK_lower_bound": float(BGK_lower_bound),
        "t_dlog_seconds": t_dlog,
        "t_prime_seconds": t_prime,
        "n_random": n_random,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--primes", type=int, nargs="+",
                        default=[1009, 10007, 100003, 1000003])
    parser.add_argument("--n_random", type=int, default=100)
    parser.add_argument("--K_factors", type=float, nargs="+",
                        default=[0.5, 1.0, 2.0])
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", type=str, default="results.json")
    args = parser.parse_args()

    results = []
    for p in args.primes:
        for kf in args.K_factors:
            print(f"\n=== p={p}, K_factor={kf} ===", flush=True)
            try:
                r = run_one_p(p, n_random=args.n_random, K_factor=kf, seed=args.seed)
                results.append(r)
            except Exception as e:
                print(f"  FAILED: {e}", flush=True)
                results.append({"p": p, "K_factor": kf, "error": str(e)})

    out_path = args.out
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\nWrote {out_path}")

    print("\n=== SUMMARY ===")
    print(f"{'p':>9} {'Kf':>5} {'|A|':>5} {'|A+A|p':>7} {'|A.A|p':>7} {'g_p':>7}  "
          f"{'control':>17} "
          f"{'z_s':>7} {'z_p':>7} {'z_g':>7}")
    for r in results:
        if r.get("skipped") or "error" in r:
            continue
        for name, c in r.get("controls", {}).items():
            if c.get("skipped"):
                continue
            print(f"{r['p']:>9} {r['K_factor']:>5.2f} {r['A_size']:>5} "
                  f"{r['sumset_prime']:>7} {r['productset_prime']:>7} "
                  f"{r['gain_prime']:>7.3f}  "
                  f"{name:>17} "
                  f"{c['z_sumset']:>7.2f} {c['z_productset']:>7.2f} "
                  f"{c['z_gain']:>7.2f}")


if __name__ == "__main__":
    main()
