"""Causal-state / excess-entropy complexity of prime parity stream.

Tests whether the binary sequence q(n) = (p(n) mod 6 - 1) / 4  in {0, 1}
(0 if p(n) ≡ 1 mod 6, 1 if p(n) ≡ 5 mod 6) admits a finite-state stochastic
generator, i.e. has finite excess entropy E.

E = sum_{L>=1} (H_L - L * h_mu)  where  h_mu = lim_L (H_L - H_{L-1}).

If E plateaus as L grows -> finite-state generator exists -> potential
polylog model. If E grows linearly -> no such generator exists.

Comparison: Berlekamp-Massey linear complexity over GF(2) is already known
to be N/2 (maximal). Excess entropy is the *stochastic* analogue and could
in principle be small even when linear complexity is maximal (an HMM is
strictly more expressive than an LFSR).

Pass/fail criterion: E_L for L = 2..18 either plateaus (pass) or grows
roughly linearly (fail). Standard random binary stream has E = 0; this is
the strict null.
"""
from __future__ import annotations

import math
import time
from collections import Counter

import numpy as np
from sympy import isprime, primerange


def primes_upto(N):
    return list(primerange(2, N))


def parity_residue_stream(primes):
    """For each p > 3, encode p mod 6 as a bit (0 or 1)."""
    bits = []
    for p in primes:
        if p in (2, 3):
            continue
        r = p % 6
        bits.append(0 if r == 1 else 1)
    return np.array(bits, dtype=np.int8)


def block_entropy(seq: np.ndarray, L: int) -> float:
    """Empirical block entropy H_L = -sum p(b) log2 p(b)."""
    if L == 0:
        return 0.0
    n = len(seq) - L + 1
    if n <= 0:
        return 0.0
    counter = Counter()
    for i in range(n):
        counter[bytes(seq[i : i + L].tobytes())] += 1
    p = np.array(list(counter.values()), dtype=np.float64) / n
    return float(-(p * np.log2(p)).sum())


def baseline_random(N, p1):
    rng = np.random.default_rng(42)
    return (rng.random(N) < p1).astype(np.int8)


def main():
    t0 = time.time()
    print("Generating primes up to N=2_000_000 ...")
    pr = primes_upto(2_000_000)
    seq = parity_residue_stream(pr)
    n = len(seq)
    p1 = float(seq.mean())
    print(f"len(seq) = {n}, P(bit=1) = {p1:.6f}")

    rand = baseline_random(n, p1)

    Lmax = 18
    H_seq = []
    H_rand = []
    for L in range(1, Lmax + 1):
        h_s = block_entropy(seq, L)
        h_r = block_entropy(rand, L)
        H_seq.append(h_s)
        H_rand.append(h_r)
        print(f"L={L:2d}  H_L(prime)={h_s:.6f}  H_L(rand)={h_r:.6f}")

    # entropy rate estimate.
    # Use Delta_L at a moderate L (least biased range L=8..10) rather than
    # the contaminated H_Lmax - H_{Lmax-1}. At large L the empirical block
    # frequencies suffer severe finite-sample bias since 2^L approaches the
    # number of L-grams.
    deltas_seq = [H_seq[L] - H_seq[L - 1] for L in range(1, Lmax)]
    deltas_rand = [H_rand[L] - H_rand[L - 1] for L in range(1, Lmax)]
    # average over L=8..10 (indices 7..9 in 0-indexed)
    h_mu_seq_estimate = float(np.mean(deltas_seq[7:10]))
    h_mu_rand_estimate = float(np.mean(deltas_rand[7:10]))
    h_mu_seq = H_seq[-1] - H_seq[-2]  # original bad estimate, kept for record
    h_mu_rand = H_rand[-1] - H_rand[-2]
    print(f"\nh_mu (prime) [L=18-17 raw, biased low] ~ {h_mu_seq:.6f}")
    print(f"h_mu (rand)  [L=18-17 raw, biased low] ~ {h_mu_rand:.6f}")
    print(f"h_mu (prime) [Delta L=8..10 mean] ~ {h_mu_seq_estimate:.6f}")
    print(f"h_mu (rand)  [Delta L=8..10 mean] ~ {h_mu_rand_estimate:.6f}")
    print(f"single-bit entropy h(p1) = {-(p1*math.log2(p1)+(1-p1)*math.log2(1-p1)):.6f}")
    # Use the moderate-L estimate for excess-entropy computation below.
    h_mu_seq = h_mu_seq_estimate
    h_mu_rand = h_mu_rand_estimate

    # Excess entropy: E_L = H_L - L*h_mu.
    # For a stationary process, H_L = E + L*h_mu + O(decaying), so E_L should
    # plateau. Plateau value > 0 means finite stochastic memory; growing
    # E_L means infinite memory (no finite-state generator).
    print("\nExcess entropy estimate  E_L = H_L - L*h_mu")
    print("(if E_L plateaus -> finite-state generator; if grows -> no)")
    E_seq, E_rand = [], []
    for L in range(1, Lmax + 1):
        e_s = H_seq[L - 1] - L * h_mu_seq
        e_r = H_rand[L - 1] - L * h_mu_rand
        E_seq.append(e_s)
        E_rand.append(e_r)
        print(f"L={L:2d}  E_L(prime)={e_s:+.6f}  E_L(rand)={e_r:+.6f}")

    # Also compute the *finite-difference* per-block entropy increase:
    # Delta_L = H_L - H_{L-1}, should -> h_mu.
    print("\nDelta_L = H_L - H_{L-1}   (should -> h_mu)")
    for L in range(2, Lmax + 1):
        d_s = H_seq[L - 1] - H_seq[L - 2]
        d_r = H_rand[L - 1] - H_rand[L - 2]
        print(f"L={L:2d}  Delta_L(prime)={d_s:.6f}  Delta_L(rand)={d_r:.6f}")

    # Verdict: compare E_seq tail to E_rand tail.
    tail_seq = np.mean(E_seq[-4:])
    tail_rand = np.mean(E_rand[-4:])
    print("\nMean E_L over last 4 L values:")
    print(f"  prime parity: {tail_seq:+.6f}")
    print(f"  random null : {tail_rand:+.6f}")
    print(f"  difference  : {tail_seq - tail_rand:+.6f}")

    # If prime parity has h_mu ~ 1 (max entropy) and E_L ~ 0 like random, then
    # it's *indistinguishable* from random in this measure -> no finite-state
    # advantage.
    # If h_mu < 1 with E_L plateauing positive, there is finite-state structure.
    elapsed = time.time() - t0
    print(f"\n[elapsed: {elapsed:.1f}s]")

    return {
        "n": n,
        "p1": p1,
        "h_mu_seq": h_mu_seq,
        "h_mu_rand": h_mu_rand,
        "E_seq": E_seq,
        "E_rand": E_rand,
        "Delta_seq": [H_seq[L] - H_seq[L - 1] for L in range(1, Lmax)],
        "Delta_rand": [H_rand[L] - H_rand[L - 1] for L in range(1, Lmax)],
    }


if __name__ == "__main__":
    main()
