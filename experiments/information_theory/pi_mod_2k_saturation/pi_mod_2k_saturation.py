"""
F2 — Information rate saturation of pi(x) mod 2^k.

Question (NOVELTY_CHALLENGES.md §2.F2):
  Is H(pi(x) mod 2^k | pi(x-1) mod 2^k) = 0.537 * k bits per step
  (linear in k), or does it saturate?

Edges touched:
  E1.5 — pi(x) mod m has invariant conditional entropy 0.537 bits
         (tested only for m in {2..30}; this script extends to 2^k, k=1..10
         and many X scales).

Sharpened prediction (this experiment):
  For every modulus m >= 2 and every X >= 100,
     H( pi(x) mod m | pi(x-1) mod m ; x in [1, X] )
       =  h_2( pi(X) / X )  +  O(1/X)
  where h_2(p) = -p log_2 p - (1-p) log_2 (1-p).

Mechanism (why this should hold):
  pi(x) - pi(x-1) = 1[x is prime] in {0,1}.
  Therefore  pi(x) mod m  =  ( pi(x-1) mod m + 1[x prime] ) mod m.
  Conditioning on pi(x-1) mod m, only the prime indicator is random.
  If 1[x prime] is asymptotically conditionally independent of pi(x-1) mod m
  (very mild assumption -- known empirically to high precision; see
  pi_modular_structure_results.md), then
     H(Y|X) = h_2( P[x prime] )  =  h_2( pi(X)/X ).
  Crucially this is INDEPENDENT of m.

What gets tested:
  (i)  H_emp(m, X) is empirically constant in m (across 2^k for k=1..10
       and across primes / primorial moduli for cross-check).
  (ii) H_emp(2, X) tracks h_2(pi(X)/X) at multiple scales X.
  (iii) The "0.537 bits" of E1.5 is X-specific (X ~ 10^4).

Falsification:
  Fail if max over m of |H_emp(m, X) - h_2(pi(X)/X)| > 0.005 at any X >= 10^4.
  Or: if H_emp(2^{k+1}, X) - H_emp(2^k, X) > 0.01 at any (k, X) -- this would
  refute saturation and indicate some bit-level information addition.

Save under: experiments/information_theory/pi_mod_2k_saturation/
"""

import math
import sys
from collections import Counter

import numpy as np


def sieve_pi_sequence(X):
    """Return numpy array pi[0..X] where pi[x] = #primes <= x."""
    if X < 2:
        return np.zeros(X + 1, dtype=np.int64)
    is_prime = np.ones(X + 1, dtype=bool)
    is_prime[:2] = False
    for p in range(2, int(math.isqrt(X)) + 1):
        if is_prime[p]:
            is_prime[p * p :: p] = False
    pi = np.cumsum(is_prime, dtype=np.int64)
    return pi


def conditional_entropy(pi_seq, m):
    """
    Compute H( Y | X ) where X = pi[x-1] mod m, Y = pi[x] mod m,
    averaged over x in [1, len(pi_seq)-1].
    Returns (H_cond, H_marginal, h2_predicted, prime_density).
    """
    seq = pi_seq % m

    # Joint counts for adjacent pairs (X, Y) = (seq[i-1], seq[i]).
    # The increment Y - X mod m is in {0, 1} since pi(x) - pi(x-1) in {0,1}.
    # We can compute H(Y|X) directly using only the increment statistics:
    #   For each x_state s (=seq[i-1]), count how often Y=s vs Y=(s+1) mod m.
    transitions_zero = np.zeros(m, dtype=np.int64)  # Y == X (no prime at x)
    transitions_one = np.zeros(m, dtype=np.int64)  # Y == X+1 (prime at x)

    diffs = (seq[1:] - seq[:-1]) % m  # 0 or 1
    x_states = seq[:-1]

    # For each m, increment is exactly 0 or 1 (pi grows by 0/1 per step).
    mask_zero = diffs == 0
    mask_one = diffs == 1

    # Count (x_state, increment) pairs.
    np.add.at(transitions_zero, x_states[mask_zero], 1)
    np.add.at(transitions_one, x_states[mask_one], 1)

    total_pairs = len(diffs)
    H_cond = 0.0
    for s in range(m):
        n_total = transitions_zero[s] + transitions_one[s]
        if n_total == 0:
            continue
        p_x = n_total / total_pairs
        p0 = transitions_zero[s] / n_total
        p1 = transitions_one[s] / n_total
        h_y_given_xs = 0.0
        if p0 > 0:
            h_y_given_xs -= p0 * math.log2(p0)
        if p1 > 0:
            h_y_given_xs -= p1 * math.log2(p1)
        H_cond += p_x * h_y_given_xs

    # Marginal entropy of Y in the sequence.
    counts_y = np.bincount(seq, minlength=m).astype(np.float64)
    counts_y = counts_y[counts_y > 0]
    p_y = counts_y / counts_y.sum()
    H_marg = float(-(p_y * np.log2(p_y)).sum())

    # Predicted by closed form: h_2(pi(X)/X) where X = len(pi_seq)-1.
    X = len(pi_seq) - 1
    n_primes = int(pi_seq[-1])
    p_dens = n_primes / X
    if 0 < p_dens < 1:
        h2_pred = -p_dens * math.log2(p_dens) - (1 - p_dens) * math.log2(1 - p_dens)
    else:
        h2_pred = 0.0

    return H_cond, H_marg, h2_pred, p_dens


def run_experiment(X_values, moduli):
    """Returns dict: results[(X, m)] = (H_cond, H_marg, h2_pred, p_dens)."""
    results = {}
    for X in X_values:
        print(f"\n[sieve] computing pi(x) for x in [0, {X}] ...", flush=True)
        pi_seq = sieve_pi_sequence(X)
        for m in moduli:
            H_c, H_m, h2_pred, p_dens = conditional_entropy(pi_seq, m)
            results[(X, m)] = (H_c, H_m, h2_pred, p_dens)
            print(
                f"  X={X:>10}  m={m:>5}  H(Y|X)={H_c:.6f}  H(Y)={H_m:.6f}  "
                f"h2(pi/X)={h2_pred:.6f}  diff={H_c - h2_pred:+.6f}",
                flush=True,
            )
    return results


def report(results, X_values, moduli):
    """Build human-readable summary tables."""
    out = []
    out.append("# Per-X table: H(Y|X) for each modulus m, with closed-form prediction\n")

    out.append("| X | h_2(pi(X)/X) (predicted) |")
    out[-1] += " ".join(f" m={m} |" for m in moduli)
    out.append("|" + "---|" * (2 + len(moduli)))
    for X in X_values:
        row = [f"{X}"]
        h2_pred = results[(X, moduli[0])][2]
        row.append(f"{h2_pred:.6f}")
        for m in moduli:
            H_c = results[(X, m)][0]
            row.append(f"{H_c:.6f}")
        out.append("| " + " | ".join(row) + " |")

    out.append("\n# Per-X table: H(Y|X) - h_2(pi(X)/X) (residual)\n")
    out.append("| X |" + " ".join(f" m={m} |" for m in moduli))
    out.append("|" + "---|" * (1 + len(moduli)))
    for X in X_values:
        row = [f"{X}"]
        h2_pred = results[(X, moduli[0])][2]
        for m in moduli:
            H_c = results[(X, m)][0]
            row.append(f"{H_c - h2_pred:+.2e}")
        out.append("| " + " | ".join(row) + " |")

    out.append("\n# Saturation check: max_m H(Y|X) - min_m H(Y|X) at fixed X\n")
    out.append("| X | max_m | min_m | spread |")
    out.append("|---|---|---|---|")
    for X in X_values:
        Hs = [results[(X, m)][0] for m in moduli]
        out.append(f"| {X} | {max(Hs):.6f} | {min(Hs):.6f} | {max(Hs) - min(Hs):.2e} |")

    out.append("\n# Density check: pi(X)/X and h_2(pi(X)/X)\n")
    out.append("| X | pi(X) | pi(X)/X | h_2(pi(X)/X) |")
    out.append("|---|---|---|---|")
    for X in X_values:
        _, _, h2_pred, p_dens = results[(X, moduli[0])]
        n_primes = int(round(p_dens * X))
        out.append(f"| {X} | {n_primes} | {p_dens:.6f} | {h2_pred:.6f} |")

    return "\n".join(out)


def main():
    # Powers of 2 (the F2 family) plus cross-check moduli.
    pow2_moduli = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
    cross_check_moduli = [3, 5, 7, 11, 13, 30, 210]
    moduli = sorted(set(pow2_moduli + cross_check_moduli))

    # Scales spanning the regime where the "0.537 bits" was originally measured
    # (X ~ 10^4) up to large X where prime density is much smaller.
    if "--quick" in sys.argv:
        X_values = [10**3, 10**4, 10**5]
    else:
        X_values = [10**3, 10**4, 10**5, 10**6, 10**7]

    print(f"# F2 experiment: H(pi(x) mod 2^k | pi(x-1) mod 2^k)")
    print(f"# Moduli (powers of 2 + cross-check): {moduli}")
    print(f"# X values: {X_values}")

    results = run_experiment(X_values, moduli)

    # Falsification logic.
    print("\n# Falsification check\n")
    fail_predict = []
    fail_satur = []
    for X in X_values:
        h2_pred = results[(X, moduli[0])][2]
        for m in moduli:
            H_c = results[(X, m)][0]
            if abs(H_c - h2_pred) > 0.005 and X >= 10**4:
                fail_predict.append((X, m, H_c, h2_pred))
        Hs = [results[(X, m)][0] for m in moduli]
        spread = max(Hs) - min(Hs)
        if spread > 0.01:
            fail_satur.append((X, spread))

    if fail_predict:
        print("FAIL (closed-form prediction):")
        for row in fail_predict:
            print("  ", row)
    else:
        print("PASS (closed-form): max |H_emp - h_2(pi/X)| <= 0.005 at all X >= 10^4")

    if fail_satur:
        print("FAIL (saturation):")
        for row in fail_satur:
            print("  ", row)
    else:
        print("PASS (saturation): max-m - min-m H_emp spread <= 0.01 at all X")

    table = report(results, X_values, moduli)
    print("\n" + table)


if __name__ == "__main__":
    main()
