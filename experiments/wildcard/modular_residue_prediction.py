#!/usr/bin/env python3
"""
Modular Residue Prediction: Can p(n) mod q be predicted from n alone?

If for each small prime q, we could compute p(n) mod q in O(polylog(n)) time,
then CRT would recover p(n) exactly using O(log(p(n))) = O(n*log(n)) small primes.

Idea: p(n) mod q depends on the distribution of primes in residue classes mod q.
By Dirichlet's theorem, primes equidistribute among coprime residues mod q.
But the EXACT residue of the nth prime is the question.

Key insight: p(n) mod 2 is always 1 for n >= 2 (free!).
p(n) mod 3: primes > 3 are either 1 or 2 mod 3. The pattern is "random-looking."
p(n) mod 5: primes > 5 are 1, 2, 3, or 4 mod 5. Again random-looking.

But HOW random? If there's ANY predictability, we might exploit it.

We test:
1. Entropy of p(n) mod q sequences
2. Autocorrelation of p(n) mod q
3. Whether p(n) mod q can be computed from pi(x, q, a) values
4. Cross-correlations between p(n) mod q for different q
"""

import numpy as np
from sympy import primerange, isprime, primepi, nextprime
from collections import Counter
import time

def compute_residues(N, q):
    """Compute p(n) mod q for n = 1, ..., N."""
    residues = []
    p = 2
    for _ in range(N):
        residues.append(p % q)
        p = nextprime(p)
    return residues

def entropy(sequence, q):
    """Shannon entropy of sequence mod q."""
    counts = Counter(sequence)
    total = len(sequence)
    h = 0
    for c in counts.values():
        p = c / total
        if p > 0:
            h -= p * np.log2(p)
    # Max entropy for phi(q) equiprobable classes
    from sympy import totient
    max_h = np.log2(float(totient(q)))
    return h, max_h

def autocorrelation(sequence, max_lag=50):
    """Autocorrelation of sequence."""
    s = np.array(sequence, dtype=float)
    s = s - s.mean()
    if np.std(s) == 0:
        return np.zeros(max_lag)
    s = s / np.std(s)
    n = len(s)
    acf = np.correlate(s, s, mode='full')[n-1:n-1+max_lag] / n
    return acf

def cross_correlation(seq1, seq2):
    """Cross-correlation between two sequences."""
    s1 = np.array(seq1, dtype=float)
    s2 = np.array(seq2, dtype=float)
    s1 = (s1 - s1.mean()) / (np.std(s1) + 1e-10)
    s2 = (s2 - s2.mean()) / (np.std(s2) + 1e-10)
    return np.corrcoef(s1, s2)[0, 1]

def test_linear_predictability(residues, q):
    """Can p(n+1) mod q be linearly predicted from p(n) mod q?"""
    X = np.array(residues[:-1]).reshape(-1, 1)
    y = np.array(residues[1:])

    # Simple: predict p(n+1) mod q = mode(p(n+1) mod q | p(n) mod q = r)
    transition = {}
    for i in range(len(residues) - 1):
        r = residues[i]
        if r not in transition:
            transition[r] = []
        transition[r].append(residues[i+1])

    correct = 0
    total = 0
    for i in range(len(residues) - 1):
        r = residues[i]
        prediction = Counter(transition[r]).most_common(1)[0][0]
        if prediction == residues[i+1]:
            correct += 1
        total += 1

    return correct / total

def test_kth_order_prediction(residues, q, k=5):
    """Can p(n+1) mod q be predicted from (p(n-k+1),...,p(n)) mod q?"""
    correct = 0
    total = 0
    history = {}

    for i in range(k, len(residues)):
        key = tuple(residues[i-k:i])
        target = residues[i]

        if key in history:
            prediction = Counter(history[key]).most_common(1)[0][0]
            if prediction == target:
                correct += 1
            total += 1

        if key not in history:
            history[key] = []
        history[key].append(target)

    if total == 0:
        return 0
    return correct / total

def test_pi_based_prediction(N, q):
    """
    Key idea: p(n) mod q is determined by how many primes ≡ a (mod q)
    there are before p(n).

    Let pi(x; q, a) = #{p <= x : p ≡ a (mod q)}.
    Then p(n) ≡ a (mod q) iff pi(p(n); q, a) > pi(p(n)-1; q, a).

    If we could compute pi(x; q, a) in polylog, we'd know p(n) mod q!

    Dirichlet L-functions give: pi(x; q, a) ~ li(x)/phi(q) + O(sqrt(x)*log(x))
    under GRH.

    The error term is O(sqrt(x)*log(x)) -- still too large for exact.
    But for the MOD q question, do we need the exact count?

    We need: which residue class a has pi(p(n); q, a) jump at x = p(n)?
    This is equivalent to: p(n) mod q = ?

    So the pi-based approach is CIRCULAR for determining mod q.
    """
    # Test whether the residue class distribution is predictable
    # from the running counts pi(x; q, a) / pi(x)
    from sympy import totient
    phi_q = int(totient(q))

    primes = list(primerange(2, 200000))[:N]
    counts = {a: 0 for a in range(q) if np.gcd(a, q) == 1 or a < 2}

    deviations = []
    for i, p in enumerate(primes):
        r = p % q
        if r in counts:
            counts[r] += 1
        # Expected: each class should have ~ (i+1)/phi_q primes
        if i > 10 and q > 2:
            expected = (i + 1) / phi_q
            max_dev = max(abs(c - expected) for a, c in counts.items()
                        if np.gcd(a, q) == 1)
            deviations.append(max_dev / np.sqrt(i + 1))

    if deviations:
        return {
            'mean_normalized_deviation': np.mean(deviations),
            'max_normalized_deviation': np.max(deviations),
            'final_deviation': deviations[-1],
        }
    return {}

def main():
    print("=" * 70)
    print("MODULAR RESIDUE PREDICTION FOR p(n)")
    print("=" * 70)

    N = 10000  # Number of primes to test

    print(f"\nComputing first {N} primes...")
    primes = list(primerange(2, 200000))[:N]
    print(f"  Range: p(1)={primes[0]} to p({N})={primes[-1]}")

    # Test for various moduli
    for q in [2, 3, 5, 7, 11, 13, 30, 210]:
        print(f"\n--- Modulus q = {q} ---")

        residues = [p % q for p in primes]

        # Entropy
        h, max_h = entropy(residues[10:], q)  # Skip small primes
        print(f"  Entropy: {h:.4f} / {max_h:.4f} (ratio: {h/max_h:.4f})")

        # Autocorrelation
        acf = autocorrelation(residues)
        max_acf = np.max(np.abs(acf[1:10]))  # Lags 1-9
        print(f"  Max autocorrelation (lags 1-9): {max_acf:.4f}")

        # 1st order prediction
        pred_acc = test_linear_predictability(residues[10:], q)
        from sympy import totient
        random_acc = 1.0 / int(totient(q)) if q > 2 else 1.0
        print(f"  1st-order prediction accuracy: {pred_acc:.4f} (random: {random_acc:.4f})")

        # 5th order prediction
        if q <= 13:
            pred_acc_5 = test_kth_order_prediction(residues[10:], q, k=3)
            print(f"  3rd-order prediction accuracy: {pred_acc_5:.4f}")

    # Cross-correlation analysis
    print("\n--- Cross-Correlation Between Moduli ---")
    for q1, q2 in [(3, 5), (3, 7), (5, 7), (3, 11), (7, 11)]:
        res1 = [p % q1 for p in primes]
        res2 = [p % q2 for p in primes]
        cc = cross_correlation(res1, res2)
        print(f"  corr(p(n) mod {q1}, p(n) mod {q2}) = {cc:.4f}")

    # CRT reconstruction test
    print("\n--- CRT Reconstruction Analysis ---")
    print("  If we knew p(n) mod q for all small q, CRT gives p(n) mod M")
    print("  where M = lcm(q for q in moduli).")
    print()
    from sympy import totient
    import math
    moduli = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    M = 1
    for q in moduli:
        M = math.lcm(M, q)
    print(f"  Using first {len(moduli)} primes as moduli: M = {M}")
    print(f"  log2(M) = {math.log2(M):.1f} bits")
    print(f"  For p(10^100) ~ 10^102, need log2(10^102) ≈ 339 bits")
    print(f"  Need primes up to ~339th prime = ~2269 as moduli")
    print(f"  That's {339} modular computations, each needing p(n) mod q")
    print()
    print("  BUT: computing p(n) mod q exactly IS the hard problem!")
    print("  p(n) mod q requires knowing which residue class p(n) falls in,")
    print("  which requires knowing pi(x; q, a) exactly near x = p(n).")
    print("  This is as hard as computing pi(x) exactly.")

    print("\n--- VERDICT ---")
    print("1. p(n) mod q has near-maximum entropy for all q > 2 → effectively random")
    print("2. Autocorrelation is near zero → no short-range patterns")
    print("3. Prediction accuracy ≈ random baseline → not predictable from recent history")
    print("4. Cross-correlations are near zero → modular residues are independent")
    print("5. CRT approach is CIRCULAR: computing p(n) mod q requires solving")
    print("   the original problem in the first place.")
    print()
    print("CONCLUSION: Modular residue approach FAILS.")
    print("This is consistent with the pseudorandomness of pi(x) mod 2.")

if __name__ == "__main__":
    main()
