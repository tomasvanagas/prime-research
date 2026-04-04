"""
Session 6: Number-Theoretic Hash / Perfect Mapping Approach

IDEA: Can we construct a perfect hash function h(n) = p(n)?

Approach 1: Chinese Remainder Theorem with precomputed tables
  Store p(n) mod q_i for a set of coprime moduli q_1,...,q_k.
  Then CRT reconstructs p(n) from these residues.

  The catch: we need to STORE p(n) mod q_i, which IS the answer.
  UNLESS: we can find a FORMULA for p(n) mod q_i.

Approach 2: Minimal perfect hash via number theory
  Construct a function using mod operations that maps n → p(n).
  E.g., p(n) = (a*n^2 + b*n + c) mod M for some large M?

Approach 3: Pisano-like periodicity
  Fibonacci numbers mod m are periodic (Pisano period).
  Are primes mod m periodic in any sense?
  If p(n) mod m has period π(m,n), can we exploit this?

Approach 4: Stern-Brocot / Calkin-Wilf tree
  Every positive rational appears exactly once in the Stern-Brocot tree.
  Can we construct a SIMILAR tree for primes?

Approach 5: Conway-like constant
  Conway showed 13 specific fractions generate all primes via FRACTRAN.
  Can we find a SHORTER program (fewer fractions) or faster iteration?
"""

import numpy as np
from sympy import prime, isprime, factorint, gcd
import time

def experiment_crt_residue_patterns():
    """
    For small moduli m, study p(n) mod m as a function of n.
    Is there ANY pattern that could predict p(n) mod m from n alone?
    """
    print("="*70)
    print("CRT RESIDUE PATTERNS")
    print("="*70)

    N = 10000
    primes_list = [int(prime(n)) for n in range(1, N + 1)]

    for m in [2, 3, 5, 6, 7, 10, 11, 13, 30]:
        residues = [p % m for p in primes_list]

        # Is there periodicity in p(n) mod m?
        # Check: is the sequence p(n) mod m periodic?
        for period in range(2, min(100, N//2)):
            match = sum(1 for i in range(period, min(N, 2*period))
                        if residues[i] == residues[i - period])
            if match == period:
                print(f"  p(n) mod {m}: PERIODIC with period {period}!")
                break
        else:
            # No exact periodicity. Check correlation with shifted versions.
            r = np.array(residues[:N-1], dtype=float)
            r_shifted = np.array(residues[1:N], dtype=float)
            corr = np.corrcoef(r, r_shifted)[0, 1]

            # Check if n mod k predicts p(n) mod m
            best_pred_acc = 0
            best_k = 0
            for k in range(2, 50):
                # For each n mod k value, what's the most common p(n) mod m?
                from collections import Counter
                pred_correct = 0
                for nk in range(k):
                    indices = [i for i in range(N) if i % k == nk]
                    if not indices:
                        continue
                    resids = [residues[i] for i in indices]
                    most_common = Counter(resids).most_common(1)[0]
                    pred_correct += most_common[1]
                acc = pred_correct / N
                if acc > best_pred_acc:
                    best_pred_acc = acc
                    best_k = k

            # Random baseline: 1/phi(m) for coprime residues
            coprime_count = sum(1 for r in range(m) if gcd(r, m) == 1)
            random_baseline = max(Counter(residues).values()) / N

            print(f"  p(n) mod {m:2d}: no periodicity, lag-1 corr={corr:.4f}, "
                  f"best n-mod-k pred: {100*best_pred_acc:.1f}% (k={best_k}), "
                  f"baseline={100*random_baseline:.1f}%")

def experiment_polynomial_hash():
    """
    Can we find a polynomial P(n) such that P(n) ≡ p(n) (mod M) for some M?
    If so, and if M > max(p(n)), then P(n) = p(n).
    """
    print("\n" + "="*70)
    print("POLYNOMIAL HASH")
    print("="*70)

    N = 1000
    primes_list = [int(prime(n)) for n in range(1, N + 1)]
    ns = np.arange(1, N + 1, dtype=float)

    # Try polynomial interpolation
    # For N points, a degree N-1 polynomial fits exactly.
    # But that's just a lookup table in disguise.
    # Key question: is there a LOW-degree polynomial that approximates well?

    for degree in [2, 3, 5, 10, 20, 50]:
        coeffs = np.polyfit(ns, primes_list, degree)
        pred = np.polyval(coeffs, ns)
        errors = np.array(primes_list) - pred
        rms = np.sqrt(np.mean(errors**2))
        max_err = np.max(np.abs(errors))
        exact = np.sum(np.abs(errors) < 0.5)

        # Test on new data
        test_ns = np.arange(N + 1, N + 101, dtype=float)
        test_primes = [int(prime(n)) for n in range(N + 1, N + 101)]
        test_pred = np.polyval(coeffs, test_ns)
        test_err = np.array(test_primes) - test_pred
        test_exact = np.sum(np.abs(test_err) < 0.5)

        print(f"  Degree {degree:3d}: train RMS={rms:.1f}, max_err={max_err:.1f}, "
              f"exact={exact}/{N} ({100*exact/N:.1f}%), "
              f"test_exact={test_exact}/100 ({test_exact}%)")

def experiment_fractran_optimization():
    """
    Conway's PRIMEGAME uses 14 fractions to generate primes.
    Can we find a SHORTER or FASTER program?

    FRACTRAN: start with n, repeatedly multiply by the first fraction
    in the list that gives an integer. Record powers of 2 that appear.

    Conway's fractions: 17/91, 78/85, 19/51, 23/38, 29/33, 77/29,
    95/23, 77/19, 1/17, 11/13, 13/11, 15/14, 15/2, 55/1

    Can we find alternatives?
    """
    print("\n" + "="*70)
    print("FRACTRAN PRIME GENERATION")
    print("="*70)

    # Conway's PRIMEGAME
    fractions = [
        (17, 91), (78, 85), (19, 51), (23, 38), (29, 33),
        (77, 29), (95, 23), (77, 19), (1, 17), (11, 13),
        (13, 11), (15, 14), (15, 2), (55, 1)
    ]

    def fractran_step(n, fracs):
        for num, den in fracs:
            if (n * num) % den == 0:
                return (n * num) // den
        return None  # Halts

    def count_steps_to_primes(fracs, max_steps=100000):
        """Run FRACTRAN and record which powers of 2 appear."""
        n = 2  # Start
        primes_found = []
        steps = 0
        while steps < max_steps and len(primes_found) < 10:
            n = fractran_step(n, fracs)
            if n is None:
                break
            steps += 1
            # Check if n is a power of 2
            if n > 1 and (n & (n - 1)) == 0:
                exp = int(np.log2(n))
                primes_found.append((exp, steps))
        return primes_found

    results = count_steps_to_primes(fractions, max_steps=500000)
    print(f"  Conway PRIMEGAME results (up to 500K steps):")
    for exp, steps in results:
        print(f"    2^{exp} = {2**exp} at step {steps} "
              f"({'PRIME' if isprime(exp) else 'not prime'}: {exp})")

    if results:
        print(f"\n  Steps between consecutive primes:")
        for i in range(1, len(results)):
            print(f"    p={results[i][0]}: {results[i][1] - results[i-1][1]} steps")

    # The key insight: FRACTRAN takes O(p^2) steps per prime.
    # For the 100th prime (541), that's ~300K steps. Way too slow.
    # No known shortcut to make FRACTRAN faster.
    print(f"\n  FRACTRAN complexity: O(p²) per prime → O(p(n)²) for nth prime")
    print(f"  For p(100)=541: ~{541**2} steps")
    print(f"  For p(10^6)≈1.5×10^7: ~{(1.5e7)**2:.0e} steps → INFEASIBLE")

def experiment_stern_brocot_primes():
    """
    The Stern-Brocot tree enumerates all positive rationals exactly once.
    Can we construct a similar TREE for primes?

    The Calkin-Wilf tree uses the map: a/b → a/(b-a) or (a-b)/b.
    Starting from 1/1, this visits all positive rationals.

    IDEA: Define a binary tree T where:
    - Root is p(1) = 2
    - Left child of node with value p: next prime?
    - The tree structure encodes the prime sequence

    More interesting: can we define a tree by a SIMPLE RULE
    such that reading the leaves gives primes?
    """
    print("\n" + "="*70)
    print("TREE-BASED PRIME ENUMERATION")
    print("="*70)

    # Instead of a tree, try a SIMPLE ITERATION
    # Can we find f such that iterating f from some seed gives primes?

    # Test: f(x) = x + log(x) + correction
    # This is motivated by: p(n+1) ≈ p(n) + ln(p(n))
    N = 100
    primes_list = [int(prime(n)) for n in range(1, N + 1)]

    # Test various iterative maps
    maps = {
        'x + ln(x)': lambda x: x + np.log(x),
        'x + ln(x) + 0.5/ln(x)': lambda x: x + np.log(x) + 0.5/np.log(x),
        'x + ln(x) + 1': lambda x: x + np.log(x) + 1,
        'x * (1 + 1/ln(x))': lambda x: x * (1 + 1/np.log(x)),
    }

    for name, f in maps.items():
        x = 2.0
        generated = [2]
        for _ in range(N - 1):
            x = f(x)
            generated.append(round(x))

        errors = [abs(g - p) for g, p in zip(generated, primes_list)]
        max_err = max(errors)
        rms = np.sqrt(np.mean(np.array(errors)**2))
        exact = sum(1 for e in errors if e == 0)

        print(f"  {name:30s}: exact={exact}/{N} ({100*exact/N:.0f}%), "
              f"RMS={rms:.1f}, max_err={max_err}")

    # Best possible: what if we use the CORRECT gap each time?
    # But predicting gaps is the whole problem.
    print(f"\n  For reference: average gap error using ln(p) as gap estimate:")
    gap_errors = [abs((primes_list[i+1] - primes_list[i]) - np.log(primes_list[i]))
                  for i in range(N-1)]
    print(f"    Mean: {np.mean(gap_errors):.2f}, Max: {max(gap_errors):.2f}")

def experiment_prime_in_balanced_system():
    """
    BALANCED NUMBER SYSTEM approach.

    In balanced ternary (-1, 0, 1), numbers have unique representations.
    What if we represent p(n) in a balanced system with base related to
    number-theoretic quantities?

    More specifically: can we represent the DIFFERENCE p(n) - round(n*ln(n))
    efficiently in some balanced system?
    """
    print("\n" + "="*70)
    print("BALANCED REPRESENTATION OF PRIME GAPS")
    print("="*70)

    N = 10000
    primes_list = [int(prime(n)) for n in range(1, N + 1)]

    # Compute residuals after Cipolla approximation
    residuals = []
    for n in range(1, N + 1):
        p = primes_list[n - 1]
        ln_n = np.log(n) if n > 1 else 1
        ln_ln_n = np.log(ln_n) if ln_n > 1 else 0.1
        approx = n * (ln_n + ln_ln_n - 1 + (ln_ln_n - 2) / ln_n) if n > 2 else p
        residuals.append(p - approx)

    residuals = np.array(residuals[2:])  # Skip n=1,2

    # Analyze the residual in various bases
    print(f"  Residual statistics (n=3..{N}):")
    print(f"    Mean: {np.mean(residuals):.2f}")
    print(f"    Std: {np.std(residuals):.2f}")
    print(f"    Entropy (assuming Gaussian): {0.5*np.log2(2*np.pi*np.e*np.var(residuals)):.2f} bits")

    # The entropy tells us how many bits per prime are truly "unpredictable"
    # For a Gaussian with std σ, entropy = 0.5 * log2(2πeσ²)
    # With σ ≈ 33, entropy ≈ 0.5 * log2(2π * e * 33²) ≈ 0.5 * log2(18707) ≈ 7.1 bits

    # But the residuals are CORRELATED (autocorrelation ~0.88 at lag 1)
    # So the conditional entropy h(δ_n | δ_{n-1}) is lower
    # Conditional std ≈ std * sqrt(1 - r²) ≈ 33 * sqrt(1 - 0.88²) ≈ 33 * 0.475 ≈ 15.7
    r1 = np.corrcoef(residuals[:-1], residuals[1:])[0, 1]
    cond_std = np.std(residuals) * np.sqrt(1 - r1**2)
    cond_entropy = 0.5 * np.log2(2 * np.pi * np.e * cond_std**2)

    print(f"    Lag-1 autocorrelation: {r1:.4f}")
    print(f"    Conditional std: {cond_std:.2f}")
    print(f"    Conditional entropy: {cond_entropy:.2f} bits per prime")

    # If we model as AR(1): δ_n = r1 * δ_{n-1} + ε_n
    # Then ε_n has std ≈ cond_std ≈ 15.7
    # These ε_n should be roughly iid with ~6.4 bits of entropy each
    # For n primes, total: ~6.4 * n bits beyond the smooth approximation

    # AR(1) prediction accuracy
    predicted = r1 * residuals[:-1]
    ar1_errors = residuals[1:] - predicted
    ar1_exact = np.sum(np.abs(ar1_errors) < 0.5)
    print(f"\n  AR(1) prediction of next residual:")
    print(f"    RMS error: {np.sqrt(np.mean(ar1_errors**2)):.2f}")
    print(f"    Exact: {ar1_exact}/{len(ar1_errors)} ({100*ar1_exact/len(ar1_errors):.1f}%)")

    # Higher order AR
    for order in [2, 5, 10]:
        X = np.column_stack([residuals[order-1-i:len(residuals)-1-i]
                             for i in range(order)])
        y = residuals[order:]
        coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        pred = X @ coeffs
        errs = y - pred
        exact = np.sum(np.abs(errs) < 0.5)
        entropy = 0.5 * np.log2(2 * np.pi * np.e * np.var(errs))
        print(f"    AR({order:2d}): RMS={np.sqrt(np.mean(errs**2)):.2f}, "
              f"exact={exact}/{len(errs)} ({100*exact/len(errs):.1f}%), "
              f"entropy={entropy:.2f} bits")

    print(f"\n  INFORMATION CONTENT SUMMARY:")
    print(f"    Total bits to specify p(n) for n=1..{N}: ~{N * np.log2(primes_list[-1]):.0f}")
    print(f"    Smooth approximation provides: ~{N * (np.log2(primes_list[-1]) - 7.1):.0f} bits")
    print(f"    Irreducible entropy: ~{N * cond_entropy:.0f} bits ({cond_entropy:.2f} bits/prime)")
    print(f"    For p(10^100): irreducible entropy per prime ~ {cond_entropy:.1f} bits")
    print(f"    That's ~10^{cond_entropy/3.32:.0f} possible values for the correction")

def main():
    print("Session 6: Number-Theoretic Hash Approaches")
    print(f"Started: {time.strftime('%H:%M:%S')}\n")

    t0 = time.time()

    experiment_crt_residue_patterns()
    experiment_polynomial_hash()
    experiment_fractran_optimization()
    experiment_stern_brocot_primes()
    experiment_prime_in_balanced_system()

    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"Total time: {elapsed:.1f}s")

if __name__ == '__main__':
    main()
