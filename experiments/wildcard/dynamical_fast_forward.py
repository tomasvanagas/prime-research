"""
Dynamical Fast-Forward: Can primes be characterized by a fast-forwardable system?

ANALOGY: Fibonacci numbers satisfy F(n) = [[1,1],[1,0]]^n, computable in O(log n)
via matrix exponentiation. Is there an analogous structure for primes?

IDEA 1: PRIME GAPS AS A DYNAMICAL SYSTEM
Let g(n) = p(n+1) - p(n). Can we find a map T such that g(n+1) = T(g(n), state)?
If T is linear (or low-degree polynomial), matrix exponentiation gives p(n) fast.

IDEA 2: BEATTY SEQUENCE CONNECTION
A Beatty sequence is B(α,β) = {⌊αn + β⌋ : n ≥ 1}.
Primes are NOT a Beatty sequence, but can they be expressed as a
COMBINATION of Beatty sequences? This connects to the three-distance theorem.

IDEA 3: SYMBOLIC DYNAMICS
Encode primes as a binary sequence: b(n) = 1 if n is prime, 0 otherwise.
This sequence has known entropy (by PNT, entropy → 0 as n → ∞).
Is there a SOFIC SHIFT or SUBSTITUTION SYSTEM that generates it?

IDEA 4: COLLATZ-LIKE ITERATION
Is there a simple recurrence that, when iterated, produces primes?
p(n+1) = f(p(n), n) for some computable f?
Wilson's theorem: p is prime iff (p-1)! ≡ -1 (mod p).
Can this be turned into an efficient iteration?

IDEA 5: MATRIX GROUP ORBIT
Consider the group SL(2, Z) acting on the upper half-plane.
Primes are connected to the spectrum of the Laplacian on Γ\H.
Can we compute the nth eigenvalue fast? This is Selberg's connection.
"""

import numpy as np
import sympy
from sympy import prime, primepi, nextprime, isprime
import time

def test_gap_recurrence():
    """
    Test: Can prime gaps be predicted from previous gaps?

    Fit g(n+1) = f(g(n), g(n-1), ..., g(n-k)) for various f.
    """
    print("=== Prime gap recurrence analysis ===\n")

    N = 10000
    primes = [sympy.prime(i) for i in range(1, N+1)]
    gaps = [primes[i+1] - primes[i] for i in range(len(primes)-1)]

    # Test linear prediction: g(n+1) ≈ Σ a_i * g(n-i)
    for order in [1, 2, 3, 5, 10, 20]:
        # Build regression matrix
        X = []
        y = []
        for i in range(order, len(gaps)-1):
            X.append(gaps[i-order:i])
            y.append(gaps[i])
        X = np.array(X)
        y = np.array(y)

        # Least squares fit
        try:
            coeffs, residuals, rank, sv = np.linalg.lstsq(X, y, rcond=None)
            y_pred = X @ coeffs
            mse = np.mean((y - y_pred)**2)
            var_y = np.var(y)
            r_squared = 1 - mse / var_y

            print(f"  Linear order {order:>2}: R² = {r_squared:.6f}, "
                  f"MSE = {mse:.2f}, σ²(gaps) = {var_y:.2f}")
        except Exception as e:
            print(f"  Linear order {order}: failed ({e})")

    # Test: are gaps MARKOV? i.e., does g(n+1) depend only on g(n)?
    print("\nMarkov analysis of prime gaps:")
    gap_transitions = {}
    for i in range(len(gaps)-1):
        g_curr = gaps[i]
        g_next = gaps[i+1]
        if g_curr not in gap_transitions:
            gap_transitions[g_curr] = []
        gap_transitions[g_curr].append(g_next)

    print(f"  Distinct gap values: {len(gap_transitions)}")
    for g in [2, 4, 6, 8, 10, 12]:
        if g in gap_transitions:
            transitions = gap_transitions[g]
            mean_next = np.mean(transitions)
            std_next = np.std(transitions)
            print(f"  g={g:>3}: mean(next)={mean_next:.2f}, std(next)={std_next:.2f}, "
                  f"count={len(transitions)}")

def test_beatty_decomposition():
    """
    Can primes be expressed as a finite union of generalized Beatty sequences?

    B(α,β) = {⌊αn + β⌋ : n ≥ 1}

    If primes = B(α₁,β₁) ∪ B(α₂,β₂) ∪ ... ∪ B(αₖ,βₖ)
    then p(n) could be computed in O(log n) time via sorting/merging.

    By Ramaré-Rumely, primes mod q are equidistributed for coprime residues.
    The "Beatty" approximation uses α ≈ log(n) which varies with n.
    """
    print("\n=== Beatty sequence decomposition ===\n")

    N = 1000
    primes_list = [sympy.prime(i) for i in range(1, N+1)]

    # Fit: p(n) ≈ α*n*ln(n) + β*n + γ
    # From PNT: p(n) ~ n*ln(n), so p(n)/n - ln(n) should be roughly constant
    residuals = []
    for i, p in enumerate(primes_list):
        n = i + 1
        if n > 1:
            expected = n * np.log(n)
            residuals.append(p - expected)

    residuals = np.array(residuals)
    ns = np.arange(2, N+1)

    print(f"p(n) - n·ln(n) statistics:")
    print(f"  Mean: {residuals.mean():.4f}")
    print(f"  Std: {residuals.std():.4f}")
    print(f"  Min: {residuals.min():.4f}, Max: {residuals.max():.4f}")

    # Better fit: p(n) ≈ n*(ln(n) + ln(ln(n)) - 1) (Cipolla)
    residuals2 = []
    for i, p in enumerate(primes_list):
        n = i + 1
        if n > 2:
            expected = n * (np.log(n) + np.log(np.log(n)) - 1)
            residuals2.append(p - expected)

    residuals2 = np.array(residuals2)
    ns2 = np.arange(3, N+1)

    print(f"\np(n) - n·(ln(n) + ln(ln(n)) - 1) statistics:")
    print(f"  Mean: {residuals2.mean():.4f}")
    print(f"  Std: {residuals2.std():.4f}")

    # Key question: is the residual (after removing smooth part) compressible?
    # Apply DFT to residual
    res_fft = np.abs(np.fft.fft(residuals2))
    sorted_fft = sorted(res_fft, reverse=True)
    cumul = np.cumsum(np.array(sorted_fft)**2)
    total = cumul[-1]

    print(f"\nFourier analysis of Cipolla residual:")
    for frac in [0.5, 0.8, 0.9, 0.95]:
        k = np.searchsorted(cumul, frac * total) + 1
        print(f"  {frac*100:.0f}% energy: top {k}/{len(residuals2)} coefficients ({k/len(residuals2)*100:.1f}%)")

def test_substitution_system():
    """
    Can the prime indicator sequence be generated by a substitution system?

    A substitution system maps symbols to strings: 0 → 01, 1 → 10 (Thue-Morse).
    If primes have such structure, iteration gives exponential speedup.

    Test: compare prime indicator with known automatic sequences.
    """
    print("\n=== Substitution system analysis ===\n")

    N = 1000
    # Prime indicator
    prime_ind = np.array([1 if isprime(n) else 0 for n in range(N)])

    # Thue-Morse sequence
    tm = np.zeros(N, dtype=int)
    for n in range(N):
        tm[n] = bin(n).count('1') % 2

    # Rudin-Shapiro sequence
    rs = np.zeros(N, dtype=int)
    for n in range(N):
        b = bin(n)
        count_11 = sum(1 for i in range(len(b)-1) if b[i:i+2] == '11')
        rs[n] = count_11 % 2

    # Correlation with automatic sequences
    def correlation(a, b):
        a = a - a.mean()
        b = b - b.mean()
        return np.abs(np.correlate(a, b)[0]) / (np.linalg.norm(a) * np.linalg.norm(b))

    print(f"Correlation of prime indicator with:")
    print(f"  Thue-Morse: {correlation(prime_ind, tm):.6f}")
    print(f"  Rudin-Shapiro: {correlation(prime_ind, rs):.6f}")

    # Is the prime indicator k-automatic for any k?
    # A sequence is k-automatic if the k-kernel is finite.
    # The k-kernel of (a_n) is {(a_{k^i*n+j})_{n≥0} : i≥0, 0≤j<k^i}.
    # Test: does the 2-kernel of the prime indicator stabilize?

    print(f"\n2-kernel analysis of prime indicator (first {N} terms):")
    kernel_sequences = set()

    for i in range(6):  # k^i up to 64
        ki = 2 ** i
        for j in range(ki):
            subseq = tuple(prime_ind[ki*n + j] for n in range(min(50, (N-j)//ki)))
            kernel_sequences.add(subseq)

    print(f"  Distinct kernel sequences (length 50): {len(kernel_sequences)}")
    print(f"  (If finite and small → sequence is 2-automatic)")
    print(f"  (Growing with depth → NOT automatic)")

    # Compare: for Thue-Morse, kernel should be small
    tm_kernel = set()
    for i in range(6):
        ki = 2 ** i
        for j in range(ki):
            subseq = tuple(tm[ki*n + j] for n in range(min(50, (N-j)//ki)))
            tm_kernel.add(subseq)

    print(f"  Thue-Morse kernel sequences: {len(tm_kernel)} (should be small)")

def test_matrix_exponentiation_prime_gaps():
    """
    WILD IDEA: What if there's a hidden linear recurrence in a LIFTED space?

    For Fibonacci, the 2D state (F_n, F_{n-1}) evolves linearly.
    For primes, maybe a high-dimensional state evolves (approximately) linearly?

    State = (p(n), p(n-1), ..., p(n-k), n, 1)
    Try: is there a matrix M such that state(n+1) ≈ M · state(n)?

    This is equivalent to: p(n+1) ≈ Σ a_i * p(n-i) + b*n + c
    i.e., p(n) satisfies a linear recurrence with polynomial coefficients.
    """
    print("\n=== Matrix exponentiation for prime sequence ===\n")

    N = 5000
    primes_list = [sympy.prime(i) for i in range(1, N+1)]

    for order in [2, 5, 10, 20, 50]:
        # Build regression: p(n+1) = Σ a_i*p(n-i+1) + b*n + c
        X = []
        y = []
        for i in range(order, len(primes_list)-1):
            features = list(primes_list[i-order+1:i+1]) + [i+1, 1]
            X.append(features)
            y.append(primes_list[i])

        X = np.array(X, dtype=float)
        y = np.array(y, dtype=float)

        coeffs, residuals, rank, sv = np.linalg.lstsq(X, y, rcond=None)
        y_pred = X @ coeffs
        mse = np.mean((y - y_pred)**2)
        rmse = np.sqrt(mse)

        # Test prediction accuracy
        max_err = np.max(np.abs(y - y_pred))

        print(f"  Order {order:>2}: RMSE = {rmse:.2f}, Max error = {max_err:.2f}")

    print("\n  (For exact computation, max error must be < gap/2 ≈ 1-5)")
    print("  (If RMSE doesn't shrink to <1, linear recurrence fails)")

def test_p_adic_structure():
    """
    p-ADIC ANGLE: In the p-adic integers Z_p, primes have specific structure.

    The sequence p(n) mod p^k for increasing k traces a path in Z_p.
    Is this path predictable?

    More specifically: does p(n) mod p^k have a pattern as a function of n?
    If so, we could compute p(n) mod p^k for many p and reconstruct via CRT.
    """
    print("\n=== p-adic structure of p(n) ===\n")

    N = 1000
    primes_list = [sympy.prime(i) for i in range(1, N+1)]

    for p in [2, 3, 5, 7]:
        for k in [1, 2, 3]:
            mod = p ** k
            residues = [pr % mod for pr in primes_list]

            # Analyze: is the sequence (p(n) mod p^k) periodic or structured?
            # Count frequency of each residue
            from collections import Counter
            freq = Counter(residues)

            # Expected: uniform over residues coprime to p
            coprime_residues = [r for r in range(mod) if r % p != 0]
            expected_count = N / len(coprime_residues) if coprime_residues else 0

            max_deviation = max(abs(freq.get(r, 0) - expected_count) for r in coprime_residues)
            rel_deviation = max_deviation / expected_count if expected_count > 0 else 0

            print(f"  p(n) mod {p}^{k} ({mod:>4}): {len(coprime_residues)} residue classes, "
                  f"max deviation from uniform: {rel_deviation:.4f}")

            # Test autocorrelation of residue sequence
            res_arr = np.array(residues, dtype=float)
            res_arr -= res_arr.mean()
            if np.std(res_arr) > 0:
                autocorr = np.correlate(res_arr[:200], res_arr[:200], mode='full')
                autocorr = autocorr[len(autocorr)//2:]
                autocorr /= autocorr[0]

                # Find significant autocorrelation peaks
                peaks = []
                for lag in range(1, min(50, len(autocorr))):
                    if abs(autocorr[lag]) > 0.1:
                        peaks.append((lag, autocorr[lag]))

                if peaks:
                    print(f"    Significant autocorrelation peaks: {peaks[:5]}")

if __name__ == "__main__":
    print("=" * 60)
    print("EXPERIMENT: Dynamical Fast-Forward Systems")
    print("=" * 60)

    test_gap_recurrence()
    test_beatty_decomposition()
    test_substitution_system()
    test_matrix_exponentiation_prime_gaps()
    test_p_adic_structure()
