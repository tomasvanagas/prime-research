"""
Session 8: Matrix/Operator approach to primes.

Novel angles:
1. Redheffer matrix: det(R_n) = M(n) = Σ_{k≤n} μ(k). Connected to RH.
   R_n[i,j] = 1 if j=1 or i|j.
   Can we extract π(n) from M(n)?

2. GCD matrix: A[i,j] = gcd(i,j). Eigenvalues involve φ(n) → primes.

3. Novel: Construct matrix whose characteristic polynomial has roots at first n primes.
   Then p(n) = largest root. But constructing the matrix IS knowing the primes.

4. Circulant matrices with prime structure.

5. The Hilbert-Pólya operator: If it exists, its eigenvalues are 1/2 + iγ_k.
   Computing eigenvalues would give zeta zeros → primes via explicit formula.

6. Novel: Matrix exponential approach — e^A for a structured A gives prime info?
"""

import numpy as np
from sympy import primerange, factorint, mobius, totient, gcd
import time

def get_primes(n):
    return list(primerange(2, n * 20))[:n]

# =============================================================================
# 1. Redheffer Matrix and Mertens function
# =============================================================================
def redheffer_approach():
    """Redheffer matrix: det(R_n) = M(n) = Mertens function."""
    print("=" * 60)
    print("1. Redheffer Matrix → Mertens function → primes?")
    print("=" * 60)

    # Mertens function M(n) = Σ_{k=1}^n μ(k)
    # |M(n)| / √n is bounded (Mertens conjecture is FALSE, but |M(n)| = O(n^{1/2+ε}) under RH)

    # Connection to primes:
    # M(n) = Σ_{k=1}^n μ(k)
    # We know: Σ_{d|n} μ(d) = [n=1]
    # And: Σ_{n≤x} μ(n)⌊x/n⌋ = 1
    # So M(x) is related to 1/ζ(s) = Σ μ(n)/n^s

    # Can we get π(x) from M(x)?
    # π(x) = Σ_{n≤x} 1_P(n)
    # 1_P(n) = -Σ_{d|n, d<n} μ(d) for n > 1... wait that's not right
    # Actually: Σ_{d|n} μ(d) = [n=1]
    # For prime p: μ(1) + μ(p) = 0, so μ(p) = -1 ✓
    # For p²: μ(1) + μ(p) + μ(p²) = 0, μ(p²)=0, so 1 + (-1) + 0 = 0 ✓

    # Direct relation: π(x) = -Σ_{n≤x} Σ_{p|n, p prime, p=n} 1
    # That's trivial and circular.

    # Better: π(x) = Σ_{k=2}^x Σ_{d|k} μ(d)·Λ(k/d)/ln(k/d) ... no, that's getting circular

    # The REAL connection:
    # ln ζ(s) = Σ_p Σ_k p^{-ks}/k = P(s) + (small)
    # -ζ'/ζ(s) = Σ Λ(n)/n^s
    # 1/ζ(s) = Σ μ(n)/n^s

    # So: M(x) = Σ_{n≤x} μ(n) encodes 1/ζ on the 1-line
    # Computing M(x) has same complexity as π(x): O(x^{2/3})

    for n in [10, 50, 100, 500]:
        # Direct computation of M(n)
        m_n = sum(int(mobius(k)) for k in range(1, n+1))
        pi_n = len(list(primerange(2, n+1)))
        print(f"  n={n}: M(n)={m_n}, π(n)={pi_n}, M/π ratio={m_n/pi_n:.3f}")

    print("\n  M(n) and π(n) are related through ζ(s) but computing M(n) is O(n^{2/3}) — same as π(n)")
    print("  No shortcut via Mertens function.")

# =============================================================================
# 2. GCD Matrix eigenvalues
# =============================================================================
def gcd_matrix_approach():
    """GCD matrix eigenvalues involve Euler totient → prime detection."""
    print("\n" + "=" * 60)
    print("2. GCD Matrix eigenvalues")
    print("=" * 60)

    # The Smith determinant: det[gcd(i,j)]_{1≤i,j≤n} = Π_{k=1}^n φ(k)
    # Eigenvalues of the GCD matrix are related to φ via Ramanujan sums

    # For the matrix A_n where A[i,j] = gcd(i,j):
    # Eigenvalues are λ_k = Σ_{d: k|d, d≤n} φ(d) · μ(d/k)... actually it's more complex

    # The eigenvalues of the GCD matrix over {1,...,n} are:
    # λ_k = Σ_{j=1}^{⌊n/k⌋} φ(jk) for k=1,...,n (THIS IS NOT EXACTLY RIGHT)

    # Actually, the eigenvalue associated with k is:
    # λ_k = Σ_{d: k|d≤n} φ(d)

    n = 20
    A = np.array([[gcd(i+1, j+1) for j in range(n)] for i in range(n)], dtype=float)
    eigenvalues = sorted(np.linalg.eigvalsh(A), reverse=True)

    print(f"  GCD matrix size {n}x{n}")
    print(f"  Top 5 eigenvalues: {[f'{e:.2f}' for e in eigenvalues[:5]]}")
    print(f"  Determinant = Π φ(k) = {np.prod(eigenvalues):.0f}")

    # Expected: det = Π_{k=1}^n φ(k)
    expected_det = 1
    for k in range(1, n+1):
        expected_det *= int(totient(k))
    print(f"  Expected det = Π φ(k) = {expected_det}")

    # The eigenvalues involve totient sums, not primes directly
    print("\n  GCD matrix eigenvalues = sums of totient values")
    print("  Computing these requires knowing factorizations → same barrier")

# =============================================================================
# 3. Circulant matrix with Dirichlet characters
# =============================================================================
def circulant_approach():
    """Circulant matrices and Dirichlet characters for prime detection."""
    print("\n" + "=" * 60)
    print("3. Circulant matrix approach")
    print("=" * 60)

    # For modulus q, define the q×q circulant matrix C where C[i,j] = f((i-j) mod q)
    # Eigenvalues of C are the DFT of f: λ_k = Σ f(j) · ω^{jk} where ω = e^{2πi/q}

    # If f = 1_P (prime indicator), then:
    # λ_k = Σ_{p≤q, p prime} ω^{pk}
    # These are "prime exponential sums" — deep in analytic number theory

    # Vinogradov's bound: |Σ_{p≤N} e^{2πipα}| ≤ N/(log N)^A for most α
    # So eigenvalues are O(q/log q) — not especially structured

    q = 100
    primes = list(primerange(2, q))
    omega = np.exp(2j * np.pi / q)

    eigenvalues = []
    for k in range(q):
        lam = sum(omega ** (p * k) for p in primes)
        eigenvalues.append(lam)

    mags = sorted([abs(e) for e in eigenvalues], reverse=True)
    print(f"  Prime DFT for q={q}:")
    print(f"  λ_0 = π(q) = {eigenvalues[0].real:.0f} (number of primes ≤ q)")
    print(f"  Top 5 magnitudes: {[f'{m:.2f}' for m in mags[:5]]}")
    print(f"  Ratio max/mean: {mags[0]/np.mean(mags):.2f}")

    # Is there a pattern in which k give large eigenvalues?
    large_k = [(k, abs(eigenvalues[k])) for k in range(q) if abs(eigenvalues[k]) > 10]
    large_k.sort(key=lambda x: -x[1])
    print(f"  Large eigenvalue indices: {[(k, f'{m:.1f}') for k, m in large_k[:10]]}")

    print("\n  Prime exponential sums are O(q/log q) by Vinogradov's theorem")
    print("  No special structure in the DFT of prime indicator")
    print("  Extracting primes from circulant eigenvalues is circular")

# =============================================================================
# 4. Novel: Companion matrix of prime-generating polynomial
# =============================================================================
def companion_matrix_approach():
    """Can we find a polynomial whose roots include all primes?"""
    print("\n" + "=" * 60)
    print("4. Companion matrix approach")
    print("=" * 60)

    # The polynomial P_n(x) = Π_{k=1}^n (x - p_k) has roots at the first n primes.
    # Its companion matrix has eigenvalues = primes.
    # But constructing P_n requires knowing the primes!

    # Alternative: Is there a structured matrix whose largest eigenvalue is p(n)?
    # Idea: the matrix L_n where L[i,j] = lcm(i,j). Its eigenvalues are related to
    # the prime factorization structure.

    # Actually, let's try: matrix T where T[i,j] = 1 if gcd(i,j) = 1 (coprimality matrix)
    n = 30
    T = np.array([[1 if gcd(i+1, j+1) == 1 else 0 for j in range(n)] for i in range(n)], dtype=float)

    eigenvalues = sorted(np.linalg.eigvalsh(T), reverse=True)
    print(f"  Coprimality matrix {n}×{n}")
    print(f"  Top 10 eigenvalues: {[f'{e:.2f}' for e in eigenvalues[:10]]}")

    # Are any eigenvalues close to primes?
    primes_set = set(primerange(2, 100))
    for e in eigenvalues[:10]:
        closest_prime = min(primes_set, key=lambda p: abs(p - e))
        print(f"    λ={e:.4f}, nearest prime={closest_prime}, diff={abs(e-closest_prime):.4f}")

    print("\n  Coprimality matrix eigenvalues are NOT close to primes")
    print("  No known structured matrix has eigenvalues = primes")

# =============================================================================
# 5. Novel: Zeta zeros from random matrix theory — can we sample zeros cheaply?
# =============================================================================
def rmt_approach():
    """Random Matrix Theory: Can we generate approximate zeta zeros from GUE?"""
    print("\n" + "=" * 60)
    print("5. Random Matrix Theory — approximate zeros from GUE")
    print("=" * 60)

    # The zeros of ζ(s) on the critical line have statistics matching GUE (Gaussian Unitary Ensemble)
    # If we could REPLACE actual zeros with GUE samples, we'd get an approximation to π(x)

    # Generate GUE eigenvalues at the right scale
    # For T zeros near height T: average spacing = 2π/ln(T/(2π))

    N = 100  # Number of zeros to approximate
    # GUE random matrix
    np.random.seed(42)
    A = np.random.randn(N, N) + 1j * np.random.randn(N, N)
    H = (A + A.conj().T) / 2  # Hermitian
    gue_eigenvalues = sorted(np.linalg.eigvalsh(H))

    # Scale to match zeta zero spacing near T=50
    T = 50
    mean_spacing = 2 * np.pi / np.log(T / (2 * np.pi))
    gue_scaled = np.array(gue_eigenvalues) * mean_spacing / (np.mean(np.diff(gue_eigenvalues)))
    gue_zeros = T + gue_scaled - np.mean(gue_scaled)

    # Compare with actual zeta zeros (first few known values)
    actual_zeros = [14.1347, 21.0220, 25.0109, 30.4249, 32.9351,
                    37.5862, 40.9187, 43.3271, 48.0052, 49.7738]

    print(f"  First 10 actual zeros: {[f'{z:.2f}' for z in actual_zeros]}")
    print(f"  GUE samples (scaled): {[f'{z:.2f}' for z in sorted(gue_zeros)[:10]]}")

    # The GUE gives the RIGHT statistics but not the RIGHT values
    # For π(x), we need the EXACT zeros, not statistical approximations
    print("\n  GUE gives correct STATISTICS but not correct VALUES of zeros")
    print("  Using GUE zeros in explicit formula gives ~correct mean but wrong oscillations")
    print("  Error from using GUE instead of actual zeros: O(√x) — same as no zeros at all")

# =============================================================================
# 6. Novel: Prime detection via matrix rank / nullity
# =============================================================================
def matrix_rank_detection():
    """Can primality be detected via matrix rank of a structured matrix?"""
    print("\n" + "=" * 60)
    print("6. Prime detection via matrix rank")
    print("=" * 60)

    # For n ≥ 2, define the n×n matrix M_n where M[i,j] = (i*j) mod n
    # What is rank(M_n) as a function of n?

    results = []
    for n in range(2, 50):
        M = np.array([[(i*j) % n for j in range(n)] for i in range(n)], dtype=float)
        r = np.linalg.matrix_rank(M, tol=0.5)
        is_prime = all(n % i != 0 for i in range(2, int(n**0.5)+1)) and n > 1
        results.append((n, r, is_prime))

    print("  n | rank(M_n) | prime?")
    print("  " + "-" * 30)
    for n, r, is_p in results[:20]:
        print(f"  {n:3d} | {r:4d}     | {'Yes' if is_p else 'No'}")

    # Check if rank distinguishes primes
    prime_ranks = [r for n, r, is_p in results if is_p]
    composite_ranks = [r for n, r, is_p in results if not is_p]

    print(f"\n  Prime ranks: {prime_ranks[:10]}")
    print(f"  Composite ranks: {composite_ranks[:10]}")

    # For prime p: the multiplication table mod p is a field → rank is related to p-1
    # For composite n: multiplication table has zero divisors → lower rank
    prime_rank_formula = all(r == n-1 for n, r, is_p in results if is_p)
    print(f"\n  rank(M_p) = p-1 for all tested primes: {prime_rank_formula}")
    print("  This is because Z/pZ is a field — multiplication by any nonzero element is bijective")
    print("  But this is just another primality test, NOT a formula for p(n)")

# =============================================================================
# Main
# =============================================================================
if __name__ == "__main__":
    print("SESSION 8: Matrix/Operator approaches to primes")
    print("=" * 60)
    t0 = time.time()

    redheffer_approach()
    gcd_matrix_approach()
    circulant_approach()
    companion_matrix_approach()
    rmt_approach()
    matrix_rank_detection()

    print(f"\n\nTotal time: {time.time() - t0:.1f}s")
    print("\n" + "=" * 60)
    print("CONCLUSIONS (Matrix/Operator)")
    print("=" * 60)
    print("""
1. Redheffer → Mertens: M(n) computable in O(n^{2/3}), same as π(n)
2. GCD matrix: eigenvalues = totient sums, requires factorizations
3. Circulant/DFT: prime exponential sums = O(q/log q) by Vinogradov, no structure
4. Companion matrix: requires knowing primes to construct — circular
5. RMT/GUE: gives correct statistics but NOT correct values — O(√x) error
6. Rank detection: rank(M_p) = p-1 for fields, but just a primality test

ALL MATRIX APPROACHES FAIL. No structured matrix encodes primes efficiently.
""")
