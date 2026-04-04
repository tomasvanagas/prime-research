#!/usr/bin/env python3
"""
Session 13 Experiment: Algebraic Normal Form of Prime Indicator over GF(2)

The prime indicator function χ_P: {0,1}^N → {0,1} maps an N-bit binary
representation of n to 1 if n is prime, 0 otherwise.

Questions:
1. What is the algebraic degree of χ_P over GF(2) (= max degree in ANF)?
2. How does the degree scale with N?
3. Does the ANF have any exploitable structure (few high-degree terms)?
4. What fraction of the 2^N ANF coefficients are nonzero?

Relevance: If algebraic degree is d, then:
- χ_P can be computed by a polynomial-size circuit of depth O(d)
- The SUM of χ_P over all inputs might be computable more efficiently
- Low degree would imply pseudorandomness properties that help counting

We know:
- PRIMES ∉ AC^0 implies degree ≥ ω(1)
- PRIMES ∉ AC^0[2] implies the function is not well-approximated by low-degree
  GF(2) polynomials (more precisely, correlation ≤ 2^{-n^{Ω(1)}} with deg-d for d = o(n^c))
"""

import numpy as np
from sympy import isprime
from itertools import combinations

def truth_table(N):
    """Generate truth table of prime indicator for N-bit numbers."""
    table = []
    for n in range(2**N):
        table.append(1 if isprime(n) else 0)
    return np.array(table, dtype=np.int64)

def compute_anf(f_table, N):
    """
    Compute the Algebraic Normal Form (ANF) over GF(2).

    The ANF of f: {0,1}^N → {0,1} is:
    f(x) = XOR_{S ⊆ [N]} a_S * PROD_{i ∈ S} x_i

    where a_S = XOR_{T ⊆ S} f(T) (Möbius inversion on the power set lattice).

    Uses the butterfly (Walsh-Hadamard-like) algorithm.
    """
    # a[i] will hold ANF coefficient for the monomial corresponding to bits of i
    a = f_table.copy() % 2

    for i in range(N):
        mask = 1 << i
        for j in range(2**N):
            if j & mask:
                a[j] = (a[j] ^ a[j ^ mask]) % 2

    return a

def analyze_anf(anf, N):
    """Analyze the structure of the ANF."""
    total_terms = np.sum(anf)

    # Degree distribution
    degree_counts = {}
    max_degree = 0

    for idx in range(2**N):
        if anf[idx]:
            deg = bin(idx).count('1')
            degree_counts[deg] = degree_counts.get(deg, 0) + 1
            max_degree = max(max_degree, deg)

    return max_degree, total_terms, degree_counts

def format_monomial(idx, N):
    """Format a monomial index as x_i * x_j * ..."""
    bits = []
    for i in range(N):
        if idx & (1 << i):
            bits.append(f"x{i}")
    if not bits:
        return "1"
    return "*".join(bits)

def experiment_anf_degree():
    """Main experiment: how does ANF degree scale with N?"""
    print("="*70)
    print("EXPERIMENT: Algebraic Normal Form of Prime Indicator over GF(2)")
    print("="*70)

    results = []

    for N in range(3, 21):
        f = truth_table(N)
        num_primes = np.sum(f)

        anf = compute_anf(f, N)
        max_deg, total_terms, deg_counts = analyze_anf(anf, N)

        # Sparsity
        total_possible = 2**N
        sparsity = total_terms / total_possible

        results.append((N, 2**N, num_primes, max_deg, total_terms, sparsity))

        print(f"\nN = {N} ({2**N} numbers, {num_primes} primes)")
        print(f"  ANF degree: {max_deg} (max possible: {N})")
        print(f"  Nonzero terms: {total_terms} / {total_possible} ({100*sparsity:.1f}%)")
        print(f"  Degree distribution: ", end="")
        for d in sorted(deg_counts.keys()):
            print(f"d={d}:{deg_counts[d]} ", end="")
        print()

        # For small N, show the actual polynomial
        if N <= 6:
            print(f"  ANF: ", end="")
            terms = []
            for idx in range(2**N):
                if anf[idx]:
                    terms.append(format_monomial(idx, N))
            if len(terms) <= 20:
                print(" ⊕ ".join(terms))
            else:
                print(f"{len(terms)} terms (too many to display)")

    # Scaling analysis
    print("\n" + "="*70)
    print("SCALING ANALYSIS")
    print("="*70)

    print(f"\n{'N':>3} | {'2^N':>8} | {'primes':>8} | {'ANF deg':>8} | {'deg/N':>6} | {'terms':>8} | {'sparsity':>10}")
    print("-"*70)
    for N, size, nprimes, deg, terms, sp in results:
        print(f"{N:3d} | {size:8d} | {nprimes:8d} | {deg:8d} | {deg/N:6.3f} | {terms:8d} | {sp:10.4f}")

    # Fit degree scaling
    Ns = np.array([r[0] for r in results if r[0] >= 5])
    degs = np.array([r[3] for r in results if r[0] >= 5])

    if len(Ns) > 2:
        # Try: degree = a * N + b
        coeffs = np.polyfit(Ns, degs, 1)
        print(f"\nLinear fit: degree ≈ {coeffs[0]:.3f} * N + {coeffs[1]:.3f}")

        # Try: degree = c * N^alpha
        log_N = np.log(Ns)
        log_d = np.log(degs)
        alpha, log_c = np.polyfit(log_N, log_d, 1)
        c = np.exp(log_c)
        print(f"Power fit: degree ≈ {c:.3f} * N^{alpha:.3f}")

    # Key question: is degree always N or N-1?
    full_degree_count = sum(1 for r in results if r[3] == r[0])
    print(f"\nFull degree (deg = N): {full_degree_count}/{len(results)} cases")

    # Check: fraction of ANF terms at each degree level
    print(f"\nFraction of ANF at highest degree:")
    for N in range(5, min(17, max(r[0] for r in results)+1)):
        f = truth_table(N)
        anf = compute_anf(f, N)
        _, _, deg_counts = analyze_anf(anf, N)
        if N in deg_counts:
            total = sum(deg_counts.values())
            print(f"  N={N}: {deg_counts.get(N,0)}/{total} terms at deg N "
                  f"({100*deg_counts.get(N,0)/total:.1f}%)")

def experiment_fourier_weight():
    """Analyze the Fourier weight distribution of prime indicator."""
    print("\n" + "="*70)
    print("EXPERIMENT: GF(2) Fourier (Walsh-Hadamard) Analysis")
    print("="*70)

    for N in range(4, 18):
        f = truth_table(N)
        # Convert to {-1,+1} basis: g = 1 - 2f
        g = 1 - 2 * f.astype(float)

        # Walsh-Hadamard transform
        wht = np.copy(g)
        for i in range(N):
            stride = 1 << i
            for j in range(2**N):
                if not (j & stride):
                    a, b = wht[j], wht[j | stride]
                    wht[j] = a + b
                    wht[j | stride] = a - b
        wht = wht / (2**N)  # Normalize

        # Fourier weight at each level
        level_weight = {}
        for idx in range(2**N):
            level = bin(idx).count('1')
            level_weight[level] = level_weight.get(level, 0.0) + wht[idx]**2

        total_weight = sum(level_weight.values())

        # Spectral concentration
        low_weight = sum(level_weight.get(d, 0) for d in range(min(5, N+1)))

        print(f"\nN={N:2d}: total Fourier weight = {total_weight:.4f}")
        print(f"  Weight at level 0 (bias): {level_weight.get(0,0):.6f}")
        print(f"  Weight at levels 0-4: {low_weight:.4f} ({100*low_weight/total_weight:.1f}%)")
        print(f"  Weight at levels ≥N/2: ", end="")
        high_weight = sum(level_weight.get(d, 0) for d in range(N//2, N+1))
        print(f"{high_weight:.4f} ({100*high_weight/total_weight:.1f}%)")

        # Average sensitivity
        sensitivity = 0
        for i in range(N):
            for x in range(2**N):
                if f[x] != f[x ^ (1 << i)]:
                    sensitivity += 1
        sensitivity /= (N * 2**N)
        print(f"  Average sensitivity: {sensitivity:.4f}")

if __name__ == "__main__":
    experiment_anf_degree()
    experiment_fourier_weight()

    print("\n" + "="*70)
    print("IMPLICATIONS FOR pi(x)")
    print("="*70)
    print("""
If ANF degree = Theta(N): prime indicator needs depth Omega(N) in GF(2) model.
This doesn't directly imply TC^0 lower bounds (different model), but suggests
the function has intrinsic complexity Omega(N).

If Fourier weight is concentrated at high levels: the function is "complex"
in the Boolean analysis sense — far from any junta or low-degree polynomial.

If average sensitivity is Theta(1): each input bit matters (primes depend
on all bits of the number, not just a few).

COMPARISON: For random functions, ANF degree = N and sensitivity = 1/2.
If primes match random, it confirms the barrier.
""")
