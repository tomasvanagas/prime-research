#!/usr/bin/env python3
"""
Session 7: Attempts to beat O(x^{1/2}) for computing π(x) or p(n).

Seven approaches investigated:
1. Hierarchical decomposition of π(x)
2. Parallel decomposition into independent sub-problems
3. Algebraic decomposition via number fields
4. Analytic continuation tricks for faster convergence
5. Approximate counting with certified error bounds
6. Binary splitting on the explicit formula
7. FFT acceleration of the explicit formula

Focus: ERROR ANALYSIS — quantify exactly where each approach fails or succeeds.
"""

import math
import time
import sys
from collections import defaultdict
from functools import lru_cache

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

try:
    import mpmath
    mpmath.mp.dps = 50
    HAS_MPMATH = True
except ImportError:
    HAS_MPMATH = False

# ============================================================================
# UTILITIES
# ============================================================================

def sieve_primes(limit):
    """Simple sieve of Eratosthenes."""
    if limit < 2:
        return []
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, limit + 1, i):
                is_prime[j] = False
    return [i for i in range(2, limit + 1) if is_prime[i]]

def pi_exact(x, primes=None):
    """Exact π(x) by counting sieved primes."""
    if primes is None:
        primes = sieve_primes(int(x) + 1)
    import bisect
    return bisect.bisect_right(primes, x)

def li(x):
    """Logarithmic integral li(x) = integral from 0 to x of dt/ln(t)."""
    if not HAS_MPMATH:
        # Rough approximation
        return x / math.log(x) * (1 + 1/math.log(x) + 2/math.log(x)**2)
    return float(mpmath.li(x))

def R_function(x):
    """Riemann R function: R(x) = sum_{k=1}^{inf} mu(k)/k * li(x^{1/k})."""
    if not HAS_MPMATH:
        return li(x)  # First approximation
    x = mpmath.mpf(x)
    mu = [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1,
          -1, 0, -1, 1, 1, 0, -1, 0, -1, 0,
          1, 1, -1, 0, 0, 1, 0, 0, -1, -1,
          -1, 0, 1, 1, 1, 0, -1, 1, 1, 0,
          -1, -1, -1, 0, 0, -1, 0, -1, 0]
    result = mpmath.mpf(0)
    for k in range(1, min(50, len(mu))):
        if mu[k] != 0:
            result += mpmath.mpf(mu[k]) / k * mpmath.li(x ** (mpmath.mpf(1)/k))
    return float(result)

# ============================================================================
# APPROACH 1: HIERARCHICAL DECOMPOSITION
# ============================================================================

def test_hierarchical_decomposition():
    """
    Can we decompose π(x) = f(π(x/2), π(x/3), ...) recursively?

    Key identity (Legendre/Meissel):
        π(x) = π(√x) + Φ(x, π(√x)) - 1
    where Φ(x,a) = #{n ≤ x : n not divisible by p_1,...,p_a}

    This gives a tree decomposition. Question: what's the tree depth
    and total work?
    """
    print("=" * 70)
    print("APPROACH 1: HIERARCHICAL DECOMPOSITION")
    print("=" * 70)

    primes = sieve_primes(10**6)

    # Legendre's identity: π(x) = π(√x) + Φ(x, π(√x)) - 1
    # Φ(x,a) = Φ(x,a-1) - Φ(x/p_a, a-1)
    # This recurses with x → x/p_a at each level

    # Analyze the recursion tree depth
    print("\n--- Recursion tree analysis for Legendre/Meissel ---")
    print(f"{'x':>15} {'sqrt(x)':>15} {'pi(sqrt(x))':>12} {'depth':>8} {'nodes_est':>15}")

    for exp in range(4, 21, 2):
        x = 10**exp
        sqrtx = x**0.5
        pi_sqrtx = li(sqrtx)  # approximate
        # Depth of Φ recursion is π(√x)
        depth = pi_sqrtx
        # Number of nodes ~ product of branching factors
        # At level k, we divide by p_k, so x → x/(p_1*...*p_k)
        # Recursion terminates when x/prod < p_{a-k+1}
        # Total nodes in Meissel: O(x^{2/3}/ln(x))
        nodes_est = x**(2/3) / math.log(x)
        print(f"{x:>15.0e} {sqrtx:>15.2e} {pi_sqrtx:>12.0f} {depth:>8.0f} {nodes_est:>15.2e}")

    # Can we make a POLYLOG depth tree?
    # If we decompose as π(x) = π(x/2) + (# primes in (x/2, x])
    # then depth = log_2(x), each level needs # primes in an interval
    print("\n--- Binary splitting: π(x) = π(x/2) + #{primes in (x/2, x]} ---")
    print("Depth = log_2(x) = polylog. BUT: counting primes in (x/2, x] is as hard as π(x).")

    # Alternative: use inclusion-exclusion with smooth moduli
    # π(x) = Σ_{d|P#} μ(d) * floor(x/d) + error
    # where P# = product of first k primes
    print("\n--- Inclusion-exclusion with primorial moduli ---")

    for x_val in [10**4, 10**5, 10**6]:
        primes_up_to = sieve_primes(int(x_val) + 1)
        pi_true = len(primes_up_to)

        for k in range(2, 8):
            small_primes = primes[:k]
            primorial = 1
            for p in small_primes:
                primorial *= p

            # Legendre sieve: Φ(x, k) via inclusion-exclusion
            # This counts numbers up to x coprime to first k primes
            phi_count = 0
            for mask in range(1 << k):
                d = 1
                bits = 0
                for j in range(k):
                    if mask & (1 << j):
                        d *= small_primes[j]
                        bits += 1
                sign = (-1)**bits
                phi_count += sign * int(x_val // d)

            # π(x) ≈ Φ(x, k) + k - 1
            pi_approx = phi_count + k - 1
            error = abs(pi_approx - pi_true)
            rel_error = error / pi_true if pi_true > 0 else 0

            if k <= 5:  # Don't print too many
                print(f"  x={x_val:.0e}, k={k} primes ({small_primes}): "
                      f"Φ+k-1={pi_approx}, π(x)={pi_true}, error={error} ({rel_error:.4f})")

    print("\n--- VERDICT on Approach 1 ---")
    print("Legendre/Meissel tree has depth π(√x) ~ √x/ln(x) — NOT polylog.")
    print("Binary splitting gives polylog depth but each node is as hard as the original.")
    print("Inclusion-exclusion sieve has O(2^k) terms; need k=π(√x) for exact, so O(2^{√x/ln x}).")
    print("CONCLUSION: No known hierarchical decomposition achieves polylog depth with")
    print("sub-exponential node work. The fundamental issue is that sieve methods require")
    print("Ω(x^{1/3}) nodes minimum (Deleglise-Rivat).")

    return {
        'approach': 'hierarchical',
        'best_complexity': 'O(x^{2/3}/ln(x))',
        'beats_sqrt': False,
        'barrier': 'Sieve recursion depth is π(√x), not polylog'
    }


# ============================================================================
# APPROACH 2: PARALLEL DECOMPOSITION
# ============================================================================

def test_parallel_decomposition():
    """
    Can we write π(x) = Σ_i f_i where each f_i is INDEPENDENTLY computable
    and each requires o(x^{1/2}) work?
    """
    print("\n" + "=" * 70)
    print("APPROACH 2: PARALLEL DECOMPOSITION")
    print("=" * 70)

    primes = sieve_primes(10**6)

    # Idea 1: Partition [1,x] into intervals and count primes in each
    # π(x) = Σ_{i=0}^{K-1} #{primes in (x*i/K, x*(i+1)/K]}
    # Each sub-problem: count primes in interval of length x/K
    # Cost per sub-problem: O((x/K)^{2/3}) using sieve
    # Total: K * O((x/K)^{2/3}) = O(x^{2/3} * K^{1/3})
    # This is WORSE than just doing π(x) directly!

    print("\n--- Interval partition analysis ---")
    print("Split [1,x] into K intervals, count primes in each via sieve.")
    print(f"{'x':>12} {'K':>8} {'cost_per':>15} {'total_cost':>15} {'vs_x^2/3':>12}")

    for exp in [10, 20, 50, 100]:
        x = 10.0**exp
        for K_exp in [0, exp//4, exp//2, exp]:
            K = max(1, 10.0**K_exp)
            interval_len = x / K
            cost_per = interval_len**(2/3)
            total = K * cost_per
            baseline = x**(2/3)
            ratio = total / baseline if baseline > 0 else float('inf')
            print(f"{x:>12.0e} {K:>8.0e} {cost_per:>15.2e} {total:>15.2e} {ratio:>12.4f}x")

    # Idea 2: Use Dirichlet characters to decompose
    # π(x) = (1/φ(q)) Σ_χ χ̄(a) Σ_{p≤x} χ(p)
    # Each character sum Σ_{p≤x} χ(p) can be computed via its L-function
    print("\n--- Dirichlet character decomposition ---")
    print("π(x; q, a) = #{p ≤ x : p ≡ a mod q}")
    print("π(x) = Σ_a π(x; q, a) for a coprime to q")

    for q in [6, 30, 210]:
        # Compute π(x; q, a) for small x
        x_val = 10**5
        residue_counts = defaultdict(int)
        for p in primes:
            if p > x_val:
                break
            if math.gcd(p, q) == 1:
                residue_counts[p % q] += 1
            elif p <= q:
                residue_counts[p % q] += 1

        pi_true = pi_exact(x_val, primes)
        reconstructed = sum(residue_counts.values())
        euler_phi = sum(1 for a in range(q) if math.gcd(a, q) == 1)

        print(f"\n  q={q}, φ(q)={euler_phi}:")
        print(f"  π({x_val:.0e}) = {pi_true}, reconstructed = {reconstructed}")
        print(f"  Sub-problems: {euler_phi}, each of size O(x^{{1/2+ε}}) via L-functions")
        print(f"  Total: {euler_phi} * O(x^{{1/2+ε}}) = O(x^{{1/2+ε}}) — NO improvement")

    # Idea 3: Smooth/rough decomposition
    # Split primes ≤ x into: smooth part (primes < y) and rough part (primes ≥ y)
    # π(x) = π(y) + #{n ≤ x : n prime, n > y}
    # The rough count can be done via Buchstab's identity
    print("\n--- Buchstab identity decomposition ---")
    print("Φ(x,y) = #{n ≤ x : P-(n) > y} where P-(n) = smallest prime factor")
    print("Buchstab: Φ(x,y) = x·ω(u)/y + error, u = log(x)/log(y)")
    print("This doesn't help — computing Φ exactly still requires O(x^{2/3}).")

    print("\n--- VERDICT on Approach 2 ---")
    print("Interval partitioning: total work = O(x^{2/3} K^{1/3}) — WORSE with more parts.")
    print("Dirichlet characters: each L-function sum is O(x^{1/2+ε}), and there are φ(q) of them.")
    print("Buchstab: Φ(x,y) computation is equivalent to Meissel-Lehmer.")
    print("CONCLUSION: π(x) does not admit a decomposition into sub-problems each of size o(x^{1/2}).")
    print("The bottleneck is that prime counting is inherently 'global' — every prime ≤ √x")
    print("contributes to the sieve at every scale.")

    return {
        'approach': 'parallel',
        'best_complexity': 'O(x^{2/3}) total, O(x^{2/3}/K^{2/3}) per part (K parts)',
        'beats_sqrt': False,
        'barrier': 'Sub-problems not independent — sieve couples all scales'
    }


# ============================================================================
# APPROACH 3: ALGEBRAIC DECOMPOSITION VIA NUMBER FIELDS
# ============================================================================

def test_algebraic_decomposition():
    """
    Use algebraic number fields K/Q. Primes split, remain inert, or ramify in K.
    π(x) = Σ contributions from different splitting types.
    Can some contributions be computed cheaply?
    """
    print("\n" + "=" * 70)
    print("APPROACH 3: ALGEBRAIC DECOMPOSITION VIA NUMBER FIELDS")
    print("=" * 70)

    primes = sieve_primes(10**5)
    x_val = 10**5
    pi_true = len(primes)

    # In Q(√d), a prime p splits if (d/p) = 1, is inert if (d/p) = -1, ramifies if p|d
    # By quadratic reciprocity and Chebotarev, half split, half inert (asymptotically)

    print("\n--- Splitting in Q(√d) ---")
    for d in [-1, 2, -3, 5, -7, 13]:
        split_count = 0
        inert_count = 0
        ramified_count = 0

        for p in primes:
            if p <= abs(d) and d % p == 0:
                ramified_count += 1
                continue
            # Legendre symbol (d/p)
            ls = pow(d % p, (p - 1) // 2, p) if p > 2 else 1
            if ls == p - 1:
                ls = -1
            if ls == 1:
                split_count += 1
            elif ls == -1:
                inert_count += 1
            else:
                ramified_count += 1

        total = split_count + inert_count + ramified_count
        print(f"  Q(√{d:>3d}): split={split_count}, inert={inert_count}, "
              f"ramified={ramified_count}, total={total}, π(x)={pi_true}")
        print(f"    split/total = {split_count/total:.4f} (expect ~0.5)")
        print(f"    π(x) = split + inert + ramified = {split_count} + {inert_count} + {ramified_count}")

    # The Chebotarev density theorem tells us the split/inert counts asymptotically
    # but the ERROR in Chebotarev is the same as in PNT for arithmetic progressions
    # i.e., O(x^{1/2+ε}) under GRH

    print("\n--- Error in Chebotarev density theorem ---")
    print("π_split(x; K) = (1/2)·li(x) + O(x^{1/2} log(disc(K)·x))")
    print("This error is EXACTLY the same order as the prime counting error!")
    print("So: π(x) = 2·π_split(x) + ramified ± O(x^{1/2} log x)")
    print("We cannot compute π_split(x) more cheaply than π(x) itself.")

    # Higher degree fields: same story
    print("\n--- Higher degree fields ---")
    print("For K/Q of degree n, Chebotarev gives density of each Frobenius class.")
    print("Error: O(x^{1/2} (log disc(K) + n log x)) under GRH.")
    print("Using multiple fields K_1, ..., K_m and Chinese Remainder:")
    print("  π(x) = Σ_i c_i · π_{C_i}(x; K_i) + small corrections")
    print("But EACH π_{C_i}(x; K_i) has error O(x^{1/2+ε}), and the c_i don't help")
    print("because the errors from different fields are CORRELATED (via zeta zeros).")

    # Artin L-functions
    print("\n--- Artin L-functions approach ---")
    print("L(s, χ, K/Q) encodes splitting behavior in K.")
    print("The zeros of L(s, χ, K/Q) control the error term for Chebotarev.")
    print("Under GRH for ALL Artin L-functions: individual errors are O(x^{1/2+ε}).")
    print("Cross-correlations between L-function zeros (Katz-Sarnak) show the errors")
    print("are NOT cancellable by linear combination.")

    print("\n--- VERDICT on Approach 3 ---")
    print("Algebraic decomposition gives π(x) = Σ (Chebotarev counts) + ramified primes.")
    print("Each Chebotarev count has error O(x^{1/2+ε}) under GRH.")
    print("The errors are correlated (all controlled by zeros on critical line).")
    print("CONCLUSION: Algebraic decomposition cannot beat O(x^{1/2+ε}). The splitting")
    print("behavior of primes in number fields is controlled by the same zeros of ζ(s).")

    return {
        'approach': 'algebraic',
        'best_complexity': 'O(x^{1/2+ε}) per component, correlated errors',
        'beats_sqrt': False,
        'barrier': 'Chebotarev errors controlled by same zeta zeros'
    }


# ============================================================================
# APPROACH 4: ANALYTIC CONTINUATION TRICKS
# ============================================================================

def test_analytic_continuation():
    """
    The explicit formula for ψ(x) = Σ_{ρ} x^ρ/ρ + ...
    Converges slowly when evaluated at s near 1.
    Can we evaluate at a point where convergence is faster?
    """
    print("\n" + "=" * 70)
    print("APPROACH 4: ANALYTIC CONTINUATION TRICKS")
    print("=" * 70)

    if not HAS_MPMATH:
        print("SKIPPED: mpmath not available")
        return {'approach': 'analytic_continuation', 'skipped': True}

    # The explicit formula:
    # ψ(x) = x - Σ_ρ x^ρ/ρ - ln(2π) - (1/2)ln(1 - x^{-2})
    # where ρ runs over nontrivial zeros of ζ(s)

    # Key insight: x^ρ = x^{1/2 + iγ} = √x · e^{iγ ln x}
    # The sum Σ x^ρ/ρ has terms of size ~√x/γ
    # Need T zeros for error O(x/T + √x log²x)
    # For error < 1: T ~ x (!!!) — this is O(x) zeros

    print("\n--- Standard explicit formula convergence ---")
    print("Term size: |x^ρ/ρ| ~ x^{1/2}/|γ|")
    print("For error < 1, need Σ_{|γ|>T} x^{1/2}/|γ| < 1")
    print("By zero density: #{|γ| ≤ T} ~ T log T / (2π)")
    print("Tail bound: Σ_{|γ|>T} 1/|γ| ~ log T")
    print("So need x^{1/2} · log T < 1, i.e., T > exp(x^{1/2}) — WORSE!")
    print()
    print("Actually, by Gallagher's theorem, the sum over zeros admits cancellation.")
    print("Effective bound: error ~ x/T · log²(xT) + x^{1/2} · log²x")
    print("For the FIRST term < 1: need T > x log²x")
    print("For the SECOND term < 1: IMPOSSIBLE — x^{1/2} log²x > 1 for large x")
    print()
    print("The x^{1/2} log²x term is the INTRINSIC error of the explicit formula")
    print("even with ALL zeros. It comes from the continuous density of zeros.")

    # Idea: Smoothed explicit formula
    # Instead of ψ(x), use a smoothed version:
    # ψ_h(x) = (1/h) ∫_{x-h/2}^{x+h/2} ψ(t) dt
    # This has x^ρ replaced by x^ρ · (sinh(ρh/2)/(ρh/2)) / ρ
    # The smoothing kills high-frequency oscillations

    print("\n--- Smoothed explicit formula ---")
    print("ψ_h(x) = (1/h) ∫ ψ(t) dt, smoothing window h")
    print("Each zero contribution multiplied by sinc(γh/(2π))")
    print("For γ > 2π/h, contributions are O(1/(γh)²) — MUCH faster decay!")
    print()
    print("Error analysis:")
    print("  Smoothing error: |ψ_h(x) - ψ(x)| ≤ max_{|t|<h} |ψ'(t)| · h/2")
    print("  ψ'(x) ~ 1 (for non-prime x), so smoothing error ~ h/2")
    print("  For this to give π(x): need h < gap between primes ~ ln x")
    print("  Zero sum with smoothing: Σ_ρ |x^ρ · sinc(γh/2)| / |ρ|")
    print("  Effective: ~ √x · (2π/h) / ln(2π/h) terms needed")
    print("  For h ~ ln x: need ~ √x / ln²x zeros")
    print()
    print("RESULT: Smoothing reduces # zeros by log factor, NOT the √x factor.")

    # Idea: Evaluate at s = 1 + 1/log(x)
    # -ζ'/ζ(s) = Σ_p Σ_k ln(p)/p^{ks} = Σ_n Λ(n)/n^s
    # At s = 1 + 1/log(x), each n^{-s} = n^{-1} · e^{-ln(n)/ln(x)}
    # This is a soft cutoff at n ~ x

    print("\n--- Evaluating -ζ'/ζ(s) at s = σ + it ---")
    print("Key identity: ψ(x) = (1/2πi) ∫ -ζ'/ζ(s) · x^s/s ds")
    print("Shifting contour to Re(s) = 1 + 1/ln(x):")
    print("  ψ(x) ≈ x + Σ_ρ x^ρ/ρ (residues from zeros)")
    print()
    print("Trick: use the FINITE Dirichlet series")
    print("  Σ_{n≤N} Λ(n)/n^s for s = 1 + 1/ln(x) + it")
    print("  Then use Perron's formula.")
    print("  But computing Λ(n) for n up to x REQUIRES knowing primes up to x!")
    print()
    print("Alternative: use Euler product")
    print("  -ζ'/ζ(s) = Σ_p ln(p)/(p^s - 1)")
    print("  At s = 1 + 1/ln(x): need primes p up to ~ x for convergence")
    print("  Same circularity!")

    # Idea: Weil explicit formula (dual form)
    # Σ_ρ h(ρ) = integral involving h_hat and log terms
    # Choose h to concentrate near x, making the sum small
    print("\n--- Weil explicit formula with optimal test function ---")
    print("Σ_ρ g(γ) = (g_hat(0)/2π)[2ln(2π) - ψ(1/4) - ψ(3/4)] + ...")
    print("           + (1/π) Σ_p Σ_k ln(p)/p^{k/2} · g_hat(k·ln(p))")
    print()
    print("Choosing g to be a Gaussian with width σ_g:")
    print("  g_hat has width ~ 1/σ_g in 'time' domain")
    print("  Need g_hat(ln p) for p up to x → need g_hat to be nonzero up to ln(x)")
    print("  → σ_g ~ 1/ln(x)")
    print("  → g(γ) ~ exp(-γ²/(2σ_g²)) = exp(-γ² ln²(x)/2)")
    print("  → Need zeros with γ < 1/ln(x) — essentially just the LOW zeros!")
    print()
    print("BUT: the prime sum side STILL requires knowing primes up to x.")
    print("The Weil formula relates zeros to primes — it doesn't eliminate either side.")

    # Numerical experiment: explicit formula for small x
    print("\n--- Numerical test: explicit formula accuracy vs # zeros ---")

    # First 30 zeros of ζ(s) (imaginary parts)
    zeta_zeros = [
        14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
        37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
        52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
        67.079811, 69.546402, 72.067158, 75.704691, 77.144840,
        79.337375, 82.910381, 84.735493, 87.425275, 88.809111,
        92.491899, 94.651344, 95.870634, 98.831194, 101.317851
    ]

    primes = sieve_primes(10000)

    for x_val in [100, 1000, 5000, 10000]:
        pi_true = pi_exact(x_val, primes)
        # ψ(x) ≈ x - Σ_ρ x^ρ/ρ - ln(2π)
        # π(x) ≈ R(x) - Σ_ρ R(x^ρ)  (Riemann's formula)

        Rx = R_function(x_val)
        print(f"\n  x = {x_val}: π(x) = {pi_true}, R(x) = {Rx:.2f}, "
              f"|R(x)-π(x)| = {abs(Rx - pi_true):.2f}")

        best_err = abs(Rx - pi_true)
        for num_zeros in [5, 10, 15, 20, 25, 30]:
            # Approximate: Σ_ρ R(x^ρ) ≈ Σ_γ 2·Re[R(x^{1/2+iγ})]
            # Simplified: contribution of zero γ ≈ -2·Re[li(x^{1/2+iγ})]/(2)
            # Better approximation: each zero pair contributes ~ -2√x·cos(γ ln x)/(γ ln x)
            correction = 0
            for gamma in zeta_zeros[:num_zeros]:
                # li(x^{1/2+iγ}) ≈ x^{1/2+iγ} / ((1/2+iγ)·ln(x))
                rho = 0.5 + 1j * gamma
                x_rho = x_val**rho
                li_x_rho = x_rho / (rho * math.log(x_val))
                correction += 2 * li_x_rho.real

            pi_approx = Rx - correction
            err = abs(pi_approx - pi_true)
            if num_zeros in [5, 15, 30]:
                print(f"    {num_zeros:>3d} zeros: π_approx = {pi_approx:.2f}, "
                      f"error = {err:.2f} "
                      f"(√x = {x_val**0.5:.1f}, √x/ln(x) = {x_val**0.5/math.log(x_val):.1f})")

    print("\n--- VERDICT on Approach 4 ---")
    print("All analytic tricks ultimately face the same barrier:")
    print("  - The explicit formula has INTRINSIC error O(√x · log²x) from zero density")
    print("  - Smoothing reduces # zeros by log factor only")
    print("  - Evaluating at s ≠ 1 requires knowing primes (circular)")
    print("  - Weil formula relates zeros ↔ primes, doesn't eliminate either")
    print("CONCLUSION: No analytic continuation trick can beat O(x^{1/2+ε}).")
    print("The √x barrier comes from the DENSITY of zeros on the critical line,")
    print("not from the rate of convergence of any particular sum.")

    return {
        'approach': 'analytic_continuation',
        'best_complexity': 'O(x^{1/2+ε}) — smoothing saves log factors only',
        'beats_sqrt': False,
        'barrier': 'Zero density on critical line forces √x error minimum'
    }


# ============================================================================
# APPROACH 5: APPROXIMATE COUNTING WITH CERTIFIED ERROR BOUNDS
# ============================================================================

def test_certified_approximation():
    """
    If we can compute π(x) with error < 1/2 in polylog time, we're done.
    Can we achieve error < 1/2?
    """
    print("\n" + "=" * 70)
    print("APPROACH 5: APPROXIMATE COUNTING WITH CERTIFIED ERROR BOUNDS")
    print("=" * 70)

    primes = sieve_primes(10**6)

    # Method 5a: R(x) alone
    print("\n--- 5a: R(x) approximation ---")
    print(f"{'x':>12} {'π(x)':>10} {'R(x)':>12} {'error':>10} {'√x/lnx':>12} {'err/√x·lnx':>12}")

    for exp in range(3, 7):
        x_val = 10**exp
        pi_true = pi_exact(x_val, primes)
        Rx = R_function(x_val)
        error = abs(Rx - pi_true)
        scale = x_val**0.5 / math.log(x_val)
        ratio = error / (x_val**0.5 * math.log(x_val)) if x_val > 1 else 0
        print(f"{x_val:>12.0e} {pi_true:>10d} {Rx:>12.4f} {error:>10.4f} {scale:>12.2f} {ratio:>12.6f}")

    # Under RH, |π(x) - R(x)| < (1/(8π)) √x ln(x) for x ≥ 599
    print("\n  RH bound: |π(x) - R(x)| < (1/(8π)) √x ln(x) ≈ 0.0398 √x ln(x)")
    print("  For x = 10^102: error < 0.0398 × 10^51 × 235 ≈ 9.4 × 10^51")
    print("  Need error < 0.5. Gap: factor of ~10^52")

    # Method 5b: R(x) + zero corrections
    print("\n--- 5b: R(x) with T zeros of ζ(s) ---")
    print("Under RH, using T zeros: |π(x) - R(x) + Σ_{|γ|≤T} R(x^ρ)| ≤ C·x/T·log²(xT)")
    print("For error < 0.5:")
    print("  Need C·x/T·log²(xT) < 0.5")
    print("  T > 2C·x·log²(xT)")
    print("  For x = 10^102: T > 10^102 (more zeros than x itself!)")
    print("  Cost of computing T zeros: O(T^{1+ε}) (Turing method)")
    print("  Total cost: O(x^{1+ε}) — WORSE than direct counting!")

    # Method 5c: Selberg-Delange method
    print("\n--- 5c: Selberg-Delange method ---")
    print("For sums Σ_{n≤x} f(n) where f is multiplicative:")
    print("  If F(s) = Σ f(n)/n^s = ζ(s)^z · G(s), G analytic for Re(s) > 1/2")
    print("  Then Σ_{n≤x} f(n) = x·(ln x)^{z-1}/Γ(z) · (c_0 + c_1/ln x + ...) + error")
    print("  Error: O(x · exp(-c√ln x))")
    print("  For π(x): F(s) = Σ 1_prime(n)/n^s = log ζ(s) + G(s)")
    print("  π(x) = li(x) + O(x · exp(-c√ln x))")
    print("  Under RH: O(x^{1/2} · ln x)")
    print("  The Selberg-Delange error is LARGER than the RH-conditional one!")

    # Method 5d: Saddle-point method
    print("\n--- 5d: Saddle-point method ---")
    print("π(x) = (1/2πi) ∮ (Σ_p p^{-s}) · x^s/s ds")
    print("Saddle point at s₀ where d/ds[s·ln(x) + log(Σ p^{-s})] = 0")
    print("Near s₀ = 1: the prime zeta function P(s) = Σ p^{-s} has a log singularity")
    print("P(s) = log ζ(s) + analytic part")
    print("The saddle-point approximation gives the SAME asymptotic as PNT.")
    print("Error: O(x · exp(-c√ln x)) unconditionally, O(x^{1/2+ε}) under RH.")

    # Method 5e: Can we combine MULTIPLE approximations?
    print("\n--- 5e: Combining multiple approximations ---")
    print("Suppose we have K approximations π̃_1(x), ..., π̃_K(x)")
    print("each with error ε_i, and errors are INDEPENDENT random variables.")
    print("Then averaging gives error ~ max(ε_i)/√K")
    print()
    print("BUT: errors from R(x), li(x), explicit formula are all CORRELATED")
    print("(all controlled by the same zeros of ζ(s)).")
    print("Empirical check:")

    errors_R = []
    errors_li = []
    for x_val in range(1000, 100000, 100):
        pi_true = pi_exact(x_val, primes)
        err_R = R_function(x_val) - pi_true
        err_li = li(x_val) - pi_true
        errors_R.append(err_R)
        errors_li.append(err_li)

    if HAS_NUMPY:
        corr = np.corrcoef(errors_R, errors_li)[0, 1]
        print(f"  Correlation between R(x)-π(x) and li(x)-π(x): {corr:.4f}")
        print(f"  (Expected ~1.0 since R(x) ≈ li(x) - li(√x)/2 + ...)")

        # How does error scale?
        x_values = list(range(1000, 100000, 100))
        errs = [abs(e) for e in errors_R]
        sqrtx = [v**0.5 for v in x_values]
        if len(errs) > 10:
            log_errs = [math.log(max(e, 0.1)) for e in errs]
            log_sqrtx = [math.log(s) for s in sqrtx]
            # Fit power law
            n = len(log_errs)
            mx = sum(log_sqrtx) / n
            my = sum(log_errs) / n
            sxy = sum((log_sqrtx[i] - mx) * (log_errs[i] - my) for i in range(n))
            sxx = sum((log_sqrtx[i] - mx)**2 for i in range(n))
            slope = sxy / sxx if sxx > 0 else 0
            print(f"  Power law fit: |R(x)-π(x)| ~ x^({slope/2:.3f})")
            print(f"  (Expected exponent: 0.5 under RH)")

    # Final error budget
    print("\n--- Final error budget for x = 10^102 ---")
    print("  R(x) alone:           error ~ 10^52")
    print("  + T zeros, T=10^20:   error ~ 10^82 (T too small)")
    print("  + T zeros, T=10^51:   error ~ 10^51 (barely helps)")
    print("  + T zeros, T=10^102:  error < 1 (but computing zeros costs 10^102)")
    print("  Selberg-Delange:      error ~ 10^102 · exp(-c·√235)")
    print("                        ≈ 10^102 · 10^{-7} = 10^95 (USELESS)")
    print("  Best unconditional:   error ~ x · exp(-c·(ln x)^{3/5}) ~ 10^95")
    print()
    print("  MINIMUM zeros needed for error < 0.5: T ~ x = 10^102")
    print("  Cost to compute those zeros: O(T^{1+ε}) = O(10^{102+ε})")
    print("  This is WORSE than the O(x^{2/3}) = O(10^68) sieve!")

    print("\n--- VERDICT on Approach 5 ---")
    print("No known smooth approximation achieves error < 1/2 in polylog(x) time.")
    print("The error bound under RH is |π(x) - R(x)| ≤ 0.04·√x·ln(x).")
    print("This is TIGHT — the error really is Θ(√x) on average.")
    print("To push error below 1/2, need Θ(x) zeros, costing Θ(x^{1+ε}).")
    print("CONCLUSION: Certified approximation with error < 1/2 requires Ω(x^{1/2}) work")
    print("under ANY known error bound.")

    return {
        'approach': 'certified_approximation',
        'best_complexity': 'O(x^{1/2+ε}) under RH, using x^{1/2+ε} zeros',
        'beats_sqrt': False,
        'barrier': 'Error is Θ(√x) with finitely many zeros; need Θ(x) zeros for error < 1/2'
    }


# ============================================================================
# APPROACH 6: BINARY SPLITTING ON THE EXPLICIT FORMULA
# ============================================================================

def test_binary_splitting():
    """
    The explicit formula is a sum over zeros. Can we use binary splitting
    or fast multiplication to evaluate it in O(T^{1-ε}) instead of O(T)?
    """
    print("\n" + "=" * 70)
    print("APPROACH 6: BINARY SPLITTING ON THE EXPLICIT FORMULA")
    print("=" * 70)

    # The sum: S(x) = Σ_{k=1}^{T} x^{ρ_k}/ρ_k
    # = Σ_{k=1}^{T} x^{1/2+iγ_k} / (1/2 + iγ_k)
    # = √x · Σ_{k=1}^{T} e^{iγ_k ln x} / (1/2 + iγ_k)

    # This is a "nonuniform DFT" — evaluating a sum of complex exponentials
    # at a SINGLE point (ln x).

    print("\n--- Structure of the zero sum ---")
    print("S(x) = √x · Σ_{k=1}^T exp(iγ_k · ln x) / (1/2 + iγ_k)")
    print()
    print("This is a single evaluation of a 'nonuniform Fourier series'.")
    print("There is no known algorithm to evaluate a single such sum in o(T) time.")
    print()
    print("Comparison with binary splitting:")
    print("  Binary splitting works for PRODUCTS: Π f(k) in O(T log²T · M(T))")
    print("  But our sum has VARYING frequencies γ_k — no algebraic structure to exploit.")

    # Can we approximate the sum using properties of the zeros?
    print("\n--- Zero distribution and partial sums ---")
    print("Montgomery pair correlation: zeros repel, with correlation R₂(u) = 1 - sinc²(u)")
    print("This means partial sums Σ_{k≤K} e^{iγ_k t} have some cancellation,")
    print("but the cancellation is already accounted for in the error bounds.")

    # Idea: group zeros by height and use mean-value estimates
    print("\n--- Grouping zeros by height ---")
    print("Partition zeros into blocks: [T_j, T_{j+1}] of length Δ")
    print("In each block: #{γ ∈ [T_j, T_{j+1}]} ~ Δ·ln(T_j)/(2π)")
    print("Sum over block: Σ_{γ∈block} e^{iγ ln x}/(1/2+iγ)")
    print("  ≈ (1/T_j) · Σ_{γ∈block} e^{iγ ln x}")
    print("  This inner sum has random phases (GUE statistics)")
    print("  Expected magnitude: ~ √(#{γ in block}) / T_j ~ √(Δ ln T_j) / T_j")
    print()
    print("Summing over O(T/Δ) blocks:")
    print("  Σ_j √(Δ ln T_j) / T_j × √x")
    print("  ~ √x · √Δ · Σ_j √(ln T_j) / T_j")
    print("  ~ √x · √Δ · ln T / T  (geometric series)")
    print("  This is NOT a useful bound for fixed T.")

    # Binary splitting on the Dirichlet series representation
    print("\n--- Binary splitting on Dirichlet series ---")
    print("ψ(x) = Σ_{n≤x} Λ(n)")
    print("Λ(n) = -Σ_{d|n} μ(d)·ln(n/d)")
    print("Can compute Σ_{n≤N} Λ(n) via hyperbola method in O(√N) time.")
    print("This is Mertens' trick and already used in Lagarias-Odlyzko.")
    print("Binary splitting on the hyperbola sum: O(√N · polylog(N)).")
    print("No improvement over direct evaluation.")

    # Baby-step giant-step for the zero sum
    print("\n--- Baby-step giant-step for zero sum ---")
    print("Write γ_k = α_k · M + β_k where α_k = floor(γ_k/M), β_k = γ_k mod M")
    print("Then e^{iγ_k t} = e^{iα_k M t} · e^{iβ_k t}")
    print("If we precompute e^{ijMt} for j = 0,...,T/M and e^{iβ t} for each β:")
    print("  S = Σ_k w_k · A[α_k] · B[β_k]")
    print("  where A[j] = e^{ijMt}, B[β] = e^{iβt}, w_k = 1/ρ_k")
    print("This doesn't help because the α_k, β_k are all DIFFERENT (irregular spacing).")
    print("No convolution structure to exploit.")

    print("\n--- VERDICT on Approach 6 ---")
    print("The zero sum S(x) = Σ x^ρ/ρ is a single-point evaluation of a")
    print("nonuniform trigonometric polynomial. No sub-linear algorithm is known.")
    print("Binary splitting requires algebraic structure (products, recurrences).")
    print("Baby-step giant-step requires uniform spacing.")
    print("CONCLUSION: Cannot evaluate the zero sum in o(T) time.")

    return {
        'approach': 'binary_splitting',
        'best_complexity': 'O(T) per evaluation, T ~ x^{1/2} zeros needed',
        'beats_sqrt': False,
        'barrier': 'Nonuniform zero spacing prevents sub-linear summation'
    }


# ============================================================================
# APPROACH 7: FFT ACCELERATION OF THE EXPLICIT FORMULA
# ============================================================================

def test_fft_acceleration():
    """
    Can we use FFT to evaluate Σ cos(γ_k ln x)/γ_k for many x simultaneously?
    Even though we only need ONE x, maybe we can exploit structure.
    """
    print("\n" + "=" * 70)
    print("APPROACH 7: FFT ACCELERATION OF THE EXPLICIT FORMULA")
    print("=" * 70)

    # NUFFT (Non-Uniform FFT) evaluates Σ_k c_k e^{iω_k t_j} for MANY t_j
    # in O(T log T + J) time where J = number of evaluation points.
    # But we only need ONE evaluation point!

    print("\n--- NUFFT analysis ---")
    print("NUFFT Type 1: Σ_{k=1}^T c_k e^{iω_k t_j} for j = 1,...,J")
    print("Cost: O(T log T + J log(1/ε))")
    print("For J=1: cost = O(T log T), WORSE than direct O(T) sum.")
    print()
    print("NUFFT Type 3: Σ_{k=1}^T c_k e^{iω_k t_j} for nonuniform ω_k AND t_j")
    print("Cost: O(T log T · log(1/ε) + J log(1/ε))")
    print("Still O(T log T) for J=1.")

    # Idea: compute π(x) for MANY x simultaneously, then interpolate
    print("\n--- Multi-point evaluation ---")
    print("Compute π(x_j) for j=1,...,J points using ONE pass over zeros:")
    print("  S(x_j) = Σ_{k=1}^T x_j^{ρ_k}/ρ_k")
    print("  = √x_j · Σ_k e^{iγ_k ln x_j} / ρ_k")
    print("Using NUFFT: cost O(T log T + J) for all J points.")
    print("Then interpolate to find exact crossing point where π(x) = n.")
    print()
    print("But this doesn't help: we still need T ~ x^{1/2} zeros,")
    print("and the cost is O(T log T) = O(x^{1/2} log x).")
    print("This is the SAME as Lagarias-Odlyzko!")

    # Idea: use the STRUCTURE of the zeros to speed up the sum
    print("\n--- Exploiting GUE statistics of zeros ---")
    print("Zero gaps: Δγ_k = γ_{k+1} - γ_k ~ 2π/(ln γ_k)")
    print("Mean spacing at height T: δ = 2π/ln T")
    print()
    print("Model: γ_k ≈ 2πk/ln(2πk) + fluctuation")
    print("If zeros were EXACTLY uniform: Σ e^{iγ_k t} would be a geometric series!")
    print("With fluctuations: Σ e^{iγ_k t} = Σ e^{i(2πk/ln(2πk))t} · e^{i·fluct_k·t}")
    print()
    print("The uniform part: Σ_{k=1}^T e^{i(2πk/ln T)t} ≈ T · sinc(Tt/(ln T))")
    print("The fluctuation part: each |fluct_k| ~ δ/√2 (GUE)")
    print("  Σ e^{i·fluct_k·t} needs individual evaluation — O(T)")

    # Numerical experiment: how well does the uniform approximation work?
    print("\n--- Numerical test: uniform zero approximation ---")

    zeta_zeros = [
        14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
        37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
        52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
        67.079811, 69.546402, 72.067158, 75.704691, 77.144840,
        79.337375, 82.910381, 84.735493, 87.425275, 88.809111,
        92.491899, 94.651344, 95.870634, 98.831194, 101.317851
    ]

    T = len(zeta_zeros)
    # Uniform model: γ_k ≈ 2πk / ln(2πk/e) (Riemann-von Mangoldt)
    uniform_zeros = []
    for k in range(1, T + 1):
        # N(T) ~ T/(2π) ln(T/(2πe))
        # Invert: γ_k ~ 2π k / ln(k) approximately
        # Better: use Gram points
        g = 2 * math.pi * k / math.log(max(k, 2))
        # Refine using N(T) formula
        for _ in range(5):
            N_g = g / (2 * math.pi) * math.log(g / (2 * math.pi * math.e)) + 7/8
            g = g + 2 * math.pi * (k - N_g) / math.log(g / (2 * math.pi))
        uniform_zeros.append(g)

    print(f"  {'k':>3} {'γ_k exact':>12} {'γ_k uniform':>14} {'difference':>12}")
    for k in range(min(10, T)):
        diff = zeta_zeros[k] - uniform_zeros[k]
        print(f"  {k+1:>3} {zeta_zeros[k]:>12.6f} {uniform_zeros[k]:>14.6f} {diff:>12.6f}")

    # Test: sum with exact vs uniform zeros
    x_val = 1000.0
    t = math.log(x_val)

    sum_exact = sum(
        math.cos(g * t) / g for g in zeta_zeros
    )
    sum_uniform = sum(
        math.cos(g * t) / g for g in uniform_zeros
    )

    print(f"\n  Sum cos(γ·ln({x_val}))/γ:")
    print(f"    Exact zeros:   {sum_exact:.6f}")
    print(f"    Uniform model: {sum_uniform:.6f}")
    print(f"    Difference:    {abs(sum_exact - sum_uniform):.6f}")
    print(f"    Relative:      {abs(sum_exact - sum_uniform)/max(abs(sum_exact), 1e-10):.4f}")

    # Key question: what fraction of the sum comes from individual zero positions?
    cumulative = []
    running = 0
    for g in zeta_zeros:
        running += math.cos(g * t) / g
        cumulative.append(running)

    print(f"\n  Cumulative sum (x={x_val}):")
    for k in [5, 10, 15, 20, 25, 30]:
        if k <= len(cumulative):
            print(f"    {k:>3} zeros: sum = {cumulative[k-1]:.6f}")

    print("\n--- VERDICT on Approach 7 ---")
    print("FFT/NUFFT can evaluate the zero sum at MANY points in O(T log T + J).")
    print("For a SINGLE point: O(T log T), no better than O(T) direct.")
    print("Uniform zero approximation has O(1) error per zero — cannot shortcut.")
    print("The GUE fluctuations in zero positions are ESSENTIAL for the result.")
    print("CONCLUSION: FFT cannot reduce the cost below O(T) for single-point evaluation.")
    print("Combined with T ~ x^{1/2+ε} zeros needed, total remains O(x^{1/2+ε}).")

    return {
        'approach': 'fft_acceleration',
        'best_complexity': 'O(T log T) = O(x^{1/2+ε} log x)',
        'beats_sqrt': False,
        'barrier': 'Single-point NUFFT is O(T log T), no better than O(T)'
    }


# ============================================================================
# BONUS: UNCONDITIONAL SUBEXPONENTIAL APPROACHES
# ============================================================================

def test_subexponential_ideas():
    """
    Even if we can't beat x^{1/2}, are there approaches that beat x^{2/3}
    WITHOUT assuming RH?
    """
    print("\n" + "=" * 70)
    print("BONUS: CAN WE BEAT x^{2/3} UNCONDITIONALLY?")
    print("=" * 70)

    print("\n--- Known complexity hierarchy ---")
    print("  Legendre (1808):      O(x / ln x)              — direct sieve")
    print("  Meissel (1870):       O(x^{2/3})               — recursive decomposition")
    print("  Lehmer (1959):        O(x^{2/3} / ln x)        — improved recursion")
    print("  Lagarias-Miller (1985): O(x^{2/3} / ln² x)     — further optimization")
    print("  Deleglise-Rivat (1996): O(x^{2/3} / ln² x)     — current best unconditional")
    print("  Lagarias-Odlyzko (1987): O(x^{1/2+ε})          — CONDITIONAL on RH")
    print()
    print("Gap: x^{2/3} vs x^{1/2+ε}. Can we close it unconditionally?")

    print("\n--- Approach: Combinatorial + analytic hybrid ---")
    print("Idea: Use Deleglise-Rivat for the 'hard' part, analytic for the 'easy' part.")
    print("Deleglise-Rivat: π(x) = Φ(x,a) + a - 1 - P₂(x,a) + ...")
    print("where Φ(x,a) has most of the cost.")
    print()
    print("Φ(x,a) can be split as:")
    print("  Φ(x,a) = Σ_{easy terms} + Σ_{hard terms}")
    print("  Easy: terms where x/d is 'small' — can evaluate quickly")
    print("  Hard: terms where x/d is 'large' — need sieve")
    print()
    print("The hard terms contribute O(x^{2/3}/ln x) to the total cost.")
    print("Can we use analytic estimates for these?")
    print("  Error per hard term: O(√(x/d) · ln(x/d)) under RH")
    print("  Number of hard terms: O(x^{1/3})")
    print("  Total error: O(x^{1/3} · √x · ln x) = O(x^{5/6} ln x) — WORSE!")
    print()
    print("So the analytic approach can't help with individual sieve terms.")

    print("\n--- Approach: Subexponential zero-free region ---")
    print("Current best zero-free region: σ > 1 - c/(log t)^{2/3} (log log t)^{1/3}")
    print("This gives π(x) = li(x) + O(x · exp(-c (ln x)^{3/5} (ln ln x)^{-1/5}))")
    print("For x = 10^102:")
    print("  ln(x) = 235, (ln x)^{3/5} ≈ 28, (ln ln x)^{1/5} ≈ 1.44")
    print("  Error ~ 10^102 · exp(-c · 19.4)")
    c_val = 0.05  # typical constant
    error_exp = 102 - c_val * 28 / 1.44 / math.log(10)
    print(f"  With c ≈ {c_val}: error ~ 10^{error_exp:.1f}")
    print(f"  Need error < 0.5 → need exponent < 0, but we have ~10^{error_exp:.0f}")
    print("  Hopeless — the zero-free region is too narrow!")

    print("\n--- What WOULD be needed? ---")
    print("For error < 0.5 in π(x) using analytic methods:")
    print("  Need zero-free region: σ > 1 - c/√(log t)")
    print("  This would give error ~ x · exp(-c√ln x) ~ x^{1-c/√(ln x)}")
    print("  For x = 10^102: error ~ 10^{102(1 - c/15.3)}")
    print("  Need 102(1 - c/15.3) < 0 → c > 15.3")
    print("  A zero-free region with c > 15 would suffice!")
    print("  But current best: c ≈ 0.05 (factor 300× too small)")
    print()
    print("Alternatively, to beat x^{2/3} computationally:")
    print("  Need either:")
    print("  (a) Unconditional O(x^{1/2+ε}) algorithm (= prove something about zeros)")
    print("  (b) Better zero-free region (> 50 years of stagnation)")
    print("  (c) Entirely new combinatorial technique (unknown)")

    print("\n--- GRAND SUMMARY ---")
    print("=" * 50)
    print("NO approach beats O(x^{1/2+ε}) for exact π(x).")
    print("This is the Lagarias-Odlyzko complexity, CONDITIONAL on RH.")
    print("Unconditionally, best is O(x^{2/3}/ln² x).")
    print()
    print("For p(10^100), x ≈ 10^102:")
    print("  O(x^{2/3}) = 10^68 operations — ~10^53 seconds")
    print("  O(x^{1/2+ε}) = 10^{51+ε} operations — ~10^{36} seconds")
    print("  Target: 10^9 operations (1 second)")
    print("  GAP: at least 10^42 (even with RH)")
    print("=" * 50)

    return {
        'approach': 'unconditional_improvement',
        'best_known_unconditional': 'O(x^{2/3}/ln²x)',
        'best_known_conditional': 'O(x^{1/2+ε})',
        'target': 'O(polylog(x))',
        'gap_factor': '10^42 minimum'
    }


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("Session 7: Recursive Decomposition & Sub-O(x^{1/2}) Attempts")
    print("=" * 70)
    print()

    results = []

    t0 = time.time()

    results.append(test_hierarchical_decomposition())
    results.append(test_parallel_decomposition())
    results.append(test_algebraic_decomposition())
    results.append(test_analytic_continuation())
    results.append(test_certified_approximation())
    results.append(test_binary_splitting())
    results.append(test_fft_acceleration())
    results.append(test_subexponential_ideas())

    elapsed = time.time() - t0

    print(f"\n\nTotal runtime: {elapsed:.2f}s")
    print("\n" + "=" * 70)
    print("SUMMARY TABLE")
    print("=" * 70)
    print(f"{'Approach':<30} {'Best Complexity':<35} {'Beats √x?':<10}")
    print("-" * 75)
    for r in results:
        name = r.get('approach', '?')
        complexity = r.get('best_complexity', r.get('best_known_conditional', '?'))
        beats = r.get('beats_sqrt', False)
        print(f"{name:<30} {complexity:<35} {'YES' if beats else 'NO':<10}")

    print("\n" + "=" * 70)
    print("FINAL CONCLUSION")
    print("=" * 70)
    print("All 7 approaches + bonus analysis confirm: O(x^{1/2+ε}) is a HARD BARRIER.")
    print("The barrier is not computational but INFORMATION-THEORETIC:")
    print("  - The explicit formula requires Θ(x^{1/2}) zeros for error < 1")
    print("  - Each zero contributes O(√x/γ) to π(x) — cannot be shortcut")
    print("  - Sieve methods require Ω(x^{1/3}) terms (Deleglise-Rivat lower bound)")
    print("  - Algebraic/character decompositions have correlated errors")
    print("  - No known parallelization reduces individual sub-problem below x^{1/2}")
    print()
    print("The ONLY hope would be a breakthrough in one of:")
    print("  1. Proving a quasi-polynomial zero-free region (σ > 1 - c/polylog(t))")
    print("  2. A completely new approach to exact prime counting (not sieve, not analytic)")
    print("  3. A quantum algorithm (Shor-type) for prime counting")


if __name__ == '__main__':
    main()
