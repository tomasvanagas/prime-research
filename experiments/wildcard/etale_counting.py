"""
Étale Cohomology / Weil Conjectures Approach

THE DEEPEST ANALOGY:
The Weil conjectures (proved by Deligne) give EXACT point counts on
varieties over finite fields:

|V(F_q)| = Σ (-1)^i Tr(Frob | H^i_et(V))

The key: this is a FINITE sum of algebraic numbers.
For a curve of genus g, the sum has 2g + 2 terms.
Each term involves an eigenvalue of Frobenius.

For PRIME COUNTING, the "variety" is Spec(Z), and the "Frobenius"
traces are the Riemann zeta zeros. But Spec(Z) is 1-dimensional,
so naively the cohomology should be simple.

CONNES' APPROACH: The "space" of primes is a noncommutative space
(the adele class space). Its "zeta function" is the Riemann zeta.
The trace formula on this NC space gives the explicit formula for π(x).

QUESTION: Is there a FINITE-DIMENSIONAL approximation to this NC space
where the trace formula gives a POLYLOG-computable approximation to π(x)?

CONCRETE TEST:
- Compute point counts on curves over F_p using Weil formula
- Compare with direct counting
- Analyze: how does the "genus" affect computation time?
- Can we construct an "arithmetic curve" whose point count equals π(x)?
"""

import numpy as np
import sympy
from sympy import primepi, prime, nextprime, isprime, factorint
import math
import time

def weil_point_count(p, coefficients):
    """
    Compute |C(F_p)| for a hyperelliptic curve y² = f(x)
    using the Weil bound and character sums.

    For y² = f(x) over F_p, the count is:
    |C(F_p)| = p + 1 + Σ_{x=0}^{p-1} (f(x)/p)
    where (·/p) is the Legendre symbol.
    """
    count = 0
    for x in range(p):
        fx = sum(c * pow(x, i, p) for i, c in enumerate(coefficients)) % p
        if fx == 0:
            count += 1  # y = 0
        elif pow(fx, (p - 1) // 2, p) == 1:
            count += 2  # two y values
    # Add point at infinity (for odd degree)
    if len(coefficients) % 2 == 0:
        count += 1  # one point at infinity

    return count

def test_weil_computation():
    """
    Test: Weil formula gives exact point counts.
    This is what we WANT for primes: a finite trace formula.
    """
    print("=== Weil point counts on curves ===\n")

    # y² = x³ + ax + b (elliptic curves, genus 1)
    curves = [
        ([0, 0, 0, 1, 1], "y²=x³+x+1"),  # coeffs of 1 + x + 0x² + 0x³
        ([0, 1, 0, 1, 0], "y²=x³+x"),
        ([1, 0, 0, 1, 0], "y²=x³+1"),
    ]

    # Actually let's use correct format: f(x) = x³ + ax + b
    # coefficients[i] = coefficient of x^i

    print("Elliptic curve point counts vs Hasse-Weil bound:")
    for p in [5, 7, 11, 13, 17, 19, 23]:
        for coeffs, name in curves:
            N = weil_point_count(p, coeffs)
            a_p = p + 1 - N  # trace of Frobenius
            hasse_bound = 2 * math.sqrt(p)
            print(f"  {name}, F_{p}: |C| = {N}, a_p = {a_p}, "
                  f"|a_p| ≤ {hasse_bound:.2f}? {abs(a_p) <= hasse_bound}")

def test_arithmetic_surface_idea():
    """
    WILD IDEA: Can we construct a "variety" V_n such that |V_n(F_2)| = p(n)?

    If V_n has genus g_n and known Frobenius eigenvalues α_i,
    then |V_n(F_2)| = 2 + 1 - Σ α_i = p(n).

    This requires: Σ α_i = 3 - p(n), with |α_i| = √2 (Weil bound).

    For p(100) = 541: need Σ α_i = -538.
    Each |α_i| ≤ √2 ≈ 1.41, so need ≥ 538/1.41 ≈ 381 eigenvalues.
    So genus ≥ 191.

    For p(10^100): need genus ≈ p(10^100)/(2√2) ≈ 10^102.
    This is a HUGE genus, but the Weil formula still has
    only O(genus) terms. Can we evaluate the trace in O(polylog(genus))?

    The answer depends on the STRUCTURE of the eigenvalues.
    If they're given by a "nice" formula, yes.
    If they're random, no.
    """
    print("\n=== Arithmetic surface construction ===\n")

    # For a curve of genus g over F_q, the zeta function is:
    # Z(C, t) = P(t) / ((1-t)(1-qt))
    # where P(t) = Π (1 - α_i t) with |α_i| = √q

    # P(1) = |Jac(C)(F_q)| = class number
    # The class number is related to the L-function

    # For our purposes: can we find a FAMILY of curves C_n such that
    # |C_n(F_q)| = p(n)?

    # Simpler question: given a target N, find a curve over F_2 with N points
    for target_N in [5, 11, 23, 41, 541]:
        # Need: 2 + 1 - Σ α_i = target_N
        # So: Σ α_i = 3 - target_N
        total_trace = 3 - target_N
        min_genus = math.ceil(abs(total_trace) / (2 * math.sqrt(2)))

        print(f"  Target N = {target_N} (p({sympy.primepi(target_N) if isprime(target_N) else '?'})): "
              f"need trace = {total_trace}, min genus = {min_genus}")

    # The computation: given a curve of genus g over F_q,
    # computing |C(F_q)| via the Weil formula requires:
    # 1. Finding the 2g Frobenius eigenvalues (hard in general)
    # 2. Summing them (easy: O(g))
    #
    # Step 1 is the bottleneck. For a RANDOM curve, finding eigenvalues
    # requires O(g² log q) via Schoof-Elkies-Atkin type algorithms.
    # For SPECIAL curves (CM curves, modular curves), eigenvalues are
    # given by Hecke operators and may be faster.

    print("\n  Complexity analysis:")
    print("  - Finding Frobenius eigenvalues for genus g: O(g² log q) [Kedlaya]")
    print("  - If eigenvalues known: summing 2g terms in O(g)")
    print("  - For p(10^100): genus ≈ 10^102, so O(g²) = O(10^204) -- WORSE")
    print("  - UNLESS eigenvalues have special structure (CM, automorphic forms)")

def test_connes_trace_formula():
    """
    Connes' trace formula (noncommutative geometry approach):

    Σ_p Σ_m log(p) · f(m·log(p)) = f̂(0) - Σ_ρ f̂(ρ) + [other terms]

    Left side: sum over prime powers (the "geometric side")
    Right side: smooth term - sum over zeta zeros (the "spectral side")

    This is the explicit formula in disguise.
    The NC geometry reformulation doesn't help computationally UNLESS
    we can find a different test function f for which:
    1. The geometric side gives π(x) directly
    2. The spectral side has few terms

    TEST: For which f does the spectral side converge fastest?
    Gaussian f: f̂ is also Gaussian, fast decay means few zero terms needed.
    But narrow f in time ↔ wide f̂ in spectral ↔ more zeros needed.
    """
    print("\n=== Connes trace formula / test function optimization ===\n")

    # The explicit formula with a smooth test function h:
    # Σ_n Λ(n) h(log n) = ĥ(1) - Σ_ρ ĥ(ρ) + ...
    # where ĥ(s) = ∫ h(u) e^{su} du (Mellin-like transform)

    # For step function h(u) = 1_{u ≤ log x}: gives ψ(x) = Σ_{n≤x} Λ(n)
    # ĥ(s) = x^s / s
    # Spectral side: -Σ_ρ x^ρ/ρ which decays as 1/|ρ| ~ 1/γ_k

    # For Gaussian h(u) = exp(-u²/(2σ²)):
    # ĥ(s) = σ√(2π) exp(σ²s²/2)
    # At ρ = 1/2 + iγ: |ĥ(ρ)| = σ√(2π) exp(σ²/8) exp(-σ²γ²/2)
    # Decay: exp(-σ²γ²/2) which is FAST (Gaussian decay)!
    # But: this gives a SMOOTHED count, not sharp.

    from mpmath import mp, zetazero
    mp.dps = 20

    # Get first 50 zeta zeros
    gammas = [float(zetazero(k).imag) for k in range(1, 51)]

    print("Test function analysis: decay of spectral contributions")
    print(f"{'σ':>8} {'terms for 99%':>15} {'smoothing width':>15}")

    for sigma in [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
        # How many zeros needed for 99% of spectral sum?
        contributions = []
        for gamma in gammas:
            c = sigma * math.sqrt(2 * math.pi) * math.exp(sigma**2 / 8) * math.exp(-sigma**2 * gamma**2 / 2)
            contributions.append(abs(c))

        total = sum(contributions)
        cumulative = 0
        terms_99 = len(contributions)
        for k, c in enumerate(contributions):
            cumulative += c
            if cumulative > 0.99 * total:
                terms_99 = k + 1
                break

        # Smoothing width in x-space: the Gaussian averages over
        # an interval of width ~ σ around the target
        # In x-space: averaging over [x·e^{-σ}, x·e^{σ}]
        smoothing = f"x·e^{{±{sigma:.1f}}}"

        print(f"{sigma:>8.1f} {terms_99:>15} {smoothing:>15}")

    print("\n  Trade-off: narrow σ → many zeros needed (sharp count)")
    print("             wide σ → few zeros (Gaussian count) but smoothed")
    print("  For exact π(x): need σ → 0, which needs ALL zeros")
    print("  For rounding: need σ small enough that smooth ± correction < 0.5")

    # KEY QUESTION: what σ gives smoothed count within 0.5 of π(x)?
    print("\nSmoothed prime counting accuracy:")
    for sigma in [0.01, 0.05, 0.1, 0.5]:
        # Smoothed π(x) = Σ_n 1_{n prime} · exp(-(log(x/n))²/(2σ²)) / (σ√(2π)·n)
        # Roughly: weighted count of primes near x
        errors = []
        for x in [100, 500, 1000, 5000]:
            # Direct computation of smoothed count
            smoothed = 0
            for p in sympy.primerange(2, int(x * math.exp(3*sigma)) + 1):
                weight = math.exp(-(math.log(x/p))**2 / (2*sigma**2)) / (sigma * math.sqrt(2*math.pi))
                smoothed += weight

            actual = float(sympy.primepi(x))
            error = abs(smoothed - actual)
            errors.append(error)

        mean_err = np.mean(errors)
        max_err = max(errors)
        print(f"  σ={sigma:.2f}: mean_error={mean_err:.4f}, max_error={max_err:.4f}")

if __name__ == "__main__":
    print("=" * 60)
    print("EXPERIMENT: Étale Cohomology / Trace Formula Approach")
    print("=" * 60)

    test_weil_computation()
    test_arithmetic_surface_idea()
    test_connes_trace_formula()
