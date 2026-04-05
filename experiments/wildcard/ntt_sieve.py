"""
Number-Theoretic Transform (NTT) Sieve

BREAKTHROUGH ANALOGY: FFT turns O(N²) convolution into O(N log N).
The sieve of Eratosthenes IS a form of convolution — can NTT speed it up?

PRECISE FORMULATION:
Let χ(n) = prime indicator function.
Let μ(n) = Möbius function.

Then: χ(n) = Σ_{d|n} μ(d) · [n/d > 1]  ... not quite

More precisely, by Möbius inversion on the "sieve":
χ_P(n) = Σ_{d|n, d|P#} μ(d)  where P# = primorial

This is a DIRICHLET CONVOLUTION: χ_P = μ_P * 1
evaluated at n, where μ_P(d) = μ(d) if d|P#, 0 otherwise.

DIRICHLET CONVOLUTION can be computed via DIRICHLET SERIES MULTIPLICATION:
Σ f(n)/n^s · Σ g(n)/n^s = Σ (f*g)(n)/n^s

In the s-domain, this is POINTWISE multiplication!
So: transform to s-domain, multiply, transform back.

The "Dirichlet transform" is essentially the Mellin transform.
Computing it at enough points and inverting gives the convolution.

KEY QUESTION: How many "s-points" do we need?
If O(polylog), we have our algorithm.

ALSO: The PERRON FORMULA inverts Dirichlet series:
Σ_{n≤x} f(n) = (1/2πi) ∫ F(s) x^s/s ds

This integral can be approximated by sampling F(s) at O(log x) points
on the critical strip. But the error depends on how fast F(s) decays...
"""

import numpy as np
import sympy
from sympy import mobius, primepi, prime, divisors, factorint
import math
import time

def dirichlet_convolution_naive(f, g, N):
    """Compute (f * g)(n) = Σ_{d|n} f(d) · g(n/d) for n = 1..N"""
    result = np.zeros(N + 1)
    for n in range(1, N + 1):
        for d in range(1, n + 1):
            if n % d == 0:
                result[n] += f[d] * g[n // d]
    return result

def dirichlet_convolution_fast(f, g, N):
    """
    Faster Dirichlet convolution using the "hyperbola method".
    Still O(N log N) but constant factor better.
    """
    result = np.zeros(N + 1)
    for d in range(1, N + 1):
        if f[d] == 0:
            continue
        for k in range(1, N // d + 1):
            result[d * k] += f[d] * g[k]
    return result

def test_dirichlet_sieve():
    """
    Test: express the sieve as Dirichlet convolution and verify.
    """
    print("=== Dirichlet convolution sieve ===\n")

    N = 500

    # Möbius function
    mu = np.zeros(N + 1)
    for n in range(1, N + 1):
        mu[n] = int(sympy.mobius(n))

    # All-ones function
    ones = np.ones(N + 1)
    ones[0] = 0

    # μ * 1 = ε (Dirichlet identity: Σ_{d|n} μ(d) = [n=1])
    epsilon = dirichlet_convolution_fast(mu, ones, N)
    assert abs(epsilon[1] - 1) < 0.01, f"ε(1) = {epsilon[1]}"
    for n in range(2, min(N, 100)):
        assert abs(epsilon[n]) < 0.01, f"ε({n}) = {epsilon[n]}"
    print("  μ * 1 = ε: VERIFIED ✓")

    # For the sieve: we want χ(n) = prime indicator
    # χ(n) = 1 if n prime, 0 otherwise
    # There's no simple Dirichlet convolution that gives χ directly.
    # But: Λ(n) = -Σ_{d|n} μ(d)·log(d) (von Mangoldt function)
    # And: log(n) = Σ_{d|n} Λ(d)
    # So: Λ = -μ * log ... let's verify

    log_arr = np.zeros(N + 1)
    for n in range(1, N + 1):
        log_arr[n] = math.log(n)

    # μ' = μ · log
    mu_log = np.zeros(N + 1)
    for n in range(1, N + 1):
        mu_log[n] = mu[n] * math.log(n)

    # -μ' * 1 should give Λ
    neg_mu_log = -mu_log
    Lambda_computed = dirichlet_convolution_fast(neg_mu_log, ones, N)

    # Verify against known von Mangoldt
    errors = 0
    for n in range(2, min(N, 100)):
        # Λ(n) = log(p) if n = p^k, 0 otherwise
        f = sympy.factorint(n)
        if len(f) == 1:
            p, k = list(f.items())[0]
            expected = math.log(p)
        else:
            expected = 0
        if abs(Lambda_computed[n] - expected) > 0.01:
            errors += 1

    print(f"  Λ = -(μ·log) * 1: {100 - errors} / 98 correct")

    # Now: can we compute Σ_{n≤x} Λ(n) = ψ(x) via Dirichlet series?
    # ψ(x) = Σ_{n≤x} Λ(n) ≈ x (PNT)
    psi_x = np.cumsum(Lambda_computed[1:])
    print(f"\n  ψ(100) = {psi_x[99]:.4f} (expected ≈ 100)")
    print(f"  ψ(500) = {psi_x[499]:.4f} (expected ≈ 500)")

def test_perron_formula():
    """
    PERRON'S FORMULA:
    Σ_{n≤x} a(n) = (1/2πi) ∫_{c-i∞}^{c+i∞} F(s) · x^s/s ds

    where F(s) = Σ a(n)/n^s is the Dirichlet series.

    For a(n) = 1 (counting function):
    F(s) = ζ(s)
    Σ_{n≤x} 1 = ⌊x⌋ = (1/2πi) ∫ ζ(s) x^s/s ds

    For a(n) = Λ(n):
    F(s) = -ζ'(s)/ζ(s)
    ψ(x) = (1/2πi) ∫ [-ζ'(s)/ζ(s)] x^s/s ds

    Residues at poles give the explicit formula!
    Pole at s=1: residue = x (main term)
    Poles at zeros ρ: residue = -x^ρ/ρ (oscillatory terms)

    COMPUTATIONAL QUESTION:
    Can we evaluate this integral numerically with O(polylog) function evaluations?

    Using Gauss-Legendre or Clenshaw-Curtis quadrature on a vertical line:
    ∫_{c-iT}^{c+iT} F(s) x^s/s ds ≈ Σ w_k F(s_k) x^{s_k}/s_k

    Need T ≈ x for full accuracy, and O(T/log T) quadrature points.
    So: O(x / log x) evaluations of ζ(s) — WORSE than direct counting!

    BUT: each ζ(s) evaluation is O(|Im(s)|^{1/2} polylog) by Riemann-Siegel.
    Total: O(x · x^{1/2}) = O(x^{3/2}) — also worse.

    HOWEVER: what if we use EXPONENTIAL convergence quadrature?
    The integrand has poles (at zeros of ζ), not just smooth decay.
    With the RIGHT contour deformation, maybe fewer points suffice?
    """
    print("\n=== Perron formula numerical test ===\n")

    from mpmath import mp, zeta, mpf, mpc, pi as mpi, exp as mexp, log as mlog
    mp.dps = 20

    # Test: evaluate ψ(x) via truncated Perron integral
    def perron_psi(x, T, num_points=100):
        """Evaluate ψ(x) via Perron's formula with integration up to height T."""
        c = mpf('1.5')  # abscissa of convergence + margin
        x = mpf(x)

        total = mpf(0)
        # Trapezoidal rule on [c - iT, c + iT]
        dt = 2 * T / num_points
        for k in range(num_points + 1):
            t = -T + k * dt
            s = mpc(c, t)
            if abs(t) < 0.01:
                t_eff = mpc(0, 0.01)
                s = mpc(c, 0.01)

            # -ζ'(s)/ζ(s) · x^s / s
            try:
                zeta_val = zeta(s)
                if abs(zeta_val) < 1e-20:
                    continue
                zeta_prime = -zeta(s, derivative=1)  # ζ'(s)
                integrand = (zeta_prime / zeta_val) * mexp(s * mlog(x)) / s
                total += integrand * dt
            except:
                continue

        # Factor: 1/(2πi) · (result is purely real for real x)
        result = float(total.real / (2 * mpi))
        return result

    print("ψ(x) via Perron formula:")
    print(f"{'x':>8} {'T':>8} {'pts':>6} {'Perron':>12} {'Exact':>12} {'Error':>10}")

    for x in [50, 100]:
        # Compute exact ψ(x)
        psi_exact = sum(math.log(p) * k
                       for n in range(2, x + 1)
                       for p, k in sympy.factorint(n).items()
                       if len(sympy.factorint(n)) == 1)

        # Actually, Λ(n) = log p if n = p^k, else 0
        psi_exact = 0
        for n in range(2, x + 1):
            f = sympy.factorint(n)
            if len(f) == 1:
                p, k = list(f.items())[0]
                psi_exact += math.log(p)

        for T in [10, 50, 100]:
            for pts in [50, 200]:
                result = perron_psi(x, T, pts)
                error = abs(result - psi_exact)
                print(f"{x:>8} {T:>8} {pts:>6} {result:>12.4f} {psi_exact:>12.4f} {error:>10.4f}")

def test_multiplicative_fourier():
    """
    MULTIPLICATIVE FOURIER TRANSFORM:
    Instead of additive characters e^{2πinx/N}, use Dirichlet characters χ(n).

    The "Fourier transform" over (Z/qZ)* using Dirichlet characters:
    f̂(χ) = Σ f(n) χ(n)

    Inversion: f(n) = (1/φ(q)) Σ_χ f̂(χ) χ̄(n)

    For f = prime indicator, f̂(χ) = Σ_p χ(p) for p ≤ x.

    This connects to Dirichlet L-functions:
    L(s, χ) = Σ χ(n)/n^s = Π_p (1 - χ(p)/p^s)^{-1}

    KEY INSIGHT: If we could compute L(1, χ) for all χ mod q in O(polylog),
    we could reconstruct π(x; q, a) for all a, and thus π(x).

    But L(1, χ) involves all primes... circular again.

    UNLESS: there's a fast formula for L(1, χ) not involving primes!
    And there IS: L(1, χ) = -(1/q) Σ_{a=1}^{q} χ(a) ψ(a/q) for χ ≠ χ₀
    where ψ is the digamma function. This is POLYLOG computable!

    BUT: this gives L(1, χ), not Σ_{p≤x} χ(p).
    L(s, χ) = Σ χ(n)/n^s for ALL n, not just primes.
    """
    print("\n=== Multiplicative Fourier analysis ===\n")

    # Test: Dirichlet character decomposition of prime counting
    from sympy.ntheory import discrete_log

    q = 12  # modulus
    x = 1000

    # Compute π(x; q, a) for each residue class
    from collections import Counter
    residues = Counter()
    for p in sympy.primerange(2, x + 1):
        residues[p % q] += 1

    print(f"π({x}; {q}, a) for each a:")
    for a in sorted(residues.keys()):
        print(f"  a = {a:>3}: π = {residues[a]}")

    # Total check
    total = sum(residues.values())
    print(f"  Total: {total} (π({x}) = {sympy.primepi(x)})")

    # The point: to get π(x) from Dirichlet L-functions,
    # we still need to evaluate L(s, χ) along the critical strip
    # and sum over zeros. Same barrier as with ζ(s).

    print("\n  Conclusion: Multiplicative Fourier doesn't bypass the zero sum barrier.")
    print("  Each L(s, χ) has its own zeros, and we'd need to sum over all of them.")
    print("  The total number of zeros is still O(x^{1/2} log x).")

def test_ntt_convolution_sieve():
    """
    Can we speed up the Meissel-Lehmer computation using NTT?

    The key operation in Meissel-Lehmer is computing Φ(x, a).
    Φ involves floor divisions x//d for various d.

    OBSERVATION: {⌊x/n⌋ : n = 1..x} has only O(√x) distinct values.
    This is the "hyperbola method" insight.

    Can we use NTT to batch-compute all these floor divisions?
    Or to convolve the Möbius function with the counting function
    at all relevant points simultaneously?
    """
    print("\n=== NTT convolution for sieve ===\n")

    # Count distinct values of ⌊x/n⌋
    for x in [100, 1000, 10000, 100000, 1000000]:
        distinct = len(set(x // n for n in range(1, x + 1)))
        print(f"  x = {x:>10}: distinct ⌊x/n⌋ values = {distinct:>6} "
              f"(2√x = {2*int(x**0.5):>6})")

    # The O(√x) distinct values is the basis of the O(x^{2/3}) algorithm.
    # Can we push this further?

    # IDEA: Among the O(√x) distinct values, many are close together.
    # Can we interpolate between them?

    x = 10000
    vals = sorted(set(x // n for n in range(1, x + 1)), reverse=True)

    # Analyze gaps between consecutive distinct floor values
    gaps = [vals[i] - vals[i+1] for i in range(len(vals)-1)]
    print(f"\n  For x = {x}:")
    print(f"  Number of distinct values: {len(vals)}")
    print(f"  Max gap: {max(gaps)}")
    print(f"  Mean gap: {np.mean(gaps):.2f}")
    print(f"  Median gap: {np.median(gaps):.1f}")

    # Large gaps are near n = 1 (where ⌊x/n⌋ changes rapidly)
    # Small gaps are near n = √x (where ⌊x/n⌋ ≈ √x and changes by 1)

    # For the sieve: we need Σ μ(d) ⌊x/d⌋ for d squarefree
    # Group d values that give the SAME ⌊x/d⌋
    # Then: Σ μ(d) ⌊x/d⌋ = Σ_v v · (Σ_{d: ⌊x/d⌋=v} μ(d))
    # This groups the sum into O(√x) blocks

    # Can we compute Σ_{d: ⌊x/d⌋=v} μ(d) faster than iterating over d?
    # This is M(b) - M(a) where M(n) = Σ_{k≤n} μ(k) is the Mertens function
    # and [a+1, b] is the range of d with ⌊x/d⌋ = v.

    # M(n) can be computed in O(n^{2/3}) via Meissel-like recursion!
    # This gives: π(x) in O(x^{2/3}) which is the known best.

    # To beat this, we'd need M(n) in o(n^{2/3}), i.e., faster Mertens.
    # The Mertens function is conjectured to be hard (related to RH).

    print("\n  The floor value grouping gives O(√x) blocks.")
    print("  Each block needs Mertens function, computable in O(x^{2/3}).")
    print("  Total: O(x^{2/3}) — matches Meissel-Lehmer.")
    print("  Beating this requires faster Mertens computation.")

if __name__ == "__main__":
    print("=" * 60)
    print("EXPERIMENT: NTT Sieve / Dirichlet Convolution")
    print("=" * 60)

    test_dirichlet_sieve()
    test_perron_formula()
    test_multiplicative_fourier()
    test_ntt_convolution_sieve()
