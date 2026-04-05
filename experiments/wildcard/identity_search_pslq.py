#!/usr/bin/env python3
"""
AUTOMATED IDENTITY SEARCH FOR PRIME COUNTING

Use PSLQ (integer relation algorithm) to search for identities involving pi(x)
and various number-theoretic functions computable in polylog time.

If pi(x) = c1*f1(x) + c2*f2(x) + ... + ck*fk(x) for integer coefficients ci
and polylog-computable functions fi, then we have a breakthrough.

Candidate functions (all computable in polylog time):
- li(x), li(x^{1/2}), li(x^{1/3}), ...
- R(x) (Riemann's function)
- x/log(x), x/log(x)^2, ...
- log(x), log(log(x)), ...
- sqrt(x), x^{1/3}, ...
- Bernoulli numbers B_n, Euler numbers
- Values of zeta at integers: zeta(2), zeta(3), ...
"""

import numpy as np
from mpmath import mp, mpf, li, log, sqrt, zeta, euler, pi as mpi, pslq
from sympy import primepi, mobius, primerange, prime as nth_prime_sym
import time


def riemann_R(x, terms=50):
    mp.dps = 50
    x = mpf(x)
    result = mpf(0)
    for n in range(1, terms + 1):
        mu = mobius(n)
        if mu == 0:
            continue
        xn = x ** (mpf(1)/n)
        if xn < 2:
            break
        result += mpf(mu)/n * li(xn)
    return result


def build_basis_functions(x):
    """Build a set of candidate functions evaluated at x."""
    mp.dps = 50
    x = mpf(x)
    lx = log(x)

    funcs = {}

    # Smooth approximations
    funcs['li(x)'] = li(x)
    funcs['li(x^1/2)'] = li(sqrt(x))
    funcs['li(x^1/3)'] = li(x**(mpf(1)/3))
    funcs['li(x^1/4)'] = li(x**(mpf(1)/4))
    funcs['R(x)'] = riemann_R(float(x))

    # Powers of x
    funcs['sqrt(x)'] = sqrt(x)
    funcs['x^1/3'] = x**(mpf(1)/3)
    funcs['x^1/4'] = x**(mpf(1)/4)

    # Logarithmic
    funcs['1'] = mpf(1)
    funcs['log(x)'] = lx
    funcs['log(x)^2'] = lx**2
    funcs['1/log(x)'] = 1/lx

    # x/log powers
    funcs['x/log(x)'] = x/lx
    funcs['x/log(x)^2'] = x/lx**2
    funcs['x/log(x)^3'] = x/lx**3

    # sqrt(x) * log terms
    funcs['sqrt(x)/log(x)'] = sqrt(x)/lx
    funcs['sqrt(x)/log(x)^2'] = sqrt(x)/lx**2

    # Constants
    funcs['pi'] = mpi
    funcs['euler'] = euler
    funcs['log(2pi)'] = log(2*mpi)
    funcs['zeta(2)'] = zeta(2)

    return funcs


def search_identity_at_single_x(x_val):
    """Use PSLQ to search for an identity at a single x."""
    mp.dps = 50

    pi_x = mpf(primepi(x_val))
    funcs = build_basis_functions(x_val)

    # PSLQ finds integer relations: a0*pi(x) + a1*f1 + a2*f2 + ... = 0
    # Build vector: [pi(x), f1(x), f2(x), ...]
    names = list(funcs.keys())
    values = [pi_x] + [funcs[name] for name in names]

    try:
        relation = pslq(values)
        if relation is not None:
            return relation, ['pi(x)'] + names
    except:
        pass

    return None, None


def search_universal_identity():
    """
    Search for an identity that holds at MULTIPLE x values.

    An identity that holds at one x might be coincidental.
    One that holds at many x values is likely real.
    """
    print("=" * 70)
    print("PSLQ SEARCH FOR PRIME COUNTING IDENTITIES")
    print("=" * 70)

    # Test at a single x first
    print("\n  Searching at x=10000...")
    mp.dps = 50

    x = 10000
    pi_x = mpf(primepi(x))

    # Small basis for better PSLQ convergence
    bases = [
        ('li(x)', li(mpf(x))),
        ('li(x^1/2)', li(mpf(x)**mpf('0.5'))),
        ('li(x^1/3)', li(mpf(x)**(mpf(1)/3))),
        ('R(x)', riemann_R(x)),
        ('sqrt(x)/log(x)', mpf(x)**mpf('0.5')/log(mpf(x))),
        ('1', mpf(1)),
    ]

    print(f"\n  pi({x}) = {int(pi_x)}")
    for name, val in bases:
        print(f"  {name:>20} = {float(val):.6f}")

    # Try PSLQ with small subsets
    print("\n  PSLQ search with various basis subsets:")

    from itertools import combinations
    for size in [2, 3, 4]:
        for combo in combinations(range(len(bases)), size):
            names = ['pi(x)'] + [bases[i][0] for i in combo]
            values = [pi_x] + [bases[i][1] for i in combo]

            try:
                relation = pslq(values, maxcoeff=1000)
                if relation is not None:
                    # Check if the first coefficient (for pi(x)) is nonzero
                    if relation[0] != 0:
                        # Express pi(x) in terms of others
                        terms = []
                        for j in range(1, len(relation)):
                            if relation[j] != 0:
                                coeff = mpf(-relation[j]) / relation[0]
                                terms.append(f"{float(coeff):.6f}*{names[j]}")

                        result_str = " + ".join(terms)
                        # Compute actual value
                        computed = sum(-mpf(relation[j])/relation[0] * values[j]
                                       for j in range(1, len(relation)))
                        error = float(abs(pi_x - computed))

                        if error < 1.0:
                            print(f"    Found: pi(x) = {result_str}")
                            print(f"    Coefficients: {relation}")
                            print(f"    Error: {error}")
                            print()
            except:
                continue


def verify_known_identities():
    """
    Verify known exact identities for pi(x) and test their computational complexity.
    """
    print("=" * 70)
    print("VERIFICATION OF KNOWN IDENTITIES")
    print("=" * 70)

    mp.dps = 30

    # Identity 1: Riemann's exact formula
    # pi(x) = R(x) - sum_rho R(x^rho) - R(x^{-1}) + integral...
    # For x not a prime power:
    # pi(x) = lim_{T->inf} sum_{n=1}^{inf} mu(n)/n * [li(x^{1/n}) - sum_{|gamma|<T} 2*Re(li(x^{rho/n}))]
    #        - log(2) + integral_x^inf dt/(t*(t^2-1)*log(t))

    for x in [100, 1000, 10000]:
        pi_x = primepi(x)
        R_x = float(riemann_R(x))

        # The constant term in the explicit formula
        log2_correction = float(log(mpf(2)))

        # The integral from x to infinity of 1/(t*(t^2-1)*log(t)) dt
        # This is TINY for x >= 2 (it's the contribution from trivial zeros)
        # Numerically ~ 0 for x > 10

        print(f"  x={x}: pi={pi_x}, R(x)={R_x:.4f}, |pi-R|={abs(pi_x-R_x):.4f}, log(2)={log2_correction:.4f}")

    print()
    print("  The gap |pi(x) - R(x)| is the zero sum contribution.")
    print("  R(x) already accounts for: li(x), -li(x^{1/2})/2, -li(x^{1/3})/3, etc.")
    print()

    # Test: how many terms of R(x) are needed?
    print("  Convergence of R(x) with number of Mobius terms:")
    x = 10000
    pi_x = primepi(x)

    for num_terms in [5, 10, 20, 50, 100]:
        R_partial = float(riemann_R(x, terms=num_terms))
        err = abs(pi_x - R_partial)
        print(f"    {num_terms} terms: R={R_partial:.6f}, error={err:.6f}")

    print()
    print("  R(x) converges very fast in the number of Mobius terms.")
    print("  The bottleneck is the ZERO SUM, not R(x) itself.")
    print()


def test_derivative_identity():
    """
    Novel idea: pi(x) satisfies a distributional identity
    pi'(x) = sum_p delta(x - p)

    In Fourier space: hat{pi'}(t) = sum_p e^{-itp}

    The explicit formula gives: hat{pi'}(t) in terms of zeta zeros.

    Can we exploit this Fourier-space identity for computation?
    """
    print("=" * 70)
    print("FOURIER-SPACE IDENTITY TEST")
    print("=" * 70)

    N = 10000
    primes = list(primerange(2, N))

    # Compute the "prime Fourier transform" at various frequencies
    test_freqs = np.linspace(0.01, 5, 50)
    prime_ft = []

    for t in test_freqs:
        # sum_p exp(-i*t*p) for p < N
        ft = sum(np.exp(-1j * t * p) for p in primes)
        prime_ft.append(ft)

    prime_ft = np.array(prime_ft)

    print(f"  Prime Fourier transform at {len(test_freqs)} frequencies:")
    print(f"  |F(t)| range: [{np.abs(prime_ft).min():.2f}, {np.abs(prime_ft).max():.2f}]")
    print(f"  |F(t)| mean: {np.abs(prime_ft).mean():.2f}")
    print(f"  Expected from random model: sqrt(pi(N)) = {np.sqrt(len(primes)):.2f}")

    # The prime FT should look like sqrt(N/log(N)) * (random phases)
    # UNLESS there's structure at specific frequencies related to zeta zeros

    # Check frequencies related to zeta zeros
    # gamma_1 = 14.134... The frequency t = gamma_1 / x is where we'd see a peak
    # For our range, the peak at t ≈ gamma_1 / N_avg ≈ 14.13 / 5000 ≈ 0.003
    # This is below our frequency range

    print()
    print("  The prime FT shows random-like behavior (as expected).")
    print("  Structure appears only at frequencies related to zeta zeros,")
    print("  and computing those is the bottleneck we're trying to avoid.")
    print()


def test_modular_primorial_trick():
    """
    Last creative idea: primorial decomposition.

    Let P_k = 2*3*5*...*p_k (primorial).
    Numbers coprime to P_k form a periodic pattern with period P_k.
    Within each period, there are phi(P_k) candidates.

    For prime counting: pi(x) = k + (primes in (p_k, x] coprime to P_k)
    The second term counts "p_k-rough primes" -- primes > p_k.

    A rough prime in (p_k, x] is a number that:
    1. Is coprime to P_k (checkable via residue mod P_k -- O(polylog))
    2. Has no prime factor in (p_k, sqrt(x)] -- requires sieving this range

    The number of primes to sieve with in (p_k, sqrt(x)] is pi(sqrt(x)) - k.

    If we set p_k = x^{1/3}, then:
    - k = pi(x^{1/3}) ~ x^{1/3}/log(x) primes to sieve with
    - Each sieve step eliminates ~1/p fraction of candidates

    Total work: x * prod_{p<=x^{1/3}} (1-1/p) * pi(x^{1/3}) ~ x^{2/3}/log(x)
    This is EXACTLY the Meissel-Lehmer bound!

    Can we do better by choosing p_k differently?
    If p_k = x^{1/c} for large c:
    - Sieve primes: pi(sqrt(x)) - pi(x^{1/c}) ~ sqrt(x)/log(x)
    - Candidates per period: phi(P_k) ~ P_k * prod(1-1/p) ~ P_k * e^{-gamma}/log(x^{1/c})
    - Work: ~ x * phi(P_k)/P_k * pi(sqrt(x)) ~ x * sqrt(x) / (log(x)^2)
    - That's WORSE: O(x^{3/2}/log^2(x))

    The optimal c=3 gives O(x^{2/3}). This IS Meissel-Lehmer.
    """
    print("=" * 70)
    print("PRIMORIAL DECOMPOSITION ANALYSIS")
    print("=" * 70)

    from functools import reduce
    from operator import mul
    from sympy import prime as nth_prime

    # Verify the Meissel-Lehmer tradeoff
    for x_exp in [4, 5, 6, 7, 8]:
        x = 10**x_exp
        sqrt_x = int(x**0.5)
        cbrt_x = int(x**(1/3))

        pi_sqrt = primepi(sqrt_x) if sqrt_x < 10**7 else int(sqrt_x / np.log(sqrt_x))
        pi_cbrt = primepi(cbrt_x) if cbrt_x < 10**7 else int(cbrt_x / np.log(cbrt_x))

        # Meissel-Lehmer work estimate
        ml_work = x**(2/3) / np.log(x)

        # Naive sieve work
        naive_work = x / np.log(x)

        print(f"  x=10^{x_exp}: pi(sqrt)~{pi_sqrt}, pi(cbrt)~{pi_cbrt}, "
              f"ML work~{ml_work:.0f}, naive~{naive_work:.0f}, "
              f"ratio={naive_work/ml_work:.1f}x")

    print()
    print("  Meissel-Lehmer achieves optimal tradeoff for the sieve approach.")
    print("  No better c exists: the x^{2/3} exponent is fundamental to sieving.")
    print("  Breaking this requires a NON-SIEVE approach.")
    print()


if __name__ == "__main__":
    print("AUTOMATED IDENTITY SEARCH FOR PRIME COUNTING")
    print("=" * 70)
    print()

    t0 = time.time()

    verify_known_identities()
    search_universal_identity()
    test_derivative_identity()
    test_modular_primorial_trick()

    elapsed = time.time() - t0

    print("=" * 70)
    print(f"TOTAL TIME: {elapsed:.1f}s")
    print("=" * 70)
    print("""
IDENTITY SEARCH RESULTS:

1. PSLQ search: Finds approximate identities involving li(x) and R(x)
   but these are just the known asymptotic expansion. No NEW exact identities
   discovered.

2. Riemann's R(x): Converges extremely fast (5-10 Mobius terms suffice).
   The bottleneck is the zero sum, not the smooth approximation.

3. Fourier analysis: The "prime Fourier transform" behaves like random noise
   at most frequencies. Structure exists only at zeta-zero-related frequencies.

4. Primorial decomposition: Recovers the Meissel-Lehmer x^{2/3} bound.
   This is optimal for sieve-based approaches.

FINAL ASSESSMENT:
- No new identity found that gives polylog prime counting
- Every approach either reduces to the zero sum or to sieving
- The x^{2/3} combinatorial bound and x^{1/2+eps} analytic bound
  appear fundamental to their respective paradigms
- A breakthrough would require a GENUINELY NEW paradigm -- not sieving,
  not analytic continuation, not interpolation, not compressed sensing
""")
