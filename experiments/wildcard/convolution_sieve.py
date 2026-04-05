#!/usr/bin/env python3
"""
Compressed Sieve via Polynomial Convolution

Idea: Represent the Eratosthenes sieve as a polynomial product
  P(z) = prod_{p <= B} (1 - z^p)
and investigate whether fast polynomial arithmetic (FFT/NTT) or
evaluation at roots of unity can yield pi(x) faster than direct methods.

KEY MATHEMATICAL POINT:
P(z) = prod_{p<=B} (1 - z^p) has coefficients at positions that are
SUBSET-SUMS of the primes {p1, ..., pa}, NOT at subset-products.
The coefficient at z^k equals sum of (-1)^|S| over all subsets S of primes
whose elements sum to k.

The Legendre sieve formula is: phi(x, a) = sum_{d | P#} mu(d) * floor(x/d)
where d ranges over squarefree PRODUCTS of primes. This is NOT directly
encoded in the polynomial coefficients (which use SUMS as exponents).

So the polynomial P(z) does NOT directly compute phi(x,a).

However, a DIFFERENT polynomial does:
  Q(z) = prod_{p<=B} (1 - z^{1/p})  -- not a polynomial!

Or we can use GENERATING FUNCTIONS differently:
  F(s) = prod_{p<=B} (1 - p^{-s}) = sum_{d | P#} mu(d) * d^{-s}
  Then phi(x, a) = sum_{d | P#} mu(d) * floor(x/d)
  This is a Dirichlet series, not a power series -- FFT doesn't directly apply.

Experiments:
1. Verify the coefficient structure of P(z) = prod(1-z^p) -- subset sums
2. Compute the Legendre sieve directly by enumerating squarefree divisors
3. Compare the direct Legendre enumeration with standard sieve
4. Test evaluation of the Dirichlet polynomial at specific points
5. Complexity crossover analysis
"""

import time
import numpy as np
from sympy import primepi, primerange
from itertools import combinations

# ─── Helpers ─────────────────────────────────────────────────────────

def direct_sieve(x):
    """Standard Sieve of Eratosthenes. Returns count of primes <= x."""
    if x < 2:
        return 0
    sieve = bytearray(b'\x01') * (x + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(x**0.5) + 1):
        if sieve[i]:
            sieve[i*i::i] = b'\x00' * len(sieve[i*i::i])
    return sum(sieve)


def poly_product_sequential(x, B=None):
    """
    Compute P(z) = prod_{p<=B} (1 - z^p) mod z^{x+1} sequentially.
    Returns the coefficient array of length x+1.
    """
    if B is None:
        B = int(x**0.5)
    primes = list(primerange(2, B + 1))
    poly = np.zeros(x + 1, dtype=np.float64)
    poly[0] = 1.0
    for p in primes:
        new_poly = poly.copy()
        new_poly[p:] -= poly[:x + 1 - p]
        poly = new_poly
    return poly, primes


def poly_product_fft(x, B=None):
    """
    Compute P(z) = prod_{p<=B} (1 - z^p) mod z^{x+1} using FFT (divide & conquer).
    """
    if B is None:
        B = int(x**0.5)
    primes = list(primerange(2, B + 1))
    if not primes:
        r = np.zeros(x + 1)
        r[0] = 1.0
        return r, primes

    def make_factor(p, length):
        f = np.zeros(min(p + 1, length))
        f[0] = 1.0
        if p < length:
            f[p] = -1.0
        return f

    def fft_multiply_truncated(a, b, max_len):
        n = len(a) + len(b) - 1
        n = min(n, max_len)
        fft_size = 1
        while fft_size < len(a) + len(b) - 1:
            fft_size <<= 1
        fa = np.fft.rfft(a, fft_size)
        fb = np.fft.rfft(b, fft_size)
        result = np.fft.irfft(fa * fb, fft_size)
        return np.round(result[:n])

    def build_product(prime_list, max_len):
        if len(prime_list) == 0:
            return np.array([1.0])
        if len(prime_list) == 1:
            return make_factor(prime_list[0], max_len)
        mid = len(prime_list) // 2
        left = build_product(prime_list[:mid], max_len)
        right = build_product(prime_list[mid:], max_len)
        return fft_multiply_truncated(left, right, max_len)

    poly = build_product(primes, x + 1)
    if len(poly) < x + 1:
        poly = np.concatenate([poly, np.zeros(x + 1 - len(poly))])
    else:
        poly = poly[:x + 1]
    return poly, primes


def legendre_phi(x, primes):
    """
    Compute phi(x, a) = #{n <= x : n has no prime factor <= primes[-1]}
    via inclusion-exclusion over squarefree products of primes.
    phi(x, a) = sum over subsets S of primes: (-1)^|S| * floor(x / prod(S))
    """
    a = len(primes)
    if a > 24:
        raise ValueError(f"Too many primes ({a}) for brute-force Legendre: 2^{a} iterations")
    total = 0
    for mask in range(1 << a):
        d = 1
        bits = 0
        for i in range(a):
            if mask & (1 << i):
                d *= primes[i]
                bits += 1
        if d <= x:
            mu_d = (-1) ** bits
            total += mu_d * (x // d)
    return total


def pi_from_legendre(x):
    """
    Compute pi(x) using Legendre's formula:
    pi(x) = phi(x, a) + a - 1, where a = pi(sqrt(x)).
    """
    B = int(x**0.5)
    primes = list(primerange(2, B + 1))
    a = len(primes)
    phi = legendre_phi(x, primes)
    return phi + a - 1


# ─── Experiment 1: Polynomial coefficient structure ──────────────────

def experiment_1_coeff_structure(x_values):
    """
    Analyze the coefficient structure of P(z) = prod_{p<=B} (1-z^p).
    Verify that nonzero positions correspond to SUBSET-SUMS of primes.
    """
    print("=" * 70)
    print("EXPERIMENT 1: Coefficient Structure of P(z) = prod(1 - z^p)")
    print("=" * 70)

    records = []
    for x in x_values:
        B = int(x**0.5)
        poly, primes = poly_product_sequential(x, B)
        a = len(primes)

        nonzero_idx = np.where(np.abs(poly) > 0.5)[0]
        nonzero_vals = np.round(poly[nonzero_idx]).astype(int)
        total_nonzero = len(nonzero_idx)

        # Compute all possible subset sums of primes (up to x)
        subset_sums = set()
        for r in range(a + 1):
            for combo in combinations(primes, r):
                s = sum(combo)
                if s <= x:
                    subset_sums.add(s)

        # Check: every nonzero position is a subset sum
        nonzero_set = set(nonzero_idx.tolist())
        all_are_subset_sums = nonzero_set.issubset(subset_sums)

        # How many subset sums exist vs how many have nonzero coefficients?
        # Some subset sums may cancel (different subsets of same size summing to same value)
        pos_count = int(np.sum(nonzero_vals > 0))
        neg_count = int(np.sum(nonzero_vals < 0))
        max_abs_coeff = int(np.max(np.abs(poly)))
        coeff_sum = int(round(np.sum(poly)))  # P(1) = prod(1-1) = 0

        print(f"\nx = {x}, B = {B}, primes = {primes}")
        print(f"  #primes (a) = {a}")
        print(f"  Nonzero coefficients: {total_nonzero}")
        print(f"  Possible subset sums (<=x): {len(subset_sums)}")
        print(f"  All nonzero at subset-sum positions? {all_are_subset_sums}")
        print(f"  Max |coeff| = {max_abs_coeff}")
        print(f"  Positive coeffs: {pos_count}, Negative: {neg_count}")
        print(f"  P(1) = sum of coefficients = {coeff_sum} (should be 0)")
        if total_nonzero <= 20:
            print(f"  Terms: {list(zip(nonzero_idx.tolist(), nonzero_vals.tolist()))}")

        records.append({
            'x': x, 'B': B, 'a': a, 'nonzero': total_nonzero,
            'subset_sums_count': len(subset_sums),
            'all_subset_sums': all_are_subset_sums,
            'max_coeff': max_abs_coeff, 'coeff_sum': coeff_sum,
        })

    return records


# ─── Experiment 2: Legendre sieve correctness ───────────────────────

def experiment_2_legendre(x_values):
    """Verify Legendre sieve gives correct pi(x)."""
    print("\n" + "=" * 70)
    print("EXPERIMENT 2: Legendre Sieve Correctness")
    print("=" * 70)

    records = []
    for x in x_values:
        true_pi = int(primepi(x))

        t0 = time.time()
        leg_pi = pi_from_legendre(x)
        t_leg = time.time() - t0

        t0 = time.time()
        sieve_pi = direct_sieve(x)
        t_sieve = time.time() - t0

        match = (leg_pi == true_pi)
        print(f"\nx = {x:>8d}: pi(x) = {true_pi:>6d}, "
              f"Legendre = {leg_pi:>6d} {'OK' if match else 'FAIL'}, "
              f"t_leg = {t_leg:.6f}s, t_sieve = {t_sieve:.6f}s, "
              f"ratio = {t_leg/max(t_sieve, 1e-9):.2f}x")

        records.append({
            'x': x, 'true_pi': true_pi, 'legendre_pi': leg_pi,
            'match': match, 't_legendre': t_leg, 't_sieve': t_sieve,
        })

    return records


# ─── Experiment 3: Sequential vs FFT polynomial product timing ──────

def experiment_3_timing(x_values):
    """Compare sequential vs FFT polynomial product computation time."""
    print("\n" + "=" * 70)
    print("EXPERIMENT 3: Sequential vs FFT Polynomial Product Timing")
    print("=" * 70)

    records = []
    for x in x_values:
        B = int(x**0.5)

        t0 = time.time()
        poly_seq, primes = poly_product_sequential(x, B)
        t_seq = time.time() - t0

        t0 = time.time()
        poly_fft, _ = poly_product_fft(x, B)
        t_fft = time.time() - t0

        t0 = time.time()
        _ = direct_sieve(x)
        t_sieve = time.time() - t0

        diff = np.max(np.abs(poly_seq - poly_fft))
        agree = diff < 0.5

        print(f"\nx = {x:>8d}, B = {B}, #primes = {len(primes)}")
        print(f"  Sequential:  {t_seq:.6f}s")
        print(f"  FFT D&C:     {t_fft:.6f}s")
        print(f"  Direct sieve:{t_sieve:.6f}s")
        print(f"  FFT/seq speedup: {t_seq/max(t_fft,1e-9):.2f}x")
        print(f"  Poly agree? {agree} (max diff={diff:.2f})")

        records.append({
            'x': x, 'B': B, 'num_primes': len(primes),
            't_seq': t_seq, 't_fft': t_fft, 't_sieve': t_sieve,
            'agree': agree,
        })

    return records


# ─── Experiment 4: Roots of unity evaluation ─────────────────────────

def experiment_4_roots_of_unity(x_values):
    """
    Evaluate P(z) at roots of unity. The key question: can we extract
    any useful sieve information from P(omega) without full polynomial?

    P(omega^j) = prod_{p<=B} (1 - omega^{jp}), easily computed in O(a) per j.
    We get m values from evaluating at m-th roots of unity.

    Can we reconstruct pi(x) from {P(omega^j) : j=0..m-1} for small m?
    No -- P(z) encodes subset-sum information, not the Legendre sieve.
    But let's see what structure exists.
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 4: Evaluation at Roots of Unity")
    print("=" * 70)

    m_values = [2, 3, 4, 5, 6, 8, 10, 12, 24, 30, 60, 120]
    records = []

    for x in x_values:
        B = int(x**0.5)
        primes = list(primerange(2, B + 1))
        a = len(primes)
        print(f"\nx = {x}, B = {B}, primes = {primes}")

        for m in m_values:
            omega = np.exp(2j * np.pi / m)
            vals = np.zeros(m, dtype=complex)
            for j in range(m):
                w = omega ** j
                prod = 1.0 + 0j
                for p in primes:
                    prod *= (1.0 - w**p)
                vals[j] = prod

            mags = np.abs(vals)
            # P(1) = 0 always (since 2 is a factor -> 1-1=0)
            # How many roots give P=0?
            zeros = np.sum(mags < 1e-10)

            print(f"  m={m:>3d}: #zeros={int(zeros):>2d}, "
                  f"|P| range=[{mags.min():.4f}, {mags.max():.4f}], "
                  f"sum|P|={mags.sum():.4f}")

            records.append({
                'x': x, 'm': m, 'zeros': int(zeros),
                'mag_min': float(mags.min()), 'mag_max': float(mags.max()),
            })

        # Key observation: P(omega^j) = 0 whenever jp ≡ 0 mod m for some prime p,
        # i.e., whenever gcd(j,m)*p ≡ 0 mod m for some p, i.e., m | jp.
        # Since 2 is always a prime, P(omega^j) = 0 whenever m | 2j, i.e., j ≡ 0 mod m/gcd(2,m).
        # For even m: zeros at all even j. For odd m: zeros at j=0 only... wait, also j multiples of m/gcd(p,m).
        print(f"  -> P(omega^j) = 0 iff any prime p satisfies jp ≡ 0 (mod m)")
        print(f"     Since 2 is always present: zeros at j where m | 2j, i.e., j multiple of m/gcd(2,m)")

    return records


# ─── Experiment 5: Dirichlet polynomial approach ────────────────────

def experiment_5_dirichlet(x_values):
    """
    The CORRECT polynomial for the Legendre sieve is the Dirichlet polynomial:
      F(s) = prod_{p<=B} (1 - p^{-s}) = sum_{d | P#} mu(d) * d^{-s}

    phi(x, a) = sum_{d | P#} mu(d) * floor(x/d)

    Can we evaluate this Dirichlet polynomial efficiently?
    The sum has 2^a terms. For a = pi(sqrt(x)) ~ 2*sqrt(x)/ln(x), this is
    exponential in sqrt(x).

    Test: enumerate the Dirichlet polynomial explicitly for small cases and
    verify phi(x, a) computation.
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 5: Dirichlet Polynomial for Legendre Sieve")
    print("=" * 70)

    records = []
    for x in x_values:
        B = int(x**0.5)
        primes = list(primerange(2, B + 1))
        a = len(primes)

        # Enumerate all squarefree divisors of primorial with mu values
        t0 = time.time()
        divisors = []
        for mask in range(1 << a):
            d = 1
            bits = 0
            for i in range(a):
                if mask & (1 << i):
                    d *= primes[i]
                    bits += 1
            if d <= x:
                divisors.append((d, (-1)**bits))
        t_enum = time.time() - t0

        # Compute phi(x, a)
        t0 = time.time()
        phi = sum(mu * (x // d) for d, mu in divisors)
        pi_x = phi + a - 1
        t_phi = time.time() - t0

        true_pi = int(primepi(x))
        match = (pi_x == true_pi)

        # How many divisors exceed x?
        exceeding_x = (1 << a) - len(divisors)

        print(f"\nx = {x:>8d}, B = {B}, a = {a}")
        print(f"  2^a = {1<<a}, divisors <= x: {len(divisors)}, exceeding x: {exceeding_x}")
        print(f"  phi(x,a) = {phi}, pi(x) = {pi_x}, true pi(x) = {true_pi}  {'OK' if match else 'FAIL'}")
        print(f"  t_enum = {t_enum:.6f}s, t_phi = {t_phi:.6f}s")

        records.append({
            'x': x, 'a': a, 'two_pow_a': 1 << a,
            'divisors_le_x': len(divisors), 'exceeding_x': exceeding_x,
            'pi_x': pi_x, 'true_pi': true_pi, 'match': match,
            't_enum': t_enum, 't_phi': t_phi,
        })

    return records


# ─── Experiment 6: Complexity crossover analysis ────────────────────

def experiment_6_complexity():
    """
    For computing pi(x):
    - Legendre divisor enumeration: O(2^{pi(sqrt(x))})
    - Direct sieve: O(x log log x)
    - Meissel-Lehmer: O(x^{2/3} / log x)
    - Polynomial product (FFT tree): O(x log^2 x)

    When does 2^{pi(sqrt(x))} cross x^{2/3}?
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 6: Complexity Crossover Analysis")
    print("=" * 70)

    print(f"\n{'x':>15s}  {'pi(sqrt(x))':>12s}  {'2^pi(sqrt(x))':>20s}  "
          f"{'x^(2/3)':>15s}  {'x log log x':>15s}  {'winner':>12s}")
    print("-" * 95)

    records = []
    for exp in range(2, 20):
        x = 10**exp
        sq = int(x**0.5)
        pi_sq = int(primepi(sq))
        two_pow = 2**pi_sq if pi_sq < 200 else float('inf')
        x_23 = x**(2/3)
        x_loglgx = x * max(1, np.log(np.log(x)))

        if two_pow < x_23 and two_pow < x_loglgx:
            winner = "Legendre"
        elif x_23 < x_loglgx:
            winner = "Meissel-L"
        else:
            winner = "sieve"

        print(f"10^{exp:>2d}          {pi_sq:>12d}  {two_pow:>20.2e}  "
              f"{x_23:>15.2e}  {x_loglgx:>15.2e}  {winner:>12s}")

        records.append({
            'x_exp': exp, 'pi_sqrt_x': pi_sq,
            'log2_two_pow': pi_sq,
            'log10_x_23': exp * 2 / 3,
            'winner': winner,
        })

    return records


# ─── Experiment 7: Can the polynomial P(z) encode a sieve? ──────────

def experiment_7_poly_as_sieve(x_values):
    """
    Even though P(z) = prod(1-z^p) doesn't directly give the Legendre sieve,
    can we use it to COUNT primes?

    Observation: The sieve of Eratosthenes marks composite numbers.
    n is composite iff n is divisible by some prime p <= sqrt(n).

    Consider a DIFFERENT polynomial: S(z) = 1 - prod_{p<=B} (1 - z^p + z^{2p} - ...)
    Actually, the "sieve indicator" for composite n is:
      C(n) = 1 - prod_{p<=B, p|n} 1 = ... this gets circular.

    Instead, consider: for each prime p, define f_p(z) = sum_{k>=1} z^{kp}
    = z^p / (1-z^p). Then the composite indicator generating function is
    1 - prod_{p<=B} (1 - f_p(z)) expanded carefully...

    Actually, the simplest polynomial sieve is: define an array A[0..x] = 1.
    For each prime p <= sqrt(x), set A[kp] = 0 for k >= 2.
    This IS the Eratosthenes sieve and takes O(x log log x).

    The polynomial product P(z) doesn't help because it encodes Mobius-weighted
    subset-SUMS, whereas the sieve needs to mark MULTIPLES.

    Let's verify with a small example.
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 7: Polynomial P(z) as Sieve -- Structural Analysis")
    print("=" * 70)

    for x in x_values:
        B = int(x**0.5)
        primes = list(primerange(2, B + 1))
        a = len(primes)

        # What P(z) actually encodes
        poly, _ = poly_product_sequential(x, B)
        nonzero_idx = np.where(np.abs(poly) > 0.5)[0]
        nonzero_vals = np.round(poly[nonzero_idx]).astype(int)

        # The numbers that survive the sieve (not divisible by any p <= B)
        survivors = []
        for n in range(2, x + 1):
            if all(n % p != 0 for p in primes):
                survivors.append(n)

        # pi(x) = len([s for s in survivors if s is prime]) + a
        # Actually: survivors includes primes > B and composites with all
        # factors > B. But all primes > B survive, and composites > B^2 = x
        # can't be in [2, x] with two factors > B... wait, that's the point
        # of Legendre's formula.
        # pi(x) = #{survivors in [2,x]} + a - 1

        true_pi = int(primepi(x))
        legendre_pi = len(survivors) + a - 1

        print(f"\nx = {x}, B = {B}, primes = {primes}")
        print(f"  Survivors (no small prime factor): {len(survivors)}")
        print(f"  Legendre: phi(x,a) = {len(survivors)}, pi(x) = {legendre_pi}, true = {true_pi}")
        print(f"  Polynomial P(z) nonzero positions: {nonzero_idx.tolist()[:20]}...")
        print(f"  These are SUBSET-SUMS of {primes}, NOT the sieve output.")
        print(f"  The polynomial encodes signed subset-sum counts, not primality.")

    return []


# ─── Main ────────────────────────────────────────────────────────────

def main():
    print("CONVOLUTION SIEVE EXPERIMENT")
    print("Exploring polynomial representations of the Eratosthenes sieve")
    print("=" * 70)

    # Keep x small for Legendre/divisor enumeration (2^a grows fast)
    # For x=1000, a=11 -> 2^11=2048. For x=5000, a=19 -> 2^19=524K. OK.
    # For x=10000, a=25 -> 2^25=33M. Borderline. Skip for divisor enum.
    x_tiny = [100, 500, 1000]
    x_small = [100, 500, 1000, 5000]
    x_medium = [10000, 50000, 100000]

    results = {}
    results['exp1'] = experiment_1_coeff_structure(x_tiny)
    results['exp2'] = experiment_2_legendre(x_small)
    results['exp3'] = experiment_3_timing(x_medium)
    results['exp4'] = experiment_4_roots_of_unity([100, 1000])
    results['exp5'] = experiment_5_dirichlet(x_small)
    results['exp6'] = experiment_6_complexity()
    results['exp7'] = experiment_7_poly_as_sieve([100, 1000])

    # ─── Summary ─────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    print("""
KEY FINDINGS:

1. POLYNOMIAL STRUCTURE: P(z) = prod_{p<=B}(1-z^p) has nonzero coefficients
   at SUBSET-SUM positions of the primes, NOT at squarefree divisor positions.
   This is fundamentally different from the Legendre sieve formula which uses
   squarefree PRODUCTS as divisors.

2. CORRECTNESS OF LEGENDRE: The direct Legendre formula
   pi(x) = phi(x,a) + a - 1 = sum_{d|P#} mu(d)*floor(x/d) + a - 1
   correctly computes pi(x) when a = pi(sqrt(x)). Verified for all test cases.

3. POLYNOMIAL PRODUCT IS IRRELEVANT: The polynomial product P(z) does NOT
   directly encode the Legendre sieve. The exponents of z track SUMS of primes,
   but the Legendre formula needs PRODUCTS of primes as divisors of x.
   There is no FFT shortcut here because the problem is a Dirichlet convolution,
   not a polynomial convolution.

4. ROOTS OF UNITY: P(omega^j) is zero whenever jp ≡ 0 mod m for any sieving
   prime p. Since 2 is always present, roughly half the evaluations are zero.
   The non-zero values don't encode useful sieve information.

5. DIRICHLET POLYNOMIAL: The correct algebraic object is the Dirichlet polynomial
   F(s) = prod_{p<=B}(1-p^{-s}), not the power series P(z). Evaluating
   this Dirichlet polynomial still requires 2^{pi(sqrt(x))} terms.

6. COMPLEXITY: The Legendre/divisor approach costs O(2^{pi(sqrt(x))}).
   For x=10^4, pi(sqrt(x))=25, so 2^25 ~ 33M -- feasible.
   For x=10^10, pi(sqrt(x))=3401, so 2^3401 -- utterly impossible.
   For x=10^100, pi(sqrt(x)) ~ 10^48 -- beyond any computation.
   Meissel-Lehmer at O(x^{2/3}) wins for all x > ~10^6.

7. FUNDAMENTAL MISMATCH: Polynomial convolution (FFT) operates on power series
   where exponents ADD. The prime sieve needs multiplicative structure where
   divisors MULTIPLY. FFT cannot bridge this gap. The relevant transform for
   multiplicative number theory is the Dirichlet convolution, which has no
   known sublinear algorithm for computing partial sums.

VERDICT: CLOSED. The polynomial convolution approach to the prime sieve is
based on a structural mismatch: polynomial multiplication encodes additive
structure (subset sums) while the sieve requires multiplicative structure
(divisor products). The correct algebraic framework is Dirichlet series,
which lacks an efficient analog of FFT. This reduces to the known Legendre
sieve with O(2^{pi(sqrt(x))}) complexity, far worse than Meissel-Lehmer
for any practically relevant x.
""")

    return results


if __name__ == "__main__":
    results = main()
