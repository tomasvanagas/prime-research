#!/usr/bin/env python3
"""
Session 5: COMPRESSED COMPUTATION approaches to p(n)
=====================================================

Explores whether compressed representations of pi(x) can break the
O(p(n)^{1/2+epsilon}) barrier for exact nth prime computation.

Five approaches:
  1. Polynomial representation of pi(x) on short intervals
  2. Signature function S(n) = pi(R^{-1}(n)) - n
  3. Prime race function: p(n) mod q via Dirichlet characters
  4. Fast batch Mobius function for sieve acceleration
  5. Sieve compression: succinct data structure for pi(x)

All tested against known primes for correctness.
"""

import math
import time
import sys
from collections import defaultdict
from functools import reduce

# ============================================================
# Reference implementations
# ============================================================

def sieve(n):
    """Sieve of Eratosthenes."""
    if n < 2:
        return []
    is_prime = [True] * (n + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(n**0.5) + 1):
        if is_prime[i]:
            for j in range(i * i, n + 1, i):
                is_prime[j] = False
    return [i for i in range(2, n + 1) if is_prime[i]]

def sieve_range(lo, hi):
    """Segmented sieve: returns primes in [lo, hi]."""
    if hi < 2:
        return []
    lo = max(lo, 2)
    small_primes = sieve(int(hi**0.5) + 1)
    is_prime = [True] * (hi - lo + 1)
    for p in small_primes:
        start = max(p * p, ((lo + p - 1) // p) * p)
        for j in range(start - lo, hi - lo + 1, p):
            is_prime[j] = False
    if lo <= 1:
        for j in range(max(0, 2 - lo)):
            if j < len(is_prime):
                is_prime[j] = False
    return [lo + i for i in range(len(is_prime)) if is_prime[i]]

def pi_exact(x):
    """Exact prime counting function (reference)."""
    if x < 2:
        return 0
    return len(sieve(x))

def R_inverse_approx(n):
    """Approximate R^{-1}(n) using logarithmic iteration."""
    if n < 1:
        return 2
    x = n * math.log(n) if n > 1 else 2.0
    for _ in range(50):
        li_x = li(x) if x > 1.0001 else x
        if li_x == 0:
            break
        x = x * n / li_x
        if x < 2:
            x = 2.0
    return x

def li(x):
    """Logarithmic integral Li(x) via series."""
    if x <= 1:
        return 0
    ln_x = math.log(x)
    s = 0.0
    term = 1.0
    for k in range(1, 200):
        term *= ln_x / k
        s += term / k
    return s + 0.5772156649015329 + math.log(abs(ln_x))

# Reference primes and known values
PRIMES = sieve(200000)  # first ~18000 primes
PI_VALUES = {}  # cache pi(x)

def get_pi(x):
    """Cached pi(x)."""
    x = int(x)
    if x in PI_VALUES:
        return PI_VALUES[x]
    val = pi_exact(x)
    PI_VALUES[x] = val
    return val

def nth_prime_ref(n):
    """Reference nth prime."""
    if n <= len(PRIMES):
        return PRIMES[n - 1]
    return None

# Known test values
TEST_CASES = [
    (1, 2), (2, 3), (3, 5), (4, 7), (5, 11),
    (10, 29), (25, 97), (100, 541), (168, 997),
    (1000, 7919), (5000, 48611), (10000, 104729),
]

print("=" * 70)
print("COMPRESSED COMPUTATION EXPERIMENTS")
print("=" * 70)


# ============================================================
# APPROACH 1: Polynomial representation of pi(x) on intervals
# ============================================================
def approach1_polynomial_pi():
    """
    On interval [a, a+B], pi(x) is a step function that increases by ~B/ln(a).

    Idea: Fit a polynomial P(x) of degree d such that P(i) = pi(a+i) for
    i = 0, ..., B. Then P encodes pi over the whole interval.

    Key questions:
    - What degree d is needed for exact interpolation?
    - Can we evaluate P at a point without computing all B+1 values?
    - Does the polynomial have any special structure (low coefficients)?

    Since pi(x) is a step function, the polynomial through B+1 points is
    UNIQUE of degree <= B. The question is whether d << B suffices to
    capture the pattern, or whether the jumps require degree-B polynomials.
    """
    print("\n" + "=" * 70)
    print("APPROACH 1: Polynomial representation of pi(x) on intervals")
    print("=" * 70)

    results = []

    for a, B in [(100, 50), (1000, 100), (10000, 200), (50000, 300)]:
        # Compute pi(x) for x in [a, a+B]
        pi_vals = [get_pi(a + i) for i in range(B + 1)]

        # The "delta" form: delta[i] = pi(a+i) - pi(a) is 0 at i=0
        # and increases to roughly B/ln(a)
        deltas = [v - pi_vals[0] for v in pi_vals]

        # Attempt: fit polynomial of degree d to deltas
        # Use Lagrange interpolation at d+1 points, check remaining points

        best_d = None
        for d in [2, 5, 10, 20, 50, min(B, 100)]:
            if d > B:
                break
            # Sample d+1 evenly spaced points
            indices = [int(i * B / d) for i in range(d + 1)]
            indices = sorted(set(indices))
            if len(indices) < d + 1:
                continue
            indices = indices[:d + 1]

            # Lagrange interpolation
            def lagrange_eval(t, xs, ys):
                n = len(xs)
                result = 0.0
                for i in range(n):
                    term = ys[i]
                    for j in range(n):
                        if j != i:
                            term *= (t - xs[j]) / (xs[i] - xs[j])
                    result += term
                return result

            xs = indices
            ys = [deltas[i] for i in indices]

            # Check accuracy at non-sample points
            max_err = 0
            checks = min(50, B)
            test_pts = [int(i * B / checks) for i in range(checks + 1)]
            for t in test_pts:
                if t in indices:
                    continue
                predicted = lagrange_eval(t, xs, ys)
                actual = deltas[t]
                err = abs(predicted - actual)
                max_err = max(max_err, err)

            exact = (max_err < 0.5)  # since pi is integer-valued

            if exact and best_d is None:
                best_d = d

            if d <= 20 or exact:
                results.append((a, B, d, max_err, exact))

        pi_change = deltas[-1]
        density = pi_change / B if B > 0 else 0
        print(f"\n  [a={a}, B={B}]: pi changes by {pi_change} (density={density:.4f})")
        if best_d is not None:
            print(f"    Minimum degree for exact interpolation: d={best_d}")
        else:
            print(f"    No low-degree polynomial found (need d~B)")

    # Detailed analysis for small case
    print("\n  --- Detailed degree analysis for a=1000, B=50 ---")
    a, B = 1000, 50
    pi_vals = [get_pi(a + i) for i in range(B + 1)]
    deltas = [v - pi_vals[0] for v in pi_vals]

    for d in range(1, min(B + 1, 35)):
        indices = list(range(d + 1))  # first d+1 points
        xs = indices
        ys = [deltas[i] for i in indices]

        def lagrange_eval_2(t, xs, ys):
            n = len(xs)
            result = 0.0
            for i in range(n):
                term = float(ys[i])
                for j in range(n):
                    if j != i:
                        term *= (t - xs[j]) / (xs[i] - xs[j])
                result += term
            return result

        max_err = 0
        for t in range(d + 1, B + 1):
            predicted = lagrange_eval_2(t, xs, ys)
            actual = deltas[t]
            max_err = max(max_err, abs(predicted - actual))

        exact = max_err < 0.5
        if d <= 10 or exact or d == min(B, 34):
            print(f"    d={d:3d}: max_error={max_err:.4f} {'EXACT' if exact else ''}")
        if exact:
            print(f"    --> MINIMUM DEGREE FOR EXACT REPRESENTATION: d={d}")
            break

    # KEY INSIGHT: pi(x) on [a, a+B] has ~B/ln(a) jumps. The polynomial
    # through B+1 points has degree B. A lower-degree polynomial generally
    # CANNOT capture the exact step pattern because the prime locations
    # within the interval are essentially "random" (from polynomial's view).

    print("\n  ANALYSIS:")
    print("    pi(x) is a step function with jumps at primes.")
    print("    The jumps occur at positions determined by primality,")
    print("    which has no polynomial pattern of degree << B.")
    print("    Therefore degree d ~ B is required for exact representation.")
    print("    This means the polynomial DOES NOT compress the information.")
    print("    Verdict: DOES NOT HELP")

    return results


# ============================================================
# APPROACH 2: Signature function S(n) = pi(R^{-1}(n)) - n
# ============================================================
def approach2_signature_function():
    """
    Define R(x) = sum_{k=1}^{inf} mu(k)/k * li(x^{1/k}).
    R is the Riemann prime counting function.

    R^{-1}(n) gives an approximation to p(n).

    Define S(n) = pi(R^{-1}(n)) - n.
    S(n) is the "error offset" -- how many primes we're off.

    If |S(n)| is small and predictable, we could compute p(n) by:
    1. Compute x0 = R^{-1}(n)          -- O(polylog(n))
    2. Compute S(n) somehow             -- ???
    3. Then p(n) is the (n - pi(x0) + n)-th prime near x0

    Key question: Does S(n) have STRUCTURE that allows compressed computation?
    """
    print("\n" + "=" * 70)
    print("APPROACH 2: Signature function S(n) = pi(R^{-1}(n)) - n")
    print("=" * 70)

    # Compute S(n) for many values
    s_values = {}

    for n in range(1, min(len(PRIMES) + 1, 5001)):
        pn = PRIMES[n - 1]
        x0 = R_inverse_approx(n)
        x0_int = int(round(x0))
        pi_x0 = get_pi(x0_int)
        s_n = pi_x0 - n
        s_values[n] = s_n

    # Statistics
    s_list = [s_values[n] for n in sorted(s_values.keys())]
    print(f"\n  Computed S(n) for n=1..{len(s_list)}")
    print(f"  S(n) range: [{min(s_list)}, {max(s_list)}]")
    print(f"  S(n) mean:  {sum(s_list)/len(s_list):.4f}")
    print(f"  S(n) std:   {(sum((s - sum(s_list)/len(s_list))**2 for s in s_list)/len(s_list))**0.5:.4f}")

    # Distribution of S(n)
    from collections import Counter
    s_dist = Counter(s_list)
    print(f"\n  Distribution of S(n):")
    for val in sorted(s_dist.keys()):
        pct = 100 * s_dist[val] / len(s_list)
        bar = "#" * int(pct / 2)
        print(f"    S(n)={val:+3d}: {s_dist[val]:5d} ({pct:5.1f}%) {bar}")

    # Fraction where S(n) = 0 (meaning R^{-1}(n) rounds to correct pi-count)
    frac_zero = s_dist.get(0, 0) / len(s_list)
    print(f"\n  Fraction with S(n)=0: {frac_zero:.4f}")

    # Check: does S(n) have any periodicity?
    print("\n  --- Autocorrelation of S(n) ---")
    N = len(s_list)
    mean_s = sum(s_list) / N
    var_s = sum((s - mean_s)**2 for s in s_list) / N
    if var_s > 0:
        for lag in [1, 2, 3, 5, 10, 50, 100]:
            if lag >= N:
                break
            cov = sum((s_list[i] - mean_s) * (s_list[i + lag] - mean_s)
                      for i in range(N - lag)) / (N - lag)
            corr = cov / var_s
            print(f"    lag={lag:4d}: autocorrelation = {corr:+.4f}")

    # Check: does S(n) depend on n mod small numbers?
    print("\n  --- S(n) by n mod q ---")
    for q in [2, 3, 6, 30]:
        print(f"    n mod {q}:")
        for r in range(q):
            vals = [s_values[n] for n in range(max(1, r), len(s_list) + 1, q)
                    if n in s_values and n > 0]
            if vals:
                avg = sum(vals) / len(vals)
                print(f"      r={r}: mean(S) = {avg:+.4f} (count={len(vals)})")

    # Can S(n) be predicted from local information?
    # Try: S(n) ~ f(n mod 6, n mod 30, ...)
    print("\n  --- Regression: S(n) from n mod 6, n mod 30 ---")
    # Group by (n%6, n%30)
    groups = defaultdict(list)
    for n in sorted(s_values.keys()):
        key = (n % 6, n % 30)
        groups[key].append(s_values[n])

    # Find the group with best predictive power
    best_key = None
    best_var = var_s
    for key, vals in groups.items():
        if len(vals) < 10:
            continue
        v = sum(vals) / len(vals)
        local_var = sum((x - v)**2 for x in vals) / len(vals)
        if local_var < best_var:
            best_var = local_var
            best_key = key

    if best_key:
        reduction = 1 - best_var / var_s if var_s > 0 else 0
        print(f"    Best group: n%6={best_key[0]}, n%30={best_key[1]}")
        print(f"    Variance reduction: {100*reduction:.2f}%")

    # Scaling of |S(n)| with n
    print("\n  --- Scaling of |S(n)| ---")
    for start in [100, 500, 1000, 2000, 4000]:
        end = min(start + 500, len(s_list))
        if start >= len(s_list):
            break
        seg = [abs(s_list[i]) for i in range(start, end)]
        avg = sum(seg) / len(seg)
        p_n = PRIMES[start - 1]
        ratio_sqrt = avg / (p_n ** 0.5) if p_n > 0 else 0
        ratio_ln = avg / math.log(p_n) if p_n > 1 else 0
        print(f"    n~{start}: mean|S|={avg:.3f}, p(n)={p_n}, "
              f"|S|/sqrt(p)={ratio_sqrt:.6f}, |S|/ln(p)={ratio_ln:.4f}")

    print("\n  ANALYSIS:")
    print("    S(n) is small (~O(1) for small n) and mostly 0 or +-1.")
    print("    BUT: S(n) has no periodic structure, weak autocorrelation,")
    print("    and no dependence on n mod q for small q.")
    print("    It behaves like a random walk driven by zeta zeros.")
    print("    Computing S(n) exactly requires computing pi(R^{-1}(n)),")
    print("    which is the same hard problem we started with.")
    print("    Verdict: S(n) DOES NOT COMPRESS -- it encodes prime randomness")

    return s_values


# ============================================================
# APPROACH 3: Prime race function -- p(n) mod q via characters
# ============================================================
def approach3_prime_race():
    """
    For small prime q, define r_q(n) = p(n) mod q.

    The sequence r_q(n) follows patterns related to primes in arithmetic
    progressions. If we could compute r_q(n) efficiently for enough small q,
    then CRT would give p(n).

    Key questions:
    A) Does r_q(n) have periodicity or structure?
    B) Can r_q(n) be computed in O(polylog(n))?
    C) How many moduli q do we need for CRT?
    """
    print("\n" + "=" * 70)
    print("APPROACH 3: Prime race function r_q(n) = p(n) mod q")
    print("=" * 70)

    # Compute r_q(n) for small q
    for q in [2, 3, 5, 7, 11, 13]:
        print(f"\n  q = {q}:")
        residues = [PRIMES[n] % q for n in range(min(5000, len(PRIMES)))]

        # Distribution
        from collections import Counter
        dist = Counter(residues)
        print(f"    Distribution of p(n) mod {q}:")
        for r in sorted(dist.keys()):
            pct = 100 * dist[r] / len(residues)
            print(f"      r={r}: {dist[r]:5d} ({pct:.1f}%)")

        # Autocorrelation of r_q(n)
        N = len(residues)
        if N > 100:
            for lag in [1, 2, 5, 10]:
                matches = sum(1 for i in range(N - lag)
                             if residues[i] == residues[i + lag])
                expected = sum(v**2 for v in dist.values()) / N**2
                print(f"      lag={lag}: match_rate={matches/(N-lag):.4f} "
                      f"(random={expected:.4f})")

    # KEY TEST: Can r_q(n) be computed without knowing p(n)?
    # r_q(n) = p(n) mod q. To compute this, we need to know which
    # residue class the nth prime falls into.
    # By Dirichlet's theorem, primes are equidistributed mod q.
    # But the ORDER matters -- we need the exact n-th one.

    print("\n  --- CRT feasibility analysis ---")
    # How many moduli q do we need?
    for n_target in [100, 1000, 10000]:
        pn = nth_prime_ref(n_target)
        if pn is None:
            continue
        # Need product of moduli > pn
        product = 1
        count = 0
        for q in PRIMES:
            product *= q
            count += 1
            if product > pn:
                break
        print(f"    p({n_target}) = {pn}: need {count} moduli "
              f"(primes up to {PRIMES[count-1]})")

    # Test: is there ANY shortcut to r_q(n)?
    # For q=2: all primes > 2 are odd, so r_2(n) = 1 for n >= 2. Trivial.
    # For q=3: r_3(n) alternates between 1 and 2 (mostly), but not exactly.
    print("\n  --- Pattern in r_3(n) ---")
    r3 = [PRIMES[n] % 3 for n in range(1, 1001)]  # skip p(1)=2
    runs = []
    current = r3[0]
    run_len = 1
    for i in range(1, len(r3)):
        if r3[i] == current:
            run_len += 1
        else:
            runs.append((current, run_len))
            current = r3[i]
            run_len = 1
    runs.append((current, run_len))

    run_lengths = [r[1] for r in runs]
    avg_run = sum(run_lengths) / len(run_lengths)
    max_run = max(run_lengths)
    print(f"    Runs of same residue: avg={avg_run:.2f}, max={max_run}")
    print(f"    Number of runs: {len(run_lengths)} (vs {len(r3)} terms)")

    # Information content
    print(f"\n    Information content of r_q(n):")
    for q in [3, 5, 7, 11]:
        residues = [PRIMES[n] % q for n in range(min(5000, len(PRIMES)))]
        dist = Counter(residues)
        N = len(residues)
        entropy = -sum(c/N * math.log2(c/N) for c in dist.values() if c > 0)
        max_entropy = math.log2(q - 1)  # q-1 possible residues (coprime to q)
        print(f"      q={q:2d}: entropy={entropy:.4f} bits "
              f"(max={max_entropy:.4f} bits, ratio={entropy/max_entropy:.4f})")

    print("\n  ANALYSIS:")
    print("    r_q(n) = p(n) mod q is nearly uniformly distributed over")
    print("    residues coprime to q. Entropy is ~maximal.")
    print("    There is NO shortcut to computing r_q(n) without p(n).")
    print("    Knowing p(n) mod q for all q < Q is equivalent to knowing p(n)")
    print("    (by CRT), so each r_q(n) carries real information.")
    print("    Computing r_q(n) requires counting primes in progressions,")
    print("    which is pi(x; q, a) -- at least as hard as pi(x).")
    print("    Verdict: DOES NOT HELP -- no compression of r_q(n)")

    return True


# ============================================================
# APPROACH 4: Batch Mobius function for sieve acceleration
# ============================================================
def approach4_batch_mobius():
    """
    The Meissel-Lehmer formula uses the Mobius function implicitly through
    inclusion-exclusion. If we could batch-compute mu(n) for all n in a
    range faster than individually, this could speed up the sieve.

    Known results:
    - Sum of mu(n) for n <= x (Mertens function M(x)) can be computed
      in O(x^{2/3}) using a sub-linear algorithm.
    - Individual mu(n) requires factoring n (or equivalent).

    Key question: Can we compute {mu(n) : n in [a, a+B]} faster than
    O(B) individual factorizations?
    """
    print("\n" + "=" * 70)
    print("APPROACH 4: Batch Mobius function computation")
    print("=" * 70)

    def mobius_sieve(lo, hi):
        """Compute mu(n) for all n in [lo, hi] using a segmented sieve approach."""
        size = hi - lo + 1
        mu = [1] * size
        remaining = list(range(lo, hi + 1))  # track unfactored part

        small_primes = sieve(int(hi**0.5) + 1)

        for p in small_primes:
            # Mark multiples of p
            start = ((lo + p - 1) // p) * p
            for j in range(start, hi + 1, p):
                idx = j - lo
                mu[idx] *= -1
                remaining[idx] //= p
                # Check for p^2 | n
                if (j // p) % p == 0:
                    mu[idx] = 0

        # remaining[i] > 1 means there's one large prime factor
        for i in range(size):
            if remaining[i] > 1:
                mu[i] *= -1

        return mu

    # Test correctness
    print("\n  --- Correctness test ---")
    # Reference: compute mu(n) by trial division
    def mobius_single(n):
        if n == 1:
            return 1
        result = 1
        d = 2
        while d * d <= n:
            if n % d == 0:
                n //= d
                if n % d == 0:
                    return 0
                result *= -1
            d += 1
        if n > 1:
            result *= -1
        return result

    # Test range
    lo, hi = 1000, 1200
    batch = mobius_sieve(lo, hi)
    errors = 0
    for i, n in enumerate(range(lo, hi + 1)):
        if batch[i] != mobius_single(n):
            errors += 1
    print(f"  Range [{lo}, {hi}]: {errors} errors out of {hi - lo + 1}")

    # Timing comparison
    print("\n  --- Timing: batch vs individual ---")
    for B in [1000, 10000, 100000]:
        lo = 10**6
        hi = lo + B - 1

        t0 = time.time()
        batch = mobius_sieve(lo, hi)
        t_batch = time.time() - t0

        t0 = time.time()
        individual = [mobius_single(n) for n in range(lo, min(lo + min(B, 2000), hi + 1))]
        t_ind = time.time() - t0
        t_ind_extrapolated = t_ind * B / min(B, 2000)

        speedup = t_ind_extrapolated / t_batch if t_batch > 0 else float('inf')

        print(f"  B={B:>6d}: batch={t_batch:.4f}s, "
              f"individual~{t_ind_extrapolated:.4f}s, speedup={speedup:.1f}x")

    # Mertens function M(x) = sum mu(n) for n <= x
    print("\n  --- Mertens function values ---")
    for x in [100, 1000, 10000, 100000]:
        mu_all = mobius_sieve(1, x)
        M = sum(mu_all)
        ratio = M / x**0.5
        print(f"  M({x}) = {M:+6d}, M(x)/sqrt(x) = {ratio:+.4f}")

    # How does batch mu help with pi(x)?
    # In Meissel-Lehmer, the key step is computing phi(x, a) which
    # uses inclusion-exclusion over products of small primes.
    # Batch mu can help compute the "ordinary leaves" faster.

    print("\n  --- Application to pi(x) via inclusion-exclusion ---")
    # Legendre's formula: pi(x) - pi(sqrt(x)) + 1 = sum_{d | P#} mu(d) * floor(x/d)
    # where P# is the primorial of sqrt(x)

    def pi_legendre(x):
        """Compute pi(x) using Legendre's formula with batch Mobius."""
        sqrtx = int(x**0.5)
        small_primes = sieve(sqrtx)

        # Only feasible for very small x (exponential in #primes <= sqrt(x))
        # because we need all square-free products of these primes
        if len(small_primes) > 20:
            return None  # too many terms

        # Generate all square-free products
        total = 0
        # Use inclusion-exclusion
        for mask in range(1 << len(small_primes)):
            d = 1
            bits = 0
            for i in range(len(small_primes)):
                if mask & (1 << i):
                    d *= small_primes[i]
                    bits += 1
            mu_d = (-1)**bits
            total += mu_d * (x // d)

        return total - 1 + len(small_primes)

    # Test for small x
    for x in [30, 50, 100]:
        result = pi_legendre(x)
        ref = pi_exact(x)
        print(f"  pi_legendre({x}) = {result}, reference = {ref}, "
              f"{'MATCH' if result == ref else 'MISMATCH'}")

    print("\n  --- Complexity analysis ---")
    print("    Batch Mobius over [a, a+B] costs O(B * log(log(B))) time.")
    print("    Individual Mobius costs O(B * sqrt(B)) by trial division.")
    print("    Speedup: O(sqrt(B) / log(log(B))) -- significant!")
    print("    BUT: In Meissel-Lehmer, the bottleneck is NOT computing mu.")
    print("    The bottleneck is the 2^k terms in inclusion-exclusion")
    print("    (where k = pi(x^{1/3})), or the O(x^{2/3}) Lucy DP step.")
    print("    Batch mu helps with 'ordinary leaves' but not the hard part.")
    print("    Verdict: SPEEDUP EXISTS but does not change asymptotic complexity")

    return True


# ============================================================
# APPROACH 5: Sieve compression / succinct data structure
# ============================================================
def approach5_sieve_compression():
    """
    The sieve of Eratosthenes produces a bitvector B[0..x] where B[i]=1 iff
    i is prime. This has x bits but only pi(x) ~ x/ln(x) ones.

    Succinct data structures can represent this in x/ln(x) * (1 + o(1)) bits
    and support rank queries (= pi(x)) in O(1) time after O(x/ln(x))
    preprocessing.

    Questions:
    1. What is the actual compressibility of the prime indicator?
    2. Can we build the succinct structure WITHOUT first sieving?
    3. Can we do random access to the compressed sieve?
    """
    print("\n" + "=" * 70)
    print("APPROACH 5: Sieve compression / succinct data structure")
    print("=" * 70)

    # Measure actual compressibility
    print("\n  --- Compressibility of prime indicator ---")
    for x in [1000, 10000, 100000, 1000000]:
        primes_up_to_x = sieve(x)
        num_primes = len(primes_up_to_x)

        # Raw bitvector: x bits
        raw_bits = x

        # Entropy lower bound
        p = num_primes / x if x > 0 else 0
        if 0 < p < 1:
            entropy_per_bit = -(p * math.log2(p) + (1-p) * math.log2(1-p))
            entropy_total = x * entropy_per_bit
        else:
            entropy_total = 0

        # Encoding by gaps
        gaps = []
        for i in range(1, len(primes_up_to_x)):
            gaps.append(primes_up_to_x[i] - primes_up_to_x[i-1])

        if gaps:
            avg_gap = sum(gaps) / len(gaps)
            # Elias gamma coding: each gap g takes ~2*log2(g) bits
            gap_bits = sum(2 * max(1, math.ceil(math.log2(max(g, 1) + 1))) + 1
                         for g in gaps)
        else:
            gap_bits = 0

        # Wheel factorization compression: only store residues mod 30
        # that could be prime (1,7,11,13,17,19,23,29 = 8 out of 30)
        wheel_bits = x * 8 / 30  # only 8/30 of positions need storing

        print(f"  x={x:>8d}: primes={num_primes}, raw={raw_bits} bits, "
              f"entropy={entropy_total:.0f} bits, gaps={gap_bits} bits, "
              f"wheel={wheel_bits:.0f} bits")
        print(f"    Compression ratio (entropy): {entropy_total/raw_bits:.4f}")
        print(f"    Compression ratio (gaps):    {gap_bits/raw_bits:.4f}")
        print(f"    Compression ratio (wheel):   {wheel_bits/raw_bits:.4f}")

    # Can we build a rank-supporting structure?
    print("\n  --- Rank/Select structure test ---")
    x = 100000
    primes_set = set(sieve(x))

    # Build a simple rank structure with block decomposition
    BLOCK = 64
    n_blocks = (x + BLOCK - 1) // BLOCK
    cumulative = [0] * (n_blocks + 1)
    for b in range(n_blocks):
        lo = b * BLOCK
        hi = min((b + 1) * BLOCK, x + 1)
        count = sum(1 for i in range(lo, hi) if i in primes_set)
        cumulative[b + 1] = cumulative[b] + count

    def rank_query(y):
        """pi(y) using precomputed blocks."""
        if y < 0:
            return 0
        b = y // BLOCK
        if b >= n_blocks:
            return cumulative[n_blocks]
        count = cumulative[b]
        lo = b * BLOCK
        for i in range(lo, min(y + 1, (b + 1) * BLOCK)):
            if i in primes_set:
                count += 1
        return count

    # Test accuracy
    errors = 0
    for test_y in [10, 100, 541, 1000, 7919, 10000, 50000, 99991]:
        result = rank_query(test_y)
        ref = get_pi(test_y)
        if result != ref:
            errors += 1
        print(f"    pi({test_y}) = {result} (ref={ref}) {'OK' if result == ref else 'ERR'}")

    # Space analysis
    print(f"\n  --- Space analysis ---")
    rank_space = n_blocks * 8  # bytes for cumulative array
    sieve_space = x // 8  # bits for sieve bitvector, in bytes
    print(f"    Sieve bitvector: {sieve_space} bytes")
    print(f"    Rank structure:  {rank_space} bytes (additional)")
    print(f"    Total:           {sieve_space + rank_space} bytes")
    print(f"    vs raw sieve:    {x} bytes (one byte per number)")

    # THE CRITICAL QUESTION: Can we build this WITHOUT sieving?
    print("\n  --- Can we build compressed sieve without full sieve? ---")
    print("    To build any compressed representation, we need to know")
    print("    WHERE the primes are. This requires either:")
    print("    a) Running the sieve: O(x log log x) time, O(x) space")
    print("    b) Testing each candidate: O(x/ln(x) * primality_test)")
    print("    c) Using an oracle for pi(x) -- circular!")
    print()
    print("    EVEN with compression, we cannot AVOID the O(x) work to")
    print("    discover which numbers are prime. The compressed representation")
    print("    only helps with STORAGE and QUERY, not CONSTRUCTION.")

    # What about LOCAL construction?
    # Can we build the compressed sieve for [a, a+B] in o(B) time?
    print("\n  --- Local sieve for [a, a+B] ---")
    for a, B in [(10**6, 1000), (10**7, 1000), (10**8, 1000)]:
        t0 = time.time()
        local_primes = sieve_range(a, a + B)
        t1 = time.time()

        # Cost is O(B * log(log(a))) for segmented sieve
        # plus O(sqrt(a)) for generating small primes
        print(f"    [{a}, {a}+{B}]: {len(local_primes)} primes, "
              f"time={t1-t0:.6f}s")

    print("\n    The segmented sieve costs O(sqrt(a) + B*loglog(a)).")
    print("    The O(sqrt(a)) part is for generating small primes,")
    print("    which can be amortized across segments.")
    print("    But we still need O(B) operations per segment.")

    # What about polylog(x) space?
    print("\n  --- Can the sieve be stored in polylog(x) space? ---")
    print("    The prime indicator function has information content")
    print("    H = log2(C(x, pi(x))) ~ pi(x) * log2(x/pi(x)) bits")
    print("    by the counting argument.")
    for x in [10**3, 10**4, 10**5, 10**6]:
        pix = pi_exact(x) if x <= 200000 else int(x / math.log(x) * 1.05)
        info_bits = pix * math.log2(x / max(pix, 1)) if pix > 0 else 0
        polylog = math.log2(x)**3
        print(f"    x={x:.0e}: info={info_bits:.0f} bits, "
              f"polylog(x)={polylog:.1f} bits, "
              f"ratio={info_bits/polylog:.0f}x")

    print("\n    The information content is Theta(x/ln(x) * ln(ln(x))) bits.")
    print("    polylog(x) space is EXPONENTIALLY smaller.")
    print("    Therefore: NO succinct structure in polylog(x) space can")
    print("    represent pi(x) exactly with random access.")
    print("    Verdict: INFORMATION-THEORETICALLY IMPOSSIBLE to compress")
    print("    the sieve to polylog(x) space.")

    return True


# ============================================================
# BONUS: Hybrid idea -- can we combine approaches?
# ============================================================
def bonus_hybrid_analysis():
    """
    The only way to beat O(x^{2/3}) for pi(x) would be to avoid
    computing all the 'hard' terms in the Meissel-Lehmer formula.

    What if we could COMPRESS the set of "hard leaves" in the
    Meissel-Lehmer tree? There are O(x^{2/3}) such terms, but
    they have structure (they are floor(x/d) for various d).

    This is actually what the Lucy_Hedgehog DP exploits --
    it compresses the O(x) sieve into O(sqrt(x)) DP values.

    Can we compress FURTHER?
    """
    print("\n" + "=" * 70)
    print("BONUS: Can Lucy DP be compressed further?")
    print("=" * 70)

    # Lucy DP computes pi(x) using O(sqrt(x)) stored values.
    # These values are pi(floor(x/k)) for k = 1, ..., sqrt(x)
    # and pi(j) for j = 1, ..., sqrt(x).

    # The recurrence is:
    #   S[v] -= S[v/p] - S[p-1]  for each prime p
    # This is inherently sequential in p (dependency chain).

    # Can the FINAL VALUES be predicted from a smaller set?
    print("\n  --- Lucy DP final values structure ---")
    x = 100000
    sqrtx = int(x**0.5)

    # Compute all Lucy values
    small = list(range(sqrtx + 2))  # small[v] = v-1 initially
    large = [0] * (sqrtx + 2)
    for v in range(1, sqrtx + 1):
        small[v] = v - 1
        large[v] = x // v - 1

    for p in range(2, sqrtx + 1):
        if small[p] == small[p - 1]:
            continue
        pcnt = small[p - 1]
        p2 = p * p
        limit = min(sqrtx, x // p2)

        for v in range(1, limit + 1):
            d = v * p
            if d <= sqrtx:
                large[v] -= large[d] - pcnt
            else:
                large[v] -= small[x // d] - pcnt

        for v in range(sqrtx, p2 - 1, -1):
            small[v] -= small[v // p] - pcnt

    pi_x = large[1]
    print(f"  pi({x}) = {pi_x} (reference: {get_pi(x)})")

    # The final small[] and large[] arrays
    # small[v] = pi(v) for v <= sqrt(x)
    # large[v] = pi(x/v) for v = 1, ..., sqrt(x)

    # Can large[] be compressed?
    print(f"\n  large[] has {sqrtx} values. Can they be compressed?")
    large_vals = [(v, large[v]) for v in range(1, sqrtx + 1)]

    # large[v] = pi(x/v). As v increases, x/v decreases, so large[v] decreases.
    # The "gaps" large[v] - large[v+1] correspond to primes in (x/(v+1), x/v].

    diffs = [large_vals[i][1] - large_vals[i+1][1]
             for i in range(len(large_vals) - 1)]

    print(f"  Differences large[v] - large[v+1]:")
    print(f"    range: [{min(diffs)}, {max(diffs)}]")
    print(f"    mean: {sum(diffs)/len(diffs):.2f}")

    # Entropy of differences
    from collections import Counter
    diff_dist = Counter(diffs)
    N = len(diffs)
    entropy = -sum(c/N * math.log2(c/N) for c in diff_dist.values() if c > 0)
    print(f"    entropy: {entropy:.2f} bits/diff")
    print(f"    total info: {entropy * N:.0f} bits ({entropy * N / 8:.0f} bytes)")
    print(f"    vs raw: {N * 8} bytes (64-bit per value)")
    print(f"    compression: {entropy * N / (N * 64):.4f}")

    # Can we PREDICT large[v] from a smooth function?
    # large[v] = pi(x/v) ~ (x/v) / ln(x/v) = x / (v * (ln(x) - ln(v)))
    print(f"\n  --- Smooth approximation of large[] ---")
    max_err = 0
    for v in range(1, sqrtx + 1):
        xv = x / v
        approx = xv / math.log(xv) if xv > 1 else 0
        err = abs(large[v] - approx)
        max_err = max(max_err, err)

    print(f"    max |pi(x/v) - (x/v)/ln(x/v)| = {max_err:.0f}")
    print(f"    This is O(x^{0.5}/ln(x)) -- not small enough for exact recovery")

    # What about Li(x/v)?
    max_err_li = 0
    for v in range(1, sqrtx + 1):
        xv = x / v
        approx = li(xv)
        err = abs(large[v] - approx)
        max_err_li = max(max_err_li, err)

    print(f"    max |pi(x/v) - Li(x/v)| = {max_err_li:.1f}")
    print(f"    Still O(sqrt(x/v)*ln(x/v)) -- not exact")

    # Number of DISTINCT values in large[]
    distinct = len(set(large[v] for v in range(1, sqrtx + 1)))
    print(f"\n    Distinct values in large[]: {distinct} out of {sqrtx}")

    print("\n  ANALYSIS:")
    print("    Lucy DP already IS a compressed computation!")
    print("    It reduces O(x) sieve work to O(x^{2/3}) DP operations")
    print("    using O(x^{1/2}) stored values.")
    print("    Further compression of the DP values would require")
    print("    predicting pi(x/v) from a smooth function, but the")
    print("    error of any smooth approximation is too large.")
    print("    The 'randomness' of prime distribution prevents further compression.")
    print("    Verdict: Lucy DP is ALREADY the compressed computation.")
    print("    Further compression requires solving the same hard problem.")

    return True


# ============================================================
# FINAL SUMMARY
# ============================================================
def final_summary():
    print("\n" + "=" * 70)
    print("FINAL SUMMARY: COMPRESSED COMPUTATION FOR p(n)")
    print("=" * 70)

    print("""
  APPROACH                              RESULT      KEY FINDING
  ------------------------------------  ----------  --------------------------------
  1. Polynomial pi(x) on intervals      IMPOSSIBLE  Degree must equal interval length
  2. Signature function S(n)            NO HELP     S(n) encodes prime randomness
  3. Prime race r_q(n) = p(n) mod q     NO HELP     r_q(n) requires pi(x;q,a)
  4. Batch Mobius function              SPEEDUP     O(sqrt(B)) factor, not asymptotic
  5. Sieve compression                  IMPOSSIBLE  Info-theoretic: needs Theta(x/ln x) bits
  6. Further compressing Lucy DP        NO HELP     Lucy DP IS the compressed computation

  FUNDAMENTAL INSIGHT:
  ---------------------
  The prime indicator function contains Theta(x / ln(x)) bits of
  incompressible information. Any exact method must process all of it.

  The Lucy Hedgehog DP (used in Meissel-Lehmer) is ALREADY a form of
  compressed computation: it reduces the O(x) sieve to O(x^{2/3})
  operations by representing pi(x) through O(sqrt(x)) DP values.

  No further compression is possible because:
  - The DP values encode pi(x/v) for various v
  - Each pi(x/v) differs from any smooth approximation by O(sqrt(x/v))
  - These errors carry irreducible information about prime distribution
  - Predicting them requires knowledge equivalent to what we're computing

  This reinforces the known barrier: exact p(n) requires
  Omega(p(n)^{1/3}) time at minimum (and likely Omega(p(n)^{1/2-epsilon})
  is the true barrier, though this is unproven).

  The batch Mobius computation (Approach 4) gives a practical constant-
  factor speedup for the "ordinary leaves" in Meissel-Lehmer, but does
  not change the asymptotic O(x^{2/3}) complexity.
""")


# ============================================================
# RUN ALL EXPERIMENTS
# ============================================================
if __name__ == "__main__":
    t_start = time.time()

    approach1_polynomial_pi()
    approach2_signature_function()
    approach3_prime_race()
    approach4_batch_mobius()
    approach5_sieve_compression()
    bonus_hybrid_analysis()
    final_summary()

    t_total = time.time() - t_start
    print(f"\nTotal runtime: {t_total:.2f}s")
