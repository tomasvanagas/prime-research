"""
Experiment C: Recursive Identity Search for pi(x)

Search for identities of the form:
  pi(x) = f(pi(x/2), pi(x/3), pi(x/5), ..., pi(x/p_k)) + correction(x)

where correction(x) is efficiently computable and the recursion has
depth O(log x), giving total work O(polylog(x)).

This is different from Meissel-Lehmer which uses:
  pi(x) = pi(sqrt(x)) + phi(x, pi(sqrt(x))) - 1
because phi requires O(x^{2/3}) evaluations.

We search for identities using ONLY pi at FRACTIONAL arguments (x/2, x/3, ...)
with O(log x)-depth recursion trees.

Date: 2026-04-04
"""

import numpy as np
from sympy import primepi, isprime, nextprime, prime
from itertools import combinations
import math
from collections import Counter

# Precompute pi values using sieve (fast)
X_MAX = 100000
print(f"Precomputing pi(x) for x=0..{X_MAX} via sieve...")
_is_prime = [False, False] + [True] * (X_MAX - 1)
for _i in range(2, int(X_MAX**0.5) + 1):
    if _is_prime[_i]:
        for _j in range(_i*_i, X_MAX + 1, _i):
            _is_prime[_j] = False
pi_cache = [0] * (X_MAX + 1)
count = 0
for n in range(1, X_MAX + 1):
    if _is_prime[n]:
        count += 1
    pi_cache[n] = count
print(f"Done. pi({X_MAX}) = {pi_cache[X_MAX]}")

def pi(x):
    x = int(x)
    if x < 0:
        return 0
    if x <= X_MAX:
        return pi_cache[x]
    return int(primepi(x))

# ============================================================
# Test 1: Legendre-type identities using only pi(x/p)
# ============================================================
print("\n" + "=" * 70)
print("TEST 1: Searching for pi(x) = sum c_i * pi(x/p_i) + correction")
print("=" * 70)

# The inclusion-exclusion principle gives:
# pi(x) - pi(sqrt(x)) = #{n <= x : n has no prime factor <= sqrt(x)}
#                      = sum_{d | P(sqrt(x))} mu(d) * floor(x/d)
# But this has 2^{pi(sqrt(x))} terms.

# Instead, look for SIMPLE linear identities.
# Test: pi(x) = a*pi(x/2) + b*pi(x/3) + c*pi(x/5) + correction

# Sample many x values and fit
from numpy.linalg import lstsq

divisors_to_try = [2, 3, 5, 7, 11, 13]

for n_divisors in range(2, 7):
    divs = divisors_to_try[:n_divisors]

    # Build system: pi(x) = sum_j c_j * pi(floor(x/d_j)) + c_0
    x_samples = list(range(100, 10001, 1))
    n_samples = len(x_samples)

    A = np.zeros((n_samples, n_divisors + 1))  # +1 for constant term
    b = np.zeros(n_samples)

    for i, x in enumerate(x_samples):
        b[i] = pi(x)
        A[i, 0] = 1  # constant
        for j, d in enumerate(divs):
            A[i, j+1] = pi(x // d)

    coeffs, residuals, rank, sv = lstsq(A, b, rcond=None)

    # Compute max and mean absolute error
    predictions = A @ coeffs
    errors = np.abs(b - predictions)
    max_err = np.max(errors)
    mean_err = np.mean(errors)
    rms_err = np.sqrt(np.mean(errors**2))

    print(f"\nDivisors {divs}: coeffs = [{', '.join(f'{c:.4f}' for c in coeffs)}]")
    print(f"  Max error: {max_err:.2f}, Mean error: {mean_err:.2f}, RMS: {rms_err:.2f}")

    # Check: is the error an integer-valued function?
    int_errors = [round(b[i] - predictions[i]) for i in range(n_samples)]
    error_is_integer = all(abs(e - round(e)) < 0.01 for e in (b - predictions))
    print(f"  Errors are integers: {error_is_integer}")

    if error_is_integer:
        err_vals = Counter(int_errors)
        print(f"  Error distribution: {dict(sorted(err_vals.items())[:20])}")

# ============================================================
# Test 2: Exact recursive identities from Legendre
# ============================================================
print("\n" + "=" * 70)
print("TEST 2: Legendre's identity and its recursive structure")
print("=" * 70)
print("""
Legendre: pi(x) = pi(sqrt(x)) + phi(x, pi(sqrt(x))) - 1
where phi(x, a) = #{n <= x : all prime factors of n > p_a}

phi(x, a) = phi(x, a-1) - phi(x/p_a, a-1)
phi(x, 0) = floor(x)

The recursion tree for phi has depth pi(sqrt(x)) and the number
of DISTINCT (floor(x/d), a) pairs is O(x^{2/3}).

Can we find a DIFFERENT recursion with fewer distinct subproblems?
""")

# Count distinct subproblems in various recursion schemes
def count_legendre_subproblems(x):
    """Count distinct (floor value, depth) pairs in Legendre recursion."""
    seen = set()

    def phi_count(x_val, a, primes):
        if a == 0:
            return
        key = (x_val, a)
        if key in seen:
            return
        seen.add(key)
        phi_count(x_val, a - 1, primes)
        phi_count(x_val // primes[a-1], a - 1, primes)

    sqrt_x = int(x**0.5)
    primes = [p for p in range(2, sqrt_x + 1) if isprime(p)]
    a = len(primes)

    phi_count(x, a, primes)
    return len(seen)

# Test 3: Alternative recursion -- pi(x) from pi at smooth points only
print("\n" + "=" * 70)
print("TEST 3: Recursion using ONLY pi(x/d) for smooth d")
print("=" * 70)
print("""
Idea: Can we express pi(x) exactly using only:
  pi(x/2), pi(x/3), pi(x/4), pi(x/5), pi(x/6), ...

i.e., pi at x divided by small numbers?

By Mobius inversion-type identities:
  sum_{d=1}^{x} mu(d) * floor(x/d) = 1  (for x >= 1)

But we want pi(x), not floor(x).

Buchstab's identity: pi(x) - pi(x^{1/u}) = sum_{x^{1/u} < p <= x^{1/(u-1)}} (1 + pi(x/p) - pi(p))
This involves pi(x/p) for PRIMES p, which is circular.
""")

# Test: the "smooth recursion" approach
# pi(x) = pi(x-1) + isprime(x) -- trivial, O(x) depth
# Can we batch this?

# The fundamental question: what is the MINIMUM number of pi evaluations
# needed to determine pi(x)?

# Test with inclusion-exclusion using the first few primes
def pi_from_ie(x, k_primes):
    """
    Compute pi(x) using inclusion-exclusion over the first k_primes primes.
    Uses: #{composites with ALL factors > p_k} = x - 1 - #{composites with some factor <= p_k}

    pi(x) = -1 + sum_{d in S} mu(d) * floor(x/d)
    where S = {products of primes <= p_k}

    Plus: pi(p_k) for the small primes themselves.
    Actually: pi(x) = pi(p_k) + sum_{d | P_k} mu(d) * floor(x/d) for P_k = prod_{i<=k} p_i

    This is Legendre's sieve formula.
    """
    primes = []
    p = 2
    for _ in range(k_primes):
        primes.append(p)
        p = int(nextprime(p))

    # Generate all squarefree products of subsets of primes
    terms = 0
    from itertools import combinations
    total_floor_evals = 0

    result = 0
    for r in range(k_primes + 1):
        for subset in combinations(primes, r):
            d = 1
            for p in subset:
                d *= p
            mu_d = (-1) ** r
            result += mu_d * (x // d)
            total_floor_evals += 1

    # result = phi(x, k) = #{n <= x : gcd(n, P_k) = 1}
    # pi(x) = phi(x, k) + k - 1
    # (subtract 1 because phi counts n=1)
    result = result + k_primes - 1

    return result, total_floor_evals

print("\nLegendre sieve (inclusion-exclusion) subproblem counts:")
for x in [100, 1000, 10000, 100000]:
    sqrt_x = int(x**0.5)
    k = pi(sqrt_x)

    if 2**k > 10**8:
        print(f"x={x:>6d}: k={k}, 2^k={2**k:.2e} -- SKIPPED (too many terms)")
        continue
    pi_ie, n_evals = pi_from_ie(x, k)
    pi_true = pi(x)

    print(f"x={x:>6d}: k={k}, 2^k={2**k:>10d}, floor evals={n_evals:>10d}, "
          f"pi_ie={pi_ie}, pi_true={pi_true}, correct={pi_ie == pi_true}")

# ============================================================
# Test 4: Search for O(log x)-depth recursion
# ============================================================
print("\n" + "=" * 70)
print("TEST 4: Testing if ANY recursion on pi(x/d) achieves O(log x) depth")
print("=" * 70)
print("""
Key question: is there a function f and divisors d_1,...,d_k such that
  pi(x) = f(pi(x/d_1), ..., pi(x/d_k), x) exactly?

If the recursion is pi(x) = F(pi(x/d), x) with d >= 2, then:
  - Depth = O(log x) (good!)
  - Width at depth j: depends on branching factor

The crucial test: what is pi(x) - 2*pi(x/2)?
(This is related to prime counting in [x/2, x])
""")

# Compute pi(x) - c*pi(x/d) for various c, d
print("\nResidual analysis: pi(x) - c*pi(x/d)")
for d in [2, 3, 5]:
    residuals = []
    for x in range(100, 10001):
        r = pi(x) - pi(x // d)
        residuals.append(r)

    residuals = np.array(residuals)

    # This residual = number of primes in (x/d, x]
    # Is this simpler than pi(x)?
    # It should be approximately pi(x) * (1 - 1/d) by PNT
    # But the EXACT correction is as hard as pi(x) itself

    # Check: can we express the residual as a function of log(x)?
    x_vals = np.array(range(100, 10001), dtype=float)
    approx = x_vals * (1 - 1/d) / np.log(x_vals)

    errors = residuals - approx
    print(f"d={d}: residual pi(x)-pi(x/{d}), mean_residual={np.mean(residuals):.1f}, "
          f"mean_approx_error={np.mean(np.abs(errors)):.2f}, max_approx_error={np.max(np.abs(errors)):.1f}")

# ============================================================
# Test 5: Exhaustive search for exact simple recursions
# ============================================================
print("\n" + "=" * 70)
print("TEST 5: Exhaustive search for exact recursions with small correction")
print("=" * 70)
print("Searching: pi(x) = sum_{i} a_i * pi(x/d_i) + g(x)")
print("where g(x) depends only on x (not on primes), and a_i are integers.")

# For each candidate set of divisors, find integer coefficients
# such that the correction g(x) = pi(x) - sum a_i*pi(x/d_i) is "simple"

divisor_sets = [
    [2],
    [2, 3],
    [2, 3, 5],
    [2, 3, 5, 7],
    [2, 4],
    [2, 3, 4],
    [2, 3, 4, 5, 6],
    [2, 3, 5, 6, 10, 15, 30],  # all divisors of 30
]

coeff_range = range(-3, 4)

best_results = []

for divs in divisor_sets:
    best_complexity = float('inf')
    best_coeffs = None

    # Try small integer coefficients
    from itertools import product as iproduct

    # Limit search space
    if len(divs) <= 4:
        coeff_space = list(iproduct(coeff_range, repeat=len(divs)))
    else:
        # Too many combinations, use least squares with integer rounding
        x_samples = list(range(100, 5001))
        A_mat = np.zeros((len(x_samples), len(divs)))
        b_vec = np.zeros(len(x_samples))
        for i, x in enumerate(x_samples):
            b_vec[i] = pi(x)
            for j, d in enumerate(divs):
                A_mat[i, j] = pi(x // d)
        c_opt, _, _, _ = lstsq(A_mat, b_vec, rcond=None)
        # Round to nearest integers and check nearby
        c_round = [round(c) for c in c_opt]
        coeff_space = [tuple(c_round)]
        # Also try +-1 on each coordinate
        for idx in range(len(divs)):
            for delta in [-1, 0, 1]:
                c2 = list(c_round)
                c2[idx] += delta
                coeff_space.append(tuple(c2))
        coeff_space = list(set(coeff_space))

    for coeffs in coeff_space:
        if all(c == 0 for c in coeffs):
            continue

        # Compute correction g(x) = pi(x) - sum a_i * pi(x/d_i)
        corrections = []
        for x in range(100, 5001):
            g = pi(x) - sum(c * pi(x // d) for c, d in zip(coeffs, divs))
            corrections.append(g)

        corrections = np.array(corrections)

        # Measure "complexity" of the correction
        # Simple = small range, low entropy, or polynomial in log(x)
        unique_vals = len(set(corrections))
        val_range = np.max(corrections) - np.min(corrections)
        diff1 = np.diff(corrections)
        n_changes = np.sum(diff1 != 0)

        # Complexity score: fewer unique values and changes = simpler
        complexity = unique_vals * n_changes

        if complexity < best_complexity:
            best_complexity = complexity
            best_coeffs = coeffs
            best_corrections = corrections

    if best_coeffs is not None:
        g = best_corrections
        print(f"\nDivisors {divs}: best coeffs = {best_coeffs}")
        print(f"  Correction range: [{int(np.min(g))}, {int(np.max(g))}], "
              f"unique values: {len(set(g))}, "
              f"changes: {np.sum(np.diff(g) != 0)}/{len(g)-1}")

        # The correction g(x) is the number of primes in certain intervals
        # Check if it grows like O(sqrt(x))
        x_vals = np.array(range(100, 5001), dtype=float)
        # Fit g(x) ~ a * x^alpha
        log_x = np.log(x_vals)
        log_g = np.log(np.maximum(np.abs(g), 1))
        alpha_fit = np.polyfit(log_x, log_g, 1)
        print(f"  Growth exponent: g(x) ~ x^{alpha_fit[0]:.3f}")

# ============================================================
# Test 6: The fundamental counting argument
# ============================================================
print("\n" + "=" * 70)
print("TEST 6: Fundamental counting -- minimum subproblem tree size")
print("=" * 70)
print("""
For any recursion pi(x) = F(pi(x/d_1), ..., pi(x/d_k)):
- The leaves of the recursion tree are pi(c) for small constants c.
- The recursion tree has O(log_d(x)) depth.
- At each node, we need EXACT values (no approximation allowed).

Key insight: the set of reachable arguments is
  S = {floor(x / (d_{i1} * d_{i2} * ... * d_{ij})) : choices of divisors}

The size of S is the number of DISTINCT floor(x/n) values for n in the
set of all products of divisors. For divisors {2,3,5}, the products are
all 5-smooth numbers up to x, giving |S| = O(x^{1/3 + epsilon}) distinct values.

This is FEWER than the O(x^{2/3}) of Meissel-Lehmer, but the correction
terms at each node are NOT efficiently computable.
""")

# Count distinct reachable values for different divisor sets
def count_reachable(x, divisors, max_depth=None):
    """Count distinct floor(x/d) values reachable by dividing by divisors."""
    if max_depth is None:
        max_depth = int(math.log(x) / math.log(min(divisors))) + 1

    reachable = {x}
    frontier = {x}

    for depth in range(max_depth):
        new_frontier = set()
        for val in frontier:
            for d in divisors:
                new_val = val // d
                if new_val >= 1 and new_val not in reachable:
                    reachable.add(new_val)
                    new_frontier.add(new_val)
        frontier = new_frontier
        if not frontier:
            break

    return len(reachable), depth + 1

for x in [1000, 10000, 100000]:
    print(f"\nx = {x}:")
    for divs in [[2], [2, 3], [2, 3, 5], [2, 3, 5, 7]]:
        n_reach, depth = count_reachable(x, divs)
        # Compare with total floor values
        n_floor = count_legendre_subproblems(x) if x <= 100000 else int(2 * x**(2/3))
        print(f"  Divisors {divs}: reachable={n_reach}, depth={depth}, "
              f"vs total floor vals~{n_floor}")

def count_legendre_subproblems(x):
    values = set()
    for n in range(1, int(x**0.5) + 2):
        values.add(x // n)
        values.add(n)
    return len(values)

print("""
CONCLUSION: Recursive identities using pi(x/d) for fixed small divisors d
produce recursion trees that reach FEWER subproblems than Meissel-Lehmer
(e.g., O(x^{1/3}) instead of O(x^{2/3})).

HOWEVER, the correction function at each node of the recursion is NOT
efficiently computable. The correction encodes the number of primes
in specific intervals, which has the same information content as pi(x) itself.

In fact, any exact identity pi(x) = F(pi(x/d_1), ..., pi(x/d_k)) + g(x)
requires g(x) to encode the primes between x/max(d_i) and x, which is
essentially pi(x) - pi(x/2) ~ x/(2*ln(x)) -- a quantity as hard to compute
as pi(x).

The recursion merely SHIFTS the problem into the correction terms.
This is Failure Mode E (Equivalence) -- the difficulty is conserved.

VERDICT: Recursive identity approach is CLOSED.
""")
