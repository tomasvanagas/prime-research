"""
Analytic and approximation-based methods for finding the nth prime number.

Implements:
  1. Prime Number Theorem approximations (simple, refined, Cipolla)
  2. Logarithmic integral li(x) and Riemann R function
  3. Gram series for computing R(x)
  4. Inverse R function to approximate p(n)
  5. Hybrid approach: analytic approximation + local sieve for exact result
"""

import math
from typing import Callable

# ---------------------------------------------------------------------------
# 1. Prime Number Theorem (PNT) approximations
# ---------------------------------------------------------------------------

def pnt_simple(n: int) -> float:
    """p(n) ~ n * ln(n).  Crude; only meaningful for large n."""
    if n < 2:
        return 2.0
    ln_n = math.log(n)
    return n * ln_n


def pnt_refined(n: int) -> float:
    """p(n) ~ n * (ln(n) + ln(ln(n)) - 1).  Much better for moderate n."""
    if n < 2:
        return 2.0
    ln_n = math.log(n)
    return n * (ln_n + math.log(ln_n) - 1)


def cipolla(n: int) -> float:
    """
    Cipolla's asymptotic expansion for p(n), including terms through O(1/L^2).

    Let L = ln(n), M = ln(ln(n)).  Then:
        p(n) ~ n * [ L + M - 1 + (M - 2)/L
                      - (M^2 - 6*M + 11) / (2*L^2) ]

    Reference: Cipolla (1902), also see Dusart's refinements.
    """
    if n < 6:
        # For tiny n the asymptotics are meaningless; fall back.
        _small = [2, 3, 5, 7, 11, 13]
        return float(_small[min(n, 5) - 1]) if n >= 1 else 2.0
    L = math.log(n)
    M = math.log(L)
    term0 = L + M - 1
    term1 = (M - 2) / L
    term2 = -(M * M - 6 * M + 11) / (2 * L * L)
    return n * (term0 + term1 + term2)


# ---------------------------------------------------------------------------
# 2. Logarithmic integral li(x) and Riemann R function
# ---------------------------------------------------------------------------

def li(x: float) -> float:
    """
    Logarithmic integral li(x) = integral from 0 to x of dt/ln(t).

    Uses the well-known series expansion (Ramanujan's):
        li(x) = gamma + ln|ln(x)| + sum_{k=1}^{inf} (ln(x))^k / (k * k!)

    where gamma is the Euler-Mascheroni constant.
    """
    if x <= 0:
        return 0.0
    if x == 1.0:
        return -math.inf  # pole at x=1
    gamma = 0.5772156649015329  # Euler-Mascheroni
    ln_x = math.log(x)
    result = gamma + math.log(abs(ln_x))
    term = 1.0
    for k in range(1, 200):
        term *= ln_x / k
        contrib = term / k
        result += contrib
        if abs(contrib) < 1e-15 * abs(result):
            break
    return result


def _mobius_values(limit: int) -> list[int]:
    """Compute Mobius function mu(k) for k = 0..limit via sieve."""
    # Start with mu[k] = 1 for all k >= 1, then adjust.
    mu = [0] * (limit + 1)
    for i in range(1, limit + 1):
        mu[i] = 1
    # Track whether any prime factor appears squared
    has_sq = [False] * (limit + 1)
    for p in range(2, limit + 1):
        # Check if p is prime: mu[p] should still be 1 at this point
        # if no smaller prime divides it with square. We need a separate
        # primality tracker.
        pass
    # Cleaner approach: standard Mobius sieve
    mu = [0] * (limit + 1)
    mu[1] = 1
    is_prime_flag = [True] * (limit + 1)
    primes = []
    for i in range(2, limit + 1):
        if is_prime_flag[i]:
            primes.append(i)
            mu[i] = -1  # prime: exactly one prime factor
        for p in primes:
            if i * p > limit:
                break
            is_prime_flag[i * p] = False
            if i % p == 0:
                mu[i * p] = 0  # p^2 divides i*p
                break
            else:
                mu[i * p] = -mu[i]
    return mu


# Pre-compute Mobius values up to a reasonable limit for R(x).
_MU_LIMIT = 100
_MU = _mobius_values(_MU_LIMIT)


def riemann_R(x: float, num_terms: int = 80) -> float:
    """
    Riemann R function using the Gram series:

        R(x) = 1 + sum_{k=1}^{inf} (ln x)^k / (k * k! * zeta(k+1))

    This is equivalent to:
        R(x) = sum_{k=1}^{inf} mu(k)/k * li(x^(1/k))

    but the Gram series converges faster and avoids computing many li() calls.

    We use pre-computed values of 1/(k * zeta(k+1)) via the identity
    zeta(s) for integer s.
    """
    return _gram_series(x, num_terms)


def riemann_R_mobius(x: float, num_terms: int = 60) -> float:
    """
    Riemann R function via the Mobius definition:
        R(x) = sum_{k=1}^{N} mu(k)/k * li(x^(1/k))

    Slower but a useful cross-check against the Gram series.
    """
    limit = min(num_terms, _MU_LIMIT)
    result = 0.0
    for k in range(1, limit + 1):
        if _MU[k] == 0:
            continue
        xk = x ** (1.0 / k)
        if xk <= 1.0:
            break
        result += _MU[k] / k * li(xk)
    return result


# ---------------------------------------------------------------------------
# 3. Gram series for R(x)
# ---------------------------------------------------------------------------

def _zeta_int(s: int) -> float:
    """
    Riemann zeta function at positive integer s >= 2.
    Computed by direct summation (converges fast for s >= 2).
    """
    if s == 2:
        return math.pi ** 2 / 6
    if s == 4:
        return math.pi ** 4 / 90
    total = 0.0
    for k in range(1, 10_000):
        term = 1.0 / (k ** s)
        total += term
        if term < 1e-15:
            break
    return total


def _gram_series(x: float, num_terms: int = 100) -> float:
    """
    Gram series for R(x):
        R(x) = 1 + sum_{k=1}^{inf} (ln x)^k / (k * k! * zeta(k+1))

    Very efficient: only needs log(x) and zeta at integer points.
    """
    ln_x = math.log(x)
    result = 1.0
    term = 1.0  # will accumulate ln_x^k / k!
    for k in range(1, num_terms + 1):
        term *= ln_x / k
        zk1 = _zeta_int(k + 1)
        contrib = term / (k * zk1)
        result += contrib
        if abs(contrib) < 1e-15 * abs(result):
            break
    return result


# ---------------------------------------------------------------------------
# 4. Inverting R(x) to approximate p(n)  -- Newton's method
# ---------------------------------------------------------------------------

def inverse_R(n: int, tol: float = 0.5) -> float:
    """
    Find x such that R(x) ~ n, i.e. approximate p(n).

    Uses Newton's method with R'(x) = 1 / (x * ln(x))  [leading term].
    Starting guess from Cipolla's formula.
    """
    # Starting guess
    x = cipolla(n)
    if x < 10:
        x = float(n * 3)  # safety for small n

    for _ in range(200):
        rx = riemann_R(x)
        err = rx - n
        if abs(err) < tol:
            break
        # R'(x) ~ 1/ln(x)  since R(x) ~ li(x) and d/dx li(x) = 1/ln(x)
        deriv = 1.0 / math.log(x)
        x -= err / deriv
        if x < 2:
            x = 2.0
    return x


# ---------------------------------------------------------------------------
# 5. Sieve utilities for exact prime finding
# ---------------------------------------------------------------------------

def segmented_sieve(low: int, high: int) -> list[int]:
    """
    Return all primes in [low, high] using a segmented sieve.
    Requires high - low to be manageable (< ~10^7 or so).
    """
    if high < 2:
        return []
    low = max(low, 2)
    size = high - low + 1

    # Generate small primes up to sqrt(high)
    sqrt_high = int(math.isqrt(high)) + 1
    is_small_prime = [True] * (sqrt_high + 1)
    is_small_prime[0] = is_small_prime[1] = False
    for i in range(2, int(math.isqrt(sqrt_high)) + 1):
        if is_small_prime[i]:
            for j in range(i * i, sqrt_high + 1, i):
                is_small_prime[j] = False
    small_primes = [i for i in range(2, sqrt_high + 1) if is_small_prime[i]]

    # Sieve the segment [low, high]
    is_prime_seg = bytearray(b'\x01') * size
    for p in small_primes:
        start = ((low + p - 1) // p) * p
        if start == p:
            start += p  # don't mark the prime itself
        for j in range(start - low, size, p):
            is_prime_seg[j] = 0

    # Handle edge: if low <= 1
    if low <= 1:
        for i in range(min(2 - low, size)):
            is_prime_seg[i] = 0

    return [low + i for i in range(size) if is_prime_seg[i]]


def count_primes_sieve(limit: int) -> int:
    """Count primes up to limit using a basic sieve. For moderate limits."""
    if limit < 2:
        return 0
    sieve = bytearray(b'\x01') * (limit + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(math.isqrt(limit)) + 1):
        if sieve[i]:
            sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
    return sum(sieve)


# ---------------------------------------------------------------------------
# 6. Hybrid approach: analytic estimate + local sieve
# ---------------------------------------------------------------------------

def nth_prime_hybrid(n: int) -> int:
    """
    Find the exact nth prime using:
      1. Inverse Riemann R to get a very close estimate x_est of p(n).
      2. Bracket p(n) within [x_est - margin, x_est + margin].
      3. Use a segmented sieve to find all primes in that bracket.
      4. Count primes below x_est - margin (via sieve or li estimate)
         to determine exactly which prime in the bracket is the nth.

    Practical for n up to ~10^7 or so on a typical machine, since we
    need to sieve a small interval plus count primes below the interval.
    """
    if n <= 0:
        raise ValueError("n must be >= 1")
    if n <= 6:
        return [2, 3, 5, 7, 11, 13][n - 1]

    x_est = inverse_R(n, tol=0.5)

    # Margin: for large n the error of R^{-1} is roughly sqrt(x)*ln(x).
    # We use a generous margin to be safe.
    margin = max(int(math.sqrt(x_est) * math.log(x_est)), 1000)
    # Cap the sieve window at something reasonable
    max_window = 50_000_000
    if 2 * margin > max_window:
        margin = max_window // 2

    lo = max(2, int(x_est) - margin)
    hi = int(x_est) + margin

    # Count primes below lo
    count_below = count_primes_sieve(lo - 1)

    # Sieve [lo, hi]
    primes_in_window = segmented_sieve(lo, hi)

    # The nth prime is at index (n - count_below - 1) in primes_in_window
    idx = n - count_below - 1
    if 0 <= idx < len(primes_in_window):
        return primes_in_window[idx]

    # If we missed, widen the window (shouldn't happen for reasonable n)
    # Fall back: extend upward
    hi2 = hi + margin
    primes_ext = segmented_sieve(hi + 1, hi2)
    primes_in_window.extend(primes_ext)
    if 0 <= idx < len(primes_in_window):
        return primes_in_window[idx]

    raise RuntimeError(
        f"Could not find prime({n}): estimate={x_est:.1f}, "
        f"window=[{lo},{hi2}], count_below={count_below}, "
        f"primes_found={len(primes_in_window)}"
    )


# ---------------------------------------------------------------------------
# 7. Error analysis helpers
# ---------------------------------------------------------------------------

def relative_error(approx: float, exact: int) -> float:
    """Relative error as a fraction."""
    return (approx - exact) / exact


def analyse_approximations(test_cases: dict[int, int] | None = None) -> list[dict]:
    """
    Run all approximation methods against known values and return
    a list of result dicts.
    """
    if test_cases is None:
        test_cases = {
            10: 29,
            100: 541,
            1_000: 7919,
            10_000: 104729,
            100_000: 1299709,
            1_000_000: 15485863,
        }

    methods: list[tuple[str, Callable]] = [
        ("n*ln(n)", pnt_simple),
        ("n*(ln n + ln ln n - 1)", pnt_refined),
        ("Cipolla", cipolla),
        ("Inverse R (Gram)", inverse_R),
    ]

    results = []
    for n, exact in sorted(test_cases.items()):
        for name, fn in methods:
            approx = fn(n)
            err = relative_error(approx, exact)
            results.append({
                "n": n,
                "exact": exact,
                "method": name,
                "approx": approx,
                "rel_error": err,
                "rel_error_pct": err * 100,
            })
    return results


# ---------------------------------------------------------------------------
# Main: run analysis and print results
# ---------------------------------------------------------------------------

def main():
    print("=" * 80)
    print("ANALYTIC APPROXIMATION METHODS FOR THE NTH PRIME")
    print("=" * 80)

    # ---- Approximation comparison table ----
    print("\n--- Approximation accuracy ---\n")
    results = analyse_approximations()
    current_n = None
    for r in results:
        if r["n"] != current_n:
            current_n = r["n"]
            print(f"\nn = {current_n:>12,}   exact p(n) = {r['exact']:>15,}")
            print(f"  {'Method':<28s} {'Approx':>18s}  {'Rel Error':>12s}")
            print(f"  {'-'*28} {'-'*18}  {'-'*12}")
        print(f"  {r['method']:<28s} {r['approx']:>18,.1f}  {r['rel_error_pct']:>+11.4f}%")

    # ---- Cross-check: Gram series vs Mobius definition of R(x) ----
    print("\n\n--- Cross-check: Gram series vs Mobius definition of R(x) ---\n")
    for x in [100, 1000, 1e6, 1e9]:
        rg = riemann_R(x)
        rm = riemann_R_mobius(x)
        print(f"  R({x:>12,.0f}):  Gram = {rg:>16,.6f}   Mobius = {rm:>16,.6f}   diff = {abs(rg-rm):.2e}")

    # ---- Hybrid exact results ----
    print("\n\n--- Hybrid method (exact) ---\n")
    hybrid_cases = {
        10: 29,
        100: 541,
        1_000: 7919,
        10_000: 104729,
        100_000: 1299709,
        1_000_000: 15485863,
    }
    all_correct = True
    for n, expected in sorted(hybrid_cases.items()):
        result = nth_prime_hybrid(n)
        ok = "OK" if result == expected else "FAIL"
        if result != expected:
            all_correct = False
        print(f"  prime({n:>10,}) = {result:>12,}  expected {expected:>12,}  [{ok}]")

    if all_correct:
        print("\n  All hybrid results correct.")
    else:
        print("\n  WARNING: some hybrid results are incorrect!")

    # ---- li(x) vs pi(x) vs R(x) comparison ----
    print("\n\n--- pi(x) vs li(x) vs R(x) for selected x ---\n")
    pi_values = {
        100: 25,
        1000: 168,
        10000: 1229,
        100000: 9592,
        1000000: 78498,
    }
    print(f"  {'x':>10s}  {'pi(x)':>10s}  {'li(x)':>14s}  {'R(x)':>14s}  {'li err%':>10s}  {'R err%':>10s}")
    for x, pi_x in sorted(pi_values.items()):
        li_val = li(x)
        r_val = riemann_R(x)
        li_err = (li_val - pi_x) / pi_x * 100
        r_err = (r_val - pi_x) / pi_x * 100
        print(f"  {x:>10,}  {pi_x:>10,}  {li_val:>14,.4f}  {r_val:>14,.4f}  {li_err:>+9.4f}%  {r_err:>+9.4f}%")


if __name__ == "__main__":
    main()
