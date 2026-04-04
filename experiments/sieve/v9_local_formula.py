"""
V9: LOCAL PRIME COUNTING — A Novel Hybrid Approach
====================================================

KEY INSIGHT: Instead of computing pi(x) GLOBALLY (O(x^{2/3})),
what if we could determine the LOCAL prime count in a small interval
around R^{-1}(n)?

The approach:
1. x0 = R^{-1}(n) — O(polylog) — initial estimate
2. Sieve a small interval [x0 - W, x0 + W] — O(W * sqrt(x0))
3. Count primes in the interval using the sieve
4. Use the ANALYTIC approximation of pi(x0 - W) to determine the offset
5. Combine to find p(n) exactly

The bottleneck is step 2: sieving requires primes up to sqrt(x0),
which themselves need to be found.

NOVEL TRICK: Use a segmented sieve with PRIMALITY TESTING instead
of a full sieve. For each candidate in [x0-W, x0+W], test primality
using Miller-Rabin. This costs O(W * log^3(x0)).

Then: the number of primes found = pi(x0+W) - pi(x0-W)
And we need: pi(x0-W) to determine the offset.

Under RH: pi(x0) ≈ R(x0) = n (by construction of x0)
           pi(x0-W) ≈ R(x0-W) ≈ n - W/ln(x0)
           Error: |pi(x0-W) - R(x0-W)| ≤ C*sqrt(x0)*ln(x0)

For this error to be < 1: C*sqrt(x0)*ln(x0) < 1
i.e., x0 < 1/C^2, which is only possible for TINY x.

SO THIS DOESN'T WORK FOR LARGE n.

BUT: For MODERATE n (up to ~10^6), R(x) approximates pi(x) well
enough that the local sieve approach gives exact results!

Let's find the LARGEST n for which this "formula" works in < 1 second.
"""

import math
import time

# ============================================================
# Fast R function and inverse
# ============================================================

_MU = [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1, 0, -1, 0, -1, 0,
       1, 1, -1, 0, 0, 1, 0, 0, -1, -1, -1, 0, 1, 1, 1, 0, -1, 1, 1, 0, -1, -1,
       -1, 0, 0, -1, -1, 0, 0, 0]

def _li(x):
    if x <= 1: return 0.0
    ln_x = math.log(x)
    r = 0.5772156649015329 + math.log(abs(ln_x))
    t = 1.0
    for k in range(1, 200):
        t *= ln_x / k
        r += t / k
        if abs(t / k) < 1e-15: break
    return r

def R(x):
    if x <= 1: return 0.0
    r = 0.0
    for k in range(1, len(_MU)):
        if _MU[k] == 0: continue
        xk = x ** (1.0 / k)
        if xk <= 1.0001: break
        r += _MU[k] / k * _li(xk)
    return r

def R_inv(n):
    if n <= 5: return [0, 2, 3, 5, 7, 11][n]
    x = float(n) * math.log(n) + float(n) * math.log(math.log(n))
    for _ in range(100):
        rx = R(x)
        dx = (n - rx) * math.log(x)
        x += dx
        if abs(dx) < 1e-10: break
    return x

# ============================================================
# Fast Miller-Rabin
# ============================================================

def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    if n < 25: return True
    witnesses = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    d = n - 1
    r = 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in witnesses:
        if a >= n: continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1: continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1: break
        else:
            return False
    return True

# ============================================================
# Lucy DP for exact pi(x)
# ============================================================

def lucy_pi(x):
    x = int(x)
    if x < 2: return 0
    if x < 3: return 1
    sqrtx = int(x**0.5)
    small = list(range(-1, sqrtx + 1))  # small[v] = v - 1
    large = [0] * (sqrtx + 2)
    for v in range(1, sqrtx + 1):
        large[v] = x // v - 1
    for p in range(2, sqrtx + 1):
        if small[p] == small[p - 1]: continue
        pcnt = small[p - 1]
        p2 = p * p
        for v in range(1, min(sqrtx, x // p2) + 1):
            d = v * p
            if d <= sqrtx:
                large[v] -= large[d] - pcnt
            else:
                large[v] -= small[x // d] - pcnt
        for v in range(sqrtx, p2 - 1, -1):
            small[v] -= small[v // p] - pcnt
    return large[1]

# ============================================================
# V9: The "Local Formula" approach
# ============================================================

def nth_prime_local(n):
    """
    Attempt to find p(n) using LOCAL information only:
    1. x0 = R^{-1}(n) — very close to p(n)
    2. Test candidates near x0 for primality
    3. Use R(x) to estimate pi(x0) and determine which nearby prime is p(n)

    This works when |R(x0) - pi(x0)| < 0.5, i.e., when R is exact for pi.
    """
    if n <= 15:
        return [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47][n-1]

    x0 = R_inv(n)
    x0_int = int(round(x0))

    # Estimate pi(x0) using R(x0)
    pi_estimate = R(x0)
    offset = round(n - pi_estimate)  # how many primes we're off

    # Search window: based on expected error of R
    # Error of R(x) is roughly O(sqrt(x) * ln(x)) but in practice much smaller
    # For x < 10^8, error < 30 typically
    window = max(int(math.sqrt(x0_int) * math.log(x0_int) * 0.1), 200)
    window = min(window, 100000)  # cap the window

    # Find all primes in [x0 - window, x0 + window]
    lo = max(2, x0_int - window)
    hi = x0_int + window

    # Collect primes in the window
    primes_in_window = []
    for candidate in range(lo, hi + 1):
        if is_prime(candidate):
            primes_in_window.append(candidate)

    if not primes_in_window:
        return None  # shouldn't happen

    # The prime p(n) should be at position (n - pi(lo-1)) in this list
    # pi(lo-1) ≈ R(lo-1)
    pi_lo = R(lo - 1)
    target_idx = round(n - pi_lo) - 1  # 0-indexed

    if 0 <= target_idx < len(primes_in_window):
        return primes_in_window[target_idx]
    else:
        return None  # need larger window or exact pi


def nth_prime_v9(n):
    """
    V9 solver: Try local formula first, fall back to full Lucy DP.
    """
    if n <= 15:
        return [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47][n-1]

    # Try local approach first
    result = nth_prime_local(n)
    if result is not None:
        return result

    # Fallback to Lucy DP approach
    x0 = R_inv(n)
    x = int(round(x0))
    if x % 2 == 0: x += 1

    for _ in range(20):
        pi_x = lucy_pi(x)
        if pi_x == n:
            if is_prime(x): return x
            while not is_prime(x): x -= 1
            return x
        elif pi_x < n:
            step = max(int((n - pi_x) * math.log(x)), 2)
            x += step
            if x % 2 == 0: x += 1
        else:
            step = max(int((pi_x - n) * math.log(x)), 2)
            x -= step
            if x < 3: x = 3
            if x % 2 == 0: x += 1

    return x


# ============================================================
# TEST
# ============================================================

def sieve_ref(limit):
    s = bytearray(b'\x01') * (limit + 1)
    s[0] = s[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if s[i]:
            s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(2, limit+1) if s[i]]

if __name__ == "__main__":
    print("=" * 60)
    print("V9: LOCAL PRIME COUNTING FORMULA")
    print("=" * 60)

    primes = sieve_ref(200000)
    print(f"Sieved {len(primes)} reference primes")

    # Test local formula accuracy
    print("\n--- Local formula accuracy ---")
    for max_n in [100, 500, 1000, 5000, 10000]:
        exact = 0
        for n in range(1, max_n + 1):
            result = nth_prime_local(n)
            if result == primes[n-1]:
                exact += 1
        print(f"  n=1..{max_n}: {exact}/{max_n} exact ({100*exact/max_n:.1f}%)")

    # How far does local formula work?
    print("\n--- Finding the boundary where local formula fails ---")
    first_fail = None
    for n in range(1, min(len(primes), 17001)):
        result = nth_prime_local(n)
        if result != primes[n-1]:
            print(f"  First failure at n={n}: got {result}, expected {primes[n-1]}")
            first_fail = n
            break
    if first_fail is None:
        print(f"  Local formula correct for all n up to {min(len(primes), 17000)}!")

    # Check if failures are systematic
    if first_fail and first_fail < 10000:
        fails = []
        for n in range(first_fail, min(first_fail + 100, len(primes))):
            result = nth_prime_local(n)
            if result != primes[n-1]:
                fails.append(n)
        print(f"  Failures in [{first_fail}, {first_fail+99}]: {len(fails)}/100")

    # Benchmark v9
    print("\n--- V9 Benchmark ---")
    for n_val in [10, 100, 1000, 10000, 100000]:
        t0 = time.time()
        result = nth_prime_v9(n_val)
        dt = time.time() - t0
        correct = result == primes[n_val - 1] if n_val <= len(primes) else "?"
        print(f"  p({n_val:>6d}) = {result:>10d}  ({dt:.4f}s) {'OK' if correct else 'FAIL'}")

    # The critical experiment: for what range does the LOCAL formula
    # (using R(x) as pi(x) proxy) give exact results?
    print("\n--- Critical test: R(x) as pi(x) proxy ---")
    for x_test in [100, 1000, 10000, 100000, 1000000]:
        pi_true = lucy_pi(x_test)
        pi_approx = R(x_test)
        error = pi_approx - pi_true
        print(f"  x={x_test:>10d}: pi={pi_true:>7d}, R(x)={pi_approx:>10.2f}, error={error:+.2f}")

    print("\n" + "=" * 60)
    print("ANALYSIS")
    print("=" * 60)
    print("""
The local formula works when |R(x) - pi(x)| < 0.5 at the relevant points.
For small x, R(x) approximates pi(x) well but not exactly.

The fundamental limitation: R(x) - pi(x) oscillates with amplitude
O(sqrt(x)*ln(x)/(8*pi)). This exceeds 0.5 around x ≈ 100-200.

So the local formula CANNOT work for all n — it fails when the
R function's approximation of pi crosses the 0.5 error threshold.

For a TRUE formula that works for all n, we need EXACT pi(x),
which requires O(x^{2/3}) computation (Lucy DP).
""")
