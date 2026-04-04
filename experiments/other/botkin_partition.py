"""
BOTKIN-DAWSEY-HEMMER-JUST-SCHNEIDER PARTITION-THEORETIC PRIME MODEL
====================================================================
arXiv:2501.00580v3 (Dec 2024, revised Oct 2025)

Key formula (for n >= 2):
  p_n = 1 + 2 * sum_{j=1}^{n-1} ceil(d(j)/2) + epsilon(n)

where d(k) is the divisor function, epsilon(n) is an error term.

Equivalent formulation (Remark 2.2, Eq 2.4):
  p_n = 1 + sum_{k=1}^{n-1} floor((n-1)/k) + floor(sqrt(n-1))

  Because: 2*ceil(d(k)/2) = d(k) when k is not a perfect square,
           2*ceil(d(k)/2) = d(k)+1 when k is a perfect square.
  And sum_{j=1}^{N} d(j) = sum_{k=1}^{N} floor(N/k).
  Number of perfect squares <= N is floor(sqrt(N)).

Models:
  Model 1: epsilon(n) = 0
  Model 2: epsilon(n) = floor(pi_2(p_{n-1}) - 2*gamma*(n-1))
           where pi_2(x) = semiprime counting function
  Model 2*: epsilon(n) = floor((n-1)*(log log(n-1) - 2*gamma))
           (uses Landau asymptotic for pi_2)

This script tests all models and measures accuracy.
"""

import math
import time
import sys

# Euler-Mascheroni constant
GAMMA = 0.5772156649015329


def sieve_primes(limit):
    """Sieve of Eratosthenes returning list of primes up to limit."""
    if limit < 2:
        return []
    is_prime = bytearray(b'\x01') * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return [i for i in range(2, limit + 1) if is_prime[i]]


def divisor_function_sieve(limit):
    """Compute d(k) for k = 0..limit using a sieve."""
    d = [0] * (limit + 1)
    for i in range(1, limit + 1):
        for j in range(i, limit + 1, i):
            d[j] += 1
    return d


def semiprime_count_sieve(limit, primes):
    """Count semiprimes (products of exactly 2 primes, with multiplicity) up to limit.
    Returns array pi2 where pi2[n] = number of semiprimes <= n."""
    is_semiprime = bytearray(limit + 1)
    for i, p in enumerate(primes):
        if p * p > limit:
            break
        for j in range(i, len(primes)):
            q = primes[j]
            val = p * q
            if val > limit:
                break
            is_semiprime[val] = 1
    # Build cumulative count
    pi2 = [0] * (limit + 1)
    for i in range(1, limit + 1):
        pi2[i] = pi2[i - 1] + is_semiprime[i]
    return pi2


# ============================================================
# MODEL 1: Main term only, epsilon = 0
# ============================================================

def model1_main_term(n, div_func=None):
    """
    M(n) = 1 + 2 * sum_{j=1}^{n-1} ceil(d(j)/2)

    Equivalent to: 1 + sum_{k=1}^{n-1} floor((n-1)/k) + floor(sqrt(n-1))
    """
    if n == 1:
        return 2
    # Use the floor-sum formulation for efficiency
    N = n - 1
    total = 0
    for k in range(1, N + 1):
        total += N // k
    total += int(math.isqrt(N))
    return 1 + total


def model1_main_term_fast(n):
    """
    Compute 1 + sum_{k=1}^{n-1} floor((n-1)/k) + floor(sqrt(n-1))
    using O(sqrt(n)) trick for the Dirichlet hyperbola method.
    """
    if n == 1:
        return 2
    N = n - 1
    sqrtN = int(math.isqrt(N))
    total = 0
    # sum_{k=1}^{N} floor(N/k) = 2*sum_{k=1}^{sqrtN} floor(N/k) - sqrtN^2
    for k in range(1, sqrtN + 1):
        total += N // k
    total = 2 * total - sqrtN * sqrtN
    total += sqrtN  # add floor(sqrt(N))
    return 1 + total


def model1_batch(max_n):
    """Compute Model 1 predictions for n=1..max_n using incremental updates."""
    predictions = [0, 2]  # index 0 unused, p_1 = 2
    # M(n) = 1 + sum_{k=1}^{n-1} floor((n-1)/k) + floor(sqrt(n-1))
    # We compute these using the fast method
    for n in range(2, max_n + 1):
        predictions.append(model1_main_term_fast(n))
    return predictions


# ============================================================
# MODEL 2*: With asymptotic error correction
# ============================================================

def model2star_batch(max_n):
    """
    Model 2*: p_n = M(n) + floor((n-1)*(log(log(n-1)) - 2*gamma))
    for n >= 3. p_1 = 2, p_2 = 3.
    """
    predictions = [0, 2, 3]  # index 0 unused, p_1=2, p_2=3

    for n in range(3, max_n + 1):
        main = model1_main_term_fast(n)
        nm1 = n - 1
        loglog = math.log(math.log(nm1))
        eps = math.floor(nm1 * (loglog - 2 * GAMMA))
        predictions.append(main + eps)
    return predictions


# ============================================================
# MODEL 2: With semiprime-based error correction
# ============================================================

def model2_batch(max_n, primes, pi2):
    """
    Model 2: p_n = M(n) + floor(pi_2(p_{n-1}) - 2*gamma*(n-1))
    Uses actual semiprime counts (requires knowing actual primes).
    """
    predictions = [0, 2, 3]  # index 0 unused
    for n in range(3, max_n + 1):
        main = model1_main_term_fast(n)
        p_prev = primes[n - 2]  # p_{n-1}, 0-indexed
        if p_prev < len(pi2):
            sp_count = pi2[p_prev]
        else:
            sp_count = 0
        eps = math.floor(sp_count - 2 * GAMMA * (n - 1))
        predictions.append(main + eps)
    return predictions


# ============================================================
# R^{-1}(n) baseline comparison
# ============================================================

def li(x):
    """Logarithmic integral li(x) = integral from 0 to x of dt/ln(t)."""
    if x <= 1:
        return 0.0
    ln_x = math.log(x)
    r = GAMMA + math.log(abs(ln_x))
    t = 1.0
    for k in range(1, 200):
        t *= ln_x / k
        r += t / k
        if abs(t / k) < 1e-15:
            break
    return r


def R_function(x):
    """Riemann R function: R(x) = sum_{k=1}^{inf} mu(k)/k * li(x^{1/k})."""
    MU = [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1, 0, -1, 0, -1, 0,
          1, 1, -1, 0, 0, 1, 0, 0, -1, -1, -1, 0, 1, 1, 1, 0, -1, 1, 1, 0, -1, -1,
          -1, 0, 0, -1, -1, 0, 0, 0]
    total = 0.0
    for k in range(1, len(MU)):
        if MU[k] == 0:
            continue
        root = x ** (1.0 / k)
        if root <= 1.001:
            break
        total += MU[k] / k * li(root)
    return total


def R_inverse(n):
    """Approximate R^{-1}(n) by bisection: find x such that R(x) = n."""
    if n <= 0:
        return 2
    # Initial estimate from PNT: p_n ~ n*ln(n)
    est = max(3, n * math.log(max(n, 2)))
    lo, hi = 2, est * 3
    for _ in range(100):
        mid = (lo + hi) / 2
        if R_function(mid) < n:
            lo = mid
        else:
            hi = mid
        if hi - lo < 0.5:
            break
    return round((lo + hi) / 2)


# ============================================================
# MAIN EXPERIMENT
# ============================================================

def run_experiments():
    MAX_N = 10000

    print("=" * 72)
    print("BOTKIN et al. PARTITION-THEORETIC PRIME MODEL")
    print("arXiv:2501.00580v3")
    print("=" * 72)

    # Generate reference primes
    print("\n[1] Generating reference primes...")
    t0 = time.time()
    # p_10000 ~ 104729, need primes up to ~120000
    primes = sieve_primes(120000)
    t_sieve = time.time() - t0
    print(f"    Sieved {len(primes)} primes in {t_sieve:.4f}s")
    print(f"    p(1)={primes[0]}, p(10000)={primes[MAX_N-1]}")

    # Generate semiprime counts for Model 2
    print("\n[2] Computing semiprime counts...")
    t0 = time.time()
    pi2 = semiprime_count_sieve(primes[MAX_N - 1] + 1, primes)
    t_pi2 = time.time() - t0
    print(f"    Semiprime sieve up to {primes[MAX_N-1]} in {t_pi2:.4f}s")

    # ---- MODEL 1 ----
    print("\n[3] Computing Model 1 (epsilon=0)...")
    t0 = time.time()
    m1 = model1_batch(MAX_N)
    t_m1 = time.time() - t0
    print(f"    Computed in {t_m1:.4f}s")

    exact_m1 = 0
    max_err_m1 = 0
    sum_rel_err_m1 = 0.0
    for i in range(1, MAX_N + 1):
        actual = primes[i - 1]
        pred = m1[i]
        if pred == actual:
            exact_m1 += 1
        err = abs(pred - actual)
        max_err_m1 = max(max_err_m1, err)
        sum_rel_err_m1 += abs(pred - actual) / actual

    print(f"    Exact matches: {exact_m1}/{MAX_N}")
    print(f"    Max absolute error: {max_err_m1}")
    print(f"    Mean relative error: {sum_rel_err_m1/MAX_N:.6f}")

    # Show first 20 comparisons
    print("\n    First 20 values (Model 1):")
    print(f"    {'n':>5} {'actual':>8} {'Model1':>8} {'error':>8}")
    for i in range(1, 21):
        actual = primes[i - 1]
        pred = m1[i]
        print(f"    {i:>5} {actual:>8} {pred:>8} {pred - actual:>8}")

    # ---- MODEL 2* ----
    print("\n[4] Computing Model 2* (asymptotic epsilon)...")
    t0 = time.time()
    m2s = model2star_batch(MAX_N)
    t_m2s = time.time() - t0
    print(f"    Computed in {t_m2s:.4f}s")

    exact_m2s = 0
    max_err_m2s = 0
    sum_rel_err_m2s = 0.0
    for i in range(1, MAX_N + 1):
        actual = primes[i - 1]
        pred = m2s[i]
        if pred == actual:
            exact_m2s += 1
        err = abs(pred - actual)
        max_err_m2s = max(max_err_m2s, err)
        sum_rel_err_m2s += abs(pred - actual) / actual

    print(f"    Exact matches: {exact_m2s}/{MAX_N}")
    print(f"    Max absolute error: {max_err_m2s}")
    print(f"    Mean relative error: {sum_rel_err_m2s/MAX_N:.6f}")

    print("\n    First 20 values (Model 2*):")
    print(f"    {'n':>5} {'actual':>8} {'Model2*':>8} {'error':>8}")
    for i in range(1, 21):
        actual = primes[i - 1]
        pred = m2s[i]
        print(f"    {i:>5} {actual:>8} {pred:>8} {pred - actual:>8}")

    # ---- MODEL 2 (with actual semiprimes) ----
    print("\n[5] Computing Model 2 (semiprime-based epsilon)...")
    t0 = time.time()
    m2 = model2_batch(MAX_N, primes, pi2)
    t_m2 = time.time() - t0
    print(f"    Computed in {t_m2:.4f}s")

    exact_m2 = 0
    max_err_m2 = 0
    sum_rel_err_m2 = 0.0
    for i in range(1, MAX_N + 1):
        actual = primes[i - 1]
        pred = m2[i]
        if pred == actual:
            exact_m2 += 1
        err = abs(pred - actual)
        max_err_m2 = max(max_err_m2, err)
        sum_rel_err_m2 += abs(pred - actual) / actual

    print(f"    Exact matches: {exact_m2}/{MAX_N}")
    print(f"    Max absolute error: {max_err_m2}")
    print(f"    Mean relative error: {sum_rel_err_m2/MAX_N:.6f}")

    print("\n    First 20 values (Model 2):")
    print(f"    {'n':>5} {'actual':>8} {'Model2':>8} {'error':>8}")
    for i in range(1, 21):
        actual = primes[i - 1]
        pred = m2[i]
        print(f"    {i:>5} {actual:>8} {pred:>8} {pred - actual:>8}")

    # ---- R^{-1}(n) COMPARISON ----
    print("\n[6] Computing R^{-1}(n) for comparison...")
    t0 = time.time()
    r_inv = [0]
    for i in range(1, MAX_N + 1):
        r_inv.append(R_inverse(i))
    t_rinv = time.time() - t0
    print(f"    Computed in {t_rinv:.4f}s")

    exact_rinv = 0
    max_err_rinv = 0
    sum_rel_err_rinv = 0.0
    for i in range(1, MAX_N + 1):
        actual = primes[i - 1]
        pred = r_inv[i]
        if pred == actual:
            exact_rinv += 1
        err = abs(pred - actual)
        max_err_rinv = max(max_err_rinv, err)
        sum_rel_err_rinv += abs(pred - actual) / actual

    print(f"    Exact matches: {exact_rinv}/{MAX_N}")
    print(f"    Max absolute error: {max_err_rinv}")
    print(f"    Mean relative error: {sum_rel_err_rinv/MAX_N:.6f}")

    # ---- COMPARISON TABLE ----
    print("\n" + "=" * 72)
    print("COMPARISON SUMMARY")
    print("=" * 72)
    print(f"{'Model':<20} {'Exact':>8} {'MaxErr':>8} {'MeanRelErr':>12} {'Time(s)':>10}")
    print("-" * 60)
    print(f"{'Model 1 (eps=0)':<20} {exact_m1:>8} {max_err_m1:>8} {sum_rel_err_m1/MAX_N:>12.6f} {t_m1:>10.4f}")
    print(f"{'Model 2* (asymp)':<20} {exact_m2s:>8} {max_err_m2s:>8} {sum_rel_err_m2s/MAX_N:>12.6f} {t_m2s:>10.4f}")
    print(f"{'Model 2 (semiprime)':<20} {exact_m2:>8} {max_err_m2:>8} {sum_rel_err_m2/MAX_N:>12.6f} {t_m2:>10.4f}")
    print(f"{'R^(-1)(n)':<20} {exact_rinv:>8} {max_err_rinv:>8} {sum_rel_err_rinv/MAX_N:>12.6f} {t_rinv:>10.4f}")

    # ---- ERROR DISTRIBUTION for Model 2* ----
    print("\n" + "=" * 72)
    print("ERROR DISTRIBUTION (Model 2*)")
    print("=" * 72)
    err_buckets = {}
    for i in range(1, MAX_N + 1):
        actual = primes[i - 1]
        pred = m2s[i]
        err = pred - actual
        bucket = err
        err_buckets[bucket] = err_buckets.get(bucket, 0) + 1

    sorted_errs = sorted(err_buckets.items())
    # Show top error buckets
    print(f"  {'Error':>8} {'Count':>8} {'Pct':>8}")
    for err, cnt in sorted_errs[:15]:
        print(f"  {err:>8} {cnt:>8} {100*cnt/MAX_N:>7.2f}%")
    if len(sorted_errs) > 30:
        print(f"  ... ({len(sorted_errs)} distinct error values total)")
    for err, cnt in sorted_errs[-15:]:
        print(f"  {err:>8} {cnt:>8} {100*cnt/MAX_N:>7.2f}%")

    # ---- PRIME GAP PREDICTION ----
    print("\n" + "=" * 72)
    print("GAP PREDICTION TEST (Model 1)")
    print("=" * 72)
    print("Model 1 predicts: p_{n+1} - p_n = 2*ceil(d(n)/2)")
    print("  = d(n) if n is not a perfect square")
    print("  = d(n)+1 if n is a perfect square")
    print()

    # Compute divisor function
    div_limit = MAX_N + 1
    d = divisor_function_sieve(div_limit)

    gap_exact = 0
    for i in range(1, min(MAX_N, len(primes) - 1)):
        actual_gap = primes[i] - primes[i - 1]
        n = i  # gap between p_n and p_{n+1}
        if n < div_limit:
            sq = int(math.isqrt(n))
            predicted_gap = d[n] + (1 if sq * sq == n else 0)
            if predicted_gap == actual_gap:
                gap_exact += 1

    print(f"    Gap predictions correct: {gap_exact}/{min(MAX_N, len(primes)-1)}")

    # ---- COMPLEXITY ANALYSIS ----
    print("\n" + "=" * 72)
    print("COMPLEXITY ANALYSIS")
    print("=" * 72)
    print()
    print("Model 1 main term: M(n) = 1 + sum_{k=1}^{n-1} floor((n-1)/k) + floor(sqrt(n-1))")
    print("  Using hyperbola method: O(sqrt(n)) per evaluation")
    print("  Total for n terms: O(n*sqrt(n))")
    print()
    print("Model 2* correction: epsilon(n) = floor((n-1)*(log log(n-1) - 2*gamma))")
    print("  This is O(1) per evaluation -- trivially computable")
    print("  BUT: this is an APPROXIMATION. Not exact.")
    print()
    print("Model 2 correction: epsilon(n) = floor(pi_2(p_{n-1}) - 2*gamma*(n-1))")
    print("  Requires knowing p_{n-1} (circular!) and computing pi_2(p_{n-1})")
    print("  pi_2(x) can be computed in O(x/log(x)) via sieve, or ~O(x^{2/3}) analytically")
    print("  NOT O(polylog(n)) -- requires knowledge of actual primes")
    print()
    print("KEY FINDING: The error term epsilon(n) CANNOT give exact p_n values")
    print("  because it depends on knowing p_{n-1} (Model 2) or is approximate (Model 2*).")
    print("  The formula is a MODEL, not a closed-form solution.")
    print()

    # ---- SCALING TEST ----
    print("=" * 72)
    print("SCALING: Time to compute M(n) for single values")
    print("=" * 72)
    for test_n in [100, 1000, 10000, 100000]:
        t0 = time.time()
        val = model1_main_term_fast(test_n)
        dt = time.time() - t0
        actual = primes[test_n - 1] if test_n <= len(primes) else "?"
        print(f"  n={test_n:>7}: M(n)={val:>12}, actual p_n={actual}, time={dt:.6f}s")

    # ---- VERIFY FORMULA EQUIVALENCE ----
    print("\n" + "=" * 72)
    print("FORMULA EQUIVALENCE CHECK")
    print("=" * 72)
    print("Verifying: 2*sum ceil(d(j)/2) == sum floor((n-1)/k) + floor(sqrt(n-1))")
    d_vals = divisor_function_sieve(100)
    all_match = True
    for n in range(2, 100):
        # Method 1: explicit divisor sum
        s1 = sum(2 * math.ceil(d_vals[j] / 2) for j in range(1, n))
        # Method 2: floor sum
        N = n - 1
        sqN = int(math.isqrt(N))
        s2 = sum(N // k for k in range(1, N + 1)) + sqN
        if s1 != s2:
            print(f"  MISMATCH at n={n}: {s1} vs {s2}")
            all_match = False
    print(f"  All match for n=2..99: {all_match}")


if __name__ == "__main__":
    run_experiments()
