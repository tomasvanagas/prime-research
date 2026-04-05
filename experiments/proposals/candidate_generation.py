#!/usr/bin/env python3
"""
PROPOSAL: Heuristic Candidate Generation + AKS/Miller-Rabin Verification

Idea: Generate a small candidate set C guaranteed to contain p(n),
then verify each candidate for primality.

Strategy:
  1. Compute x = R^{-1}(n) via Newton's method on R(x) - n
  2. The error |p(n) - x| < C * sqrt(x) * log(x)  (under RH)
  3. Apply progressive filters to narrow the candidate set:
     a. Residue class constraints (mod 30 coprime residues)
     b. Zeta zero oscillatory correction to narrow the interval
     c. Sieving out composites with small primes
  4. Measure: how small can the candidate set get while always containing p(n)?

Test for n up to 10000.
"""

import math
import time
import numpy as np
from sympy import prime, isprime, mobius, primerange
from scipy.special import expi as li_scipy

# ---------------------------------------------------------------------------
# Logarithmic integral li(x) using scipy's Ei
# ---------------------------------------------------------------------------
def li(x):
    """Offset logarithmic integral: li(x) = Ei(ln(x))."""
    if x <= 1:
        return 0.0
    return li_scipy(math.log(x))

# ---------------------------------------------------------------------------
# Riemann R(x) = sum_{k=1}^{K} mu(k)/k * li(x^{1/k})
# ---------------------------------------------------------------------------
def riemann_R(x, K=100):
    """Riemann's prime counting function approximation."""
    if x <= 1:
        return 0.0
    total = 0.0
    for k in range(1, K + 1):
        mu_k = int(mobius(k))
        if mu_k == 0:
            continue
        xk = x ** (1.0 / k)
        if xk <= 1.0:
            break
        total += mu_k / k * li(xk)
    return total

# ---------------------------------------------------------------------------
# R^{-1}(n) via Newton's method
# ---------------------------------------------------------------------------
def riemann_R_inverse(n, tol=1e-6, max_iter=200):
    """Compute R^{-1}(n): the x such that R(x) = n."""
    if n <= 1:
        return 2.0
    # Initial guess: n * ln(n)
    x = n * math.log(n)
    if x < 2:
        x = 2.0
    for _ in range(max_iter):
        rx = riemann_R(x)
        err = rx - n
        if abs(err) < tol:
            break
        # Derivative of R(x) approx 1/ln(x) for the dominant term
        drdx = 1.0 / math.log(x) if x > 1 else 1.0
        step = err / drdx
        x -= step
        if x < 2:
            x = 2.0
    return x

# ---------------------------------------------------------------------------
# First K nontrivial zeta zeros (imaginary parts)
# ---------------------------------------------------------------------------
# Precomputed first 30 zeta zeros
ZETA_ZEROS = [
    14.134725142, 21.022039639, 25.010857580, 30.424876126, 32.935061588,
    37.586178159, 40.918719012, 43.327073281, 48.005150881, 49.773832478,
    52.970321478, 56.446247697, 59.347044003, 60.831778525, 65.112544048,
    67.079810529, 69.546401711, 72.067157674, 75.704690699, 77.144840069,
    79.337375020, 82.910380854, 84.735492981, 87.425274613, 88.809111208,
    92.491899271, 94.651344041, 95.870634228, 98.831194218, 101.317851006,
]

def zeta_oscillation(x, K_zeros=20):
    """
    Estimate the oscillatory correction to pi(x) from the first K zeta zeros.
    pi(x) ~ R(x) - sum_{rho} R(x^rho)
    For each zero gamma_k, the contribution is approximately:
      -2 * Re[ R(x^{1/2 + i*gamma_k}) ]
    We approximate this as:
      -2 * x^{1/2} / log(x) * cos(gamma_k * log(x))  (dominant term)
    """
    if x <= 4:
        return 0.0
    lnx = math.log(x)
    sqrtx = math.sqrt(x)
    correction = 0.0
    for k in range(min(K_zeros, len(ZETA_ZEROS))):
        gamma = ZETA_ZEROS[k]
        correction -= 2.0 * sqrtx / lnx * math.cos(gamma * lnx)
    return correction

# ---------------------------------------------------------------------------
# Coprime residues mod 30 (the primes > 5 lie in these classes)
# ---------------------------------------------------------------------------
COPRIME_30 = {1, 7, 11, 13, 17, 19, 23, 29}  # 8 out of 30

# ---------------------------------------------------------------------------
# Small prime sieve: given an interval [lo, hi], sieve out multiples of
# small primes to get candidate primes
# ---------------------------------------------------------------------------
def sieve_interval(lo, hi, sieve_limit=1000):
    """
    Sieve the interval [lo, hi] removing multiples of primes up to sieve_limit.
    Returns sorted list of surviving candidates.
    """
    length = hi - lo + 1
    if length <= 0:
        return []
    is_candidate = bytearray(b'\x01') * length
    for p in primerange(2, sieve_limit + 1):
        # First multiple of p >= lo
        start = ((lo + p - 1) // p) * p
        if start == p:
            start += p  # don't sieve out p itself
        for j in range(start - lo, length, p):
            is_candidate[j] = 0
    # Also remove 1 if in range
    if lo <= 1:
        is_candidate[1 - lo] = 0
    candidates = []
    for i in range(length):
        if is_candidate[i]:
            val = lo + i
            if val >= 2:
                candidates.append(val)
    return candidates

# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------
def run_experiment(n_max=10000, sample_points=None, verbose=True):
    """
    For each n in sample_points (or 1..n_max), compute:
      - R^{-1}(n) = smooth approximation
      - Actual p(n)
      - Error |p(n) - R^{-1}(n)|
      - Interval width under RH bound
      - Candidate set sizes after progressive filtering
      - Whether p(n) is always in the surviving set
    """
    if sample_points is None:
        # Sample at key points for speed
        sample_points = list(range(1, min(n_max + 1, 101)))  # all of 1..100
        sample_points += list(range(200, min(n_max + 1, 1001), 100))  # 200..1000 step 100
        sample_points += list(range(2000, min(n_max + 1, 10001), 1000))  # 2000..10000 step 1000
        sample_points = sorted(set(sample_points))

    results = {
        'n': [],
        'pn': [],
        'R_inv': [],
        'error': [],
        'rh_interval_width': [],
        'cand_after_interval': [],
        'cand_after_residue': [],
        'cand_after_zeta': [],
        'cand_after_sieve': [],
        'pn_in_set': [],
        'density_ratio': [],  # candidate_set_size / interval_width
    }

    failures = []  # cases where p(n) is NOT in the candidate set
    t0 = time.time()

    print(f"Running candidate generation experiment for {len(sample_points)} values of n (max={max(sample_points)})")
    print("=" * 100)
    print(f"{'n':>6} | {'p(n)':>8} | {'R_inv':>10} | {'error':>8} | {'RH_width':>8} | "
          f"{'#interval':>9} | {'#residue':>8} | {'#zeta':>8} | {'#sieve':>8} | {'ok?':>4} | {'ratio':>6}")
    print("-" * 100)

    for n in sample_points:
        pn = prime(n)

        # Step 1: Compute R^{-1}(n)
        x = riemann_R_inverse(n)

        # Step 2: Error and interval width
        error = pn - x
        # RH bound: |pi(x) - R(x)| < (1/(8*pi)) * sqrt(x) * log(x)
        # Translating to p(n): the error in p(n) ~ sqrt(pn) * log(pn)
        # Use a generous constant to ensure containment
        if pn > 10:
            rh_width = 3.0 * math.sqrt(pn) * math.log(pn)
        else:
            rh_width = max(20, 3.0 * math.sqrt(pn) * math.log(pn))

        lo = max(2, int(math.floor(x - rh_width)))
        hi = int(math.ceil(x + rh_width))
        interval_size = hi - lo + 1

        # Step 3a: Count integers in interval (baseline)
        cand_interval = interval_size

        # Step 3b: Residue class filter (mod 30)
        # Fraction surviving: 8/30 ~ 0.267, but for small numbers include 2,3,5
        cand_residue = 0
        for v in range(lo, hi + 1):
            if v <= 5:
                if v in (2, 3, 5):
                    cand_residue += 1
            elif v % 30 in COPRIME_30:
                cand_residue += 1

        # Step 3c: Zeta zero narrowing
        # Use zeta oscillation to get a better center estimate
        zeta_corr = zeta_oscillation(x, K_zeros=20)
        x_corrected = x + 0.0  # R^{-1} already targets pi(x)=n
        # The zeta correction helps us estimate pi(x) more precisely,
        # so we can narrow the interval. The correction to p(n) is roughly:
        # delta_p ~ -zeta_corr * log(x)  (by inverting the density 1/log(x))
        x_zeta = x - zeta_corr * math.log(x) if x > 4 else x

        # Narrower interval using zeta-corrected center
        if pn > 10:
            # After zeta correction, residual error is smaller
            # Empirically scale: with K zeros, error reduces by ~1/sqrt(K) heuristically
            zeta_width = rh_width / math.sqrt(5)  # conservative
        else:
            zeta_width = rh_width

        lo_z = max(2, int(math.floor(x_zeta - zeta_width)))
        hi_z = int(math.ceil(x_zeta + zeta_width))

        cand_zeta = 0
        for v in range(lo_z, hi_z + 1):
            if v <= 5:
                if v in (2, 3, 5):
                    cand_zeta += 1
            elif v % 30 in COPRIME_30:
                cand_zeta += 1

        # Step 3d: Full sieve of the zeta-narrowed interval
        sieve_limit = min(1000, int(math.sqrt(hi_z)) + 1)
        sieved = sieve_interval(lo_z, hi_z, sieve_limit=sieve_limit)
        cand_sieve = len(sieved)

        # Check if p(n) is in the final candidate set
        pn_in_set = pn in sieved

        # Density ratio: final candidates / interval width
        density = cand_sieve / max(1, interval_size)

        results['n'].append(n)
        results['pn'].append(pn)
        results['R_inv'].append(x)
        results['error'].append(error)
        results['rh_interval_width'].append(rh_width)
        results['cand_after_interval'].append(cand_interval)
        results['cand_after_residue'].append(cand_residue)
        results['cand_after_zeta'].append(cand_zeta)
        results['cand_after_sieve'].append(cand_sieve)
        results['pn_in_set'].append(pn_in_set)
        results['density_ratio'].append(density)

        if not pn_in_set:
            failures.append(n)

        # Print selected rows
        if n <= 20 or n % 100 == 0 or n in sample_points[-5:] or not pn_in_set:
            print(f"{n:>6} | {pn:>8} | {x:>10.1f} | {error:>8.1f} | {rh_width:>8.0f} | "
                  f"{cand_interval:>9} | {cand_residue:>8} | {cand_zeta:>8} | {cand_sieve:>8} | "
                  f"{'YES' if pn_in_set else 'FAIL':>4} | {density:>6.4f}")

    elapsed = time.time() - t0

    # ---------------------------------------------------------------------------
    # Summary statistics
    # ---------------------------------------------------------------------------
    print("\n" + "=" * 100)
    print("SUMMARY")
    print("=" * 100)
    n_tested = len(sample_points)
    n_ok = sum(results['pn_in_set'])
    n_fail = n_tested - n_ok

    print(f"Values tested:     {n_tested}")
    print(f"p(n) in final set: {n_ok}/{n_tested} ({100*n_ok/n_tested:.1f}%)")
    print(f"Failures:          {n_fail}")
    if failures:
        print(f"Failed at n =      {failures[:20]}{'...' if len(failures)>20 else ''}")

    # Candidate set size statistics
    sieve_sizes = np.array(results['cand_after_sieve'])
    interval_sizes = np.array(results['cand_after_interval'])
    pn_arr = np.array(results['pn'])
    n_arr = np.array(results['n'])

    print(f"\nCandidate set sizes (after full sieve):")
    print(f"  Min:    {sieve_sizes.min()}")
    print(f"  Max:    {sieve_sizes.max()}")
    print(f"  Mean:   {sieve_sizes.mean():.1f}")
    print(f"  Median: {np.median(sieve_sizes):.1f}")

    print(f"\nInterval sizes (RH bound):")
    print(f"  Min:    {interval_sizes.min()}")
    print(f"  Max:    {interval_sizes.max()}")
    print(f"  Mean:   {interval_sizes.mean():.1f}")
    print(f"  Median: {np.median(interval_sizes):.1f}")

    # Reduction ratio at each stage
    res_ratio = np.array(results['cand_after_residue']) / np.maximum(interval_sizes, 1)
    zeta_ratio = np.array(results['cand_after_zeta']) / np.maximum(interval_sizes, 1)
    sieve_ratio = sieve_sizes / np.maximum(interval_sizes, 1)

    print(f"\nAverage reduction ratios (candidates / interval_size):")
    print(f"  After residue filter (mod 30):  {res_ratio.mean():.4f}")
    print(f"  After zeta narrowing + residue: {zeta_ratio.mean():.4f}")
    print(f"  After full sieve:               {sieve_ratio.mean():.4f}")

    # Key question: scaling of candidate set with n
    print(f"\nScaling analysis (candidate set size vs n):")
    # Fit log(candidates) vs log(n) for large n
    mask = n_arr >= 100
    if mask.sum() >= 3:
        log_n = np.log(n_arr[mask])
        log_c = np.log(sieve_sizes[mask].astype(float))
        # Linear fit in log-log
        coeffs = np.polyfit(log_n, log_c, 1)
        print(f"  log(#candidates) ~ {coeffs[0]:.3f} * log(n) + {coeffs[1]:.3f}")
        print(f"  => #candidates ~ n^{coeffs[0]:.3f}")
        print(f"  (polylog would need exponent -> 0)")

    # Also fit against sqrt(p(n)) * log(p(n)) / log(p(n))^2 to see if
    # candidates scale as expected from prime density in the interval
    if mask.sum() >= 3:
        expected_primes = 2.0 * np.array(results['rh_interval_width'])[mask] / np.log(pn_arr[mask].astype(float))
        actual = sieve_sizes[mask].astype(float)
        ratio = actual / expected_primes
        print(f"\n  Ratio actual_candidates / expected_primes_in_interval:")
        print(f"    Mean:   {ratio.mean():.3f}")
        print(f"    (should be ~1.0 if sieve works well, >1 means composites survive)")

    # The real question: is the candidate set O(polylog(n))?
    print(f"\nCRITICAL ASSESSMENT:")
    if mask.sum() >= 3:
        # candidates ~ sqrt(p(n)) * log(p(n)) / log(p(n)) = sqrt(p(n))
        # p(n) ~ n*log(n), so candidates ~ sqrt(n*log(n)) ~ sqrt(n)
        print(f"  The candidate set after sieving has ~pi(interval) primes")
        print(f"  Interval width ~ sqrt(p(n)) * log(p(n)) ~ sqrt(n*log(n)) * log(n*log(n))")
        print(f"  Number of primes in interval ~ sqrt(p(n)) ~ sqrt(n*log(n))")
        print(f"  This is O(sqrt(n)), NOT O(polylog(n))")
        print(f"  The bottleneck: we cannot narrow the interval below O(sqrt(p(n)))")
        print(f"  without computing the zeta zeros' contribution, which is the SAME")
        print(f"  as computing p(n) directly.")

    print(f"\nTotal runtime: {elapsed:.2f}s")

    return results, failures


def run_containment_stress_test(n_range=(1, 1000)):
    """
    For every n in range, check if p(n) is always in the candidate set.
    Uses tighter intervals to find the MINIMUM width needed.
    """
    print(f"\n{'='*80}")
    print("STRESS TEST: Finding minimum interval width that always contains p(n)")
    print(f"{'='*80}")

    multipliers = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]
    containment = {m: 0 for m in multipliers}
    total = 0

    for n in range(max(1, n_range[0]), n_range[1] + 1):
        pn = prime(n)
        x = riemann_R_inverse(n)

        # Zeta-corrected center
        zeta_corr = zeta_oscillation(x, K_zeros=20)
        x_zeta = x - zeta_corr * math.log(max(x, 2.1)) if x > 4 else x

        for m in multipliers:
            if pn > 10:
                width = m * math.sqrt(pn) * math.log(pn)
            else:
                width = max(10, m * math.sqrt(pn) * math.log(pn))

            # Check if p(n) is within [x_zeta - width, x_zeta + width]
            if abs(pn - x_zeta) <= width:
                containment[m] += 1

        total += 1

    print(f"\nTested n = {n_range[0]} to {n_range[1]} ({total} values)")
    print(f"\n{'Multiplier C':>15} | {'Contained':>10} | {'Rate':>8} | Width formula: C*sqrt(p(n))*log(p(n))")
    print("-" * 65)
    for m in multipliers:
        rate = containment[m] / total * 100
        marker = " <-- minimum for 100%" if rate == 100 and (m == multipliers[0] or containment[multipliers[multipliers.index(m)-1]] < total) else ""
        print(f"{m:>15.1f} | {containment[m]:>10} | {rate:>7.1f}% |{marker}")

    # Also test: minimum absolute error
    print(f"\nAbsolute error analysis (|p(n) - R_inv_zeta(n)|):")
    errors = []
    rel_errors = []
    for n in range(max(10, n_range[0]), n_range[1] + 1):
        pn = prime(n)
        x = riemann_R_inverse(n)
        zeta_corr = zeta_oscillation(x, K_zeros=20)
        x_zeta = x - zeta_corr * math.log(max(x, 2.1)) if x > 4 else x
        errors.append(abs(pn - x_zeta))
        rel_errors.append(abs(pn - x_zeta) / (math.sqrt(pn) * math.log(pn)))

    errors = np.array(errors)
    rel_errors = np.array(rel_errors)
    print(f"  Max |error|:                    {errors.max():.1f}")
    print(f"  Mean |error|:                   {errors.mean():.1f}")
    print(f"  Max |error|/sqrt(p)*log(p):     {rel_errors.max():.4f}")
    print(f"  Mean |error|/sqrt(p)*log(p):    {rel_errors.mean():.4f}")
    print(f"  99th percentile ratio:          {np.percentile(rel_errors, 99):.4f}")


if __name__ == "__main__":
    # Main experiment
    results, failures = run_experiment(n_max=10000)

    # Stress test for containment
    run_containment_stress_test(n_range=(1, 1000))
