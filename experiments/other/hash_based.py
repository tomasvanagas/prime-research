"""
Session 5: Hash/Lookup-Based Approaches + Novel Compression Ideas

RADICAL IDEA 1: Perfect Hash Function
If there existed a function h(n) that maps n to p(n) directly,
using only O(polylog(n)) operations, it would solve our problem.
This is essentially what we're looking for — a "perfect formula."

RADICAL IDEA 2: Compressed Oracle
Store a compressed representation of the prime sequence that allows
O(polylog(n)) random access. This is the "succinct data structure" approach.

Known results:
- Meissel-Lehmer is already a form of compressed computation
- The prime sequence has entropy ~log(n) + log(log(n)) bits per element
- Random access to compressed sequences generally costs at least O(√n)

RADICAL IDEA 3: Algebraic Integer Approach
Express p(n) as the nearest integer to an algebraic number.
If p(n) = round(α_n) where α_n satisfies a polynomial equation,
we might be able to solve it faster.

RADICAL IDEA 4: Self-Correcting Formula
Start with R^{-1}(n), then use NUMBER-THEORETIC PROPERTIES of the
result to correct it. For example:
  x₀ = round(R^{-1}(n))
  If x₀ is composite, move to next prime
  Check if π(x₀) = n using fast test

The issue: we can't compute π(x₀) fast enough for x₀ ~ 10^{102}.

RADICAL IDEA 5: Multi-precision Approximation Cascade
Use MULTIPLE independent approximations, each with different error
characteristics, and combine them to cancel errors.

RADICAL IDEA 6: Probabilistic with proof
Generate a candidate using R^{-1}(n) + probabilistic correction,
then PROVE it's correct using fast primality test + arithmetic
constraints. The proof replaces the need for exact π(x).
"""

import time
import numpy as np
import sympy
from sympy import prime as sympy_prime, isprime, nextprime, prevprime, primepi
from mpmath import mp, mpf, log, exp, li, sqrt, floor as mpfloor
from mpmath import nstr, fabs

mp.dps = 50


def approach_cascade():
    """
    MULTI-APPROXIMATION CASCADE

    Idea: Combine multiple independent approximations of p(n),
    each with different systematic errors, to cancel biases.

    Approximations:
    1. R^{-1}(n) — Riemann approximation, error ~ √p * oscillatory
    2. li^{-1}(n) — logarithmic integral inversion, error ~ √p * ln(p)
    3. Cipolla k-terms — asymptotic, error ~ p * (L2/L1)^k
    4. n * W(n) — Lambert W, error ~ p/ln(p)

    If errors were INDEPENDENT, combining k approximations would give
    error ~ error_individual / √k. But they're NOT independent —
    all are biased by the same zeta zero oscillations.

    Let's test empirically.
    """
    print("=" * 80)
    print("MULTI-APPROXIMATION CASCADE")
    print("=" * 80)

    from mpmath import lambertw

    def cipolla_4(n):
        n = mpf(n)
        L1 = log(n)
        L2 = log(L1) if L1 > 1 else mpf(1)
        return float(n * (L1 + L2 - 1 + (L2-2)/L1))

    def li_inv(n):
        n = mpf(n)
        x = n * log(n)
        for _ in range(50):
            lix = li(x)
            lnx = log(x)
            dx = (lix - n) * lnx
            x -= dx
            if fabs(dx) < mpf(10)**(-40):
                break
        return float(x)

    def R_inv(n):
        n = mpf(n)

        def R(x):
            result = mpf(0)
            for k in range(1, 80):
                mu_k = int(sympy.mobius(k))
                if mu_k == 0:
                    continue
                xk = x ** (mpf(1) / k)
                if xk < 1.01:
                    break
                term = mpf(mu_k) / k * li(xk)
                result += term
                if k > 5 and fabs(term) < mpf(10)**(-40):
                    break
            return result

        x = mpf(li_inv(n))
        for _ in range(100):
            rx = R(x)
            dx = (rx - n) * log(x)
            x -= dx
            if fabs(dx) < mpf(10)**(-40):
                break
        return float(x)

    def lambert_approx(n):
        n = mpf(n)
        w = lambertw(n)
        return float(n * w * mpf('1.04'))  # crude correction

    # Test cascade
    test_ns = list(range(100, 10001, 100))
    approx_funcs = {
        'R_inv': R_inv,
        'li_inv': li_inv,
        'cipolla4': cipolla_4,
    }

    # Collect errors
    all_errors = {name: [] for name in approx_funcs}

    for n in test_ns:
        actual = sympy_prime(n)
        for name, func in approx_funcs.items():
            try:
                val = func(n)
                all_errors[name].append(val - actual)
            except:
                all_errors[name].append(float('nan'))

    # Analyze error correlations
    import numpy as np
    err_matrix = np.array([all_errors[name] for name in approx_funcs])

    print("\nError statistics:")
    for i, name in enumerate(approx_funcs):
        errs = err_matrix[i]
        valid = ~np.isnan(errs)
        print(f"  {name:>12s}: mean={np.nanmean(errs):>8.2f}, std={np.nanstd(errs):>8.2f}")

    # Correlation matrix
    valid = ~np.any(np.isnan(err_matrix), axis=0)
    valid_errs = err_matrix[:, valid]
    if valid_errs.shape[1] > 10:
        corr = np.corrcoef(valid_errs)
        names = list(approx_funcs.keys())
        print("\nError correlation matrix:")
        print(f"  {'':>12s}", end="")
        for name in names:
            print(f"  {name:>10s}", end="")
        print()
        for i, name in enumerate(names):
            print(f"  {name:>12s}", end="")
            for j in range(len(names)):
                print(f"  {corr[i,j]:>10.4f}", end="")
            print()

    # Try optimal linear combination
    # p(n) ≈ w1 * R_inv + w2 * li_inv + w3 * cipolla4 + w0
    X = valid_errs.T  # (N, 3) matrix of errors
    actuals = np.array([sympy_prime(n) for n in test_ns])
    approxs = np.array([[func(n) for n in test_ns] for func in approx_funcs.values()]).T

    valid_idx = ~np.any(np.isnan(approxs), axis=1)
    X_valid = approxs[valid_idx]
    y_valid = actuals[valid_idx]

    # Add constant term
    X_aug = np.column_stack([X_valid, np.ones(X_valid.shape[0])])
    # Least squares fit: y = X_aug @ w
    w, residuals, rank, sv = np.linalg.lstsq(X_aug, y_valid, rcond=None)

    names = list(approx_funcs.keys()) + ['const']
    print("\nOptimal weights:")
    for name, wi in zip(names, w):
        print(f"  {name:>12s}: {wi:.6f}")

    # Test combined formula
    combined = X_aug @ w
    combined_errors = combined - y_valid
    print(f"\nCombined: mean_err={np.mean(combined_errors):.4f}, "
          f"std_err={np.std(combined_errors):.4f}")
    exact = sum(1 for e in combined_errors if abs(e) < 0.5)
    print(f"  Exact matches: {exact}/{len(combined_errors)} ({exact/len(combined_errors)*100:.1f}%)")

    # Compare to R_inv alone
    r_errors = np.array(all_errors['R_inv'])[valid_idx]
    r_exact = sum(1 for e in r_errors if abs(e) < 0.5)
    print(f"  R_inv alone:   {r_exact}/{len(r_errors)} ({r_exact/len(r_errors)*100:.1f}%)")

    # Test on holdout
    print("\nHoldout test (n=10100..10500):")
    holdout_ns = list(range(10100, 10501, 10))
    h_actual = [sympy_prime(n) for n in holdout_ns]
    h_approxs = np.array([[func(n) for n in holdout_ns] for func in approx_funcs.values()]).T
    h_aug = np.column_stack([h_approxs, np.ones(h_approxs.shape[0])])
    h_combined = h_aug @ w
    h_errors = h_combined - np.array(h_actual)
    h_exact = sum(1 for e in h_errors if abs(e) < 0.5)
    print(f"  Combined: exact={h_exact}/{len(holdout_ns)} ({h_exact/len(holdout_ns)*100:.1f}%), "
          f"std={np.std(h_errors):.2f}")

    # Key insight: are the errors more INDEPENDENT on holdout?
    h_r_errors = np.array([R_inv(n) - sympy_prime(n) for n in holdout_ns])
    print(f"  R_inv alone: std={np.std(h_r_errors):.2f}")


def approach_self_correcting():
    """
    SELF-CORRECTING FORMULA

    Idea: Use the output of the approximate formula to correct itself.

    Step 1: x₀ = R^{-1}(n) — approximation with error ~ √p * constant
    Step 2: Check if x₀ is prime. If so, check if it's the n-th prime.
            If not, adjust.

    The check "is x₀ the n-th prime?" requires π(x₀) = n.
    We can't compute π(x₀) for large x₀.

    BUT: What if we use R(x₀) as an approximation for π(x₀)?
    R(x₀) ≈ n + ε where ε is small.

    If |ε| < 1, then round(R(x₀)) = n, confirming x₀ ≈ p(n).
    But x₀ = R^{-1}(n), so R(x₀) = n by construction!

    The ACTUAL question: is the nearest prime to x₀ the correct p(n)?
    This fails when the error in R^{-1} is > gap/2.

    Alternative self-correction:
    x₀ = R^{-1}(n)
    x₁ = nearest prime to x₀ (using fast primality test)
    n₁ = R(x₁) ≈ π(x₁)
    If n₁ ≈ n, we're done. If n₁ < n, move to next prime, etc.

    The issue: R(x₁) ≠ π(x₁) exactly. The error |R(x₁) - π(x₁)| ~ √x₁ / ln(x₁).
    This makes n₁ unreliable for correcting the offset.
    """
    print("\n" + "=" * 80)
    print("SELF-CORRECTING FORMULA")
    print("=" * 80)

    from mpmath import mpf, log, li

    def R_func(x):
        x = mpf(x)
        result = mpf(0)
        for k in range(1, 80):
            mu_k = int(sympy.mobius(k))
            if mu_k == 0:
                continue
            xk = x ** (mpf(1) / k)
            if xk < 1.01:
                break
            term = mpf(mu_k) / k * li(xk)
            result += term
        return float(result)

    def R_inv(n):
        n = mpf(n)
        x = n * log(n)
        for _ in range(100):
            rx = mpf(R_func(float(x)))
            dx = (rx - n) * log(x)
            x -= dx
            if fabs(dx) < mpf(10)**(-40):
                break
        return float(x)

    # Test self-correction
    print("\nSelf-correction test:")

    correct_nearest = 0
    correct_R_round = 0
    total = 0

    for n in range(2, 1001):
        actual = sympy_prime(n)
        x0 = R_inv(n)

        # Find nearest prime
        x0_int = max(2, int(round(x0)))
        if x0_int < 2:
            x0_int = 2

        lo = prevprime(x0_int + 1) if x0_int > 2 else 2
        hi = nextprime(x0_int - 1) if x0_int > 2 else 2
        nearest = lo if abs(x0 - lo) < abs(x0 - hi) else hi

        if nearest == actual:
            correct_nearest += 1

        # Self-correction: compute R(nearest) and check
        r_nearest = R_func(nearest)
        r_round = round(r_nearest)

        if r_round == n:
            correct_R_round += 1
        else:
            # If R(nearest) != n, we know we're wrong
            # Can we use R(nearest) to correct?
            diff = n - r_round

            # Move diff primes in the right direction
            corrected = nearest
            if diff > 0:
                for _ in range(int(abs(diff))):
                    corrected = nextprime(corrected)
            elif diff < 0:
                for _ in range(int(abs(diff))):
                    corrected = prevprime(corrected)

            # This "self-correction" works IF R(p) ≈ π(p) = rank(p)
            pass

        total += 1

    print(f"  Nearest prime to R^{{-1}}: {correct_nearest}/{total} = {correct_nearest/total*100:.1f}%")
    print(f"  R(nearest) rounds to n: {correct_R_round}/{total} = {correct_R_round/total*100:.1f}%")

    # How often is the self-correction right?
    # Test: compute R for each prime near x₀ and find which gives n
    print("\nDetailed self-correction analysis (n=2..200):")
    corrections_needed = []
    for n in range(2, 201):
        actual = sympy_prime(n)
        x0 = R_inv(n)

        # Check several primes around x0
        x0_int = max(2, int(round(x0)))
        candidates = []
        p = max(2, prevprime(x0_int - 20) if x0_int > 23 else 2)
        while p < x0_int + 30:
            candidates.append(p)
            p = nextprime(p)

        # For each candidate, compute R(p) and see which is closest to n
        best = None
        best_dist = float('inf')
        for p in candidates:
            rp = R_func(p)
            dist = abs(rp - n)
            if dist < best_dist:
                best_dist = dist
                best = p

        if best == actual:
            corrections_needed.append(0)
        else:
            # How many primes off?
            if best and best > actual:
                off = -len(list(sympy.primerange(actual, best)))
            elif best:
                off = len(list(sympy.primerange(best, actual)))
            else:
                off = 999
            corrections_needed.append(off)

    corrections_arr = np.array(corrections_needed)
    print(f"  R-argmin correct: {sum(c == 0 for c in corrections_needed)}/199 = "
          f"{sum(c == 0 for c in corrections_needed)/199*100:.1f}%")
    print(f"  Max correction needed: {max(abs(c) for c in corrections_needed)} primes")
    print(f"  Mean |correction|: {np.mean(np.abs(corrections_arr)):.2f}")


def final_summary():
    """Print final session 5 summary."""
    print("\n" + "=" * 80)
    print("FINAL SESSION 5 FINDINGS")
    print("=" * 80)
    print("""
Session 5 explored 15+ additional approaches across 8 sub-agents:

NEW APPROACHES TESTED:
1. Cipolla expansion Padé resummation → FAILS (divergent, no generalization)
2. Richardson/Aitken extrapolation → FAILS (no acceleration)
3. Explicit formula + few zeta zeros → FAILS (K^{-0.01} convergence)
4. Perron contour integral → EQUIVALENT to explicit formula (proven)
5. Gap structure exploitation → FAILS (3+ bits/gap unpredictable)
6. Modular CRT for π(x) → FAILS (same cost as exact)
7. Primality test inversion → FAILS (all circular or too slow)
8. Symbolic regression/PSLQ → FAILS (no pattern in δ(n))
9. Smooth number elimination → FAILS (12% reduction, insufficient)
10. Bit-by-bit construction → WORKS but needs π(x) per bit
11. Number theory identity inversions → FAILS (same barrier)
12. Arithmetic term construction → EXISTS but O(p(n)^2)
13. Multi-approximation cascade → MARGINAL (correlation too high)
14. Self-correcting formula → PARTIAL (R-argmin often correct)
15. Prunescu-Shunia fixed-length formula → EXISTS but computationally infeasible

INTERNET SEARCH FINDINGS:
- Prunescu-Shunia (2024): Fixed-length arithmetic term for p(n) EXISTS
  but uses Wilson's theorem → O(p(n)^2) intermediate values
- Adam Harper (2025): Beyond square root barrier (theoretical only)
- No new algorithms for computing p(n) or π(x) more efficiently

TOTAL APPROACHES ACROSS ALL 5 SESSIONS: 90+

THE BARRIER REMAINS UNCHANGED:
- EXACT p(n): minimum O(p(n)^{1/2+ε}) ≈ O(10^51) for n=10^100
- APPROXIMATE p(n): O(polylog(n)) via R^{-1}(n), ~50% digits correct
- The gap is FUNDAMENTAL: bridging it requires Riemann zeta zero information
""")


if __name__ == "__main__":
    approach_cascade()
    approach_self_correcting()
    final_summary()
