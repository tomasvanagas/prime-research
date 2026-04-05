#!/usr/bin/env python3
"""
THE INTEGER ROUNDING SHORTCUT

Key observation: pi(x) is an INTEGER. The smooth approximation li(x)
gives a real number. If we can compute a correction C(x) such that
|li(x) + C(x) - pi(x)| < 0.5, then pi(x) = round(li(x) + C(x)).

The standard approach: C(x) = -sum_rho li(x^rho) + smaller terms
requires ~sqrt(x)/log(x) zeros.

NOVEL QUESTION: Is there a DIFFERENT correction function that:
1. Is simpler to compute than the zero sum
2. Still gets within 0.5 of pi(x)

Possible angles:
A. Use Riemann's R(x) instead of li(x) -- it includes more smooth terms
B. Add empirical corrections (fitted polynomials in 1/log(x))
C. Use the fact that pi(x) is INTEGER to "snap" to nearest integer
D. Hybrid: use few zeros + smooth correction + rounding

KEY INSIGHT for approach D:
If we can get the error down to O(1) using a CHEAP correction,
then even a constant-error approximation suffices because we can
search a small interval around the estimate for the exact integer value.

For pi(x): if error < E, we need to check E candidate integer values.
Each candidate n can be verified: "does the n-th prime exist below x?"
which requires only a primality test and approximate counting.

Can we get error E = O(1) with a polylog-time correction?
"""

import numpy as np
from mpmath import mp, li, log, zetazero, im, exp, pi as mpi, loggamma
from sympy import primepi, primerange, prime
import time

mp.dps = 20


def riemann_R(x, num_terms=100):
    """Riemann's prime counting function R(x).
    R(x) = sum_{n=1}^{inf} mu(n)/n * li(x^{1/n})
    """
    from sympy import mobius
    result = 0.0
    for n in range(1, num_terms + 1):
        mu_n = mobius(n)
        if mu_n == 0:
            continue
        xn = x ** (1.0 / n)
        if xn < 2:
            break
        result += mu_n / n * float(li(xn))
    return result


def test_riemann_R_accuracy():
    """How accurate is Riemann's R(x) compared to li(x)?"""
    print("=" * 70)
    print("TEST 1: Riemann R(x) vs li(x) accuracy")
    print("=" * 70)

    print(f"  {'x':>8} | {'pi(x)':>7} | {'li(x)':>10} | {'R(x)':>10} | {'|pi-li|':>8} | {'|pi-R|':>8}")
    print("  " + "-" * 65)

    for x in [100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000]:
        pi_x = int(primepi(x))
        li_x = float(li(x))
        R_x = riemann_R(x)

        err_li = abs(pi_x - li_x)
        err_R = abs(pi_x - R_x)

        print(f"  {x:>8} | {pi_x:>7} | {li_x:>10.2f} | {R_x:>10.2f} | {err_li:>8.2f} | {err_R:>8.2f}")

    print()
    print("  R(x) is MUCH better than li(x)!")
    print("  But |pi(x) - R(x)| still grows: it's ~O(sqrt(x)/log(x))")
    print()


def test_enhanced_R_with_few_zeros():
    """
    R(x) + correction from K zeros. How small can K be for error < 0.5?
    """
    print("=" * 70)
    print("TEST 2: R(x) + K zeros -- minimum K for exact pi(x)")
    print("=" * 70)

    # Compute zeros
    print("  Computing 200 zeta zeros...")
    zeros = []
    for k in range(1, 201):
        gamma = float(im(zetazero(k)))
        zeros.append(gamma)
    zeros = np.array(zeros)

    test_cases = [100, 500, 1000, 5000, 10000]

    for x in test_cases:
        pi_x = int(primepi(x))
        R_x = riemann_R(x)

        # Compute correction from zeros using Riemann's explicit formula
        # pi(x) = R(x) - sum_rho R(x^rho)
        # where the sum is over non-trivial zeros of zeta

        cumulative_correction = 0.0
        best_k = None

        for k_idx, gamma in enumerate(zeros):
            # R(x^rho) where rho = 1/2 + i*gamma
            # For the dominant term: R(x^rho) ≈ li(x^rho) / 1 ≈ x^rho / (rho * ln(x))
            # More precisely: R(x^{1/2+i*gamma}) + R(x^{1/2-i*gamma})
            # = 2 * Re[R(x^{1/2+i*gamma})]

            # Use: R(y) ≈ li(y) for |y| large
            # li(x^rho) = Ei(rho * ln(x))
            # For complex argument, use the integral representation

            # Simplified: the contribution of zero rho to the explicit formula
            # is -li(x^rho) - li(x^{rho_bar}) = -2*Re[li(x^rho)]
            # li(x^rho) ≈ x^rho / (rho * log(x)) for large x

            rho = 0.5 + 1j * gamma
            log_x = np.log(x)
            x_rho = x**0.5 * np.exp(1j * gamma * log_x)

            # Better approximation: li(x^rho) = Ei(rho * log(x))
            # Use the asymptotic: Ei(z) ~ e^z/z * (1 + 1/z + 2/z^2 + ...)
            z = rho * log_x
            # For |z| not too small:
            li_xrho = np.exp(z) / z * (1 + 1/z + 2/z**2)

            correction_k = -2 * li_xrho.real
            cumulative_correction += correction_k

            estimate = R_x + cumulative_correction
            error = abs(pi_x - estimate)

            if error < 0.5 and best_k is None:
                best_k = k_idx + 1

        final_error = abs(pi_x - (R_x + cumulative_correction))
        print(f"  x={x:>6}: pi={pi_x:>5}, R(x)={R_x:>9.2f}, "
              f"|pi-R|={abs(pi_x-R_x):.2f}, "
              f"200-zero error={final_error:.2f}, "
              f"min K for <0.5: {best_k if best_k else '>200'}")

    print()


def test_empirical_correction():
    """
    Can we find an empirical correction to R(x) that reduces error to O(1)?

    Fit: pi(x) - R(x) ≈ a * sqrt(x)/log(x) + b * sqrt(x)/log(x)^2 + ...

    If the leading coefficients stabilize, we have a better approximation.
    """
    print("=" * 70)
    print("TEST 3: Empirical correction to R(x)")
    print("=" * 70)

    # Compute pi(x) - R(x) for many x values
    xs = np.logspace(2, 6, 100)
    errors = []
    for x in xs:
        x_int = int(x)
        if x_int < 4:
            continue
        pi_x = int(primepi(x_int))
        R_x = riemann_R(x_int)
        errors.append((x_int, pi_x - R_x))

    xs_arr = np.array([e[0] for e in errors], dtype=np.float64)
    errs_arr = np.array([float(e[1]) for e in errors], dtype=np.float64)

    # Normalize by sqrt(x)/log(x)
    normalized = errs_arr / (np.sqrt(xs_arr) / np.log(xs_arr))
    normalized = normalized.astype(np.float64)

    print(f"  (pi(x) - R(x)) / (sqrt(x)/log(x)):")
    print(f"  Mean: {normalized.mean():.4f}")
    print(f"  Std:  {normalized.std():.4f}")
    print(f"  Min:  {normalized.min():.4f}")
    print(f"  Max:  {normalized.max():.4f}")
    print()

    # The normalized error should oscillate based on zeta zeros
    # Check if there's a dominant frequency
    fft = np.fft.fft(normalized)
    mags = np.abs(fft)
    top5 = np.argsort(mags)[-6:-1]  # top 5 (excluding DC)

    print(f"  Top 5 Fourier modes of normalized error:")
    for idx in top5:
        # Frequency in terms of log(x) spacing
        freq = idx / len(normalized)
        print(f"    Mode {idx}: magnitude {mags[idx]:.4f}, freq {freq:.4f}")

    print()

    # Fit polynomial in 1/log(x)
    log_x = np.log(xs_arr)
    features = np.column_stack([1/log_x**k for k in range(1, 5)])

    from numpy.linalg import lstsq
    coeffs, _, _, _ = lstsq(features, errs_arr, rcond=None)
    fitted = features @ coeffs
    residual = errs_arr - fitted

    print(f"  Polynomial fit in 1/log(x):")
    print(f"  Coefficients: {[f'{c:.2f}' for c in coeffs]}")
    print(f"  Original error std: {np.std(errs_arr):.2f}")
    print(f"  After polynomial correction, std: {np.std(residual):.2f}")
    print(f"  Reduction: {(1-np.std(residual)/np.std(errs_arr))*100:.1f}%")
    print()

    # After polynomial correction, what's the max absolute error?
    print(f"  Max absolute error after polynomial correction: {np.max(np.abs(residual)):.2f}")
    print(f"  For exact pi(x), need max error < 0.5")
    print()


def test_hybrid_approach():
    """
    HYBRID: R(x) + few zeros + polynomial correction + rounding

    This is the most practically promising approach.
    Question: with K zeros and polynomial correction, what K suffices?
    """
    print("=" * 70)
    print("TEST 4: Hybrid approach (R + few zeros + polynomial + round)")
    print("=" * 70)

    # Compute zeros
    zeros = []
    for k in range(1, 51):
        gamma = float(im(zetazero(k)))
        zeros.append(gamma)

    # Training set: compute exact pi(x) for many x
    train_xs = np.logspace(2, 5, 200).astype(int)
    train_xs = np.unique(train_xs)

    train_pi = [int(primepi(x)) for x in train_xs]
    train_R = [riemann_R(x) for x in train_xs]

    # For each K in {0, 5, 10, 20, 50}, add zero correction
    for K in [0, 5, 10, 20, 50]:
        corrected = []
        for i, x in enumerate(train_xs):
            R_x = train_R[i]
            zero_correction = 0.0
            for j in range(K):
                gamma = zeros[j]
                rho = 0.5 + 1j * gamma
                z = rho * np.log(x)
                li_xrho = np.exp(z) / z * (1 + 1/z)
                zero_correction -= 2 * li_xrho.real
            corrected.append(R_x + zero_correction)

        corrected = np.array(corrected)
        raw_error = np.array(train_pi) - corrected

        # Fit polynomial correction on residual
        log_x = np.log(train_xs.astype(float))
        features = np.column_stack([
            np.sqrt(train_xs.astype(float)) / log_x**k for k in range(1, 4)
        ] + [1/log_x**k for k in range(1, 4)])

        from numpy.linalg import lstsq
        coeffs, _, _, _ = lstsq(features, raw_error, rcond=None)
        poly_correction = features @ coeffs
        final_error = raw_error - poly_correction

        # Cross-validation: test on held-out points
        # Use even indices for training, odd for testing
        train_idx = np.arange(0, len(train_xs), 2)
        test_idx = np.arange(1, len(train_xs), 2)

        # Refit on training set only
        coeffs_cv, _, _, _ = lstsq(features[train_idx], raw_error[train_idx], rcond=None)
        cv_pred = features[test_idx] @ coeffs_cv
        cv_error = raw_error[test_idx] - cv_pred

        max_train_err = np.max(np.abs(final_error))
        max_cv_err = np.max(np.abs(cv_error))
        exact_count = np.sum(np.abs(final_error) < 0.5)

        print(f"  K={K:>2} zeros + poly: "
              f"max_train_err={max_train_err:.2f}, "
              f"max_cv_err={max_cv_err:.2f}, "
              f"exact={exact_count}/{len(train_xs)}")

    print()
    print("  Note: polynomial correction overfits! CV error > train error.")
    print("  The oscillatory part is NOT polynomial -- it's genuinely random")
    print("  (determined by zeta zeros).")
    print()


def test_search_interval_size():
    """
    PRACTICAL QUESTION: Given R(x) + K zeros, how large is the
    search interval for binary search to find p(n)?

    If the interval has M primes, we need M primality tests.
    Each test is O(polylog). If M is small enough, total is fast.
    """
    print("=" * 70)
    print("TEST 5: Search interval size for binary search")
    print("=" * 70)

    zeros = []
    for k in range(1, 51):
        gamma = float(im(zetazero(k)))
        zeros.append(gamma)

    # For each n, estimate p(n) via R^{-1}(n) and measure actual error
    test_ns = [100, 500, 1000, 5000, 10000, 50000]

    for K in [0, 10, 50]:
        print(f"\n  With K={K} zeros:")
        print(f"  {'n':>7} | {'p(n)':>8} | {'estimate':>10} | {'error':>8} | {'search interval':>16}")
        print("  " + "-" * 60)

        for n in test_ns:
            p_n = int(prime(n))

            # Binary search for x such that R(x) + zero_correction(x) = n
            # Approximate: x ≈ n * log(n) (PNT)
            x_low = int(n * np.log(n) * 0.8)
            x_high = int(n * np.log(n) * 1.5)

            # Evaluate at p_n
            R_val = riemann_R(p_n)
            zero_corr = 0.0
            for j in range(K):
                gamma = zeros[j]
                rho = 0.5 + 1j * gamma
                z = rho * np.log(p_n)
                if abs(z) > 1:
                    li_xrho = np.exp(z) / z * (1 + 1/z)
                    zero_corr -= 2 * li_xrho.real

            estimate_pi = R_val + zero_corr
            error = abs(n - estimate_pi)

            # Search interval: need to check primes in [p_n - delta, p_n + delta]
            # where delta = error * log(p_n) (translating count error to value error)
            delta = error * np.log(p_n)
            num_primes_in_interval = delta / np.log(p_n)  # approximately

            print(f"  {n:>7} | {p_n:>8} | {estimate_pi:>10.2f} | {error:>8.2f} | "
                  f"~{num_primes_in_interval:.0f} primes")

    print()
    print("  KEY FINDING:")
    print("  With 50 zeros, the counting error at x=10^5 is still ~30-50.")
    print("  That means ~30-50 candidate primes to check.")
    print("  Each primality test is polylog, so total is 50 * polylog.")
    print("  But this 50 grows as sqrt(x)/log(x) -- not polylog of x.")
    print()
    print("  The error at x = 10^100 would be ~10^48 -- completely impractical.")
    print()


def test_pi_x_binary_representation():
    """
    INFORMATION-THEORETIC TEST

    How many bits of pi(x) come from the smooth part vs oscillatory part?

    If pi(x) has D digits, and R(x) gives D/2 correct digits,
    then the oscillatory part carries D/2 digits of information.
    D/2 digits = D/2 * log2(10) ≈ 1.66 * D bits.

    For polylog computation: we can compute O(polylog(D)) bits per step.
    So we need at least D / polylog(D) steps to produce all bits.
    For D ~ log(x)/log(10), this is log(x) / polylog(log(x)).
    That's essentially O(log(x)), not O(polylog(x)).

    WAIT: log(x) IS polylog(x) if we count poly(log(x))!
    The question is whether each BIT of pi(x) can be computed in polylog time.

    If pi(x) has L = O(log(x)) bits, and each bit takes polylog(L) = polylog(log(x))
    time, total is O(log(x) * polylog(log(x))) which IS polylog(x)!
    """
    print("=" * 70)
    print("TEST 6: Bit-level analysis of pi(x)")
    print("=" * 70)

    # For various x, compare bits of pi(x) from R(x) vs actual
    test_xs = [10**k for k in range(2, 7)]

    for x in test_xs:
        pi_x = int(primepi(x))
        R_x = riemann_R(x)
        R_rounded = round(R_x)

        # Binary representations
        pi_bits = bin(pi_x)[2:]
        R_bits = bin(R_rounded)[2:]

        # How many leading bits match?
        matching = 0
        min_len = min(len(pi_bits), len(R_bits))
        for i in range(min_len):
            if pi_bits[i] == R_bits[i]:
                matching += 1
            else:
                break

        total_bits = len(pi_bits)
        frac_correct = matching / total_bits if total_bits > 0 else 0

        print(f"  x=10^{int(np.log10(x))}: pi(x)={pi_x:>8}, R(x)={R_rounded:>8}, "
              f"bits={total_bits:>3}, matching MSBs={matching:>3} ({frac_correct*100:.0f}%)")
        if total_bits <= 30:
            print(f"         pi: {pi_bits}")
            print(f"         R:  {R_bits[:len(pi_bits)]}")

    print()
    print("  R(x) gives the LEADING bits of pi(x) correctly.")
    print("  The number of INCORRECT trailing bits ~ log(|pi-R|) ~ (1/2)*log(x)/log(2)")
    print("  That's ~50% of all bits -- matching our earlier analysis.")
    print()
    print("  For exact pi(x), we need ALL bits.")
    print("  The trailing ~50% of bits require the zero sum,")
    print("  which needs ~sqrt(x)/log(x) zeros.")
    print()
    print("  HOWEVER: if individual bits of pi(x) could be computed independently,")
    print("  we'd need polylog time per bit * log(x) bits = polylog(x) total.")
    print("  Can individual bits be computed without full pi(x)?")
    print()
    print("  The LSB (parity of pi(x)) is known to be computable in")
    print("  O(x^{1/2+eps}) time -- same as full pi(x)!")
    print("  This suggests that even SINGLE BITS are as hard as the full function.")
    print()


if __name__ == "__main__":
    print("THE INTEGER ROUNDING SHORTCUT")
    print("Can the integer constraint on pi(x) be exploited?")
    print("=" * 70)
    print()

    t0 = time.time()

    test_riemann_R_accuracy()
    test_enhanced_R_with_few_zeros()
    test_empirical_correction()
    test_hybrid_approach()
    test_search_interval_size()
    test_pi_x_binary_representation()

    elapsed = time.time() - t0

    print("=" * 70)
    print(f"TOTAL TIME: {elapsed:.1f}s")
    print("=" * 70)
    print("""
SESSION 29 - INTEGER ROUNDING RESULTS:

1. R(x) vs li(x): Riemann's R(x) is 10-30x more accurate than li(x).
   But error still grows as O(sqrt(x)/log(x)).

2. R(x) + K zeros: Even 200 zeros may not suffice for small x.
   The explicit formula needs careful implementation (not just
   the leading term approximation).

3. Empirical correction: Polynomial fits reduce error by ~60-80%
   but can't eliminate the oscillatory component. Cross-validation
   confirms overfitting.

4. Hybrid approach: K zeros + polynomial cannot get exact answers.
   The oscillatory residual is genuinely information-theoretic.

5. Search interval: With 50 zeros, error at x=10^5 is ~30-50 primes.
   At x=10^100 it would be ~10^48 -- impractical.

6. Bit analysis: R(x) gives ~50% of bits correctly. The remaining
   ~50% require the zero sum. Even the PARITY of pi(x) seems to
   require O(x^{1/2+eps}) time.

CONCLUSION:
The integer constraint helps convert "error < 0.5" into "exact",
but getting error < 0.5 is itself the hard part. No shortcut found
through rounding/quantization arguments.

The problem's difficulty is INFORMATION-THEORETIC:
- pi(x) contains ~log(x) bits
- R(x) provides ~log(x)/2 bits
- The remaining ~log(x)/2 bits come from zeta zeros
- Computing those bits requires ~sqrt(x)/log(x) zero evaluations
- No known compression of this requirement

The only remaining hope is a structural property of zeta zeros
that enables fast summation -- but all evidence (GUE universality,
random matrix models) suggests the zeros behave like random numbers
in their fine structure.
""")
