"""
Session 5: Cipolla Expansion + Padé Resummation for p(n)

The Cipolla (1902) asymptotic expansion gives:
  p(n) ~ n * (L1 + L2 - 1 + (L2-2)/L1 - ((L2^2 - 6*L2 + 11))/(2*L1^2) + ...)
where L1 = ln(n), L2 = ln(ln(n))

This is a DIVERGENT asymptotic series. We try:
1. Computing many terms of the expansion
2. Padé approximant resummation
3. Borel summation
4. Richardson extrapolation
5. Optimal truncation

Goal: Can we get error < prime gap (~ln(p(n))) for large n?
"""

from mpmath import mp, mpf, log, exp, fac, matrix, lu_solve, nstr
import sympy
from sympy import prime as sympy_prime, primepi, nextprime
import time

mp.dps = 150  # high precision

def cipolla_coefficients(max_order=30):
    """
    Compute coefficients of the Cipolla expansion to high order.

    p(n) = n * (L1 + L2 - 1 + sum_{k=1}^{inf} P_k(L2) / L1^k)

    where P_k is a polynomial in L2 of degree k.
    We compute these polynomials.

    Using the inversion of li(p) = n, where li is the logarithmic integral.
    """
    # We use the known expansion from inverting li(x) = n
    # Let w = W(n/e) (Lambert W), then p(n) ~ n*w*(1 + 1/w + ...)
    #
    # More practically, use the recursive inversion:
    # If li(x) = n, and x = n*(ln(n) + ln(ln(n)) - 1 + sum c_k / ln(n)^k * poly(ln(ln(n))))
    #
    # We compute terms iteratively using the expansion of li^{-1}

    # For now, use the known first several terms and try to extend
    # The expansion is: p(n) ≈ n*(L1 + L2 - 1 + (L2-2)/L1 + ((L2^2 - 6*L2 + 11)/2)/L1^2 + ...)

    # Actually, let's compute this numerically by fitting
    # For a range of n, compute p(n) exactly and fit the coefficients

    return None  # We'll use numerical approach below


def cipolla_approx(n, num_terms=6):
    """
    Standard Cipolla approximation with known terms.
    """
    n = mpf(n)
    L1 = log(n)
    L2 = log(L1) if L1 > 0 else mpf(0)

    if n <= 1:
        return mpf(2)

    # Term 0: n * L1
    # Term 1: + n * L2
    # Term 2: - n
    # Higher terms involve L2/L1, L2^2/L1^2, etc.

    result = n * L1
    if num_terms >= 1:
        result += n * L2
    if num_terms >= 2:
        result -= n
    if num_terms >= 3:
        result += n * (L2 - 2) / L1
    if num_terms >= 4:
        result += n * (L2**2 - 6*L2 + 11) / (2 * L1**2)
    if num_terms >= 5:
        result += n * (L2**3 - 9*L2**2 + 39*L2 - 62) / (6 * L1**3)  # was negative? check
    if num_terms >= 6:
        # 4th order term
        result += n * (L2**4 - 12*L2**3 + 78*L2**2 - 264*L2 + 371) / (24 * L1**4) # check sign

    return result


def inverse_li(n):
    """
    Compute li^{-1}(n) using Newton's method on li(x) = n.
    This is the Riemann approximation R^{-1}(n) simplified.
    """
    from mpmath import li, mpf
    n = mpf(n)

    # Initial guess from Cipolla
    x = cipolla_approx(int(n), 4)
    if x < 2:
        x = mpf(2)

    # Newton iterations: li(x) = n => x_{k+1} = x_k - (li(x_k) - n) * ln(x_k)
    for _ in range(100):
        lix = li(x)
        lnx = log(x)
        if lnx == 0:
            break
        dx = (lix - n) * lnx
        x -= dx
        if abs(dx) < mpf(10)**(-mp.dps + 10):
            break

    return x


def riemann_R_inverse(n):
    """
    Compute R^{-1}(n) where R(x) = sum_{k=1}^{inf} mu(k)/k * li(x^{1/k})
    Uses Newton's method.
    """
    from mpmath import li, mobius, mpf, sqrt, power

    n = mpf(n)

    def R(x):
        """Riemann's R function"""
        result = mpf(0)
        for k in range(1, 100):
            mu_k = int(sympy.mobius(k))
            if mu_k == 0:
                continue
            term = mpf(mu_k) / k * li(power(x, mpf(1)/k))
            result += term
            if k > 1 and abs(term) < mpf(10)**(-mp.dps + 10):
                break
        return result

    def R_prime(x):
        """Derivative of R: R'(x) = 1/ln(x)"""
        return 1 / log(x)

    # Initial guess
    x = inverse_li(n)

    # Newton iterations
    for _ in range(100):
        rx = R(x)
        dx = (rx - n) / R_prime(x)
        x -= dx
        if abs(dx) < mpf(10)**(-mp.dps + 10):
            break

    return x


def test_cipolla_accuracy():
    """Test how accurate Cipolla is at various orders."""
    print("=" * 80)
    print("CIPOLLA EXPANSION ACCURACY TEST")
    print("=" * 80)

    test_ns = [100, 1000, 10000, 100000, 1000000]

    for n in test_ns:
        actual = sympy_prime(n)
        print(f"\nn = {n}, p(n) = {actual}")

        for terms in range(1, 7):
            approx = float(cipolla_approx(n, terms))
            error = approx - actual
            rel_error = abs(error) / actual * 100
            print(f"  {terms} terms: {approx:.1f}, error = {error:+.1f} ({rel_error:.4f}%)")

        # Also test li^{-1}
        li_inv = float(inverse_li(n))
        error = li_inv - actual
        print(f"  li^-1:  {li_inv:.1f}, error = {error:+.1f} ({abs(error)/actual*100:.4f}%)")


def pade_approximant(coeffs, m, n):
    """
    Compute [m/n] Padé approximant from power series coefficients.
    Given f(x) = sum c_k x^k, find P(x)/Q(x) where deg(P)=m, deg(Q)=n.

    Returns (P_coeffs, Q_coeffs) as lists.
    """
    # Total coefficients needed: m + n + 1
    N = m + n + 1
    if len(coeffs) < N:
        coeffs = list(coeffs) + [mpf(0)] * (N - len(coeffs))

    # Solve for Q coefficients using the Padé equations
    # c_{m+1} + c_m*q_1 + ... + c_{m+1-n}*q_n = 0
    # ...
    # c_{m+n} + c_{m+n-1}*q_1 + ... + c_m*q_n = 0

    if n == 0:
        return coeffs[:m+1], [mpf(1)]

    # Build the system for q_1, ..., q_n
    A = matrix(n, n)
    b_vec = matrix(n, 1)

    for i in range(n):
        b_vec[i, 0] = -coeffs[m + 1 + i] if m + 1 + i < len(coeffs) else mpf(0)
        for j in range(n):
            idx = m - j + i
            A[i, j] = coeffs[idx] if 0 <= idx < len(coeffs) else mpf(0)

    try:
        q = lu_solve(A, b_vec)
    except:
        return coeffs[:m+1], [mpf(1)]

    Q = [mpf(1)] + [q[i, 0] for i in range(n)]

    # Compute P coefficients: p_k = c_k + sum_{j=1}^{min(k,n)} q_j * c_{k-j}
    P = []
    for k in range(m + 1):
        pk = coeffs[k]
        for j in range(1, min(k, n) + 1):
            pk += Q[j] * coeffs[k - j]
        P.append(pk)

    return P, Q


def eval_pade(P, Q, x):
    """Evaluate Padé approximant P(x)/Q(x)."""
    num = sum(p * x**k for k, p in enumerate(P))
    den = sum(q * x**k for k, q in enumerate(Q))
    if den == 0:
        return mpf('inf')
    return num / den


def test_pade_resummation():
    """
    Try Padé resummation of the Cipolla expansion.

    p(n) = n * (L1 + L2 - 1 + f(L2/L1))
    where f(t) = sum a_k * t^k (some series in t = L2/L1)

    We fit the coefficients of f numerically, then apply Padé.
    """
    print("\n" + "=" * 80)
    print("PADÉ RESUMMATION OF CIPOLLA")
    print("=" * 80)

    # Strategy: compute exact p(n) for many n, extract the "correction" series
    # Express as power series in t = L2/L1

    # For each n, compute: correction = p(n)/n - L1 - L2 + 1
    # This should be expressible as a function of t = L2/L1

    # Collect data points
    data_n = list(range(100, 10001, 10))
    corrections = []
    t_values = []

    for n in data_n:
        pn = sympy_prime(n)
        L1 = float(log(mpf(n)))
        L2 = float(log(mpf(L1)))
        t = L2 / L1
        corr = pn / n - L1 - L2 + 1
        corrections.append(corr)
        t_values.append(t)

    # Fit polynomial in t: correction ≈ sum a_k * t^k
    import numpy as np

    t_arr = np.array(t_values)
    c_arr = np.array(corrections)

    # Try different polynomial degrees
    best_deg = 0
    best_rmse = float('inf')

    for deg in range(1, 20):
        coeffs = np.polyfit(t_arr, c_arr, deg)
        pred = np.polyval(coeffs, t_arr)
        rmse = np.sqrt(np.mean((pred - c_arr)**2))
        if rmse < best_rmse:
            best_rmse = rmse
            best_deg = deg

    print(f"\nBest polynomial fit: degree {best_deg}, RMSE = {best_rmse:.6f}")

    # Now try Padé on the polynomial coefficients
    # First, get coefficients in ascending order of t
    coeffs_fit = np.polyfit(t_arr, c_arr, best_deg)[::-1]  # ascending order

    mp_coeffs = [mpf(c) for c in coeffs_fit]

    # Test different Padé orders
    print("\nPadé approximant tests:")

    # Test on larger n values
    test_ns = [50000, 100000, 500000, 1000000]

    for n in test_ns:
        actual = sympy_prime(n)
        L1 = float(log(mpf(n)))
        L2 = float(log(mpf(L1)))
        t = L2 / L1

        # Direct polynomial
        poly_corr = sum(float(mp_coeffs[k]) * t**k for k in range(len(mp_coeffs)))
        poly_pred = n * (L1 + L2 - 1 + poly_corr)

        # Various Padé approximants
        results = {}
        for m in range(1, min(8, len(mp_coeffs))):
            for k in range(1, min(8, len(mp_coeffs) - m)):
                try:
                    P, Q = pade_approximant(mp_coeffs, m, k)
                    pade_corr = float(eval_pade(P, Q, mpf(t)))
                    pade_pred = n * (L1 + L2 - 1 + pade_corr)
                    results[(m,k)] = abs(pade_pred - actual)
                except:
                    pass

        if results:
            best_pade = min(results, key=results.get)
            print(f"\nn={n}: actual={actual}")
            print(f"  Poly[{best_deg}]: error = {abs(poly_pred - actual):.1f}")
            print(f"  Best Padé [{best_pade[0]}/{best_pade[1]}]: error = {results[best_pade]:.1f}")
            print(f"  Gap ~= {float(log(mpf(actual))):.1f}")
            print(f"  Error/gap = {results[best_pade] / float(log(mpf(actual))):.3f}")


def test_richardson_extrapolation():
    """
    Apply Richardson extrapolation to the Cipolla expansion.
    If the error has the form E(h) = a*h + b*h^2 + ... where h = 1/ln(n),
    then Richardson extrapolation can eliminate terms systematically.
    """
    print("\n" + "=" * 80)
    print("RICHARDSON EXTRAPOLATION ON CIPOLLA")
    print("=" * 80)

    # For a sequence of n values, compute Cipolla at various orders
    # Then apply Richardson extrapolation

    import numpy as np

    # Test: for n=10^6, use multiple evaluations at nearby points
    target_n = 1000000
    actual = sympy_prime(target_n)

    # Compute li^{-1} for several nearby values and extrapolate
    # The idea: li^{-1}(n) has error that's a function of n
    # If we evaluate at n, 2n, 4n, ... we can extrapolate

    # Actually, let's try a different Richardson:
    # Compute p(n) estimates using different numbers of Cipolla terms
    # and extrapolate the limit

    estimates = []
    for terms in range(2, 7):
        est = float(cipolla_approx(target_n, terms))
        estimates.append(est)

    print(f"\nTarget: p({target_n}) = {actual}")
    print(f"\nCipolla estimates:")
    for i, est in enumerate(estimates):
        print(f"  {i+2} terms: {est:.2f} (error: {est - actual:+.2f})")

    # Richardson table
    R = [estimates[:]]
    for j in range(1, len(estimates)):
        new_row = []
        for i in range(len(R[-1]) - 1):
            # Standard Richardson: R_{j} = (2^j * R_{j-1}[i+1] - R_{j-1}[i]) / (2^j - 1)
            # But our "h" isn't halving. Let's try Aitken's delta-squared
            pass
        if len(new_row) > 0:
            R.append(new_row)

    # Try Aitken's delta-squared process
    if len(estimates) >= 3:
        print(f"\nAitken's Δ² acceleration:")
        for i in range(len(estimates) - 2):
            s0, s1, s2 = estimates[i], estimates[i+1], estimates[i+2]
            denom = s2 - 2*s1 + s0
            if abs(denom) > 1e-10:
                aitken = s0 - (s1 - s0)**2 / denom
                print(f"  Aitken[{i}]: {aitken:.2f} (error: {aitken - actual:+.2f})")


def test_li_inverse_correction():
    """
    Test if the correction δ(n) = p(n) - li^{-1}(n) has any exploitable pattern
    when viewed as a function of n.
    """
    print("\n" + "=" * 80)
    print("li^{-1} CORRECTION PATTERN ANALYSIS")
    print("=" * 80)

    from mpmath import li
    import numpy as np

    # Compute corrections for a range of n
    ns = list(range(10, 5001))
    deltas = []

    for n in ns:
        pn = sympy_prime(n)
        li_inv = float(inverse_li(n))
        delta = pn - li_inv
        deltas.append(delta)

    deltas = np.array(deltas)
    ns_arr = np.array(ns, dtype=float)

    print(f"\nδ(n) = p(n) - li^{{-1}}(n) statistics:")
    print(f"  Mean: {np.mean(deltas):.4f}")
    print(f"  Std:  {np.std(deltas):.4f}")
    print(f"  Min:  {np.min(deltas):.4f}")
    print(f"  Max:  {np.max(deltas):.4f}")

    # Check if δ(n) / sqrt(p(n)) has a pattern
    pn_arr = np.array([sympy_prime(n) for n in ns], dtype=float)
    normalized = deltas / np.sqrt(pn_arr)

    print(f"\nδ(n)/√p(n) statistics:")
    print(f"  Mean: {np.mean(normalized):.6f}")
    print(f"  Std:  {np.std(normalized):.6f}")

    # Autocorrelation
    d_centered = deltas - np.mean(deltas)
    autocorr = np.correlate(d_centered, d_centered, 'full')
    autocorr = autocorr[len(autocorr)//2:]
    autocorr /= autocorr[0]

    print(f"\nAutocorrelation of δ(n):")
    for lag in [1, 2, 3, 5, 10, 50, 100]:
        if lag < len(autocorr):
            print(f"  lag {lag}: {autocorr[lag]:.4f}")

    # NEW: Check if δ(n) has periodicity related to primes
    from numpy.fft import fft
    spectrum = np.abs(fft(d_centered))[:len(d_centered)//2]
    top_freqs = np.argsort(spectrum)[-10:][::-1]

    print(f"\nTop frequency components:")
    for f in top_freqs:
        period = len(d_centered) / f if f > 0 else float('inf')
        print(f"  freq={f}, period={period:.1f}, amplitude={spectrum[f]:.4f}")

    # NEW: Check if δ has structure in PRIME indices vs composite indices
    prime_idx = [i for i, n in enumerate(ns) if sympy.isprime(n)]
    comp_idx = [i for i, n in enumerate(ns) if not sympy.isprime(n) and i > 0]

    delta_at_prime_n = deltas[prime_idx]
    delta_at_comp_n = deltas[comp_idx]

    print(f"\nδ at prime n: mean={np.mean(delta_at_prime_n):.4f}, std={np.std(delta_at_prime_n):.4f}")
    print(f"δ at composite n: mean={np.mean(delta_at_comp_n):.4f}, std={np.std(delta_at_comp_n):.4f}")

    # NEW: Check if δ correlates with n mod small numbers
    for mod in [2, 3, 6, 12, 30]:
        residues = ns_arr % mod
        for r in range(mod):
            mask = residues == r
            if np.sum(mask) > 10:
                mean_d = np.mean(deltas[mask])
                # Only print if significantly different from overall mean
                if abs(mean_d - np.mean(deltas)) > 0.5 * np.std(deltas):
                    print(f"  n ≡ {r} (mod {mod}): mean δ = {mean_d:.4f}")


def test_R_inverse_high_precision():
    """
    Test R^{-1}(n) with very high precision.
    How close can we get with just R^{-1}?
    """
    print("\n" + "=" * 80)
    print("R^{-1}(n) HIGH PRECISION TEST")
    print("=" * 80)

    mp.dps = 50

    test_cases = [
        (1000, None),
        (10000, None),
        (100000, None),
        (1000000, None),
    ]

    for n, _ in test_cases:
        actual = sympy_prime(n)

        t0 = time.time()
        r_inv = riemann_R_inverse(n)
        t1 = time.time()

        error = float(r_inv - actual)
        gap = sympy_prime(n+1) - actual if n < 1000000 else nextprime(actual) - actual

        print(f"\nn={n}: p(n)={actual}")
        print(f"  R^{{-1}}(n) = {nstr(r_inv, 15)}")
        print(f"  Error: {error:+.6f}")
        print(f"  Gap to next prime: {gap}")
        print(f"  |error| / gap: {abs(error)/gap:.4f}")
        print(f"  |error| < gap/2: {abs(error) < gap/2}")
        print(f"  Time: {t1-t0:.3f}s")

    # Key question: what fraction of the time is |error| < gap/2?
    # (i.e., rounding R^{-1}(n) to nearest prime gives the correct answer)
    print("\n\nFraction where round(R^{-1}(n)) = p(n):")
    correct = 0
    total = 0
    for n in range(2, 10001):
        actual = sympy_prime(n)
        r_inv = float(inverse_li(n))  # Use simpler li^{-1} for speed

        # Find nearest prime
        if r_inv >= 2:
            lower = sympy.prevprime(int(r_inv) + 1)
            upper = nextprime(int(r_inv))
            nearest = lower if abs(r_inv - lower) < abs(r_inv - upper) else upper
            if nearest == actual:
                correct += 1
        total += 1

    print(f"  Correct: {correct}/{total} = {correct/total*100:.1f}%")
    mp.dps = 150


if __name__ == "__main__":
    test_cipolla_accuracy()
    test_pade_resummation()
    test_richardson_extrapolation()
    test_li_inverse_correction()
    test_R_inverse_high_precision()
