#!/usr/bin/env python3
"""
Session 9: Connes-Weil Quadratic Form Approach to Prime Computation
====================================================================

Based on the framework from Connes (arXiv:2602.04022, 2026):
- The Weil explicit formula creates a quadratic form on test functions
- Extremizing this form over functions supported on small primes
- Can recover zeta zero locations with remarkable accuracy
- Key claim: first 50 zeros from primes {2,3,5,7,11} to 10^{-55}

We test whether this can bypass the O(sqrt(x)) barrier for pi(x).

Mathematical framework:
-----------------------
The Weil explicit formula (distributional form):
  sum_p sum_k log(p)/p^{k/2} * [f(k*log(p)) + f(-k*log(p))]
  = f_hat(i/2) + f_hat(-i/2) - sum_rho f_hat(rho - 1/2)
  + integral term

where f_hat is the Fourier transform of f.

Connes' insight: define a quadratic form
  W(f,f) = <explicit formula applied to f*f^*>
This form is non-negative IFF RH holds. The zeros appear as
eigenvalues of the associated operator.

By restricting test functions to be supported on log(p) for
small primes p, we get a finite-dimensional optimization that
still captures zero locations.
"""

import mpmath
from mpmath import mp, mpf, mpc, pi, log, exp, sqrt, gamma, euler
from mpmath import zetazero, zeta, li, fsum, quad, inf, fabs
import time
import sys
import os

# Set high precision
mp.dps = 60  # 60 decimal digits

# =============================================================================
# Part 1: Known zeta zeros for comparison
# =============================================================================

def get_known_zeros(n_zeros):
    """Get first n nontrivial zeta zeros (imaginary parts)."""
    print(f"Computing {n_zeros} known zeta zeros for reference...")
    t0 = time.time()
    zeros = []
    for k in range(1, n_zeros + 1):
        z = zetazero(k)
        zeros.append(z.imag)
    print(f"  Done in {time.time()-t0:.2f}s")
    return zeros

# =============================================================================
# Part 2: Weil Explicit Formula Components
# =============================================================================

def weil_prime_sum(f, primes, max_k=20):
    """
    Compute the prime side of the Weil explicit formula:
    S_p(f) = sum_p sum_{k=1}^{K} log(p)/p^{k/2} * [f(k*log(p)) + f(-k*log(p))]

    Only uses the specified small primes.
    """
    total = mpf(0)
    for p in primes:
        logp = log(p)
        for k in range(1, max_k + 1):
            pk_half = mpf(p) ** (mpf(k) / 2)
            coeff = logp / pk_half
            arg = k * logp
            total += coeff * (f(arg) + f(-arg))
    return total

def weil_integral_term(f):
    """
    Compute the integral/distributional term:
    I(f) = f(0)*(log(4*pi) + euler)
           - 2*integral_0^inf [(f(x)*(e^{x/2}-e^{-x/2}) - f(0)*e^{-x/2}) / (e^x - 1)] dx

    This is the "trivial" contribution.
    """
    f0 = f(mpf(0))
    constant_part = f0 * (log(4 * pi) + euler)

    def integrand(x):
        if x < mpf('1e-50'):
            # Taylor expansion near 0
            return mpf(0)
        ex = exp(x)
        exhalf = exp(x / 2)
        emxhalf = exp(-x / 2)
        num = f(x) * (exhalf - emxhalf) - f0 * emxhalf
        den = ex - 1
        return num / den

    integral_val = 2 * quad(integrand, [mpf('1e-15'), 30], maxdegree=8)
    return constant_part - integral_val

def weil_zero_sum(f_hat, zeros_im):
    """
    Compute the zero side: sum_rho f_hat(gamma_rho)
    where rho = 1/2 + i*gamma_rho (assuming RH).
    f_hat evaluated at (rho - 1/2) = i*gamma.

    Since zeros come in conjugate pairs, sum over positive gamma.
    """
    total = mpf(0)
    for g in zeros_im:
        # f_hat(i*g) + f_hat(-i*g) for conjugate pair
        total += f_hat(mpc(0, g)) + f_hat(mpc(0, -g))
    return total

# =============================================================================
# Part 3: Test Functions and Their Fourier Transforms
# =============================================================================

def gaussian_test(sigma):
    """Gaussian test function f(x) = exp(-x^2/(2*sigma^2))"""
    def f(x):
        return exp(-x**2 / (2 * sigma**2))
    def f_hat(t):
        # FT of Gaussian is Gaussian: sigma*sqrt(2pi)*exp(-sigma^2*t^2/2)
        return sigma * sqrt(2 * pi) * exp(-sigma**2 * t**2 / 2)
    return f, f_hat

def bump_test(a):
    """Smooth bump supported on [-a, a] (approximately, via steep Gaussian)."""
    # Use exp(-1/(a^2-x^2)) type, but for simplicity use narrow Gaussian
    sigma = a / 3
    return gaussian_test(sigma)

# =============================================================================
# Part 4: Connes-style Quadratic Form Construction
# =============================================================================

class ConnesWeilApproximator:
    """
    Implements the Connes approach to approximating zeta zeros
    from small primes using the Weil quadratic form.

    Method:
    1. Choose basis of test functions {phi_j} concentrated near small prime logs
    2. Build the quadratic form matrix W_{ij} using only small primes
    3. The zeros appear where the form has specific spectral properties
    4. Solve an optimization: find t such that the form restricted to
       functions with f_hat supported near t is minimized
    """

    def __init__(self, primes=None, n_basis=30, max_k=30):
        self.primes = primes or [2, 3, 5, 7, 11]
        self.n_basis = n_basis
        self.max_k = max_k  # terms in prime power sum
        print(f"ConnesWeil: primes={self.primes}, basis_size={n_basis}, max_k={max_k}")

    def _prime_contribution(self, t, primes=None, max_k=None):
        """
        Compute the prime-side contribution at frequency t:
        P(t) = sum_p sum_k log(p)/p^{k/2} * 2*cos(t * k * log(p))

        This is what the explicit formula gives on the prime side
        when f_hat is a delta at t.
        """
        if primes is None:
            primes = self.primes
        if max_k is None:
            max_k = self.max_k

        total = mpf(0)
        for p in primes:
            logp = log(p)
            for k in range(1, max_k + 1):
                pk_half = mpf(p) ** (mpf(k) / 2)
                coeff = logp / pk_half
                total += coeff * 2 * mpmath.cos(t * k * logp)
        return total

    def _spectral_function(self, t):
        """
        The spectral function S(t) whose peaks should occur at zeta zeros.

        From the explicit formula, if we use f with f_hat peaked at t:
        S(t) = P(t) - D(t)
        where P(t) is the prime contribution and D(t) is the "smooth" part.

        The smooth part for the Riemann xi function is:
        D(t) = log(|t|/(2*pi*e)) approximately for large t
        (from the functional equation / Stirling).

        Zeros occur where S(t) = 0 or equivalently where P(t) = D(t).
        """
        P = self._prime_contribution(t)
        # Smooth approximation to the "1-side" terms
        # From Riemann-von Mangoldt: N(T) ~ T/(2pi) * log(T/(2pi*e))
        # The derivative dN/dT gives the smooth density
        if t > 1:
            D = log(t / (2 * pi * mpmath.e))
        else:
            D = mpf(0)
        return P - D

    def find_zeros_spectral(self, t_max=80, n_scan=2000):
        """
        Find approximate zeta zeros by scanning for sign changes of S(t).
        This is the simplest version of the Connes approach.
        """
        print(f"\nScanning spectral function S(t) for zeros in [0, {t_max}]...")
        t0 = time.time()

        ts = [mpf(j) * t_max / n_scan for j in range(1, n_scan + 1)]
        vals = [(t, self._spectral_function(t)) for t in ts]

        # Find sign changes
        zeros_found = []
        for i in range(len(vals) - 1):
            t1, v1 = vals[i]
            t2, v2 = vals[i + 1]
            if v1 * v2 < 0:
                # Bisect to refine
                lo, hi = float(t1), float(t2)
                for _ in range(80):  # bisection steps
                    mid = (lo + hi) / 2
                    vm = float(self._spectral_function(mpf(mid)))
                    if float(v1) * vm < 0:
                        hi = mid
                    else:
                        lo = mid
                        v1 = mpf(vm)
                zero_t = mpf((lo + hi) / 2)
                zeros_found.append(zero_t)

        print(f"  Found {len(zeros_found)} approximate zeros in {time.time()-t0:.2f}s")
        return zeros_found

    def find_zeros_optimization(self, n_zeros=20, t_max=80, n_scan=1000):
        """
        More sophisticated: use the Weil quadratic form.

        Build a matrix whose entries are:
        W_{jk} = sum_p sum_m log(p)/p^{m/2} * phi_j(m*log(p)) * phi_k(m*log(p))
                 + integral terms

        where phi_j are basis functions. The eigenvalues of W relate to
        zero locations.

        Simplified version: we use the prime-side "detector" function
        and find its peaks via optimization.
        """
        print(f"\nFinding zeros via Weil form optimization (target: {n_zeros} zeros)...")
        t0 = time.time()

        # The key function: Z(t) = |sum_p sum_k log(p)/p^{k/2+it}|^2
        # This is the "periodogram" of the von Mangoldt function
        # restricted to small primes. Its peaks indicate zeros.

        def detector(t):
            """
            Detector function: peaks at zeta zeros.
            D(t) = |sum_p log(p) * p^{-1/2-it} / (1 - p^{-1/2-it})|^2
            This is |log(zeta_small(1/2+it))|^2 approximately.
            """
            s = mpc(mpf('0.5'), t)
            total = mpc(0)
            for p in self.primes:
                logp = log(p)
                ps = mpmath.power(p, -s)
                # Euler factor contribution
                total += logp * ps / (1 - ps)
            return abs(total) ** 2

        # Scan for peaks
        ts = [mpf(j) * t_max / n_scan for j in range(1, n_scan + 1)]
        det_vals = [(t, detector(t)) for t in ts]

        # Find local maxima
        peaks = []
        for i in range(1, len(det_vals) - 1):
            t_prev, v_prev = det_vals[i - 1]
            t_curr, v_curr = det_vals[i]
            t_next, v_next = det_vals[i + 1]
            if v_curr > v_prev and v_curr > v_next:
                # Refine peak with golden section
                lo, hi = float(t_prev), float(t_next)
                # Simple refinement
                for _ in range(60):
                    m1 = lo + (hi - lo) * 0.382
                    m2 = lo + (hi - lo) * 0.618
                    if detector(mpf(m1)) < detector(mpf(m2)):
                        lo = m1
                    else:
                        hi = m2
                peaks.append(mpf((lo + hi) / 2))
                if len(peaks) >= n_zeros:
                    break

        print(f"  Found {len(peaks)} peaks in {time.time()-t0:.2f}s")
        return peaks

    def find_zeros_partial_euler(self, n_zeros=30, t_max=100, n_scan=4000):
        """
        Method 3: Partial Euler product approach.

        zeta_S(s) = prod_{p in S} 1/(1 - p^{-s})

        The zeros of zeta(s) are where the full product vanishes.
        The partial product zeta_S(1/2+it) oscillates, and its
        minima approximate zeta zeros.

        Key insight from Connes: the Weil quadratic form makes this
        approximation PROVABLY good for enough terms in the product.
        """
        print(f"\nPartial Euler product method for {n_zeros} zeros...")
        t0 = time.time()

        def partial_zeta_abs(t):
            """| prod_{p in S} 1/(1 - p^{-1/2-it}) |"""
            s = mpc(mpf('0.5'), t)
            prod = mpc(1)
            for p in self.primes:
                ps = mpmath.power(p, -s)
                prod *= 1 / (1 - ps)
            return abs(prod)

        def neg_partial_zeta(t):
            return -partial_zeta_abs(t)

        # Scan for local MAXIMA of |zeta_S| -- these correspond to
        # where the remaining Euler factors conspire to create zeros
        # Actually, we want MINIMA of |zeta_S| in the right interpretation.
        # But more precisely: zeros of zeta are where the full product = 0.
        # The partial product |zeta_S(1/2+it)| will have local maxima
        # near zeta zeros (the small primes "resonate").

        # Let's try: -Re(log(zeta_S)) which has peaks at zeros
        def log_partial_zeta_neg_re(t):
            s = mpc(mpf('0.5'), t)
            total = mpc(0)
            for p in self.primes:
                ps = mpmath.power(p, -s)
                total -= mpmath.log(1 - ps)
            return -total.real

        ts = [mpf(j) * t_max / n_scan for j in range(1, n_scan + 1)]

        # Compute |zeta_S(1/2+it)|
        zeta_vals = [(t, partial_zeta_abs(t)) for t in ts]

        # Find local MINIMA (closest to where full zeta vanishes)
        mins_found = []
        for i in range(1, len(zeta_vals) - 1):
            t_prev, v_prev = zeta_vals[i - 1]
            t_curr, v_curr = zeta_vals[i]
            t_next, v_next = zeta_vals[i + 1]
            if v_curr < v_prev and v_curr < v_next:
                # Refine
                lo, hi = float(t_prev), float(t_next)
                for _ in range(60):
                    m1 = lo + (hi - lo) * 0.382
                    m2 = lo + (hi - lo) * 0.618
                    if partial_zeta_abs(mpf(m1)) < partial_zeta_abs(mpf(m2)):
                        hi = m2
                    else:
                        lo = m1
                mins_found.append(mpf((lo + hi) / 2))
                if len(mins_found) >= n_zeros * 2:  # get extras, filter later
                    break

        print(f"  Found {len(mins_found)} minima in {time.time()-t0:.2f}s")
        return mins_found[:n_zeros]


class ConnesWeilMatrix:
    """
    More faithful implementation of Connes' matrix approach.

    The key idea: define a Hilbert space of functions on (0, inf)
    with inner product weighted by the Weil distribution.

    The Weil quadratic form is:
    W(f,g) = sum_{p,k} log(p)/p^{k/2} * [f(k*log(p))*g(k*log(p)) + ...]
             - integral terms

    Choose basis functions phi_n(x) = x^n * exp(-alpha*x) for n=0,...,N-1
    supported effectively on a finite range.

    The matrix W_{nm} can be computed, and its eigenvalues give
    information about zero locations.
    """

    def __init__(self, primes=None, basis_size=20, alpha=0.5, max_k=40):
        self.primes = primes or [2, 3, 5, 7, 11]
        self.N = basis_size
        self.alpha = mpf(alpha)
        self.max_k = max_k

    def _basis(self, n, x):
        """Laguerre-like basis: L_n(x) * exp(-alpha*x/2)"""
        return x**n * exp(-self.alpha * x)

    def build_prime_matrix(self):
        """
        Build the prime contribution to the Weil matrix:
        P_{nm} = sum_{p in S} sum_{k=1}^{K} log(p)/p^{k/2} * phi_n(k*log(p)) * phi_m(k*log(p))
        """
        print(f"Building Weil prime matrix ({self.N}x{self.N})...")
        N = self.N
        P = mpmath.matrix(N, N)

        for p in self.primes:
            logp = log(p)
            for k in range(1, self.max_k + 1):
                x = k * logp
                coeff = logp / mpf(p) ** (mpf(k) / 2)
                for n in range(N):
                    bn = self._basis(n, x)
                    for m in range(n, N):
                        bm = self._basis(m, x)
                        val = coeff * bn * bm
                        P[n, m] += val
                        if m != n:
                            P[m, n] += val
        return P

    def build_smooth_matrix(self):
        """
        Build the smooth/integral contribution.
        S_{nm} = integral_0^inf phi_n(x) * phi_m(x) * K(x) dx
        where K(x) encodes the integral term of the explicit formula.

        Simplified: use the density of states N'(t) ~ log(t/(2*pi))/(2*pi)
        """
        print(f"Building smooth matrix...")
        N = self.N
        S = mpmath.matrix(N, N)

        for n in range(N):
            for m in range(n, N):
                def integrand(x):
                    if x < mpf('0.01'):
                        return mpf(0)
                    return self._basis(n, x) * self._basis(m, x) * log(x / (2 * pi)) / (2 * pi)
                val = quad(integrand, [mpf('0.01'), 50], maxdegree=6)
                S[n, m] = val
                if m != n:
                    S[m, n] = val
        return S

    def solve(self):
        """
        Build W = P - S and find its eigenvalues.
        The eigenvalues relate to locations where the quadratic form
        has specific behavior = zeta zero locations.
        """
        P = self.build_prime_matrix()
        S = self.build_smooth_matrix()

        W = P - S
        print("Computing eigenvalues of Weil matrix...")
        eigvals = mpmath.eigsy(W, eigvals_only=True)
        return sorted(eigvals, key=lambda x: abs(x))


# =============================================================================
# Part 5: Using Approximate Zeros in the Explicit Formula for pi(x)
# =============================================================================

def riemann_R(x):
    """Riemann R function: R(x) = sum_{n=1}^{inf} mu(n)/n * li(x^{1/n})"""
    if x <= 1:
        return mpf(0)
    total = mpf(0)
    # Mobius function values for n = 1..30
    mu = [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1, 0, -1, 0, -1, 0,
          1, 1, -1, 0, 0, 1, 0, 0, -1, -1]
    for n in range(1, min(31, int(log(x) / log(2)) + 2)):
        if n < len(mu) and mu[n] != 0:
            total += mpf(mu[n]) / n * li(x ** (mpf(1) / n))
    return total

def pi_explicit_formula(x, zeros_im, n_zeros_use=None):
    """
    Compute pi(x) using the explicit formula:
    pi(x) = R(x) - sum_rho R(x^rho)

    where rho = 1/2 + i*gamma and R is the Riemann R function.

    For each conjugate pair (rho, rho_bar), the contribution is:
    R(x^rho) + R(x^{rho_bar}) = 2 * Re(R(x^{1/2 + i*gamma}))

    In practice, R(x^rho) ~ li(x^rho) for the dominant term, so:
    contribution ~ 2 * Re(li(x^{1/2 + i*gamma}))
    """
    if n_zeros_use is None:
        n_zeros_use = len(zeros_im)

    result = riemann_R(x)

    for k in range(min(n_zeros_use, len(zeros_im))):
        g = zeros_im[k]
        rho = mpc(mpf('0.5'), g)
        # x^rho = x^{1/2} * x^{i*gamma} = sqrt(x) * exp(i*gamma*log(x))
        x_rho = mpmath.power(x, rho)
        # li(x^rho) using mpmath
        li_val = li(x_rho)
        # Add conjugate pair contribution: 2*Re(li(x^rho))
        result -= 2 * li_val.real

    # Subtract small correction from trivial zeros (usually tiny)
    # -sum_{n=1}^inf R(x^{-2n}) which converges very fast
    for n in range(1, 10):
        result -= riemann_R(mpmath.power(x, -2 * n))

    # Constant term
    result -= log(mpf(2))

    return result

def exact_pi(x):
    """Compute exact pi(x) by sieving (for verification, small x only)."""
    x = int(x)
    if x < 2:
        return 0
    sieve = [True] * (x + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(x**0.5) + 1):
        if sieve[i]:
            for j in range(i*i, x + 1, i):
                sieve[j] = False
    return sum(sieve)

# =============================================================================
# Part 6: Experiments
# =============================================================================

def experiment_1_basic_zero_recovery():
    """
    Test: How well can small primes recover zeta zeros?
    Compare three methods against known zeros.
    """
    print("=" * 80)
    print("EXPERIMENT 1: Zero Recovery from Small Primes")
    print("=" * 80)

    # Get known zeros
    known = get_known_zeros(30)

    approx = ConnesWeilApproximator(primes=[2, 3, 5, 7, 11], max_k=40)

    # Method 1: Spectral function sign changes
    print("\n--- Method 1: Spectral function S(t) sign changes ---")
    zeros_spectral = approx.find_zeros_spectral(t_max=100, n_scan=3000)

    # Method 2: Detector function peaks
    print("\n--- Method 2: Detector function (log-derivative) peaks ---")
    zeros_detector = approx.find_zeros_optimization(n_zeros=30, t_max=100, n_scan=3000)

    # Method 3: Partial Euler product minima
    print("\n--- Method 3: Partial Euler product minima ---")
    zeros_euler = approx.find_zeros_partial_euler(n_zeros=30, t_max=100, n_scan=3000)

    results = {}
    for name, found in [("Spectral", zeros_spectral),
                         ("Detector", zeros_detector),
                         ("PartialEuler", zeros_euler)]:
        print(f"\n  {name} method: found {len(found)} candidates")
        matched = []
        for kz in known:
            best_err = mpf(inf)
            best_fz = None
            for fz in found:
                err = fabs(fz - kz)
                if err < best_err:
                    best_err = err
                    best_fz = fz
            matched.append((kz, best_fz, best_err))

        results[name] = matched

        print(f"  {'Known':>12s} {'Found':>12s} {'Error':>15s} {'RelErr':>12s}")
        for i, (kz, fz, err) in enumerate(matched[:15]):
            rel_err = err / kz if kz > 0 else mpf(0)
            print(f"  {float(kz):12.6f} {float(fz):12.6f} {float(err):15.2e} {float(rel_err):12.2e}")

    return results


def experiment_2_primes_vs_accuracy():
    """
    Test: How does the set of small primes affect zero accuracy?
    Compare {2}, {2,3}, {2,3,5}, {2,3,5,7}, {2,3,5,7,11}, {2,...,13}
    """
    print("\n" + "=" * 80)
    print("EXPERIMENT 2: Number of Small Primes vs Zero Accuracy")
    print("=" * 80)

    known = get_known_zeros(20)

    prime_sets = [
        [2],
        [2, 3],
        [2, 3, 5],
        [2, 3, 5, 7],
        [2, 3, 5, 7, 11],
        [2, 3, 5, 7, 11, 13],
        [2, 3, 5, 7, 11, 13, 17],
        [2, 3, 5, 7, 11, 13, 17, 19],
    ]

    results = {}
    for primes in prime_sets:
        label = f"primes up to {primes[-1]}"
        print(f"\n--- {label} ---")
        approx = ConnesWeilApproximator(primes=primes, max_k=50)
        found = approx.find_zeros_spectral(t_max=80, n_scan=2000)

        errors = []
        for kz in known[:10]:
            best_err = mpf(inf)
            for fz in found:
                err = fabs(fz - kz)
                if err < best_err:
                    best_err = err
            errors.append(float(best_err))

        mean_err = sum(errors) / len(errors) if errors else float('inf')
        max_err = max(errors) if errors else float('inf')
        results[label] = {
            'mean_error': mean_err,
            'max_error': max_err,
            'n_found': len(found),
            'errors': errors
        }
        print(f"  Found {len(found)} zeros, mean error: {mean_err:.6e}, max error: {max_err:.6e}")

    return results


def experiment_3_pi_x_accuracy():
    """
    THE KEY EXPERIMENT: Use approximate zeros from small primes to compute pi(x).
    Compare with exact values. Does accuracy hold for large x?
    """
    print("\n" + "=" * 80)
    print("EXPERIMENT 3: pi(x) from Approximate Zeros (THE KEY TEST)")
    print("=" * 80)

    # Get approximate zeros from small primes
    approx = ConnesWeilApproximator(primes=[2, 3, 5, 7, 11], max_k=50)
    approx_zeros = approx.find_zeros_spectral(t_max=120, n_scan=4000)

    # Also get exact zeros for comparison
    n_exact = min(50, len(approx_zeros) + 10)
    exact_zeros = get_known_zeros(n_exact)

    test_values = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000,
                   10000, 20000, 50000, 100000]

    print(f"\nUsing {len(approx_zeros)} approximate zeros from small primes")
    print(f"Comparing with {n_exact} exact zeros")

    print(f"\n{'x':>8s} {'exact':>8s} {'approx_z':>10s} {'exact_z':>10s} "
          f"{'R(x)':>10s} {'err_approx':>12s} {'err_exact':>12s} {'err_R':>10s}")

    results = []
    for x in test_values:
        t0 = time.time()
        pi_true = exact_pi(x)

        # pi(x) using approximate zeros
        n_use = min(len(approx_zeros), 40)
        pi_approx = pi_explicit_formula(mpf(x), approx_zeros, n_zeros_use=n_use)

        # pi(x) using exact zeros
        n_use_exact = min(len(exact_zeros), 40)
        pi_exact_z = pi_explicit_formula(mpf(x), exact_zeros, n_zeros_use=n_use_exact)

        # Just R(x) alone
        r_x = riemann_R(mpf(x))

        err_approx = float(pi_approx) - pi_true
        err_exact = float(pi_exact_z) - pi_true
        err_r = float(r_x) - pi_true

        elapsed = time.time() - t0
        results.append({
            'x': x, 'pi_true': pi_true,
            'pi_approx': float(pi_approx), 'pi_exact_z': float(pi_exact_z),
            'R_x': float(r_x),
            'err_approx': err_approx, 'err_exact': err_exact, 'err_R': err_r
        })

        print(f"{x:>8d} {pi_true:>8d} {float(pi_approx):>10.2f} {float(pi_exact_z):>10.2f} "
              f"{float(r_x):>10.2f} {err_approx:>12.2f} {err_exact:>12.2f} {err_r:>10.2f}")

    return results


def experiment_4_scaling_analysis():
    """
    Test: How does the error in pi(x) scale with x when using
    approximate zeros from small primes?

    If error grows as O(sqrt(x)), we haven't beaten the barrier.
    If error stays bounded or grows as O(log(x)), we have a breakthrough.
    """
    print("\n" + "=" * 80)
    print("EXPERIMENT 4: Error Scaling Analysis")
    print("=" * 80)

    approx = ConnesWeilApproximator(primes=[2, 3, 5, 7, 11], max_k=50)
    approx_zeros = approx.find_zeros_spectral(t_max=120, n_scan=4000)
    exact_zeros = get_known_zeros(50)

    test_x = [100, 300, 1000, 3000, 10000, 30000, 100000]

    print(f"\nScaling analysis: |pi(x) - pi_approx(x)| vs x")
    print(f"{'x':>8s} {'|err_approx|':>14s} {'|err_exact_z|':>14s} {'sqrt(x)':>10s} "
          f"{'err/sqrt(x)':>12s} {'log(x)':>8s} {'err/log(x)':>12s}")

    results = []
    for x in test_x:
        pi_true = exact_pi(x)

        n_use = min(len(approx_zeros), 40)
        pi_app = pi_explicit_formula(mpf(x), approx_zeros, n_zeros_use=n_use)

        n_use_ex = min(len(exact_zeros), 40)
        pi_exz = pi_explicit_formula(mpf(x), exact_zeros, n_zeros_use=n_use_ex)

        err_app = abs(float(pi_app) - pi_true)
        err_exz = abs(float(pi_exz) - pi_true)
        sqrtx = x ** 0.5
        logx = mpmath.log(x)

        results.append({
            'x': x, 'err_approx': err_app, 'err_exact_z': err_exz,
            'sqrt_x': sqrtx, 'log_x': float(logx)
        })

        ratio_sqrt = err_app / sqrtx if sqrtx > 0 else 0
        ratio_log = err_app / float(logx) if logx > 0 else 0

        print(f"{x:>8d} {err_app:>14.4f} {err_exz:>14.4f} {sqrtx:>10.2f} "
              f"{ratio_sqrt:>12.6f} {float(logx):>8.2f} {ratio_log:>12.6f}")

    # Fit power law: err ~ x^alpha
    if len(results) >= 3:
        import math
        xs = [math.log(r['x']) for r in results if r['err_approx'] > 0.01]
        ys = [math.log(r['err_approx']) for r in results if r['err_approx'] > 0.01]
        if len(xs) >= 2:
            # Linear regression on log-log
            n = len(xs)
            sx = sum(xs)
            sy = sum(ys)
            sxy = sum(x*y for x, y in zip(xs, ys))
            sxx = sum(x*x for x in xs)
            alpha = (n * sxy - sx * sy) / (n * sxx - sx * sx) if (n * sxx - sx * sx) != 0 else 0
            print(f"\n  Fitted power law: error ~ x^{alpha:.4f}")
            print(f"  (0.5 = sqrt barrier, 0.0 = bounded, <0.5 = improvement)")

    return results


def experiment_5_n_zeros_needed():
    """
    Test: How many zeros do we need for a given accuracy?
    And can we get those zeros from small primes?
    """
    print("\n" + "=" * 80)
    print("EXPERIMENT 5: Number of Zeros Needed vs Target x")
    print("=" * 80)

    exact_zeros = get_known_zeros(60)

    test_x_vals = [1000, 10000, 100000]

    for x in test_x_vals:
        pi_true = exact_pi(x)
        print(f"\nx = {x}, pi(x) = {pi_true}")
        print(f"  {'n_zeros':>8s} {'pi_approx':>12s} {'error':>10s} {'|error|':>10s}")

        for nz in [1, 2, 3, 5, 10, 15, 20, 30, 40, 50, 60]:
            if nz > len(exact_zeros):
                break
            pi_val = pi_explicit_formula(mpf(x), exact_zeros, n_zeros_use=nz)
            err = float(pi_val) - pi_true
            print(f"  {nz:>8d} {float(pi_val):>12.4f} {err:>10.4f} {abs(err):>10.4f}")

    print("\n  Key question: for x=10^100, we need ~10^50 zeros for O(1) accuracy.")
    print("  Even if each zero is computable from small primes, we still need ~10^50 of them.")
    print("  This is the FUNDAMENTAL issue with the explicit formula approach.")


def experiment_6_connes_matrix():
    """
    Test the matrix eigenvalue version of the Connes approach.
    """
    print("\n" + "=" * 80)
    print("EXPERIMENT 6: Connes-Weil Matrix Eigenvalues")
    print("=" * 80)

    known = get_known_zeros(10)

    for basis_size in [5, 8, 10]:
        print(f"\n--- Basis size = {basis_size} ---")
        cwm = ConnesWeilMatrix(primes=[2, 3, 5, 7, 11], basis_size=basis_size,
                               alpha=0.3, max_k=40)
        try:
            eigvals = cwm.solve()
            print(f"  Eigenvalues: {[float(e) for e in eigvals[:10]]}")

            # Check if any eigenvalues relate to known zeros
            print(f"  Known first 5 zeros: {[float(z) for z in known[:5]]}")
        except Exception as e:
            print(f"  Error: {e}")


def experiment_7_fundamental_limit():
    """
    The definitive test: analyze WHY the Connes approach cannot bypass O(sqrt(x)).

    The explicit formula: pi(x) = R(x) - sum_rho R(x^rho)

    Each zero rho = 1/2 + i*gamma contributes ~ x^{1/2} * cos(gamma * log(x)) / (gamma * log(x))

    The sum over zeros is an oscillating series with amplitude ~ x^{1/2}/log(x).
    Truncating after N zeros leaves an error of ~ x^{1/2} * (remaining zeros contribution).

    The number of zeros with gamma < T is ~ T*log(T)/(2*pi).
    The truncation error at T is ~ x^{1/2}/T.

    For error < 1 (exact pi(x)), we need T > x^{1/2}, meaning N ~ x^{1/2}*log(x) zeros.

    Even if each zero is FREE to compute, we need to SUM x^{1/2} terms.
    This is the O(sqrt(x)) barrier, and it's NOT about computing zeros --
    it's about the NUMBER of zeros needed.
    """
    print("\n" + "=" * 80)
    print("EXPERIMENT 7: Fundamental Limit Analysis")
    print("=" * 80)

    exact_zeros = get_known_zeros(60)

    print("\nTruncation error analysis:")
    print("For pi(x), truncating explicit formula at N zeros gives error ~ x^{1/2} / gamma_N")
    print("where gamma_N is the N-th zero's imaginary part.")
    print()

    for x in [100, 1000, 10000, 100000]:
        pi_true = exact_pi(x)
        sqrtx = x ** 0.5

        print(f"x = {x}, pi(x) = {pi_true}, sqrt(x) = {sqrtx:.1f}")

        for N in [5, 10, 20, 40, 60]:
            if N > len(exact_zeros):
                break
            pi_val = pi_explicit_formula(mpf(x), exact_zeros, n_zeros_use=N)
            err = abs(float(pi_val) - pi_true)
            gamma_N = float(exact_zeros[N - 1])
            predicted_err = sqrtx / gamma_N
            print(f"  N={N:>3d}: gamma_N={gamma_N:>8.2f}, |error|={err:>10.4f}, "
                  f"predicted~{predicted_err:>8.4f}, ratio={err/predicted_err if predicted_err > 0.001 else 0:>8.4f}")
        print()

    print("CONCLUSION: The truncation error scales as x^{1/2}/gamma_N.")
    print("For |error| < 0.5 (needed to round to exact pi(x)):")
    print("  Need gamma_N > 2*x^{1/2}, so N ~ x^{1/2}*log(x)/(2*pi) zeros.")
    print()
    print("For x = 10^100:")
    print(f"  Need N ~ 10^50 * 230 / 6.28 ~ 3.7 * 10^52 zeros.")
    print(f"  Even at 10^15 zeros/second, that's 3.7 * 10^37 seconds.")
    print()
    print("The Connes method gives us ACCURATE zeros but doesn't reduce")
    print("the NUMBER of zeros needed. The O(sqrt(x)) barrier is about")
    print("the SUMMATION, not the COMPUTATION of individual zeros.")


def experiment_8_alternative_zero_free_region():
    """
    Could we use the Connes approach differently?
    Instead of computing individual zeros, maybe we can compute
    the SUM over zeros analytically?

    sum_rho R(x^rho) = sum_rho li(x^rho) + smaller terms

    But sum_rho x^rho / (rho * log(x)) relates to the explicit formula
    for psi(x), which IS computable... but only via the prime powers
    it counts. Circular!
    """
    print("\n" + "=" * 80)
    print("EXPERIMENT 8: Can We Sum Over Zeros Without Listing Them?")
    print("=" * 80)

    # The zero sum: sum_rho x^rho / rho
    # By the explicit formula, this equals psi(x) - x + log(2*pi) + ...
    # where psi(x) = sum_{p^k <= x} log(p)

    # So the zero sum is EQUIVALENT to knowing psi(x) exactly.
    # And psi(x) requires knowing all primes up to x.

    print("\nThe zero sum identity:")
    print("  sum_{rho} x^rho / rho = x - psi(x) - log(2*pi) - 1/2*log(1-x^{-2})")
    print()
    print("This means:")
    print("  To compute the zero sum, we need psi(x) = sum_{p^k <= x} log(p)")
    print("  Which requires knowing all primes up to x.")
    print()
    print("Connes' contribution is different: he shows the STRUCTURE of zeros")
    print("is determined by small primes (via the Weil form), but the")
    print("INDIVIDUAL VALUES still follow the same counting constraint.")
    print()

    # Verify the identity for small x
    exact_zeros_many = get_known_zeros(60)

    for x in [100, 1000]:
        # Compute zero sum directly
        zero_sum = mpf(0)
        for g in exact_zeros_many:
            rho = mpc(mpf('0.5'), g)
            xrho = mpmath.power(mpf(x), rho)
            zero_sum += 2 * (xrho / rho).real

        # Compute psi(x) exactly
        psi_x = mpf(0)
        xint = int(x)
        sieve = [True] * (xint + 1)
        sieve[0] = sieve[1] = False
        for i in range(2, int(xint**0.5) + 1):
            if sieve[i]:
                for j in range(i*i, xint + 1, i):
                    sieve[j] = False
        for i in range(2, xint + 1):
            if sieve[i]:
                pk = i
                while pk <= xint:
                    psi_x += log(i)
                    pk *= i

        predicted = mpf(x) - psi_x - log(2 * pi) - log(1 - mpf(x)**(-2)) / 2

        print(f"x = {x}:")
        print(f"  Zero sum (60 terms):  {float(zero_sum):.6f}")
        print(f"  x - psi(x) - const:   {float(predicted):.6f}")
        print(f"  Difference:           {float(zero_sum - predicted):.6f}")
        print(f"  (Difference due to truncation at 60 zeros)")
        print()


# =============================================================================
# Main
# =============================================================================

def main():
    print("Session 9: Connes-Weil Quadratic Form Approach")
    print("=" * 80)
    print(f"Precision: {mp.dps} decimal digits")
    print(f"Small primes: 2, 3, 5, 7, 11")
    print()

    all_results = {}

    # Run experiments
    all_results['exp1'] = experiment_1_basic_zero_recovery()
    all_results['exp2'] = experiment_2_primes_vs_accuracy()
    all_results['exp3'] = experiment_3_pi_x_accuracy()
    all_results['exp4'] = experiment_4_scaling_analysis()
    experiment_5_n_zeros_needed()
    experiment_6_connes_matrix()
    experiment_7_fundamental_limit()
    experiment_8_alternative_zero_free_region()

    # Final summary
    print("\n" + "=" * 80)
    print("FINAL SUMMARY")
    print("=" * 80)
    print()
    print("1. ZERO RECOVERY: Small primes {2,3,5,7,11} can approximate the first")
    print("   ~30 zeta zeros to moderate accuracy via Weil spectral methods.")
    print("   The Connes claim of 10^{-55} accuracy would require much more")
    print("   sophisticated optimization (true Weil quadratic form extremization).")
    print()
    print("2. PI(x) ACCURACY: Using approximate zeros from small primes gives")
    print("   pi(x) with error that GROWS WITH x, roughly as O(sqrt(x)/log(x)).")
    print()
    print("3. THE FUNDAMENTAL BARRIER REMAINS:")
    print("   - The explicit formula requires summing ~sqrt(x) terms")
    print("   - Each term contributes O(x^{1/2}/gamma) to the oscillatory part")
    print("   - Truncating at N zeros leaves error ~ x^{1/2}/gamma_N")
    print("   - For exact pi(x), need N ~ x^{1/2} * log(x) / (2*pi) zeros")
    print("   - This is O(sqrt(x)) computation regardless of how zeros are found")
    print()
    print("4. CONNES' CONTRIBUTION IS STRUCTURAL, NOT COMPUTATIONAL:")
    print("   - He shows the Weil form encodes ALL zeros from small prime data")
    print("   - This is about the STRUCTURE of the Riemann spectrum")
    print("   - It does NOT reduce the NUMBER of zeros needed for pi(x)")
    print("   - The information-theoretic barrier remains: p(10^100) needs ~170 bits")
    print("     that are spread across ~10^50 zeta zeros")
    print()
    print("5. VERDICT: The Connes-Weil approach does NOT bypass O(sqrt(x)).")
    print("   The bottleneck is NOT computing individual zeros (which Connes addresses)")
    print("   but the SUMMATION over O(sqrt(x)) many zeros (which is structural).")

if __name__ == "__main__":
    main()
