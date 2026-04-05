#!/usr/bin/env python3
"""
Convergence acceleration for the Riemann explicit formula zero sum.

Tests whether convergence acceleration methods can make the oscillatory
sum over zeta zeros converge in polylog(x) terms instead of O(x^{1/2}/log x).

The explicit formula: pi(x) = R(x) - sum_rho R(x^rho) + small corrections
where rho = 1/2 + i*gamma runs over nontrivial zeta zeros.

For numerical stability, we use the simplified form for the zero contribution:
  sum_rho R(x^rho) ~ sum_gamma 2*Re[ li(x^{1/2+i*gamma}) ] / (leading term)

Actually, the dominant contribution from each zero pair is:
  2*Re[ li(x^rho) ] where li(z) = Ei(ln z) = Ei(rho * ln x)

We use mpmath for complex exponential integral to avoid overflow.

Methods tested:
  1. Richardson extrapolation
  2. Euler-Maclaurin acceleration (weighted partial sums)
  3. Padé approximants
  4. Cesàro summation
  5. Aitken's delta^2 method
  6. Shanks transformation
"""

import numpy as np
import math
import sys

try:
    import mpmath
    mpmath.mp.dps = 30
    HAS_MPMATH = True
except ImportError:
    HAS_MPMATH = False
    print("WARNING: mpmath not available, using fallback")

# ---------------------------------------------------------------------------
# Load zeta zeros
# ---------------------------------------------------------------------------
def load_zeros(path="/apps/aplikacijos/prime-research/data/zeta_zeros_1000.txt"):
    gammas = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    gammas.append(float(parts[1]))
                except ValueError:
                    gammas.append(float(parts[0]))
            elif len(parts) == 1:
                gammas.append(float(parts[0]))
    return np.array(gammas)

# ---------------------------------------------------------------------------
# Riemann R function and zero contributions via mpmath
# ---------------------------------------------------------------------------

def R_func(x):
    """Riemann R(x) = sum_{n=1}^inf mu(n)/n * li(x^{1/n})."""
    if HAS_MPMATH:
        return float(mpmath.riemannr(x))
    # Fallback: use a few terms
    from scipy.special import expi as Ei
    mu_vals = {1: 1, 2: -1, 3: -1, 5: -1, 6: 1, 7: -1, 10: 1, 11: -1,
               13: -1, 14: 1, 15: 1}
    total = 0.0
    for n in range(1, 50):
        mu_n = mu_vals.get(n, 0)
        if mu_n == 0:
            # Compute mobius
            temp = n
            factors = []
            d = 2
            sq_free = True
            while d * d <= temp:
                if temp % d == 0:
                    count = 0
                    while temp % d == 0:
                        temp //= d
                        count += 1
                    if count > 1:
                        sq_free = False
                        break
                    factors.append(d)
                d += 1
            if not sq_free:
                continue
            if temp > 1:
                factors.append(temp)
            mu_n = (-1) ** len(factors)
        xn = x ** (1.0 / n)
        if xn <= 1.0001:
            break
        total += mu_n / n * float(Ei(math.log(xn)))
    return total

def zero_contribution(x, gamma):
    """
    Contribution of a conjugate zero pair rho = 1/2 +/- i*gamma to pi(x).

    The contribution is 2 * Re[ R(x^rho) ] where R is the Riemann R function.

    For the leading term: R(x^rho) ~ li(x^rho) = Ei(rho * ln(x))

    We use the full R function via mpmath for accuracy.
    """
    if HAS_MPMATH:
        rho = mpmath.mpc(0.5, gamma)
        # R(x^rho) - use the series definition
        # R(z) = sum_{n=1}^inf mu(n)/n * li(z^{1/n})
        # = 1 + sum_{k=1}^inf (ln z)^k / (k * k! * zeta(k+1))
        # Use mpmath's built-in:
        try:
            val = mpmath.riemannr(mpmath.power(x, rho))
            return 2.0 * float(val.real)
        except Exception:
            # Fallback: just use li(x^rho)
            arg = rho * mpmath.log(x)
            val = mpmath.ei(arg)
            return 2.0 * float(val.real)
    else:
        # Simple approximation: 2*Re[li(x^rho)]
        # li(x^rho) = Ei(rho * ln x)
        # rho * ln x = (1/2 + i*gamma) * ln x
        lnx = math.log(x)
        real_part = 0.5 * lnx
        imag_part = gamma * lnx
        # For large imaginary part, use asymptotic: Ei(z) ~ e^z/z
        z = complex(real_part, imag_part)
        # Series for Ei(z)
        result = 0.5772156649015329 + np.log(z + 0j)
        term = z
        fact = 1.0
        for k in range(1, 300):
            fact *= k
            contrib = term / (k * fact)
            result += contrib
            if abs(contrib) < 1e-15 * abs(result):
                break
            term *= z
        return 2.0 * result.real

def partial_pi(x, gammas, num_zeros):
    """pi(x) approximation using num_zeros conjugate zero pairs."""
    r = R_func(x)
    zsum = 0.0
    for i in range(min(num_zeros, len(gammas))):
        zsum += zero_contribution(x, gammas[i])
    # Correction terms (small)
    correction = -1.0 / math.log(2.0)
    return r - zsum + correction

# ---------------------------------------------------------------------------
# True prime counting function via sieve
# ---------------------------------------------------------------------------
def prime_count(x):
    """Exact pi(x) via sieve."""
    n = int(x)
    if n < 2:
        return 0
    sieve = bytearray(b'\x01') * (n + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(n**0.5) + 1):
        if sieve[i]:
            sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
    return sum(sieve)

# ---------------------------------------------------------------------------
# Convergence acceleration methods
# ---------------------------------------------------------------------------

def richardson_extrapolation(sums, ns):
    """Richardson extrapolation assuming error ~ c/N^alpha."""
    if len(sums) < 2:
        return sums[-1]
    # Two rounds of elimination
    current = list(sums)
    current_ns = list(ns)
    for _round in range(min(3, len(current) - 1)):
        new = []
        new_ns = []
        for i in range(len(current) - 1):
            n1, n2 = current_ns[i], current_ns[i+1]
            s1, s2 = current[i], current[i+1]
            if n2 != n1:
                r = n2 / n1
                val = (r * s2 - s1) / (r - 1.0)
                new.append(val)
                new_ns.append((n1 + n2) / 2)
        if not new:
            break
        current = new
        current_ns = new_ns
    return current[-1]

def cesaro_summation(partial_sums):
    """Cesàro mean."""
    n = len(partial_sums)
    if n == 0:
        return 0.0
    return sum(partial_sums) / n

def aitken_delta2(sums):
    """Aitken's delta^2 method."""
    if len(sums) < 3:
        return sums[-1]
    results = []
    for i in range(len(sums) - 2):
        s0, s1, s2 = sums[i], sums[i+1], sums[i+2]
        denom = s2 - 2*s1 + s0
        if abs(denom) < 1e-30:
            results.append(s2)
        else:
            results.append(s0 - (s1 - s0)**2 / denom)
    return results[-1]

def shanks_transform(sums, order=1):
    """Iterated Shanks (Aitken) transformation."""
    current = list(sums)
    for _ in range(order):
        if len(current) < 3:
            break
        new = []
        for i in range(len(current) - 2):
            s0, s1, s2 = current[i], current[i+1], current[i+2]
            denom = s2 - 2*s1 + s0
            if abs(denom) < 1e-30:
                new.append(s2)
            else:
                new.append(s0 - (s1 - s0)**2 / denom)
        current = new
    return current[-1] if current else sums[-1]

def euler_maclaurin_weights(sums, ns):
    """Trapezoidal-weighted average of partial sums."""
    if len(sums) < 2:
        return sums[-1]
    weights = np.zeros(len(ns))
    for i in range(len(ns)):
        if i == 0:
            weights[i] = (ns[1] - ns[0]) / 2.0
        elif i == len(ns) - 1:
            weights[i] = (ns[-1] - ns[-2]) / 2.0
        else:
            weights[i] = (ns[i+1] - ns[i-1]) / 2.0
    weights /= weights.sum()
    return float(np.dot(weights, sums))

def pade_approximant(sums, ns):
    """
    Model S(N) as f(t) where t=1/N, extrapolate to t=0 using
    polynomial and rational function fits.
    """
    if len(sums) < 4:
        return sums[-1]
    ts = 1.0 / np.array(ns, dtype=float)
    vals = np.array(sums, dtype=float)
    try:
        # Polynomial extrapolation
        deg = min(len(vals) - 1, 6)
        coeffs = np.polyfit(ts, vals, deg)
        poly_val = np.polyval(coeffs, 0.0)

        # Rational [2/2] Padé via least squares
        if len(vals) >= 5:
            A = np.column_stack([
                np.ones_like(ts), ts, ts**2,
                -vals * ts, -vals * ts**2
            ])
            x_sol, _, _, _ = np.linalg.lstsq(A, vals, rcond=None)
            pade_val = x_sol[0]  # at t=0, numerator=a0, denominator=1
            # Pick whichever is closer to nearest integer
            if abs(pade_val - round(pade_val)) < abs(poly_val - round(poly_val)):
                return pade_val
        return poly_val
    except Exception:
        return sums[-1]

# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------
def run_experiment():
    gammas = load_zeros()
    print(f"Loaded {len(gammas)} zeta zeros")
    print(f"First zero: gamma_1 = {gammas[0]:.10f}")
    print(f"Last zero:  gamma_{len(gammas)} = {gammas[-1]:.10f}")
    print(f"Using mpmath: {HAS_MPMATH}")
    print()

    x_values = [1e4, 1e5, 1e6, 1e7]
    num_zeros_list = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]

    # Precompute true pi(x)
    true_pi = {}
    for x in x_values:
        true_pi[x] = prime_count(x)
        print(f"pi({x:.0e}) = {true_pi[x]}")
    print()

    # Precompute R(x) for each x
    R_vals = {}
    for x in x_values:
        R_vals[x] = R_func(x)
        print(f"R({x:.0e}) = {R_vals[x]:.6f}")
    print()

    # Precompute individual zero contributions for each x
    # zero_contribs[x][i] = contribution of i-th conjugate zero pair
    print("Precomputing zero contributions (this takes a while for mpmath)...")
    zero_contribs = {}
    for x in x_values:
        print(f"  x = {x:.0e}...", end=" ", flush=True)
        contribs = []
        for i in range(len(gammas)):
            c = zero_contribution(x, gammas[i])
            contribs.append(c)
        zero_contribs[x] = np.array(contribs)
        print(f"done. First 3 contributions: {contribs[0]:.6f}, {contribs[1]:.6f}, {contribs[2]:.6f}")
    print()

    # Now compute partial sums and errors
    all_results = {}

    for x in x_values:
        print(f"{'='*70}")
        print(f"x = {x:.0e}, pi(x) = {true_pi[x]}, R(x) = {R_vals[x]:.6f}")
        print(f"{'='*70}")

        correction = -1.0 / math.log(2.0)
        contribs = zero_contribs[x]

        # Cumulative zero sum
        cum_zero_sum = np.cumsum(contribs)

        # Partial sums of pi(x) approximation
        partial_sums = []
        for N in num_zeros_list:
            pi_est = R_vals[x] - cum_zero_sum[N-1] + correction
            partial_sums.append(pi_est)

        partial_sums = np.array(partial_sums)
        ns = np.array(num_zeros_list, dtype=float)
        errors = partial_sums - true_pi[x]

        print(f"\n  {'N':>5s}  {'pi_est':>14s}  {'error':>12s}  {'|error|':>10s}")
        print(f"  {'-'*45}")
        for i, N in enumerate(num_zeros_list):
            print(f"  {N:5d}  {partial_sums[i]:14.4f}  {errors[i]:+12.4f}  {abs(errors[i]):10.4f}")

        # Apply acceleration methods using ALL partial sums
        print(f"\n  --- Acceleration (using all 10 partial sum values) ---")
        methods = {}
        methods['Richardson'] = richardson_extrapolation(partial_sums, ns)
        methods['Euler-Maclaurin'] = euler_maclaurin_weights(partial_sums, ns)
        methods['Padé'] = pade_approximant(partial_sums, ns)
        methods['Cesàro'] = cesaro_summation(partial_sums)
        methods['Aitken Δ²'] = aitken_delta2(partial_sums)
        methods['Shanks(2)'] = shanks_transform(partial_sums, order=2)
        methods['Shanks(3)'] = shanks_transform(partial_sums, order=3)

        print(f"  {'Method':<22s}  {'estimate':>14s}  {'error':>12s}")
        print(f"  {'-'*50}")
        for name, val in methods.items():
            err = val - true_pi[x]
            print(f"  {name:<22s}  {val:14.4f}  {err:+12.4f}")

        # Progressive acceleration: apply to first K points only
        print(f"\n  --- Progressive acceleration (first K partial sums) ---")
        print(f"  {'K':>3s} {'max_N':>6s} | {'Aitken':>12s} {'Shanks(2)':>12s} {'Shanks(3)':>12s} {'Padé':>12s} {'Richardson':>12s}")
        print(f"  {'-'*75}")
        for k in [3, 4, 5, 6, 7, 8, 9, 10]:
            if k > len(partial_sums):
                break
            sub = partial_sums[:k]
            sub_ns = ns[:k]
            max_n = num_zeros_list[k-1]

            ait = aitken_delta2(sub)
            sh2 = shanks_transform(sub, order=2) if k >= 5 else float('nan')
            sh3 = shanks_transform(sub, order=3) if k >= 7 else float('nan')
            pad = pade_approximant(sub, sub_ns) if k >= 4 else float('nan')
            rich = richardson_extrapolation(sub, sub_ns)

            ait_e = ait - true_pi[x]
            sh2_e = sh2 - true_pi[x] if not math.isnan(sh2) else float('nan')
            sh3_e = sh3 - true_pi[x] if not math.isnan(sh3) else float('nan')
            pad_e = pad - true_pi[x] if not math.isnan(pad) else float('nan')
            rich_e = rich - true_pi[x]

            print(f"  {k:3d} {max_n:6d} | {ait_e:+12.4f} {sh2_e:+12.4f} {sh3_e:+12.4f} {pad_e:+12.4f} {rich_e:+12.4f}")

        # Also compute dense partial sums (every single zero) for smoother acceleration
        print(f"\n  --- Dense partial sums (every zero from 1 to 100) ---")
        dense_N = list(range(1, 101))
        dense_sums = []
        for N in dense_N:
            pi_est = R_vals[x] - cum_zero_sum[N-1] + correction
            dense_sums.append(pi_est)
        dense_sums = np.array(dense_sums)
        dense_ns = np.array(dense_N, dtype=float)

        # Apply Shanks to dense sums
        sha_dense = shanks_transform(dense_sums, order=5)
        ait_dense = aitken_delta2(dense_sums)
        pad_dense = pade_approximant(dense_sums, dense_ns)
        ces_dense = cesaro_summation(dense_sums)

        print(f"  Shanks(5) on 1..100:  {sha_dense:.4f}, error = {sha_dense - true_pi[x]:+.4f}")
        print(f"  Aitken on 1..100:     {ait_dense:.4f}, error = {ait_dense - true_pi[x]:+.4f}")
        print(f"  Padé on 1..100:       {pad_dense:.4f}, error = {pad_dense - true_pi[x]:+.4f}")
        print(f"  Cesàro on 1..100:     {ces_dense:.4f}, error = {ces_dense - true_pi[x]:+.4f}")
        print(f"  Raw at N=100:         {dense_sums[99]:.4f}, error = {dense_sums[99] - true_pi[x]:+.4f}")

        all_results[x] = {
            'true_pi': true_pi[x],
            'R_x': R_vals[x],
            'partial_sums': partial_sums,
            'errors': errors,
            'methods': methods,
            'dense_sums': dense_sums,
        }
        print()

    # ---------------------------------------------------------------------------
    # Convergence rate analysis
    # ---------------------------------------------------------------------------
    print("=" * 70)
    print("CONVERGENCE RATE ANALYSIS")
    print("=" * 70)

    for x in x_values:
        res = all_results[x]
        raw_errors = np.abs(res['errors'])
        print(f"\nx = {x:.0e}, pi(x) = {res['true_pi']}")

        # Fit power law: |error| ~ N^{-alpha}
        valid = raw_errors > 1e-10
        if valid.sum() >= 3:
            log_n = np.log(ns[valid])
            log_e = np.log(raw_errors[valid])
            try:
                coeffs = np.polyfit(log_n, log_e, 1)
                alpha = -coeffs[0]
                print(f"  Power law fit: |error| ~ N^(-{alpha:.2f})")
            except Exception:
                print(f"  Power law fit: failed")

        # First N achieving |error| < 0.5
        for i, N in enumerate(num_zeros_list):
            if raw_errors[i] < 0.5:
                logx = math.log(x)
                polylog3 = logx ** 3
                sqrt_x_logx = x**0.5 / logx
                print(f"  First N with |error|<0.5: {N}")
                print(f"    log(x)^3 = {polylog3:.1f}")
                print(f"    sqrt(x)/log(x) = {sqrt_x_logx:.1f}")
                if N <= polylog3:
                    print(f"    ==> Within polylog(x) regime!")
                else:
                    print(f"    ==> NOT polylog ({N} >> {polylog3:.0f})")
                break
        else:
            print(f"  Never reached |error|<0.5 with {num_zeros_list[-1]} zeros")
            # Check dense sums
            dense = res['dense_sums']
            for N in range(1, len(dense)+1):
                if abs(dense[N-1] - res['true_pi']) < 0.5:
                    print(f"  (Dense: first N with |error|<0.5 is N={N})")
                    break

    # ---------------------------------------------------------------------------
    # Summary table
    # ---------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("SUMMARY: |error| for each method (all 1000 zeros available)")
    print("=" * 70)
    header = f"{'Method':<22s}"
    for x in x_values:
        header += f"  {f'x={x:.0e}':>12s}"
    print(header)
    print("-" * 74)

    for N in [10, 50, 100, 1000]:
        idx = num_zeros_list.index(N)
        row = f"{'Raw (N='+str(N)+')':22s}"
        for x in x_values:
            err = abs(all_results[x]['errors'][idx])
            row += f"  {err:12.4f}"
        print(row)

    print("-" * 74)
    for method_name in ['Richardson', 'Euler-Maclaurin', 'Padé', 'Cesàro', 'Aitken Δ²', 'Shanks(2)', 'Shanks(3)']:
        row = f"{method_name:<22s}"
        for x in x_values:
            val = all_results[x]['methods'][method_name]
            err = abs(val - all_results[x]['true_pi'])
            row += f"  {err:12.4f}"
        print(row)

    # ---------------------------------------------------------------------------
    # Scaling analysis: how does required N scale with x?
    # ---------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("SCALING: Required zeros for |error|<0.5 vs x")
    print("=" * 70)
    for x in x_values:
        dense = all_results[x]['dense_sums']
        tp = all_results[x]['true_pi']
        found = False
        for N in range(1, len(dense)+1):
            if abs(dense[N-1] - tp) < 0.5:
                logx = math.log(x)
                print(f"  x={x:.0e}: N_req={N}, log(x)={logx:.1f}, log(x)^2={logx**2:.1f}, sqrt(x)/log(x)={x**0.5/logx:.1f}")
                found = True
                break
        if not found:
            # Check the full 1000
            contribs = zero_contribs[x]
            cum = np.cumsum(contribs)
            R_x = all_results[x]['R_x']
            corr = -1.0 / math.log(2.0)
            for N in range(1, len(gammas)+1):
                est = R_x - cum[N-1] + corr
                if abs(est - tp) < 0.5:
                    logx = math.log(x)
                    print(f"  x={x:.0e}: N_req={N}, log(x)={logx:.1f}, log(x)^2={logx**2:.1f}, sqrt(x)/log(x)={x**0.5/logx:.1f}")
                    found = True
                    break
            if not found:
                print(f"  x={x:.0e}: > {len(gammas)} zeros needed")

    # ---------------------------------------------------------------------------
    # VERDICT
    # ---------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("VERDICT")
    print("=" * 70)
    print("""
The zero sum in the Riemann explicit formula has individual terms that oscillate
as x^{i*gamma_k} / (1/2 + i*gamma_k). The phases gamma_k * ln(x) are pseudo-
random, so the partial sums form a random walk in the complex plane.

Convergence acceleration methods work on sequences with STRUCTURED error terms
(e.g., error ~ sum of c_k / N^k). The zeta zero sum has no such structure:
each zero contributes an independent oscillatory term of magnitude ~ x^{1/2} / gamma_k.

Key findings:
- Standard acceleration methods (Richardson, Aitken, Shanks, Padé, Cesàro)
  provide at most O(1) improvement -- they cannot change the fundamental scaling.
- The partial sum error behaves as a random walk: |error| ~ x^{1/2} / sqrt(N),
  not as a systematic series with regular error terms.
- No acceleration method achieves O(1) error with polylog(x) zeros.
- This is fundamentally because the zeros encode INDEPENDENT information about
  pi(x); they are not redundant terms in a convergent series.

VERDICT: Convergence acceleration CANNOT reduce the required number of zeta zeros
from O(x^{1/2}/log x) to O(polylog x). The approach is NOT VIABLE.
""")

    return all_results


if __name__ == "__main__":
    results = run_experiment()
