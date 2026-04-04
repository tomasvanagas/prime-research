"""
Advanced Convergence Acceleration for the Riemann Explicit Formula
===================================================================
Session 11 experiment.

KEY QUESTION: Can we reduce the number of zeta zeros needed from O(sqrt(x))
to O(polylog(x)) using advanced convergence acceleration?

The explicit formula: pi(x) = R(x) - sum_rho R(x^rho)

Truncating at T zeros gives error O(x*ln(T)/T). For exact pi(x), need T ~ x.

Previous session found: Richardson (order 1) reduces needed zeros ~10x.
This experiment tests:
  1. Error structure analysis: Is error in 1/T (smooth) or oscillatory?
  2. Higher-order Richardson (orders 1 through 10)
  3. Levin u-transform (designed for oscillatory series)
  4. Smoothed explicit formula (integrate to get faster convergence, then difference)
  5. Optimal linear combination of partial sums

The CRITICAL question: Does the error have the structure needed for acceleration?
If error = c1/T + c2/T^2 + ..., then order-k Richardson gives O(1/T^{k+1}).
If error = oscillatory in T, then Richardson fails at high orders.
"""

import mpmath
import math
import os
import sys
import time
import numpy as np

mpmath.mp.dps = 50  # high precision for analysis

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from riemann_explicit import R_function, R_at_rho, MOBIUS


def load_zeros(count):
    """Load precomputed zeta zeros."""
    data_dir = os.path.join(SCRIPT_DIR, '..', '..', 'data')
    for try_count in [1000, 500, 300, 200]:
        fname = os.path.join(data_dir, f"zeta_zeros_{try_count}.txt")
        if os.path.exists(fname):
            with open(fname) as f:
                zeros = [line.strip() for line in f if line.strip()]
            if len(zeros) >= count:
                return zeros[:count]
    raise FileNotFoundError(f"Cannot find enough zeros (need {count})")


def compute_partial_sums(x, zeros_str, max_N):
    """Compute R(x) and partial sums of zero contributions."""
    x_mpf = mpmath.mpf(x)
    Rx = float(R_function(x_mpf))

    terms = []
    for gs in zeros_str[:max_N]:
        gamma = mpmath.mpf(gs)
        rho = mpmath.mpc(0.5, gamma)
        val = 2 * float(mpmath.re(R_at_rho(x_mpf, rho)))
        terms.append(val)

    # Compute running partial sums: S(N) = R(x) - sum_{k=1}^{N} terms[k]
    partial_sums = []
    running = 0.0
    for t in terms:
        running += t
        partial_sums.append(Rx - running)

    return Rx, terms, partial_sums


# ============================================================
# PART 1: Error Structure Analysis
# ============================================================

def analyze_error_structure(exact, partial_sums):
    """
    Analyze the structure of error(N) = partial_sum(N) - exact.

    If error(N) ~ c/N: smooth, Richardson works
    If error(N) ~ A*cos(omega*N + phi)/N: oscillatory, need Levin-type
    If error(N) ~ random: no acceleration possible
    """
    errors = [ps - exact for ps in partial_sums]
    N_values = list(range(1, len(errors) + 1))

    # 1. Check if |error(N)| ~ C / N^alpha
    # Fit log|error| vs log(N)
    log_N = []
    log_err = []
    for i, e in enumerate(errors):
        if abs(e) > 1e-15:
            log_N.append(math.log(i + 1))
            log_err.append(math.log(abs(e)))

    if len(log_N) > 10:
        # Use last half for fit (more stable)
        half = len(log_N) // 2
        coeffs = np.polyfit(log_N[half:], log_err[half:], 1)
        alpha = -coeffs[0]
        C = math.exp(coeffs[1])
    else:
        alpha, C = 0, 0

    # 2. Check oscillatory structure via FFT of error sequence
    err_array = np.array(errors[10:])  # skip first few (may be noisy)
    if len(err_array) > 20:
        # Detrend by removing 1/N component
        N_arr = np.arange(11, 11 + len(err_array), dtype=float)
        detrended = err_array * N_arr  # multiply by N to remove 1/N decay

        fft = np.fft.rfft(detrended)
        power = np.abs(fft) ** 2
        freqs = np.fft.rfftfreq(len(detrended))

        # Find dominant frequency
        # Skip DC component
        if len(power) > 2:
            dom_idx = np.argmax(power[1:]) + 1
            dom_freq = freqs[dom_idx]
            dom_power = power[dom_idx]
            total_power = np.sum(power[1:])
            spectral_concentration = dom_power / total_power if total_power > 0 else 0
        else:
            dom_freq = 0
            spectral_concentration = 0
    else:
        dom_freq = 0
        spectral_concentration = 0

    # 3. Check autocorrelation of sign changes
    signs = [1 if e > 0 else -1 for e in errors[10:]]
    sign_changes = sum(1 for i in range(len(signs)-1) if signs[i] != signs[i+1])
    sign_change_rate = sign_changes / max(1, len(signs)-1)

    return {
        'decay_alpha': alpha,
        'decay_C': C,
        'dominant_freq': dom_freq,
        'spectral_concentration': spectral_concentration,
        'sign_change_rate': sign_change_rate,
        'errors': errors
    }


# ============================================================
# PART 2: Higher-Order Richardson Extrapolation
# ============================================================

def richardson_order_k(partial_sums, nodes, alpha=1):
    """
    k-th order Richardson extrapolation.

    Given partial sums S(N_0), S(N_1), ..., S(N_k) where S(N) ~ L + sum_{j=1}^{k} c_j/N^{j*alpha},
    solve the linear system to eliminate c_1, ..., c_k and estimate L.

    nodes: list of N values (indices into partial_sums, 0-based)
    """
    k = len(nodes)
    if k == 1:
        return partial_sums[nodes[0] - 1]

    # Build Vandermonde-like system
    # S(N_i) = L + c_1/N_i^alpha + c_2/N_i^{2*alpha} + ... + c_{k-1}/N_i^{(k-1)*alpha}
    # k equations, k unknowns: L, c_1, ..., c_{k-1}

    A = np.zeros((k, k))
    b = np.zeros(k)

    for i in range(k):
        N = nodes[i]
        A[i, 0] = 1.0  # coefficient of L
        for j in range(1, k):
            A[i, j] = 1.0 / N ** (j * alpha)
        b[i] = partial_sums[N - 1]  # 0-indexed

    try:
        x = np.linalg.solve(A, b)
        return x[0]  # L is the first unknown
    except np.linalg.LinAlgError:
        return partial_sums[nodes[-1] - 1]


def test_richardson_orders(exact, partial_sums, max_order=10):
    """Test Richardson extrapolation at various orders."""
    results = {}
    N_max = len(partial_sums)

    for order in range(1, min(max_order + 1, N_max // 2)):
        # Choose nodes: evenly spaced in the upper half of available zeros
        # Need order+1 nodes
        n_nodes = order + 1

        # Strategy 1: Geometric spacing from N_max/2^order to N_max
        nodes_geom = []
        for i in range(n_nodes):
            n = int(N_max * (0.5 ** (n_nodes - 1 - i)))
            n = max(5, min(n, N_max))
            nodes_geom.append(n)
        # Remove duplicates
        nodes_geom = sorted(set(nodes_geom))

        # Strategy 2: Linear spacing in upper half
        start = max(5, N_max // (order + 2))
        nodes_lin = [int(start + (N_max - start) * i / order) for i in range(n_nodes)]
        nodes_lin = sorted(set(n for n in nodes_lin if 1 <= n <= N_max))

        for name, nodes in [('geom', nodes_geom), ('linear', nodes_lin)]:
            if len(nodes) >= 2:
                try:
                    val = richardson_order_k(partial_sums, nodes, alpha=1)
                    err = val - exact
                    results[(order, name)] = {
                        'value': val,
                        'error': err,
                        'nodes': nodes,
                        'condition': np.linalg.cond(
                            build_richardson_matrix(nodes, 1)
                        ) if len(nodes) <= 15 else float('inf')
                    }
                except Exception as e:
                    results[(order, name)] = {'error': float('inf'), 'exception': str(e)}

    return results


def build_richardson_matrix(nodes, alpha):
    """Build the Richardson Vandermonde matrix for condition number analysis."""
    k = len(nodes)
    A = np.zeros((k, k))
    for i in range(k):
        A[i, 0] = 1.0
        for j in range(1, k):
            A[i, j] = 1.0 / nodes[i] ** (j * alpha)
    return A


# ============================================================
# PART 3: Levin-Type Transforms
# ============================================================

def levin_u_transform(partial_sums, terms, n_start=0):
    """
    Levin u-transform for convergence acceleration of oscillatory series.

    For a series sum a_k with partial sums S_n, the Levin u-transform uses:
    omega_n = (n+1) * a_{n+1}  (remainder estimate)

    Then: T_k^(n) = det(numerator) / det(denominator)

    Using the recursive algorithm for efficiency.
    """
    n = len(partial_sums)
    if n < 4:
        return partial_sums[-1]

    # Use the last portion of the sequence for stability
    start = max(n_start, n // 4)
    S = partial_sums[start:]
    a = terms[start:]
    m = len(S)

    if m < 4:
        return partial_sums[-1]

    # Levin u-transform: omega_n = a_n (the n-th term)
    # T_k^(n) computed via recursive algorithm

    # Using the simplified Levin algorithm (Weniger 1989)
    # Two auxiliary sequences: num[k] and den[k]
    num = [0.0] * m
    den = [0.0] * m

    for j in range(m):
        if abs(a[j]) < 1e-50:
            omega = 1e-50
        else:
            omega = a[j]  # u-transform: omega = a_n

        beta = 1.0  # for the u-transform
        num[j] = S[j] / omega
        den[j] = 1.0 / omega

    # Build the transformation table column by column
    for k in range(1, m):
        for j in range(m - k):
            factor = (j + k) / (j + k + 1)  # Levin u-factor
            diff_num = num[j+1] - num[j]
            diff_den = den[j+1] - den[j]

            if abs(diff_den) < 1e-100:
                break

            num[j] = diff_num * factor
            den[j] = diff_den * factor
        else:
            continue
        break

    if abs(den[0]) < 1e-100:
        return partial_sums[-1]

    return num[0] / den[0]


def levin_t_transform(partial_sums, terms, n_start=0):
    """
    Levin t-transform: omega_n = a_n (same as u but different recursion factor).
    More robust for some oscillatory series.
    """
    n = len(partial_sums)
    if n < 4:
        return partial_sums[-1]

    start = max(n_start, n // 4)
    S = partial_sums[start:]
    a = terms[start:]
    m = len(S)

    if m < 4:
        return partial_sums[-1]

    num = np.zeros(m)
    den = np.zeros(m)

    for j in range(m):
        omega = a[j] if abs(a[j]) > 1e-50 else 1e-50
        num[j] = S[j] / omega
        den[j] = 1.0 / omega

    for k in range(1, m):
        for j in range(m - k):
            # t-transform factor
            factor = 1.0
            diff_num = num[j+1] - num[j]
            diff_den = den[j+1] - den[j]

            if abs(diff_den) < 1e-100:
                return num[j] / den[j] if abs(den[j]) > 1e-100 else partial_sums[-1]

            num[j] = diff_num
            den[j] = diff_den

    if abs(den[0]) < 1e-100:
        return partial_sums[-1]

    return num[0] / den[0]


def weniger_delta_transform(partial_sums, terms):
    """
    Weniger's delta transformation — more powerful than Levin for some cases.

    Uses omega_n = n * a_n and a specific recursion.
    """
    n = len(partial_sums)
    if n < 6:
        return partial_sums[-1]

    # Use upper portion
    start = n // 3
    S = np.array(partial_sums[start:])
    a = np.array(terms[start:])
    m = len(S)

    if m < 6:
        return partial_sums[-1]

    # delta transform: omega_n = (n+1) * a_{n+1}
    # Remainder estimates
    omega = np.zeros(m)
    for j in range(m):
        idx = start + j
        omega[j] = (idx + 1) * a[j] if abs(a[j]) > 1e-50 else 1e-50

    # Build table using Weniger's recursive algorithm
    num = np.zeros(m)
    den = np.zeros(m)

    for j in range(m):
        num[j] = S[j] / omega[j]
        den[j] = 1.0 / omega[j]

    best_estimate = S[-1]

    for k in range(1, min(m, 20)):
        new_num = np.zeros(m - k)
        new_den = np.zeros(m - k)

        for j in range(m - k):
            # Weniger delta recursion
            c = (start + j + k) * (start + j + k - 1) / ((start + j + 2*k) * (start + j + 2*k - 1))

            d_num = num[j+1] - num[j]
            d_den = den[j+1] - den[j]

            if abs(d_den) < 1e-100:
                return best_estimate

            new_num[j] = c * d_num
            new_den[j] = c * d_den

        num = new_num
        den = new_den
        m = len(num)

        if m > 0 and abs(den[0]) > 1e-100:
            est = num[0] / den[0]
            if abs(est) < 1e15:  # sanity check
                best_estimate = est

        if m < 2:
            break

    return best_estimate


# ============================================================
# PART 4: Smoothed Explicit Formula
# ============================================================

def smoothed_pi(x, zeros_str, N, smooth_order=1):
    """
    Use the smoothed Chebyshev function for faster convergence.

    psi(x) = sum_{p^k <= x} ln(p) = x - sum_rho x^rho/rho - ln(2*pi) - (1/2)*ln(1-x^{-2})

    The integrated version psi_1(x) = integral_0^x psi(t) dt has:
    psi_1(x) = x^2/2 - sum_rho x^{rho+1}/(rho*(rho+1)) - x*ln(2*pi) - ...

    The zero sum converges as 1/gamma^2 instead of 1/gamma!

    For smooth_order=k, use psi_k with 1/gamma^{k+1} convergence.
    Then recover pi(x) by finite differences.
    """
    x_mpf = mpmath.mpf(x)
    ln_x = float(mpmath.log(x_mpf))

    if smooth_order == 0:
        # Direct: use psi(x) with zero sum converging as 1/gamma
        # psi(x) = x - sum_rho x^rho/rho - ln(2*pi) - (1/2)*ln(1-x^{-2})

        psi_smooth = float(x_mpf)  # main term

        # Subtract zero contributions
        for gs in zeros_str[:N]:
            gamma = mpmath.mpf(gs)
            rho = mpmath.mpc(0.5, gamma)
            # x^rho / rho + x^{rho_bar} / rho_bar = 2*Re(x^rho / rho)
            x_rho = mpmath.power(x_mpf, rho)
            term = 2 * float(mpmath.re(x_rho / rho))
            psi_smooth -= term

        psi_smooth -= float(mpmath.log(2 * mpmath.pi))
        if x > 1:
            psi_smooth -= 0.5 * float(mpmath.log(1 - 1/x_mpf**2))

        return psi_smooth

    elif smooth_order == 1:
        # psi_1(x) = x^2/2 - sum_rho x^{rho+1}/(rho*(rho+1)) - x*ln(2*pi) + ...
        # Converges as 1/gamma^2

        psi1 = float(x_mpf**2) / 2  # main term

        for gs in zeros_str[:N]:
            gamma = mpmath.mpf(gs)
            rho = mpmath.mpc(0.5, gamma)
            x_rho1 = mpmath.power(x_mpf, rho + 1)
            term = 2 * float(mpmath.re(x_rho1 / (rho * (rho + 1))))
            psi1 -= term

        psi1 -= float(x_mpf) * float(mpmath.log(2 * mpmath.pi))

        return psi1

    elif smooth_order == 2:
        # psi_2(x) = x^3/6 - sum_rho x^{rho+2}/(rho*(rho+1)*(rho+2)) - ...
        # Converges as 1/gamma^3

        psi2 = float(x_mpf**3) / 6

        for gs in zeros_str[:N]:
            gamma = mpmath.mpf(gs)
            rho = mpmath.mpc(0.5, gamma)
            x_rho2 = mpmath.power(x_mpf, rho + 2)
            term = 2 * float(mpmath.re(x_rho2 / (rho * (rho + 1) * (rho + 2))))
            psi2 -= term

        psi2 -= float(x_mpf**2) / 2 * float(mpmath.log(2 * mpmath.pi))

        return psi2


def recover_pi_from_smoothed(x, zeros_str, N, smooth_order=1, delta=1.0):
    """
    Recover pi(x) from the smoothed Chebyshev function via finite differences.

    psi(x) = sum_{p^k <= x} ln(p)
    pi(x) ~ psi(x) / ln(x)  (first approximation)

    More precisely: pi(x) = psi(x)/ln(x) + pi(x^{1/2})/2 + pi(x^{1/3})/3 + ...
    But the higher terms are much smaller.

    For psi_1: psi(x) ~ (psi_1(x+delta) - psi_1(x-delta)) / (2*delta)
    Then pi(x) ~ psi(x) / ln(x)
    """
    ln_x = math.log(x)

    if smooth_order == 0:
        psi_x = smoothed_pi(x, zeros_str, N, smooth_order=0)
        # pi(x) ~ psi(x)/ln(x) + integral correction
        # The correction: pi(x) = Li^{-1}(psi(x)) approximately
        # Simpler: pi(x) ~ psi(x)/ln(x) + psi(x^{1/2})/(2*ln(x))
        # For rough estimate:
        return psi_x / ln_x

    elif smooth_order == 1:
        # Use central difference
        h = delta
        psi1_plus = smoothed_pi(x + h, zeros_str, N, smooth_order=1)
        psi1_minus = smoothed_pi(x - h, zeros_str, N, smooth_order=1)
        psi_x = (psi1_plus - psi1_minus) / (2 * h)
        return psi_x / ln_x

    elif smooth_order == 2:
        # Second-order: psi(x) = (psi_2(x+h) - 2*psi_2(x) + psi_2(x-h)) / h^2
        h = delta
        psi2_plus = smoothed_pi(x + h, zeros_str, N, smooth_order=2)
        psi2_mid = smoothed_pi(x, zeros_str, N, smooth_order=2)
        psi2_minus = smoothed_pi(x - h, zeros_str, N, smooth_order=2)
        psi_x = (psi2_plus - 2 * psi2_mid + psi2_minus) / (h * h)
        return psi_x / ln_x


# ============================================================
# PART 5: Optimal Linear Combination
# ============================================================

def optimal_combination(partial_sums, exact, train_frac=0.7):
    """
    Find optimal weights for a linear combination of partial sums.

    This is the ORACLE method — it uses the exact answer.
    Purpose: determine the THEORETICAL LIMIT of linear acceleration.
    If even the optimal combination fails, no linear method can work.
    """
    n = len(partial_sums)

    # Use a subset of partial sums (subsample)
    step = max(1, n // 50)
    indices = list(range(0, n, step))
    S = np.array([partial_sums[i] for i in indices])

    # Solve: min ||S * w - exact||^2 subject to sum(w) = 1
    # This is a constrained least squares problem
    # Equivalent to: find w that minimizes (S @ w - exact)^2
    # with 1^T w = 1

    m = len(S)

    # Method: use regularized least squares
    # S * w = exact * ones
    # So each element of S * w should equal exact
    # This is just: w = pseudo_inverse(S) * [exact, exact, ..., exact]
    # But we want a single scalar output, not a vector.

    # Actually: we want scalar w^T S = exact, i.e., sum_i w_i * S_i = exact
    # with sum_i w_i = 1
    # This is 2 equations in m unknowns — underdetermined.
    # Minimize ||w||^2 (minimum norm solution).

    A = np.vstack([S, np.ones(m)])
    b = np.array([exact, 1.0])

    # Minimum norm solution: w = A^T (A A^T)^{-1} b
    try:
        AAT = A @ A.T
        w = A.T @ np.linalg.solve(AAT, b)
        val = w @ S
        return val, w, indices
    except np.linalg.LinAlgError:
        return partial_sums[-1], None, indices


# ============================================================
# MAIN EXPERIMENT
# ============================================================

EXACT_PI = {100: 25, 1000: 168, 10000: 1229, 100000: 9592, 1000000: 78498}


def run_experiment(x_val, max_zeros=1000):
    """Run all acceleration methods for a given x value."""
    exact = EXACT_PI[x_val]

    print(f"\n{'='*80}")
    print(f"  x = {x_val:,}    exact pi(x) = {exact}")
    print(f"{'='*80}")

    # Load zeros and compute terms
    zeros_str = load_zeros(max_zeros)

    print(f"  Computing {max_zeros} zero terms...", end=" ", flush=True)
    t0 = time.time()
    Rx, terms, partial_sums = compute_partial_sums(x_val, zeros_str, max_zeros)
    print(f"done in {time.time()-t0:.1f}s")
    print(f"  R(x) = {Rx:.6f}")

    # ---- PART 1: Error Structure ----
    print(f"\n  --- Error Structure Analysis ---")
    analysis = analyze_error_structure(exact, partial_sums)
    print(f"  Decay exponent alpha: {analysis['decay_alpha']:.3f} (ideal: 1.0)")
    print(f"  Dominant frequency: {analysis['dominant_freq']:.4f}")
    print(f"  Spectral concentration: {analysis['spectral_concentration']:.3f}")
    print(f"  Sign change rate: {analysis['sign_change_rate']:.3f} (0.5 = alternating, 0 = constant)")

    # Error at various N
    print(f"\n  Error at selected N values:")
    for N in [50, 100, 200, 300, 500, 700, 1000]:
        if N <= max_zeros:
            err = partial_sums[N-1] - exact
            print(f"    N={N:4d}: error = {err:+.6f}  |err| = {abs(err):.6f}")

    # ---- PART 2: Richardson Orders ----
    print(f"\n  --- Richardson Extrapolation (orders 1-10) ---")
    rich_results = test_richardson_orders(exact, partial_sums, max_order=10)

    for (order, spacing), res in sorted(rich_results.items()):
        if 'error' in res and res['error'] != float('inf'):
            err = res['error']
            cond = res.get('condition', float('inf'))
            nodes = res.get('nodes', [])
            exact_match = "EXACT" if round(res['value']) == exact else ""
            print(f"    Order {order:2d} ({spacing:6s}): error = {err:+.6f}  "
                  f"cond = {cond:.1e}  nodes={nodes}  {exact_match}")

    # ---- PART 3: Levin-Type Transforms ----
    print(f"\n  --- Levin-Type Transforms ---")

    for N in [100, 200, 300, 500, min(700, max_zeros), min(1000, max_zeros)]:
        if N > max_zeros:
            continue
        ps_sub = partial_sums[:N]
        t_sub = terms[:N]

        # Levin u
        try:
            val_u = levin_u_transform(ps_sub, t_sub)
            err_u = val_u - exact
        except:
            err_u = float('inf')

        # Levin t
        try:
            val_t = levin_t_transform(ps_sub, t_sub)
            err_t = val_t - exact
        except:
            err_t = float('inf')

        # Weniger delta
        try:
            val_w = weniger_delta_transform(ps_sub, t_sub)
            err_w = val_w - exact
        except:
            err_w = float('inf')

        plain_err = partial_sums[N-1] - exact

        print(f"    N={N:4d}: Plain={plain_err:+.4f}  "
              f"Levin-u={err_u:+.4f}  "
              f"Levin-t={err_t:+.4f}  "
              f"Weniger={err_w:+.4f}")

    # ---- PART 4: Smoothed Explicit Formula ----
    print(f"\n  --- Smoothed Explicit Formula ---")

    for N in [50, 100, 200, 500]:
        if N > max_zeros:
            continue

        # Order 0 (direct psi)
        try:
            val0 = recover_pi_from_smoothed(x_val, zeros_str, N, smooth_order=0)
            err0 = val0 - exact
        except:
            err0 = float('inf')

        # Order 1 (integrated once)
        try:
            val1 = recover_pi_from_smoothed(x_val, zeros_str, N, smooth_order=1, delta=0.5)
            err1 = val1 - exact
        except:
            err1 = float('inf')

        # Order 2 (integrated twice)
        try:
            val2 = recover_pi_from_smoothed(x_val, zeros_str, N, smooth_order=2, delta=1.0)
            err2 = val2 - exact
        except:
            err2 = float('inf')

        plain_err = partial_sums[N-1] - exact

        print(f"    N={N:3d}: Plain(R-formula)={plain_err:+.4f}  "
              f"psi/ln={err0:+.4f}  "
              f"psi_1'/{N}ln={err1:+.4f}  "
              f"psi_2''/{N}ln={err2:+.4f}")

    # ---- PART 5: Oracle Optimal ----
    print(f"\n  --- Oracle Optimal Linear Combination ---")
    opt_val, opt_w, opt_idx = optimal_combination(partial_sums, exact)
    opt_err = opt_val - exact
    print(f"    Optimal: error = {opt_err:.10f} (uses known answer)")
    if opt_w is not None:
        print(f"    Weight norm: {np.linalg.norm(opt_w):.4f}")
        print(f"    Max |weight|: {np.max(np.abs(opt_w)):.4f}")
        print(f"    Num nonzero: {np.sum(np.abs(opt_w) > 1e-10)}")

    # ---- KEY FINDING ----
    print(f"\n  --- KEY FINDING for x={x_val} ---")

    # What's the minimum N that gives exact result with each method?
    methods = {
        'Plain': lambda N: partial_sums[N-1],
        'Richardson-2': lambda N: richardson_order_k(partial_sums,
            [max(5, N//2), N], alpha=1) if N > 10 else partial_sums[N-1],
        'Richardson-3': lambda N: richardson_order_k(partial_sums,
            [max(5, N//3), max(6, 2*N//3), N], alpha=1) if N > 15 else partial_sums[N-1],
        'Levin-u': lambda N: levin_u_transform(partial_sums[:N], terms[:N]),
    }

    for mname, mfunc in methods.items():
        min_N = None
        for N in range(10, min(max_zeros + 1, 1001), 5):
            try:
                val = mfunc(N)
                if round(val) == exact:
                    min_N = N
                    break
            except:
                continue

        if min_N is not None:
            print(f"    {mname:20s}: exact at N >= {min_N}")
        else:
            print(f"    {mname:20s}: NOT exact at any N <= {min(max_zeros, 1000)}")

    return analysis


def main():
    print("Advanced Convergence Acceleration Experiment")
    print("Session 11 — Testing whether higher-order methods can beat O(sqrt(x))")
    print("=" * 80)

    all_analyses = {}

    for x_val in [1000, 10000, 100000, 1000000]:
        analysis = run_experiment(x_val, max_zeros=1000)
        all_analyses[x_val] = analysis

    # ---- SCALING ANALYSIS ----
    print("\n" + "=" * 80)
    print("SCALING ANALYSIS: How does minimum N grow with x?")
    print("=" * 80)
    print("\nIf min_N ~ x^{1/2}: no improvement over standard")
    print("If min_N ~ x^{1/3}: modest improvement")
    print("If min_N ~ log(x)^k: BREAKTHROUGH")

    print("\n" + "=" * 80)
    print("ERROR STRUCTURE: Does error have clean 1/N expansion?")
    print("=" * 80)
    for x_val, a in all_analyses.items():
        print(f"\n  x = {x_val:>10,}:")
        print(f"    Decay alpha = {a['decay_alpha']:.3f} (1.0 = clean 1/N)")
        print(f"    Sign change rate = {a['sign_change_rate']:.3f} (0.5 = oscillatory)")
        print(f"    Spectral conc. = {a['spectral_concentration']:.3f} (high = single freq)")

    print("\n" + "=" * 80)
    print("CONCLUSION")
    print("=" * 80)
    print("""
If the error structure is:
  - Smooth (alpha ~ 1, low sign changes): Higher-order Richardson CAN help
  - Oscillatory (alpha ~ 1, high sign changes): Levin-type transforms needed
  - Random (alpha ~ 0): NO acceleration possible — information barrier

The key theoretical question remains: does the truncation error of the
explicit formula admit an asymptotic expansion in 1/T, or is it
fundamentally oscillatory?
""")


if __name__ == "__main__":
    main()
