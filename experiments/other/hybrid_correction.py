"""
Hybrid Correction for Riemann Explicit Formula
================================================

Explores whether the truncation error E_N(x) = pi_exact(x) - pi_N(x)
can be modeled by a smooth correction function, enabling exact pi(x)
by rounding after correction.

Approaches:
  1. Error characterization: tabulate E_N(x) for various x and N
  2. Smooth correction fitting: fit E_N(x) ~ f(x) using known pi(x) values
  3. Self-calibrating: calibrate at known points, predict at unknown ones
  4. Differential approach: compute pi(x)-pi(y) via Riemann differences
"""

import os
import sys
import time
import math

# Ensure we can import from the same directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

import mpmath
from riemann_explicit import R_function, R_at_rho, load_or_compute_zeros, pi_riemann
from nth_prime import _prime_count

mpmath.mp.dps = 30

# ============================================================================
# Utility: fast batch Riemann evaluation (reuse R(x) across different N)
# ============================================================================

def pi_riemann_incremental(x, zeros, max_n=500):
    """
    Compute pi_N(x) for all N at once, returning a dict {N: value}.
    This avoids recomputing R(x) and earlier zero contributions.
    """
    if x < 2:
        return {n: 0.0 for n in range(0, max_n+1)}

    Rx = float(R_function(x))

    cumulative_sum = 0.0
    results = {0: Rx}

    for i, gs in enumerate(zeros[:max_n]):
        gamma = mpmath.mpf(gs)
        rho = mpmath.mpc(0.5, gamma)
        contribution = 2 * float(mpmath.re(R_at_rho(x, rho)))
        cumulative_sum += contribution
        results[i + 1] = Rx - cumulative_sum

    return results


# ============================================================================
# Part 1: Error Characterization
# ============================================================================

def error_characterization():
    """
    Study E_N(x) = pi_exact(x) - pi_N(x) for various x and N.
    """
    print("=" * 90)
    print("PART 1: ERROR CHARACTERIZATION")
    print("=" * 90)

    zeros = load_or_compute_zeros(500)

    # Known exact values at powers of 10
    known_pow10 = {
        10: 4, 100: 25, 1000: 168, 10000: 1229,
        100000: 9592, 1000000: 78498, 10000000: 664579,
    }

    print("\n--- E_N(x) at powers of 10 for various N ---")
    print(f"{'x':>12s}", end="")
    Ns = [50, 100, 200, 300, 400, 500]
    for N in Ns:
        print(f"  {'E_'+str(N):>10s}", end="")
    print()
    print("-" * (12 + 12 * len(Ns)))

    for x_val, exact in sorted(known_pow10.items()):
        if x_val > 1000000:
            continue  # skip 10^7, too slow
        incremental = pi_riemann_incremental(x_val, zeros, max_n=500)
        print(f"{x_val:>12d}", end="")
        for N in Ns:
            pi_N = incremental[N]
            error = exact - pi_N
            print(f"  {error:>+10.4f}", end="")
        print()

    # Now tabulate E_300(x) for x = 10000, 20000, ..., 100000
    print("\n--- E_300(x) for x = 10000, 20000, ..., 100000 ---")
    print(f"{'x':>10s}  {'pi_exact':>8s}  {'pi_300':>12s}  {'E_300':>10s}  {'E_300/sqrt(x)*log(x)':>22s}")
    print("-" * 68)

    errors_300 = []
    for x_val in range(10000, 100001, 10000):
        exact = _prime_count(x_val)
        pi_300 = pi_riemann(x_val, num_zeros=300, zeros=zeros)
        error = exact - pi_300
        normalized = error / (math.sqrt(x_val) / math.log(x_val))
        errors_300.append((x_val, exact, pi_300, error, normalized))
        print(f"{x_val:>10d}  {exact:>8d}  {pi_300:>12.4f}  {error:>+10.4f}  {normalized:>+22.6f}")

    return errors_300


# ============================================================================
# Part 2: Smooth Correction Function
# ============================================================================

def smooth_correction():
    """
    Compute E_300(x) for x = 1000, 2000, ..., 100000 and try fitting models.
    """
    print("\n" + "=" * 90)
    print("PART 2: SMOOTH CORRECTION FUNCTION")
    print("=" * 90)

    zeros = load_or_compute_zeros(500)

    # Collect error data (use coarser grid for speed: every 5000 instead of 1000)
    xs = list(range(5000, 100001, 5000))
    errors = []

    print(f"\nComputing E_300(x) for {len(xs)} points (x = 5000, 10000, ..., 100000)...")
    t0 = time.time()

    data = []
    for i, x_val in enumerate(xs):
        exact = _prime_count(x_val)
        pi_300 = pi_riemann(x_val, num_zeros=300, zeros=zeros)
        error = exact - pi_300
        data.append((x_val, exact, pi_300, error))
        if (i + 1) % 5 == 0:
            print(f"  ...computed {i+1}/{len(xs)} ({time.time()-t0:.1f}s)")

    elapsed = time.time() - t0
    print(f"Done in {elapsed:.1f}s")

    xs_arr = [d[0] for d in data]
    errors_arr = [d[3] for d in data]

    # Print sample
    print(f"\n{'x':>8s}  {'pi_exact':>8s}  {'pi_300':>12s}  {'E_300':>10s}")
    print("-" * 44)
    for x_val, exact, pi300, err in data[::4]:
        print(f"{x_val:>8d}  {exact:>8d}  {pi300:>12.4f}  {err:>+10.4f}")

    # ---- Fit model 1: E(x) = a*sqrt(x)/log(x) + b*sqrt(x)/log(x)^2 + c ----
    print("\n--- Fitting Model 1: E(x) = a*sqrt(x)/log(x) + b*sqrt(x)/log(x)^2 + c ---")

    # Least squares via normal equations (3x3 system)
    # basis functions: f1 = sqrt(x)/log(x), f2 = sqrt(x)/log(x)^2, f3 = 1
    n = len(xs_arr)
    F = []
    for x_val in xs_arr:
        sq = math.sqrt(x_val)
        lg = math.log(x_val)
        F.append([sq / lg, sq / (lg * lg), 1.0])

    # Build A^T A and A^T b
    ATA = [[0.0]*3 for _ in range(3)]
    ATb = [0.0]*3
    for i in range(n):
        for j in range(3):
            for k in range(3):
                ATA[j][k] += F[i][j] * F[i][k]
            ATb[j] += F[i][j] * errors_arr[i]

    # Solve 3x3 system using Cramer's rule
    def det3(m):
        return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
               -m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
               +m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))

    def solve3(A, b):
        d = det3(A)
        if abs(d) < 1e-30:
            return None
        coeffs = []
        for col in range(3):
            M = [row[:] for row in A]
            for row in range(3):
                M[row][col] = b[row]
            coeffs.append(det3(M) / d)
        return coeffs

    coeffs1 = solve3(ATA, ATb)
    if coeffs1:
        a, b, c = coeffs1
        print(f"  Coefficients: a={a:.8f}, b={b:.8f}, c={c:.8f}")

        # Evaluate fit quality
        max_resid = 0
        sum_sq = 0
        for i, (x_val, exact, pi300, err) in enumerate(data):
            sq = math.sqrt(x_val)
            lg = math.log(x_val)
            fitted = a * sq / lg + b * sq / (lg * lg) + c
            resid = err - fitted
            sum_sq += resid ** 2
            if abs(resid) > abs(max_resid):
                max_resid = resid

        rmse = math.sqrt(sum_sq / n)
        print(f"  RMSE = {rmse:.6f}")
        print(f"  Max residual = {max_resid:+.6f}")
        print(f"  Can round to exact? {'YES' if abs(max_resid) < 0.5 else 'NO'} (need max |resid| < 0.5)")

    # ---- Fit model 2: E(x) = a * x^alpha ----
    print("\n--- Fitting Model 2: E(x) = a * x^alpha (log-linear fit) ---")

    # Take log: log|E| = log|a| + alpha * log(x)
    # Filter out zero/negative errors for log fit
    log_xs = []
    log_es = []
    signs = []
    for x_val, exact, pi300, err in data:
        if abs(err) > 0.01:  # skip near-zero errors
            log_xs.append(math.log(x_val))
            log_es.append(math.log(abs(err)))
            signs.append(1 if err > 0 else -1)

    if len(log_xs) > 2:
        n_log = len(log_xs)
        sum_x = sum(log_xs)
        sum_y = sum(log_es)
        sum_xy = sum(x*y for x, y in zip(log_xs, log_es))
        sum_xx = sum(x*x for x in log_xs)

        alpha = (n_log * sum_xy - sum_x * sum_y) / (n_log * sum_xx - sum_x ** 2)
        log_a = (sum_y - alpha * sum_x) / n_log
        a_pow = math.exp(log_a)

        print(f"  alpha = {alpha:.6f}")
        print(f"  a = {a_pow:.8f}")
        print(f"  Model: E(x) ~ {a_pow:.6f} * x^{alpha:.4f}")

        # Quality check
        max_resid2 = 0
        for x_val, exact, pi300, err in data:
            fitted = a_pow * (x_val ** alpha)
            resid = err - fitted
            if abs(resid) > abs(max_resid2):
                max_resid2 = resid
        print(f"  Max residual = {max_resid2:+.6f}")
        print(f"  Can round to exact? {'YES' if abs(max_resid2) < 0.5 else 'NO'}")

    # ---- Fit model 3: E(x) = a*sqrt(x)/log(x) + b*sqrt(x)/log(x)^2 + c*sqrt(x)/log(x)^3 ----
    print("\n--- Fitting Model 3: E(x) = a*sqrt(x)/log(x) + b*sqrt(x)/log(x)^2 + c*sqrt(x)/log(x)^3 ---")

    F3 = []
    for x_val in xs_arr:
        sq = math.sqrt(x_val)
        lg = math.log(x_val)
        F3.append([sq / lg, sq / (lg * lg), sq / (lg ** 3)])

    ATA3 = [[0.0]*3 for _ in range(3)]
    ATb3 = [0.0]*3
    for i in range(len(xs_arr)):
        for j in range(3):
            for k in range(3):
                ATA3[j][k] += F3[i][j] * F3[i][k]
            ATb3[j] += F3[i][j] * errors_arr[i]

    coeffs3 = solve3(ATA3, ATb3)
    if coeffs3:
        a3, b3, c3 = coeffs3
        print(f"  Coefficients: a={a3:.8f}, b={b3:.8f}, c={c3:.8f}")

        max_resid3 = 0
        sum_sq3 = 0
        for i, (x_val, exact, pi300, err) in enumerate(data):
            sq = math.sqrt(x_val)
            lg = math.log(x_val)
            fitted = a3 * sq / lg + b3 * sq / (lg ** 2) + c3 * sq / (lg ** 3)
            resid = err - fitted
            sum_sq3 += resid ** 2
            if abs(resid) > abs(max_resid3):
                max_resid3 = resid

        rmse3 = math.sqrt(sum_sq3 / len(xs_arr))
        print(f"  RMSE = {rmse3:.6f}")
        print(f"  Max residual = {max_resid3:+.6f}")
        print(f"  Can round to exact? {'YES' if abs(max_resid3) < 0.5 else 'NO'}")

    return data, coeffs1, coeffs3


# ============================================================================
# Part 3: Self-Calibrating Approach
# ============================================================================

def self_calibrating():
    """
    Calibrate at known points, predict at unknown ones.
    """
    print("\n" + "=" * 90)
    print("PART 3: SELF-CALIBRATING APPROACH")
    print("=" * 90)

    zeros = load_or_compute_zeros(500)

    # Calibration points — compute exact pi(x) and Riemann pi(x)
    cal_points = [5000, 10000, 20000, 50000, 75000, 100000]
    print(f"\nCalibrating at: {cal_points}")

    cal_data = []
    for x_val in cal_points:
        exact = _prime_count(x_val)
        pi_500 = pi_riemann(x_val, num_zeros=500, zeros=zeros)
        error = exact - pi_500
        cal_data.append((x_val, exact, pi_500, error))
        print(f"  x={x_val:>8d}  pi_exact={exact:>6d}  pi_500={pi_500:>12.4f}  E_500={error:>+8.4f}")

    # Fit correction: E(x) = a*sqrt(x)/log(x) + b*sqrt(x)/log(x)^2 + c
    def basis(x):
        sq = math.sqrt(x)
        lg = math.log(x)
        return [sq / lg, sq / (lg * lg), 1.0]

    def det3(m):
        return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
               -m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
               +m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))

    def solve3(A, b):
        d = det3(A)
        if abs(d) < 1e-30:
            return None
        coeffs = []
        for col in range(3):
            M = [row[:] for row in A]
            for row in range(3):
                M[row][col] = b[row]
            coeffs.append(det3(M) / d)
        return coeffs

    n = len(cal_data)
    ATA = [[0.0]*3 for _ in range(3)]
    ATb = [0.0]*3
    for x_val, exact, pi_est, error in cal_data:
        f = basis(x_val)
        for j in range(3):
            for k in range(3):
                ATA[j][k] += f[j] * f[k]
            ATb[j] += f[j] * error

    coeffs = solve3(ATA, ATb)
    if not coeffs:
        print("  Failed to fit correction!")
        return

    a, b, c = coeffs
    print(f"\n  Fitted correction: E(x) = {a:.6f}*sqrt(x)/log(x) + {b:.6f}*sqrt(x)/log(x)^2 + {c:.6f}")

    # Test calibration residuals
    print(f"\n  Calibration residuals:")
    for x_val, exact, pi_est, error in cal_data:
        f = basis(x_val)
        predicted_error = sum(ci * fi for ci, fi in zip(coeffs, f))
        corrected = pi_est + predicted_error
        resid = exact - corrected
        print(f"    x={x_val:>8d}  corrected={corrected:>10.4f}  round={round(corrected):>6d}  "
              f"exact={exact:>6d}  resid={resid:>+8.4f}  {'OK' if round(corrected)==exact else 'MISS'}")

    # Predict at test points
    test_points = [150000, 200000, 500000, 1000000]
    print(f"\n  Predictions at test points (using 500 zeros + correction):")
    print(f"  {'x':>10s}  {'pi_500':>12s}  {'correction':>12s}  {'corrected':>12s}  {'rounded':>8s}  {'exact':>8s}  {'ok?':>5s}")
    print("  " + "-" * 75)

    for x_val in test_points:
        t0 = time.time()
        pi_est = pi_riemann(x_val, num_zeros=500, zeros=zeros)
        f = basis(x_val)
        correction = sum(ci * fi for ci, fi in zip(coeffs, f))
        corrected = pi_est + correction
        rounded = round(corrected)
        exact = _prime_count(x_val)
        elapsed = time.time() - t0
        ok = "YES" if rounded == exact else "NO"
        print(f"  {x_val:>10d}  {pi_est:>12.4f}  {correction:>+12.4f}  {corrected:>12.4f}  "
              f"{rounded:>8d}  {exact:>8d}  {ok:>5s}  ({elapsed:.1f}s)")

    return coeffs


# ============================================================================
# Part 4: Differential Approach
# ============================================================================

def differential_approach():
    """
    Compute pi(x) - pi(y) via Riemann differences.
    Error may cancel better for nearby x, y.
    """
    print("\n" + "=" * 90)
    print("PART 4: DIFFERENTIAL APPROACH")
    print("=" * 90)

    zeros = load_or_compute_zeros(500)

    # Test: how well does pi_N(x) - pi_N(y) approximate pi(x) - pi(y)?
    print("\n--- Differential accuracy: pi_500(x) - pi_500(y) vs pi(x) - pi(y) ---")

    base_points = [10000, 50000, 100000]
    offsets = [5000, 10000, 50000]

    print(f"{'y':>10s}  {'x':>10s}  {'exact diff':>10s}  {'Riemann diff':>14s}  {'error':>10s}  {'rounded ok?':>12s}")
    print("-" * 72)

    # Cache Riemann values
    riemann_cache = {}

    for y in base_points:
        for dx in offsets:
            x = y + dx
            if x > 1000000:
                continue

            if y not in riemann_cache:
                riemann_cache[y] = pi_riemann(y, num_zeros=500, zeros=zeros)
            if x not in riemann_cache:
                riemann_cache[x] = pi_riemann(x, num_zeros=500, zeros=zeros)

            exact_y = _prime_count(y)
            exact_x = _prime_count(x)
            exact_diff = exact_x - exact_y

            riemann_diff = riemann_cache[x] - riemann_cache[y]
            error = exact_diff - riemann_diff
            rounded_ok = round(riemann_diff) == exact_diff

            print(f"{y:>10d}  {x:>10d}  {exact_diff:>10d}  {riemann_diff:>14.4f}  "
                  f"{error:>+10.4f}  {'YES' if rounded_ok else 'NO':>12s}")

    # "Walking" approach: start from known pi(10^5) = 9592, walk in steps
    print("\n--- Walking approach: start from pi(100000) = 9592, walk in steps ---")

    start_x = 100000
    start_pi = 9592
    step = 100000
    targets = list(range(start_x + step, 1000001, step))

    print(f"\nStep size: {step}")
    print(f"{'x':>10s}  {'walked pi(x)':>12s}  {'exact pi(x)':>12s}  {'walk error':>12s}  {'direct 500 err':>15s}")
    print("-" * 65)

    current_pi = start_pi
    prev_x = start_x
    prev_riemann = pi_riemann(start_x, num_zeros=500, zeros=zeros)

    for target_x in targets:
        target_riemann = pi_riemann(target_x, num_zeros=500, zeros=zeros)
        riemann_diff = target_riemann - prev_riemann
        current_pi += round(riemann_diff)

        exact = _prime_count(target_x)
        walk_error = current_pi - exact
        direct_error = round(target_riemann) - exact

        print(f"{target_x:>10d}  {current_pi:>12d}  {exact:>12d}  {walk_error:>+12d}  {direct_error:>+15d}")

        prev_x = target_x
        prev_riemann = target_riemann

    # Finer steps
    print(f"\n--- Walking with step=10000 from pi(100000) = 9592 ---")
    step = 10000
    targets = list(range(start_x + step, 160001, step))

    current_pi = start_pi
    prev_riemann = pi_riemann(start_x, num_zeros=500, zeros=zeros)

    print(f"{'x':>10s}  {'walked pi(x)':>12s}  {'exact pi(x)':>12s}  {'cumul error':>12s}")
    print("-" * 50)

    for target_x in targets:
        target_riemann = pi_riemann(target_x, num_zeros=500, zeros=zeros)
        riemann_diff = target_riemann - prev_riemann
        current_pi += round(riemann_diff)

        exact = _prime_count(target_x)
        walk_error = current_pi - exact

        print(f"{target_x:>10d}  {current_pi:>12d}  {exact:>12d}  {walk_error:>+12d}")

        prev_riemann = target_riemann


# ============================================================================
# Part 5: Multi-N correction (use multiple zero counts to extrapolate)
# ============================================================================

def multi_N_correction():
    """
    Use pi_N(x) for several N values to extrapolate toward the true value.
    If the error E_N(x) ~ C(x)/N^alpha, then Richardson extrapolation applies.
    """
    print("\n" + "=" * 90)
    print("PART 5: MULTI-N EXTRAPOLATION (Richardson-style)")
    print("=" * 90)

    zeros = load_or_compute_zeros(500)

    test_xs = [10000, 50000, 100000, 200000]
    N_values = [100, 200, 300, 400, 500]

    print(f"\n{'x':>10s}", end="")
    for N in N_values:
        print(f"  {'pi_'+str(N):>12s}", end="")
    print(f"  {'exact':>8s}  {'Richardson':>12s}  {'Rich. err':>10s}")
    print("-" * (10 + 14 * len(N_values) + 34))

    for x_val in test_xs:
        exact = _prime_count(x_val)
        incremental = pi_riemann_incremental(x_val, zeros, max_n=500)

        pi_vals = {}
        print(f"{x_val:>10d}", end="")
        for N in N_values:
            pi_vals[N] = incremental[N]
            print(f"  {pi_vals[N]:>12.4f}", end="")

        # Richardson: if E_N ~ C/N, then pi_true ~ (N2*pi_N2 - N1*pi_N1)/(N2-N1)
        # Try with N=300 and N=500
        N1, N2 = 300, 500
        rich = (N2 * pi_vals[N2] - N1 * pi_vals[N1]) / (N2 - N1)
        rich_err = exact - rich

        print(f"  {exact:>8d}  {rich:>12.4f}  {rich_err:>+10.4f}")


# ============================================================================
# Main: run all experiments
# ============================================================================

if __name__ == "__main__":
    overall_start = time.time()

    print("HYBRID CORRECTION FOR RIEMANN EXPLICIT FORMULA")
    print(f"mpmath precision: {mpmath.mp.dps} decimal places")
    print()

    errors_300 = error_characterization()
    data, coeffs1, coeffs3 = smooth_correction()
    cal_coeffs = self_calibrating()
    differential_approach()
    multi_N_correction()

    total = time.time() - overall_start
    print(f"\n\nTotal runtime: {total:.1f}s")
