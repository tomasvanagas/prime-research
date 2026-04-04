"""
Session 14: Can pi(x) be expressed as a matrix power entry or trace?

If there exists a fixed-size matrix A over Z such that (A^x)_{ij} = pi(x)
for some i,j, then pi(x) would be computable in TC^0 (scalar powering = TC^0,
fixed-dim matrix powering = TC^0 by Mereghetti-Palano 2000).

Key insight: (A^n)_{ij} = sum over paths of length n from i to j, weighted
by edge products. This is a LINEAR RECURRENCE SEQUENCE (LRS).

OBSTRUCTION: Mauduit-Rivat (2010) proved the prime indicator is NOT k-automatic.
And Skolem-Mahler-Lech: an integer LRS is zero iff it's eventually periodic in
its zero set. So (A^n)_{ij} cannot equal 1_prime(n) (the prime indicator isn't
eventually periodic in its support).

BUT: pi(x) = sum_{k<=x} 1_prime(k) is the CUMULATIVE SUM. A cumulative sum
of a non-LRS might still be an LRS? No — if pi(x) were LRS, then
pi(x) - pi(x-1) = 1_prime(x) would also be LRS. Contradiction.

HOWEVER: what about (A^x)_{ij} mod m for carefully chosen A and m? Or
working in a quotient ring?

Let's also test: can pi(x) mod m be a linear recurrence mod m for small m?

And explore: perhaps not (A^x)_{ij} = pi(x), but rather a MORE COMPLEX
extraction like Tr(A^x * B) or det(I - z*A^x) at z=1?
"""

import numpy as np
from sympy import primepi, isprime, Matrix
from itertools import product as cartesian_product

def check_lrs_modm(seq, m, max_order=20):
    """
    Check if sequence mod m satisfies a linear recurrence of order <= max_order.
    Returns (True, order, coefficients) or (False, None, None).
    """
    seq_mod = [s % m for s in seq]
    n = len(seq_mod)

    for order in range(1, min(max_order + 1, n // 2)):
        # Build system: c[0]*a[k] + c[1]*a[k-1] + ... + c[order-1]*a[k-order+1] = a[k+1] (mod m)
        # for k = order-1, ..., n-2
        rows = n - order - 1
        if rows < order:
            continue

        # Check if the recurrence holds for all windows
        valid = True
        # Try to find recurrence from first 'order+1' windows
        A_mat = []
        b_vec = []
        for k in range(order, min(2 * order + 5, n)):
            row = [seq_mod[k - j - 1] % m for j in range(order)]
            A_mat.append(row)
            b_vec.append(seq_mod[k] % m)

        # Try all possible coefficient vectors mod m (brute force for small m, order)
        if m <= 5 and order <= 4:
            found = False
            for coeffs in cartesian_product(range(m), repeat=order):
                ok = True
                for k in range(order, n):
                    pred = sum(coeffs[j] * seq_mod[k - j - 1] for j in range(order)) % m
                    if pred != seq_mod[k]:
                        ok = False
                        break
                if ok:
                    found = True
                    return True, order, coeffs
            if not found:
                continue
        else:
            # Heuristic: use numpy
            A_np = np.array(A_mat, dtype=float)
            b_np = np.array(b_vec, dtype=float)
            if A_np.shape[0] >= order:
                try:
                    coeffs, res, _, _ = np.linalg.lstsq(A_np[:order], b_np[:order], rcond=None)
                    coeffs_int = [round(c) % m for c in coeffs]
                    # Verify
                    ok = True
                    for k in range(order, n):
                        pred = sum(coeffs_int[j] * seq_mod[k - j - 1] for j in range(order)) % m
                        if pred != seq_mod[k]:
                            ok = False
                            break
                    if ok:
                        return True, order, coeffs_int
                except:
                    pass

    return False, None, None


def test_pi_lrs():
    """Test whether pi(x) mod m is a linear recurrence for various m."""
    N = 200
    pi_values = [primepi(x) for x in range(N)]

    print("Testing pi(x) mod m as linear recurrence sequence")
    print("=" * 60)

    for m in [2, 3, 4, 5, 6, 7, 8, 10, 12, 16]:
        is_lrs, order, coeffs = check_lrs_modm(pi_values, m, max_order=8)
        if is_lrs:
            print(f"pi(x) mod {m}: LRS of order {order}, coeffs = {coeffs}")
        else:
            print(f"pi(x) mod {m}: NOT LRS (order <= 8)")


def search_matrix_encoding():
    """
    Brute-force search: find a small integer matrix A such that
    Tr(A^x) or (A^x)_{0,0} approximates pi(x).

    Since pi(x) grows like x/ln(x), and eigenvalue powers grow exponentially
    (or polynomially if eigenvalues on unit circle), we'd need eigenvalues
    that somehow encode prime distribution.

    We know this can't work exactly (LRS argument), but let's see how close
    we can get for small matrices.
    """
    print("\n" + "=" * 60)
    print("Searching for matrix encodings of pi(x)")
    print("=" * 60)

    N = 50  # test range
    target = [int(primepi(x)) for x in range(N)]

    # For 2x2 matrices over small integers
    best_err = float('inf')
    best_A = None
    best_extract = None

    tested = 0
    for a in range(-3, 4):
        for b in range(-3, 4):
            for c in range(-3, 4):
                for d in range(-3, 4):
                    A = np.array([[a, b], [c, d]], dtype=float)
                    # Try (A^x)_{0,0} and Tr(A^x)
                    powers = [np.linalg.matrix_power(A.astype(int), x) for x in range(N)]

                    for extract_name, extract_fn in [("(0,0)", lambda P: int(P[0,0])),
                                                      ("trace", lambda P: int(np.trace(P)))]:
                        try:
                            vals = [extract_fn(P) for P in powers]
                            err = sum(abs(vals[x] - target[x]) for x in range(2, N))
                            if err < best_err:
                                best_err = err
                                best_A = A.copy()
                                best_extract = extract_name
                        except (OverflowError, ValueError):
                            pass
                    tested += 1

    print(f"Tested {tested} 2x2 matrices")
    print(f"Best: A = {best_A.astype(int).tolist()}, extract = {best_extract}, total error = {best_err}")

    if best_A is not None:
        A = best_A.astype(int)
        print("\nComparison (x, pi(x), A^x value):")
        for x in range(2, min(20, N)):
            P = np.linalg.matrix_power(A, x)
            if best_extract == "(0,0)":
                val = int(P[0,0])
            else:
                val = int(np.trace(P))
            print(f"  x={x:3d}: pi(x)={target[x]:4d}, A^x={val:10d}, diff={val-target[x]:+d}")


def analyze_pi_differences():
    """
    Analyze the difference sequence delta(x) = pi(x) - pi(x-1) = 1_prime(x).
    Also higher-order differences.
    If pi(x) satisfied any polynomial recurrence, the differences would satisfy
    a lower-order one.
    """
    print("\n" + "=" * 60)
    print("Analyzing pi(x) difference structure")
    print("=" * 60)

    N = 100
    pi_vals = [primepi(x) for x in range(N)]

    # 1st differences: delta_1(x) = pi(x) - pi(x-1) = 1_prime(x)
    d1 = [pi_vals[x] - pi_vals[x-1] for x in range(1, N)]

    # 2nd differences: delta_2(x) = delta_1(x) - delta_1(x-1)
    d2 = [d1[x] - d1[x-1] for x in range(1, len(d1))]

    print(f"1st differences (prime indicator): {d1[:30]}")
    print(f"2nd differences: {d2[:30]}")

    # Autocorrelation of 1st differences
    d1_arr = np.array(d1[:80], dtype=float)
    d1_centered = d1_arr - d1_arr.mean()
    autocorr = np.correlate(d1_centered, d1_centered, mode='full')
    autocorr = autocorr[len(d1_centered)-1:]  # Keep only non-negative lags
    autocorr /= autocorr[0]
    print(f"\nAutocorrelation of prime indicator (lags 0-10): {autocorr[:11].round(4)}")

    # Key question: is there ANY linear recurrence in the differences?
    # Test: correlation between d1[x] and d1[x-k] for various k
    print("\nLag correlations of prime indicator:")
    for k in range(1, 15):
        corr = np.corrcoef(d1_arr[k:], d1_arr[:-k])[0,1]
        print(f"  lag {k:2d}: r = {corr:+.4f}")


def companion_matrix_search():
    """
    The companion matrix approach: if pi(x) satisfies a recurrence
    pi(x) = c_1*pi(x-1) + c_2*pi(x-2) + ... + c_k*pi(x-k) + f(x)
    then pi(x) = (C^x * v)_0 where C is the companion matrix and v is
    the initial conditions.

    Test whether pi(x) is CLOSE to satisfying such a recurrence for
    small k, and characterize the residual.
    """
    print("\n" + "=" * 60)
    print("Companion matrix / linear recurrence fitting")
    print("=" * 60)

    N = 200
    pi_vals = np.array([primepi(x) for x in range(N)], dtype=float)

    for order in [2, 3, 4, 5, 10, 20]:
        if order + 10 > N:
            continue
        # Fit recurrence from training data
        train = pi_vals[:N//2]
        test = pi_vals[N//2:]

        # Build system A * c = b
        A_rows = []
        b_vec = []
        for x in range(order, len(train)):
            row = [train[x - j - 1] for j in range(order)]
            A_rows.append(row)
            b_vec.append(train[x])

        A_mat = np.array(A_rows)
        b_arr = np.array(b_vec)

        coeffs, residuals, _, _ = np.linalg.lstsq(A_mat, b_arr, rcond=None)

        # Test on training data
        train_errs = []
        for x in range(order, len(train)):
            pred = sum(coeffs[j] * train[x - j - 1] for j in range(order))
            train_errs.append(abs(pred - train[x]))

        # Test on test data (extrapolate)
        test_preds = list(pi_vals[:N//2])
        test_errs = []
        for x in range(N//2, N):
            pred = sum(coeffs[j] * test_preds[x - j - 1] for j in range(order))
            test_preds.append(pred)
            test_errs.append(abs(pred - pi_vals[x]))

        print(f"\nOrder {order}: coeffs = {coeffs.round(6)[:5]}{'...' if order > 5 else ''}")
        print(f"  Train MAE: {np.mean(train_errs):.4f}, max: {max(train_errs):.4f}")
        print(f"  Test  MAE: {np.mean(test_errs):.4f}, max: {max(test_errs):.4f}")
        print(f"  Test exact matches: {sum(1 for e in test_errs if e < 0.5)}/{len(test_errs)}")


if __name__ == '__main__':
    test_pi_lrs()
    search_matrix_encoding()
    analyze_pi_differences()
    companion_matrix_search()
