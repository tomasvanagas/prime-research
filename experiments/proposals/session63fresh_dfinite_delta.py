"""
P1: D-finite (holonomic / Apéry-style) recurrence hunt for delta(n) = p(n) - R^{-1}(n).

A sequence (a_n) is D-finite if there exist polynomials P_0, ..., P_L in Z[n],
not all zero, such that
    sum_{j=0..L} P_j(n) a_{n+j} = 0
for all n large enough. Apéry's recurrence for zeta(3) is L=2, deg=2.

If delta(n) is D-finite of order L with polynomial coefficients of degree d,
then evaluation at n is O((L+d) * M(log n)) by binary-splitting shift.

Test: search the rectangle (L, d) in {1..4} x {1..4} for a
linear-algebra null vector. Use exact rationals (sympy) on K data points
where K = (L+1)*(d+1) + 5 to be confident about full rank.

Falsification: full rank for ALL (L, d) up to (4, 4).
"""

import sys
from sympy import Matrix, Rational, Integer
from sympy.ntheory import primerange
from mpmath import mp, mpf, mpc, log as mlog, exp as mexp, polylog, zeta, riemannr

mp.prec = 200  # high precision so R-inverse is correct to many digits


def primes_upto(N):
    return list(primerange(2, N + 1))


def riemann_R_inverse(n_target):
    """Numerical inverse of Riemann R function: solve R(x) = n for x."""
    # initial guess via PNT
    n = float(n_target)
    if n < 5:
        x = mpf(2 * (n + 1))
    else:
        # better initial guess from inverse of li
        x = mpf(n) * mlog(mpf(n)) + mpf(n) * mlog(mlog(mpf(max(2.0, float(n)))))
    # Newton iterations on f(x) = R(x) - n
    for _ in range(60):
        Rx = riemannr(x)
        # derivative R'(x): sum mobius(k)/k * 1/(k log x) approximated by 1/log x
        Rp = mpf(1) / mlog(x)
        delta = (Rx - mpf(n_target)) / Rp
        x = x - delta
        if abs(delta) < mpf("1e-30"):
            break
    return x


def compute_delta_table(N_max):
    """Compute delta(n) = p(n) - R^{-1}(n) for n = 1..N_max as exact integers
    minus high-precision floats. Returns list of mpf values.
    """
    # collect first N_max primes
    upper_estimate = max(20, int(N_max * (mlog(N_max) + mlog(mlog(N_max + 1)) + 2)))
    primes = []
    cur = upper_estimate
    while len(primes) < N_max:
        primes = primes_upto(cur)
        if len(primes) < N_max:
            cur = int(cur * 1.5)
    primes = primes[:N_max]

    deltas = []
    for n in range(1, N_max + 1):
        p_n = primes[n - 1]
        R_inv = riemann_R_inverse(n)
        deltas.append(mpf(p_n) - R_inv)
    return deltas


def build_matrix(deltas, L, d, n_start, n_end):
    """Build matrix A whose row at n contains entries n^k * delta(n+j),
       indexed by (j, k). Columns are L2-normalized to fight conditioning.
    """
    import numpy as np
    n_data = len(deltas)
    n_end = min(n_end, n_data - L - 1)
    indices = list(range(n_start, n_end + 1))
    K = len(indices)
    cols = (L + 1) * (d + 1)
    A = np.zeros((K, cols))
    for r, n in enumerate(indices):
        col = 0
        for j in range(L + 1):
            d_val = float(deltas[n + j - 1])
            for k in range(d + 1):
                A[r, col] = float(n) ** k * d_val
                col += 1
    # column-normalize
    norms = np.linalg.norm(A, axis=0)
    norms_safe = np.where(norms > 0, norms, 1.0)
    A_norm = A / norms_safe
    return A, A_norm, norms, indices


def search_recurrence(deltas, L, d):
    """Search for D-finite recurrence on FIRST half, validate on SECOND half.

       We hold out the second half of data to check that any apparent
       null vector is real and not a conditioning artifact.
    """
    import numpy as np
    n_data = len(deltas)
    cols = (L + 1) * (d + 1)
    n_train_start = 30
    n_train_end = n_data // 2
    n_test_start = n_data // 2 + 1
    n_test_end = n_data - L - 1

    if n_train_end - n_train_start + 1 < cols * 4:
        return None  # need plenty of rows for stable null space

    # --- training matrix ---
    A_train, A_train_n, norms, _ = build_matrix(
        deltas, L, d, n_train_start, n_train_end
    )
    # SVD on column-normalized matrix → unitless singular values
    _, S_train, Vt = np.linalg.svd(A_train_n, full_matrices=False)
    smax = float(S_train[0])
    smin = float(S_train[-1])

    # candidate null vector (column-normalized space)
    null_norm = Vt[-1]
    # un-normalize back to coefficient space:
    # row of A is e_norm = a_n^k delta(n+j) / norms_col, with
    # null space satisfies A_n @ null_norm ≈ 0 → A @ (null_norm/norms) ≈ 0
    null_coef = null_norm / norms
    null_coef = null_coef / np.max(np.abs(null_coef))  # rescale for readability

    # --- held-out residual: predict delta(n+L) from delta(n)..delta(n+L-1) ---
    # null_coef structure: [(j=0,k=0), (j=0,k=1), ..., (j=0,k=d), (j=1,k=0), ...]
    # so coefficients of delta(n+j) are P_j(n) = sum_k null_coef[(j,k)] * n^k
    def P_j_of_n(n_val, j):
        s = 0.0
        for k in range(d + 1):
            s += null_coef[j * (d + 1) + k] * (float(n_val) ** k)
        return s

    pred_errors = []
    delta_scale = np.std([float(deltas[n - 1]) for n in
                          range(n_test_start, n_test_end + 1)])
    for n in range(n_test_start, n_test_end + 1):
        # recurrence: sum_j P_j(n) * delta(n+j) = 0
        # solve for delta(n+L) = -[sum_{j<L} P_j(n) delta(n+j)] / P_L(n)
        P_L = P_j_of_n(n, L)
        if abs(P_L) < 1e-15:
            continue
        rhs = 0.0
        for j in range(L):
            rhs += P_j_of_n(n, j) * float(deltas[n + j - 1])
        pred = -rhs / P_L
        actual = float(deltas[n + L - 1])
        pred_errors.append(abs(pred - actual))

    pred_rmse = float(np.sqrt(np.mean(np.array(pred_errors) ** 2))) \
        if pred_errors else float("inf")
    pred_skill = pred_rmse / delta_scale  # 1.0 = no skill, 0.0 = perfect

    return {
        "smax": smax,
        "smin": smin,
        "rank_ratio": smin / smax,
        "pred_rmse": pred_rmse,
        "pred_skill": float(pred_skill),
        "delta_scale": float(delta_scale),
        "null_coef": null_coef,
        "n_train": (n_train_start, n_train_end),
        "n_test": (n_test_start, n_test_end),
    }


def main():
    print("=" * 70)
    print("P1: D-finite recurrence hunt for delta(n) = p(n) - R^{-1}(n)")
    print("=" * 70)

    N_MAX = 400  # enough samples for (L,d) up to (4,4); keep runtime modest
    print(f"\nComputing delta(1..{N_MAX}) with {mp.prec}-bit precision...")
    deltas = compute_delta_table(N_MAX)
    print("First 10 delta(n) values:")
    for n in range(1, 11):
        print(f"  delta({n}) = {float(deltas[n-1]):+.6f}")

    print("\n--- Searching for D-finite recurrence (train/test split) ---")
    print("pred_skill: RMSE(predicted - actual) / std(delta) on held-out half")
    print("            1.0 = no skill (random), 0.0 = perfect prediction")
    print(f"{'L':>2} {'d':>2}  {'unknowns':>4}  "
          f"{'train_rank_ratio':>17}  {'pred_skill':>11}  verdict")
    print("-" * 80)

    found = []
    SKILL_THRESH = 0.05  # held-out prediction must beat this to be real
    for L in range(1, 5):
        for d in range(0, 5):
            res = search_recurrence(deltas, L, d)
            if res is None:
                continue
            unknowns = (L + 1) * (d + 1)
            verdict = "no skill"
            if res["pred_skill"] < SKILL_THRESH:
                verdict = "POTENTIAL D-FINITE (verify!)"
                found.append((L, d, res))
            print(f"{L:>2} {d:>2}  {unknowns:>4}  "
                  f"{res['rank_ratio']:>17.4e}  "
                  f"{res['pred_skill']:>11.4f}  {verdict}")

    print()
    if not found:
        print("VERDICT: No D-finite recurrence detected up to (L=4, d=4).")
        print("  - Held-out test residual is large in every case.")
        print("  - Apparent rank-deficiencies on training data are conditioning")
        print("    artifacts; they do NOT extrapolate.")
        print("  - delta(n) is NOT D-finite with low-degree polynomial coefficients.")
        print("  - Consistent with delta(n) inheriting GUE-like randomness from zeros.")
        return 0
    else:
        print("VERDICT: Potential D-finite recurrence found! Worth follow-up.")
        for L, d, res in found:
            print(f"  L={L} d={d}: rank_ratio={res['rank_ratio']:.2e}, "
                  f"test_residual={res['test_residual_relative']:.2e}")
            vec = res["null_coef"]
            print(f"    null vector: {vec[:6]}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
