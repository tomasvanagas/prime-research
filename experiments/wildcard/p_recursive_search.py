"""
P-recursive (holonomic) structure search for pi(n) and related sequences.

IDEA
----
A sequence a_n is *P-recursive* (a.k.a. D-finite or holonomic) if there exist
polynomials p_0, ..., p_k of bounded degree d such that

    sum_{j=0}^k p_j(n) a_{n+j} = 0      for all n >= n_0.

P-recursive sequences include n!, Catalan numbers, Apery numbers, and any
sequence that's a polynomial-coefficient finite linear combination of its
shifts. Crucially:

    Holonomic sequence -> evaluable at index N in O(sqrt(N) log^c N) ops
    (Bostan-Gaudry-Schost factorial trick), and at all indices up to N in
    O(N) ops with O(d * k) memory.

If pi(n) — or pi(n) - Li(n), or some simple transform — were P-recursive
with low order k and low degree d, we'd have an immediate sub-linear
algorithm. The search for such a recurrence is a closed numerical procedure:
solve a linear system over n_0+M equations in (k+1)(d+1) unknowns.

This is the FRESH angle: no entry in the wildcard list specifically tests
holonomic / P-recursive structure. The closest is `p_recursive_search` ...
which doesn't exist yet. Other tests (LFSR, transfer matrix, automata)
restrict to *constant-coefficient* recurrences, which are strictly weaker.

EXPECTATION
-----------
Failure mode is almost certainly **C (circularity) or I (information)**:
- If pi(n) were P-recursive of order k and degree d, the residual after
  fitting on [1, M] would extend correctly to [M, 2M]. Random-like sequences
  do NOT admit such extension at low (k, d).

But the QUANTITATIVE answer matters: we want to plot the smallest residual
(order * degree) across (k, d) grid and confirm it grows linearly with M.
"""

import numpy as np
from sympy import primepi
from mpmath import li
import time


def fit_p_recursive(a: np.ndarray, k: int, d: int):
    """
    Fit a P-recurrence of order k, polynomial degree d:
        sum_{j=0}^k p_j(n) a[n+j] = 0
    where p_j(n) = sum_{r=0}^d c_{j,r} n^r.

    Set up linear system: for each test row n in [n0, n0+M),
        sum_{j=0}^k sum_{r=0}^d c_{j,r} n^r a[n+j] = 0

    Solve homogeneous system in unknowns c_{j,r}; report smallest singular
    value (= residual of best non-trivial recurrence).
    """
    N = len(a)
    n_unknowns = (k + 1) * (d + 1)
    n0 = 0
    n_eqns = N - k - n0
    if n_eqns < n_unknowns + 5:
        return None
    A_mat = np.zeros((n_eqns, n_unknowns))
    for row, n in enumerate(range(n0, n0 + n_eqns)):
        col = 0
        for j in range(k + 1):
            for r in range(d + 1):
                A_mat[row, col] = (n ** r) * a[n + j]
                col += 1
    # Solve via SVD; smallest singular value = residual
    # Normalize per-row (so coefficient magnitudes don't dominate)
    row_norms = np.linalg.norm(A_mat, axis=1, keepdims=True)
    row_norms[row_norms == 0] = 1.0
    A_norm = A_mat / row_norms
    s = np.linalg.svd(A_norm, compute_uv=False)
    return s[-1], s[0]  # min and max singular value


def fit_recurrence(a_train: np.ndarray, k: int, d: int):
    """Fit homogeneous P-recurrence on a_train; return polynomial coefficients."""
    n_unknowns = (k + 1) * (d + 1)
    n_eqns = len(a_train) - k
    if n_eqns < n_unknowns + 5:
        return None
    A_train = np.zeros((n_eqns, n_unknowns))
    for row, n in enumerate(range(n_eqns)):
        col = 0
        for j in range(k + 1):
            for r in range(d + 1):
                A_train[row, col] = (n ** r) * a_train[n + j]
                col += 1
    _, _, Vt = np.linalg.svd(A_train, full_matrices=False)
    return Vt[-1].reshape((k + 1, d + 1))


def predict_one_step(a: np.ndarray, C: np.ndarray, n: int):
    """
    Use recurrence with coefficients C (shape (k+1, d+1)) at index n to predict
    a[n + k] from a[n], a[n+1], ..., a[n+k-1].
    """
    k = C.shape[0] - 1
    d = C.shape[1] - 1
    n_powers = np.array([n ** r for r in range(d + 1)])
    p_vals = C @ n_powers
    if abs(p_vals[k]) < 1e-12:
        return None
    return -sum(p_vals[j] * a[n + j] for j in range(k)) / p_vals[k]


def predictive_residual(a: np.ndarray, k: int, d: int, train_frac: float = 0.6):
    """One-step prediction error in the held-out tail region."""
    N = len(a)
    M = int(train_frac * N)
    a_train = a[:M]
    C = fit_recurrence(a_train, k, d)
    if C is None:
        return None
    test_errs = []
    for n in range(M, N - k):
        a_pred = predict_one_step(a, C, n)
        if a_pred is None:
            continue
        a_true = a[n + k]
        denom = max(abs(a_true), 1.0)
        test_errs.append(abs(a_pred - a_true) / denom)
    if not test_errs:
        return None
    return np.median(test_errs), np.mean(test_errs)


def extrapolation_test(a: np.ndarray, k: int, d: int, train_end: int):
    """
    Fit on [0, train_end). For each test index n in [train_end, N-k), use
    recurrence iteratively (rolling forward, replacing a[n+k] with prediction)
    to assess if the recurrence captures real long-range structure.

    Returns (one-step errors, rollout abs errors) keyed by n.
    """
    a_train = a[:train_end]
    C = fit_recurrence(a_train, k, d)
    if C is None:
        return None
    one_step_errs = {}
    rollout_errs = {}
    a_rollout = a.copy()
    for n in range(train_end - k, len(a) - k):
        a_pred = predict_one_step(a, C, n)
        if a_pred is None:
            continue
        a_true = a[n + k]
        one_step_errs[n + k] = abs(a_pred - a_true)
        # rollout: use predicted values going forward
        a_rollout_pred = predict_one_step(a_rollout, C, n)
        if a_rollout_pred is None:
            a_rollout[n + k] = a_true  # fallback
        else:
            a_rollout[n + k] = a_rollout_pred
        rollout_errs[n + k] = abs(a_rollout[n + k] - a_true)
    return one_step_errs, rollout_errs


def main():
    print("P-RECURSIVE STRUCTURE SEARCH for pi(n) and related sequences")
    print("=" * 72)

    # Build sequences
    N_max = 600
    print(f"Building sequences up to n = {N_max} ...")
    pi_seq = np.array([int(primepi(n)) for n in range(1, N_max + 1)], dtype=float)
    li_seq = np.array([float(li(n)) for n in range(2, N_max + 2)], dtype=float)
    delta_seq = pi_seq - li_seq[:len(pi_seq)]  # pi(n) - li(n)
    diff_seq = np.diff(pi_seq)  # 1 if n+1 prime else 0

    sequences = {
        "pi(n)": pi_seq,
        "li(n)": li_seq[: len(pi_seq)],
        "delta(n) = pi(n)-li(n)": delta_seq,
        "indicator chi_P(n+1)": diff_seq,
    }

    # Reference: fit a known holonomic sequence (n!) for sanity (small range
    # to avoid overflow). n! satisfies a_{n+1} = n a_n, i.e., p_1(n)=1, p_0(n)=-n.
    import math as _math
    fact_small = np.array([float(_math.factorial(n)) for n in range(20)])
    sequences["n! small (control)"] = fact_small

    print()
    print("PART A: Fit residual on a window, then in-window prediction error")
    print(f"{'sequence':<28} {'k':>2} {'d':>2} {'fit_resid':>12} {'cond':>10} {'pred_rel_err':>14}")
    print("-" * 72)

    for name, seq in sequences.items():
        # Use a smaller window for clarity
        if "small" in name:
            seq_use = seq[:20]
        else:
            seq_use = seq[:200]

        for (k, d) in [(1, 1), (2, 2), (3, 3), (4, 3)]:
            try:
                fit = fit_p_recursive(seq_use, k, d)
                if fit is None:
                    continue
                s_min, s_max = fit
                cond = s_max / max(s_min, 1e-300)
                pred = predictive_residual(seq_use, k, d, train_frac=0.6)
                if pred is None:
                    pred_str = "n/a"
                else:
                    pred_str = f"{pred[0]:.4g}"
                print(f"{name:<28} {k:>2} {d:>2} {s_min:>12.4g} {cond:>10.2e} {pred_str:>14}")
            except Exception as e:
                print(f"{name:<28} {k:>2} {d:>2}  ERROR {e}")
        print()

    print()
    print("PART B: Extrapolation test — fit on [0,200), predict on [200, 600)")
    print("        (the real test of structural recurrence)")
    print(f"{'sequence':<28} {'k':>2} {'d':>2} {'mean_1step':>14} {'max_1step':>14} {'mean_rollout':>14}")
    print("-" * 90)

    for name, seq in sequences.items():
        if "small" in name:
            continue  # only 20 elems
        seq_use = seq[:600]
        for (k, d) in [(2, 2), (3, 3), (4, 3)]:
            res = extrapolation_test(seq_use, k, d, train_end=200)
            if res is None:
                continue
            one_step, rollout = res
            if not one_step:
                continue
            m1 = np.mean(list(one_step.values()))
            x1 = np.max(list(one_step.values()))
            mr = np.mean(list(rollout.values()))
            print(f"{name:<28} {k:>2} {d:>2} {m1:>14.4g} {x1:>14.4g} {mr:>14.4g}")
        print()

    print()
    print("PART C: Round-to-nearest-integer recovery on the extrapolated region")
    print("        (does the recurrence give pi(n) exactly?)")
    print(f"{'sequence':<28} {'k':>2} {'d':>2} {'frac_exact':>14} {'mean_abs_err':>14}")
    print("-" * 90)

    for name, seq in [("pi(n)", pi_seq), ("delta(n) = pi(n)-li(n)", delta_seq)]:
        seq_use = seq[:600]
        for (k, d) in [(2, 2), (3, 3), (4, 3), (5, 4)]:
            res = extrapolation_test(seq_use, k, d, train_end=200)
            if res is None:
                continue
            one_step, _ = res
            if not one_step:
                continue
            ns = list(one_step.keys())
            errs_int = []
            n_exact = 0
            for n in ns:
                C = fit_recurrence(seq_use[:200], k, d)
                a_pred = predict_one_step(seq_use, C, n - k)
                if a_pred is None:
                    continue
                a_true = seq_use[n]
                if name == "pi(n)":
                    if int(round(a_pred)) == int(a_true):
                        n_exact += 1
                else:  # delta -> reconstruct pi
                    pi_pred = a_pred + float(li(n + 1))
                    if int(round(pi_pred)) == int(pi_seq[n]):
                        n_exact += 1
                errs_int.append(abs(a_pred - a_true))
            if errs_int:
                frac = n_exact / len(errs_int)
                print(f"{name:<28} {k:>2} {d:>2} {frac:>14.4g} {np.mean(errs_int):>14.4g}")
        print()


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"\nTotal time: {time.time() - t0:.1f}s")
    print("\nINTERPRETATION")
    print("-" * 72)
    print("- For a true holonomic sequence (n!), small (k,d) should give tiny")
    print("  fit residual AND tiny one-step predictive error on hold-out.")
    print("- For pi(n) / delta(n): if predictive error remains O(1) regardless of")
    print("  (k,d), the sequence is NOT P-recursive at the orders tested.")
    print("- A surprise win would be predictive_rel_err << 1 at some moderate (k,d).")
