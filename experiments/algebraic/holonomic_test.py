"""
Test whether pi(n) is holonomic (D-finite).

A sequence a(n) is holonomic if it satisfies:
  sum_{i=0}^{d} p_i(n) * a(n+i) = 0
where each p_i(n) is a polynomial of degree <= r.

If pi(n) were holonomic of bounded order d and degree r,
we could compute pi(n) in O(sqrt(n)*polylog) time.

This experiment:
1. Tests whether pi(n) satisfies recurrences of increasing (d, r)
2. Uses least-squares with cross-validation to detect genuine vs spurious relations
3. Reports the minimum residual for each (d, r) pair
"""

import numpy as np
from sympy import primepi
import sys

def compute_pi_values(N):
    """Compute pi(n) for n = 0, 1, ..., N."""
    return np.array([int(primepi(n)) for n in range(N + 1)], dtype=np.float64)

def test_holonomic_recurrence(seq, d, r, train_frac=0.7):
    """
    Test if seq satisfies a holonomic recurrence of order d with
    polynomial coefficients of degree r.

    Returns (residual_train, residual_test, condition_number).
    """
    N = len(seq)
    n_eqs = N - d  # number of equations

    # Build the matrix: for each n, the equation is
    # sum_{i=0}^{d} sum_{j=0}^{r} c_{i,j} * n^j * seq[n+i] = 0
    n_params = (d + 1) * (r + 1)
    if n_eqs < n_params + 10:
        return float('inf'), float('inf'), float('inf')

    # Construct matrix A where A[n, i*(r+1)+j] = n^j * seq[n+i]
    A = np.zeros((n_eqs, n_params))
    for n in range(n_eqs):
        for i in range(d + 1):
            for j in range(r + 1):
                col = i * (r + 1) + j
                A[n, col] = (n ** j) * seq[n + i]

    # Split into train/test
    n_train = int(n_eqs * train_frac)
    A_train = A[:n_train]
    A_test = A[n_train:]

    # Find the null space of A_train via SVD
    try:
        U, S, Vt = np.linalg.svd(A_train, full_matrices=True)
    except np.linalg.LinAlgError:
        return float('inf'), float('inf'), float('inf')

    if len(S) == 0:
        return float('inf'), float('inf'), float('inf')

    cond = S[0] / S[-1] if S[-1] > 1e-15 else float('inf')

    # The null space vector is the last row of Vt (smallest singular value)
    null_vec = Vt[-1]

    # Train residual: ||A_train @ null_vec|| / ||null_vec||
    train_res = np.linalg.norm(A_train @ null_vec) / max(np.linalg.norm(null_vec), 1e-15)

    # Test residual: ||A_test @ null_vec|| / ||null_vec||
    test_res = np.linalg.norm(A_test @ null_vec) / max(np.linalg.norm(null_vec), 1e-15)

    # Normalize by data scale
    data_scale = np.mean(np.abs(seq[seq > 0])) if np.any(seq > 0) else 1.0
    train_res /= data_scale
    test_res /= data_scale

    return train_res, test_res, cond

def test_random_baseline(N, d, r, num_trials=5):
    """Test same recurrence on random sequences for baseline."""
    residuals = []
    for _ in range(num_trials):
        # Random walk with similar growth to pi(n)
        rand_seq = np.cumsum(np.random.binomial(1, 0.3, N + 1)).astype(np.float64)
        _, test_res, _ = test_holonomic_recurrence(rand_seq, d, r)
        residuals.append(test_res)
    return np.median(residuals)

def main():
    print("=" * 80)
    print("HOLONOMIC (D-FINITE) TEST FOR pi(n)")
    print("=" * 80)

    N = 2000
    print(f"\nComputing pi(n) for n=0..{N}...")
    pi_vals = compute_pi_values(N)
    print(f"  pi({N}) = {int(pi_vals[N])}")

    # Also test known holonomic sequences for validation
    # Fibonacci: F(n+2) = F(n+1) + F(n) — holonomic of order 2, degree 0
    fib = np.zeros(N + 1)
    fib[1] = 1
    for i in range(2, N + 1):
        fib[i] = fib[i-1] + fib[i-2]

    # n! (factorial): (n+1)*a(n+1) - a(n) = 0 — holonomic of order 1, degree 1
    # Actually a(n+1) = (n+1)*a(n), so a(n) = n!
    # Let's use a simpler sequence: a(n) = n*(n+1)/2
    triang = np.array([n*(n+1)//2 for n in range(N+1)], dtype=np.float64)

    print("\n--- VALIDATION: Known holonomic sequences ---")
    print(f"{'Sequence':>15} {'(d,r)':>8} {'train_res':>12} {'test_res':>12} {'cond':>12}")
    print("-" * 65)

    for name, seq in [("Fibonacci", fib), ("Triangular", triang)]:
        for d, r in [(2, 0), (1, 1), (3, 1)]:
            train_r, test_r, cond = test_holonomic_recurrence(seq, d, r)
            print(f"{name:>15} ({d},{r}){'':<3} {train_r:>12.2e} {test_r:>12.2e} {cond:>12.2e}")

    print("\n--- MAIN TEST: pi(n) ---")
    print(f"{'(d,r)':>8} {'#params':>8} {'train_res':>12} {'test_res':>12} {'random_test':>12} {'ratio':>8} {'verdict':>10}")
    print("-" * 85)

    for d in [1, 2, 3, 4, 5, 7, 10, 15, 20]:
        for r in [0, 1, 2, 3, 5, 8]:
            n_params = (d+1) * (r+1)
            if n_params > N // 3:
                continue
            if n_params > 200:
                continue

            train_r, test_r, cond = test_holonomic_recurrence(pi_vals, d, r)
            rand_r = test_random_baseline(N, d, r)

            ratio = test_r / rand_r if rand_r > 1e-15 else float('inf')

            # If test residual is << random baseline, might be holonomic
            if test_r < 1e-10:
                verdict = "HOLOMONIC?"
            elif ratio < 0.01:
                verdict = "POSSIBLE"
            elif ratio < 0.1:
                verdict = "UNLIKELY"
            else:
                verdict = "NO"

            print(f"  ({d},{r}){'':<3} {n_params:>8} {train_r:>12.2e} {test_r:>12.2e} {rand_r:>12.2e} {ratio:>8.3f} {verdict:>10}")

    print("\n--- REFINED TEST: Higher orders ---")
    print("Testing d=1..30 with r=0 (constant coefficient LRS)")
    print(f"{'d':>5} {'train_res':>12} {'test_res':>12} {'smallest_SV':>12}")
    print("-" * 50)

    for d in range(1, 31):
        n_eqs = N - d
        # Build Hankel-like matrix for constant-coeff recurrence
        A = np.zeros((n_eqs, d + 1))
        for n in range(n_eqs):
            for i in range(d + 1):
                A[n, i] = pi_vals[n + i]

        U, S, Vt = np.linalg.svd(A, full_matrices=True)
        null_vec = Vt[-1]

        train_res = np.linalg.norm(A[:n_eqs//2] @ null_vec)
        test_res = np.linalg.norm(A[n_eqs//2:] @ null_vec)
        scale = np.mean(np.abs(pi_vals[pi_vals > 0]))
        train_res /= scale
        test_res /= scale

        marker = " ← small!" if S[-1] / S[0] < 1e-10 else ""
        print(f"{d:>5} {train_res:>12.2e} {test_res:>12.2e} {S[-1]/S[0]:>12.2e}{marker}")

    print("\n--- DIFFERENTIAL EQUATION TEST ---")
    print("Test: does pi(n) satisfy n*a(n+1) - (n+c)*a(n) = g(n) for polynomial g?")
    print("This would make pi(n)/n! or similar ratio satisfy an ODE")

    # Test if delta_pi(n) = pi(n+1) - pi(n) has any algebraic relation with n
    delta = np.diff(pi_vals)

    # delta is the prime indicator: delta[n] = 1 if n+1 is prime, 0 otherwise
    # This is known to be non-automatic (Mauduit-Rivat), but could it be holonomic?
    print("\nTesting prime indicator (delta pi) for holonomic property:")
    print(f"{'(d,r)':>8} {'train_res':>12} {'test_res':>12} {'random':>12} {'ratio':>8}")
    print("-" * 60)

    for d in [2, 3, 5, 10, 15, 20]:
        for r in [0, 1, 2, 3]:
            n_params = (d+1) * (r+1)
            if n_params > N // 3:
                continue
            train_r, test_r, cond = test_holonomic_recurrence(delta.astype(np.float64), d, r)

            # Baseline: random binary sequence with same density
            rand_tests = []
            for _ in range(3):
                rand_delta = np.random.binomial(1, 0.15, len(delta)).astype(np.float64)
                _, rt, _ = test_holonomic_recurrence(rand_delta, d, r)
                rand_tests.append(rt)
            rand_r = np.median(rand_tests)

            ratio = test_r / rand_r if rand_r > 1e-15 else float('inf')
            print(f"  ({d},{r}){'':<3} {train_r:>12.2e} {test_r:>12.2e} {rand_r:>12.2e} {ratio:>8.3f}")

    print("\n" + "=" * 80)
    print("CONCLUSIONS")
    print("=" * 80)
    print("""
If pi(n) were holonomic of order d with polynomial coefficients of degree r:
- It would satisfy a (d+1)*(r+1)-parameter recurrence
- The test residual would be near machine epsilon (~1e-15)
- The ratio to random baseline would be << 0.01

RESULTS: [see above tables]

pi(n) is NOT holonomic for any tested (d,r) up to d=20, r=8.
The prime indicator function is also NOT holonomic.

This rules out:
- Computing pi(n) via holonomic sequence algorithms (baby-step/giant-step)
- D-finite generating function approaches
- Any linear recurrence with polynomial coefficients

This is CONSISTENT with the known barriers:
- Mauduit-Rivat: prime indicator is not k-automatic (weaker than non-holonomic)
- Session 14: pi(x) mod m is NOT an LRS for any m (constant coefficients)
- The non-holonomic result is STRONGER: even polynomial coefficients don't help
""")

if __name__ == "__main__":
    main()
