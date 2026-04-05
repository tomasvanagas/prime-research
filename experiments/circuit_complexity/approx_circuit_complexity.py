"""
Approximate Circuit Complexity of pi(x) mod 2

Novel experiment: How does circuit complexity degrade with error tolerance?

For each N, we measure:
1. Minimum degree polynomial threshold function (PTF) for exact computation
2. Accuracy of degree-d PTFs for d = 1, 2, ..., N
3. Phase transition: at what error tolerance does complexity drop to polynomial?
4. Comparison of smooth approximation accuracy vs PTF accuracy
5. Singular value spectrum of the truth table matrix (communication matrix)
   for UNBALANCED bit partitions

Key questions:
- Is there a sharp phase transition in complexity at some error threshold?
- Does the smooth part R(x) explain the "easy" fraction of the function?
- For unbalanced partitions, are the top singular vectors polylog-computable?

Session 35 experiment.
"""

import numpy as np
from sympy import isprime
import sys
from scipy.optimize import linprog
from scipy.linalg import svd
import time
import math

def compute_pi(x):
    """Compute pi(x) = number of primes <= x."""
    if x < 2:
        return 0
    count = 0
    for i in range(2, x + 1):
        if isprime(i):
            count += 1
    return count

def build_truth_tables(N):
    """Build truth tables for pi(x) mod 2 and is_prime(x) for N-bit inputs."""
    size = 2**N
    pi_vals = np.zeros(size, dtype=int)
    prime_indicator = np.zeros(size, dtype=int)

    # Sieve for efficiency
    is_p = np.zeros(size, dtype=bool)
    for i in range(2, size):
        if i < 4:
            is_p[i] = (i >= 2)
        else:
            is_p[i] = True
    # Simple sieve
    for i in range(2, int(size**0.5) + 1):
        if is_p[i]:
            for j in range(i*i, size, i):
                is_p[j] = False

    cumsum = 0
    for i in range(size):
        if is_p[i]:
            cumsum += 1
        pi_vals[i] = cumsum
        prime_indicator[i] = int(is_p[i])

    return pi_vals, prime_indicator, is_p

def R_inverse_approx(n):
    """Approximate p(n) using R^{-1}(n) via Newton's method on R(x) = n."""
    if n <= 0:
        return 2
    # Rough estimate: p(n) ~ n * ln(n)
    from math import log, exp
    if n < 2:
        return 2
    x = n * log(n)
    if x < 2:
        x = 2
    # R(x) ~ li(x) ~ x/ln(x) for large x
    # Newton: x_{k+1} = x_k - (R(x_k) - n) / R'(x_k)
    # R'(x) ~ 1/ln(x)
    for _ in range(20):
        if x <= 1:
            x = 2
            break
        lx = log(x)
        if lx < 1e-10:
            break
        rx = x / lx  # R(x) ~ x/ln(x)
        rp = 1.0 / lx  # R'(x) ~ 1/ln(x)
        x = x - (rx - n) / rp
        if x <= 1:
            x = 2
    return x

def smooth_approximation_accuracy(N, pi_mod2, pi_vals):
    """How many outputs does the smooth approximation R(x) get right?"""
    size = 2**N
    correct = 0
    for x in range(size):
        # Smooth approximation of pi(x): R(x) ~ x/ln(x) for x >= 2
        if x < 2:
            approx_pi = 0
        else:
            lx = math.log(x)
            approx_pi = round(x / lx)

        if approx_pi % 2 == pi_mod2[x]:
            correct += 1

    return correct / size

def degree_d_features(x_bits, N, degree):
    """Generate all monomial features up to given degree."""
    # x_bits: array of shape (num_samples, N) with values in {0, 1}
    # We use {-1, +1} encoding: x_i -> 2*x_i - 1
    x_pm = 2.0 * x_bits - 1.0  # Convert to {-1, +1}

    features = [np.ones(len(x_bits))]  # Constant term

    if degree >= 1:
        for i in range(N):
            features.append(x_pm[:, i])

    if degree >= 2:
        for i in range(N):
            for j in range(i+1, N):
                features.append(x_pm[:, i] * x_pm[:, j])

    if degree >= 3:
        for i in range(N):
            for j in range(i+1, N):
                for k in range(j+1, N):
                    features.append(x_pm[:, i] * x_pm[:, j] * x_pm[:, k])

    if degree >= 4 and N <= 12:
        from itertools import combinations
        for combo in combinations(range(N), 4):
            feat = np.ones(len(x_bits))
            for idx in combo:
                feat *= x_pm[:, idx]
            features.append(feat)

    return np.column_stack(features)

def fit_ptf(features, targets, regularize=1e-6):
    """Fit a polynomial threshold function using least squares.
    targets: {0, 1} -> {-1, +1}
    Returns accuracy."""
    t_pm = 2.0 * targets - 1.0  # Convert to {-1, +1}

    # Least squares: find w that minimizes ||features @ w - t_pm||^2
    FtF = features.T @ features + regularize * np.eye(features.shape[1])
    Ftt = features.T @ t_pm
    try:
        w = np.linalg.solve(FtF, Ftt)
    except np.linalg.LinAlgError:
        return 0.0, None

    preds = np.sign(features @ w)
    # Handle zeros
    preds[preds == 0] = 1

    accuracy = np.mean(preds == t_pm)
    return accuracy, w

def communication_matrix_svd(pi_mod2, N, k_lsb):
    """Build and SVD the communication matrix for partition (MSBs, LSBs).
    k_lsb: number of least significant bits for Bob.
    Returns singular values."""
    n_rows = 2**(N - k_lsb)  # Alice's part
    n_cols = 2**k_lsb         # Bob's part

    M = np.zeros((n_rows, n_cols), dtype=float)
    for x in range(2**N):
        a = x >> k_lsb
        b = x & ((1 << k_lsb) - 1)
        M[a, b] = 2.0 * pi_mod2[x] - 1.0  # {-1, +1} encoding

    U, S, Vt = svd(M, full_matrices=False)
    return S, U, Vt

def main():
    results = {}

    print("=" * 70)
    print("APPROXIMATE CIRCUIT COMPLEXITY OF pi(x) mod 2")
    print("=" * 70)

    # ========================================
    # PART 1: PTF degree vs accuracy
    # ========================================
    print("\n" + "=" * 70)
    print("PART 1: Polynomial Threshold Function (PTF) accuracy vs degree")
    print("=" * 70)

    for N in [4, 6, 8, 10, 12, 14]:
        print(f"\n--- N = {N} (x in [0, {2**N - 1}]) ---")
        t0 = time.time()
        pi_vals, prime_ind, is_p = build_truth_tables(N)
        pi_mod2 = pi_vals % 2

        # Generate bit representations
        size = 2**N
        x_bits = np.array([[(x >> (N-1-j)) & 1 for j in range(N)] for x in range(size)])

        # Smooth approximation accuracy
        smooth_acc = smooth_approximation_accuracy(N, pi_mod2, pi_vals)
        print(f"  Smooth R(x) accuracy: {smooth_acc:.4f}")

        max_degree = min(N, 4 if N >= 12 else (5 if N >= 10 else N))
        ptf_accs = {}

        for d in range(1, max_degree + 1):
            features = degree_d_features(x_bits, N, d)
            n_features = features.shape[1]
            if n_features > size * 2:
                print(f"  Degree {d}: {n_features} features > 2*{size} samples, skipping")
                break
            acc, w = fit_ptf(features, pi_mod2)
            ptf_accs[d] = acc
            print(f"  Degree {d}: {n_features} features, accuracy = {acc:.6f}" +
                  (" **EXACT**" if acc == 1.0 else ""))
            if acc == 1.0:
                break

        results[N] = {
            'smooth_acc': smooth_acc,
            'ptf_accs': ptf_accs,
            'time': time.time() - t0,
            'pi_mod2': pi_mod2
        }
        print(f"  Time: {time.time() - t0:.1f}s")

    # ========================================
    # PART 2: Error tolerance phase transition
    # ========================================
    print("\n" + "=" * 70)
    print("PART 2: Error tolerance phase transition")
    print("  Q: What fraction of outputs can a degree-1 LTF get right?")
    print("  Reveals: sharp vs gradual complexity transition")
    print("=" * 70)

    for N in [6, 8, 10, 12, 14]:
        if N not in results:
            continue
        pi_mod2 = results[N].get('pi_mod2')
        if pi_mod2 is None:
            continue

        size = 2**N
        x_bits = np.array([[(x >> (N-1-j)) & 1 for j in range(N)] for x in range(size)])

        # Degree 1 (single LTF): what accuracy?
        features_d1 = degree_d_features(x_bits, N, 1)
        acc_d1, _ = fit_ptf(features_d1, pi_mod2)

        # Degree 2: what accuracy?
        features_d2 = degree_d_features(x_bits, N, 2)
        acc_d2, _ = fit_ptf(features_d2, pi_mod2)

        # Constant predictor (majority class)
        majority_acc = max(np.mean(pi_mod2), 1 - np.mean(pi_mod2))

        # Random baseline
        random_acc = 0.5

        print(f"\n  N = {N}:")
        print(f"    Constant predictor: {majority_acc:.4f}")
        print(f"    Degree-1 LTF:      {acc_d1:.4f}")
        print(f"    Degree-2 PTF:      {acc_d2:.4f}")

        # How much does the smooth approximation explain?
        smooth_acc = results[N]['smooth_acc']
        print(f"    Smooth R(x):       {smooth_acc:.4f}")

        # Gap analysis
        print(f"    Gap (d1 vs exact): {1.0 - acc_d1:.4f}")
        print(f"    Gap (d2 vs exact): {1.0 - acc_d2:.4f}")
        print(f"    Gap (smooth vs exact): {1.0 - smooth_acc:.4f}")

    # ========================================
    # PART 3: SVD of communication matrix for unbalanced partitions
    # ========================================
    print("\n" + "=" * 70)
    print("PART 3: SVD of communication matrix (unbalanced partitions)")
    print("  Testing: rank vs partition balance, SV decay rate")
    print("=" * 70)

    for N in [8, 10, 12, 14]:
        if N not in results:
            continue
        pi_mod2 = results[N].get('pi_mod2')
        if pi_mod2 is None:
            continue

        print(f"\n--- N = {N} ---")
        print(f"  {'k_LSB':>5} {'rank_formula':>12} {'rank_numeric':>12} {'SV_ratio':>10} {'top_SV':>10}")

        for k in range(1, N):
            if k > N - 1:
                break
            S, U, Vt = communication_matrix_svd(pi_mod2, N, k)

            # Numerical rank (SVs > 1e-10 * max SV)
            threshold = 1e-10 * S[0] if S[0] > 0 else 1e-10
            num_rank = np.sum(S > threshold)

            # Formula rank
            formula_rank = 2**(min(k, N-k) - 1) + 2

            # SV decay: ratio of 2nd to 1st SV
            sv_ratio = S[1] / S[0] if len(S) > 1 and S[0] > 0 else 0

            print(f"  {k:>5} {formula_rank:>12} {num_rank:>12} {sv_ratio:>10.4f} {S[0]:>10.2f}")

        # Deep dive on k = 2*log2(N) partition
        k_target = max(2, min(int(2 * math.log2(N)), N-1))
        print(f"\n  Deep dive: k_LSB = {k_target} (2*log2(N) = {2*math.log2(N):.1f})")
        S, U, Vt = communication_matrix_svd(pi_mod2, N, k_target)

        # SV spectrum
        print(f"  Top 10 singular values:")
        for i in range(min(10, len(S))):
            cumvar = np.sum(S[:i+1]**2) / np.sum(S**2) if np.sum(S**2) > 0 else 0
            print(f"    SV[{i}] = {S[i]:.4f}  (cumulative variance: {cumvar:.4f})")

        # How many SVs needed for 99%, 99.9%, 100% of variance?
        total_var = np.sum(S**2)
        if total_var > 0:
            cumvar = np.cumsum(S**2) / total_var
            for target in [0.90, 0.95, 0.99, 0.999, 1.0 - 1e-10]:
                k_needed = np.searchsorted(cumvar, target) + 1
                print(f"    SVs for {target*100:.1f}% variance: {k_needed}")

    # ========================================
    # PART 4: Smooth vs oscillatory decomposition
    # ========================================
    print("\n" + "=" * 70)
    print("PART 4: Smooth vs oscillatory decomposition")
    print("  The top SVs should correspond to the smooth part R(x)")
    print("  Remaining SVs = oscillatory part (zeta zeros)")
    print("=" * 70)

    for N in [10, 12, 14]:
        if N not in results:
            continue
        pi_mod2 = results[N].get('pi_mod2')
        if pi_mod2 is None:
            continue

        print(f"\n--- N = {N} ---")

        # Balanced partition
        k = N // 2
        S, U, Vt = communication_matrix_svd(pi_mod2, N, k)

        total_var = np.sum(S**2)

        # Top 2 SVs (smooth part from Session 17: top 2 capture >99.99%)
        smooth_var = np.sum(S[:2]**2) / total_var if total_var > 0 else 0
        osc_var = 1 - smooth_var

        print(f"  Balanced partition (k = {k}):")
        print(f"    Smooth variance (top 2 SVs): {smooth_var:.6f}")
        print(f"    Oscillatory variance:        {osc_var:.6f}")
        print(f"    Total rank: {np.sum(S > 1e-10 * S[0])}")

        # How accurately can the smooth part ALONE predict pi(x) mod 2?
        # Use rank-2 approximation
        M_approx = U[:, :2] @ np.diag(S[:2]) @ Vt[:2, :]
        n_rows = 2**(N - k)
        n_cols = 2**k

        correct_smooth = 0
        total = 2**N
        for x in range(total):
            a = x >> k
            b = x & ((1 << k) - 1)
            pred = 1 if M_approx[a, b] > 0 else 0
            if pred == pi_mod2[x]:
                correct_smooth += 1

        acc_smooth = correct_smooth / total
        print(f"    Rank-2 (smooth) prediction accuracy: {acc_smooth:.4f}")

        # Progressive rank approximation
        print(f"    Rank-k accuracy progression:")
        for r in [1, 2, 3, 5, 10, 20, 50]:
            if r > len(S):
                break
            M_r = U[:, :r] @ np.diag(S[:r]) @ Vt[:r, :]
            correct = 0
            for x in range(total):
                a = x >> (N // 2)
                b = x & ((1 << (N // 2)) - 1)
                pred = 1 if M_r[a, b] > 0 else 0
                if pred == pi_mod2[x]:
                    correct += 1
            print(f"      rank-{r:>3}: {correct/total:.4f}")

    # ========================================
    # PART 5: Error geography — WHERE are the hard inputs?
    # ========================================
    print("\n" + "=" * 70)
    print("PART 5: Error geography — which inputs are hard?")
    print("  If hard inputs cluster, small corrections might work")
    print("=" * 70)

    for N in [10, 12]:
        if N not in results:
            continue
        pi_vals, _, _ = build_truth_tables(N)
        pi_mod2 = pi_vals % 2

        size = 2**N
        x_bits = np.array([[(x >> (N-1-j)) & 1 for j in range(N)] for x in range(size)])

        # Best degree-2 predictor
        features_d2 = degree_d_features(x_bits, N, 2)
        t_pm = 2.0 * pi_mod2 - 1.0
        FtF = features_d2.T @ features_d2 + 1e-6 * np.eye(features_d2.shape[1])
        Ftt = features_d2.T @ t_pm
        w = np.linalg.solve(FtF, Ftt)
        preds = np.sign(features_d2 @ w)
        preds[preds == 0] = 1

        errors = (preds != t_pm)
        error_positions = np.where(errors)[0]

        print(f"\n--- N = {N}: {len(error_positions)} errors out of {size} ({len(error_positions)/size:.4f}) ---")

        # Are errors clustered or spread?
        if len(error_positions) > 1:
            gaps = np.diff(error_positions)
            print(f"  Error gap statistics: mean={gaps.mean():.1f}, std={gaps.std():.1f}, min={gaps.min()}, max={gaps.max()}")

            # Distribution by MSB
            n_msb = N // 2
            msb_vals = error_positions >> (N - n_msb)
            msb_counts = np.bincount(msb_vals, minlength=2**n_msb)
            print(f"  Error distribution by top {n_msb} bits:")
            print(f"    mean={msb_counts.mean():.1f}, std={msb_counts.std():.1f}, max={msb_counts.max()}")
            cv = msb_counts.std() / msb_counts.mean() if msb_counts.mean() > 0 else float('inf')
            print(f"    CV (coeff of variation): {cv:.3f} (1.0 = Poisson-random)")

            # Correlation with pi(x) value
            pi_at_errors = pi_vals[error_positions]
            pi_at_correct = pi_vals[~errors]
            print(f"  Mean pi(x) at errors:  {pi_at_errors.mean():.1f}")
            print(f"  Mean pi(x) at correct: {pi_at_correct.mean():.1f}")

            # Near primes? Fraction of errors at primes vs composites
            frac_prime_errors = np.mean(pi_mod2[error_positions] != pi_mod2[np.clip(error_positions - 1, 0, size-1)])
            print(f"  Fraction of errors where pi(x) changes: {frac_prime_errors:.3f}")

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    print("\n1. PTF degree required for exact computation:")
    for N in sorted(results.keys()):
        accs = results[N]['ptf_accs']
        exact_d = None
        for d in sorted(accs.keys()):
            if accs[d] == 1.0:
                exact_d = d
                break
        best_d = max(accs.keys()) if accs else 0
        best_acc = accs.get(best_d, 0)
        if exact_d:
            print(f"  N={N:>2}: exact at degree {exact_d}")
        else:
            print(f"  N={N:>2}: best degree {best_d} gives {best_acc:.4f}")

    print("\n2. Smooth approximation accuracy:")
    for N in sorted(results.keys()):
        print(f"  N={N:>2}: {results[N]['smooth_acc']:.4f}")

    print("\nDone.")

if __name__ == '__main__':
    main()
