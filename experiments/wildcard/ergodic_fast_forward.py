#!/usr/bin/env python3
"""
Ergodic Fast-Forwarding: Can prime gaps be modeled as a fast-forwardable dynamical system?

Analogy: Linear recurrences like Fibonacci can be computed in O(log n) via matrix
exponentiation. Can we find a dynamical system for prime gaps with similar structure?

The sequence g(n) = p(n+1) - p(n) has been studied extensively.
Cramér's model: g(n) ~ Exponential(ln(p(n))) for "random" primes.
But actual gaps have more structure (divisibility constraints, etc.)

Key question: Is there a MAP T such that:
  (p(n+1), state) = T(p(n), state)
where state is finite-dimensional and T can be iterated via repeated squaring?

If state dimension is O(polylog), then T^n gives p(n) in O(polylog) time.

This experiment tests:
1. Is the gap sequence g(n) a function of finitely many previous gaps?
2. Can gaps be predicted from a finite-dimensional hidden state?
3. Does the system have a Lyapunov exponent (chaos test)?
4. Can a recurrence relation be fit to g(n)?
"""

import numpy as np
from sympy import primerange, nextprime
import time
import math

def compute_gaps(N):
    """Compute first N prime gaps."""
    primes = list(primerange(2, 2 * N * int(np.log(N + 10) + 5)))[:N+1]
    gaps = [primes[i+1] - primes[i] for i in range(len(primes)-1)]
    return gaps, primes

def test_markov_property(gaps, order=1):
    """
    Test if g(n+1) is determined by (g(n), g(n-1), ..., g(n-order+1)).
    Measure: conditional entropy H(g(n+1) | g(n),...,g(n-order+1)).
    """
    from collections import Counter

    # Build transition table
    transitions = {}
    for i in range(order, len(gaps)):
        state = tuple(gaps[i-order:i])
        next_gap = gaps[i]
        if state not in transitions:
            transitions[state] = []
        transitions[state].append(next_gap)

    # Compute conditional entropy
    total_entropy = 0
    total_count = 0
    deterministic = 0
    total_states = 0

    for state, nexts in transitions.items():
        counts = Counter(nexts)
        n = len(nexts)
        total_states += 1

        if len(counts) == 1:
            deterministic += 1

        h = 0
        for c in counts.values():
            p = c / n
            if p > 0:
                h -= p * np.log2(p)
        total_entropy += h * n
        total_count += n

    cond_entropy = total_entropy / total_count if total_count > 0 else 0

    # Unconditional entropy of gaps
    gap_counts = Counter(gaps)
    unc_entropy = 0
    for c in gap_counts.values():
        p = c / len(gaps)
        if p > 0:
            unc_entropy -= p * np.log2(p)

    return {
        'order': order,
        'cond_entropy': cond_entropy,
        'unconditional_entropy': unc_entropy,
        'entropy_reduction': 1 - cond_entropy / unc_entropy if unc_entropy > 0 else 0,
        'unique_states': total_states,
        'deterministic_states': deterministic,
        'det_fraction': deterministic / total_states if total_states > 0 else 0,
    }

def test_recurrence_fit(gaps, max_order=10):
    """
    Fit g(n) = a1*g(n-1) + a2*g(n-2) + ... + ak*g(n-k) + c
    and measure R^2.
    """
    results = []
    gaps_arr = np.array(gaps, dtype=float)

    for order in range(1, max_order + 1):
        # Build design matrix
        X = np.zeros((len(gaps) - order, order))
        y = gaps_arr[order:]

        for k in range(order):
            X[:, k] = gaps_arr[order - k - 1: len(gaps) - k - 1]

        # Add constant
        X = np.column_stack([X, np.ones(len(y))])

        # Least squares
        try:
            coeffs, residuals, _, _ = np.linalg.lstsq(X, y, rcond=None)
            y_pred = X @ coeffs
            ss_res = np.sum((y - y_pred)**2)
            ss_tot = np.sum((y - y.mean())**2)
            r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        except:
            r_squared = 0
            coeffs = np.zeros(order + 1)

        results.append({
            'order': order,
            'r_squared': r_squared,
            'coeffs': coeffs[:order],
            'constant': coeffs[-1] if len(coeffs) > order else 0,
        })

    return results

def test_embedding_dimension(gaps, max_dim=15):
    """
    False nearest neighbors test for embedding dimension.
    If gaps come from a low-dimensional dynamical system,
    the false neighbor fraction drops to ~0 at the true dimension.
    """
    gaps_arr = np.array(gaps, dtype=float)
    n = len(gaps_arr)

    results = []
    for dim in range(1, max_dim + 1):
        # Build delay embedding
        if dim > n - 1:
            break

        embedded = np.zeros((n - dim, dim))
        for k in range(dim):
            embedded[:, k] = gaps_arr[k:n - dim + k]

        # For each point, find nearest neighbor
        # (subsample for speed)
        sample_size = min(500, len(embedded))
        indices = np.random.choice(len(embedded), sample_size, replace=False)

        false_nn = 0
        total = 0

        for idx in indices:
            point = embedded[idx]
            dists = np.linalg.norm(embedded - point, axis=1)
            dists[idx] = np.inf
            nn_idx = np.argmin(dists)
            nn_dist = dists[nn_idx]

            if nn_dist < 1e-10:
                continue

            # Check if neighbor is still close in next dimension
            if dim < len(gaps_arr) and idx + dim < len(gaps_arr) and nn_idx + dim < len(gaps_arr):
                extra_dist = abs(gaps_arr[idx + dim] - gaps_arr[nn_idx + dim])
                if extra_dist / nn_dist > 10:  # Threshold for "false"
                    false_nn += 1
            total += 1

        fnn_fraction = false_nn / total if total > 0 else 0
        results.append({
            'dim': dim,
            'fnn_fraction': fnn_fraction,
        })

    return results

def test_lyapunov_exponent(gaps, embedding_dim=3):
    """
    Estimate largest Lyapunov exponent of the gap sequence.
    Positive → chaotic → not fast-forwardable.
    Zero → quasi-periodic → potentially fast-forwardable.
    Negative → contracting → trivially fast-forwardable (but wrong model).
    """
    gaps_arr = np.array(gaps, dtype=float)
    n = len(gaps_arr) - embedding_dim

    # Build delay embedding
    embedded = np.zeros((n, embedding_dim))
    for k in range(embedding_dim):
        embedded[:, k] = gaps_arr[k:n + k]

    # Estimate Lyapunov exponent via nearest neighbor divergence
    lyap_estimates = []

    sample_size = min(200, n - 10)
    indices = np.random.choice(n - 10, sample_size, replace=False)

    for idx in indices:
        point = embedded[idx]
        dists = np.linalg.norm(embedded[:n-10] - point, axis=1)
        dists[idx] = np.inf

        nn_idx = np.argmin(dists)
        initial_dist = dists[nn_idx]

        if initial_dist < 1e-10:
            continue

        # Track divergence over several steps
        for steps in [1, 2, 5, 10]:
            if idx + steps < n and nn_idx + steps < n:
                new_dist = np.linalg.norm(embedded[idx + steps] - embedded[nn_idx + steps])
                if new_dist > 0 and initial_dist > 0:
                    lyap = np.log(new_dist / initial_dist) / steps
                    lyap_estimates.append((steps, lyap))

    # Average by step count
    from collections import defaultdict
    by_steps = defaultdict(list)
    for steps, lyap in lyap_estimates:
        by_steps[steps].append(lyap)

    return {steps: np.mean(vals) for steps, vals in sorted(by_steps.items())}

def main():
    print("=" * 70)
    print("ERGODIC FAST-FORWARDING OF PRIME GAPS")
    print("=" * 70)

    N = 20000
    print(f"\nComputing first {N} prime gaps...")
    t0 = time.time()
    gaps, primes = compute_gaps(N)
    print(f"  Done in {time.time()-t0:.1f}s")
    print(f"  Range: p(1)={primes[0]} to p({N})={primes[-1]}")
    print(f"  Gap statistics: mean={np.mean(gaps):.2f}, std={np.std(gaps):.2f}, "
          f"max={max(gaps)}")

    # Test 1: Markov property
    print("\n--- Test 1: Conditional Entropy (Markov Order) ---")
    for order in [1, 2, 3, 4, 5, 8, 10]:
        result = test_markov_property(gaps, order)
        print(f"  order={order:2d}: H(g|history)={result['cond_entropy']:.3f} "
              f"(uncond: {result['unconditional_entropy']:.3f}, "
              f"reduction: {result['entropy_reduction']*100:.1f}%, "
              f"det_states: {result['det_fraction']*100:.1f}%)")

    # Test 2: Linear recurrence fit
    print("\n--- Test 2: Linear Recurrence R^2 ---")
    rec_results = test_recurrence_fit(gaps)
    for r in rec_results:
        print(f"  order={r['order']:2d}: R^2 = {r['r_squared']:.6f}")

    # Test 3: Embedding dimension
    print("\n--- Test 3: False Nearest Neighbors (Embedding Dimension) ---")
    np.random.seed(42)
    fnn_results = test_embedding_dimension(gaps)
    for r in fnn_results:
        bar = '#' * int(r['fnn_fraction'] * 50)
        print(f"  dim={r['dim']:2d}: FNN = {r['fnn_fraction']:.3f} {bar}")

    # Test 4: Lyapunov exponent
    print("\n--- Test 4: Lyapunov Exponent ---")
    np.random.seed(42)
    lyap = test_lyapunov_exponent(gaps)
    for steps, val in lyap.items():
        sign = "CHAOTIC" if val > 0.1 else "NEUTRAL" if val > -0.1 else "CONTRACTING"
        print(f"  steps={steps:2d}: lambda = {val:.4f} ({sign})")

    # Test 5: Compare with random model
    print("\n--- Test 5: Comparison with Cramér Random Model ---")
    # Generate Cramér-model gaps: Exp(ln(p)) for each p
    cramer_gaps = []
    np.random.seed(42)
    for p in primes[:-1]:
        lp = np.log(p)
        cramer_gaps.append(max(1, int(np.round(np.random.exponential(lp)))))

    # Compare entropies
    real_markov = test_markov_property(gaps, 3)
    cramer_markov = test_markov_property(cramer_gaps[:len(gaps)], 3)
    print(f"  Real gaps:   H(g|3-history) = {real_markov['cond_entropy']:.3f}")
    print(f"  Cramér gaps: H(g|3-history) = {cramer_markov['cond_entropy']:.3f}")
    print(f"  Difference: {abs(real_markov['cond_entropy'] - cramer_markov['cond_entropy']):.3f}")

    rec_real = test_recurrence_fit(gaps, max_order=5)
    rec_cramer = test_recurrence_fit(cramer_gaps[:len(gaps)], max_order=5)
    for r, c in zip(rec_real, rec_cramer):
        print(f"  order {r['order']}: R^2 real={r['r_squared']:.6f}, "
              f"Cramér={c['r_squared']:.6f}")

    print("\n--- VERDICT ---")
    print("For fast-forwarding to work, we need:")
    print("  1. Low conditional entropy (predictable) → FAILS if entropy stays high")
    print("  2. Good linear recurrence fit (R^2 > 0.5) → FAILS if R^2 ≈ 0")
    print("  3. Low embedding dimension → FAILS if FNN never reaches 0")
    print("  4. Non-positive Lyapunov exponent → FAILS if chaotic")
    print()
    print("If gaps behave like Cramér random model, no fast-forwarding is possible")
    print("because random sequences have no exploitable dynamical structure.")

if __name__ == "__main__":
    main()
