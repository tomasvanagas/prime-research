#!/usr/bin/env python3
"""
Threshold Circuit Depth Profile for Prime Indicator and pi(x) mod 2

Measures the minimum depth of threshold circuits (LTF circuits) needed to
compute is_prime(x) and pi(x) mod 2 on N-bit inputs, for depths 1-5 and
polynomial gate budgets.

For each (N, depth d, gate budget k):
  - Construct depth-d threshold circuits with k gates per layer
  - Use LP to find optimal weights at each layer
  - Report whether exact computation is achieved

This extends S35's depth-2 analysis to depths 3-5, answering:
  "Does depth help? Does poly(N) gates at depth 3-4 suffice?"

Key comparison: if BPSW is in TC^0 (depth 4), then is_prime should be
computable at low depth. We measure this empirically for small N.
"""

import numpy as np
from itertools import product as iproduct
import sys
import time

def sieve_primes(limit):
    """Simple sieve of Eratosthenes."""
    if limit < 2:
        return set()
    is_p = [True] * (limit + 1)
    is_p[0] = is_p[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if is_p[i]:
            for j in range(i*i, limit + 1, i):
                is_p[j] = False
    return set(i for i in range(limit + 1) if is_p[i])

def build_truth_tables(N):
    """Build truth tables for is_prime and pi(x) mod 2 on N-bit inputs."""
    domain = range(2**N)
    primes = sieve_primes(2**N - 1)

    is_prime_tt = np.array([1 if x in primes else 0 for x in domain], dtype=np.int8)

    # pi(x) mod 2: running parity of prime count
    pi_mod2_tt = np.zeros(2**N, dtype=np.int8)
    count = 0
    for x in domain:
        if x in primes:
            count += 1
        pi_mod2_tt[x] = count % 2

    # Random function with same weight as pi_mod2
    rng = np.random.RandomState(42)
    w = int(pi_mod2_tt.sum())
    rand_tt = np.zeros(2**N, dtype=np.int8)
    indices = rng.choice(2**N, size=w, replace=False)
    rand_tt[indices] = 1

    return is_prime_tt, pi_mod2_tt, rand_tt

def build_input_matrix(N):
    """Build the 2^N x N binary input matrix."""
    X = np.zeros((2**N, N), dtype=np.float64)
    for i in range(2**N):
        for b in range(N):
            X[i, b] = (i >> b) & 1
    return X

def find_ltf(features, target):
    """
    Find a linear threshold function (LTF) of the given features that computes target.
    Uses LP: find w such that sign(features @ w) = 2*target - 1.

    Returns (success, weights) where success is True if exact match found.
    """
    n_samples, n_features = features.shape
    labels = 2 * target.astype(np.float64) - 1  # {-1, +1}

    # We want: labels[i] * (features[i] @ w) > 0 for all i
    # Equivalently: (labels * features) @ w > 0
    # LP: maximize min_i { labels[i] * (features[i] @ w) } subject to ||w|| bounded

    try:
        from scipy.optimize import linprog

        # Maximize gamma subject to: labels[i] * (features[i] @ w) >= gamma, ||w||_inf <= 1
        # Rewrite: minimize -gamma
        # Variables: [w_1, ..., w_n_features, gamma]

        c = np.zeros(n_features + 1)
        c[-1] = -1  # minimize -gamma

        # Constraints: labels[i] * features[i] @ w - gamma >= 0
        # i.e.: -(labels[i] * features[i]) @ w + gamma <= 0  (for >= constraint as <=)
        A_ub = np.zeros((n_samples + 2 * n_features, n_features + 1))
        b_ub = np.zeros(n_samples + 2 * n_features)

        for i in range(n_samples):
            A_ub[i, :n_features] = -labels[i] * features[i]
            A_ub[i, -1] = 1  # +gamma on wrong side

        # Bound weights: w_j <= 1 and -w_j <= 1
        for j in range(n_features):
            A_ub[n_samples + 2*j, j] = 1
            b_ub[n_samples + 2*j] = 1
            A_ub[n_samples + 2*j + 1, j] = -1
            b_ub[n_samples + 2*j + 1] = 1

        # Gamma unbounded below
        bounds = [(None, None)] * (n_features + 1)

        result = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method='highs')

        if result.success:
            w = result.x[:n_features]
            gamma = result.x[-1]

            # Verify
            predictions = (features @ w > 0).astype(np.int8)
            accuracy = np.mean(predictions == target)

            if accuracy == 1.0:
                return True, w
            else:
                return False, w
        else:
            return False, None

    except ImportError:
        # Fallback: simple perceptron-like approach
        return False, None

def random_ltf_features(X, k, rng):
    """
    Generate k random LTF features from input matrix X.
    Each feature is sign(X @ w + b) for random w, b.
    """
    N = X.shape[1]
    W = rng.randn(N, k)
    b = rng.randn(k) * 0.5
    features = (X @ W + b > 0).astype(np.float64)
    return features

def measure_depth(N, target_tt, target_name, max_depth=5, gates_per_layer_options=None,
                  n_random_trials=20):
    """
    For a given target function, measure minimum depth needed with poly(N) gates.

    At each depth d, tries increasing gate counts until exact computation is achieved
    or budget is exhausted.
    """
    X = build_input_matrix(N)
    n_points = 2**N

    if gates_per_layer_options is None:
        gates_per_layer_options = [N, 2*N, N**2, 2*N**2, N**3]
        # Keep reasonable for small N
        gates_per_layer_options = [k for k in gates_per_layer_options if k <= 500]

    results = {}
    rng = np.random.RandomState(123)

    for depth in range(1, max_depth + 1):
        best_accuracy = 0.0
        best_k = None
        exact = False

        for k in gates_per_layer_options:
            if exact:
                break

            trial_accuracies = []

            for trial in range(n_random_trials):
                # Build depth-d circuit
                # Layer 1: random LTFs of inputs
                current_features = X.copy()

                for layer in range(depth - 1):
                    n_in = current_features.shape[1]
                    if layer == 0:
                        layer_features = random_ltf_features(current_features, k, rng)
                    else:
                        layer_features = random_ltf_features(current_features, min(k, current_features.shape[1] * 2), rng)
                    current_features = layer_features

                # Final layer: single LTF via LP
                success, w = find_ltf(current_features, target_tt)

                if success:
                    trial_accuracies.append(1.0)
                    exact = True
                    best_k = k
                    break
                else:
                    if w is not None:
                        preds = (current_features @ w > 0).astype(np.int8)
                        acc = np.mean(preds == target_tt)
                    else:
                        acc = 0.5
                    trial_accuracies.append(acc)

            if trial_accuracies:
                max_acc = max(trial_accuracies)
                if max_acc > best_accuracy:
                    best_accuracy = max_acc
                    if not exact:
                        best_k = k

        results[depth] = {
            'exact': exact,
            'best_accuracy': best_accuracy if not exact else 1.0,
            'gates_per_layer': best_k,
            'total_gates': best_k * depth if best_k else None
        }

        status = "EXACT" if exact else f"{best_accuracy:.4f}"
        print(f"  depth={depth}, k={best_k}: {status}")

        if exact:
            # No need to try deeper
            break

    return results

def main():
    print("=" * 70)
    print("Threshold Circuit Depth Profile for Prime-Related Functions")
    print("Measuring minimum depth with poly(N) gates per layer")
    print("=" * 70)

    # Test range
    N_values = [4, 6, 8, 10]

    # For larger N, reduce trials
    trial_counts = {4: 50, 6: 30, 8: 20, 10: 10}

    all_results = {}

    for N in N_values:
        print(f"\n{'='*60}")
        print(f"N = {N}  (domain size = {2**N})")
        print(f"{'='*60}")

        is_prime_tt, pi_mod2_tt, rand_tt = build_truth_tables(N)

        print(f"  is_prime weight: {is_prime_tt.sum()}/{2**N}")
        print(f"  pi_mod2 weight: {pi_mod2_tt.sum()}/{2**N}")
        print(f"  random weight: {rand_tt.sum()}/{2**N}")

        gates = [N, 2*N, N**2, 2*N**2]
        gates = [k for k in gates if k <= 300]

        n_trials = trial_counts.get(N, 10)

        for name, tt in [("is_prime", is_prime_tt), ("pi_mod2", pi_mod2_tt), ("random", rand_tt)]:
            print(f"\n  --- {name} (N={N}) ---")
            t0 = time.time()
            res = measure_depth(N, tt, name, max_depth=5,
                              gates_per_layer_options=gates,
                              n_random_trials=n_trials)
            elapsed = time.time() - t0
            print(f"  ({elapsed:.1f}s)")
            all_results[(N, name)] = res

    # Summary table
    print("\n" + "=" * 70)
    print("SUMMARY: Minimum depth achieving exact computation")
    print("(with poly(N) gates per layer, heuristic search)")
    print("=" * 70)

    print(f"\n{'N':>4} | {'is_prime':>20} | {'pi_mod2':>20} | {'random':>20}")
    print("-" * 70)

    for N in N_values:
        row = f"{N:>4} |"
        for name in ["is_prime", "pi_mod2", "random"]:
            res = all_results.get((N, name), {})
            min_depth = None
            best_acc = 0.0
            for d in sorted(res.keys()):
                if res[d]['exact']:
                    min_depth = d
                    break
                best_acc = max(best_acc, res[d]['best_accuracy'])

            if min_depth is not None:
                cell = f"depth={min_depth} (k={res[min_depth]['gates_per_layer']})"
            else:
                cell = f">{max(res.keys()) if res else '?'} (best={best_acc:.3f})"
            row += f" {cell:>20} |"
        print(row)

    print("\n" + "=" * 70)
    print("DEPTH SCALING ANALYSIS")
    print("=" * 70)

    for name in ["is_prime", "pi_mod2", "random"]:
        print(f"\n  {name}:")
        for N in N_values:
            res = all_results.get((N, name), {})
            depths_str = []
            for d in sorted(res.keys()):
                r = res[d]
                if r['exact']:
                    depths_str.append(f"d={d}:EXACT(k={r['gates_per_layer']})")
                else:
                    depths_str.append(f"d={d}:{r['best_accuracy']:.3f}")
            print(f"    N={N}: {', '.join(depths_str)}")

    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    print("""
If is_prime achieves exact computation at constant depth with poly(N) gates
as N grows, this provides empirical evidence for PRIMES in TC^0.

If depth must grow (e.g., log N), this suggests NC^1 \ TC^0 behavior.

If is_prime behaves identically to random functions, this extends the
pseudorandomness evidence from S35 to the depth dimension.

CAVEATS:
- Heuristic search (random features + LP) gives UPPER bounds on depth
- True minimum depth could be lower with optimal feature selection
- Small N (4-10) may not reflect asymptotic behavior
- LP may fail to find separating hyperplanes that exist
""")

if __name__ == "__main__":
    main()
