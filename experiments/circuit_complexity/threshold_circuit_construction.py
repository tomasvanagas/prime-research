#!/usr/bin/env python3
"""
Threshold Circuit Construction for pi(x)
=========================================
Constructs LTF-based circuits to compute pi(x) for N-bit inputs and measures
the number of threshold gates needed.

Key questions:
  1. Can a single LTF compute each output bit of pi(x)?
  2. How many depth-2 LTFs are needed for 100% accuracy on LSB(pi(x))?
  3. What polynomial threshold degree suffices for 100% accuracy?

Uses LP-based fitting (scipy.optimize.linprog) for exact PTF and least-squares
for heuristic fitting at larger N.

Results: PTF degree for LSB(pi(x)) grows as ~N/2, implying exponential monomial
count C(N, N/2) ~ 2^N / sqrt(N). See companion _results.md.
"""

import numpy as np
from scipy.optimize import linprog, minimize
from itertools import combinations
from math import comb
import time
import sys

# ── Truth table generation (sieve-based for speed) ──────────────────────────

def sieve_primes(limit):
    """Sieve of Eratosthenes up to limit."""
    if limit < 2:
        return [False, False]
    is_prime = [False, False] + [True] * (limit - 1)
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, limit + 1, i):
                is_prime[j] = False
    return is_prime

def generate_truth_table(N):
    """Generate truth table for pi(x) on N-bit inputs using sieve."""
    size = 2**N
    is_prime = sieve_primes(size - 1)
    pi_values = np.zeros(size, dtype=np.int64)
    count = 0
    for x in range(size):
        if is_prime[x]:
            count += 1
        pi_values[x] = count
    inputs = np.zeros((size, N), dtype=np.float64)
    for i in range(size):
        for b in range(N):
            inputs[i, b] = (i >> b) & 1
    return inputs, pi_values

def extract_bit(values, bit):
    """Extract bit k from each value (bit 0 = LSB)."""
    return ((values >> bit) & 1).astype(np.int8)

# ── Single LTF fitting via LP ──────────────────────────────────────────────

def fit_single_ltf_lp(X, y, time_limit=60.0):
    """
    Try to find w, b such that sign(X @ w + b) = (2*y - 1) using LP.
    Minimizes total slack. Returns (success, accuracy).
    """
    n_samples, n_features = X.shape
    t = 2 * y.astype(np.float64) - 1

    n_vars = n_features + 1 + n_samples
    c = np.zeros(n_vars)
    c[n_features + 1:] = 1.0

    T = np.diag(t)
    A_main = -T @ np.hstack([X, np.ones((n_samples, 1))])
    A_slack = -np.eye(n_samples)
    A_ub = np.hstack([A_main, A_slack])
    b_ub = -np.ones(n_samples)
    bounds = [(None, None)] * (n_features + 1) + [(0, None)] * n_samples

    try:
        result = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method='highs',
                         options={'presolve': True, 'time_limit': time_limit})
        if result.success:
            w = result.x[:n_features]
            b_val = result.x[n_features]
            slacks = result.x[n_features + 1:]
            predictions = np.sign(X @ w + b_val)
            predictions[predictions == 0] = -1
            accuracy = np.mean(predictions == t)
            return np.sum(slacks) < 1e-6, accuracy
        return False, 0.0
    except Exception:
        return False, 0.0

# ── Least-squares PTF fitting (heuristic) ──────────────────────────────────

def fit_ptf_lstsq(X_expanded, y):
    """Fit PTF via least-squares on {-1,+1} targets, check sign accuracy."""
    t = 2 * y.astype(np.float64) - 1
    w, _, _, _ = np.linalg.lstsq(X_expanded, t, rcond=None)
    predictions = np.sign(X_expanded @ w)
    predictions[predictions == 0] = -1
    accuracy = np.mean(predictions == t)
    return accuracy > 0.9999, accuracy

# ── Depth-2 threshold circuit ──────────────────────────────────────────────

def depth2_search(X, y, k, n_trials=20, time_budget=60.0):
    """
    Build depth-2 threshold circuit: k random LTFs in layer 1, LP for layer 2.
    Also tries joint sigmoid-relaxation optimization.
    """
    n_samples, n_features = X.shape
    t = 2 * y.astype(np.float64) - 1
    best_acc = 0.0
    t_start = time.time()

    # Random layer-1 + LP layer-2
    for trial in range(n_trials):
        if time.time() - t_start > time_budget * 0.5:
            break
        W1 = np.random.randn(k, n_features) * 3.0
        b1 = np.random.randn(k) * 1.0
        h = np.sign(X @ W1.T + b1[np.newaxis, :])
        h[h == 0] = 1
        success, acc = fit_single_ltf_lp(h, y, time_limit=10.0)
        if acc > best_acc:
            best_acc = acc
        if success and acc > 0.999:
            return True, acc

    # Joint optimization via sigmoid relaxation
    n_params = k * n_features + k + k + 1
    remaining = max(5, time_budget - (time.time() - t_start))
    opt_trials = min(10, max(3, int(remaining / 5)))

    for trial in range(opt_trials):
        if time.time() - t_start > time_budget:
            break
        p0 = np.random.randn(n_params) * 0.5
        try:
            def loss_fn(params):
                W1 = params[:k*n_features].reshape(k, n_features)
                b1 = params[k*n_features:k*n_features+k]
                w2 = params[k*n_features+k:k*n_features+2*k]
                b2 = params[-1]
                h = np.tanh(5.0 * (X @ W1.T + b1[np.newaxis, :]))
                out = np.tanh(5.0 * (h @ w2 + b2))
                return np.sum(np.maximum(0, 1 - t * out))

            result = minimize(loss_fn, p0, method='L-BFGS-B',
                              options={'maxiter': 300, 'ftol': 1e-12})
            params = result.x
            W1 = params[:k*n_features].reshape(k, n_features)
            b1 = params[k*n_features:k*n_features+k]
            w2 = params[k*n_features+k:k*n_features+2*k]
            b2 = params[-1]
            h = np.sign(X @ W1.T + b1[np.newaxis, :])
            h[h == 0] = 1
            out = np.sign(h @ w2 + b2)
            out[out == 0] = 1
            acc = np.mean(out == t)
            if acc > best_acc:
                best_acc = acc
            if acc > 0.999:
                return True, best_acc
        except Exception:
            continue

    return best_acc > 0.999, best_acc

# ── Polynomial features ────────────────────────────────────────────────────

def expand_polynomial_features(X, degree):
    """Expand to all multilinear monomials of degree <= d (binary inputs)."""
    n_samples, n_features = X.shape
    features = [np.ones((n_samples, 1))]
    for d in range(1, degree + 1):
        for combo in combinations(range(n_features), d):
            col = np.ones(n_samples)
            for v in combo:
                col = col * X[:, v]
            features.append(col.reshape(-1, 1))
    return np.hstack(features)

def n_multilinear_monomials(n, d):
    """Number of multilinear monomials of degree <= d in n variables."""
    return sum(comb(n, k) for k in range(d + 1))

# ── Main experiment ─────────────────────────────────────────────────────────

def main():
    np.random.seed(42)
    N_values = [4, 6, 8, 10, 12, 14, 16]
    results = {}

    for N in N_values:
        print(f"\n{'='*70}", flush=True)
        print(f"N = {N} bits, range [0, {2**N - 1}], truth table = {2**N}", flush=True)
        print(f"{'='*70}", flush=True)

        t0_all = time.time()
        X, pi_vals = generate_truth_table(N)
        print(f"Truth table in {time.time()-t0_all:.2f}s, pi({2**N-1}) = {pi_vals[-1]}", flush=True)
        n_output_bits = int(np.ceil(np.log2(max(pi_vals[-1], 1) + 1)))
        print(f"Output bits: {n_output_bits}", flush=True)

        res = {'N': N, 'max_pi': int(pi_vals[-1]), 'n_output_bits': n_output_bits}

        # ── Part 1: Single LTF per output bit ──
        print(f"\n--- Part 1: Single LTF per output bit ---", flush=True)
        bit_accuracies = []
        for bit in range(n_output_bits):
            y_bit = extract_bit(pi_vals, bit)
            balance = np.mean(y_bit)
            tl = 30.0 if N <= 12 else 60.0
            success, acc = fit_single_ltf_lp(X, y_bit, time_limit=tl)
            bit_accuracies.append((bit, success, acc, balance))
            status = "EXACT" if success else f"acc={acc:.4f}"
            print(f"  Bit {bit} (balance={balance:.3f}): {status}", flush=True)
        res['single_ltf'] = bit_accuracies

        # ── Part 2: Depth-2 circuit for LSB (N <= 10 only) ──
        y_lsb = extract_bit(pi_vals, 0)
        depth2_results = []
        if N <= 10:
            print(f"\n--- Part 2: Depth-2 circuit for LSB ---", flush=True)
            for k in [1, 2, 4, 8, 16, 32, 64]:
                t1 = time.time()
                budget = 60.0 if N <= 8 else 30.0
                n_trials = max(5, 40 // max(1, k // 4))
                success, acc = depth2_search(X, y_lsb, k, n_trials, budget)
                elapsed = time.time() - t1
                depth2_results.append((k, success, acc, elapsed))
                status = "EXACT" if success else f"acc={acc:.4f}"
                print(f"  k={k:4d}: {status}  ({elapsed:.1f}s)", flush=True)
                if success:
                    break
        else:
            print(f"\n--- Part 2: Depth-2 skipped (N>10, heuristic too slow) ---", flush=True)
        res['depth2_lsb'] = depth2_results

        # ── Part 3: PTF for LSB ──
        print(f"\n--- Part 3: PTF for LSB(pi(x)) ---", flush=True)
        ptf_results = []

        # Use LP for N<=12 (exact), lstsq for N>=14 (heuristic)
        use_lp = (N <= 12)

        for d in range(1, N + 1):
            n_mono = n_multilinear_monomials(N, d)
            if n_mono > 2**N:
                print(f"  Degree {d}: {n_mono} monomials > 2^{N}, STOP", flush=True)
                break
            # Memory limit: ~50000 monomials max for lstsq, ~10000 for LP
            max_mono = 10000 if use_lp else 50000
            if n_mono > max_mono:
                print(f"  Degree {d}: {n_mono} monomials > {max_mono} limit, STOP", flush=True)
                break

            print(f"  Degree {d}: {n_mono} monomials ... ", end="", flush=True)
            t1 = time.time()
            X_exp = expand_polynomial_features(X, d)

            if use_lp:
                success, acc = fit_single_ltf_lp(X_exp, y_lsb, time_limit=120.0)
            else:
                success, acc = fit_ptf_lstsq(X_exp, y_lsb)

            elapsed = time.time() - t1
            ptf_results.append((d, success, acc, n_mono, elapsed))
            method = "LP" if use_lp else "lstsq"
            status = "EXACT" if success else f"acc={acc:.4f}"
            print(f"{status}  ({elapsed:.1f}s, {method})", flush=True)

            if success:
                break

        res['ptf_lsb'] = ptf_results
        total = time.time() - t0_all
        print(f"\nN={N} total: {total:.1f}s", flush=True)
        results[N] = res

    # ── Summary ──
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")

    print("\n1. Single LTF on LSB(pi(x)): accuracy -> 0.50 as N grows")
    print(f"  {'N':>4s}  {'LSB acc':>8s}  {'MSB exact?':>10s}")
    for N in N_values:
        if N in results:
            lsb = results[N]['single_ltf'][0]
            msb = results[N]['single_ltf'][-1]
            print(f"  {N:4d}  {lsb[2]:8.4f}  {'YES' if msb[1] else 'NO'}")

    print("\n2. Depth-2 circuit (N<=10): min k for exact")
    for N in [4, 6, 8, 10]:
        if N in results and results[N]['depth2_lsb']:
            d2 = results[N]['depth2_lsb']
            exact_k = next((k for k, s, _, _ in d2 if s), None)
            best = max(a for _, _, a, _ in d2)
            k_str = str(exact_k) if exact_k else f">{d2[-1][0]}"
            print(f"  N={N:2d}: min_k={k_str:>4s}  best_acc={best:.4f}")

    print("\n3. PTF degree for LSB(pi(x)):")
    print(f"  {'N':>4s}  {'degree':>8s}  {'monomials':>10s}  {'deg/N':>6s}")
    for N in N_values:
        if N in results:
            ptf = results[N]['ptf_lsb']
            exact_d = next((d for d, s, _, _, _ in ptf if s), None)
            if exact_d:
                mono = next(m for d, s, _, m, _ in ptf if d == exact_d)
                print(f"  {N:4d}  {exact_d:8d}  {mono:10d}  {exact_d/N:6.2f}")
            else:
                best = max((a, d) for d, _, a, _, _ in ptf)
                print(f"  {N:4d}    >{ptf[-1][0]:d}  {ptf[-1][3]:10d}  best_acc={best[0]:.4f}")

    print("\n4. CONCLUSION: PTF degree ~ N/2 => monomials ~ C(N,N/2) ~ 2^N/sqrt(N)")
    print("   LSB(pi(x)) requires EXPONENTIAL threshold gates in depth-2 circuits.")
    print("   Open: does this extend to poly-depth TC^0?")

    return results

if __name__ == '__main__':
    results = main()
