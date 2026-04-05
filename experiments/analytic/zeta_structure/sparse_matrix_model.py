#!/usr/bin/env python3
"""
Sparse Matrix Model for Riemann Zeta Zeros
===========================================
Test whether the zeta zero sequence can be modeled as eigenvalues of a SPARSE
Hermitian matrix with O(N) entries instead of O(N^2).

Approaches:
  A) Tridiagonal (Jacobi) matrix -- 2N-1 parameters
  B) Banded matrix with bandwidth k -- O(kN) parameters
  C) Structured sparse: circulant+diagonal, Toeplitz, prime-based

Key question: minimum bandwidth for <1% per-eigenvalue error, and whether
that bandwidth grows with N.
"""

import numpy as np
from scipy.optimize import minimize
from scipy.linalg import eigvalsh, toeplitz
from scipy.sparse import diags
import time
import os
import sys

DATA_PATH = "/apps/aplikacijos/prime-research/data/zeta_zeros_1000.txt"
RESULTS_PATH = "/apps/aplikacijos/prime-research/experiments/analytic/zeta_structure/sparse_matrix_results.md"

# ── Load zeros ──────────────────────────────────────────────────────────────
def load_zeros(path, n=None):
    zeros = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                zeros.append(float(line))
    zeros = np.array(zeros)
    if n is not None:
        zeros = zeros[:n]
    return zeros

ALL_ZEROS = load_zeros(DATA_PATH)
print(f"Loaded {len(ALL_ZEROS)} zeros. Range: [{ALL_ZEROS[0]:.4f}, {ALL_ZEROS[-1]:.4f}]")

# ── Helpers ─────────────────────────────────────────────────────────────────
def relative_error(eigs, targets):
    """Mean relative error per eigenvalue."""
    return np.mean(np.abs(eigs - targets) / np.abs(targets))

def rms_error(eigs, targets):
    """RMS absolute error."""
    return np.sqrt(np.mean((eigs - targets)**2))

def max_relative_error(eigs, targets):
    return np.max(np.abs(eigs - targets) / np.abs(targets))

# ── Approach A: Tridiagonal (Jacobi) matrix ─────────────────────────────────
def jacobi_eigenvalues(params, N):
    """Build symmetric tridiagonal matrix and return sorted eigenvalues."""
    a = params[:N]       # diagonal
    b = params[N:2*N-1]  # off-diagonal
    M = np.diag(a) + np.diag(b, 1) + np.diag(b, -1)
    return np.sort(np.linalg.eigvalsh(M))

def jacobi_loss(params, targets):
    N = len(targets)
    eigs = jacobi_eigenvalues(params, N)
    # Weight by 1/gamma_i^2 so relative error matters
    return np.sum(((eigs - targets) / targets)**2)

def fit_jacobi(targets, maxiter=2000):
    N = len(targets)
    # Initialize: diagonal ~ targets, off-diagonal ~ mean spacing / 2
    spacing = np.diff(targets)
    a0 = targets.copy()
    b0 = np.concatenate([spacing / 2.0])
    x0 = np.concatenate([a0, b0])

    res = minimize(jacobi_loss, x0, args=(targets,),
                   method='L-BFGS-B', options={'maxiter': maxiter, 'ftol': 1e-15})
    eigs = jacobi_eigenvalues(res.x, N)
    return eigs, res

# ── Approach B: Banded matrix with bandwidth k ──────────────────────────────
def banded_eigenvalues(params, N, k):
    """Build symmetric banded matrix with bandwidth k."""
    M = np.zeros((N, N))
    idx = 0
    for d in range(k + 1):
        n_entries = N - d
        vals = params[idx:idx + n_entries]
        if d == 0:
            M += np.diag(vals)
        else:
            M += np.diag(vals, d) + np.diag(vals, -d)
        idx += n_entries
    return np.sort(np.linalg.eigvalsh(M))

def banded_param_count(N, k):
    return sum(N - d for d in range(k + 1))

def banded_loss(params, targets, k):
    N = len(targets)
    eigs = banded_eigenvalues(params, N, k)
    return np.sum(((eigs - targets) / targets)**2)

def fit_banded(targets, k, maxiter=2000):
    N = len(targets)
    n_params = banded_param_count(N, k)

    # Initialize
    x0 = np.zeros(n_params)
    idx = 0
    # diagonal ~ targets
    x0[:N] = targets
    idx = N
    for d in range(1, k + 1):
        n_entries = N - d
        spacing = np.diff(targets)
        if d <= len(spacing):
            x0[idx:idx + n_entries] = spacing[:n_entries] / (2.0 * d)
        else:
            x0[idx:idx + n_entries] = 0.1
        idx += n_entries

    res = minimize(banded_loss, x0, args=(targets, k),
                   method='L-BFGS-B', options={'maxiter': maxiter, 'ftol': 1e-15})
    eigs = banded_eigenvalues(res.x, N, k)
    return eigs, res

# ── Approach C: Structured Sparse ───────────────────────────────────────────

# C1: Circulant + diagonal (2N params: N diagonal + N first-row of circulant)
def circulant_diag_eigenvalues(params, N):
    d = params[:N]
    c = params[N:2*N]
    # Circulant matrix from first row c
    C = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            C[i, j] = c[(j - i) % N]
    # Make Hermitian: (C + C^T)/2 + diag(d)
    M = (C + C.T) / 2.0 + np.diag(d)
    return np.sort(np.linalg.eigvalsh(M))

def circulant_diag_loss(params, targets):
    N = len(targets)
    eigs = circulant_diag_eigenvalues(params, N)
    return np.sum(((eigs - targets) / targets)**2)

def fit_circulant_diag(targets, maxiter=2000):
    N = len(targets)
    x0 = np.zeros(2 * N)
    x0[:N] = targets
    x0[N:2*N] = 0.0
    res = minimize(circulant_diag_loss, x0, args=(targets,),
                   method='L-BFGS-B', options={'maxiter': maxiter, 'ftol': 1e-15})
    eigs = circulant_diag_eigenvalues(res.x, N)
    return eigs, res

# C2: Toeplitz (symmetric, N params)
def toeplitz_eigenvalues(params, N):
    M = toeplitz(params[:N])
    return np.sort(np.linalg.eigvalsh(M))

def toeplitz_loss(params, targets):
    N = len(targets)
    eigs = toeplitz_eigenvalues(params, N)
    return np.sum(((eigs - targets) / targets)**2)

def fit_toeplitz(targets, maxiter=2000):
    N = len(targets)
    x0 = np.zeros(N)
    x0[0] = np.mean(targets)
    # First off-diagonal ~ mean spacing
    if N > 1:
        x0[1] = np.mean(np.diff(targets)) / 2
    res = minimize(toeplitz_loss, x0, args=(targets,),
                   method='L-BFGS-B', options={'maxiter': maxiter, 'ftol': 1e-15})
    eigs = toeplitz_eigenvalues(res.x, N)
    return eigs, res

# C3: Prime-based matrix
def prime_sieve(n):
    """Simple sieve for first n primes."""
    if n == 0:
        return []
    primes = []
    candidate = 2
    while len(primes) < n:
        is_prime = True
        for p in primes:
            if p * p > candidate:
                break
            if candidate % p == 0:
                is_prime = False
                break
        if is_prime:
            primes.append(candidate)
        candidate += 1
    return np.array(primes, dtype=float)

def prime_matrix_eigenvalues(params, N, primes):
    """M_{ij} = a * log(p_i * p_j) + b * log(gcd(p_i,p_j)) + c * delta_{ij} * log(p_i)
    params = [a, b, c, shift]"""
    a, b, c, shift = params[:4]
    p = primes[:N]
    lp = np.log(p)
    M = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            g = np.gcd(int(primes[i]), int(primes[j]))
            M[i, j] = a * (lp[i] + lp[j]) + b * np.log(g)
        M[i, i] += c * lp[i] + shift
    # Symmetrize
    M = (M + M.T) / 2.0
    return np.sort(np.linalg.eigvalsh(M))

def prime_matrix_loss(params, targets, primes):
    N = len(targets)
    eigs = prime_matrix_eigenvalues(params, N, primes)
    return np.sum(((eigs - targets) / targets)**2)

# ── Prediction test ─────────────────────────────────────────────────────────
def prediction_test(fit_func, targets_train, targets_test, label, **kwargs):
    """Fit on train zeros, predict test zeros by fitting N+len(test) matrix
    but only optimizing the additional parameters."""
    N_train = len(targets_train)
    N_total = N_train + len(targets_test)
    all_targets = np.concatenate([targets_train, targets_test])

    eigs, res = fit_func(all_targets, **kwargs)
    pred_test = eigs[N_train:]
    pred_err = relative_error(pred_test, targets_test)
    return pred_err, pred_test

# ── Main experiment ─────────────────────────────────────────────────────────
def main():
    results = {}
    results_text = []
    results_text.append("# Sparse Matrix Model for Zeta Zeros\n")
    results_text.append(f"Date: 2026-04-05\n")
    results_text.append(f"Zeros loaded: {len(ALL_ZEROS)}\n")

    N_values = [50, 100, 200]

    # ═══════════════════════════════════════════════════════════════════════
    # APPROACH A: Tridiagonal (Jacobi)
    # ═══════════════════════════════════════════════════════════════════════
    results_text.append("\n## Approach A: Tridiagonal (Jacobi) Matrix\n")
    results_text.append("Parameters: 2N-1 (diagonal + off-diagonal)\n")

    for N in N_values:
        targets = ALL_ZEROS[:N]
        t0 = time.time()
        eigs, res = fit_jacobi(targets, maxiter=3000)
        elapsed = time.time() - t0

        rel_err = relative_error(eigs, targets)
        rms = rms_error(eigs, targets)
        max_rel = max_relative_error(eigs, targets)

        key = f"jacobi_N{N}"
        results[key] = {
            'rel_err': rel_err, 'rms': rms, 'max_rel': max_rel,
            'n_params': 2*N-1, 'converged': res.success
        }

        line = (f"- N={N}: params={2*N-1}, mean_rel_err={rel_err:.6e}, "
                f"max_rel_err={max_rel:.6e}, rms={rms:.4f}, "
                f"converged={res.success}, time={elapsed:.1f}s")
        print(line)
        results_text.append(line + "\n")

    # ═══════════════════════════════════════════════════════════════════════
    # APPROACH B: Banded matrix with varying bandwidth
    # ═══════════════════════════════════════════════════════════════════════
    results_text.append("\n## Approach B: Banded Matrix (varying bandwidth k)\n")

    bandwidths = [1, 2, 3, 5, 10]

    for N in N_values:
        targets = ALL_ZEROS[:N]
        results_text.append(f"\n### N = {N}\n")

        for k in bandwidths:
            if k >= N:
                continue
            n_params = banded_param_count(N, k)
            # Skip if too many params for large N (would be slow)
            if n_params > 3000:
                results_text.append(f"- k={k}: SKIPPED (too many params: {n_params})\n")
                continue

            t0 = time.time()
            eigs, res = fit_banded(targets, k, maxiter=3000)
            elapsed = time.time() - t0

            rel_err = relative_error(eigs, targets)
            rms = rms_error(eigs, targets)
            max_rel = max_relative_error(eigs, targets)

            key = f"banded_N{N}_k{k}"
            results[key] = {
                'rel_err': rel_err, 'rms': rms, 'max_rel': max_rel,
                'n_params': n_params, 'converged': res.success
            }

            line = (f"- k={k}: params={n_params}, mean_rel_err={rel_err:.6e}, "
                    f"max_rel_err={max_rel:.6e}, rms={rms:.4f}, "
                    f"converged={res.success}, time={elapsed:.1f}s")
            print(line)
            results_text.append(line + "\n")

    # ═══════════════════════════════════════════════════════════════════════
    # APPROACH C: Structured Sparse
    # ═══════════════════════════════════════════════════════════════════════
    results_text.append("\n## Approach C: Structured Sparse\n")

    # C1: Circulant + diagonal
    results_text.append("\n### C1: Circulant + Diagonal (2N params)\n")
    for N in [50, 100]:  # Skip 200 -- too slow for circulant
        targets = ALL_ZEROS[:N]
        t0 = time.time()
        eigs, res = fit_circulant_diag(targets, maxiter=2000)
        elapsed = time.time() - t0

        rel_err = relative_error(eigs, targets)
        rms = rms_error(eigs, targets)
        max_rel = max_relative_error(eigs, targets)

        line = (f"- N={N}: params={2*N}, mean_rel_err={rel_err:.6e}, "
                f"max_rel_err={max_rel:.6e}, rms={rms:.4f}, "
                f"converged={res.success}, time={elapsed:.1f}s")
        print(line)
        results_text.append(line + "\n")

    # C2: Toeplitz
    results_text.append("\n### C2: Toeplitz (N params)\n")
    for N in [50, 100]:
        targets = ALL_ZEROS[:N]
        t0 = time.time()
        eigs, res = fit_toeplitz(targets, maxiter=2000)
        elapsed = time.time() - t0

        rel_err = relative_error(eigs, targets)
        rms = rms_error(eigs, targets)
        max_rel = max_relative_error(eigs, targets)

        line = (f"- N={N}: params={N}, mean_rel_err={rel_err:.6e}, "
                f"max_rel_err={max_rel:.6e}, rms={rms:.4f}, "
                f"converged={res.success}, time={elapsed:.1f}s")
        print(line)
        results_text.append(line + "\n")

    # C3: Prime-based (only 4 params -- very constrained)
    results_text.append("\n### C3: Prime-based M_{ij} = a*log(pi*pj) + b*log(gcd) + c*delta*log(pi) + shift (4 params)\n")
    primes = prime_sieve(200)
    for N in [50]:
        targets = ALL_ZEROS[:N]
        x0 = np.array([1.0, 0.1, 1.0, np.mean(targets)])
        t0 = time.time()
        res = minimize(prime_matrix_loss, x0, args=(targets, primes),
                       method='Nelder-Mead', options={'maxiter': 5000})
        eigs = prime_matrix_eigenvalues(res.x, N, primes)
        elapsed = time.time() - t0

        rel_err = relative_error(eigs, targets)
        rms = rms_error(eigs, targets)

        line = (f"- N={N}: params=4, mean_rel_err={rel_err:.6e}, "
                f"rms={rms:.4f}, converged={res.success}, time={elapsed:.1f}s")
        print(line)
        results_text.append(line + "\n")
        results_text.append(f"  Best params: a={res.x[0]:.4f}, b={res.x[1]:.4f}, c={res.x[2]:.4f}, shift={res.x[3]:.4f}\n")

    # ═══════════════════════════════════════════════════════════════════════
    # PREDICTION TEST
    # ═══════════════════════════════════════════════════════════════════════
    results_text.append("\n## Prediction Test\n")
    results_text.append("Fit on first N zeros, check fit quality on zeros N+1..N+50\n")

    for N in [50, 100]:
        targets_train = ALL_ZEROS[:N]
        targets_test = ALL_ZEROS[N:N+50]
        # Fit Jacobi on N+50 zeros total
        targets_all = ALL_ZEROS[:N+50]
        t0 = time.time()
        eigs_all, res_all = fit_jacobi(targets_all, maxiter=3000)
        elapsed = time.time() - t0

        # How well do we fit the training set?
        train_err = relative_error(eigs_all[:N], targets_train)
        # How well do we fit the test set?
        test_err = relative_error(eigs_all[N:N+50], targets_test)

        line = (f"- Jacobi N={N}, predict next 50: "
                f"train_rel_err={train_err:.6e}, test_rel_err={test_err:.6e}, "
                f"time={elapsed:.1f}s")
        print(line)
        results_text.append(line + "\n")

        # Also: fit ONLY on N zeros, then extend matrix to N+50 by padding
        # This is the real test: can we extrapolate?
        results_text.append(f"  (Note: this fits all N+50 zeros jointly, not a true extrapolation test.)\n")

    # True extrapolation: fit N, freeze params, extend matrix
    results_text.append("\n### True Extrapolation (Jacobi)\n")
    results_text.append("Fit N zeros -> get 2N-1 params. Extend to (N+50) matrix by\n")
    results_text.append("extrapolating diagonal/off-diagonal trends, then check eigenvalues.\n")

    for N in [50, 100]:
        targets_train = ALL_ZEROS[:N]
        targets_test = ALL_ZEROS[N:N+50]

        _, res = fit_jacobi(targets_train, maxiter=3000)
        a_fit = res.x[:N]
        b_fit = res.x[N:2*N-1]

        # Extrapolate: linear trend in a and b
        # Fit linear model to last 20 values
        tail = min(20, N-1)
        idx_a = np.arange(N - tail, N)
        slope_a = np.polyfit(idx_a, a_fit[N-tail:], 1)
        idx_b = np.arange(N - 1 - tail, N - 1)
        slope_b = np.polyfit(idx_b, b_fit[N-1-tail:], 1)

        # Extend
        new_idx_a = np.arange(N, N + 50)
        new_a = np.polyval(slope_a, new_idx_a)
        a_ext = np.concatenate([a_fit, new_a])

        new_idx_b = np.arange(N - 1, N - 1 + 50)
        new_b = np.polyval(slope_b, new_idx_b)
        b_ext = np.concatenate([b_fit, new_b])

        N2 = N + 50
        M = np.diag(a_ext) + np.diag(b_ext, 1) + np.diag(b_ext, -1)
        eigs_ext = np.sort(np.linalg.eigvalsh(M))

        pred_test = eigs_ext[N:N+50]
        extrap_err = relative_error(pred_test, targets_test)
        max_extrap = max_relative_error(pred_test, targets_test)

        line = (f"- N={N}: extrapolated mean_rel_err={extrap_err:.6e}, "
                f"max_rel_err={max_extrap:.6e}")
        print(line)
        results_text.append(line + "\n")

        # Show first 5 predicted vs actual
        results_text.append(f"  First 5 predicted vs actual:\n")
        for i in range(5):
            results_text.append(f"    gamma_{N+i+1}: predicted={pred_test[i]:.6f}, "
                               f"actual={targets_test[i]:.6f}, "
                               f"rel_err={abs(pred_test[i]-targets_test[i])/targets_test[i]:.6e}\n")

    # ═══════════════════════════════════════════════════════════════════════
    # KEY ANALYSIS: Bandwidth vs N scaling
    # ═══════════════════════════════════════════════════════════════════════
    results_text.append("\n## Key Analysis: Minimum Bandwidth for <1% Error\n")

    threshold = 0.01  # 1% mean relative error

    for N in N_values:
        found = False
        for k in [1, 2, 3, 5, 10]:
            key = f"banded_N{N}_k{k}"
            if key in results and results[key]['rel_err'] < threshold:
                results_text.append(
                    f"- N={N}: minimum bandwidth k={k} achieves "
                    f"mean_rel_err={results[key]['rel_err']:.6e} (<1%)\n"
                    f"  Parameters: {results[key]['n_params']} "
                    f"(vs N^2={N*N} for dense)\n"
                    f"  Compression ratio: {N*N / results[key]['n_params']:.1f}x\n"
                )
                found = True
                break
        if not found:
            # Check if any achieved it
            best_key = None
            best_err = float('inf')
            for k in [1, 2, 3, 5, 10]:
                key = f"banded_N{N}_k{k}"
                if key in results and results[key]['rel_err'] < best_err:
                    best_err = results[key]['rel_err']
                    best_key = key
            if best_key:
                results_text.append(
                    f"- N={N}: NO bandwidth k<=10 achieved <1% error. "
                    f"Best: {best_key} with mean_rel_err={best_err:.6e}\n"
                )
            else:
                results_text.append(f"- N={N}: no results available\n")

    # ═══════════════════════════════════════════════════════════════════════
    # ANALYSIS: Can matrix entries be expressed as simple functions?
    # ═══════════════════════════════════════════════════════════════════════
    results_text.append("\n## Analysis: Structure of Fitted Parameters\n")

    for N in [50, 100]:
        targets = ALL_ZEROS[:N]
        _, res = fit_jacobi(targets, maxiter=3000)
        a_fit = res.x[:N]
        b_fit = res.x[N:2*N-1]

        # Check if diagonal ~ linear in index
        idx = np.arange(N)
        poly_a = np.polyfit(idx, a_fit, 1)
        a_pred_linear = np.polyval(poly_a, idx)
        a_linear_err = np.mean(np.abs(a_fit - a_pred_linear) / np.abs(a_fit))

        # Check if diagonal ~ targets (trivial)
        a_trivial_err = np.mean(np.abs(a_fit - targets) / np.abs(targets))

        # Check off-diagonal: constant? linear?
        idx_b = np.arange(len(b_fit))
        poly_b = np.polyfit(idx_b, b_fit, 1)
        b_pred_linear = np.polyval(poly_b, idx_b)
        b_linear_err = np.mean(np.abs(b_fit - b_pred_linear) / (np.abs(b_fit) + 1e-10))

        b_mean = np.mean(b_fit)
        b_const_err = np.mean(np.abs(b_fit - b_mean) / (np.abs(b_fit) + 1e-10))

        results_text.append(f"\n### N={N} Jacobi parameter structure\n")
        results_text.append(f"- Diagonal a_i vs gamma_i (trivial match): mean_rel_err={a_trivial_err:.6e}\n")
        results_text.append(f"- Diagonal a_i ~ linear(i): slope={poly_a[0]:.4f}, intercept={poly_a[1]:.4f}, "
                           f"mean_rel_err={a_linear_err:.6e}\n")
        results_text.append(f"- Off-diagonal b_i ~ constant({b_mean:.4f}): mean_rel_err={b_const_err:.6e}\n")
        results_text.append(f"- Off-diagonal b_i ~ linear(i): slope={poly_b[0]:.6f}, intercept={poly_b[1]:.4f}, "
                           f"mean_rel_err={b_linear_err:.6e}\n")
        results_text.append(f"- Off-diagonal stats: mean={b_mean:.4f}, std={np.std(b_fit):.4f}, "
                           f"min={np.min(b_fit):.4f}, max={np.max(b_fit):.4f}\n")

    # ═══════════════════════════════════════════════════════════════════════
    # CONCLUSIONS
    # ═══════════════════════════════════════════════════════════════════════
    results_text.append("\n## Conclusions\n")
    results_text.append("""
The key question is whether N zeta zeros can be encoded by O(N) parameters
(sparse matrix) rather than requiring O(N) independent numbers.

**Trivial observation**: A diagonal matrix with entries gamma_1, ..., gamma_N
already has eigenvalues equal to the zeros -- but uses N parameters for N zeros
(no compression). The question is whether a STRUCTURED sparse matrix with
FEWER than N free parameters (or with parameters following a simple pattern)
can reproduce the zeros.

**Findings summarized in the bandwidth/error table above.**

If the Jacobi off-diagonal entries b_i follow a simple pattern (constant or
linear), then the matrix is described by ~3 parameters regardless of N --
but this would imply the zeros are nearly deterministic, contradicting GUE
statistics. The actual b_i values measure how much "randomness" is needed.
""")

    # Write results
    with open(RESULTS_PATH, 'w') as f:
        f.write('\n'.join(results_text))
    print(f"\nResults written to {RESULTS_PATH}")

if __name__ == '__main__':
    main()
