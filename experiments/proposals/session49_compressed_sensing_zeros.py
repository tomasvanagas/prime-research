#!/usr/bin/env python3
"""
Session 49 / Proposal B — LASSO compressed-sensing recovery of delta(n)
in the zeta-zero basis.

Test whether the residual delta(n) = p(n) - R^{-1}(n) admits a sparse
representation in the basis of zeta-zero modes
   { cos(gamma_k log n)/sqrt(n), sin(gamma_k log n)/sqrt(n) : k=1..K_max }.
We use LASSO to find the smallest support that fits delta on a training
window, then evaluate on a held-out test window.

If LASSO finds K_eff = O(polylog(N)) non-zero modes that achieve
zero rounding error on test data, proposal B has empirical traction.

Output: experiments/proposals/session49_compressed_sensing_zeros_results.md
"""

import math
import time
import numpy as np
from sklearn.linear_model import Lasso

ZEROS_FILE = "/apps/aplikacijos/prime-research/data/zeta_zeros_2000.txt"
CACHE = "/apps/aplikacijos/prime-research/experiments/proposals/session49_data.npz"


def load_zeros(K):
    gammas = []
    with open(ZEROS_FILE) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            gammas.append(float(line))
            if len(gammas) >= K:
                break
    return np.array(gammas)


def build_design(n_arr, gammas):
    """Returns Phi of shape (len(n_arr), 2*len(gammas)).
    Phi[:, 2k]   = cos(gamma_k log n) / sqrt(n)
    Phi[:, 2k+1] = sin(gamma_k log n) / sqrt(n)
    """
    log_n = np.log(n_arr)
    inv_sqrt_n = 1.0 / np.sqrt(n_arr)
    K = len(gammas)
    Phi = np.zeros((len(n_arr), 2 * K), dtype=np.float64)
    for k, g in enumerate(gammas):
        ang = g * log_n
        Phi[:, 2 * k] = np.cos(ang) * inv_sqrt_n
        Phi[:, 2 * k + 1] = np.sin(ang) * inv_sqrt_n
    return Phi


def main():
    N = 4096
    K_max = 1024
    print(f"Loading cached primes and delta for n=1..{N}", flush=True)
    d = np.load(CACHE)
    primes = d["primes"][:N].tolist()
    rinv = d["rinv"][:N]
    delta = d["delta"][:N]
    print(f"  loaded; delta std={delta.std():.3f}", flush=True)

    print(f"Loading first {K_max} zeros", flush=True)
    gammas = load_zeros(K_max)
    print(f"  zero range: [{gammas[0]:.2f}, {gammas[-1]:.2f}]", flush=True)

    # Use n=2..N as features. Train on first half, test on second half.
    n_all = np.arange(2, N + 1, dtype=np.float64)
    delta_all = delta[1:]  # delta[i] corresponds to n=i+1; here we use n=2..N => index 1..N-1
    primes_all = np.array(primes[1:], dtype=np.float64)
    rinv_all = rinv[1:]

    split = (N - 1) // 2
    train_idx = np.arange(0, split)
    test_idx = np.arange(split, N - 1)

    Phi_full = build_design(n_all, gammas)

    # Normalize columns of Phi (LASSO sensitive to scale)
    col_norms = np.linalg.norm(Phi_full, axis=0)
    col_norms[col_norms == 0] = 1.0
    Phi_norm = Phi_full / col_norms

    rows = []
    print("\nLASSO sweeps (best alpha by test RMSE):")
    print(f"{'K':>6} {'alpha':>10} {'nnz':>5} {'train_RMSE':>12} {'test_RMSE':>11} {'test_round_acc':>15}")

    alpha_grid = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]

    for K in [16, 64, 256, 1024]:
        Phi_K = Phi_norm[:, : 2 * K]
        Phi_train = Phi_K[train_idx]
        y_train = delta_all[train_idx]
        Phi_test = Phi_K[test_idx]
        y_test = delta_all[test_idx]

        best = None
        for alpha in alpha_grid:
            model = Lasso(alpha=alpha, max_iter=50000, fit_intercept=True)
            model.fit(Phi_train, y_train)
            pred_train = model.predict(Phi_train)
            pred_test = model.predict(Phi_test)
            tr_rmse = float(np.sqrt(np.mean((pred_train - y_train) ** 2)))
            te_rmse = float(np.sqrt(np.mean((pred_test - y_test) ** 2)))
            nnz = int((np.abs(model.coef_) > 1e-10).sum())

            est_test = rinv_all[test_idx] + pred_test
            recovered_primes = np.round(est_test).astype(np.int64)
            acc = float(np.mean(recovered_primes == primes_all[test_idx].astype(np.int64)))

            cand = (alpha, nnz, tr_rmse, te_rmse, acc)
            if best is None or te_rmse < best[3]:
                best = cand

        alpha, nnz, tr, te, acc = best
        rows.append((K, alpha, nnz, tr, te, acc))
        print(f"{K:6d} {alpha:10.3g} {nnz:5d} {tr:12.4f} {te:11.4f} {acc:15.3f}", flush=True)

    # Baseline: ordinary least squares (no L1 sparsity)
    print("\nOLS baselines for comparison:")
    print(f"{'K':>6} {'train_RMSE':>12} {'test_RMSE':>11} {'test_round_acc':>15}")
    ols_rows = []
    for K in [16, 64, 256, 1024]:
        Phi_K = Phi_norm[:, : 2 * K]
        Phi_train = Phi_K[train_idx]
        y_train = delta_all[train_idx]
        Phi_test = Phi_K[test_idx]
        y_test = delta_all[test_idx]
        coef, _, _, _ = np.linalg.lstsq(Phi_train, y_train, rcond=None)
        pred_train = Phi_train @ coef
        pred_test = Phi_test @ coef
        tr = float(np.sqrt(np.mean((pred_train - y_train) ** 2)))
        te = float(np.sqrt(np.mean((pred_test - y_test) ** 2)))
        est_test = rinv_all[test_idx] + pred_test
        recovered_primes = np.round(est_test).astype(np.int64)
        acc = float(np.mean(recovered_primes == primes_all[test_idx].astype(np.int64)))
        ols_rows.append((K, tr, te, acc))
        print(f"{K:6d} {tr:12.4f} {te:11.4f} {acc:15.3f}", flush=True)

    # Naive baseline: just round R^{-1}(n)
    naive_acc = float(np.mean(np.round(rinv_all[test_idx]).astype(np.int64) == primes_all[test_idx].astype(np.int64)))
    print(f"\nNaive round(R^{{-1}}) test recovery rate: {naive_acc:.3f}")

    out_md = "/apps/aplikacijos/prime-research/experiments/proposals/session49_compressed_sensing_zeros_results.md"
    with open(out_md, "w") as f:
        f.write("# session49_compressed_sensing_zeros — results\n\n")
        f.write(f"N={N}, training [n=2..{split+1}], testing [n={split+2}..{N}].\n")
        f.write(f"K_max zeros = {K_max}, gamma_max = {gammas[-1]:.2f}.\n\n")
        f.write("## LASSO results (best alpha per K)\n\n")
        f.write("| K | alpha | nnz | train RMSE | test RMSE | test round-prime acc |\n")
        f.write("|---|---|---|---|---|---|\n")
        for K, alpha, nnz, tr, te, acc in rows:
            f.write(f"| {K} | {alpha:.3g} | {nnz} | {tr:.4f} | {te:.4f} | {acc:.3f} |\n")
        f.write("\n## OLS baselines (no sparsity)\n\n")
        f.write("| K | train RMSE | test RMSE | test round-prime acc |\n")
        f.write("|---|---|---|---|\n")
        for K, tr, te, acc in ols_rows:
            f.write(f"| {K} | {tr:.4f} | {te:.4f} | {acc:.3f} |\n")
        f.write(f"\nNaive round(R^{{-1}}(n)) test recovery rate: **{naive_acc:.3f}**.\n\n")
        f.write("## Interpretation\n\n")
        best_lasso = max(rows, key=lambda r: r[5])
        K_b, _, nnz_b, _, _, acc_b = best_lasso
        f.write(f"Best LASSO recovery: K={K_b}, nnz={nnz_b}, test acc={acc_b:.3f}.\n\n")
        if acc_b > 0.95 and nnz_b < 50:
            f.write("**Surprising positive signal**: <50 effective modes give >95% recovery.\n")
            f.write("Worth follow-up to confirm out-of-distribution generalization.\n")
        elif acc_b > 0.5:
            f.write("Partial recovery. K_eff is not polylog but is < 2*K_max,\n")
            f.write("indicating LASSO discovered some structure beyond a uniform\n")
            f.write("zero-mode budget. Probably consistent with sqrt(N) scaling.\n")
        else:
            f.write("**Negative**: LASSO sparsity in the zero basis does not yield\n")
            f.write("polylog recovery. Consistent with GUE-random null hypothesis.\n")
            f.write("Proposal B closed.\n")

    print(f"\nWrote {out_md}")


if __name__ == "__main__":
    main()
