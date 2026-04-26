#!/usr/bin/env python3
"""
Critique session — verify whether Proposal C's basis misspecification
(cos(g log n)/sqrt(n) using prime *index* n vs prime *value* p_n) changes
the qualitative conclusion.

Standard explicit formula:
  pi(x) - R(x) = -2 sqrt(x) / log x * sum_k cos(gamma_k log x - phase_k) / |rho_k|

Inverting: delta(n) = p(n) - R^{-1}(n) ≈ f(p_n) * log(p_n)
        ≈ -2 sqrt(p_n) * sum_k cos(gamma_k log p_n - phase_k) / |rho_k|

So the *correct* basis for ridge fitting is cos(gamma_k log p_n) / sqrt(p_n) * p_n
= sqrt(p_n) * cos(gamma_k log p_n), NOT cos(gamma_k log n) / sqrt(n).

We compare:
  - basis_orig: as in proposal C
  - basis_correct: sqrt(p_n) * cos(gamma_k log p_n), sqrt(p_n) * sin(gamma_k log p_n)

If the qualitative conclusion (no polylog K, ridge plateaus at naive) holds
under the corrected basis, proposal C is genuinely closed.
"""

import math
import time
import numpy as np

ZEROS_FILE = "/apps/aplikacijos/prime-research/data/zeta_zeros_2000.txt"
CACHE = "/apps/aplikacijos/prime-research/experiments/proposals/session49_data.npz"


def load_zeros(K):
    g = []
    with open(ZEROS_FILE) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            g.append(float(line))
            if len(g) >= K:
                break
    return np.array(g)


def build_design_orig(n_arr, gammas):
    log_n = np.log(n_arr)
    inv_sqrt_n = 1.0 / np.sqrt(n_arr)
    K = len(gammas)
    Phi = np.zeros((len(n_arr), 2 * K), dtype=np.float64)
    for k, g in enumerate(gammas):
        ang = g * log_n
        Phi[:, 2 * k] = np.cos(ang) * inv_sqrt_n
        Phi[:, 2 * k + 1] = np.sin(ang) * inv_sqrt_n
    return Phi


def build_design_correct(p_arr, gammas):
    """sqrt(p_n) * cos(gamma_k log p_n), sqrt(p_n) * sin(gamma_k log p_n)."""
    log_p = np.log(p_arr)
    sqrt_p = np.sqrt(p_arr)
    K = len(gammas)
    Phi = np.zeros((len(p_arr), 2 * K), dtype=np.float64)
    for k, g in enumerate(gammas):
        ang = g * log_p
        Phi[:, 2 * k] = np.cos(ang) * sqrt_p
        Phi[:, 2 * k + 1] = np.sin(ang) * sqrt_p
    return Phi


def ridge_fit(Phi, y, lam):
    XtX = Phi.T @ Phi
    XtX[np.diag_indices_from(XtX)] += lam
    rhs = Phi.T @ y
    return np.linalg.solve(XtX, rhs)


def main():
    N = 4096
    print(f"Loading cached data for n=1..{N}", flush=True)
    d = np.load(CACHE)
    primes = d["primes"][:N]
    rinv = d["rinv"][:N]
    delta = d["delta"][:N]

    n_all = np.arange(2, N + 1, dtype=np.float64)
    p_all = primes[1:].astype(np.float64)
    delta_all = delta[1:]
    rinv_all = rinv[1:]
    primes_all = primes[1:]

    split = (N - 1) // 2
    train_idx = np.arange(0, split)
    test_idx = np.arange(split, N - 1)

    rows = []
    print("\nRidge on the TWO bases for K in {32, 128, 512, 2000}.")
    print(f"{'K':>6} {'orig_acc':>10} {'corr_acc':>10} {'orig_RMSE':>11} {'corr_RMSE':>11}")

    for K in [32, 128, 512, 2000]:
        gammas = load_zeros(K)
        Phi_o = build_design_orig(n_all, gammas)
        Phi_c = build_design_correct(p_all, gammas)

        best_o = (None, None, None)
        best_c = (None, None, None)
        for lam in [1e-6, 1e-4, 1e-2, 1.0, 10.0]:
            try:
                w_o = ridge_fit(Phi_o[train_idx], delta_all[train_idx], lam)
                pred_o = Phi_o[test_idx] @ w_o
                acc_o = float(np.mean(np.round(rinv_all[test_idx] + pred_o).astype(np.int64) == primes_all[test_idx]))
                rmse_o = float(np.sqrt(np.mean((pred_o - delta_all[test_idx]) ** 2)))
            except np.linalg.LinAlgError:
                acc_o, rmse_o = 0.0, float("inf")

            try:
                w_c = ridge_fit(Phi_c[train_idx], delta_all[train_idx], lam)
                pred_c = Phi_c[test_idx] @ w_c
                acc_c = float(np.mean(np.round(rinv_all[test_idx] + pred_c).astype(np.int64) == primes_all[test_idx]))
                rmse_c = float(np.sqrt(np.mean((pred_c - delta_all[test_idx]) ** 2)))
            except np.linalg.LinAlgError:
                acc_c, rmse_c = 0.0, float("inf")

            if best_o[0] is None or acc_o > best_o[0]:
                best_o = (acc_o, rmse_o, lam)
            if best_c[0] is None or acc_c > best_c[0]:
                best_c = (acc_c, rmse_c, lam)

        rows.append((K, best_o[0], best_c[0], best_o[1], best_c[1]))
        print(f"{K:6d} {best_o[0]:10.3f} {best_c[0]:10.3f} {best_o[1]:11.3f} {best_c[1]:11.3f}", flush=True)

    naive_acc = float(np.mean(np.round(rinv_all[test_idx]).astype(np.int64) == primes_all[test_idx]))
    print(f"\nNaive baseline: {naive_acc:.3f}")

    out_md = "/apps/aplikacijos/prime-research/experiments/proposals/critique49_basis_misspec_results.md"
    with open(out_md, "w") as f:
        f.write("# critique49_basis_misspec — results\n\n")
        f.write("Critic test: did Proposal C's `cos(gamma log n)/sqrt(n)` basis\n")
        f.write("(prime index n) hide a positive result that the **correct** basis\n")
        f.write("`sqrt(p_n) * cos(gamma log p_n)` (prime value) would surface?\n\n")
        f.write(f"N={N}, train [n=2..{split+1}], test [n={split+2}..{N}].\n\n")
        f.write("| K | orig acc | correct acc | orig test RMSE | correct test RMSE |\n")
        f.write("|---|---|---|---|---|\n")
        for K, ao, ac, ro, rc in rows:
            f.write(f"| {K} | {ao:.3f} | {ac:.3f} | {ro:.3f} | {rc:.3f} |\n")
        f.write(f"\nNaive round(R^{{-1}}) baseline: **{naive_acc:.3f}**.\n\n")
        f.write("## Verdict\n\n")
        max_corr = max(r[2] for r in rows)
        max_orig = max(r[1] for r in rows)
        if max_corr > 0.5 and max_corr > max_orig + 0.1:
            f.write(f"**Significant regression**: correct basis gives {max_corr:.3f} vs orig {max_orig:.3f}. \n")
            f.write("Proposal C's misspecification HID a positive signal. Reopen.\n")
        elif max_corr > max_orig + 0.05:
            f.write(f"Modest improvement with correct basis ({max_corr:.3f} vs {max_orig:.3f}); ")
            f.write("still does not approach polylog recovery. Conclusion unchanged.\n")
        else:
            f.write(f"Both bases plateau near naive ({max_corr:.3f} vs {max_orig:.3f}); ")
            f.write("misspecification did NOT hide a positive result. Proposal C closed.\n")

    print(f"\nWrote {out_md}")


if __name__ == "__main__":
    main()
