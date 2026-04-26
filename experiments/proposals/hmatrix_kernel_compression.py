"""
Proposal B1: Hierarchical-matrix (H-matrix) compression of the
explicit-formula kernel.

The Riemann explicit formula evaluated on a vector of test points x
amounts to a matrix-vector product

    [psi_0(x_i)]_i  =  K v

with K[i, k] = x_i^{1/2 + i gamma_k} / (1/2 + i gamma_k), and v_k = -1.

If K admits hierarchical low-rank structure (off-diagonal blocks have
small numerical rank), then a one-time HSS / HODLR / H-matrix factorisation
costs O(N r^2) and each new x evaluation costs O(N r log N), with r the
typical off-diagonal rank.

This script builds K for a range of test windows and zero-truncations,
measures off-diagonal SVD ranks at successive bisection levels, and
reports the rank-vs-block-size growth law.

Win condition: ranks plateau at O(log) or grow only sub-linearly.
Failure: ranks ~ block-size (matrix is essentially full-rank).

Usage:
    python3 hmatrix_kernel_compression.py
    python3 hmatrix_kernel_compression.py --x_lo 1e3 --x_hi 1e4 --N 64 --M 64
"""
from __future__ import annotations
import argparse
import json
import time
import numpy as np
from mpmath import zetazero, mp


def kernel_matrix(x_grid: np.ndarray, gammas: np.ndarray) -> np.ndarray:
    """K[i, k] = x_i^{1/2 + i gamma_k} / (1/2 + i gamma_k)."""
    s = 0.5 + 1j * gammas
    log_x = np.log(x_grid)
    return np.exp(np.outer(log_x, s)) / s[None, :]


def off_diagonal_ranks(K: np.ndarray, eps: float = 1e-10) -> dict:
    """
    For each dyadic level, look at off-diagonal blocks of size
    (rows) x (block-size 2^level), report epsilon-rank.

    Specifically: for level l, partition columns into blocks of width
    w = 2^l; for each block, take the row-strip *not* indexed by the
    near-diagonal rows (we don't have a column->row alignment, so
    instead we treat the whole row-set --- this is a column-side rank).
    """
    M, N = K.shape
    out: dict[int, dict[int, list[int]]] = {}
    max_level = int(np.floor(np.log2(N)))
    for level in range(max_level + 1):
        w = 2 ** level
        ranks: list[int] = []
        for j0 in range(0, N - w + 1, w):
            block = K[:, j0:j0 + w]
            if block.size == 0:
                continue
            s = np.linalg.svd(block, compute_uv=False)
            top = s[0] if s[0] > 0 else 1.0
            r = int(np.sum(s > eps * top))
            ranks.append(r)
        if ranks:
            out[level] = {
                "block_width": w,
                "n_blocks": len(ranks),
                "ranks": ranks,
                "rank_max": max(ranks),
                "rank_mean": float(np.mean(ranks)),
            }
    return out


def hodlr_offdiag_ranks(K: np.ndarray, eps: float = 1e-10) -> list[dict]:
    """
    Proper HODLR-style off-diagonal block ranks: at each bisection,
    split rows AND columns and take both off-diagonal corner blocks.
    """
    out: list[dict] = []

    def recurse(rs: int, re: int, cs: int, ce: int, depth: int) -> None:
        m = re - rs
        n = ce - cs
        if m < 2 or n < 2:
            return
        rmid = rs + m // 2
        cmid = cs + n // 2
        # off-diagonal corner blocks
        TR = K[rs:rmid, cmid:ce]
        BL = K[rmid:re, cs:cmid]
        for label, B in (("TR", TR), ("BL", BL)):
            if B.size == 0:
                continue
            s = np.linalg.svd(B, compute_uv=False)
            top = s[0] if s[0] > 0 else 1.0
            r = int(np.sum(s > eps * top))
            out.append({
                "depth": depth,
                "block_shape": list(B.shape),
                "block": label,
                "rank": r,
                "norm": float(np.linalg.norm(B)),
            })
        recurse(rs, rmid, cs, cmid, depth + 1)
        recurse(rmid, re, cmid, ce, depth + 1)

    recurse(0, K.shape[0], 0, K.shape[1], 0)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--x_lo", type=float, default=1e3)
    parser.add_argument("--x_hi", type=float, default=1e4)
    parser.add_argument("--N", type=int, default=64,
                        help="number of zeros (gammas)")
    parser.add_argument("--M", type=int, default=64,
                        help="number of x test points")
    parser.add_argument("--eps", type=float, default=1e-10)
    parser.add_argument("--mpmath_dps", type=int, default=30)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    mp.dps = args.mpmath_dps
    np.random.seed(args.seed)

    print(f"[setup] computing first {args.N} zeros at dps={args.mpmath_dps}")
    t0 = time.time()
    gammas = np.array([float(zetazero(k).imag) for k in range(1, args.N + 1)])
    print(f"[setup] zeros computed in {time.time()-t0:.2f}s")

    xs = np.linspace(args.x_lo, args.x_hi, args.M)

    print(f"[build] kernel shape = ({args.M}, {args.N})")
    K = kernel_matrix(xs, gammas)
    print(f"[build] |K|_F = {np.linalg.norm(K):.3e}, "
          f"max|K| = {np.max(np.abs(K)):.3e}")

    full_svd = np.linalg.svd(K, compute_uv=False)
    full_rank = int(np.sum(full_svd > args.eps * full_svd[0]))
    print(f"[svd] full numerical rank (eps={args.eps}) = {full_rank} "
          f"out of min(M,N)={min(args.M,args.N)}")

    print("\n[column-bisection ranks]")
    col_ranks = off_diagonal_ranks(K, eps=args.eps)
    for level, info in col_ranks.items():
        print(f"  level {level}: block_width={info['block_width']:>4d}  "
              f"n_blocks={info['n_blocks']:>3d}  "
              f"rank_max={info['rank_max']:>3d}  "
              f"rank_mean={info['rank_mean']:.2f}")

    print("\n[HODLR off-diagonal corner block ranks]")
    hodlr = hodlr_offdiag_ranks(K, eps=args.eps)
    by_depth: dict[int, list[int]] = {}
    for entry in hodlr:
        by_depth.setdefault(entry["depth"], []).append(entry["rank"])
    for depth in sorted(by_depth.keys()):
        rs = by_depth[depth]
        print(f"  depth {depth}: n_blocks={len(rs):>3d}  "
              f"rank_max={max(rs):>3d}  rank_mean={np.mean(rs):.2f}")

    # rank growth law analysis
    print("\n[rank growth law]")
    levels = sorted(col_ranks.keys())
    if levels:
        log_w = np.log2([col_ranks[l]['block_width'] for l in levels])
        rmax = np.array([col_ranks[l]['rank_max'] for l in levels])
        # fit r ~ a * log_w + b  vs  r ~ c * w
        # linear-in-log fit
        A = np.column_stack([log_w, np.ones_like(log_w)])
        coef_log, *_ = np.linalg.lstsq(A, rmax, rcond=None)
        a_log, b_log = coef_log
        # linear-in-w fit
        w = 2.0 ** log_w
        A2 = np.column_stack([w, np.ones_like(w)])
        coef_lin, *_ = np.linalg.lstsq(A2, rmax, rcond=None)
        a_lin, b_lin = coef_lin
        # residuals
        res_log = float(np.linalg.norm(rmax - (a_log * log_w + b_log)))
        res_lin = float(np.linalg.norm(rmax - (a_lin * w + b_lin)))
        print(f"  fit  rank_max ~ {a_log:.3f} * log2(w) + {b_log:.3f}  "
              f"residual={res_log:.3f}")
        print(f"  fit  rank_max ~ {a_lin:.5f} * w + {b_lin:.3f}  "
              f"residual={res_lin:.3f}")

    out = {
        "args": vars(args),
        "n_gammas": int(args.N),
        "gammas_first5": gammas[:5].tolist(),
        "full_rank_eps": full_rank,
        "K_frobenius": float(np.linalg.norm(K)),
        "col_bisection_ranks": col_ranks,
        "hodlr_summary": {str(d): {
            "n_blocks": len(rs),
            "rank_max": max(rs),
            "rank_mean": float(np.mean(rs)),
        } for d, rs in by_depth.items()},
    }

    out_json = (
        "/apps/aplikacijos/prime-research/experiments/proposals/"
        "hmatrix_kernel_compression_data.json"
    )
    with open(out_json, "w") as f:
        json.dump(out, f, indent=2, default=float)
    print(f"\n[save] wrote {out_json}")


if __name__ == "__main__":
    main()
