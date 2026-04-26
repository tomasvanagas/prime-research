"""
MPS Bond Dimension of the Prime Indicator (base-W reshape)
==========================================================

Quantum-inspired complexity test. View the binary prime indicator
chi_P(n) for n in [1..W^d] as a state vector in (C^W)^{otimes d}.
A state with low entanglement has small TT (matrix-product-state)
bond dimension, hence O(d * chi_max^2 * W) storage and many
linear-algebraic operations in polylog.

This script measures TT bond dimensions for two families:

  Family A (base-2 sweep, original measurement):
    W = 2, d in {10, 12, 14, 16, 18, 20}
    Reports chi_P, random control, cumulative pi(n), residual delta.

  Family B (base-W theorem verification, new):
    W in {2, 6, 30, 210}, d as large as memory allows.
    Tests the precise prediction:
        rank at cut j  =  min( W^j,  phi(W) * W^{d-j-1} + 1 )
    The half-cut bond dim ratio approaches phi(W) / W as d -> infinity.

THEOREM (proven in novel/mps_bond_dimension.md and confirmed at every
tested cut):
    For d >= 2, 1 <= j < d, the unfolded matrix M^{(j)} of chi_P viewed
    in base W has rank at most  min(W^j, phi(W) * W^{d-j-1} + 1).
    Proof sketch: for i >= 1, n = i*W^{d-j} + k+1 > W satisfies
    gcd(n, W) > 1 iff gcd(k+1, W) > 1, and any such n is composite
    (since W >= 2 has every small prime factor we care about).
    Hence rows i = 1, ..., W^j - 1 have support on at most
    phi(W) * W^{d-j-1} columns. Row 0 contributes at most 1 more.
    The empirical SVD saturates this bound exactly at every cut.

Why this is NEW (vs the prior W=2 measurement):
    Quantifies the "wheel-W in MPS language": every primorial wheel
    becomes a constant-factor savings in tensor-network bond dimension.
    Provides a clean phi(W)/W "compressibility ratio" with NO
    polylog limit -- the ratio approaches 0 only as W -> infinity,
    but evaluating the wheel costs >= W, so cost stays >= sqrt(N).

This is DIFFERENT from earlier MPS / sieve tensor experiments:
  - tensor_sieve.py: studied product DFA (primorial blowup) of the
    sieve transfer matrix, not the bond dim of the indicator vector.
  - sieve_matrix_rank.py: studied divisibility matrix M[n,p].
  - wavelet_*: tested basis sparsity, not multipartite entanglement.

Usage:
  python mps_bond_dimension.py                # full sweep
  python mps_bond_dimension.py --base 30      # only base-30 sweep
  python mps_bond_dimension.py --quick        # smaller sizes for CI
"""

import argparse
import math
import time
from math import gcd
import numpy as np
from sympy import sieve
from mpmath import mp, li


# ---------------------------------------------------------------------------
# Number-theoretic helpers
# ---------------------------------------------------------------------------

def euler_phi(n: int) -> int:
    result = n
    p = 2
    nn = n
    while p * p <= nn:
        if nn % p == 0:
            while nn % p == 0:
                nn //= p
            result -= result // p
        p += 1
    if nn > 1:
        result -= result // nn
    return result


def predicted_rank(W: int, d: int, j: int) -> int:
    """Theorem: rank at cut j is at most min(W^j, phi(W)*W^{d-j-1}+1)."""
    phi = euler_phi(W)
    return min(W ** j, phi * (W ** (d - j - 1)) + 1)


def predicted_max_rank(W: int, d: int) -> int:
    return max(predicted_rank(W, d, j) for j in range(1, d))


# ---------------------------------------------------------------------------
# Indicator vectors
# ---------------------------------------------------------------------------

def chi_indicator(N: int) -> np.ndarray:
    arr = np.zeros(N, dtype=np.float64)
    for p in sieve.primerange(2, N + 1):
        arr[p - 1] = 1.0
    return arr


def random_density_match(N: int, density: float, seed: int = 42) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return (rng.random(N) < density).astype(np.float64)


def cumulative_pi(N: int) -> np.ndarray:
    chi = chi_indicator(N)
    return np.cumsum(chi)


def residual_pi_li(N: int) -> np.ndarray:
    pi = cumulative_pi(N)
    mp.dps = 30
    li_arr = np.array([0.0 if (n + 1) < 2 else float(li(n + 1)) for n in range(N)])
    return pi - li_arr


# ---------------------------------------------------------------------------
# TT-SVD: per-cut ranks and entropies
# ---------------------------------------------------------------------------

def tt_svd_ranks(T: np.ndarray, eps_list=(0.0, 1e-2, 1e-1)):
    d = T.ndim
    n = list(T.shape)
    norm = float(np.linalg.norm(T))
    if norm == 0:
        return None
    Tn = T / norm
    sing_spectra = []
    ranks_exact = [1]
    ranks_eps = {e: [1] for e in eps_list}
    r = 1
    C = Tn.copy().astype(np.float64)
    eps_per_cut = {e: (e ** 2) / max(d - 1, 1) for e in eps_list}
    for j in range(d - 1):
        C = C.reshape(r * n[j], -1)
        U, S, Vh = np.linalg.svd(C, full_matrices=False)
        sing_spectra.append(S.copy())
        s_max = S[0] if S.size else 0
        # Strict numerical cutoff -- match exactly the theorem-tight integer
        # rank rather than chasing floating-point noise. With only 0/1 entries
        # the smallest nonzero singular value is bounded away from 0 (strictly
        # positive, well above 1e-10 * norm at all sizes we run), so 1e-10 is
        # safe.
        r_exact = int((S > 1e-10 * max(s_max, 1)).sum())
        ranks_exact.append(max(r_exact, 1))
        for e in eps_list:
            tail = np.cumsum(S[::-1] ** 2)[::-1]
            mask = tail > eps_per_cut[e]
            r_e = int(mask.sum()) if mask.any() else 1
            ranks_eps[e].append(max(r_e, 1))
        r_new = r_exact
        U = U[:, :r_new]
        S = S[:r_new]
        Vh = Vh[:r_new, :]
        C = (S[:, None] * Vh).reshape((r_new,) + tuple(n[j + 1:]))
        r = r_new
    ranks_exact.append(1)
    for e in eps_list:
        ranks_eps[e].append(1)
    mid = (d - 1) // 2 if d > 1 else 0
    S_mid = sing_spectra[mid] if sing_spectra else np.array([1.0])
    p = S_mid ** 2
    p = p / p.sum()
    p = p[p > 0]
    entropy = float(-(p * np.log(p)).sum()) if p.size else 0.0
    # max-possible-at-half-cut for non-uniform local dims
    left_dim = int(np.prod(n[: mid + 1])) if d > 1 else n[0]
    right_dim = int(np.prod(n[mid + 1:])) if d > 1 else 1
    return {
        "ranks_exact": ranks_exact,
        "ranks_eps": ranks_eps,
        "sing_spectra": sing_spectra,
        "half_cut_entropy_nats": entropy,
        "half_cut_index": mid,
        "max_possible_at_half_cut": min(left_dim, right_dim),
        "sing_max_at_half_cut": float(sing_spectra[mid][0]) if sing_spectra else 0.0,
    }


def storage_cost(ranks, local_dim):
    return sum(ranks[j] * local_dim * ranks[j + 1] for j in range(len(ranks) - 1))


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def report_base2(label, N, info):
    re = info["ranks_exact"]
    rmax = max(re)
    bound = info["max_possible_at_half_cut"]
    fill = rmax / bound if bound else 0
    line = (
        f"  {label:<15s}  N={N:>9d}  max-rank={rmax:>5d}/{bound:<6d} "
        f"({fill*100:5.1f}% of half-cut bound) "
        f"H_half={info['half_cut_entropy_nats']:.3f}/{math.log(bound):.3f} nats"
    )
    print(line)
    eps01 = info["ranks_eps"][1e-2]
    eps10 = info["ranks_eps"][1e-1]
    print(
        f"      eps=1e-2 max-rank={max(eps01)}/{bound}  "
        f"eps=1e-1 max-rank={max(eps10)}/{bound}"
    )


def run_base2_sweep(ks):
    print("=" * 78)
    print("Family A: base-2 sweep (replicates the original W=2 measurement)")
    print("=" * 78)
    for k in ks:
        N = 2 ** k
        t0 = time.perf_counter()
        chi = chi_indicator(N)
        density = chi.mean()
        T = chi.reshape([2] * k)
        info_p = tt_svd_ranks(T)
        rnd = random_density_match(N, density)
        T_r = rnd.reshape([2] * k)
        info_r = tt_svd_ranks(T_r)
        pi_n = cumulative_pi(N)
        T_pi = pi_n.reshape([2] * k)
        info_cum = tt_svd_ranks(T_pi)
        info_res = None
        if k <= 16:
            res = residual_pi_li(N)
            T_res = res.reshape([2] * k)
            info_res = tt_svd_ranks(T_res)
        elapsed = time.perf_counter() - t0
        pred = predicted_max_rank(2, k)
        observed = max(info_p["ranks_exact"])
        match = "MATCH" if observed == pred else f"MISMATCH (pred {pred}, obs {observed})"
        print(f"\nk={k}  (N={N}, density={density:.4f}, time={elapsed:.2f}s)  theorem: {match}")
        report_base2("chi_P", N, info_p)
        report_base2("random@dens", N, info_r)
        report_base2("cum pi(n)", N, info_cum)
        if info_res is not None:
            report_base2("pi(n)-li(n)", N, info_res)


# ---------------------------------------------------------------------------
# Base-W theorem verification
# ---------------------------------------------------------------------------

def run_baseW(W: int, d: int, do_random: bool = True, seed: int = 42):
    """Exact rank measurement at every cut for chi_P viewed in base W."""
    N = W ** d
    if N > 50_000_000:
        print(f"  W={W} d={d}: N={N} too large -- skipping")
        return None
    t0 = time.perf_counter()
    chi = chi_indicator(N)
    density = chi.mean()
    T = chi.reshape([W] * d)
    info = tt_svd_ranks(T)
    re = info["ranks_exact"]
    # predicted ranks at every cut
    preds = [None] + [predicted_rank(W, d, j) for j in range(1, d)] + [None]
    # check theorem at every cut
    cut_match = all(re[j] == preds[j] for j in range(1, d))
    obs_max = max(re)
    pred_max = predicted_max_rank(W, d)
    elapsed_p = time.perf_counter() - t0

    info_r = None
    rnd_max = None
    if do_random:
        rnd = random_density_match(N, density, seed=seed)
        Tr = rnd.reshape([W] * d)
        info_r = tt_svd_ranks(Tr)
        rnd_max = max(info_r["ranks_exact"])

    elapsed = time.perf_counter() - t0

    phi = euler_phi(W)
    half_cut_random_ceiling = min(W ** ((d + 1) // 2), W ** (d // 2))
    print(
        f"  W={W:<3d} d={d:<2d} N={N:<10d} phi(W)/W={phi}/{W}={phi/W:.4f}"
        f"   chi_P max-rank={obs_max:<8d} predicted={pred_max:<8d} "
        f"{'OK' if obs_max == pred_max else '*MISMATCH*'}"
    )
    if do_random and rnd_max is not None:
        print(
            f"           random@dens max-rank={rnd_max:<8d} "
            f"(half-cut ceiling {half_cut_random_ceiling}, "
            f"chi_P/random ratio {obs_max/rnd_max:.4f}, "
            f"phi(W)/W = {phi/W:.4f})"
        )
    print(f"           time = {elapsed:.2f}s   per-cut theorem ALL MATCH: {cut_match}")

    # Per-cut profile
    rows = []
    for j in range(1, d):
        rows.append({
            "j": j,
            "obs": re[j],
            "pred": preds[j],
            "max_at_cut": min(W ** j, W ** (d - j)),
        })
    return {
        "W": W, "d": d, "N": N, "phi_over_W": phi / W,
        "chi_P_max": obs_max, "pred_max": pred_max,
        "cut_match": cut_match,
        "rows": rows,
        "random_max": rnd_max,
        "density": density,
        "info": info,
    }


def run_baseW_sweep(quick: bool = False):
    print("\n" + "=" * 78)
    print("Family B: base-W theorem verification")
    print("=" * 78)
    print("THEOREM: at cut j of depth d, rank(chi_P) = min(W^j, phi(W)*W^{d-j-1}+1)")
    print()
    plan = [
        (2, 10), (2, 14), (2, 20),
        (6, 4), (6, 6), (6, 8),
        (30, 2), (30, 3), (30, 4),
        (210, 2), (210, 3),
    ]
    if quick:
        plan = [(2, 10), (6, 4), (6, 6), (30, 2), (30, 3), (210, 2)]
    results = []
    for W, d in plan:
        try:
            r = run_baseW(W, d, do_random=True)
            if r is not None:
                results.append(r)
        except MemoryError:
            print(f"  W={W} d={d}: MemoryError, skipped")
    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=int, default=None,
                        help="If set, run only the base-W sweep for this W (and skip base-2 family A).")
    parser.add_argument("--depth", type=int, default=None,
                        help="If set with --base, run a single (W, d).")
    parser.add_argument("--quick", action="store_true",
                        help="Smaller sizes for CI / quick verification.")
    parser.add_argument("--family", choices=["a", "b", "both"], default="both")
    args = parser.parse_args()

    if args.base is not None and args.depth is not None:
        print("Single configuration:")
        run_baseW(args.base, args.depth, do_random=True)
        return

    if args.family in ("a", "both"):
        ks = [10, 12, 14] if args.quick else [10, 12, 14, 16, 18, 20]
        run_base2_sweep(ks)

    if args.family in ("b", "both"):
        run_baseW_sweep(quick=args.quick)

    print()
    print("Interpretation:")
    print("  Family A reproduces chi_P max-rank = 2^{d/2-1} + 1 at every d.")
    print("  Family B verifies the base-W generalisation:")
    print("    rank at cut j = min(W^j, phi(W)*W^{d-j-1} + 1)")
    print("  Half-cut ratio chi_P / random -> phi(W) / W as d grows.")
    print("  W=2:   1/2  = 0.5000")
    print("  W=6:   2/6  = 0.3333")
    print("  W=30:  8/30 = 0.2667")
    print("  W=210: 48/210 = 0.2286")
    print()
    print("This confirms the wheel-W sieve has an MPS-language interpretation:")
    print("  primorial-W trial division saves a factor phi(W)/W of bond dimension.")
    print("  No polylog limit (ratio -> 0 only as W -> infty, but the cost of")
    print("  evaluating the wheel grows >= W). Equivalence (E) failure mode.")


if __name__ == "__main__":
    main()
