#!/usr/bin/env python3
"""C2 — Free cumulants of chi_P operator from MPS unfolding.

Composes E2.1 (MPS bond-dim identity at every primorial cut) with the
free probability framework partially explored in S53
(experiments/wildcard/free_probability_delta.py).

Construction. View chi_P : [1, W^d] -> {0,1} as a vector in (C^W)^{otimes d}
via the natural base-W reshape. For a cut 1 <= j < d the unfolding
M^(j) is a rectangular matrix of shape (W^j, W^{d-j}). Treat the
spectrum of M^(j) M^(j)^T as the spectral distribution of a positive
self-adjoint operator A_chi^pos in a finite tracial von Neumann algebra.
Compute the first 4 free cumulants and compare against the canonical
free distributions (Marchenko-Pastur / free-Poisson, Wigner semicircle).

Cross-domain source. Voiculescu free-probability framework. Reference:
Mingo & Speicher, "Free Probability and Random Matrices" (Fields Inst.
Monographs 35, 2017), Ch. 1-2 for free cumulants and Ch. 4 for MP/free-
Poisson relation. Hiai & Petz, "Semicircle Law, Free Random Variables
and Entropy" for moment-cumulant transforms on non-crossing partitions.

Falsification (pre-stated):
- PR1: For W=2 at half-cut, the standardized free cumulants
  (kappa_1, kappa_2, kappa_3, kappa_4) of chi_P's full sv2 spectrum match
  MP(c = phi(W)/W = 0.5) within 10% relative for all four.
- PR2: For random Bernoulli baseline of matched density and shape,
  the same cumulants match MP(c = density-corrected ratio) within 10%.
- PR3: Across W in {2, 3, 5, 6, 30}, the chi_P free cumulants under the
  "drop zeros" normalization are constant in r at MP-rate cumulants
  match phi(W)/W (free-Poisson signature of the wheel-W sieve).

Outcomes measured under three normalizations:
  A_full:   include rank-deficit zeros, total dim min(W^j, W^{d-j}).
  B_nonzero: drop zeros, normalize over rank-many sv2.
  C_bulk:   drop zeros and the leading (rank-1 mean) sv2.
"""

import argparse
import json
import time
from math import comb
from pathlib import Path

import numpy as np


def sieve(N):
    is_p = np.ones(N + 1, dtype=bool)
    is_p[:2] = False
    for i in range(2, int(N**0.5) + 1):
        if is_p[i]:
            is_p[i * i :: i] = False
    return is_p


def chi_P_vec(N):
    """chi_P_vec[k] = 1 if (k+1) is prime, 0 otherwise; length N (so domain [1, N])."""
    return sieve(N)[1 : N + 1].astype(np.float64)


def euler_phi(n):
    r = n
    p = 2
    nn = n
    while p * p <= nn:
        if nn % p == 0:
            while nn % p == 0:
                nn //= p
            r -= r // p
        p += 1
    if nn > 1:
        r -= r // nn
    return r


def free_cumulants(moments):
    """First 4 free cumulants from raw moments m_1..m_4."""
    m1, m2, m3, m4 = moments
    return np.array(
        [
            m1,
            m2 - m1**2,
            m3 - 3 * m1 * m2 + 2 * m1**3,
            m4 - 4 * m1 * m3 - 2 * m2**2 + 10 * m1**2 * m2 - 5 * m1**4,
        ]
    )


def classical_cumulants(moments):
    """First 4 classical cumulants from raw moments m_1..m_4."""
    m1, m2, m3, m4 = moments
    return np.array(
        [
            m1,
            m2 - m1**2,
            m3 - 3 * m1 * m2 + 2 * m1**3,
            m4 - 4 * m1 * m3 - 3 * m2**2 + 12 * m1**2 * m2 - 6 * m1**4,
        ]
    )


def cumulants_from_sv2(sv2, m_total):
    """Return (moments, free_cumulants, classical_cumulants) under standardization mean=1.

    sv2 : array of (squared singular) values, possibly with zeros dropped.
    m_total : total number of positions (>= len(sv2)). Zeros are absorbed implicitly.
    """
    if len(sv2) == 0:
        return (np.array([np.nan] * 4),) * 3
    s = sv2.sum()
    if s == 0:
        return (np.array([np.nan] * 4),) * 3
    mean = s / m_total
    moments = np.array([np.sum(sv2**k) / m_total / mean**k for k in range(1, 5)])
    return moments, free_cumulants(moments), classical_cumulants(moments)


def mp_moments(c, kmax=4):
    """Moments of MP(c) in the convention with mean = c."""
    out = []
    for K in range(1, kmax + 1):
        s = 0.0
        for r in range(K):
            s += comb(K, r) * comb(K - 1, r) * c ** (r + 1) / (r + 1)
        out.append(s)
    return np.array(out)


def mp_normalized_free_cumulants(c, kmax=4):
    """Free cumulants of MP(c) under standard random-matrix-theory convention.

    For X iid m x n (m > n) with variance 1 and c = n/m, eigenvalues of (1/m) X X^T
    have mean 1 and free cumulants kappa_r = c^{r-1}. So kappa_1 = 1, kappa_2 = c,
    kappa_3 = c^2, kappa_4 = c^3. Verified empirically against random Gaussian.
    """
    return np.array([c ** (r - 1) for r in range(1, kmax + 1)])


def mp_normalized_moments(c, kmax=4):
    """Moments of MP(c) after rescaling to mean=1."""
    return mp_moments(c, kmax) / c ** np.arange(1, kmax + 1)


def squared_singular_values(M):
    """Sorted descending squared singular values."""
    s = np.linalg.svd(M, compute_uv=False)
    return (s**2)[::-1][::-1]  # ensure descending sort


DROP_LEVELS = [0, 1, 2, 5, 10, 20, 50]


def analyze_matrix(M, label, threshold_factor=1e-10, theory_rank=None,
                   drop_top_phi=None):
    """Compute free cumulants under several drop levels.

    A_full        : full spectrum including rank-deficit zeros.
    B_nonzero     : drop zeros, normalize over rank-many sv2.
    Drop_<k>      : drop zeros and the leading k sv2 (k in DROP_LEVELS).
    """
    sv2 = squared_singular_values(M)
    m, n = M.shape
    cap = min(m, n)
    smax = sv2.max() if sv2.size else 0.0
    threshold = threshold_factor * smax if smax > 0 else 0.0

    if theory_rank is not None:
        sv2_nz = sv2[:theory_rank]
    else:
        sv2_nz = sv2[sv2 > threshold]

    out = {
        "label": label,
        "shape": [m, n],
        "rank_numeric": int(np.sum(sv2 > threshold)),
        "rank_used": int(len(sv2_nz)),
        "leading_sv2": float(sv2[0]) if sv2.size else float("nan"),
        "sv2_top10": sv2[: min(10, sv2.size)].tolist(),
        "sv2_min_nz": float(sv2_nz.min()) if len(sv2_nz) else float("nan"),
        "drop_levels": DROP_LEVELS,
    }

    moments, fc, cc = cumulants_from_sv2(sv2_nz, cap)
    out["A_full"] = {"moments": moments.tolist(), "free_cumulants": fc.tolist(),
                     "classical_cumulants": cc.tolist()}
    moments, fc, cc = cumulants_from_sv2(sv2_nz, len(sv2_nz))
    out["B_nonzero"] = {"moments": moments.tolist(), "free_cumulants": fc.tolist(),
                        "classical_cumulants": cc.tolist()}

    for k in DROP_LEVELS:
        if len(sv2_nz) > k:
            sub = sv2_nz[k:]
            moments, fc, cc = cumulants_from_sv2(sub, len(sub))
        else:
            moments, fc, cc = (np.array([np.nan]*4),) * 3
        out[f"drop_{k}"] = {"moments": moments.tolist(),
                            "free_cumulants": fc.tolist(),
                            "classical_cumulants": cc.tolist(),
                            "n_used": int(max(0, len(sv2_nz) - k))}
    return out


def adaptive_outlier_count(M, theory_rank, c_pred, tol=0.1, kmax=None):
    """Smallest k such that drop_k yields free cumulants within tol of MP(c_pred)."""
    sv2 = squared_singular_values(M)
    sv2_nz = sv2[:theory_rank]
    pred = mp_normalized_free_cumulants(c_pred)
    if kmax is None:
        kmax = max(1, len(sv2_nz) - 5)
    for k in range(0, kmax):
        sub = sv2_nz[k:]
        if len(sub) < 5:
            return None
        moments, fc, _ = cumulants_from_sv2(sub, len(sub))
        # Compare kappa_2, kappa_3, kappa_4 to MP prediction.
        rel = np.array([abs(fc[i] - pred[i]) / max(abs(pred[i]), 1e-3) for i in (1, 2, 3)])
        if np.all(rel < tol):
            return k
    return None


def analyze_chi_p(W, d, j=None):
    if j is None:
        j = d // 2
    N = W**d
    chi = chi_P_vec(N)
    M = chi.reshape(W**j, W ** (d - j))
    phiW = euler_phi(W)
    theory_rank = min(W**j, phiW * W ** (d - j - 1) + 1)
    info = analyze_matrix(M, label=f"chi_P W={W} d={d} j={j}",
                          theory_rank=theory_rank)
    info["W"] = W
    info["d"] = d
    info["j"] = j
    info["N"] = N
    info["theory_rank"] = theory_rank
    info["pi_N"] = int(chi.sum())
    info["density"] = info["pi_N"] / N
    info["phi_W_over_W"] = phiW / W
    info["mp_c_active"] = phiW / W
    info["mp_kappa_active_pred"] = mp_normalized_free_cumulants(phiW / W).tolist()
    info["outlier_count_tol_0p1"] = adaptive_outlier_count(M, theory_rank, phiW / W, tol=0.1)
    info["outlier_count_tol_0p2"] = adaptive_outlier_count(M, theory_rank, phiW / W, tol=0.2)
    info["outlier_count_tol_0p5"] = adaptive_outlier_count(M, theory_rank, phiW / W, tol=0.5)
    return info


def _accum_keys():
    return ["A_full", "B_nonzero"] + [f"drop_{k}" for k in DROP_LEVELS]


def analyze_bernoulli(shape, density, seed=0, n_trials=3, theory_rank=None):
    rng = np.random.default_rng(seed)
    accum_keys = _accum_keys()
    sums = {k: {"moments": np.zeros(4), "free_cumulants": np.zeros(4), "classical_cumulants": np.zeros(4)} for k in accum_keys}
    for _ in range(n_trials):
        M = rng.binomial(1, density, size=shape).astype(np.float64)
        info = analyze_matrix(M, label="bernoulli", theory_rank=theory_rank)
        for k in accum_keys:
            for q in ("moments", "free_cumulants", "classical_cumulants"):
                sums[k][q] += np.array(info[k][q])
    out = {"shape": list(shape), "density": density, "n_trials": n_trials}
    for k in accum_keys:
        out[k] = {q: (v / n_trials).tolist() for q, v in sums[k].items()}
    return out


def analyze_active_bernoulli(W, d, density, seed=0, n_trials=3, j=None):
    """Bernoulli baseline matched to the ACTIVE submatrix shape."""
    if j is None:
        j = d // 2
    phiW = euler_phi(W)
    m_act = W**j - 1
    n_act = phiW * W ** (d - j - 1)
    if m_act <= 0 or n_act <= 0:
        return None
    rng = np.random.default_rng(seed)
    accum_keys = _accum_keys()
    sums = {k: {"moments": np.zeros(4), "free_cumulants": np.zeros(4), "classical_cumulants": np.zeros(4)} for k in accum_keys}
    for _ in range(n_trials):
        M = rng.binomial(1, density, size=(m_act, n_act)).astype(np.float64)
        info = analyze_matrix(M, label="active_bernoulli", theory_rank=None)
        for k in accum_keys:
            for q in ("moments", "free_cumulants", "classical_cumulants"):
                sums[k][q] += np.array(info[k][q])
    out = {"shape": [m_act, n_act], "density": density, "n_trials": n_trials}
    for k in accum_keys:
        out[k] = {q: (v / n_trials).tolist() for q, v in sums[k].items()}
    return out


def fmt(arr, prec=4):
    return "(" + ", ".join(f"{x:+.{prec}f}" for x in arr) + ")"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ws", type=str, default="2,3,5,6,30",
                        help="comma-separated W values")
    parser.add_argument("--dmax_w2", type=int, default=18)
    parser.add_argument("--dmax_w3", type=int, default=12)
    parser.add_argument("--dmax_w5", type=int, default=8)
    parser.add_argument("--dmax_w6", type=int, default=8)
    parser.add_argument("--dmax_w30", type=int, default=4)
    parser.add_argument("--with_bernoulli", action="store_true")
    parser.add_argument("--out", type=str, default="free_cumulants_chi_p_results.json")
    parser.add_argument("--quick", action="store_true",
                        help="smallest configs only")
    args = parser.parse_args()

    Ws = [int(x) for x in args.ws.split(",")]
    dmax = {2: args.dmax_w2, 3: args.dmax_w3, 5: args.dmax_w5, 6: args.dmax_w6, 30: args.dmax_w30}

    configs = []
    for W in Ws:
        if W == 2:
            d_iter = range(8, dmax[W] + 1, 2) if not args.quick else [10, 12]
        elif W == 3:
            d_iter = range(6, dmax[W] + 1, 2) if not args.quick else [6, 8]
        elif W == 5:
            d_iter = range(4, dmax[W] + 1, 2) if not args.quick else [4, 6]
        elif W == 6:
            d_iter = range(4, dmax[W] + 1, 2) if not args.quick else [4, 6]
        elif W == 30:
            d_iter = range(2, dmax[W] + 1) if not args.quick else [2, 3]
        else:
            d_iter = range(2, 6)
        for d in d_iter:
            if W**d > 6e7:
                continue
            configs.append((W, d))

    print(f"# Free cumulants of chi_P operator (C2)")
    print(f"# {len(configs)} configurations")
    print(f"# Cross-domain reference: Mingo-Speicher 2017, free probability.")
    print()

    results = []
    for W, d in configs:
        t0 = time.time()
        try:
            info = analyze_chi_p(W, d)
        except MemoryError:
            print(f"# W={W} d={d}: MemoryError, skipping")
            continue
        elapsed = time.time() - t0
        info["elapsed_sec_chi"] = elapsed

        if args.with_bernoulli:
            t1 = time.time()
            phiW = euler_phi(W)
            info["bernoulli"] = analyze_bernoulli(
                shape=tuple(info["shape"]),
                density=info["density"],
                seed=hash((W, d)) & 0xFFFF,
                n_trials=3,
                theory_rank=None,
            )
            # Active-submatrix density: primality density on coprime residues only.
            phi_density = info["pi_N"] / (info["N"] * phiW / W) if phiW > 0 else 0.0
            info["active_bernoulli"] = analyze_active_bernoulli(
                W=W, d=d, density=phi_density,
                seed=(hash((W, d)) ^ 0x12345) & 0xFFFF,
                n_trials=3,
            )
            info["elapsed_sec_bern"] = time.time() - t1

        # console summary
        print(f"W={W:2d} d={d:2d}  N={info['N']:>10d}  shape={info['shape']}  "
              f"rank={info['theory_rank']}/{info['rank_numeric']}  "
              f"density={info['density']:.6f}  pi(N)={info['pi_N']}")
        c_pred = info["mp_c_active"]
        print(f"   MP(c={c_pred:.4f}) prediction kappa = {fmt(info['mp_kappa_active_pred'])}")
        print(f"   outlier counts (k s.t. drop_k matches MP within rel-tol): "
              f"tol=0.1 -> {info['outlier_count_tol_0p1']}  "
              f"tol=0.2 -> {info['outlier_count_tol_0p2']}  "
              f"tol=0.5 -> {info['outlier_count_tol_0p5']}")
        for k in DROP_LEVELS:
            v = info[f"drop_{k}"]
            line = f"   chi_P drop_{k:<2d} (n={v.get('n_used','?'):>4})  kappa={fmt(v['free_cumulants'])}"
            if args.with_bernoulli:
                ab = info.get("active_bernoulli")
                if ab is not None:
                    av = ab[f"drop_{k}"]
                    line += f"   ActB kappa={fmt(av['free_cumulants'])}"
            print(line)

        print(f"   elapsed: chi={elapsed:.2f}s"
              + (f" bern={info.get('elapsed_sec_bern', 0):.2f}s" if args.with_bernoulli else ""))
        print()
        results.append(info)

    out_path = Path(__file__).parent / args.out
    out_path.write_text(json.dumps(results, indent=2))
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
