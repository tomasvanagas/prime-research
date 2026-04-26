"""
Construction C5 (composing E1.4 + E2.5):
    Test the "N/2 universality" of pi(x) mod 2 against a class of natural
    number-theoretic Boolean functions.

E1.4 records that for the prime indicator chi_P : {0,...,2^N - 1} -> {0,1},
six independent complexity measures all land at the boundary N/2:

    * approximate degree at eps=0.49 :  ceil(N/2)
    * communication matrix rank      :  2^{N/2-1} + 2
    * PTF degree                     :  N/2
    * GF(2) LFSR linear complexity   :  ~ 2^{N-1} on length-2^N sequence
    * per-bit R-correlation crossover:  bit ~ N/2
    * 3-way real tensor rank         :  ~ 2^{N/2}

E1.4 was measured on chi_P only. This experiment runs the cheap subset of
the battery on five Boolean functions:

    f_pi      : f(n) = 1 iff n is prime           (baseline)
    f_sqfree  : f(n) = 1 iff n is squarefree
    f_mu_pos  : f(n) = 1 iff Mobius(n) = +1
    f_lam_pos : f(n) = 1 iff Liouville(n) = +1
    f_sqfree3 : f(n) = 1 iff n is squarefree AND n mod 3 == 1
    f_aes     : pseudorandom control: AES counter-mode output bit

Measurements per function (small N because LP is exponential):

    M1 communication-matrix rank       N in {6, 8, 10, 12, 14}
    M2 GF(2) Berlekamp-Massey LFSR     N in {10, 12, 14}
    M3 approximate degree at eps=0.49  N in {4, 6, 8, 10}
    M4 PTF degree                      N in {4, 6, 8}

Falsification statements:

    PR1  M1: comm-rank universal class
        For every f and every N, comm-matrix rank of M[a,b] = f(a*2^k + b)
        is in {2^{N/2-1}, 2^{N/2-1}+1, 2^{N/2-1}+2, 2^{N/2}}.  In particular,
        all five non-pi functions deviate from f_pi's "+2" anomaly: they
        either saturate at 2^{N/2} (no anomaly) or exhibit their own +c
        anomaly tied to the function's own monotonicity / parity structure.

    PR2  M2: BM-LFSR linear complexity ratio L/length is in [0.4, 0.6]
        for all five functions on length-2^N sequences (N >= 10).

    PR3  M3: approximate degree at eps=0.49 equals ceil(N/2) for every
        Boolean function in the class with density bounded away from 0,1.
        I.e., the N/2 universality holds for the WHOLE class, not just chi_P.

    PR4  M4: PTF degree equals ceil(N/2) for every function in the class.

Outcome interpretation:
    * PR1-PR4 all PASS for all five functions => N/2 universality is a
      META-THEOREM of natural NT Boolean functions, not a property of
      chi_P specifically. New negative-shape edge candidate.
    * Some function deviates => chi_P is special at that measure;
      records the deviation as a structural property unique to primes.

Save under: experiments/constructions/n_over_2_universality_class/
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from itertools import combinations

import numpy as np
from scipy.optimize import linprog


# -----------------------------------------------------------------------
# Boolean function generators
# -----------------------------------------------------------------------


def sieve_smallest_prime_factor(N: int) -> np.ndarray:
    """Return spf[0..N] : smallest prime factor of n for 2 <= n <= N."""
    spf = np.zeros(N + 1, dtype=np.int32)
    for p in range(2, N + 1):
        if spf[p] == 0:
            slc = spf[p::p]
            slc[slc == 0] = p
            spf[p::p] = slc
    return spf


def function_truth_tables(N: int) -> dict:
    """Compute Boolean truth tables on {0, 1, ..., 2^N - 1}.

    Returns a dict of np.float64 arrays of length 2^N (entries in {0, 1})."""
    M = 2 ** N
    spf = sieve_smallest_prime_factor(M - 1)

    # Omega(n) = number of prime factors with multiplicity
    Omega = np.zeros(M, dtype=np.int32)
    is_squarefree = np.ones(M, dtype=bool)
    is_squarefree[0] = False
    for n in range(2, M):
        m = n
        last_p = -1
        repeated = False
        omega = 0
        while m > 1:
            p = spf[m]
            if p == last_p:
                repeated = True
            omega += 1
            last_p = p
            m //= p
        Omega[n] = omega
        if repeated:
            is_squarefree[n] = False
    is_squarefree[1] = True
    Omega[0] = 0
    Omega[1] = 0

    is_prime = np.zeros(M, dtype=bool)
    for n in range(2, M):
        is_prime[n] = (spf[n] == n)

    # Mobius(n): 0 if not squarefree, else (-1)^Omega(n)
    mobius = np.zeros(M, dtype=np.int32)
    for n in range(1, M):
        if is_squarefree[n]:
            mobius[n] = 1 if Omega[n] % 2 == 0 else -1
    mobius[1] = 1  # mu(1) = 1

    # Liouville(n) = (-1)^Omega(n)
    liouville = np.where(Omega % 2 == 0, 1, -1).astype(np.int32)
    liouville[0] = 0  # n = 0 excluded

    # Boolean functions
    f_pi = is_prime.astype(np.float64)
    f_sqfree = is_squarefree.astype(np.float64)
    f_mu_pos = (mobius == 1).astype(np.float64)
    f_lam_pos = (liouville == 1).astype(np.float64)
    n_arr = np.arange(M)
    f_sqfree3 = (is_squarefree & (n_arr % 3 == 1)).astype(np.float64)

    # Pseudorandom control: AES-CTR(seed=0xc5) bits, sampled at the same
    # density as f_pi for fair comparison. We just use SHA256 for ease of
    # implementation - good enough as a TC^0 / cryptographic PRNG witness.
    def prf_bit(n):
        h = hashlib.sha256(b"c5_universality_class:" + n.to_bytes(8, "big")).digest()
        return h[0] & 1

    f_prf = np.array([prf_bit(n) for n in range(M)], dtype=np.float64)

    # Density-matched PRF: pick the lowest-hash MSBs to give density ~ pi(M)/M
    # so that approximate-degree comparisons are at the same density as f_pi.
    target_density = float(f_pi.mean())
    if target_density > 0:
        hash_vals = np.array([
            int.from_bytes(
                hashlib.sha256(b"density:" + n.to_bytes(8, "big")).digest()[:4],
                "big",
            )
            for n in range(M)
        ], dtype=np.uint32)
        threshold = np.quantile(hash_vals, target_density)
        f_prf_matched = (hash_vals < threshold).astype(np.float64)
    else:
        f_prf_matched = np.zeros(M, dtype=np.float64)

    return {
        "f_pi": f_pi,
        "f_sqfree": f_sqfree,
        "f_mu_pos": f_mu_pos,
        "f_lam_pos": f_lam_pos,
        "f_sqfree3": f_sqfree3,
        "f_prf": f_prf,
        "f_prf_matched": f_prf_matched,
    }


# -----------------------------------------------------------------------
# M1: Communication-matrix rank
# -----------------------------------------------------------------------


def comm_matrix_rank(f: np.ndarray, N: int) -> tuple[int, int, int]:
    """Reshape f into a balanced-split communication matrix and return
    (rank, n_rows, n_cols) using exact rational rank via SVD threshold.
    Balanced split: M[a, b] = f(a * 2^{N-k} + b), k = N // 2."""
    k = N // 2
    rows = 2 ** k
    cols = 2 ** (N - k)
    M_mat = f.reshape(rows, cols)
    # Exact rank: count singular values above 1e-10 * max
    s = np.linalg.svd(M_mat, compute_uv=False)
    if s[0] == 0:
        return 0, rows, cols
    rank = int(np.sum(s > 1e-9 * s[0]))
    return rank, rows, cols


# -----------------------------------------------------------------------
# M2: Berlekamp-Massey GF(2) linear complexity
# -----------------------------------------------------------------------


def berlekamp_massey_gf2(seq: np.ndarray) -> int:
    """Linear complexity of binary sequence over GF(2)."""
    n = len(seq)
    s = [int(x) & 1 for x in seq]

    c = [1]
    b = [1]
    L = 0
    m = 1

    for i in range(n):
        d = s[i]
        for j in range(1, L + 1):
            if j < len(c) and c[j]:
                d ^= s[i - j]

        if d == 0:
            m += 1
        elif 2 * L <= i:
            t = c[:]
            while len(c) < len(b) + m:
                c.append(0)
            for k in range(len(b)):
                c[k + m] ^= b[k]
            L = i + 1 - L
            b = t
            m = 1
        else:
            while len(c) < len(b) + m:
                c.append(0)
            for k in range(len(b)):
                c[k + m] ^= b[k]
            m += 1
    return L


# -----------------------------------------------------------------------
# M3: Approximate degree LP at eps = 0.49
# -----------------------------------------------------------------------


def monomial_matrix(N: int, max_degree: int) -> np.ndarray:
    """Build matrix where row x is the bitmask values prod_{i in S} x_i for
    each S of size <= max_degree.  x_i is bit i of x."""
    n_points = 2 ** N
    masks = []
    for d in range(max_degree + 1):
        for S in combinations(range(N), d):
            mask = 0
            for bit in S:
                mask |= 1 << bit
            masks.append(mask)
    n_mono = len(masks)
    M_mat = np.zeros((n_points, n_mono), dtype=np.float64)
    for j, mask in enumerate(masks):
        for x in range(n_points):
            if (x & mask) == mask:
                M_mat[x, j] = 1.0
    return M_mat


def n_monomials(N: int, d: int) -> int:
    return sum(math.comb(N, k) for k in range(d + 1))


def min_approx_error_lp(M_mat: np.ndarray, f: np.ndarray, timeout: float = 120.0):
    n_points, n_mono = M_mat.shape
    n_vars = n_mono + 1

    c_obj = np.zeros(n_vars)
    c_obj[-1] = 1.0

    ones_col = -np.ones((n_points, 1))
    A_upper = np.hstack([M_mat, ones_col])
    A_lower = np.hstack([-M_mat, ones_col])
    A_ub = np.vstack([A_upper, A_lower])
    b_ub = np.concatenate([f, -f])

    bounds = [(None, None)] * n_mono + [(0, None)]

    result = linprog(
        c_obj,
        A_ub=A_ub,
        b_ub=b_ub,
        bounds=bounds,
        method="highs",
        options={"presolve": True, "time_limit": timeout},
    )
    if result.success:
        return float(result.x[-1])
    return None


def approx_degree_eps49(f: np.ndarray, N: int, max_mono: int = 4500):
    """Return the smallest d such that min eps approximating f by a degree-d
    polynomial is < 0.49.  None if not found within max_mono budget."""
    eps_by_d = {}
    for d in range(0, N + 1):
        nm = n_monomials(N, d)
        if nm > max_mono:
            return None, eps_by_d
        M_mat = monomial_matrix(N, d)
        eps = min_approx_error_lp(M_mat, f)
        if eps is None:
            return None, eps_by_d
        eps_by_d[d] = eps
        if eps < 0.49:
            return d, eps_by_d
    return None, eps_by_d


# -----------------------------------------------------------------------
# M4: PTF degree (sign-representation)
# -----------------------------------------------------------------------


def ptf_degree(f: np.ndarray, N: int, max_mono: int = 4500):
    """Return the smallest d such that there exists a polynomial p of degree
    <= d with p(x) > 0 iff f(x) = 1 and p(x) < 0 iff f(x) = 0, in {0,1}^N.

    LP: find c such that  M @ c  has the right sign with margin >= 1.
    For x with f=1:  M[x] @ c  >= 1  ->   - M[x] @ c  <= -1
    For x with f=0:  M[x] @ c  <= -1 ->     M[x] @ c  <= -1
    No objective; feasibility check.
    """
    n_points = 2 ** N
    pos = f > 0.5
    neg = ~pos

    for d in range(0, N + 1):
        nm = n_monomials(N, d)
        if nm > max_mono:
            return None
        M_mat = monomial_matrix(N, d)
        # Build A_ub * c <= b_ub with margin
        A_pos = -M_mat[pos]  # need M_mat[pos] @ c >= 1
        b_pos = -np.ones(np.sum(pos))
        A_neg = M_mat[neg]  # need M_mat[neg] @ c <= -1
        b_neg = -np.ones(np.sum(neg))
        A_ub = np.vstack([A_pos, A_neg])
        b_ub = np.concatenate([b_pos, b_neg])
        # Feasibility: minimize 0
        c_obj = np.zeros(nm)
        bounds = [(None, None)] * nm
        result = linprog(
            c_obj,
            A_ub=A_ub,
            b_ub=b_ub,
            bounds=bounds,
            method="highs",
            options={"presolve": True, "time_limit": 60},
        )
        if result.status == 0:
            return d
    return None


# -----------------------------------------------------------------------
# Driver
# -----------------------------------------------------------------------


def run_battery(maxN_m1: int, maxN_m2: int, maxN_m3: int, maxN_m4: int,
                quick: bool = False) -> dict:
    """Run all four measurements on all functions across the requested N range."""

    # Truth tables sized for the largest N needed
    biggest_N = max(maxN_m1, maxN_m2, maxN_m3, maxN_m4)
    print(f"Computing truth tables on {{0,...,2^{biggest_N}-1}}...")
    t0 = time.time()
    tables = function_truth_tables(biggest_N)
    print(f"  done in {time.time() - t0:.2f}s")
    print()

    func_names = [
        "f_pi",
        "f_sqfree",
        "f_mu_pos",
        "f_lam_pos",
        "f_sqfree3",
        "f_prf_matched",
    ]

    densities = {name: float(tables[name].mean()) for name in func_names}
    print("Densities (at largest N):")
    for name in func_names:
        print(f"  {name:20s}  density = {densities[name]:.4f}")
    print()

    out = {"densities_largest_N": densities, "biggest_N": biggest_N}

    # ---- M1 communication-matrix rank ----
    print("=" * 72)
    print("M1: Communication-matrix rank (balanced split)")
    print("=" * 72)
    m1 = {}
    Ns_m1 = list(range(6, maxN_m1 + 1, 2))
    print(f"Ns: {Ns_m1}")
    for name in func_names:
        m1[name] = {}
        for N in Ns_m1:
            f = tables[name][:2 ** N]
            rank, rows, cols = comm_matrix_rank(f, N)
            half = N // 2
            anomaly = rank - 2 ** half  # for chi_P this is +2 - 2 = 0 actually...
            # NB: 2^{N/2-1} + 2 vs 2^{N/2} : difference = 2 - 2^{N/2-1}
            # The "+2" anomaly in E1.4 is rank = 2^{N/2-1} + 2 when N is EVEN
            # and we use a balanced split.  Easier to just record rank vs 2^k.
            m1[name][N] = {
                "rank": rank,
                "rows": rows,
                "cols": cols,
                "vs_full": rank - min(rows, cols),
                "vs_half": rank - 2 ** (half - 1),
            }
            print(f"  {name:20s}  N={N:>2}  rank={rank:>6}  full={min(rows,cols):>6}  "
                  f"rank - 2^(N/2-1) = {rank - 2**(half-1):>6}")
        print()
    out["M1_comm_rank"] = m1

    # ---- M2 Berlekamp-Massey ----
    print("=" * 72)
    print("M2: GF(2) Berlekamp-Massey linear complexity (length 2^N)")
    print("=" * 72)
    m2 = {}
    Ns_m2 = list(range(10, maxN_m2 + 1, 2))
    if quick:
        Ns_m2 = [n for n in Ns_m2 if n <= 12]
    print(f"Ns: {Ns_m2}")
    for name in func_names:
        m2[name] = {}
        for N in Ns_m2:
            seq = tables[name][:2 ** N].astype(np.int8)
            t0 = time.time()
            L = berlekamp_massey_gf2(seq)
            dt = time.time() - t0
            ratio = L / len(seq)
            m2[name][N] = {"L": L, "M": len(seq), "ratio": ratio, "time": dt}
            print(f"  {name:20s}  N={N:>2}  L={L:>7}  M=2^N={len(seq):>7}  "
                  f"L/M={ratio:.4f}  ({dt:.1f}s)")
        print()
    out["M2_lfsr"] = m2

    # ---- M3 approximate degree ----
    print("=" * 72)
    print("M3: Approximate degree at eps = 0.49")
    print("=" * 72)
    m3 = {}
    Ns_m3 = list(range(4, maxN_m3 + 1, 2))
    print(f"Ns: {Ns_m3}")
    for name in func_names:
        m3[name] = {}
        for N in Ns_m3:
            f = tables[name][:2 ** N]
            t0 = time.time()
            adeg, eps_by_d = approx_degree_eps49(f, N, max_mono=4500)
            dt = time.time() - t0
            m3[name][N] = {
                "adeg": adeg,
                "eps_by_degree": {str(k): v for k, v in eps_by_d.items()},
                "time": dt,
            }
            adeg_str = str(adeg) if adeg is not None else "?"
            print(f"  {name:20s}  N={N:>2}  adeg(0.49)={adeg_str}  ({dt:.1f}s)")
        print()
    out["M3_approx_degree"] = m3

    # ---- M4 PTF degree ----
    print("=" * 72)
    print("M4: PTF degree (sign-representation)")
    print("=" * 72)
    m4 = {}
    Ns_m4 = list(range(4, maxN_m4 + 1, 2))
    print(f"Ns: {Ns_m4}")
    for name in func_names:
        m4[name] = {}
        for N in Ns_m4:
            f = tables[name][:2 ** N]
            t0 = time.time()
            d = ptf_degree(f, N, max_mono=4500)
            dt = time.time() - t0
            m4[name][N] = {"ptf": d, "time": dt}
            print(f"  {name:20s}  N={N:>2}  ptf={d}  ({dt:.1f}s)")
        print()
    out["M4_ptf"] = m4

    return out


def evaluate_predictions(out: dict) -> dict:
    """Apply PR1-PR4 falsification criteria."""
    func_names = list(out["M1_comm_rank"].keys())
    verdicts = {}

    # PR1: comm rank in {2^{k-1}, 2^{k-1}+1, +2, 2^k}
    pr1_pass = True
    pr1_details = {}
    for name in func_names:
        for N, rec in out["M1_comm_rank"][name].items():
            half = N // 2
            band_lo = 2 ** (half - 1)
            band_hi = 2 ** half
            r = rec["rank"]
            ok = band_lo <= r <= band_hi
            if not ok:
                pr1_pass = False
            pr1_details.setdefault(name, {})[N] = {
                "rank": r,
                "in_band": ok,
                "band": [band_lo, band_hi],
            }
    verdicts["PR1"] = {"pass": pr1_pass, "details": pr1_details}

    # PR2: BM ratio in [0.4, 0.6]
    pr2_pass = True
    pr2_details = {}
    for name in func_names:
        if name not in out["M2_lfsr"]:
            continue
        for N, rec in out["M2_lfsr"][name].items():
            ratio = rec["ratio"]
            ok = 0.4 <= ratio <= 0.6
            if not ok:
                pr2_pass = False
            pr2_details.setdefault(name, {})[N] = {
                "ratio": ratio,
                "in_band": ok,
            }
    verdicts["PR2"] = {"pass": pr2_pass, "details": pr2_details}

    # PR3: adeg(eps=0.49) == ceil(N/2) for all functions and N
    pr3_pass = True
    pr3_details = {}
    for name in func_names:
        for N, rec in out["M3_approx_degree"][name].items():
            target = math.ceil(N / 2)
            adeg = rec["adeg"]
            ok = adeg == target
            if not ok:
                pr3_pass = False
            pr3_details.setdefault(name, {})[N] = {
                "adeg": adeg,
                "target": target,
                "match": ok,
            }
    verdicts["PR3"] = {"pass": pr3_pass, "details": pr3_details}

    # PR4: ptf degree == ceil(N/2)
    pr4_pass = True
    pr4_details = {}
    for name in func_names:
        if name not in out["M4_ptf"]:
            continue
        for N, rec in out["M4_ptf"][name].items():
            target = math.ceil(N / 2)
            ptf = rec["ptf"]
            ok = ptf == target
            if not ok:
                pr4_pass = False
            pr4_details.setdefault(name, {})[N] = {
                "ptf": ptf,
                "target": target,
                "match": ok,
            }
    verdicts["PR4"] = {"pass": pr4_pass, "details": pr4_details}

    return verdicts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--maxN-m1", type=int, default=14)
    parser.add_argument("--maxN-m2", type=int, default=14)
    parser.add_argument("--maxN-m3", type=int, default=10)
    parser.add_argument("--maxN-m4", type=int, default=8)
    parser.add_argument("--quick", action="store_true",
                        help="Smaller N for fast iteration.")
    parser.add_argument("--out", default="n_over_2_universality_class_data.json")
    args = parser.parse_args()

    if args.quick:
        args.maxN_m1 = min(args.maxN_m1, 12)
        args.maxN_m2 = min(args.maxN_m2, 12)
        args.maxN_m3 = min(args.maxN_m3, 8)
        args.maxN_m4 = min(args.maxN_m4, 6)

    out = run_battery(args.maxN_m1, args.maxN_m2, args.maxN_m3, args.maxN_m4,
                      quick=args.quick)
    out["verdicts"] = evaluate_predictions(out)

    print("=" * 72)
    print("FINAL VERDICTS")
    print("=" * 72)
    for k, v in out["verdicts"].items():
        print(f"  {k}: {'PASS' if v['pass'] else 'FAIL'}")

    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=2, default=str)
    print(f"\nData written to {args.out}")


if __name__ == "__main__":
    main()
