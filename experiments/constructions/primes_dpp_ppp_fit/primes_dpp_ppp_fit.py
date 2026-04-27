"""
ATTACK_VECTORS §D7 — DPP / Permanantal-PP fit to the integer prime sequence.

Question: is chi_P consistent with being a translation-invariant
determinantal point process (DPP) or its permanantal counterpart (PPP)?

Mathematical setup
------------------
A translation-invariant point process on Z with intensity rho and a
real-symmetric kernel K(t) = K(-t) has 2-point inclusion probability

    R_2(t) = P(n in X, n+t in X) = K(0,0) K(t,t) - K(0,t)^2     (DPP, det)
           = K(0,0) K(t,t) + K(0,t)^2                           (PPP, perm)

Setting K(0) = rho (since R_1(n) = K(n, n) = rho), this gives:

    R_2(t) = rho^2 - K(t)^2     (DPP -- repulsive: R_2 <= rho^2)
    R_2(t) = rho^2 + K(t)^2     (PPP -- attractive: R_2 >= rho^2)

For a DPP, K(t)^2 = rho^2 - R_2(t) >= 0 is REQUIRED at every t.
For a PPP, K(t)^2 = R_2(t) - rho^2 >= 0 is REQUIRED at every t.

Hardy-Littlewood prime k-tuple conjecture predicts
    R_2(t) = rho^2 * S(0, t) / [scaling],   (asymptotically rho^2 S(t))
where S(0, t) is the singular series.

  - For t even: S(0, t) > 1 (twin / cousin / sexy primes have positive
    density). PPP is locally feasible, DPP is INFEASIBLE.
  - For t odd:  S(0, t) = 0 (only one even prime). DPP locally feasible,
    PPP INFEASIBLE.

Therefore primes cannot be a translation-invariant real DPP nor PPP
because R_2(t) - rho^2 changes sign across the odd/even offset partition.

Three-point test (deeper):
For PPP, restricted to admissible all-even triples (0, t_1, t_2):
    R_3^PPP(t_1, t_2) = perm[K(0), K(t_1), K(t_2); ...]
                      = rho^3 + rho [K(t_1)^2 + K(t_2)^2 + K(t_2-t_1)^2]
                              + 2 K(t_1) K(t_2) K(t_2-t_1)

with K(t) = rho * sqrt(S(0, t) - 1) extracted from the pair correlations.

Compare R_3^PPP to R_3^HL = rho^3 * S(0, t_1, t_2). These match iff a
specific factorisation identity holds among S(0, t_1, t_2) and the
pairwise singular series. If the empirical 3-point R_3^emp matches
R_3^HL but DEVIATES from R_3^PPP by > 5% at any single triple, primes
are quantitatively NOT a PPP.

Falsification statements (pre-stated)
-------------------------------------
F1 (DPP pair-level): for every admissible even t > 0, K^2_DPP(t)
   = rho^2 - R_2(t) is NEGATIVE -- DPP fails at the pair level.

F2 (PPP pair-level): for every odd t > 1, K^2_PPP(t) = R_2(t) - rho^2
   is NEGATIVE -- PPP fails at the pair level.

F3 (PPP 3-point identity): for some admissible triple (0, t_1, t_2)
   with t_1, t_2 both even, the relative gap
       |R_3^PPP - R_3^HL| / R_3^HL  > 0.10
   i.e., the permanent identity does NOT reproduce HL even when we
   plug in the HL-correct pair K(t).

If F1, F2, F3 all hold: B-grade structural negative -- "primes are not
a translation-invariant signed DPP/PPP, in any of the standard formal
senses, even at the 3-point level on the all-even sub-lattice".

A-grade success (NOT EXPECTED): R_3^PPP matches R_3^HL within sample
noise at every tested triple, and matches empirical R_3^emp -- primes
admit a non-trivial PPP structure with kernel K(t); kernel becomes a
new arithmetic invariant.

Author: session95.
Cross-domain ref: Hough-Krishnapur-Peres-Virag 2009 *Zeros of Gaussian
Analytic Functions and DPP* (AMS ULect 51); Soshnikov 2000
arXiv:math/0002099.
"""

from __future__ import annotations
import json
import os
import time
import argparse

import numpy as np


# ---- Sieve & chi_P --------------------------------------------------------

def build_chi_P(N: int) -> np.ndarray:
    """Return uint8 array chi[0..N-1] with chi[n] = 1 iff n is prime."""
    is_prime = np.ones(N, dtype=bool)
    if N > 0:
        is_prime[0] = False
    if N > 1:
        is_prime[1] = False
    pmax = int(np.sqrt(N)) + 1
    for p in range(2, pmax + 1):
        if is_prime[p]:
            is_prime[p * p::p] = False
    return is_prime.astype(np.uint8)


# ---- Pair / triple correlations ------------------------------------------

def pair_correlations(chi: np.ndarray, t_max: int) -> np.ndarray:
    """R_2[t] = (1/N) sum_n chi[n] chi[n+t] for t in [0, t_max].
    Includes t=0 (gives rho exactly).
    """
    N = len(chi)
    out = np.zeros(t_max + 1, dtype=np.float64)
    out[0] = float(chi.sum()) / N
    chi_int = chi.astype(np.int64)
    for t in range(1, t_max + 1):
        out[t] = float((chi_int[: N - t] * chi_int[t:]).sum()) / N
    return out


def triple_correlations(chi: np.ndarray, triples: list[tuple[int, int]]) -> dict:
    """R_3[(t1, t2)] = (1/N) sum_n chi[n] chi[n+t1] chi[n+t2]."""
    N = len(chi)
    chi_int = chi.astype(np.int64)
    res = {}
    for (t1, t2) in triples:
        if not (0 < t1 < t2):
            raise ValueError(f"require 0 < t1 < t2, got {(t1, t2)}")
        n_max = N - t2
        v = chi_int[:n_max] * chi_int[t1: t1 + n_max] * chi_int[t2: t2 + n_max]
        res[(t1, t2)] = float(v.sum()) / N
    return res


# ---- Hardy-Littlewood singular series -------------------------------------

def hl_singular_series(h_tuple: tuple[int, ...], primes: list[int]) -> float:
    """Singular series for the prime k-tuple h_tuple = (h_1, ..., h_k).

    S(h) = prod_p (1 - nu_p(h)/p) / (1 - 1/p)^k

    where nu_p(h) = # distinct residues of h_i mod p.
    Convention: h must be admissible (no p | n+h_i for all i means nu_p < p).
    Returns 0 if any prime saturates nu_p = p (inadmissible).
    """
    k = len(h_tuple)
    val = 1.0
    for p in primes:
        residues = set(x % p for x in h_tuple)
        nu = len(residues)
        if nu == p:
            return 0.0
        # factor = (p - nu) p^(k-1) / (p - 1)^k
        factor = (1.0 - nu / p) / (1.0 - 1.0 / p) ** k
        val *= factor
    return val


# ---- DPP / PPP prediction formulas ---------------------------------------

def perm_3x3_K(rho: float, K_t1: float, K_t2: float, K_dt: float) -> float:
    """Permanent of the 3x3 K-matrix for translation-invariant K(0)=rho."""
    return (
        rho**3
        + rho * (K_t1**2 + K_t2**2 + K_dt**2)
        + 2.0 * K_t1 * K_t2 * K_dt
    )


def det_3x3_K(rho: float, K_t1: float, K_t2: float, K_dt: float) -> float:
    """Determinant of the 3x3 K-matrix (real symmetric K)."""
    return (
        rho**3
        - rho * (K_t1**2 + K_t2**2 + K_dt**2)
        + 2.0 * K_t1 * K_t2 * K_dt
    )


# ---- Helpers --------------------------------------------------------------

def primes_up_to(P: int) -> list[int]:
    is_prime = [True] * (P + 1)
    is_prime[0] = is_prime[1] = False
    for p in range(2, int(P**0.5) + 1):
        if is_prime[p]:
            for k in range(p * p, P + 1, p):
                is_prime[k] = False
    return [p for p in range(P + 1) if is_prime[p]]


# ---- Main -----------------------------------------------------------------

def main(N: int, t_pair_max: int, hl_pmax: int, out_dir: str) -> dict:
    print(f"[1/4] Build chi_P up to N = {N:,}")
    t0 = time.time()
    chi = build_chi_P(N)
    pi_N = int(chi.sum())
    rho = pi_N / N
    print(f"      pi(N) = {pi_N}, rho = {rho:.6e}, log(N) = {np.log(N):.4f}")
    print(f"      built in {time.time() - t0:.2f}s")

    # Pre-compute triples to test:
    # all-even admissible 3-tuples (t1, t2) with t1 < t2 <= 30:
    even_triples = [
        (2, 6), (2, 8), (4, 10), (6, 12), (4, 14), (6, 16), (2, 12),
        (4, 6), (4, 16), (6, 8), (6, 18), (8, 14), (8, 18), (10, 16),
        (12, 18), (10, 22), (14, 20), (16, 22), (18, 24), (20, 26),
    ]
    # parity-mixed triples (one even one odd) -- expected R_3 ~ 0:
    mixed_triples = [(1, 3), (1, 5), (3, 5), (5, 7), (1, 7), (3, 7)]
    all_triples = sorted(set(even_triples + mixed_triples))
    needed = max(t2 for (_, t2) in all_triples)
    if t_pair_max < needed:
        print(f"      [INFO] extending t_pair_max from {t_pair_max} to {needed}")
        t_pair_max = needed

    print(f"[2/4] Compute R_2(t) for t in [0, {t_pair_max}]")
    t0 = time.time()
    R2 = pair_correlations(chi, t_pair_max)
    print(f"      done in {time.time() - t0:.2f}s")

    print(f"[3/4] Compute R_3 for {len(all_triples)} triples")
    t0 = time.time()
    R3 = triple_correlations(chi, all_triples)
    print(f"      done in {time.time() - t0:.2f}s")

    print(f"[4/4] HL singular series with primes up to {hl_pmax}")
    primes = primes_up_to(hl_pmax)
    print(f"      using {len(primes)} primes")

    # ----- Pair-level analysis: F1 (DPP) and F2 (PPP) ----------------------
    pair_rows = []
    f1_dpp_violated_at_even = []
    f2_ppp_violated_at_odd = []
    for t in range(1, t_pair_max + 1):
        S_t = hl_singular_series((0, t), primes)
        R2_HL = rho**2 * S_t
        K2_DPP = rho**2 - R2[t]
        K2_PPP = R2[t] - rho**2
        # ratio
        ratio = (R2[t] / rho**2) if rho > 0 else 0.0
        row = {
            "t": t,
            "R2_emp": R2[t],
            "S_HL": S_t,
            "R2_HL": R2_HL,
            "ratio_R2_over_rho2": ratio,
            "K2_DPP": K2_DPP,        # = rho^2 - R2 ; DPP needs >= 0
            "K2_PPP": K2_PPP,        # = R2 - rho^2 ; PPP needs >= 0
        }
        pair_rows.append(row)
        if t % 2 == 0 and S_t > 0 and K2_DPP < 0:
            f1_dpp_violated_at_even.append(t)
        if t % 2 == 1 and t > 1 and K2_PPP < 0:
            f2_ppp_violated_at_odd.append(t)

    # ----- Triple-level analysis: F3 (PPP permanent identity) -------------
    triple_rows = []
    max_ppp_hl_gap_pct = 0.0
    f3_ppp_violated_triples = []
    for (t1, t2) in all_triples:
        dt = t2 - t1
        # HL prediction
        S3 = hl_singular_series((0, t1, t2), primes)
        R3_HL = rho**3 * S3
        # PPP prediction
        K2_t1 = R2[t1] - rho**2
        K2_t2 = R2[t2] - rho**2
        K2_dt = R2[dt] - rho**2
        # If any K^2 < 0, PPP is locally infeasible
        ppp_feasible = (K2_t1 >= 0) and (K2_t2 >= 0) and (K2_dt >= 0)
        if ppp_feasible:
            K_t1 = np.sqrt(K2_t1)
            K_t2 = np.sqrt(K2_t2)
            K_dt = np.sqrt(K2_dt)
            R3_PPP = perm_3x3_K(rho, K_t1, K_t2, K_dt)
        else:
            R3_PPP = None
        # DPP prediction: only feasible if all K^2_DPP >= 0
        K2D_t1 = rho**2 - R2[t1]
        K2D_t2 = rho**2 - R2[t2]
        K2D_dt = rho**2 - R2[dt]
        dpp_feasible = (K2D_t1 >= 0) and (K2D_t2 >= 0) and (K2D_dt >= 0)
        if dpp_feasible:
            R3_DPP = det_3x3_K(rho, np.sqrt(K2D_t1), np.sqrt(K2D_t2), np.sqrt(K2D_dt))
        else:
            R3_DPP = None
        # Empirical
        R3_emp = R3[(t1, t2)]
        # Sigma estimate (Poisson-like): sigma ~ sqrt(R3_emp / N)
        sigma = (R3_emp / N) ** 0.5

        def pct(num, den):
            if den is None or den == 0 or num is None:
                return None
            return 100.0 * (num - den) / den

        row = {
            "t1": t1,
            "t2": t2,
            "tuple": (0, t1, t2),
            "S3_HL": S3,
            "R3_HL": R3_HL,
            "R3_emp": R3_emp,
            "R3_PPP": R3_PPP,
            "R3_DPP": R3_DPP,
            "ppp_feasible": ppp_feasible,
            "dpp_feasible": dpp_feasible,
            "emp_minus_HL_pct": pct(R3_emp, R3_HL),
            "PPP_minus_HL_pct": pct(R3_PPP, R3_HL),
            "DPP_minus_HL_pct": pct(R3_DPP, R3_HL),
            "emp_minus_PPP_pct": pct(R3_emp, R3_PPP),
            "sigma_emp": sigma,
            "rel_sigma_pct": (100.0 * sigma / R3_emp) if R3_emp > 0 else None,
        }
        triple_rows.append(row)

        # Signed-K analysis: solve for the sign cross-term sigma s.t.
        #   R_3^HL = rho^3 + rho [K(t1)^2 + K(t2)^2 + K(dt)^2] + 2 sigma |K(t1)||K(t2)||K(dt)|
        # Required sigma for HL match:
        if ppp_feasible and R3_HL > 0 and (K2_t1 * K2_t2 * K2_dt) > 0:
            num = R3_HL - rho**3 - rho * (K2_t1 + K2_t2 + K2_dt)
            den = 2.0 * np.sqrt(K2_t1 * K2_t2 * K2_dt)
            sigma_req = num / den
        else:
            sigma_req = None
        row["sigma_HL_required"] = sigma_req
        # For real K with magnitudes |K(t)| = rho sqrt(S(t)-1) (fixed by pair):
        #   real-signed K matches HL <=> sigma_req in {-1, +1} (with tolerance).
        # For complex Hermitian K of same magnitudes:
        #   complex K LOCALLY matches HL <=> |sigma_req| <= 1 (a phase exists).
        TOL = 0.05
        row["sigma_real_signed_match"] = (sigma_req is not None) and (
            abs(abs(sigma_req) - 1.0) <= TOL
        )
        row["sigma_complex_feasible"] = (sigma_req is not None) and (
            abs(sigma_req) <= 1.0 + TOL
        )

        if ppp_feasible and R3_HL > 0:
            gap_pct = abs(row["PPP_minus_HL_pct"])
            if gap_pct > max_ppp_hl_gap_pct:
                max_ppp_hl_gap_pct = gap_pct
            if gap_pct > 10.0:
                f3_ppp_violated_triples.append((t1, t2, gap_pct))

    # ----- Summary --------------------------------------------------------
    summary = {
        "config": {
            "N": N,
            "t_pair_max": t_pair_max,
            "hl_pmax": hl_pmax,
        },
        "pi_N": pi_N,
        "rho": rho,
        "log_N": float(np.log(N)),
        "rho_log_N": rho * np.log(N),
        "twin_C2_HL_S_at_t2": pair_rows[1]["S_HL"],
        "F1_dpp_violated_count_even_t": len(f1_dpp_violated_at_even),
        "F1_dpp_violated_offsets_first10": f1_dpp_violated_at_even[:10],
        "F2_ppp_violated_count_odd_t": len(f2_ppp_violated_at_odd),
        "F2_ppp_violated_offsets_first10": f2_ppp_violated_at_odd[:10],
        "F3_max_PPP_minus_HL_pct_abs": max_ppp_hl_gap_pct,
        "F3_violated_triples_at_10pct": f3_ppp_violated_triples,
        "pair_rows_first40": pair_rows[:40],
        "triple_rows": triple_rows,
    }

    out_path = os.path.join(out_dir, "main_run.json")
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\nSaved summary to {out_path}")

    # ----- Pretty-print key tables ---------------------------------------
    print("\n" + "=" * 78)
    print("PAIR LEVEL: small offsets (showing R_2(t)/rho^2 and HL singular series)")
    print("=" * 78)
    print(f"{'t':>3}  {'parity':>6}  {'R2_emp':>12}  {'rho^2':>12}  "
          f"{'R2/rho^2':>9}  {'S_HL':>9}  {'K2_DPP':>13}  {'K2_PPP':>13}")
    for row in pair_rows[:30]:
        parity = "odd" if row["t"] % 2 else "even"
        print(f"{row['t']:>3}  {parity:>6}  {row['R2_emp']:.6e}  {rho**2:.6e}  "
              f"{row['ratio_R2_over_rho2']:>9.5f}  {row['S_HL']:>9.5f}  "
              f"{row['K2_DPP']:>13.4e}  {row['K2_PPP']:>13.4e}")

    print("\n" + "=" * 78)
    print("TRIPLE LEVEL: PPP permanent prediction vs HL prediction vs empirical")
    print("=" * 78)
    print(f"{'tuple':>14}  {'R3_emp':>10}  {'R3_HL':>10}  {'R3_PPP':>10}  "
          f"{'PPP-HL%':>8}  {'emp-HL%':>8}  {'emp-PPP%':>9}  {'feas':>4}")
    for row in triple_rows:
        ppp = f"{row['R3_PPP']:.3e}" if row['R3_PPP'] is not None else "    --   "
        ph = f"{row['PPP_minus_HL_pct']:>7.2f}" if row['PPP_minus_HL_pct'] is not None else "  --"
        eh = f"{row['emp_minus_HL_pct']:>7.2f}" if row['emp_minus_HL_pct'] is not None else "  --"
        ep = f"{row['emp_minus_PPP_pct']:>8.2f}" if row['emp_minus_PPP_pct'] is not None else "  --"
        feas = "yes" if row['ppp_feasible'] else "no"
        print(f"{str(row['tuple']):>14}  {row['R3_emp']:.3e}  {row['R3_HL']:.3e}  {ppp:>10}  "
              f"{ph:>8}  {eh:>8}  {ep:>9}  {feas:>4}")

    print("\n" + "=" * 78)
    print("FALSIFICATION OUTCOMES")
    print("=" * 78)
    print(f"  F1 (DPP fails at pair level for even t): "
          f"{len(f1_dpp_violated_at_even)} violations of {sum(1 for r in pair_rows if r['t']%2==0 and r['S_HL']>0)} admissible even t")
    print(f"     First 10 offending offsets: {f1_dpp_violated_at_even[:10]}")
    print(f"  F2 (PPP fails at pair level for odd t > 1): "
          f"{len(f2_ppp_violated_at_odd)} violations")
    print(f"     First 10 offending offsets: {f2_ppp_violated_at_odd[:10]}")
    print(f"  F3 (PPP permanent vs HL gap > 10% on even triples): "
          f"max |gap| = {max_ppp_hl_gap_pct:.2f}%")
    print(f"     Violating triples: {f3_ppp_violated_triples}")

    # ----- Signed-K test (could a signed real / complex Hermitian K fit?) -
    print("\n" + "=" * 78)
    print("SIGNED-K TEST: required cross-term ratio sigma to match HL at each triple")
    print("  sigma_req := [R_3^HL - rho^3 - rho * sum K(t)^2] / [2 |K(t1)||K(t2)||K(dt)|]")
    print("  Real-signed K (signs ±1) matches HL only if sigma_req ∈ {-1, +1}")
    print("  Complex Hermitian K can locally match if |sigma_req| ≤ 1 (phase choice)")
    print("=" * 78)
    print(f"{'tuple':>14}  {'sigma_req':>10}  {'real-signed?':>13}  "
          f"{'complex-feas?':>14}")
    n_real_sgn = 0
    n_real_fail = 0
    n_complex_feas = 0
    real_signed_failures = []
    for row in triple_rows:
        sr = row["sigma_HL_required"]
        if sr is None:
            continue
        rl_real = "YES" if row["sigma_real_signed_match"] else "no"
        rl_cplx = "YES" if row["sigma_complex_feasible"] else "no"
        print(f"{str(row['tuple']):>14}  {sr:>10.6f}  {rl_real:>13}  {rl_cplx:>14}")
        if row["sigma_real_signed_match"]:
            n_real_sgn += 1
        else:
            n_real_fail += 1
            real_signed_failures.append((row["tuple"], sr))
        if row["sigma_complex_feasible"]:
            n_complex_feas += 1
    print(f"  Real-signed K matches HL at {n_real_sgn} triples, FAILS at {n_real_fail}.")
    print(f"  Complex Hermitian K: locally feasible at {n_complex_feas} of "
          f"{n_real_sgn + n_real_fail} feasible triples.")
    summary["real_signed_K_failures"] = real_signed_failures
    summary["real_signed_match_count"] = n_real_sgn
    summary["real_signed_fail_count"] = n_real_fail
    summary["complex_locally_feasible_count"] = n_complex_feas

    # ----- Complex-K phase consistency probe ------------------------------
    # For each triple, sigma = cos(phi(t1) + phi(dt) - phi(t2)) where dt = t2 - t1.
    # We test: pick phi(t) ∈ [0, 2π) for offsets {2, 4, 6, ..., dt_max}.
    # Set phi(2) = 0 (gauge fix). For each triple write equation
    #   cos(phi(t1) + phi(t2-t1) - phi(t2)) = sigma_req
    # Try simple greedy assignment; report max residual.
    print("\n" + "=" * 78)
    print("COMPLEX-K PHASE CONSISTENCY PROBE")
    print("=" * 78)

    # Collect distinct offsets:
    offsets = sorted({t for row in triple_rows for t in row["tuple"][1:]} |
                     {row["tuple"][2] - row["tuple"][1] for row in triple_rows})
    offsets = [t for t in offsets if t > 0]
    # Use only offsets where K^2 > 0 (PPP-feasible)
    feasible_offsets = [t for t in offsets if t < len(R2) and (R2[t] - rho**2) > 0]
    print(f"  Distinct PPP-feasible offsets used: {feasible_offsets}")

    # Try least-squares fit of phi via scipy.
    try:
        from scipy.optimize import least_squares

        # Map offset -> index in phi vector
        idx = {t: i for i, t in enumerate(feasible_offsets)}

        # Build residuals
        triple_eqs = []
        for row in triple_rows:
            sr = row["sigma_HL_required"]
            if sr is None or not row["sigma_complex_feasible"]:
                continue
            t1, t2 = row["t1"], row["t2"]
            dt = t2 - t1
            if t1 in idx and t2 in idx and dt in idx:
                triple_eqs.append((idx[t1], idx[dt], idx[t2], sr))

        def residuals(phi_vec):
            res = []
            for (i1, idt, i2, sr) in triple_eqs:
                res.append(np.cos(phi_vec[i1] + phi_vec[idt] - phi_vec[i2]) - sr)
            return np.array(res)

        # Try many random starts; multistart over both LM and trf methods
        best_max_res = np.inf
        best_phi = None
        rng = np.random.default_rng(42)
        for trial in range(200):
            phi0 = rng.uniform(-np.pi, np.pi, size=len(feasible_offsets))
            phi0[0] = 0.0  # gauge fix the smallest offset
            method = "lm" if trial % 2 == 0 else "trf"
            try:
                sol = least_squares(residuals, phi0, max_nfev=4000, method=method)
                rs = residuals(sol.x)
                mr = np.max(np.abs(rs))
                if mr < best_max_res:
                    best_max_res = mr
                    best_phi = sol.x
            except Exception:
                pass

        print(f"  Triples fed to LS fit: {len(triple_eqs)}")
        print(f"  Best max residual |cos(...) - sigma_req|: {best_max_res:.6f}")
        if best_phi is not None and best_max_res < 0.01:
            print("  -> Complex Hermitian K phase assignment EXISTS within tolerance.")
            print(f"     Phi(t) for t in {feasible_offsets}: "
                  f"{[round(float(p) % (2 * np.pi), 4) for p in best_phi]}")
        else:
            print("  -> No complex K phase assignment passes within 0.01 tolerance.")
            if best_phi is not None:
                rs = residuals(best_phi)
                worst = max(zip(rs, triple_eqs), key=lambda x: abs(x[0]))
                print(f"     Worst residual: {worst[0]:.4f} at indices {worst[1][:3]}")
        summary["complex_K_max_residual"] = float(best_max_res)
    except ImportError:
        print("  (scipy not available; skipping LS fit)")

    # Re-save with sigma fields
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)

    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=10_000_000)
    parser.add_argument("--t_pair_max", type=int, default=30)
    parser.add_argument("--hl_pmax", type=int, default=5000)
    parser.add_argument("--out_dir", type=str,
                        default=os.path.dirname(os.path.abspath(__file__)) or ".")
    args = parser.parse_args()
    main(args.N, args.t_pair_max, args.hl_pmax, args.out_dir)
