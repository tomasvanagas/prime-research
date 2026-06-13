"""
slot4_alpha_derivation.py — Thread 9 (P4) commit slot 4

Goal: derive (or pin down) the closed-form leading-log coefficient alpha in

    Delta(w) = T(w) - 2 C_2^2 w^2  ~  -alpha * w * log(w) + beta * w + o(w)

where T(w) = sum_{m in A, m <= w-1} (w-m) S_4(0,2,m,m+2) is the HL
4-tuple singular-series triangular sum (slot 3 result F9 / F10).

Slot 3's empirical fit: alpha ~ 5.72, beta ~ 10.50 (intercept matches
24 C_2^2 = 10.4596 to 0.35%). Slot 3 proposed structural candidate
alpha = 24 C_2^2 / log(6) = 5.838 (within 2% of empirical fit).

Slot 4 strategy:

  1. Compute S_4(0, 2, m, m+2) for all admissible m up to M_MAX ~ 200,000
     (~33,000 admissible terms) using one-shot tail-cached singular series
     evaluation. Cache the per-m values for ALL downstream sums.

  2. Compute H(N) = sum_{m in A, m <= N} S_4(m) for N = 60, 120, ..., M_MAX.
     Tabulate H(N) - 4 C_2^2 * N (the "leading-cancellation residual").

  3. Compute T(w) = sum_{m=1}^{w-1} H(m) for w = 100, 200, ..., M_MAX.
     Tabulate Delta(w) = T(w) - 2 C_2^2 w^2.

  4. Extract alpha and beta via the integral relation:
        Delta(w) = ∫_0^w (H(u) - 4 C_2^2 u) du  +  edge.
     If H(u) - 4 C_2^2 u ~ a log(u) + b + o(1), then
        Delta(w) ~ a (w log w - w) + b w
                 = a * w log w + (b - a) * w.
     So alpha = -a, and beta = b - a.

  5. Test structural candidates for alpha:
        (i)   24 C_2^2 / log(6)         = 5.8378
        (ii)  4 C_2^2 * log(constant)   = various
        (iii) related Gallagher / Goldston-Yildirim closed forms.
     Compare to high-precision empirical alpha at largest w.

  6. Also extract (a, b) directly from H residuals (linear fit in log u),
     which is cleaner than fitting Delta (no Cesaro accumulation noise).

Output:
  - slot4_h_residual.csv        H(N), 4C_2^2*N, residual, log(N)
  - slot4_t_delta.csv           T(w), 2C_2^2*w^2, Delta, log(w)
  - slot4_alpha_fits.csv        alpha(M), beta(M) from rolling fits
  - slot4_run.log               full run log
  - slot4_alpha_derivation_results.md      results writeup with theorem statement
"""

from __future__ import annotations

import csv
import math
import os
import time
from pathlib import Path

# Sieve primes up to PRIME_BOUND for singular series evaluation; primes
# above this are folded into the asymptotic tail factor.
PRIME_BOUND = 200_000

# Total admissible m's to compute (covers w up to ~6 * M_TERMS in practice).
# 200,000 admissible m's = m up to 1,200,000.
M_TERMS = 200_000


# ---------------------------------------------------------------------------
# Sieve and constants
# ---------------------------------------------------------------------------

def sieve_primes(N: int) -> list[int]:
    is_p = bytearray(b"\x01") * (N + 1)
    is_p[0] = is_p[1] = 0
    bound = int(N ** 0.5) + 1
    for i in range(2, bound):
        if is_p[i]:
            is_p[i * i : N + 1 : i] = b"\x00" * (((N - i * i) // i) + 1)
    return [i for i in range(2, N + 1) if is_p[i]]


PRIMES = sieve_primes(PRIME_BOUND)

# Twin-prime constant
def twin_prime_constant() -> float:
    val = 1.0
    for p in PRIMES:
        if p < 3:
            continue
        val *= (1.0 - 1.0 / (p - 1) ** 2)
    return val


C_2 = twin_prime_constant()
C2_SQ = C_2 * C_2
TARGET_S4_AVG = 24 * C2_SQ           # asymptotic <S_4>_admissible
LEADING_H_PER_M = 4 * C2_SQ          # H(N) ~ 4 C_2^2 * N to leading order


# ---------------------------------------------------------------------------
# Fast singular-series for tuple {0, 2, m, m+2}, m admissible (m % 6 == 0)
# ---------------------------------------------------------------------------
#
# At each prime p:
#   p = 2: factor = 8 (m even, ν=1)
#   p = 3: factor = 27/16 (m ≡ 0 mod 3, ν=2)
#   p >= 5: factor depends on m mod p as
#       m mod p == 0:        f_p(0) = (1 - 2/p) / (1-1/p)^4
#       m mod p == 2 or p-2: f_p(2) = (1 - 3/p) / (1-1/p)^4
#       otherwise:           f_p(other) = (1 - 4/p) / (1-1/p)^4
#
# For large p, f_p(other) ~ 1 - 4/p^2 + O(1/p^3).  We cache f_p(other) for
# every prime up to PRIME_BOUND, and only update factors for primes that
# divide m, m-2, or m+2 (the "special" primes for that m).  The "other"
# product at every prime is a fixed constant for all m; it is precomputed
# once.
#
# Concretely we write
#     S_4(m) = OTHER_PRODUCT * prod_{p | m, p>=5} (f_p(0)/f_p(other))
#                            * prod_{p | m-2 or p | m+2, p>=5}
#                                            (f_p(2)/f_p(other))
# for m admissible with m >= 6.
#
# OTHER_PRODUCT = factor at p=2 (=8) * factor at p=3 (=27/16) *
#                 prod_{p>=5} f_p(other).

# Precompute per-prime factors and ratios.
F_OTHER: dict[int, float] = {}       # f_p(other), p >= 5
RATIO_F0_OVER_FOTHER: dict[int, float] = {}   # f_p(0)/f_p(other)
RATIO_F2_OVER_FOTHER: dict[int, float] = {}   # f_p(2)/f_p(other)

for p in PRIMES:
    if p < 5:
        continue
    factor_other = (1.0 - 4.0 / p) / ((1.0 - 1.0 / p) ** 4)
    factor_0     = (1.0 - 2.0 / p) / ((1.0 - 1.0 / p) ** 4)
    factor_2     = (1.0 - 3.0 / p) / ((1.0 - 1.0 / p) ** 4)
    F_OTHER[p] = factor_other
    RATIO_F0_OVER_FOTHER[p] = factor_0 / factor_other
    RATIO_F2_OVER_FOTHER[p] = factor_2 / factor_other

# Precompute log of OTHER_PRODUCT for numerical stability.
LOG_OTHER_PRODUCT = math.log(8.0) + math.log(27.0 / 16.0)
for p, f in F_OTHER.items():
    LOG_OTHER_PRODUCT += math.log(f)
OTHER_PRODUCT = math.exp(LOG_OTHER_PRODUCT)


def small_prime_factors_ge5(n: int) -> list[int]:
    """Return distinct prime divisors of n that are >= 5."""
    out = []
    if n <= 0:
        return out
    # strip 2 and 3
    while n % 2 == 0:
        n //= 2
    while n % 3 == 0:
        n //= 3
    if n == 1:
        return out
    # trial divide by primes >= 5
    i = 5
    while i * i <= n:
        if n % i == 0:
            out.append(i)
            while n % i == 0:
                n //= i
        i += 2
        if i % 3 == 0:
            i += 2
    if n > 1:
        out.append(n)
    return out


def get_ratio_f0_other(p: int) -> float:
    """f_p(0)/f_p(other) for prime p >= 5. Cached."""
    if p in RATIO_F0_OVER_FOTHER:
        return RATIO_F0_OVER_FOTHER[p]
    fo = (1.0 - 4.0 / p) / ((1.0 - 1.0 / p) ** 4)
    f0 = (1.0 - 2.0 / p) / ((1.0 - 1.0 / p) ** 4)
    r = f0 / fo
    RATIO_F0_OVER_FOTHER[p] = r
    return r


def get_ratio_f2_other(p: int) -> float:
    """f_p(2)/f_p(other) for prime p >= 5. Cached."""
    if p in RATIO_F2_OVER_FOTHER:
        return RATIO_F2_OVER_FOTHER[p]
    fo = (1.0 - 4.0 / p) / ((1.0 - 1.0 / p) ** 4)
    f2 = (1.0 - 3.0 / p) / ((1.0 - 1.0 / p) ** 4)
    r = f2 / fo
    RATIO_F2_OVER_FOTHER[p] = r
    return r


def S4_admissible(m: int) -> float:
    """Compute S_4({0, 2, m, m+2}) for admissible m (m % 6 == 0, m >= 6).

    Uses the OTHER_PRODUCT and per-special-prime correction factors.
    """
    val = OTHER_PRODUCT

    # Primes p >= 5 dividing m: contribute f_p(0) / f_p(other)
    for p in small_prime_factors_ge5(m):
        val *= get_ratio_f0_other(p)

    # Primes p >= 5 dividing m-2 OR m+2: contribute f_p(2) / f_p(other)
    # Note: m-2 and m+2 share no prime factors >= 5 (they differ by 4, and
    # any common prime factor would divide 4, i.e., only p=2). So we just
    # union the divisor sets.
    p_minus = small_prime_factors_ge5(m - 2)
    p_plus = small_prime_factors_ge5(m + 2)
    for p in p_minus:
        val *= get_ratio_f2_other(p)
    for p in p_plus:
        val *= get_ratio_f2_other(p)

    return val


# ---------------------------------------------------------------------------
# Validate against slot 3's slow computation on the first 30 admissible m's
# ---------------------------------------------------------------------------

def slow_singular_series(m: int, P_max: int = PRIME_BOUND) -> float:
    """Reference (slow) singular series evaluator. Includes p's own factor
    when transitioning to tail mode (slot 3 had an off-by-one bug that
    skipped the boundary prime; the form below is correct)."""
    H = (0, 2, m, m + 2)
    val = 1.0
    k = len(H)
    diam = max(H) - min(H)
    for p in PRIMES:
        if p > P_max:
            break
        residues = set(h % p for h in H)
        nu = len(residues)
        if nu == p:
            return 0.0
        # Apply factor at p
        val *= (1.0 - nu / p) / ((1.0 - 1.0 / p) ** k)
        # If p > diam and nu == k, all subsequent primes have nu == k.
        # Apply tail product over q > p directly.
        if p > diam and nu == k:
            for q in PRIMES:
                if q <= p:
                    continue
                if q > P_max:
                    break
                val *= (1.0 - k / q) / ((1.0 - 1.0 / q) ** k)
            return val
    return val


def validate_against_slot3():
    """Compare fast S4 to slow S4 for first 30 admissible m's."""
    max_rel_err = 0.0
    rows = []
    for k in range(1, 31):
        m = 6 * k
        s_fast = S4_admissible(m)
        s_slow = slow_singular_series(m)
        rel = abs(s_fast - s_slow) / max(s_slow, 1e-12)
        max_rel_err = max(max_rel_err, rel)
        rows.append((m, s_fast, s_slow, rel))
    return max_rel_err, rows


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def main():
    out_dir = Path(__file__).parent
    log_lines: list[str] = []

    def log(msg: str = ""):
        print(msg)
        log_lines.append(str(msg))

    log("=" * 78)
    log("Thread 9 / P4 slot 4: Closed-form alpha derivation for Delta(w)")
    log("=" * 78)
    log(f"  C_2          = {C_2:.12f}")
    log(f"  C_2^2        = {C2_SQ:.12f}")
    log(f"  4  C_2^2     = {LEADING_H_PER_M:.12f}")
    log(f"  24 C_2^2     = {TARGET_S4_AVG:.12f}")
    log(f"  log(6)       = {math.log(6):.12f}")
    log(f"  24 C_2^2/log(6) = {TARGET_S4_AVG / math.log(6):.12f}")
    log(f"  PRIME_BOUND  = {PRIME_BOUND}")
    log(f"  M_TERMS      = {M_TERMS} (m up to {6 * M_TERMS})")
    log(f"  OTHER_PRODUCT = {OTHER_PRODUCT:.6f}")
    log("")

    # ------------------------------------------------------------------
    # Step 0: validate fast S_4 vs slot 3 slow evaluator
    # ------------------------------------------------------------------
    log("-" * 78)
    log("Step 0: validate fast S_4 against slot 3 slow evaluator (30 m's)")
    log("-" * 78)
    t0 = time.perf_counter()
    max_rel, rows = validate_against_slot3()
    log(f"  max relative error (fast vs slow) over m=6,...,180: {max_rel:.3e}")
    log(f"  time: {time.perf_counter()-t0:.2f}s")
    if max_rel > 1e-9:
        log(f"  *** WARNING: rel err exceeds 1e-9; review fast evaluator")
        for m, fa, sl, r in rows[:5]:
            log(f"     m={m}: fast={fa:.10f} slow={sl:.10f} rel={r:.3e}")
    else:
        log(f"  Fast and slow agree to {max_rel:.0e} (numerical roundoff).")
    log("")

    # ------------------------------------------------------------------
    # Step 1: precompute S_4 for all admissible m up to M_TERMS terms
    # ------------------------------------------------------------------
    log("-" * 78)
    log(f"Step 1: precompute S_4 for {M_TERMS} admissible m's (m up to {6*M_TERMS})")
    log("-" * 78)
    t0 = time.perf_counter()
    S4_values = [S4_admissible(6 * (k + 1)) for k in range(M_TERMS)]
    elapsed = time.perf_counter() - t0
    log(f"  computed {M_TERMS} S_4 values in {elapsed:.2f}s "
        f"({M_TERMS/elapsed:.0f} evals/sec)")
    log("")

    # Quick sanity: average should converge to 24 C_2^2
    cum_sum = 0.0
    cum_avg_log = []
    for k, s in enumerate(S4_values, start=1):
        cum_sum += s
        if k in [1, 10, 100, 1000, 10000, M_TERMS]:
            cum_avg_log.append((k, cum_sum / k, abs(cum_sum/k - TARGET_S4_AVG)/TARGET_S4_AVG))
    log(f"  Cumulative <S_4>_admissible -> 24 C_2^2 = {TARGET_S4_AVG:.6f}")
    for k, avg, rel in cum_avg_log:
        log(f"    k={k:>6d}: <S_4>={avg:.6f}  rel_err={rel:.3e}")
    log("")

    # ------------------------------------------------------------------
    # Step 2: compute H(N) - 4 C_2^2 * N for many N
    # ------------------------------------------------------------------
    log("-" * 78)
    log("Step 2: H(N) = sum_{admissible m <= N} S_4(m); tabulate residual")
    log("-" * 78)
    H_residual_rows = []
    cum = 0.0
    sample_K = sorted(set(
        list(range(10, 100, 10))
        + list(range(100, 1000, 50))
        + list(range(1000, 10000, 500))
        + list(range(10000, 50000, 2500))
        + list(range(50000, M_TERMS + 1, 10000))
    ))
    H_at_K: dict[int, float] = {}
    for k_idx, s in enumerate(S4_values, start=1):
        cum += s
        if k_idx in sample_K:
            H_at_K[k_idx] = cum
    # tabulate at sampled K's
    log(f"{'K':>6s} {'m_max':>7s} {'H(N)':>14s} {'4C2^2*N':>14s} {'residual':>11s} {'log(N)':>9s}")
    for K in sample_K:
        H_N = H_at_K[K]
        N = 6 * K
        leading = LEADING_H_PER_M * N
        residual = H_N - leading
        H_residual_rows.append({
            "K": K, "N": N, "H_N": H_N, "leading": leading,
            "residual": residual, "log_N": math.log(N),
        })
    # show every Nth row
    for row in H_residual_rows[::4]:
        log(f"{row['K']:>6d} {row['N']:>7d} {row['H_N']:>14.4f} {row['leading']:>14.4f} "
            f"{row['residual']:>11.4f} {row['log_N']:>9.4f}")
    log(f"  ... total {len(H_residual_rows)} sampled K's")
    log("")

    # Save to CSV
    h_path = out_dir / "slot4_h_residual.csv"
    with h_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["K", "N", "H_N", "leading", "residual", "log_N"])
        w.writeheader()
        for r in H_residual_rows:
            w.writerow({k: (f"{v:.6f}" if isinstance(v, float) else v) for k, v in r.items()})
    log(f"  H residual table -> {h_path}")
    log("")

    # ------------------------------------------------------------------
    # Step 3: linear fit on H residual vs log(N) to extract (a, b)
    #   H(N) - 4 C_2^2 * N = a log(N) + b + o(1)
    # ------------------------------------------------------------------
    log("-" * 78)
    log("Step 3: linear fit H_residual = a * log(N) + b, on tail (K >= 1000)")
    log("-" * 78)
    tail_rows = [r for r in H_residual_rows if r["K"] >= 1000]
    if len(tail_rows) >= 3:
        n = len(tail_rows)
        sx = sum(r["log_N"] for r in tail_rows)
        sy = sum(r["residual"] for r in tail_rows)
        sxx = sum(r["log_N"] ** 2 for r in tail_rows)
        sxy = sum(r["log_N"] * r["residual"] for r in tail_rows)
        a_h = (n * sxy - sx * sy) / (n * sxx - sx * sx)
        b_h = (sy - a_h * sx) / n
        log(f"  Fit: H_resid(N) ~ a * log(N) + b")
        log(f"     a = {a_h:.6f}    (=> alpha = -a = {-a_h:.6f})")
        log(f"     b = {b_h:.6f}")
        log(f"  fit data: K = {tail_rows[0]['K']}..{tail_rows[-1]['K']}, n={n}")
        # residuals
        max_dev = 0.0
        for r in tail_rows:
            pred = a_h * r["log_N"] + b_h
            dev = abs(r["residual"] - pred)
            max_dev = max(max_dev, dev)
        log(f"  max|residual - fit| over tail = {max_dev:.4f}")
    else:
        a_h = b_h = float("nan")
        log("  (insufficient tail data for fit)")
    log("")

    # ------------------------------------------------------------------
    # Step 4: compute T(w) for many w; fit alpha, beta
    # ------------------------------------------------------------------
    log("-" * 78)
    log("Step 4: T(w) = sum_{admissible m <= w-1} (w-m) S_4(m); Delta = T - 2 C_2^2 w^2")
    log("-" * 78)
    # Pick w's covering several decades
    # We'll only use w-aligned to admissible m boundaries: w = 6K + 5 so all
    # admissible m's up to 6K are included strictly.
    w_K_choices = sorted(set([
        50, 100, 200, 300, 500, 700, 1000, 1500, 2000, 3000, 5000,
        7500, 10000, 15000, 20000, 30000, 40000, 50000,
        70000, 100000, 150000, 200000,
    ]))

    # Also compute T(w) at the slot 3 w values, using exact triangular sum
    # (so we can directly compare slot 3's "F_pred" with the corrected
    # singular series).
    slot3_ws = [(190, 31), (259, 43), (339, 56), (2290, 381), (3117, 519), (4071, 678)]
    T_rows = []
    log(f"{'w':>8s} {'K':>7s} {'T(w)':>20s} {'2C_2^2*w^2':>20s} {'Delta':>14s} "
        f"{'log(w)':>9s} {'alpha_emp':>10s}")
    for K in w_K_choices:
        if K > M_TERMS:
            continue
        w = 6 * K + 5
        # T(w) = sum_{j=0}^{K-1} (w - 6(j+1)) * S_4(6(j+1))
        T = 0.0
        for j in range(K):
            T += (w - 6 * (j + 1)) * S4_values[j]
        leading = 2 * C2_SQ * w * w
        Delta = T - leading
        log_w = math.log(w)
        # alpha_emp from -Delta/(w log w) (still finite-w biased)
        alpha_emp_simple = -Delta / (w * log_w)
        T_rows.append({
            "w": w, "K": K, "T_w": T, "leading": leading, "Delta": Delta,
            "log_w": log_w, "alpha_emp": alpha_emp_simple,
        })
    for r in T_rows:
        log(f"{r['w']:>8d} {r['K']:>7d} {r['T_w']:>20.4f} {r['leading']:>20.4f} "
            f"{r['Delta']:>14.4f} {r['log_w']:>9.4f} {r['alpha_emp']:>10.4f}")
    log("")

    # Slot 3 comparison cells
    log("-" * 78)
    log("Step 4b: recompute T(w) at slot 3 w values with corrected S_4")
    log("-" * 78)
    log(f"  (slot 3 had off-by-one bug in singular-series tail; corrected here)")
    slot3_T_rows = []
    log(f"  {'w':>6s} {'K':>4s} {'T(w) [corrected]':>20s} {'2C_2^2*w^2':>20s} {'Delta_corr':>14s}")
    for w_s, K_s in slot3_ws:
        if K_s > M_TERMS:
            continue
        # T(w) = sum_{j=0}^{K_s-1}(w_s - 6(j+1)) S4
        Ts = 0.0
        for j in range(K_s):
            m = 6*(j+1)
            if m >= w_s:
                break
            Ts += (w_s - m) * S4_values[j]
        leading_s = 2 * C2_SQ * w_s * w_s
        Delta_s = Ts - leading_s
        slot3_T_rows.append({
            "w_slot3": w_s, "K": K_s, "T_corr": Ts, "leading": leading_s,
            "Delta_corr": Delta_s,
        })
        log(f"  {w_s:>6d} {K_s:>4d} {Ts:>20.4f} {leading_s:>20.4f} {Delta_s:>14.4f}")
    log("")

    s3cmp_path = out_dir / "slot4_slot3_comparison.csv"
    with s3cmp_path.open("w", newline="") as f:
        wcsv = csv.DictWriter(f, fieldnames=["w_slot3", "K", "T_corr", "leading", "Delta_corr"])
        wcsv.writeheader()
        for r in slot3_T_rows:
            wcsv.writerow({k: (f"{v:.6f}" if isinstance(v, float) else v) for k, v in r.items()})
    log(f"  slot3-comparison table -> {s3cmp_path}")
    log("")

    t_path = out_dir / "slot4_t_delta.csv"
    with t_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["w", "K", "T_w", "leading", "Delta", "log_w", "alpha_emp"])
        w.writeheader()
        for r in T_rows:
            w.writerow({k: (f"{v:.6f}" if isinstance(v, float) else v) for k, v in r.items()})
    log(f"  T/Delta table -> {t_path}")
    log("")

    # ------------------------------------------------------------------
    # Step 5: extract alpha, beta from Delta(w) = -alpha * w log w + beta * w
    # i.e., Delta/w = -alpha * log(w) + beta.
    # Linear fit on tail (large w).
    # ------------------------------------------------------------------
    log("-" * 78)
    log("Step 5: linear fit Delta(w)/w = -alpha * log(w) + beta, on tail (K >= 5000)")
    log("-" * 78)
    tail_T = [r for r in T_rows if r["K"] >= 5000]
    if len(tail_T) >= 3:
        n = len(tail_T)
        xs = [r["log_w"] for r in tail_T]
        ys = [r["Delta"] / r["w"] for r in tail_T]
        sx = sum(xs); sy = sum(ys)
        sxx = sum(x*x for x in xs)
        sxy = sum(x*y for x, y in zip(xs, ys))
        slope = (n * sxy - sx * sy) / (n * sxx - sx * sx)
        intercept = (sy - slope * sx) / n
        alpha_T = -slope
        beta_T = intercept
        log(f"  Fit: Delta/w ~ -alpha * log(w) + beta")
        log(f"     alpha = {alpha_T:.6f}   (slot 3 fit was 5.7165)")
        log(f"     beta  = {beta_T:.6f}    (slot 3 fit was 10.4958; 24 C_2^2 = {TARGET_S4_AVG:.4f})")
        log(f"  fit data: K = {tail_T[0]['K']}..{tail_T[-1]['K']}, n={n}")
        # residuals
        max_dev = 0.0
        for r in tail_T:
            pred = -alpha_T * r["log_w"] + beta_T
            dev = abs(r["Delta"]/r["w"] - pred)
            max_dev = max(max_dev, dev)
        log(f"  max|Delta/w - fit| over tail = {max_dev:.4e}")
    else:
        alpha_T = beta_T = float("nan")
        log("  (insufficient tail data for fit)")
    log("")

    # ------------------------------------------------------------------
    # Step 5b: log^2 fit. Delta(w) = w * (A * log^2 w + B * log w + C)
    # ------------------------------------------------------------------
    log("-" * 78)
    log("Step 5b: 3-parameter fit Delta(w) = w (A log^2 w + B log w + C)")
    log("-" * 78)
    A_fit = B_fit = C_fit = float("nan")
    for K_min in [100, 1000, 5000, 20000, 50000]:
        tail2 = [r for r in T_rows if r["K"] >= K_min]
        if len(tail2) < 4:
            continue
        n = len(tail2)
        # 3-var least squares via normal equations
        sx = [r["log_w"] for r in tail2]
        sy = [r["Delta"] / r["w"] for r in tail2]
        # Form: y = A x^2 + B x + C
        x1 = [s for s in sx]
        x2 = [s*s for s in sx]
        x0 = [1.0] * n
        # Normal eq matrix
        a11 = sum(xi*xi for xi in x2)
        a12 = sum(xi*xj for xi, xj in zip(x2, x1))
        a13 = sum(x2)
        a22 = sum(xi*xi for xi in x1)
        a23 = sum(x1)
        a33 = float(n)
        b1 = sum(xi*yi for xi, yi in zip(x2, sy))
        b2 = sum(xi*yi for xi, yi in zip(x1, sy))
        b3 = sum(sy)
        # Cramer's rule (3x3)
        def det3(M):
            return (M[0][0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1])
                    - M[0][1]*(M[1][0]*M[2][2] - M[1][2]*M[2][0])
                    + M[0][2]*(M[1][0]*M[2][1] - M[1][1]*M[2][0]))
        M = [[a11, a12, a13], [a12, a22, a23], [a13, a23, a33]]
        D = det3(M)
        Mx = [[b1, a12, a13], [b2, a22, a23], [b3, a23, a33]]
        My = [[a11, b1, a13], [a12, b2, a23], [a13, b3, a33]]
        Mz = [[a11, a12, b1], [a12, a22, b2], [a13, a23, b3]]
        A_v = det3(Mx) / D
        B_v = det3(My) / D
        C_v = det3(Mz) / D
        # Residuals
        max_res = 0.0
        for r, yi in zip(tail2, sy):
            pred = A_v * r["log_w"]**2 + B_v * r["log_w"] + C_v
            max_res = max(max_res, abs(yi - pred))
        log(f"  K_min={K_min:>6d}, n={n}: A={A_v:+.6f}  B={B_v:+.6f}  C={C_v:+.6f}  max_res={max_res:.4e}")
        if K_min == 1000:
            A_fit, B_fit, C_fit = A_v, B_v, C_v
    log("")

    # ------------------------------------------------------------------
    # Step 6: structural candidate test
    # ------------------------------------------------------------------
    log("-" * 78)
    log("Step 6: test structural candidates for log^2 coefficient A")
    log("-" * 78)
    log(f"  Best A estimate (3-param fit, K>=1000): A = {A_fit:.6f}")
    log("")
    candidates_A = [
        ("-C_2^2 / 2",        -C2_SQ / 2.0),
        ("-C_2^2",            -C2_SQ),
        ("-2 C_2^2",          -2.0 * C2_SQ),
        ("-C_2 / 3",          -C_2 / 3.0),
        ("-C_2 / 2",          -C_2 / 2.0),
        ("-C_2 / pi",         -C_2 / math.pi),
        ("-1/4",              -0.25),
        ("-1/(2 log 6)",      -1.0 / (2*math.log(6))),
        ("-log(6) C_2^2 / pi", -math.log(6) * C2_SQ / math.pi),
        ("-1/pi",             -1.0 / math.pi),
    ]
    if not math.isnan(A_fit):
        log(f"  {'Candidate (for A)':<35s} {'value':>12s} {'rel_diff':>12s}")
        for name, val in candidates_A:
            rel = abs(val - A_fit) / abs(A_fit) if A_fit != 0 else float("nan")
            log(f"  {name:<35s} {val:>12.6f} {rel:>12.2%}")
    log("")
    log("")

    # ------------------------------------------------------------------
    # Step 7: alpha "rolling" estimate over disjoint w bands
    # ------------------------------------------------------------------
    log("-" * 78)
    log("Step 7: rolling alpha estimate alpha(K_lo, K_hi) using band fits")
    log("-" * 78)
    bands = [(50, 200), (200, 500), (500, 1000), (1000, 5000),
             (5000, 20000), (20000, 50000), (50000, 200000)]
    alpha_band_rows = []
    for K_lo, K_hi in bands:
        band = [r for r in T_rows if K_lo <= r["K"] <= K_hi]
        if len(band) < 2:
            continue
        n = len(band)
        xs = [r["log_w"] for r in band]
        ys = [r["Delta"] / r["w"] for r in band]
        sx = sum(xs); sy = sum(ys)
        sxx = sum(x*x for x in xs)
        sxy = sum(x*y for x, y in zip(xs, ys))
        slope = (n * sxy - sx * sy) / (n * sxx - sx * sx)
        intercept = (sy - slope * sx) / n
        alpha_band = -slope
        alpha_band_rows.append({
            "K_lo": K_lo, "K_hi": K_hi, "n": n,
            "alpha": alpha_band, "beta": intercept,
        })
    log(f"  {'K range':<20s} {'n':>3s} {'alpha':>10s} {'beta':>10s}")
    for r in alpha_band_rows:
        rng = f"{r['K_lo']}..{r['K_hi']}"
        log(f"  {rng:<20s} {r['n']:>3d} {r['alpha']:>10.4f} {r['beta']:>10.4f}")
    log("")
    fit_path = out_dir / "slot4_alpha_fits.csv"
    with fit_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["K_lo", "K_hi", "n", "alpha", "beta"])
        w.writeheader()
        for r in alpha_band_rows:
            w.writerow({k: (f"{v:.6f}" if isinstance(v, float) else v) for k, v in r.items()})
    log(f"  alpha-band fits -> {fit_path}")
    log("")

    # ------------------------------------------------------------------
    # Final summary
    # ------------------------------------------------------------------
    log("=" * 78)
    log("Slot 4 summary")
    log("=" * 78)
    log(f"  H residual fit (clean):  alpha = {-a_h:.4f}, beta_H = {b_h:.4f}")
    log(f"  Delta fit (large w):     alpha = {alpha_T:.4f}, beta_Delta = {beta_T:.4f}")
    log(f"  Slot 3 fit (small w):    alpha = 5.7165, beta_Delta = 10.4958")
    log(f"  Candidate 24 C_2^2/log(6) = {TARGET_S4_AVG/math.log(6):.4f}")
    log("")

    log_path = out_dir / "slot4_run.log"
    log_path.write_text("\n".join(log_lines) + "\n")
    print(f"\nLog -> {log_path}")


if __name__ == "__main__":
    main()
