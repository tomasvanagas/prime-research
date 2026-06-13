"""
slot3_pair_correlation.py — Thread 9 (P4) commit slot 3

Goldston-Montgomery-style pair-correlation analysis of the cross-x HL
residual variance for windowed twin-prime counts.

The slot 2 empirical finding was:
  F_HL_kt = sigma_eff(N_2(x; w)) / sqrt(E[N_2(x; w)]) = 0.87 +/- 0.03
across three decades x in {10^6, 10^7, 10^8} and two w-regimes.

This slot derives the structure of F_HL_kt under HL conjecture, separating:
  - LEADING-ORDER cancellation (proved: F^2 -> 1 as x -> infty)
  - SUB-LEADING correction (computed: F^2 = 1 + Delta(w)/(C_2 * w * log^2 x))

Setup. Let Z_i = 1[i prime AND i+2 prime] for i in [x, x+w-1] and
N_2(x; w) = sum Z_i. Under HL, p := E[Z_i] = 2 C_2 / log^2(i) (constant
to leading order across the window).

Variance decomposition:
  Var[N_2] = E[N_2] + 2 sum_{m=1}^{w-1} (w-m) [ E[Z_i Z_{i+m}] - p^2 ]

The cross-term E[Z_i Z_{i+m}] = Pr(i, i+2, i+m, i+m+2 all prime). Under
HL 4-tuple, this equals S_4({0,2,m,m+2}) / log^4(x) for admissible
4-tuples, and 0 otherwise.

Define
  T(w) := sum_{m=1}^{w-1} (w-m) * S_4({0,2,m,m+2})
       = (1/log^4 x) * sum_{m} (w-m) E[Z_i Z_{i+m}]  scaled.

KEY IDENTITY (proved in this slot, via prime-by-prime averaging):
  Average of S_4({0,2,m,m+2}) / S_2({0,2})^2 over m mod p (each p >= 5)
  equals exactly 1.

This identity (combined with the p=2, p=3 admissibility factors of 2*3=6)
gives <S_4>_admissible = 6 * S_2^2 = 24 C_2^2, and the asymptotic

  T(w) ~ 2 C_2^2 * w^2  +  Delta(w)

where the leading 2 C_2^2 w^2 EXACTLY cancels the -p^2 * w^2 term in the
variance, yielding Var = E[N_2] * (1 + o(1)) to leading order.

The empirical F_HL_kt = 0.87 < 1 is NOT a contradiction; it is the
sub-leading Delta(w) term, which we compute numerically here.

Three deliverables:

  (1) The prime-by-prime cancellation identity (computed and verified
      symbolically/numerically for each p in [5, ..., P_MAX]).

  (2) Direct numerical evaluation of T(w) via singular-series
      computation for each m up to w. We then compute Delta(w) explicitly
      and predict F^2_pred(x; w) = 1 + Delta(w) / (C_2 * w * log^2 x).

  (3) Comparison of F^2_pred to the slot 2 empirical F^2_eff = (sigma_eff
      / sigma_pois)^2 at six (anchor, w_regime) cells.

If F^2_pred matches F^2_emp within tolerance across cells, slot 3 is
A-grade-shaped: a *new conditional theorem* for the cross-x k=2 windowed
variance factor under HL.

If the predicted Delta(w) does NOT match empirical (e.g., wrong sign,
wrong order of magnitude), then there is a missing term in our pair-
correlation analysis (probably from the off-diagonal contribution that
GM 1987 needs RH+pair-correlation conjecture for); we report that as
an honest negative.

Cross-domain ingredient: pair-correlation analysis a la
Goldston-Montgomery 1987 ("Pair correlation of zeros and primes in
short intervals"), applied here to the 4-tuple singular series.
NEW USE in CROSS_DOMAIN_TECHNIQUES.md (per-x windowed counts).

Edges cited:
  - Hardy-Littlewood 1923 k-tuple conjecture (provides S_4)
  - Gallagher 1976 (primes in short intervals are Poisson under HL)
  - Goldston-Montgomery 1987 (pair correlation, primes in short intervals)
  - S429, S430 (Thread 9 slots 1, 2)
  - S195/S224 (Thread 5 reference; pair correlation cancellation analog)
"""

from __future__ import annotations

import csv
import math
import os
import time
from pathlib import Path

# Pre-sieved primes up to PRIME_BOUND for singular-series evaluation.
PRIME_BOUND = 50_000


def sieve_primes(N: int) -> list[int]:
    is_p = bytearray(b"\x01") * (N + 1)
    is_p[0] = is_p[1] = 0
    bound = int(N ** 0.5) + 1
    for i in range(2, bound):
        if is_p[i]:
            is_p[i * i : N + 1 : i] = b"\x00" * (((N - i * i) // i) + 1)
    return [i for i in range(2, N + 1) if is_p[i]]


_PRIMES_CACHE: list[int] | None = None


def primes_up_to() -> list[int]:
    global _PRIMES_CACHE
    if _PRIMES_CACHE is None:
        _PRIMES_CACHE = sieve_primes(PRIME_BOUND)
    return _PRIMES_CACHE


# ---------------------------------------------------------------------------
# Singular series
# ---------------------------------------------------------------------------

def singular_series(H: tuple[int, ...], P_max: int = PRIME_BOUND) -> float:
    """S_k(H) = prod_p (1 - nu_H(p)/p) / (1-1/p)^k."""
    H = list(H)
    k = len(H)
    val = 1.0
    diam = max(H) - min(H)
    for p in primes_up_to():
        if p > P_max:
            break
        residues = set(h % p for h in H)
        nu = len(residues)
        if nu == p:
            return 0.0
        if p > diam and nu == k:
            # All remaining primes have nu = k (all residues distinct);
            # contribute factor 1 - k/p / (1-1/p)^k each. Compute tail
            # explicitly to save time.
            remaining = [q for q in primes_up_to() if q > p and q <= P_max]
            tail_log = 0.0
            for q in remaining:
                # Same factor for all q with nu = k:
                factor = (1.0 - k / q) / ((1.0 - 1.0 / q) ** k)
                tail_log += math.log(factor)
            val *= math.exp(tail_log)
            return val
        val *= (1.0 - nu / p) / ((1.0 - 1.0 / p) ** k)
    return val


# Twin-prime constant C_2 = prod_{p>=3} (1 - 1/(p-1)^2)
def twin_prime_constant() -> float:
    val = 1.0
    for p in primes_up_to():
        if p < 3:
            continue
        val *= (1.0 - 1.0 / (p - 1) ** 2)
    return val


C_2 = twin_prime_constant()
S_2_TWIN = 2.0 * C_2  # S_2({0, 2})


# ---------------------------------------------------------------------------
# Prime-by-prime cancellation identity
# ---------------------------------------------------------------------------

def factor_S4_at_p(p: int, m_mod_p: int) -> float:
    """Factor of S_4({0,2,m,m+2}) at prime p, given m mod p."""
    residues = {0 % p, 2 % p, m_mod_p % p, (m_mod_p + 2) % p}
    nu = len(residues)
    if nu == p:
        return 0.0
    return (1.0 - nu / p) / ((1.0 - 1.0 / p) ** 4)


def factor_S2_squared_at_p(p: int) -> float:
    """Factor of S_2({0,2})^2 at prime p."""
    residues = {0 % p, 2 % p}
    nu = len(residues)
    factor_S2 = (1.0 - nu / p) / ((1.0 - 1.0 / p) ** 2)
    return factor_S2 * factor_S2


def verify_cancellation_identity(P_max: int = 1000) -> list[dict]:
    """Verify that <S_4 factor at p> averaged over m mod p equals
    S_2^2 factor at p, for each prime p >= 5.

    For p = 2: S_4 factor (m even) / S_2^2 factor = 8 / 4 = 2.
    For p = 3: S_4 factor (m == 0 mod 3) / S_2^2 factor = (27/16)/(9/16) = 3.
    For p >= 5: <S_4 factor>_uniform / S_2^2 factor = 1 EXACTLY.

    Combined factor 6 = 2 * 3 (from admissibility constraint at small
    primes) * 1 (from each p >= 5).

    Hence <S_4>_admissible = 6 * S_2^2 = 24 C_2^2.

    We verify this identity numerically for p in [5, 7, ..., up to P_max].
    """
    out = []
    for p in primes_up_to():
        if p > P_max:
            break
        if p < 5:
            # Special cases
            if p == 2:
                # m even: ratio = 2
                f4 = factor_S4_at_p(p, 0)
                f2sq = factor_S2_squared_at_p(p)
                out.append({
                    "p": p, "case": "m even", "f4": f4, "f2sq": f2sq,
                    "ratio": f4 / f2sq, "expected": 2.0,
                })
            elif p == 3:
                # m == 0 mod 3: ratio = 3
                f4 = factor_S4_at_p(p, 0)
                f2sq = factor_S2_squared_at_p(p)
                out.append({
                    "p": p, "case": "m == 0 mod 3", "f4": f4, "f2sq": f2sq,
                    "ratio": f4 / f2sq, "expected": 3.0,
                })
            continue
        # p >= 5: average over m mod p uniformly
        avg_f4 = sum(factor_S4_at_p(p, mm) for mm in range(p)) / p
        f2sq = factor_S2_squared_at_p(p)
        out.append({
            "p": p, "case": "avg over m mod p",
            "f4": avg_f4, "f2sq": f2sq,
            "ratio": avg_f4 / f2sq, "expected": 1.0,
        })
    return out


# ---------------------------------------------------------------------------
# Numerical T(w) for twin-prime tuple
# ---------------------------------------------------------------------------

def admissible_4tuple(m: int) -> bool:
    """{0, 2, m, m+2} is admissible iff its residues mod every prime
    p <= 4 don't cover all residues mod p. For 4-tuples involving 0,2,
    only p = 3 is relevant: need m ≡ 0 (mod 3) (i.e., m ≡ 0 mod 6 since
    m must also be even for nontrivial pair). Strictly, m ≡ 0 mod 6.
    Also m must be even (otherwise some element is odd != 2)."""
    if m % 2 != 0 or m % 3 != 0:
        return False
    return True


def T_for_w(w: int, P_max: int = 5000, log_fn=print) -> tuple[float, list[float]]:
    """Compute T(w) = sum_{m=1}^{w-1} (w-m) S_4({0,2,m,m+2}).

    For inadmissible m, S_4 = 0. So we only sum over m ≡ 0 mod 6 (and we
    need m >= 6, since for m in {0, 2, 4} the tuple has duplicates or is
    inadmissible mod 3 or mod 2).

    Returns (T, S4_values) where S4_values[m // 6] is S_4 at m = 6k for
    k = 1, 2, ..., (w-1)//6.
    """
    # Iterate over m = 6, 12, 18, ..., up to w-1
    m_max = ((w - 1) // 6) * 6
    S4_values = []
    T = 0.0
    K = m_max // 6
    t0 = time.perf_counter()
    for k in range(1, K + 1):
        m = 6 * k
        s4 = singular_series((0, 2, m, m + 2), P_max=P_max)
        S4_values.append(s4)
        T += (w - m) * s4
    elapsed = time.perf_counter() - t0
    log_fn(f"    T(w={w}): K={K} admissible m's, T={T:.4f}, elapsed={elapsed:.2f}s")
    return T, S4_values


# ---------------------------------------------------------------------------
# Variance prediction
# ---------------------------------------------------------------------------

def predict_F_squared(x: int, w: int, T_w: float) -> dict:
    """Given T(w) computed above, predict F^2 = Var/E[N] where
    N = N_2(x; w).

    Formula:
      Var[N] ≈ E[N] + 2 * [T(w)/log^4(x) - 2 C_2^2 w^2 / log^4(x)]
             = E[N] + 2 * Delta(w) / log^4(x)
      Where Delta(w) = T(w) - 2 C_2^2 w^2.

      F^2 = Var/E[N] = 1 + 2 Delta(w) / (E[N] * log^4 x)
                    = 1 + Delta(w) / (C_2 * w * log^2 x).
    """
    log_x = math.log(x + w / 2.0)
    leading = 2.0 * C_2 * C_2 * w * w
    Delta = T_w - leading
    E_N = S_2_TWIN * w / (log_x ** 2)  # = 2 C_2 w / log^2(x)
    # F^2 = Var/E[N] = 1 + 2 Delta / (E[N] * log^4 x)
    #               = 1 + Delta / (C_2 * w * log^2 x).
    F_sq_pred = 1.0 + 2.0 * Delta / (E_N * log_x ** 4)
    F_sq_alt = 1.0 + Delta / (C_2 * w * log_x ** 2)  # equivalent form
    return {
        "x": x,
        "w": w,
        "log_x_mid": log_x,
        "T_w": T_w,
        "leading": leading,
        "Delta": Delta,
        "E_N": E_N,
        "F_sq_pred": F_sq_pred,
        "F_sq_check": F_sq_alt,  # should equal F_sq_pred
        "F_pred": math.sqrt(max(F_sq_pred, 0.0)),
    }


# ---------------------------------------------------------------------------
# Slot 2 empirical loader
# ---------------------------------------------------------------------------

def load_slot2_summary() -> list[dict]:
    """Load slot 2 per-cell summary."""
    path = Path(__file__).parent / "slot2_summary.csv"
    out = []
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            out.append({
                "anchor_log10": int(row["anchor_log10"]),
                "k_label": row["k_label"],
                "k": int(row["k"]),
                "C_H": float(row["C_H"]),
                "w_label": row["w_label"],
                "w": int(row["w"]),
                "N": int(row["N"]),
                "mean_HL": float(row["mean_HL"]),
                "sigma_eff": float(row["sigma_eff"]),
                "sigma_pois": float(row["sigma_pois"]),
                "F_emp": float(row["sigma_eff_over_pois"]),
            })
    return out


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def main():
    out_dir = Path(__file__).parent
    log_path = out_dir / "slot3_run.log"
    log_lines: list[str] = []

    def log_fn(msg: str = ""):
        print(msg)
        log_lines.append(str(msg))

    log_fn("=" * 78)
    log_fn("Thread 9 / P4 slot 3: Pair-correlation derivation of F_HL_kt")
    log_fn("=" * 78)
    log_fn(f"C_2 (twin-prime constant) = {C_2:.10f}")
    log_fn(f"S_2(0,2) = 2 C_2 = {S_2_TWIN:.10f}")
    log_fn(f"S_2(0,2)^2 = {S_2_TWIN**2:.10f}")
    log_fn(f"24 C_2^2 (predicted <S_4>_admissible) = {24 * C_2 * C_2:.10f}")
    log_fn("")

    # ------------------------------------------------------------------
    # 1. Verify cancellation identity prime-by-prime
    # ------------------------------------------------------------------
    log_fn("-" * 78)
    log_fn("1. Verify <S_4 at p>_uniform_m = S_2^2 at p, for each prime p")
    log_fn("-" * 78)
    log_fn("   At p=2 (m even):       ratio should be 2")
    log_fn("   At p=3 (m ≡ 0 mod 3):  ratio should be 3")
    log_fn("   At p>=5 (uniform m):   ratio should be 1 exactly")
    log_fn("")
    identity_rows = verify_cancellation_identity(P_max=1000)
    max_dev_p_ge_5 = 0.0
    for r in identity_rows[:8]:  # show first few
        log_fn(f"   p={r['p']:>4d}  case='{r['case']}'  f4={r['f4']:.6f}  f2sq={r['f2sq']:.6f}  "
               f"ratio={r['ratio']:.6f}  expected={r['expected']:.4f}  "
               f"dev={abs(r['ratio'] - r['expected']):.2e}")
    # check overall max deviation for p >= 5
    for r in identity_rows:
        if r["p"] >= 5:
            dev = abs(r["ratio"] - r["expected"])
            if dev > max_dev_p_ge_5:
                max_dev_p_ge_5 = dev
    log_fn(f"   ...")
    log_fn(f"   Max deviation from expected for p >= 5: {max_dev_p_ge_5:.2e}")
    log_fn(f"   (numerical roundoff only; cancellation identity holds EXACTLY)")
    log_fn("")

    # Save identity verification to CSV
    id_path = out_dir / "slot3_identity.csv"
    with id_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=identity_rows[0].keys())
        w.writeheader()
        for r in identity_rows:
            w.writerow({k: (f"{v:.10f}" if isinstance(v, float) else v) for k, v in r.items()})
    log_fn(f"   Identity table -> {id_path} ({len(identity_rows)} primes)")
    log_fn("")

    # ------------------------------------------------------------------
    # 2. Compute T(w) and predict F^2 for slot 2 cells
    # ------------------------------------------------------------------
    log_fn("-" * 78)
    log_fn("2. Direct evaluation of T(w) = sum_{admissible m=6,12,...,<w} (w-m) S_4(0,2,m,m+2)")
    log_fn("-" * 78)
    cells = load_slot2_summary()
    twin_cells = [c for c in cells if c["k_label"] == "k2_twin"]
    log_fn(f"   Loaded {len(twin_cells)} k=2 (twin) cells from slot2_summary.csv")
    log_fn("")

    # Compute T(w) once per unique w (not per cell, since T only depends on w)
    unique_ws = sorted(set(c["w"] for c in twin_cells))
    log_fn(f"   Unique w values: {unique_ws}")
    log_fn("")
    T_by_w: dict[int, float] = {}
    S4_by_w: dict[int, list[float]] = {}
    for w in unique_ws:
        log_fn(f"   Computing T(w={w}) ...")
        Tw, S4s = T_for_w(w, P_max=2000, log_fn=log_fn)
        T_by_w[w] = Tw
        S4_by_w[w] = S4s
    log_fn("")

    # Predict F^2 per cell
    log_fn("-" * 78)
    log_fn("3. F^2 prediction vs slot 2 empirical (k=2 twin cells)")
    log_fn("-" * 78)
    predict_rows = []
    for c in twin_cells:
        x = c["anchor_log10"]
        x_val = 10 ** x
        w = c["w"]
        Tw = T_by_w[w]
        pred = predict_F_squared(x_val, w, Tw)
        F_emp = c["F_emp"]
        F_sq_emp = F_emp * F_emp
        row = {
            "anchor_log10": x,
            "x": x_val,
            "k_label": c["k_label"],
            "w_label": c["w_label"],
            "w": w,
            "T_w": Tw,
            "leading_2C2sq_wsq": pred["leading"],
            "Delta_w": pred["Delta"],
            "Delta_over_leading": pred["Delta"] / pred["leading"],
            "E_N": pred["E_N"],
            "F_sq_pred": pred["F_sq_pred"],
            "F_pred": pred["F_pred"],
            "F_sq_emp": F_sq_emp,
            "F_emp": F_emp,
            "F_diff": pred["F_pred"] - F_emp,
            "F_sq_diff": pred["F_sq_pred"] - F_sq_emp,
        }
        predict_rows.append(row)

    # Print table
    log_fn(f"{'anc':>4s} {'w_lab':>7s} {'w':>5s} {'T(w)':>12s} {'leading':>12s} "
           f"{'Delta':>11s} {'Δ/lead':>9s} {'F_pred':>8s} {'F_emp':>8s} {'F_diff':>9s}")
    log_fn("-" * 110)
    for r in predict_rows:
        log_fn(f"{r['anchor_log10']:>4d} {r['w_label']:>7s} {r['w']:>5d} "
               f"{r['T_w']:>12.4f} {r['leading_2C2sq_wsq']:>12.4f} "
               f"{r['Delta_w']:>11.4f} {r['Delta_over_leading']:>9.4f} "
               f"{r['F_pred']:>8.4f} {r['F_emp']:>8.4f} {r['F_diff']:>9.4f}")
    log_fn("=" * 110)
    log_fn("")

    # Compute aggregate diagnostics
    F_diff_abs = [abs(r["F_diff"]) for r in predict_rows]
    log_fn(f"  |F_pred - F_emp|: max={max(F_diff_abs):.4f}, mean={sum(F_diff_abs)/len(F_diff_abs):.4f}")
    F_diffs = [r["F_diff"] for r in predict_rows]
    sign_pos = sum(1 for d in F_diffs if d > 0)
    log_fn(f"  F_pred > F_emp in {sign_pos}/{len(F_diffs)} cells "
           f"(positive => prediction overshoots empirical reduction)")
    log_fn("")

    # Save predictions
    pred_path = out_dir / "slot3_predictions.csv"
    with pred_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=predict_rows[0].keys())
        w.writeheader()
        for r in predict_rows:
            w.writerow({k: (f"{v:.6f}" if isinstance(v, float) else v) for k, v in r.items()})
    log_fn(f"  Predictions -> {pred_path}")

    # ------------------------------------------------------------------
    # 3. Per-m profile of S_4 at finest resolution (for w=339, smallest)
    # ------------------------------------------------------------------
    log_fn("")
    log_fn("-" * 78)
    log_fn("4. Per-m profile of S_4(0,2,m,m+2) for first 30 admissible m's")
    log_fn("    (mean of S_4 over ALL admissible m -> 24 C_2^2 = "
           f"{24*C_2*C_2:.6f})")
    log_fn("-" * 78)
    smallest_w = min(unique_ws)
    S4s = S4_by_w[smallest_w][:30]
    log_fn(f"{'k':>3s} {'m':>5s} {'S_4(0,2,m,m+2)':>17s} {'cum_avg':>10s}")
    cum_sum = 0.0
    for k, s4 in enumerate(S4s, start=1):
        cum_sum += s4
        cum_avg = cum_sum / k
        log_fn(f"{k:>3d} {6*k:>5d} {s4:>17.6f} {cum_avg:>10.6f}")
    log_fn("")

    # Save S_4 profile (full)
    s4_path = out_dir / "slot3_s4_profile.csv"
    with s4_path.open("w", newline="") as f:
        ww = csv.writer(f)
        ww.writerow(["k", "m", "S_4(0,2,m,m+2)"])
        for k, s4 in enumerate(S4_by_w[max(unique_ws)], start=1):
            ww.writerow([k, 6 * k, f"{s4:.6f}"])
    log_fn(f"  S_4 profile -> {s4_path} (full, w_max={max(unique_ws)})")

    # ------------------------------------------------------------------
    # 4. Cumulative <S_4> convergence diagnostic
    # ------------------------------------------------------------------
    log_fn("")
    log_fn("-" * 78)
    log_fn("5. Cumulative <S_4> over admissible m: convergence to 24 C_2^2")
    log_fn("-" * 78)
    target = 24 * C_2 * C_2
    log_fn(f"   Target <S_4>_admissible = 24 C_2^2 = {target:.6f}")
    log_fn(f"{'M_max':>8s} {'<S_4>':>12s} {'rel_err':>10s}")
    s4_full = S4_by_w[max(unique_ws)]
    cum = 0.0
    for k, s4 in enumerate(s4_full, start=1):
        cum += s4
        if k in [1, 5, 10, 20, 50, 100, 200, 500, 1000] or k == len(s4_full):
            avg = cum / k
            rel = abs(avg - target) / target
            log_fn(f"{k:>8d} {avg:>12.6f} {rel:>10.4f}")
    log_fn("")

    # Write log
    log_path.write_text("\n".join(log_lines) + "\n")
    log_fn(f"\nLog -> {log_path}")
    print(f"\nLog -> {log_path}")


if __name__ == "__main__":
    main()
