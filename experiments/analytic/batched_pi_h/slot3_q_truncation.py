"""
slot3_q_truncation.py — Thread 8 (P2) commit slot 3

Q-truncated singular series cost-vs-error tradeoff for the Hardy-
Littlewood approximation HL_h(x) = S_h * li_2(x) of the prime-gap
function pi_h(x) = #{p <= x : p, p+h both prime}.

Define
   S_Q(h) = 2 C_2 * Pi_{p odd, p|h, p<=Q} (p-1)/(p-2)
i.e. include only the odd prime factors of h that are <= Q. The full
singular series is S_inf(h) (no truncation).

Slot 1 (S418): structural dichotomy — APPROX regime is polylog-per-h,
mean |err|/sqrt(x) <= 0.10.
Slot 2 (S419): cross-h ensemble at fixed (anchor, x) is the natural
sampling for HL residual analysis; sigma_eff/sigma_pois in [0.36, 0.70];
half-Gaussian shape passes 90/90 cells.

Slot 3 mission: characterise the Q-truncation tradeoff.
  - Per-h truncation error: eps_Q(h) = |S_inf(h) - S_Q(h)|.
  - Cross-h aggregate: sigma_HL_Q(x) = RMS over h-ensemble of
    r_Q(x, h) = pi_h(x) - S_Q(h) * li_2(x).
  - Knee: smallest Q* s.t. sigma_HL_Q ~ sigma_HL_inf within tolerance.
  - Cost: trial division up to min(sqrt(h), Q) -> ~ Q operations per h
    when Q < sqrt(h_max). For h_max <= polylog(x), full S_inf computation
    is already polylog. Q-truncation is interesting when h is large.

A-grade target (per .commit_state recommended_next_action): a precise
tradeoff curve Q -> per-h cost vs total error sigma_HL_Q, with named-
exponent error formula
  sigma_HL_Q(x, h-ensemble) = sqrt[ sigma_eff_inf^2 + (li_2(x) * RMS[eps_Q])^2 ]
i.e. quadrature of intrinsic (Q-independent) and truncation contributions.

If knee is at Q << sqrt(h_max), polylog HL evaluation with named-exponent
error is established (Thread-7 Corollary B analogue).

Edges cited:
  - E1.5 (information density of pi)
  - S195 (variance formula machinery)
  - S224 (Correlation Dichotomy template)
  - S418 / S419 slot-1, slot-2 baselines
  - S240 / S244 Thread 7 named-exponent corollary (analogue)

Falsification:
  - If sigma_HL_Q monotone-decreases with Q across all measured Q (no
    plateau), then no knee exists in the tested range and the bound is
    not Q-truncation-friendly.
  - If sigma_HL_Q at Q = max_p_h(h_set) exceeds sigma_HL_inf by more
    than 10%, then the cross-h decomposition fails (truncation error
    dominates intrinsic noise even for the largest tested Q).
"""

from __future__ import annotations

import csv
import math
import time
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
C_2 = 0.6601618158468695739278121100145557784326233628

# h-values with diverse max prime factors spanning [3, 2003]:
H_VALUES = [
    2,           # max_p = - (no odd primes)
    6,           # max_p = 3
    14,          # max_p = 7
    22,          # max_p = 11
    30,          # max_p = 5  (= 2 * 3 * 5)
    38,          # max_p = 19 (= 2 * 19)
    42,          # max_p = 7  (= 2 * 3 * 7)
    46,          # max_p = 23 (= 2 * 23)
    58,          # max_p = 29 (= 2 * 29)
    62,          # max_p = 31 (= 2 * 31)
    78,          # max_p = 13 (= 2 * 3 * 13)
    82,          # max_p = 41 (= 2 * 41)
    106,         # max_p = 53 (= 2 * 53)
    158,         # max_p = 79 (= 2 * 79)
    198,         # max_p = 11 (= 2 * 3^2 * 11)
    218,         # max_p = 109
    302,         # max_p = 151
    398,         # max_p = 199
    502,         # max_p = 251
    606,         # max_p = 101 (= 2 * 3 * 101)
    802,         # max_p = 401
    1006,        # max_p = 503
    1198,        # max_p = 599
    2018,        # max_p = 1009
    3030,        # max_p = 101 (= 2 * 3 * 5 * 101)
    4006,        # max_p = 2003
]

Q_VALUES_FINITE = [30, 50, 100, 200, 500, 1000, 2000, 5000]
Q_INF_LABEL = "inf"

ANCHORS_LOG10 = [7, 8]
N_X_PER_ANCHOR = 5  # log-uniform in [x_anchor, x_anchor * 10^0.25]
WINDOW_LOG = 0.25 * math.log(10.0)


# ---------------------------------------------------------------------------
# Sieve + factorisation primitives
# ---------------------------------------------------------------------------


def sieve_eratosthenes(N: int) -> bytearray:
    is_p = bytearray(b"\x01") * (N + 1)
    is_p[0] = 0
    is_p[1] = 0
    bound = int(N ** 0.5) + 1
    for i in range(2, bound):
        if is_p[i]:
            start = i * i
            is_p[start : N + 1 : i] = b"\x00" * (((N - start) // i) + 1)
    return is_p


def li2(x: float) -> float:
    from scipy.integrate import quad
    v, _ = quad(lambda t: 1.0 / math.log(t) ** 2, 2.0, float(x), limit=400)
    return v


def factor_h_odd_primes(h: int) -> list[int]:
    """Sorted list of distinct odd prime factors of h (returns [] if h has
    only the factor 2)."""
    primes: list[int] = []
    n = h
    while n % 2 == 0:
        n //= 2
    p = 3
    while p * p <= n:
        if n % p == 0:
            primes.append(p)
            while n % p == 0:
                n //= p
        p += 2
    if n > 1:
        primes.append(n)
    return primes


def s_q(h: int, Q: float, primes_h: list[int] | None = None) -> float:
    """Truncated singular series.

    S_Q(h) = 2 C_2 * Pi_{p in primes_h, p <= Q} (p - 1)/(p - 2).

    Q can be math.inf for the full singular series.
    """
    if h % 2 != 0:
        return 0.0
    if primes_h is None:
        primes_h = factor_h_odd_primes(h)
    factor = 1.0
    for p in primes_h:
        if p <= Q:
            factor *= (p - 1) / (p - 2)
    return 2.0 * C_2 * factor


def s_q_cost_h(h: int, Q: float) -> int:
    """Number of trial-division steps to compute S_Q(h).

    Approximates the cost as # odd integers tested up to min(Q, sqrt(h)).
    For h <= Q^2: cost <= sqrt(h) / 2.
    For Q < sqrt(h): cost <= Q / 2.
    """
    bound = min(int(math.sqrt(h)), int(Q if Q != math.inf else math.sqrt(h))) + 1
    if bound < 3:
        return 0
    # # odd integers in [3, bound]
    return max(0, (bound - 3) // 2 + 1)


# ---------------------------------------------------------------------------
# Sample design
# ---------------------------------------------------------------------------


def log_uniform_xs(x_anchor: int, N: int) -> list[int]:
    log_lo = math.log(float(x_anchor))
    log_hi = log_lo + WINDOW_LOG
    raw = [int(round(math.exp(log_lo + j * (log_hi - log_lo) / (N - 1)))) for j in range(N)]
    seen, out = set(), []
    for x in raw:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return sorted(out)


# ---------------------------------------------------------------------------
# Per-anchor evaluator
# ---------------------------------------------------------------------------


def evaluate_anchor(anchor_log10: int) -> tuple[dict[tuple[int, int], int], dict[int, float]]:
    """Returns (pi_h_dict, li2_xs).

    pi_h_dict[(x, h)] = exact pi_h(x) for all x in xs, h in H_VALUES.
    li2_xs[x] = li_2(x) value.
    """
    x_anchor = 10 ** anchor_log10
    xs = log_uniform_xs(x_anchor, N_X_PER_ANCHOR)
    x_max = max(xs)
    h_max = max(H_VALUES)

    print(f"\n[anchor x = 10^{anchor_log10}]  N_x={len(xs)}  x_max={x_max:,}  h_max={h_max}")

    t0 = time.perf_counter()
    sieve = sieve_eratosthenes(x_max + h_max)
    t_sieve = time.perf_counter() - t0
    print(f"  sieve [1, {x_max + h_max:,}] in {t_sieve:.2f}s")
    is_p = np.frombuffer(sieve, dtype=np.uint8)

    t0 = time.perf_counter()
    li2_xs = {x: li2(x) for x in xs}
    print(f"  li_2 evaluations: {time.perf_counter() - t0:.2f}s")

    pi_h_dict: dict[tuple[int, int], int] = {}
    t0 = time.perf_counter()
    for h in H_VALUES:
        pair = is_p[: x_max + 1] & is_p[h : x_max + 1 + h]
        # Walk through sorted xs accumulating the prefix sum.
        cum = 0
        last_idx = 0
        for x in xs:
            cum += int(pair[last_idx : x + 1].sum())
            last_idx = x + 1
            pi_h_dict[(x, h)] = cum
        del pair
    print(f"  pair-count over {len(H_VALUES)} h-values: {time.perf_counter() - t0:.2f}s")
    return pi_h_dict, li2_xs


# ---------------------------------------------------------------------------
# Statistics: KS test against half-Gaussian
# ---------------------------------------------------------------------------


def half_normal_cdf(z: float) -> float:
    return math.erf(z / math.sqrt(2)) if z >= 0 else 0.0


def ks_stat(samples: list[float]) -> float:
    n = len(samples)
    if n == 0:
        return float("nan")
    s = sorted(samples)
    D = 0.0
    for i, z in enumerate(s):
        F_emp_lo = i / n
        F_emp_hi = (i + 1) / n
        F_th = half_normal_cdf(z)
        D = max(D, abs(F_emp_lo - F_th), abs(F_emp_hi - F_th))
    return D


def ks_pvalue(D: float, n: int) -> float:
    if n < 2:
        return float("nan")
    en = math.sqrt(n) * D + 0.12 * D + 0.11 * D / math.sqrt(n)
    s = 0.0
    for j in range(1, 50):
        term = 2 * (-1) ** (j - 1) * math.exp(-2 * j * j * en * en)
        s += term
        if abs(term) < 1e-12:
            break
    return min(1.0, max(0.0, s))


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------


def main() -> None:
    out_dir = Path(__file__).parent

    # Pre-compute factorisations and S_Q for each (h, Q).
    primes_by_h = {h: factor_h_odd_primes(h) for h in H_VALUES}
    max_p_by_h = {h: (max(primes_by_h[h]) if primes_by_h[h] else 0) for h in H_VALUES}
    s_q_table: dict[tuple[int, float], float] = {}
    for h in H_VALUES:
        s_q_table[(h, math.inf)] = s_q(h, math.inf, primes_by_h[h])
        for Q in Q_VALUES_FINITE:
            s_q_table[(h, float(Q))] = s_q(h, float(Q), primes_by_h[h])

    # Truncation table: per-(h, Q), |S_inf(h) - S_Q(h)| and relative.
    trunc_path = out_dir / "slot3_truncation.csv"
    with trunc_path.open("w", newline="") as f:
        w = csv.writer(f)
        header = ["h", "max_p", "S_inf"]
        for Q in Q_VALUES_FINITE:
            header += [f"S_{Q}", f"abs_eps_{Q}", f"rel_eps_{Q}"]
        w.writerow(header)
        for h in H_VALUES:
            S_inf = s_q_table[(h, math.inf)]
            row: list = [h, max_p_by_h[h], f"{S_inf:.6f}"]
            for Q in Q_VALUES_FINITE:
                S_Q_val = s_q_table[(h, float(Q))]
                eps = abs(S_inf - S_Q_val)
                rel = eps / S_inf if S_inf > 0 else 0.0
                row += [f"{S_Q_val:.6f}", f"{eps:.6f}", f"{rel:.6e}"]
            w.writerow(row)
    print(f"Truncation table -> {trunc_path}")

    # Evaluate at each anchor, build full data row list.
    all_data: list[dict] = []
    for anc_log10 in ANCHORS_LOG10:
        pi_h, li2_xs = evaluate_anchor(anc_log10)
        for x_j, li2_xj in sorted(li2_xs.items()):
            for h in H_VALUES:
                pi_h_xj = pi_h[(x_j, h)]
                S_inf_h = s_q_table[(h, math.inf)]
                row = {
                    "anchor_log10": anc_log10,
                    "x": x_j,
                    "h": h,
                    "max_p_h": max_p_by_h[h],
                    "S_inf": S_inf_h,
                    "pi_h": pi_h_xj,
                    "li2_x": li2_xj,
                    "sqrt_x": math.sqrt(x_j),
                    "sigma_pred_pois": math.sqrt(S_inf_h * li2_xj) if S_inf_h > 0 else 0.0,
                }
                # Per-Q residual columns
                for Q in Q_VALUES_FINITE:
                    S_Q_val = s_q_table[(h, float(Q))]
                    HL_Q = S_Q_val * li2_xj
                    r_Q = pi_h_xj - HL_Q
                    row[f"S_{Q}"] = S_Q_val
                    row[f"r_{Q}"] = r_Q
                # Q = inf baseline
                HL_inf = S_inf_h * li2_xj
                row["r_inf"] = pi_h_xj - HL_inf
                all_data.append(row)

    # Per-row CSV
    raw_path = out_dir / "slot3_data.csv"
    with raw_path.open("w", newline="") as f:
        w = csv.writer(f)
        keys = list(all_data[0].keys())
        w.writerow(keys)
        for r in all_data:
            w.writerow([
                r[k] if isinstance(r[k], int)
                else f"{r[k]:.6f}" if isinstance(r[k], float)
                else r[k]
                for k in keys
            ])
    print(f"\nRaw rows ({len(all_data)}) -> {raw_path}")

    # ========== Cross-h aggregates per (anchor, x_j, Q) ==========
    cross_h_groups: dict[tuple[int, int], list[dict]] = {}
    for r in all_data:
        key = (r["anchor_log10"], r["x"])
        cross_h_groups.setdefault(key, []).append(r)

    cross_h_path = out_dir / "slot3_cross_h.csv"
    cross_h_rows = []
    for (anc, x_j), g in sorted(cross_h_groups.items()):
        N = len(g)
        for Q_label, Q_key in [(str(Q), f"r_{Q}") for Q in Q_VALUES_FINITE] + [("inf", "r_inf")]:
            residuals = [row[Q_key] for row in g]
            abs_residuals = [abs(r) for r in residuals]
            sigma_HL_Q = math.sqrt(sum(r * r for r in residuals) / N)
            mean_abs = sum(abs_residuals) / N
            median_abs = sorted(abs_residuals)[N // 2]

            # Normalised by sigma_pred_pois (per-h)
            norm_abs = [
                abs(row[Q_key]) / row["sigma_pred_pois"] if row["sigma_pred_pois"] > 0 else 0.0
                for row in g
            ]
            norm_signed = [
                row[Q_key] / row["sigma_pred_pois"] if row["sigma_pred_pois"] > 0 else 0.0
                for row in g
            ]
            sigma_norm = math.sqrt(sum(r * r for r in norm_signed) / N)
            # KS of |r|/sigma_HL_Q (cross-h ensemble, re-normalised by ensemble RMS)
            if sigma_HL_Q > 0:
                ratios_eff = [a / sigma_HL_Q for a in abs_residuals]
                D_eff = ks_stat(ratios_eff)
                p_eff = ks_pvalue(D_eff, N)
            else:
                D_eff, p_eff = float("nan"), float("nan")

            cross_h_rows.append({
                "anchor_log10": anc,
                "x": x_j,
                "Q": Q_label,
                "N_h": N,
                "sigma_HL_Q": sigma_HL_Q,
                "mean_abs": mean_abs,
                "median_abs": median_abs,
                "sigma_norm_pois": sigma_norm,
                "ks_D_eff": D_eff,
                "ks_p_eff": p_eff,
            })

    with cross_h_path.open("w", newline="") as f:
        w = csv.writer(f)
        keys = list(cross_h_rows[0].keys())
        w.writerow(keys)
        for r in cross_h_rows:
            w.writerow([
                r[k] if isinstance(r[k], (int, str))
                else f"{r[k]:.6f}"
                for k in keys
            ])
    print(f"Cross-h aggregates ({len(cross_h_rows)}) -> {cross_h_path}")

    # ========== Knee summary per (anchor, x_j) ==========
    print("\n" + "=" * 130)
    print("Cross-h sigma_HL_Q vs Q  (knee = smallest Q with sigma_HL_Q within 5% of sigma_HL_inf)")
    print("=" * 130)
    print(f"{'anc':>5s} {'x':>13s} | "
          + " ".join([f"Q={q:<5d}" for q in Q_VALUES_FINITE])
          + f"  Q=inf       knee_Q  knee_max_p")
    knee_summary = []
    for (anc, x_j), g in sorted(cross_h_groups.items()):
        sigmas_by_Q = {}
        for Q in Q_VALUES_FINITE + [math.inf]:
            r_key = "r_inf" if Q == math.inf else f"r_{Q}"
            residuals = [row[r_key] for row in g]
            sigmas_by_Q[Q] = math.sqrt(sum(r * r for r in residuals) / len(residuals))
        sigma_inf = sigmas_by_Q[math.inf]
        # Find smallest Q s.t. sigma_HL_Q <= 1.05 * sigma_inf
        knee = None
        for Q in Q_VALUES_FINITE:
            if sigmas_by_Q[Q] <= 1.05 * sigma_inf:
                knee = Q
                break
        # Also note the max max_p_h that is <= knee (most informative)
        knee_max_p = (
            max((max_p_by_h[h] for h in H_VALUES if max_p_by_h[h] <= knee), default=0)
            if knee else 0
        )
        knee_summary.append({
            "anchor_log10": anc,
            "x": x_j,
            "sigma_inf": sigma_inf,
            "knee_Q": knee,
            "knee_max_p": knee_max_p,
            **{f"sigma_Q_{Q}": sigmas_by_Q[Q] for Q in Q_VALUES_FINITE},
            "sigma_Q_inf": sigma_inf,
        })
        sigma_str = " ".join([f"{sigmas_by_Q[Q]:>7.0f}" for Q in Q_VALUES_FINITE])
        print(f"{anc:>5d} {x_j:>13d} | {sigma_str}  {sigma_inf:>7.0f}      {str(knee):>5s}    {knee_max_p:>3d}")

    knee_path = out_dir / "slot3_knee.csv"
    with knee_path.open("w", newline="") as f:
        w = csv.writer(f)
        keys = list(knee_summary[0].keys())
        w.writerow(keys)
        for r in knee_summary:
            w.writerow([
                r[k] if isinstance(r[k], (int, str))
                else (f"{r[k]:.6f}" if r[k] is not None and isinstance(r[k], float) else r[k])
                for k in keys
            ])
    print(f"\nKnee summary -> {knee_path}")

    # ========== Cost vs accuracy summary ==========
    # Per-h cost at Q for the largest h in our set (worst-case in our ensemble).
    h_max = max(H_VALUES)
    print("\n" + "=" * 80)
    print(f"Per-h cost at h_max = {h_max} (max in ensemble), trial-division steps:")
    print(f"{'Q':>10s} {'cost(h_max, Q)':>18s} {'cost(h, Q) avg over h':>26s}")
    for Q in Q_VALUES_FINITE + [math.inf]:
        c_max = s_q_cost_h(h_max, Q)
        c_avg = sum(s_q_cost_h(h, Q) for h in H_VALUES) / len(H_VALUES)
        Q_str = "inf" if Q == math.inf else str(int(Q))
        print(f"{Q_str:>10s} {c_max:>18d} {c_avg:>26.1f}")

    # ========== KS summary across all (anchor, x_j) cells per Q ==========
    print("\n" + "=" * 80)
    print("Cross-h KS (vs half-Gaussian) per Q across all (anchor, x_j) cells:")
    print(f"{'Q':>10s} {'cells':>8s} {'med p_eff':>12s} {'cells p_eff>0.1':>17s} {'min p_eff':>12s}")
    by_Q: dict[str, list[dict]] = {}
    for r in cross_h_rows:
        by_Q.setdefault(r["Q"], []).append(r)
    for Q_label in [str(Q) for Q in Q_VALUES_FINITE] + ["inf"]:
        cells = by_Q[Q_label]
        ps = sorted([c["ks_p_eff"] for c in cells])
        med_p = ps[len(ps) // 2] if ps else float("nan")
        n_pass = sum(1 for p in ps if p > 0.1)
        min_p = min(ps) if ps else float("nan")
        print(f"{Q_label:>10s} {len(cells):>8d} {med_p:>12.4f} {n_pass:>17d} {min_p:>12.4f}")


if __name__ == "__main__":
    main()
