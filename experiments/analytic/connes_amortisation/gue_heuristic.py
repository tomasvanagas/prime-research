"""
GUE random-phase heuristic for the explicit-formula error variance.
Slot 3/5 of commit thread "connes_amortisation" (S195).

Hypothesis (random-phase model for {gamma_j log x mod 2pi}):
   Var(pi(x) - pi_K(x)) ~ x * log^2(K) / (2 pi^2 K log^2(x))   (*)

Derivation sketch (in gue_heuristic_results.md):
   Each conjugate pair contributes 2 x^{1/2} sin(gamma_j log x) / (gamma_j log x)
   (asymptotic li expansion); under iid uniform phases Var(sin) = 1/2;
   sum 1/gamma_j^2 over j > K, using gamma_j ~ 2 pi j / log(j),
   evaluates to log^2(K) / (4 pi^2 K) + lower order.

This script:
  1. Loads the S194 fluctuation CSVs (40 samples each at x ~ 1e5, x ~ 1e6).
  2. For each (x, K) computes predicted sigma per (*).
  3. Compares predicted Gaussian-quantile median |err| against empirical
     mean / median |err| across the 40 samples.
  4. Compares predicted erf-hit-rate against empirical hit-rate.
  5. Sweeps the K-policy axis to identify the K* where hit-rate -> 1.

The script is self-contained (no zeros needed beyond the CSV).
"""

from __future__ import annotations

import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path

# ----------------------------------------------------------------------------
# The heuristic predictor
# ----------------------------------------------------------------------------


def sigma_predicted(x: float, K: int) -> float:
    """Predicted sigma for pi(x) - pi_K(x) under random-phase model.

    sigma^2 ~ x * log^2(K) / (2 pi^2 K log^2(x)).
    Returns 0 for K <= 1 (degenerate).
    """
    if K <= 1 or x <= 1:
        return float("inf")
    return (math.sqrt(x) * math.log(K)) / (math.pi * math.sqrt(2 * K) * math.log(x))


def median_predicted(x: float, K: int) -> float:
    """Predicted median |err| under Gaussian assumption: 0.6745 * sigma."""
    return 0.67448975 * sigma_predicted(x, K)


def hit_rate_predicted(x: float, K: int, threshold: float = 0.5) -> float:
    """Predicted P(|err| <= threshold) under N(0, sigma^2)."""
    sigma = sigma_predicted(x, K)
    if sigma == 0:
        return 1.0
    if math.isinf(sigma):
        return 0.0
    return math.erf(threshold / (math.sqrt(2) * sigma))


# ----------------------------------------------------------------------------
# CSV ingestion
# ----------------------------------------------------------------------------


def load_fluctuation_csv(path: Path):
    """Yield rows as dicts of floats."""
    with path.open() as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            yield {k: float(v) for k, v in row.items()}


# ----------------------------------------------------------------------------
# Comparison with S194 data
# ----------------------------------------------------------------------------


POLICIES = [
    ("log2x",       "K_log2x",       "abs_err_at_K_log2x",       "hit_at_K_log2x"),
    ("log3x",       "K_log3x",       "abs_err_at_K_log3x",       "hit_at_K_log3x"),
    ("5log2x",      "K_5log2x",      "abs_err_at_K_5log2x",      "hit_at_K_5log2x"),
    ("sqrtx",       "K_sqrtx",       "abs_err_at_K_sqrtx",       "hit_at_K_sqrtx"),
    ("halfsqrtx",   "K_halfsqrtx",   "abs_err_at_K_halfsqrtx",   "hit_at_K_halfsqrtx"),
]


def compare_csv(path: Path, label: str):
    rows = list(load_fluctuation_csv(path))
    print(f"\n=== {label}: {len(rows)} samples ===")
    print(f"x range: [{min(r['x'] for r in rows):.3g}, {max(r['x'] for r in rows):.3g}]")
    print()
    print(f"{'policy':10s} {'K_med':>7s} {'pred_sig':>9s} {'pred_med':>9s} "
          f"{'emp_med':>9s} {'pred_hit':>9s} {'emp_hit':>9s} {'ratio_sig':>10s}")
    print("-" * 80)

    for name, k_col, err_col, hit_col in POLICIES:
        Ks = [r[k_col] for r in rows]
        errs = [r[err_col] for r in rows]
        hits = [r[hit_col] for r in rows]
        K_med = statistics.median(Ks)
        x_med = statistics.median([r["x"] for r in rows])

        pred_sig = sigma_predicted(x_med, int(round(K_med)))
        pred_med = median_predicted(x_med, int(round(K_med)))
        pred_hit = hit_rate_predicted(x_med, int(round(K_med)))
        emp_med = statistics.median(errs)
        emp_hit = sum(hits) / len(hits)

        # Empirical sigma estimate from data: std of signed err. We only have
        # |err|, so use 1.4826 * MAD-of-|err| as a robust scale, or use rms.
        rms = math.sqrt(sum(e * e for e in errs) / len(errs))
        ratio_sig = rms / pred_sig if pred_sig > 0 else float("nan")

        print(f"{name:10s} {K_med:>7.0f} {pred_sig:>9.3f} {pred_med:>9.3f} "
              f"{emp_med:>9.3f} {pred_hit:>9.3f} {emp_hit:>9.3f} {ratio_sig:>10.3f}")


def per_sample_comparison(path: Path, label: str, policy: str = "log2x"):
    """Per-sample scatter of empirical |err| vs predicted sigma for one policy."""
    p = next(p for p in POLICIES if p[0] == policy)
    _, k_col, err_col, _ = p
    rows = list(load_fluctuation_csv(path))
    print(f"\n--- per-sample scatter ({label}, policy {policy}) ---")
    print(f"{'x':>12s} {'K':>5s} {'sigma_pred':>10s} {'|err_emp|':>10s} {'ratio':>8s}")
    for r in rows[:8]:
        K = int(r[k_col])
        sig = sigma_predicted(r["x"], K)
        emp = r[err_col]
        print(f"{r['x']:>12.4g} {K:>5d} {sig:>10.3f} {emp:>10.3f} {emp/sig:>8.3f}")


def aggregate_correlation(paths_labels):
    """Pool all 80 samples x 5 policies = 400 (x, K, |err|) triples and report:
       - Pearson correlation of |err| vs sigma_predicted
       - Slope of |err| vs sigma (best-fit through origin)
       - Half-Gaussian comparison: under N(0,sigma^2), E[|err|] = sigma * sqrt(2/pi)
    """
    triples = []  # (sigma_pred, |err_emp|, label, policy)
    for path, label in paths_labels:
        rows = list(load_fluctuation_csv(path))
        for r in rows:
            for name, k_col, err_col, _ in POLICIES:
                K = int(r[k_col])
                if K <= 1:
                    continue
                sig = sigma_predicted(r["x"], K)
                triples.append((sig, r[err_col], label, name))

    n = len(triples)
    print(f"\n=== Aggregate validation: {n} (x, K, |err|) triples ===")
    if n < 4:
        return

    sigmas = [t[0] for t in triples]
    errs = [t[1] for t in triples]
    # Pearson r:
    mean_s = sum(sigmas) / n
    mean_e = sum(errs) / n
    cov = sum((s - mean_s) * (e - mean_e) for s, e in zip(sigmas, errs)) / n
    var_s = sum((s - mean_s) ** 2 for s in sigmas) / n
    var_e = sum((e - mean_e) ** 2 for e in errs) / n
    r = cov / math.sqrt(var_s * var_e) if var_s * var_e > 0 else float("nan")
    # Best-fit through origin: slope = sum(sigma*err)/sum(sigma^2)
    slope = sum(s * e for s, e in zip(sigmas, errs)) / sum(s * s for s in sigmas)
    expected_slope_half_gauss = math.sqrt(2 / math.pi)  # ~0.7979
    print(f"Pearson r(sigma_pred, |err|)      = {r:.4f}")
    print(f"slope-through-origin              = {slope:.4f}")
    print(f"expected slope (half-Gaussian)    = {expected_slope_half_gauss:.4f}")
    print(f"observed/expected slope ratio     = "
          f"{slope/expected_slope_half_gauss:.4f}")
    # Hit-rate aggregate
    n_hits = sum(1 for s, e, _, _ in triples if e <= 0.5)
    pooled_hit = n_hits / n
    pooled_pred_hit = sum(hit_rate_predicted_for_sigma(t[0]) for t in triples) / n
    print(f"pooled empirical hit-rate         = {pooled_hit:.4f}")
    print(f"pooled predicted hit-rate         = {pooled_pred_hit:.4f}")


def hit_rate_predicted_for_sigma(sigma: float, threshold: float = 0.5) -> float:
    if sigma <= 0:
        return 1.0
    if math.isinf(sigma):
        return 0.0
    return math.erf(threshold / (math.sqrt(2) * sigma))


# ----------------------------------------------------------------------------
# Asymptotic K* characterisation
# ----------------------------------------------------------------------------


def k_star_for_hit_rate(x: float, target_hit: float = 0.99,
                        threshold: float = 0.5) -> int:
    """Find smallest K with predicted hit_rate >= target_hit.

    Uses sigma <= threshold / (sqrt(2) * erfinv(target_hit)).
    """
    # invert: erf(threshold / (sqrt(2) sigma)) >= target_hit
    # sigma <= threshold / (sqrt(2) erfinv(target_hit))
    z = inv_erf(target_hit)
    if z <= 0:
        return 0
    sigma_max = threshold / (math.sqrt(2) * z)

    # sigma^2 = x log^2(K) / (2 pi^2 K log^2(x)) <= sigma_max^2
    # K / log^2(K) >= x / (2 pi^2 log^2(x) sigma_max^2)
    if sigma_max <= 0:
        return -1
    target_ratio = x / (2 * math.pi ** 2 * math.log(x) ** 2 * sigma_max ** 2)

    def f(K: int) -> float:
        if K <= 1:
            return 0.0
        return K / (math.log(K) ** 2)

    # Find an upper bracket where f(hi) >= target_ratio.
    hi = 4
    while f(hi) < target_ratio:
        hi *= 2
        if hi > 10 ** 200:
            return -1
    lo = max(2, hi // 2)
    while lo < hi:
        mid = (lo + hi) // 2
        if f(mid) >= target_ratio:
            hi = mid
        else:
            lo = mid + 1
    return lo


def inv_erf(y: float) -> float:
    """Approximate inverse error function (Winitzki, decent precision)."""
    if y <= -1 or y >= 1:
        return float("inf") * (1 if y > 0 else -1)
    # Use scipy if available; otherwise Winitzki approximation.
    try:
        from scipy.special import erfinv  # type: ignore
        return float(erfinv(y))
    except ImportError:
        a = 0.147
        ln = math.log(1 - y * y)
        first = (2 / (math.pi * a) + ln / 2)
        return math.copysign(math.sqrt(math.sqrt(first ** 2 - ln / a) - first), y)


def k_star_table():
    print("\n=== Predicted K* threshold for hit-rate >= 99% (|err| <= 0.5) ===")
    print(f"{'x':>10s} {'K*_99':>14s} {'K*/x':>10s} {'sqrt(x)':>10s} {'log^2(x)':>10s}")
    for log10x in [3, 4, 5, 6, 7, 8, 10, 15, 20, 30, 50, 100]:
        x = 10.0 ** log10x
        K = k_star_for_hit_rate(x, target_hit=0.99)
        K_str = f"{K:.3g}" if K > 0 else "n/a"
        ratio = K / x if K > 0 else float("nan")
        print(f"10^{log10x:<7d} {K_str:>14s} {ratio:>10.3g} "
              f"{math.sqrt(x):>10.2g} {math.log(x) ** 2:>10.2g}")


def k_star_table_50pct():
    """Hit rate threshold at 50% (median rounding succeeds)."""
    print("\n=== Predicted K* for hit-rate >= 50% (|err| <= 0.5) ===")
    print(f"{'x':>10s} {'K*_50':>14s} {'K*/x':>10s} {'K*/sqrt(x)':>13s} "
          f"{'sqrt(x)*log^2':>14s}")
    for log10x in [3, 4, 5, 6, 7, 8, 10, 15, 20, 30, 50, 100]:
        x = 10.0 ** log10x
        K = k_star_for_hit_rate(x, target_hit=0.50)
        K_str = f"{K:.3g}" if K > 0 else "n/a"
        ratio = K / x if K > 0 else float("nan")
        rsx = K / math.sqrt(x) if K > 0 else float("nan")
        sl2 = math.sqrt(x) * math.log(x) ** 2
        print(f"10^{log10x:<7d} {K_str:>14s} {ratio:>10.3g} "
              f"{rsx:>13.3g} {sl2:>14.2g}")


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------


def main():
    here = Path(__file__).parent
    csv_1e5 = here / "connes_amortisation_fluctuation_1e5.csv"
    csv_1e6 = here / "connes_amortisation_fluctuation_1e6.csv"

    csv_1e7 = here / "connes_amortisation_fluctuation_1e7.csv"

    compare_csv(csv_1e5, "x ~ 1e5..3e5")
    compare_csv(csv_1e6, "x ~ 1e6..3e6")
    if csv_1e7.exists():
        compare_csv(csv_1e7, "x ~ 1e7..3e7")

    per_sample_comparison(csv_1e5, "x ~ 1e5..3e5", policy="log2x")
    per_sample_comparison(csv_1e6, "x ~ 1e6..3e6", policy="log2x")
    if csv_1e7.exists():
        per_sample_comparison(csv_1e7, "x ~ 1e7..3e7", policy="log2x")

    paths_labels = [(csv_1e5, "1e5"), (csv_1e6, "1e6")]
    if csv_1e7.exists():
        paths_labels.append((csv_1e7, "1e7"))
    aggregate_correlation(paths_labels)

    k_star_table_50pct()
    k_star_table()

    print("\n=== Cross-decade scaling check (K = log^2(x) policy) ===")
    print(f"{'center':>8s} {'pred_med':>10s} {'emp_med':>10s} {'ratio':>10s} "
          f"{'pred/emp':>10s}")
    for path, label in paths_labels:
        rows = list(load_fluctuation_csv(path))
        x_med = statistics.median([r["x"] for r in rows])
        K_med = statistics.median([r["K_log2x"] for r in rows])
        pred_med = median_predicted(x_med, int(K_med))
        emp_med = statistics.median([r["abs_err_at_K_log2x"] for r in rows])
        print(f"{label:>8s} {pred_med:>10.3f} {emp_med:>10.3f} "
              f"{emp_med/pred_med:>10.3f} {pred_med/emp_med:>10.3f}")


if __name__ == "__main__":
    main()
