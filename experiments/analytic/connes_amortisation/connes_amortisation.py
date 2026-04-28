"""
Connes-Consani-Moscovici operator amortisation: per-query lower bound.

Thread 2 (commit S193, Connes amortisation).

S53 closed Chain A on a per-query argument: building the operator costs
at least O(K^3) for K zeros, exceeding any plausible saving over O(1)
queries. S190 asked: is the construction actually amortisable across
k = polylog(N) queries with shared spectral data? If so, per-query
cost = T_setup/k + T_query and the polylog regime opens IF
T_setup <= N^c * k for small c.

S193 adversarial finding: the answer is NO, but the structural reason
is sharper than S53's. Even with FREE amortised access to any zero,
the per-query cost is bounded below by the EXPLICIT FORMULA truncation
T(x), which is K_min(x) >= x^{1/2-eps} / log x by Riemann-von Mangoldt.

This script empirically confirms K_min(x) scales as x^{1/2} for the
worst-case error |pi_explicit(x; K) - pi(x)| <= 0.5 threshold.

The bottleneck is NOT operator construction. It is the bit-content of
the explicit-formula truncation. Thread 2's amortisation question
therefore reduces to Thread 3's (Galway frontier) "is K = polylog
sufficient in distribution?" question.

Falsifiable: if K_min(x) is observed to scale slower than x^{1/2-eps}
across the tested x values, the structural argument fails.
"""

from __future__ import annotations

import math
import time
from typing import List, Tuple

import mpmath as mp
import numpy as np

# -----------------------------------------------------------------------------
# Riemann zeros and explicit formula utilities.
# -----------------------------------------------------------------------------


def get_zeros(K: int, dps: int = 30) -> List[float]:
    """Load the imaginary parts gamma_1..gamma_K from the cached zeros files.
    Falls back to mpmath.zetazero for uncached counts.
    """
    import os as _os
    mp.mp.dps = dps
    here = _os.path.dirname(_os.path.abspath(__file__))
    data_dir = _os.path.normpath(_os.path.join(here, "..", "..", "..", "data"))
    for size in (300, 500, 1000, 2000, 8000):
        if size >= K:
            fname = _os.path.join(data_dir, f"zeta_zeros_{size}.txt")
            if _os.path.exists(fname):
                with open(fname) as f:
                    zs = [line.strip() for line in f if line.strip()]
                if len(zs) >= K:
                    return [float(z) for z in zs[:K]]
    return [float(mp.zetazero(k).imag) for k in range(1, K + 1)]


def riemann_R(x: float, M: int = 20) -> float:
    """Riemann's R(x) = sum_{n>=1} mu(n)/n * li(x^(1/n)).
    Truncated at n <= M. Uses mpmath for accuracy.
    """
    mp.mp.dps = 30
    total = mp.mpf(0)
    # mobius
    mobius = mobius_sieve(M)
    for n in range(1, M + 1):
        mu = mobius[n]
        if mu == 0:
            continue
        # li(x^(1/n))
        try:
            val = mp.li(mp.power(x, mp.mpf(1) / n))
            total += mp.mpf(mu) / n * val
        except Exception:
            continue
    return float(total)


def mobius_sieve(N: int) -> List[int]:
    mu = [1] * (N + 1)
    is_prime = [True] * (N + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, N + 1):
        if is_prime[i]:
            for j in range(i, N + 1, i):
                if j != i:
                    is_prime[j] = False
                mu[j] *= -1
            sq = i * i
            for j in range(sq, N + 1, sq):
                mu[j] = 0
    mu[0] = 0
    return mu


def R_at_rho(x: float, rho_imag: float, M: int = 12) -> complex:
    """R(x^rho) for rho = 1/2 + i*gamma_imag, computed via direct branch:
        R(x^rho) = sum_n mu(n)/n * Ei(rho * ln(x) / n).
    Avoids the principal-log phase-wrap bug of computing z = x^rho first.
    """
    mp.mp.dps = 30
    if x <= 1:
        return complex(0)
    ln_x = mp.log(mp.mpf(x))
    rho = mp.mpc(mp.mpf("0.5"), mp.mpf(rho_imag))
    mobius = mobius_sieve(M)
    total = mp.mpc(0)
    for n in range(1, M + 1):
        mu = mobius[n]
        if mu == 0:
            continue
        arg = rho * ln_x / n
        if abs(arg) < 1e-20:
            break
        total += mp.mpc(mu) / n * mp.ei(arg)
    return complex(total)


def pi_explicit(x: float, K: int, gammas: List[float]) -> float:
    """Compute pi(x) via explicit formula with first K nontrivial zeros.
        pi(x) ~ R(x) - 2 Re sum_{k<=K} R(x^rho_k).
    """
    if K > len(gammas):
        K = len(gammas)
    R_x = riemann_R(x)
    correction = 0.0
    for k in range(K):
        Rrho = R_at_rho(x, gammas[k])
        correction += 2 * Rrho.real
    return R_x - correction


def pi_brute(x: float) -> int:
    """Compute pi(x) by sieving (small x only)."""
    n = int(x)
    if n < 2:
        return 0
    sieve = bytearray(b"\x01" * (n + 1))
    sieve[0] = 0
    sieve[1] = 0
    for i in range(2, int(math.isqrt(n)) + 1):
        if sieve[i]:
            sieve[i * i::i] = b"\x00" * (((n - i * i) // i) + 1)
    return sum(sieve)


_PI_CACHE: dict = {}


def build_prime_count_array(N: int) -> np.ndarray:
    """Return cumulative prime-count array c such that c[k] = pi(k) for k <= N."""
    sieve = bytearray(b"\x01" * (N + 1))
    sieve[0] = 0
    sieve[1] = 0
    for i in range(2, int(math.isqrt(N)) + 1):
        if sieve[i]:
            sieve[i * i::i] = b"\x00" * (((N - i * i) // i) + 1)
    arr = np.frombuffer(bytes(sieve), dtype=np.uint8)
    return np.cumsum(arr, dtype=np.int64)


def pi_from_array(cum: np.ndarray, x: float) -> int:
    n = int(x)
    if n < 0:
        return 0
    if n >= len(cum):
        raise ValueError(f"pi_from_array: x={x} exceeds precomputed limit {len(cum)-1}")
    return int(cum[n])


def find_k_min(x: float, gammas: List[float], threshold: float = 0.5,
               k_max: int = None) -> Tuple[int, List[float]]:
    """Find the minimal K such that |pi_explicit(x; K) - pi(x)| <= threshold.

    Computes the running sum incrementally so total cost is O(K_max), not
    O(K_max^2).
    """
    if k_max is None:
        k_max = len(gammas)
    pi_x = pi_brute(x)
    R_x = riemann_R(x)
    correction = 0.0
    errors: List[float] = []
    K_min = -1
    for K in range(1, k_max + 1):
        Rrho = R_at_rho(x, gammas[K - 1])
        correction += 2 * Rrho.real
        approx = R_x - correction
        err = approx - pi_x
        errors.append(err)
        if K_min == -1 and abs(err) <= threshold:
            K_min = K
    return K_min, errors


# -----------------------------------------------------------------------------
# K_sustained measurement and reusable runner.
# -----------------------------------------------------------------------------


def measure_K_sustained_for_xs(xs, K_max, gammas, pi_lookup, tag=""):
    """Compute K_sustained(x) and K_first_eps05(x) for each x in xs.

    pi_lookup(x) -> int (true value of pi(x)).
    Returns list of dicts.
    """
    rows = []
    print(f"\n[{tag}] x        pi(x)    Ksust   K_eps05  sqrt(x)   K/sqrt(x)  time")
    print(f"[{tag}] " + "-" * 72)
    for x in xs:
        t1 = time.time()
        pi_x = pi_lookup(x)
        R_x = riemann_R(x)
        correction = 0.0
        rounded_eq: List[bool] = []
        first_eps05 = -1
        for K in range(1, K_max + 1):
            Rrho = R_at_rho(x, gammas[K - 1])
            correction += 2 * Rrho.real
            approx = R_x - correction
            err = approx - pi_x
            rounded_eq.append(round(approx) == pi_x)
            if first_eps05 == -1 and abs(err) <= 0.5:
                first_eps05 = K
        K_sust = -1
        for K in range(1, K_max + 1):
            if all(rounded_eq[K - 1:]):
                K_sust = K
                break
        sqrt_x = math.sqrt(x)
        ratio = K_sust / sqrt_x if K_sust > 0 else float('nan')
        elapsed = time.time() - t1
        print(f"[{tag}] %-9g %7d %7d %8d %8.1f %10.4f  (%.1fs)"
              % (x, pi_x, K_sust, first_eps05, sqrt_x, ratio, elapsed))
        rows.append(dict(x=x, pi_x=pi_x, K_sustained=K_sust,
                         K_first_eps05=first_eps05, sqrt_x=sqrt_x,
                         Ksust_over_sqrt_x=ratio, time=elapsed))
    return rows


# -----------------------------------------------------------------------------
# Main experiment.
# -----------------------------------------------------------------------------


def run_legacy_baseline():
    """The original S193 5-point baseline. Kept for regression sanity."""
    xs = [100.0, 1000.0, 5000.0, 10000.0, 50000.0]
    K_max = 500

    print("Computing first %d Riemann zeros..." % K_max)
    t0 = time.time()
    gammas = get_zeros(K_max, dps=30)
    print("  done in %.1fs (gamma_1 = %.4f, gamma_K = %.4f)"
          % (time.time() - t0, gammas[0], gammas[-1]))

    rows = []
    print("\nSUSTAINED K_min: smallest K such that round(pi_K') = pi(x) for")
    print("all K' in [K, K_max]. This is the rigour measure used by S53.")
    print("\nx        pi(x)  Ksust(x)  K_min0.5  sqrt(x)  K/sqrt(x)  time")
    print("-" * 70)
    for x in xs:
        t1 = time.time()
        pi_x = pi_brute(x)
        R_x = riemann_R(x)
        correction = 0.0
        errors: List[float] = []
        rounded_eq: List[bool] = []
        K_min_05 = -1
        for K in range(1, K_max + 1):
            Rrho = R_at_rho(x, gammas[K - 1])
            correction += 2 * Rrho.real
            approx = R_x - correction
            err = approx - pi_x
            errors.append(err)
            rounded_eq.append(round(approx) == pi_x)
            if K_min_05 == -1 and abs(err) <= 0.5:
                K_min_05 = K
        # sustained K_min: smallest K such that rounded_eq[K-1:] is all True.
        K_sust = -1
        for K in range(1, K_max + 1):
            if all(rounded_eq[K - 1:]):
                K_sust = K
                break
        sqrt_x = math.sqrt(x)
        ratio = K_sust / sqrt_x if K_sust > 0 else float('nan')
        print("%-8g %5d  %8d  %8d  %7.1f  %8.3f  (%.1fs)"
              % (x, pi_x, K_sust, K_min_05, sqrt_x, ratio, time.time() - t1))
        rows.append((x, pi_x, K_sust, K_min_05, sqrt_x, ratio,
                     errors[:50]))

    # Save CSV for plotting.
    import csv
    with open(
        "experiments/analytic/connes_amortisation/connes_amortisation_data.csv",
        "w",
    ) as f:
        w = csv.writer(f)
        w.writerow(["x", "pi_x", "K_sustained", "K_first_eps05", "sqrt_x",
                    "Ksust_over_sqrt_x"])
        for row in rows:
            w.writerow(list(row[:6]))

    # Linear regression in log-log: log K_sust ~ alpha * log(x) + beta.
    xs_log = np.log([r[0] for r in rows if r[2] > 0])
    Ks_log = np.log([r[2] for r in rows if r[2] > 0])
    if len(xs_log) >= 3:
        slope, intercept = np.polyfit(xs_log, Ks_log, 1)
        print("\nLog-log fit: K_sustained(x) ~ %.3f * x^%.3f"
              % (math.exp(intercept), slope))
        print("Predicted slope (Riemann-von Mangoldt worst case): 0.5")
        print("Empirical slope (limited range): %.3f" % slope)

    # Per-query cost statement.
    print("\n" + "=" * 70)
    print("AMORTISATION CONCLUSION (Thread 2 / S193):")
    print("Per-query cost = K_min(x; eps) zeros summed via explicit formula.")
    print("Even with FREE amortised access to all zeros, per-query is")
    print("Omega(K_min(x)) = Omega(sqrt(x)) up to log factors.")

    # ----------------------------------------------------------------
    # Section B — Connes vs Galway setup-cost comparison.
    # The adversarial finding: under FULL amortisation, S53's argument 2
    # (O(K^3) per-query) collapses, but the Connes setup cost is still
    # O(K^3) one-shot. Galway 2004 (Hiary t^(4/13) per zero) gives
    # O(K^(17/13)) one-shot for the SAME K zeros. Connes is dominated.
    print("\n" + "=" * 70)
    print("Setup-cost comparison: Connes vs Galway 2004 (Hiary)")
    print("=" * 70)
    print("Under amortisation across Q queries up to x = X:")
    print("  Connes (rank-one + dense eigensolver):")
    print("    setup     = O(K^3)        (N x N diag, N >= K)")
    print("    per-query = O(K)")
    print("  Galway 2004 + Hiary:")
    print("    setup     = O(K^(17/13))  (K zeros at O(t^(4/13)) each, t ~ K)")
    print("    per-query = O(K)")
    print("  Setup ratio = K^3 / K^(17/13) = K^(22/13) ~ K^1.692.")
    print()
    print("  %-12s %10s %14s %14s %12s" %
          ("X = 10^a", "K~sqrt(X)", "Connes K^3", "Galway K^17/13", "ratio"))
    print("  " + "-" * 65)
    for a in [4, 6, 10, 50, 100]:
        K = 10 ** (a / 2)
        c = K ** 3
        g = K ** (17.0 / 13.0)
        r = c / g
        print("  %-12s %10s %14s %14s %12s" % (
            f"10^{a}",
            f"10^{a/2:.1f}",
            f"10^{math.log10(c):.2f}",
            f"10^{math.log10(g):.2f}",
            f"10^{math.log10(r):.2f}",
        ))
    print()
    print("Conclusion: Connes setup is dominated by Galway/Hiary by")
    print("  10^(22a/26) ~ 10^(0.85*a) at x = 10^a (e.g., 10^85 worse at x=10^100).")
    print("  Even with infinite amortisation, Connes cannot beat Galway.")
    print("  The CCM operator's amortisation question therefore reduces to")
    print("  Thread 3 (Galway frontier): does K(x) = polylog(x) suffice for")
    print("  pi(x) +/- 1 in distribution? If yes, BOTH routes go polylog;")
    print("  if no, NEITHER does. Connes adds nothing distinct to amortise.")


# -----------------------------------------------------------------------------
# Slot 2 (S194) extension: wide-range and fluctuation modes.
# -----------------------------------------------------------------------------


def run_wide_range(K_max_cap: int, sieve_limit: int, output_csv: str):
    """Decade-spaced K_sustained(x) measurement on x = 1e3 .. 1e7.

    Uses the 8000-zero cache. K_max_cap caps the search depth.
    """
    xs = [1e3, 1e4, 5e4, 1e5, 5e5, 1e6, 5e6, 1e7]
    xs = [x for x in xs if x <= sieve_limit]

    print(f"WIDE-RANGE: xs = {xs}, K_max = {K_max_cap}")
    print(f"Loading {K_max_cap} cached zeros...")
    t0 = time.time()
    gammas = get_zeros(K_max_cap, dps=30)
    print(f"  loaded in {time.time() - t0:.1f}s "
          f"(gamma_1={gammas[0]:.4f}, gamma_K={gammas[-1]:.4f})")

    print(f"Building cumulative pi(.) array up to N={sieve_limit}...")
    t0 = time.time()
    cum = build_prime_count_array(sieve_limit)
    print(f"  done in {time.time() - t0:.1f}s. pi({sieve_limit})={int(cum[-1])}.")
    pi_lookup = lambda x: pi_from_array(cum, x)

    rows = measure_K_sustained_for_xs(xs, K_max_cap, gammas, pi_lookup,
                                       tag="wide")

    import csv
    with open(output_csv, "w") as f:
        w = csv.writer(f)
        w.writerow(["x", "pi_x", "K_sustained", "K_first_eps05",
                    "sqrt_x", "Ksust_over_sqrt_x", "time_seconds"])
        for r in rows:
            w.writerow([r["x"], r["pi_x"], r["K_sustained"],
                        r["K_first_eps05"], r["sqrt_x"],
                        r["Ksust_over_sqrt_x"], r["time"]])
    print(f"Wrote {output_csv}")

    # Log-log fit on rows where K_sust was found.
    valid = [r for r in rows if r["K_sustained"] > 0]
    if len(valid) >= 3:
        xs_log = np.log([r["x"] for r in valid])
        Ks_log = np.log([r["K_sustained"] for r in valid])
        slope, intercept = np.polyfit(xs_log, Ks_log, 1)
        print(f"\nLOG-LOG FIT (wide range, {len(valid)} points):")
        print(f"  K_sustained(x) ~ {math.exp(intercept):.4f} * x^{slope:.4f}")
        print(f"  Riemann-von Mangoldt baseline slope: 0.5")
    return rows


def run_fluctuation(center: float, n_samples: int, K_max_cap: int,
                    output_csv: str):
    """Dense sampling of x in a small log-band around `center`.

    For each sampled x we record the full error trajectory
        e_K(x) = R_K(x) - pi(x), for K = 1 .. K_max_cap
    where R_K(x) = R(x) - 2 sum_{k=1}^{K} Re R(x^rho_k).

    The distributional question (Thread 3 / Galway frontier):
    *fix K and ask what fraction of x in this band have |e_K(x)| <= 0.5.*
    If polylog K suffices in distribution, that fraction tends to 1
    as K grows polylogarithmically.

    Output CSV has one row per (x, K) pair for a small set of
    diagnostic K values (K = round(c * (log x)^p) for several (c, p)
    plus a sqrt(x) baseline). A second CSV stores the K_sustained
    summary per x.
    """
    lo = center
    hi = center * math.sqrt(10.0)  # half a decade
    xs = list(np.geomspace(lo, hi, n_samples))
    sieve_limit = int(hi) + 100

    print(f"FLUCTUATION: center={center}, n={n_samples}, "
          f"x in [{lo:.0f}, {hi:.0f}], K_max={K_max_cap}")
    print(f"Loading {K_max_cap} zeros...")
    t0 = time.time()
    gammas = get_zeros(K_max_cap, dps=30)
    print(f"  loaded in {time.time() - t0:.1f}s")

    print(f"Sieve up to {sieve_limit}...")
    t0 = time.time()
    cum = build_prime_count_array(sieve_limit)
    print(f"  done in {time.time() - t0:.1f}s")

    # diagnostic K-policies as functions of x.
    K_policies = {
        "log2x":   lambda x: int(round((math.log(x)) ** 2)),
        "log3x":   lambda x: int(round((math.log(x)) ** 3)),
        "5log2x":  lambda x: int(round(5.0 * (math.log(x)) ** 2)),
        "sqrtx":   lambda x: int(round(math.sqrt(x))),
        "halfsqrtx": lambda x: int(round(0.5 * math.sqrt(x))),
    }

    # Collect, for each x, the full e_K trajectory (only |error| values).
    per_x_records = []
    print(f"\n[fluc] Sweeping {len(xs)} x values, K={K_max_cap} each:")
    for i, x in enumerate(xs):
        t1 = time.time()
        pi_x = int(cum[int(x)])
        R_x = riemann_R(x)
        correction = 0.0
        errs = np.empty(K_max_cap)
        rounded_eq = np.empty(K_max_cap, dtype=bool)
        for K in range(1, K_max_cap + 1):
            Rrho = R_at_rho(x, gammas[K - 1])
            correction += 2 * Rrho.real
            approx = R_x - correction
            errs[K - 1] = approx - pi_x
            rounded_eq[K - 1] = round(approx) == pi_x
        # K_sustained: smallest K such that rounded_eq[K-1:] is all True.
        # = (last false index) + 2; if none false, K_sust = 1.
        false_idx = np.where(~rounded_eq)[0]
        if len(false_idx) == 0:
            K_sust = 1
        elif false_idx[-1] == K_max_cap - 1:
            K_sust = -1  # last point still false; not stabilised
        else:
            K_sust = int(false_idx[-1]) + 2
        per_x_records.append({
            "x": x, "pi_x": pi_x, "errs": errs,
            "K_sustained": K_sust, "rounded_eq": rounded_eq,
            "time": time.time() - t1,
        })
        if i % max(1, len(xs) // 10) == 0:
            print(f"  [{i+1}/{len(xs)}] x={x:.0f}, pi(x)={pi_x}, "
                  f"K_sust={K_sust}, t={time.time() - t1:.1f}s")

    # Summary CSV: one row per x.
    import csv
    with open(output_csv, "w") as f:
        w = csv.writer(f)
        cols = ["x", "pi_x", "K_sustained", "sqrt_x"]
        for name, fn in K_policies.items():
            cols += [f"K_{name}", f"abs_err_at_K_{name}",
                     f"hit_at_K_{name}"]
        w.writerow(cols)
        for rec in per_x_records:
            x = rec["x"]
            row = [x, rec["pi_x"], rec["K_sustained"], math.sqrt(x)]
            for name, fn in K_policies.items():
                K = max(1, min(fn(x), K_max_cap))
                err_K = rec["errs"][K - 1]
                hit = abs(err_K) <= 0.5
                row += [K, abs(err_K), int(hit)]
            w.writerow(row)
    print(f"Wrote {output_csv}")

    # Distributional summary.
    print(f"\nFLUCTUATION SUMMARY ({len(per_x_records)} x in [{lo:.0f}, "
          f"{hi:.0f}]):")
    valid = [r for r in per_x_records if r["K_sustained"] > 0]
    if valid:
        ratios = np.array([r["K_sustained"] / math.sqrt(r["x"])
                           for r in valid])
        print(f"  K_sust/sqrt(x):  min={ratios.min():.3f}, "
              f"med={np.median(ratios):.3f}, "
              f"max={ratios.max():.3f}, "
              f"std={ratios.std():.3f}")
    else:
        print("  K_sustained: no x stabilised within K_max")

    # Per-policy hit rates.
    print(f"\n  Policy           K(x_med)   median|err|   hit-rate (|err|<=0.5)")
    print(f"  " + "-" * 60)
    x_med = float(np.median([r["x"] for r in per_x_records]))
    for name, fn in K_policies.items():
        K_med = fn(x_med)
        hits = 0
        errs_at_K = []
        for rec in per_x_records:
            K = max(1, min(fn(rec["x"]), K_max_cap))
            err_K = abs(rec["errs"][K - 1])
            errs_at_K.append(err_K)
            if err_K <= 0.5:
                hits += 1
        med_err = np.median(errs_at_K)
        print(f"  {name:14s}  {K_med:8d}  {med_err:11.3f}   "
              f"{hits}/{len(per_x_records)} = "
              f"{hits / len(per_x_records):.2f}")

    return per_x_records


def main():
    import argparse
    import os as _os
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode",
                        choices=["legacy", "wide", "fluctuation", "both"],
                        default="both")
    parser.add_argument("--K-max", type=int, default=8000,
                        help="cap on zeros used (max 8000 from cache)")
    parser.add_argument("--sieve-limit", type=int, default=int(1.05e7),
                        help="upper sieve limit for wide-range mode")
    parser.add_argument("--fluc-center", type=float, default=1e5)
    parser.add_argument("--fluc-n", type=int, default=40,
                        help="number of fluctuation samples")
    parser.add_argument("--fluc-K-max", type=int, default=600,
                        help="K_max for fluctuation mode (each x is small)")
    args = parser.parse_args()

    here = _os.path.dirname(_os.path.abspath(__file__))

    if args.mode in ("legacy",):
        run_legacy_baseline()
        return

    if args.mode in ("wide", "both"):
        out = _os.path.join(here, "connes_amortisation_wide.csv")
        run_wide_range(args.K_max, args.sieve_limit, out)

    if args.mode in ("fluctuation", "both"):
        # Encode center in output name so re-runs at different centers
        # don't clobber each other.
        c = args.fluc_center
        if c >= 1.0:
            tag = f"{int(round(math.log10(c)))}"
        else:
            tag = "frac"
        out = _os.path.join(here, f"connes_amortisation_fluctuation_1e{tag}.csv")
        run_fluctuation(args.fluc_center, args.fluc_n,
                        args.fluc_K_max, out)


if __name__ == "__main__":
    main()
