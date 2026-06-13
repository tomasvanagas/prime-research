"""
Multi-sample empirical validation of the sigma-formula for pi(x) - pi_K(x).

Thread 7 (commit, slot 2, S241). Extends slot-1's single-anchor measurements
in `polylog_approx_pi.py` with multi-sample averaging at x near {10^7, 10^8,
10^9}. Tests the half-Gaussian hypothesis under random phases: if the
S195 sigma-formula

    sigma(K, x) ~ sqrt(x) * log K / (pi * sqrt(2 K) * log x)        (*)

is correct and {gamma_j log x mod 2pi} is iid uniform, then |err| =
|pi(x) - pi_K(x)| follows a half-normal with scale sigma. The empirical
distribution should have:
  - mean / sigma -> sqrt(2/pi) ~= 0.7979
  - median / sigma -> 0.6745
  - cdf(|err|/sigma <= z) = erf(z/sqrt(2))

The slot-2 contribution: tighten this empirically beyond S195's three-decade
multi-sample data (x ~ 1e5, 1e6, 1e7, 40 samples each) to x ~ 1e8, 1e9
(30 samples each) using a segmented sieve for pi(x_j) anchored at known
pi(10^k) values.

Multi-sample method:
  1. Anchor x_0 in {1e7, 1e8, 1e9} with known pi(x_0).
  2. Sample x_j = x_0 * (1 + j * delta), j = 0..N-1, delta chosen so the
     window spans many decoherence lengths (2*pi / gamma_K_max).
  3. Segmented Eratosthenes sieve over [x_0 + 1, x_max] gives prime
     indicators; cumulative sum + anchor gives pi(x_j) for all j.
  4. For each x_j compute pi_K(x_j) for K = 1..8000 using cached zeros.
  5. Record |err| at K-policy values {1, log x, log^2 x, log^3 x, K_max}.

CLI args:
  --anchors  comma list of log10 anchors  (default: 7,8,9)
  --N        samples per anchor           (default: 30)
  --K-max    hard zero cap                (default: 8000)
  --dps      mpmath precision             (default: 25)
  --M        Mobius truncation in R_at_rho (default: 8)
  --csv      output CSV path              (default: multi_sample_data.csv)
  --quiet    suppress per-K progress prints

Run time: ~5-15 minutes per anchor at N=30, K_max=8000 (mostly R_at_rho calls).
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
import time
from dataclasses import dataclass
from typing import Dict, List, Tuple

import mpmath as mp

# Reuse R(x), R_at_rho, sigma_predicted, zero loader from slot 1.
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from polylog_approx_pi import (  # type: ignore
    PI_KNOWN,
    get_zeros,
    riemann_R,
    R_at_rho,
    sigma_predicted,
)


# -----------------------------------------------------------------------------
# Segmented prime sieve (counts primes in (lo, hi]; expects lo >= 2).
# -----------------------------------------------------------------------------

def base_primes_up_to(N: int) -> List[int]:
    """Eratosthenes sieve for primes <= N."""
    if N < 2:
        return []
    sieve = bytearray(b"\x01" * (N + 1))
    sieve[0] = sieve[1] = 0
    for i in range(2, int(N ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = b"\x00" * (((N - i * i) // i) + 1)
    return [i for i in range(2, N + 1) if sieve[i]]


def segmented_prime_indicator(lo: int, hi: int) -> bytearray:
    """Return a bytearray seg of length hi - lo + 1 with seg[k] = 1 iff
    lo + k is prime. Requires lo >= 2."""
    assert lo >= 2 and hi >= lo
    sqrt_hi = int(hi ** 0.5) + 1
    bp = base_primes_up_to(sqrt_hi)
    n = hi - lo + 1
    seg = bytearray(b"\x01" * n)
    for p in bp:
        # First multiple of p in [lo, hi] that is >= p^2.
        start = max(p * p, ((lo + p - 1) // p) * p)
        if start > hi:
            continue
        a = start - lo
        count = (n - a + p - 1) // p
        seg[a::p] = b"\x00" * count
    return seg


def pi_array_for_samples(
    x_anchor: int, pi_anchor: int, sample_xs: List[int],
    chunk_size: int = 10 ** 8,
) -> List[int]:
    """For each sample x in sample_xs (assumed > x_anchor and sorted ascending),
    return pi(x). Method: chunked segmented sieve over (x_anchor, max(sample_xs)],
    cumulative prime counts via numpy, anchor-corrected at sample positions.

    chunk_size caps the bytearray length so this works for x_anchor up to 10^10
    or beyond without OOM. Default 10^8 ~ 100MB.
    """
    import numpy as np
    assert all(x > x_anchor for x in sample_xs)
    assert sample_xs == sorted(sample_xs)
    pi_xs: List[int] = []
    cum = 0
    j = 0
    lo = x_anchor + 1
    last_target = sample_xs[-1]
    while lo <= last_target and j < len(sample_xs):
        hi = min(lo + chunk_size - 1, last_target)
        seg = segmented_prime_indicator(lo, hi)
        seg_np = np.frombuffer(seg, dtype=np.uint8)
        cum_arr = np.cumsum(seg_np, dtype=np.int64)
        while j < len(sample_xs) and sample_xs[j] <= hi:
            k = sample_xs[j] - lo
            pi_xs.append(pi_anchor + cum + int(cum_arr[k]))
            j += 1
        cum += int(cum_arr[-1])
        lo = hi + 1
    assert len(pi_xs) == len(sample_xs)
    return pi_xs


# -----------------------------------------------------------------------------
# Multi-sample driver per anchor.
# -----------------------------------------------------------------------------

@dataclass
class SampleRow:
    anchor_log10: int
    x_anchor: int
    x: int
    pi_x: int
    R_x: float
    K: int
    pi_K: float
    err: float
    abs_err: float
    sigma_pred: float
    ratio_emp_over_pred: float
    policy_name: str


def policy_K_values(x: float, K_max: int) -> Dict[str, int]:
    lx = math.log(x)
    return {
        "K=1": 1,
        "logx": min(K_max, max(1, int(round(lx)))),
        "log2x": min(K_max, max(1, int(round(lx ** 2)))),
        "log3x": min(K_max, max(1, int(round(lx ** 3)))),
        "Kmax": K_max,
    }


def evaluate_anchor(
    anchor_log10: int,
    N_samples: int,
    K_max: int,
    gammas: List[mp.mpf],
    dps: int,
    M_mobius: int,
    quiet: bool,
) -> List[SampleRow]:
    x_anchor = 10 ** anchor_log10
    if x_anchor not in PI_KNOWN:
        raise SystemExit(f"pi(10^{anchor_log10}) not in PI_KNOWN.")
    pi_anchor = PI_KNOWN[x_anchor]

    # Sample at x_j = x_anchor * (1 + j * delta).
    # delta chosen so {gamma_j log x_j mod 2pi} are independent across j.
    # Decoherence length in log x: 2pi / gamma_{K_max}; gamma_{K_max} ~= 11000
    # at K_max = 8000.  Window log diff = N * delta should be >> N * 2pi / gamma.
    # Use window of half-decade in log -> delta_log = 0.5 * log(10) / N ~ 0.038.
    # Equivalently delta = exp(0.038 * j) - 1.
    log_lo = math.log(float(x_anchor))
    log_hi = log_lo + 0.5 * math.log(10.0)  # half decade
    sample_xs = sorted({int(round(math.exp(log_lo + j * (log_hi - log_lo) / (N_samples - 1))))
                        for j in range(N_samples)})
    # ensure all > x_anchor
    sample_xs = [x for x in sample_xs if x > x_anchor]
    # if dedup dropped some, top up by stepping linearly past hi.
    if len(sample_xs) < N_samples:
        extra_needed = N_samples - len(sample_xs)
        last = sample_xs[-1]
        gap = max(1, int(0.001 * x_anchor))
        for j in range(1, extra_needed + 1):
            sample_xs.append(last + j * gap)
    sample_xs = sorted(sample_xs[:N_samples])

    print(f"\n[anchor 10^{anchor_log10}]  N={len(sample_xs)} samples in "
          f"[{sample_xs[0]:.4g}, {sample_xs[-1]:.4g}]; pi(10^{anchor_log10})={pi_anchor}")

    # Compute pi(x) for each sample via segmented sieve.
    t0 = time.time()
    pi_xs = pi_array_for_samples(x_anchor, pi_anchor, sample_xs)
    print(f"  segmented-sieve pi(x_j): {time.time() - t0:.1f}s "
          f"(pi range [{pi_xs[0]}, {pi_xs[-1]}])")

    # Compute pi_K(x_j) for each sample.
    rows: List[SampleRow] = []
    for s, (x, pi_x) in enumerate(zip(sample_xs, pi_xs)):
        t0 = time.time()
        R_x = riemann_R(x, M=20, dps=dps)
        correction = mp.mpf(0)
        # Cumulative loop K = 1..K_max.
        # Record at policy points.
        policies = policy_K_values(x, K_max)
        K_targets = sorted(set(policies.values()))
        # K_target -> list of policy names hitting it.
        ktp = {k: [name for name, kv in policies.items() if kv == k] for k in K_targets}
        K_done = 0
        for K_target in K_targets:
            for k in range(K_done, K_target):
                Rrho = R_at_rho(x, gammas[k], M=M_mobius, dps=dps)
                correction += 2 * Rrho.real
            K_done = K_target

            pi_K = R_x - correction
            err = float(pi_K) - pi_x
            sig = sigma_predicted(x, K_target)
            ratio = abs(err) / sig if (sig > 0 and not math.isinf(sig)) else float("nan")
            for name in ktp[K_target]:
                rows.append(SampleRow(
                    anchor_log10=anchor_log10,
                    x_anchor=x_anchor,
                    x=x,
                    pi_x=pi_x,
                    R_x=float(R_x),
                    K=K_target,
                    pi_K=float(pi_K),
                    err=err,
                    abs_err=abs(err),
                    sigma_pred=sig,
                    ratio_emp_over_pred=ratio,
                    policy_name=name,
                ))
        if not quiet:
            print(f"  sample {s + 1}/{len(sample_xs)}: x={x}, "
                  f"|err@Kmax|={abs(rows[-1].err):.3f}, "
                  f"{time.time() - t0:.1f}s")

    return rows


# -----------------------------------------------------------------------------
# Aggregate statistics & half-normal fit.
# -----------------------------------------------------------------------------

def half_normal_cdf(z: float) -> float:
    """P(|Z| <= z) for Z ~ N(0,1) = erf(z/sqrt(2))."""
    if z < 0:
        return 0.0
    return math.erf(z / math.sqrt(2))


def ks_statistic(ratios: List[float]) -> float:
    """One-sample Kolmogorov-Smirnov statistic against the half-normal CDF."""
    n = len(ratios)
    if n == 0:
        return float("nan")
    s = sorted(ratios)
    D = 0.0
    for i, z in enumerate(s):
        F_emp_lo = i / n
        F_emp_hi = (i + 1) / n
        F_th = half_normal_cdf(z)
        D = max(D, abs(F_emp_lo - F_th), abs(F_emp_hi - F_th))
    return D


def ks_pvalue_approx(D: float, n: int) -> float:
    """Approximate two-sided KS p-value (Kolmogorov asymptotic)."""
    if n < 2:
        return float("nan")
    en = math.sqrt(n) * D + 0.12 * D + 0.11 * D / math.sqrt(n)
    # Sum of alternating series of exp(-2 j^2 en^2)
    s = 0.0
    for j in range(1, 50):
        term = 2 * (-1) ** (j - 1) * math.exp(-2 * j * j * en * en)
        s += term
        if abs(term) < 1e-12:
            break
    return min(1.0, max(0.0, s))


def aggregate_stats_per_policy(rows: List[SampleRow]):
    """Group by (anchor_log10, policy_name) and report distribution stats.

    Two scale assumptions for the KS fit are reported:
      - vs predicted: ratios = |err| / sigma_pred (tests S195's exact formula)
      - vs effective: ratios = |err| / sigma_eff, where sigma_eff is the
        rms of |err| (tests whether distribution shape is half-Gaussian
        regardless of scale; absorbs the GUE pair-correlation 0.74 factor).
    """
    groups: Dict[Tuple[int, str], List[SampleRow]] = {}
    for r in rows:
        groups.setdefault((r.anchor_log10, r.policy_name), []).append(r)

    print("\n" + "=" * 130)
    print(f"{'anchor':>7s} {'policy':>8s} {'K':>6s} {'N':>4s} "
          f"{'sig_pred':>10s} {'med|err|':>9s} {'mean|err|':>10s} "
          f"{'med/sig':>9s} {'mean/sig':>9s} "
          f"{'KSp_pred':>9s} {'KSp_eff':>9s} {'sig_eff/sig_pred':>17s}")
    print("-" * 130)

    summary = []
    for (anc, pol), grp in sorted(groups.items()):
        if pol == "K=1":
            continue  # K=1 has infinite sigma_pred; skip
        Ks = [r.K for r in grp]
        sigs = [r.sigma_pred for r in grp]
        errs = [r.abs_err for r in grp]
        ratios = [r.abs_err / r.sigma_pred for r in grp
                  if r.sigma_pred > 0 and not math.isinf(r.sigma_pred)]
        if not ratios:
            continue
        med = sorted(errs)[len(errs) // 2]
        mean = sum(errs) / len(errs)
        sig_med = sorted(sigs)[len(sigs) // 2]
        ratio_med = sorted(ratios)[len(ratios) // 2]
        ratio_mean = sum(ratios) / len(ratios)
        D_pred = ks_statistic(ratios)
        p_pred = ks_pvalue_approx(D_pred, len(ratios))
        # Scale-free KS: use sigma_eff = rms(|err|), divide and refit.
        sig_eff = math.sqrt(sum(e * e for e in errs) / len(errs))
        ratios_eff = [e / sig_eff for e in errs] if sig_eff > 0 else []
        D_eff = ks_statistic(ratios_eff)
        p_eff = ks_pvalue_approx(D_eff, len(ratios_eff))
        scale_ratio = sig_eff / sig_med if sig_med > 0 else float("nan")
        K_used = sorted(set(Ks))[len(set(Ks)) // 2]
        print(f"{'10^' + str(anc):>7s} {pol:>8s} {K_used:>6d} {len(grp):>4d} "
              f"{sig_med:>10.3f} {med:>9.3f} {mean:>10.3f} "
              f"{ratio_med:>9.3f} {ratio_mean:>9.3f} "
              f"{p_pred:>9.3f} {p_eff:>9.3f} {scale_ratio:>17.3f}")
        summary.append({
            "anchor_log10": anc,
            "policy": pol,
            "K": K_used,
            "N": len(grp),
            "sigma_pred_median": sig_med,
            "sigma_eff_rms": sig_eff,
            "median_abs_err": med,
            "mean_abs_err": mean,
            "ratio_median": ratio_med,
            "ratio_mean": ratio_mean,
            "ks_D_pred": D_pred,
            "ks_p_pred": p_pred,
            "ks_D_eff": D_eff,
            "ks_p_eff": p_eff,
            "sig_eff_over_sig_pred": scale_ratio,
        })
    return summary


# -----------------------------------------------------------------------------
# Theoretical extrapolation table.
# -----------------------------------------------------------------------------

def extrapolation_table(log10_xs: List[int], alphas: List[float]):
    """For each (x, alpha) print sigma_pred(K = log^alpha x) and sigma/sqrt(x).
    Use the named-exponent corollary."""
    print("\n" + "=" * 110)
    print("EXTRAPOLATION (sigma_pred for K = log^alpha x); rows = alpha, cols = log10 x")
    print("- raw sigma below diagonal, sigma / sqrt(x) above")
    print()
    header_x = " ".join(f"{'1e' + str(k):>11s}" for k in log10_xs)
    print(f"{'alpha':>5s}  {'K_at(1e10)':>11s}   {'-- sigma(x, K=log^alpha x) --':>20s}")
    print(f"{'':>5s}  {'':>11s}   {header_x}")
    for alpha in alphas:
        K_at_10 = max(1, int(round(math.log(1e10) ** alpha)))
        sigs = []
        for k in log10_xs:
            x = 10.0 ** k
            lx = math.log(x)
            K = max(1, int(round(lx ** alpha)))
            sig = sigma_predicted(x, K)
            sigs.append(sig)
        print(f"{alpha:>5.1f}  {K_at_10:>11d}   "
              + " ".join(f"{s:>11.3g}" for s in sigs))
    print()
    print(f"{'alpha':>5s}  {'':>11s}  {'-- sigma / sqrt(x) --':>20s}")
    print(f"{'':>5s}  {'':>11s}   {header_x}")
    for alpha in alphas:
        sigs_norm = []
        for k in log10_xs:
            x = 10.0 ** k
            lx = math.log(x)
            K = max(1, int(round(lx ** alpha)))
            sig = sigma_predicted(x, K)
            sigs_norm.append(sig / math.sqrt(x))
        print(f"{alpha:>5.1f}  {'':>11s}   "
              + " ".join(f"{s:>11.3e}" for s in sigs_norm))


# -----------------------------------------------------------------------------
# Driver.
# -----------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--anchors", default="7,8,9",
                   help="Comma list of log10 anchors.")
    p.add_argument("--N", type=int, default=30, help="Samples per anchor.")
    p.add_argument("--K-max", dest="K_max", type=int, default=8000,
                   help="Hard cap on K (upper bound for the cumulative R_at_rho loop).")
    p.add_argument("--dps", type=int, default=25)
    p.add_argument("--M", type=int, default=8)
    p.add_argument("--csv", default="multi_sample_data.csv")
    p.add_argument("--summary-csv", default="multi_sample_summary.csv")
    p.add_argument("--quiet", action="store_true")
    p.add_argument("--no-empirical", action="store_true",
                   help="Skip the per-anchor empirical run (theoretical extrapolation only).")
    return p.parse_args()


def main():
    args = parse_args()
    anchors = [int(s) for s in args.anchors.split(",")]

    rows: List[SampleRow] = []

    if not args.no_empirical:
        print(f"Loading {args.K_max} zeta zeros (mp dps={args.dps})...", flush=True)
        t0 = time.time()
        gammas = get_zeros(args.K_max, dps=args.dps)
        print(f"  loaded {len(gammas)} zeros in {time.time() - t0:.1f}s "
              f"(gamma_1={float(gammas[0]):.4f}, gamma_K={float(gammas[-1]):.4f})")

        for anc in anchors:
            sub = evaluate_anchor(anc, args.N, args.K_max, gammas,
                                  args.dps, args.M, args.quiet)
            rows.extend(sub)

        # ---------- summary ----------
        summary = aggregate_stats_per_policy(rows)

        # ---------- CSV outputs ----------
        out_rows = args.csv if os.path.isabs(args.csv) else os.path.join(HERE, args.csv)
        with open(out_rows, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow([
                "anchor_log10", "x_anchor", "x", "pi_x", "R_x", "policy", "K",
                "pi_K", "err", "abs_err", "sigma_pred", "ratio_emp_over_pred",
            ])
            for r in rows:
                w.writerow([
                    r.anchor_log10, r.x_anchor, r.x, r.pi_x,
                    f"{r.R_x:.6f}", r.policy_name, r.K,
                    f"{r.pi_K:.6f}", f"{r.err:.6f}", f"{r.abs_err:.6f}",
                    f"{r.sigma_pred:.6f}", f"{r.ratio_emp_over_pred:.6f}",
                ])
        print(f"\nWrote per-sample CSV: {out_rows}  ({len(rows)} rows)")

        out_sum = args.summary_csv if os.path.isabs(args.summary_csv) else os.path.join(HERE, args.summary_csv)
        with open(out_sum, "w", newline="") as fh:
            w = csv.writer(fh)
            keys = ["anchor_log10", "policy", "K", "N",
                    "sigma_pred_median", "sigma_eff_rms",
                    "median_abs_err", "mean_abs_err",
                    "ratio_median", "ratio_mean",
                    "ks_D_pred", "ks_p_pred", "ks_D_eff", "ks_p_eff",
                    "sig_eff_over_sig_pred"]
            w.writerow(keys)
            for s in summary:
                w.writerow([s[k] for k in keys])
        print(f"Wrote summary CSV: {out_sum}  ({len(summary)} rows)")

    # ---------- theoretical extrapolation ----------
    extrapolation_table(
        log10_xs=[6, 8, 10, 12, 15, 18, 21, 24],
        alphas=[2.0, 3.0, 4.0, 6.0, 8.0],
    )


if __name__ == "__main__":
    main()
