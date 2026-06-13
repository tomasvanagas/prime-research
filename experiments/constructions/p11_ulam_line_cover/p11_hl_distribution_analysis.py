"""
Thread 11 / Slot 5 — distributional analysis of HL singular series
across Ulam slope-+/-1 diagonals.

The LP gap LP_p / LP_r = 0.776 is NOT explained by the AVG HL boost
(prime-weighted average HL singular series over the |c| <= C_max
range is ~2.3, line-uniform ~2.0; if these were the LP boost factor,
LP_p / LP_r would be ~0.5, not 0.776). This is because:
1. The avg over ALL slope-+/-1 diagonals (including length-1 corners)
   is exactly 1 by HL theory (mean HL singular series over admissible
   quadratics = 1).
2. The diagonals in |c| <= 250 we enumerate are the LONG central ones
   (length ~sqrt(N)), which are also the HL-richest ones (small
   discriminants).
3. LP gap depends on the DISTRIBUTION of HL across diagonals (variance),
   not on the average alone.
4. The LP also has axis lines (~31% of LP weight) where HL = 1.

This script:
1. Reads the per-line HL data from p11_hl_singular_series.py.
2. Computes line-LENGTH-weighted average HL (the "honest" avg).
3. Computes weighted variance and distribution moments.
4. Models the LP gap under simple distribution models and compares.

CLI:
  python p11_hl_distribution_analysis.py --csv hl_singular_series_results.csv
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from collections import Counter

import numpy as np


def load_data(path: str) -> list[dict]:
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append({
                "direction": r["direction"],
                "c": int(r["c"]),
                "n_on_line": int(r["n_on_line"]),
                "n_primes": int(r["n_primes"]),
                "density": float(r["density"]),
                "n_quadratics": int(r["n_quadratics"]),
                "hl_mean": float(r["hl_mean"]),
                "hl_max": float(r["hl_max"]),
            })
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", default="hl_singular_series_results.csv")
    args = parser.parse_args()

    path = os.path.join(os.path.dirname(__file__), args.csv)
    rows = load_data(path)
    print(f"# Loaded {len(rows)} slope-+/-1 lines")

    # Line-uniform statistics
    hl = np.array([r["hl_mean"] for r in rows])
    n_on = np.array([r["n_on_line"] for r in rows])
    n_pr = np.array([r["n_primes"] for r in rows])
    densities = np.array([r["density"] for r in rows])

    print(f"\n## Line-uniform statistics")
    print(f"   mean HL: {hl.mean():.4f}")
    print(f"   std HL:  {hl.std():.4f}")
    print(f"   var HL:  {hl.var():.4f}")
    print(f"   min HL:  {hl.min():.4f}")
    print(f"   max HL:  {hl.max():.4f}")
    print(f"   median:  {np.median(hl):.4f}")

    # Length-weighted statistics
    w_len = n_on / n_on.sum()
    hl_len_mean = (hl * w_len).sum()
    hl_len_var = ((hl - hl_len_mean) ** 2 * w_len).sum()
    print(f"\n## Length-weighted statistics (line length n_on as weight)")
    print(f"   mean HL: {hl_len_mean:.4f}")
    print(f"   std HL:  {math.sqrt(hl_len_var):.4f}")

    # Prime-weighted statistics
    w_pr = n_pr / n_pr.sum()
    hl_pr_mean = (hl * w_pr).sum()
    hl_pr_var = ((hl - hl_pr_mean) ** 2 * w_pr).sum()
    print(f"\n## Prime-weighted statistics (n_primes as weight)")
    print(f"   mean HL: {hl_pr_mean:.4f}")
    print(f"   std HL:  {math.sqrt(hl_pr_var):.4f}")

    # Distribution histogram
    print(f"\n## HL distribution histogram (line-uniform)")
    bins = [0, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0]
    hist, _ = np.histogram(hl, bins=bins)
    for i, count in enumerate(hist):
        print(f"   HL in [{bins[i]:.1f}, {bins[i+1]:.1f}): {count:4d} lines  ({100*count/len(hl):.1f}%)")

    # Top-30 most prime-rich
    sorted_rows = sorted(rows, key=lambda r: -r["n_primes"])
    print(f"\n## Top-15 prime-richest slope-+/-1 lines (HL-confirmed)")
    print(f"{'dir':>10} {'c':>5} {'pr':>5} {'dens':>6} {'HL_avg':>8} {'HL_max':>8}")
    for r in sorted_rows[:15]:
        print(f"{r['direction']:>10} {r['c']:>5} {r['n_primes']:>5} {r['density']:>6.3f} "
              f"{r['hl_mean']:>8.4f} {r['hl_max']:>8.4f}")

    # Now estimate the LP gap under simple model:
    # Suppose LP picks the densest diagonals first. Roughly, if it puts
    # weight w_d on each diagonal d with HL constant theta_d, then
    # primes covered on d = theta_d * sqrt(N) / log N. Each prime needs
    # total weight 1.
    #
    # Toy model: cover all primes via slope-+/-1 only. For each prime p,
    # set w_+(d_+(p)) + w_-(d_-(p)) = 1. The LP min is then:
    #   sum_d w_+(d) + sum_d w_-(d)
    # where (1) on each diagonal d, sum of (1/2 contribution from each
    # prime there) = #primes_on_d / 2 = theta_d * sqrt(N) / (2 log N).
    # So at LP optimum (under symmetric model), the weight per direction
    # = #primes / (avg primes per diagonal) / 2.
    #
    # For random: avg primes per diag = sqrt(N)/log N (uniform), so
    # LP_random_diag_only = pi(N) / (sqrt(N)/log N) / 2 * 2 directions
    #                     = 2 pi(N) log N / sqrt(N) / 2 ... messy.
    #
    # Simpler: empirically compare.

    log_N = math.log(100000)  # if data is from N=10^5
    sqrt_N = math.sqrt(100000)
    pi_N = 9592

    print(f"\n## LP gap interpretation")
    print(f"   N = 10^5, sqrt(N) = {sqrt_N:.2f}, pi(N) = {pi_N}, log(N) = {log_N:.4f}")
    print(f"   Empirical LP_p / LP_r at N=10^5: 0.7807")
    print(f"   Ratio target 1/c = 1/0.7807 = 1.2809")
    print(f"   Line-uniform avg HL on |c|<=250 sampled diagonals: {hl.mean():.4f}")
    print(f"   Length-weighted avg HL: {hl_len_mean:.4f}")
    print(f"   (Should equal 1 in true infinite limit; current sample biased high)")
    print(f"   Variance of HL (line-uniform): {hl.var():.4f}")
    print(f"   Std-dev of HL (line-uniform):  {hl.std():.4f}")

    # Coverage check: total prime-line incidences from our enumerated set
    # Each prime is on 1 slope-+1 diagonal and 1 slope-(-1) diagonal,
    # so total incidences = 2 * pi(N) = 19184. We found:
    incidences_found = n_pr.sum()
    print(f"\n   2 * pi(N) (expected incidences): {2 * pi_N}")
    print(f"   Incidences found in |c|<=250: {incidences_found}")
    print(f"   Coverage fraction: {incidences_found / (2 * pi_N):.4f}")
    print(f"   The {1 - incidences_found/(2*pi_N):.0%} of incidences in unenumerated")
    print(f"   short |c|>250 diagonals likely have HL ~ 1 (uniform corners),")
    print(f"   pulling the global avg from {hl.mean():.2f} down to ~1.")


if __name__ == "__main__":
    main()
