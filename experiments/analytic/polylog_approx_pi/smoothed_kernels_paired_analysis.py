"""Paired sign-test analysis of smoothed_kernels_data.csv.

For each (anchor, K_compute, kernel != hard), test whether
err_kernel(x_i) < err_hard(x_i) more often than chance under the
null hypothesis of no difference. Outputs:

  - ratio_paired = sqrt(mean(err_kernel^2) / mean(err_hard^2))   (== sigma ratio)
  - wins / N = #{i : |err_kernel(x_i)| < |err_hard(x_i)|} / N
  - sign_pval = P(Binom(N, 0.5) >= wins) (one-sided test for kernel-beats-hard)

Also prints the inverse "hard-beats-kernel" sign p-value for cells where ratio > 1.
"""

import csv
import math
import sys
from collections import defaultdict
from math import comb

HERE = "experiments/analytic/polylog_approx_pi"
data_csv = f"{HERE}/smoothed_kernels_data.csv"


def main():
    errs = defaultdict(dict)  # (anchor, x, K) -> kernel -> abs_err
    with open(data_csv) as f:
        reader = csv.DictReader(f)
        for r in reader:
            key = (int(r['anchor_log10']), int(r['x']), int(r['K_compute']))
            errs[key][r['kernel']] = float(r['abs_err'])

    agg = defaultdict(lambda: defaultdict(list))
    for key, kerr in errs.items():
        if 'hard' not in kerr:
            continue
        for kn, e in kerr.items():
            agg[(key[0], key[2])][kn].append((e, kerr['hard']))

    print(f"{'anchor':>6s} {'K':>5s} {'kernel':>9s}  "
          f"{'ratio':>7s} {'wins':>7s} {'p_kernel_beats':>15s} {'p_hard_beats':>13s}")
    print('-' * 70)
    for (anc, K) in sorted(agg.keys()):
        for kn in ['hard', 'tukey25', 'tukey50', 'cosine', 'riesz', 'triangle',
                   'riesz4', 'hamming', 'hann']:
            if kn not in agg[(anc, K)]:
                continue
            pairs = agg[(anc, K)][kn]
            N = len(pairs)
            sumK = sum(e ** 2 for e, _ in pairs)
            sumH = sum(h ** 2 for _, h in pairs)
            ratio = math.sqrt(sumK / sumH) if sumH > 0 else float('nan')
            wins = sum(1 for e, h in pairs if e < h)
            if kn == 'hard':
                p_k = float('nan')
                p_h = float('nan')
            else:
                # one-sided "kernel beats hard"
                p_k = sum(comb(N, k) for k in range(wins, N + 1)) / (2 ** N)
                # one-sided "hard beats kernel"
                p_h = sum(comb(N, k) for k in range(N - wins, N + 1)) / (2 ** N)
            wins_str = f"{wins:>2d}/{N}"
            print(f"{anc:>6d} {K:>5d} {kn:>9s}  "
                  f"{ratio:>7.3f} {wins_str:>7s} {p_k:>15.4f} {p_h:>13.4f}")
        print()

    # ----- Aggregate per kernel: average ratio and median paired p_hard_beats -----
    print('=' * 70)
    print("Per-kernel summary across all 12 (anchor, K) cells:")
    print(f"{'kernel':>9s}  {'mean(ratio)':>12s} {'#cells_p<0.05':>15s} "
          f"{'#cells_ratio>=1.05':>20s}")
    by_kernel = defaultdict(list)
    for (anc, K) in agg:
        for kn, pairs in agg[(anc, K)].items():
            if kn == 'hard':
                continue
            N = len(pairs)
            sumK = sum(e ** 2 for e, _ in pairs)
            sumH = sum(h ** 2 for _, h in pairs)
            ratio = math.sqrt(sumK / sumH) if sumH > 0 else float('nan')
            wins = sum(1 for e, h in pairs if e < h)
            p_h = sum(comb(N, k) for k in range(N - wins, N + 1)) / (2 ** N)
            by_kernel[kn].append((ratio, p_h))
    for kn in ['tukey25', 'tukey50', 'cosine', 'riesz', 'triangle',
               'riesz4', 'hamming', 'hann']:
        cells = by_kernel[kn]
        ratios = [r for r, _ in cells]
        ph = [p for _, p in cells]
        mean_r = sum(ratios) / len(ratios)
        n_sig = sum(1 for p in ph if p < 0.05)
        n_big = sum(1 for r in ratios if r >= 1.05)
        print(f"{kn:>9s}  {mean_r:>12.4f} {n_sig:>15d} {n_big:>20d}")


if __name__ == '__main__':
    main()
