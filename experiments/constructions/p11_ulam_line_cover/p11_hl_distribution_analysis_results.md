# p11_hl_distribution_analysis — distributional moments of HL singular series

Companion analysis to `p11_hl_singular_series.py`. Reads the per-line
HL CSV and computes line-uniform / length-weighted / prime-weighted
moments of the HL distribution across slope-±1 Ulam diagonals.

## Output (N=10⁵, |c|≤250, 377 lines)

```
## Line-uniform statistics
   mean HL: 2.0220
   std HL:  0.7419
   min HL:  ~0.04
   max HL:  4.7395
   median:  2.0125

## Length-weighted statistics (line length n_on as weight)
   mean HL: 2.0318
   std HL:  0.7394

## Prime-weighted statistics (n_primes as weight)
   mean HL: 2.2975
   std HL:  0.7315

## HL distribution histogram (line-uniform)
   HL in [0.0, 1.0):   34 lines  ( 9.0%)
   HL in [1.0, 1.5):   63 lines  (16.7%)
   HL in [1.5, 2.0):   89 lines  (23.6%)
   HL in [2.0, 2.5):   97 lines  (25.7%)
   HL in [2.5, 3.0):   66 lines  (17.5%)
   HL in [3.0, 4.0):   23 lines  ( 6.1%)
   HL in [4.0, 5.0):    5 lines  ( 1.3%)
   HL in [5.0, 7.0):    0 lines  ( 0.0%)

## Coverage
   2 * pi(N) (expected incidences): 19184
   Incidences found in |c|<=250:    13908
   Coverage fraction:               0.7250
```

## Interpretation

The slope-±1 diagonals at |c|≤250 cover 72.5% of slope-±1 prime-line
incidences (each prime is on 2 slope-±1 lines, one per direction).
Remaining 27.5% of incidences are on shorter |c|>250 corner
diagonals, which presumably have HL closer to 1. The sample is
biased toward long, central, HL-rich diagonals.

## Why the naive prediction LP_p / LP_r = 1/⟨θ⟩ ≈ 0.495 fails

See `p11_hl_singular_series_results.md` § "Structural identification
of the LP gap" for the empirical decomposition

```
c = 1 / (w_axis + w_diag * ⟨θ⟩_LP) = 1 / (0.31 + 0.69 * 1.42) = 0.776
```

where ⟨θ⟩_LP ≈ 1.42 is much smaller than the line-uniform avg HL on
the |c|≤250 sample (2.02), because (a) LP weight spreads across many
diagonals not just the densest, and (b) the LP must distribute weight
between slope-+1 and slope-(-1) directions (cannot put unit weight
on every diagonal).

## Falsifiers

The 1.42 figure is INFERRED from one decomposition; the LP-active
direction split (0.31 axis / 0.69 diagonals) is from S486 N=10⁴
only. Confirming asymptotic stability of these LP-equilibrium
quantities at N=10⁶ would require dual-inspection at scale, not
performed in this slot.
