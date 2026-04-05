========================================================================
Kt(delta(n)|n) ESTIMATION
========================================================================

--- (a) Arithmetic program search (n=1..500) ---
  BFS search done in 0.0s. Found programs for 500/500 values.
  n_range -> avg Kt bits (BFS)
    n=[1,10]: mean=4.64 bits
    n=[10,50]: mean=4.64 bits
    n=[50,100]: mean=4.64 bits
    n=[100,200]: mean=4.64 bits
    n=[200,500]: mean=4.86 bits

--- (b) Description complexity via expression fitting (n=1..10000) ---
  Model: constant -> residual entropy = 7.47 bits/value
  Model: linear(n) -> residual entropy = 7.48 bits/value
  Model: poly(5) -> residual entropy = 7.43 bits/value
  Model: sqrt(n)*sin(b*log(n)+c) -> residual entropy = 7.48 bits/value (4 params)

  Sliding-window entropy vs log(n): Kt ~ 0.3456 * log(n) + 3.2909
    R^2 residual sum: 11.8046
  Power-law fit: Kt ~ n^0.0603

--- (c) Conditional compression (block-level) ---
  k=   100: raw=9.84 b/v, diff=10.67 b/v, resid=9.44 b/v
  k=   500: raw=7.44 b/v, diff=6.49 b/v, resid=7.41 b/v
  k=  1000: raw=6.76 b/v, diff=5.96 b/v, resid=6.78 b/v
  k=  5000: raw=6.04 b/v, diff=5.44 b/v, resid=6.14 b/v
  k= 10000: raw=5.89 b/v, diff=5.35 b/v, resid=5.87 b/v
  k= 50000: raw=5.86 b/v, diff=5.43 b/v, resid=5.85 b/v
  k=100000: raw=5.91 b/v, diff=5.44 b/v, resid=5.92 b/v

  Conditional test: compress delta[i] vs compress (i, delta[i]):
  k=1000: delta_only=6.76 b/v, pairs=21.67 b/v, idx_only=6.86 b/v, saving=-8.048 b/v
  k=10000: delta_only=5.89 b/v, pairs=22.46 b/v, idx_only=2.66 b/v, saving=-13.904 b/v
  k=100000: delta_only=5.91 b/v, pairs=25.67 b/v, idx_only=3.08 b/v, saving=-16.687 b/v

--- (d) Incremental complexity (prediction residual entropy) ---

  Window W=10:
    n~500: H(resid) = 5.18 bits/value
    n~2000: H(resid) = 5.51 bits/value
    n~10000: H(resid) = 5.91 bits/value
    n~50000: H(resid) = 6.16 bits/value
    n~100000: H(resid) = 6.25 bits/value

  Window W=50:
    n~500: H(resid) = 5.66 bits/value
    n~2000: H(resid) = 6.34 bits/value
    n~10000: H(resid) = 6.73 bits/value
    n~50000: H(resid) = 7.05 bits/value
    n~100000: H(resid) = 7.14 bits/value

  Window W=200:
    n~2000: H(resid) = 6.86 bits/value
    n~10000: H(resid) = 7.15 bits/value
    n~50000: H(resid) = 7.59 bits/value
    n~100000: H(resid) = 7.65 bits/value

  Growth rate analysis of incremental entropy (W=50):
  Fit H ~ 0.2196 * log(n) + 4.6485
  Fit H ~ n^0.032787  (power law exponent)

  Ratio H(n) / log(n) at selected n:
    n=2000: H=6.34, H/log2(n)=0.578, H/sqrt(n)=0.14169, H/n^(1/3)=0.5029
    n=10000: H=6.77, H/log2(n)=0.509, H/sqrt(n)=0.06768, H/n^(1/3)=0.3141
    n=50000: H=7.04, H/log2(n)=0.451, H/sqrt(n)=0.03148, H/n^(1/3)=0.1911
    n=90000: H=7.14, H/log2(n)=0.434, H/sqrt(n)=0.02380, H/n^(1/3)=0.1594

========================================================================
SUMMARY: Kt(delta(n)|n) Growth Classification
========================================================================

Method (a) BFS program search (n=1..500):
  n=[1,50]: Kt ~ 4.6 bits
  n=[50,200]: Kt ~ 4.6 bits
  n=[200,500]: Kt ~ 4.9 bits

Method (b) Expression fitting:
  Sliding entropy growth: 0.3456 * log(n) + 3.2909
  Power-law exponent: 0.0603

Method (d) Incremental entropy growth:
  Log fit: 0.2196 * log(n)
  Power-law exponent: 0.032787

CLASSIFICATION: SLOW GROWTH, possibly O(log log n)

Interpretation:
  delta(n) has nearly constant conditional entropy per value.
  This suggests the 'random' part is FIXED complexity, NOT growing with n.
  However, this measures STATISTICAL entropy, not Kt directly.
  The per-value entropy being ~constant means each delta(n) carries a fixed
  number of 'surprise bits', but the TOTAL information grows as O(n).

--- Generating plots ---
Plots saved to /apps/aplikacijos/prime-research/experiments/information_theory/kt_complexity/kt_estimation_plots.png
Growth rate plot saved to /apps/aplikacijos/prime-research/experiments/information_theory/kt_complexity/kt_growth_rate.png