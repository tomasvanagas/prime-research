BINARY STRING STRUCTURE AND CLUSTERING ANALYSIS OF DELTA(n)
delta(n) = p(n) - round(R^{-1}(n))
Data: /apps/aplikacijos/prime-research/experiments/information_theory/kt_delta_values.npy
N = 100000
Range: [-664, 654]
Mean: -9.6985, Std: 176.0207
Median: -12.0
Fraction positive: 0.457070
Fraction zero: 0.003380
Fraction negative: 0.539550

========================================================================
  EXPERIMENT 1: BINARY STRING ANALYSIS
========================================================================

Delta range: [-664, 654]
Binary width (two's complement): 11 bits

--- Bit bias (fraction with bit=1) ---
  Bit  0 (2^10): 0.539550  (bias from 0.5: +0.039550)
  Bit  1 (2^ 9): 0.535000  (bias from 0.5: +0.035000)
  Bit  2 (2^ 8): 0.541360  (bias from 0.5: +0.041360)
  Bit  3 (2^ 7): 0.515940  (bias from 0.5: +0.015940)
  Bit  4 (2^ 6): 0.508780  (bias from 0.5: +0.008780)
  Bit  5 (2^ 5): 0.505390  (bias from 0.5: +0.005390)
  Bit  6 (2^ 4): 0.501170  (bias from 0.5: +0.001170)
  Bit  7 (2^ 3): 0.500540  (bias from 0.5: +0.000540)
  Bit  8 (2^ 2): 0.499160  (bias from 0.5: -0.000840)
  Bit  9 (2^ 1): 0.500180  (bias from 0.5: +0.000180)
  Bit 10 (2^ 0): 0.497780  (bias from 0.5: -0.002220)

--- Comparison to Gaussian random with same mean/std ---
  Bit  0: delta=0.539550  random=0.524860  diff=+0.014690
  Bit  1: delta=0.535000  random=0.524050  diff=+0.010950
  Bit  2: delta=0.541360  random=0.510030  diff=+0.031330
  Bit  3: delta=0.515940  random=0.503790  diff=+0.012150
  Bit  4: delta=0.508780  random=0.503060  diff=+0.005720
  Bit  5: delta=0.505390  random=0.500000  diff=+0.005390
  Bit  6: delta=0.501170  random=0.498130  diff=+0.003040
  Bit  7: delta=0.500540  random=0.498380  diff=+0.002160
  Bit  8: delta=0.499160  random=0.497890  diff=+0.001270
  Bit  9: delta=0.500180  random=0.499870  diff=+0.000310
  Bit 10: delta=0.497780  random=0.501060  diff=-0.003280

--- Bit-level entropy per position (bits, max=1.0) ---
  Bit  0: H = 0.995482
  Bit  1: H = 0.996463
  Bit  2: H = 0.995058
  Bit  3: H = 0.999267
  Bit  4: H = 0.999778
  Bit  5: H = 0.999916
  Bit  6: H = 0.999996
  Bit  7: H = 0.999999
  Bit  8: H = 0.999998
  Bit  9: H = 1.000000
  Bit 10: H = 0.999986

--- Bit-pair correlations (top 15 by |r|) ---
  Bits ( 0, 1): r = +0.975887
  Bits ( 0, 2): r = +0.742057
  Bits ( 1, 2): r = +0.718146
  Bits ( 0, 3): r = +0.390319
  Bits ( 1, 3): r = +0.367135
  Bits ( 2, 3): r = +0.250990
  Bits ( 0, 4): r = +0.212570
  Bits ( 1, 4): r = +0.201656
  Bits ( 2, 4): r = +0.139166
  Bits ( 0, 5): r = +0.107089
  Bits ( 1, 5): r = +0.102483
  Bits ( 3, 4): r = +0.093942
  Bits ( 2, 5): r = +0.067564
  Bits ( 0, 6): r = +0.049490
  Bits ( 1, 6): r = +0.046450

--- Joint entropy for adjacent bit pairs ---
  Bits (0,1): H_joint = 1.088505 bits (indep would be ~2.0)
  Bits (1,2): H_joint = 1.579499 bits (indep would be ~2.0)
  Bits (2,3): H_joint = 1.948421 bits (indep would be ~2.0)
  Bits (3,4): H_joint = 1.992669 bits (indep would be ~2.0)
  Bits (4,5): H_joint = 1.999654 bits (indep would be ~2.0)
  Bits (5,6): H_joint = 1.999876 bits (indep would be ~2.0)
  Bits (6,7): H_joint = 1.999995 bits (indep would be ~2.0)
  Bits (7,8): H_joint = 1.999997 bits (indep would be ~2.0)
  Bits (8,9): H_joint = 1.999993 bits (indep would be ~2.0)
  Bits (9,10): H_joint = 1.999974 bits (indep would be ~2.0)

========================================================================
  EXPERIMENT 2: DELTA CLUSTERING
========================================================================

k= 2: silhouette = 0.383590
       cluster sizes: [np.int64(51181), np.int64(48819)]
       n mod 6 distribution per cluster:
         cluster 0: 0:0.167 1:0.167 2:0.167 3:0.167 4:0.167 5:0.167
         cluster 1: 0:0.167 1:0.167 2:0.167 3:0.167 4:0.167 5:0.167

k= 4: silhouette = 0.389060
       cluster sizes: [np.int64(40364), np.int64(26522), np.int64(17156), np.int64(15958)]
       n mod 6 distribution per cluster:
         cluster 0: 0:0.167 1:0.166 2:0.167 3:0.167 4:0.166 5:0.167
         cluster 1: 0:0.166 1:0.167 2:0.167 3:0.167 4:0.167 5:0.167
         cluster 2: 0:0.166 1:0.167 2:0.167 3:0.167 4:0.167 5:0.167
         cluster 3: 0:0.167 1:0.167 2:0.166 3:0.166 4:0.167 5:0.167
       n mod 30 distribution per cluster (top 5 residues):
         cluster 0: 20:584 18:582 17:581 19:580 27:576
         cluster 1: 1:1354 4:1352 3:1352 2:1349 26:1348
         cluster 2: 29:545 26:545 28:543 27:541 7:538
         cluster 3: 16:896 14:896 8:894 10:893 6:893

k= 8: silhouette = 0.410385
       cluster sizes: [np.int64(24336), np.int64(15122), np.int64(13581), np.int64(13451), np.int64(12048), np.int64(8984), np.int64(6423), np.int64(6055)]
       n mod 6 distribution per cluster:
         cluster 0: 0:0.167 1:0.167 2:0.166 3:0.168 4:0.167 5:0.165
         cluster 1: 0:0.167 1:0.166 2:0.166 3:0.167 4:0.167 5:0.167
         cluster 2: 0:0.167 1:0.167 2:0.167 3:0.165 4:0.166 5:0.168
         cluster 3: 0:0.166 1:0.166 2:0.167 3:0.166 4:0.167 5:0.168
         cluster 4: 0:0.166 1:0.167 2:0.167 3:0.167 4:0.167 5:0.165
         cluster 5: 0:0.166 1:0.166 2:0.166 3:0.167 4:0.167 5:0.167
         cluster 6: 0:0.167 1:0.167 2:0.169 3:0.165 4:0.166 5:0.166
         cluster 7: 0:0.166 1:0.167 2:0.167 3:0.168 4:0.167 5:0.164

k=16: silhouette = 0.394179
       cluster sizes: [np.int64(13086), np.int64(9108), np.int64(8359), np.int64(8093), np.int64(7333), np.int64(7258), np.int64(7153), np.int64(6074), np.int64(5630), np.int64(5003), np.int64(4988), np.int64(4581), np.int64(4460), np.int64(3514), np.int64(3474), np.int64(1886)]

--- Chi-squared test: cluster labels (k=4) vs n mod 6 ---
  chi2 = 0.46, dof = 15, p-value = 1.000000e+00


========================================================================
  EXPERIMENT 3: TRANSITION PATTERNS
========================================================================

Transition matrix: 30 x 30 bins (quantile-based)
Row sums: mean=1.000000, std=0.000000
Col sums: mean=1.000000, std=0.017171, min=0.959246, max=1.037587
Is approximately doubly stochastic? YES (col_sum std=0.017171)

--- Eigenvalue spectrum (top 10 by magnitude) ---
  lambda_0: |lambda| = 1.000000
  lambda_1: |lambda| = 0.996658
  lambda_2: |lambda| = 0.991356
  lambda_3: |lambda| = 0.980952
  lambda_4: |lambda| = 0.969326
  lambda_5: |lambda| = 0.953931
  lambda_6: |lambda| = 0.935697
  lambda_7: |lambda| = 0.911907
  lambda_8: |lambda| = 0.887786
  lambda_9: |lambda| = 0.862318
Spectral gap: 0.003342

--- Shuffled (iid baseline) eigenvalues (top 5) ---
  lambda_0: |lambda| = 1.000000
  lambda_1: |lambda| = 0.016177
  lambda_2: |lambda| = 0.016054
  lambda_3: |lambda| = 0.016054
  lambda_4: |lambda| = 0.015391
Shuffled spectral gap: 0.983823

Frobenius distance from uniform: actual=3.734090, shuffled=0.093314
Mean diagonal of P: 0.634621 (uniform would be 0.033333)
Diagonal enhancement ratio: 19.0386x

========================================================================
  EXPERIMENT 4: INFORMATION-THEORETIC MEASURES
========================================================================

--- Lempel-Ziv complexity ---
  Sign sequence (N=50000): LZ = 1298
  Random baseline:         LZ = 4813
  Ratio (delta/random):    0.269686
  Quantized (16 bins, N=50000): LZ = 5698
  Random 16-symbol baseline:    LZ = 13456
  Ratio: 0.423454

--- Approximate Entropy (ApEn) ---
  m=2: ApEn(delta) = 0.710533,  ApEn(random) = 2.671473
  m=3: ApEn(delta) = 0.691542,  ApEn(random) = 1.493815

--- Sample Entropy (SampEn, binned approx) ---
  m=2, r=0.2*std: SampEn(delta) = 0.892607, SampEn(random) = 2.873445
  m=3, r=0.2*std: SampEn(delta) = 0.874372, SampEn(random) = 2.881551

--- Permutation Entropy ---
  Order 3: H = 2.429096 / 2.584963 (normalized: 0.939703), observed 6/6 permutations
  Order 4: H = 4.122062 / 4.584963 (normalized: 0.899039), observed 24/24 permutations
  Order 5: H = 6.040770 / 6.906891 (normalized: 0.874600), observed 120/120 permutations
  Order 6: H = 8.158733 / 9.491853 (normalized: 0.859551), observed 711/720 permutations
  Order 7: H = 10.341671 / 12.299208 (normalized: 0.840840), observed 3862/5040 permutations

========================================================================
  EXPERIMENT 5: RUNS AND PATTERNS
========================================================================

--- Run length distribution (same sign) ---
  Total runs: 2269
  Mean run length: 43.9233
  Std run length: 245.5589
  Max run length: 4524
  Median run length: 3.0

  Run length distribution (top 15):
    length   1:    657 (28.96%)
    length   2:    314 (13.84%)
    length   3:    184 (8.11%)
    length   4:    121 (5.33%)
    length   5:    112 (4.94%)
    length   6:     65 (2.86%)
    length   7:     73 (3.22%)
    length   8:     45 (1.98%)
    length   9:     40 (1.76%)
    length  10:     29 (1.28%)
    length  11:     24 (1.06%)
    length  12:     24 (1.06%)
    length  13:     26 (1.15%)
    length  14:     22 (0.97%)
    length  15:     18 (0.79%)

  Geometric fit: p = 0.022767
  Comparison (observed vs geometric):
    length 1: observed=0.289555, geometric=0.022767, ratio=12.7182
    length 2: observed=0.138387, geometric=0.022249, ratio=6.2200
    length 3: observed=0.081093, geometric=0.021742, ratio=3.7298
    length 4: observed=0.053327, geometric=0.021247, ratio=2.5099
    length 5: observed=0.049361, geometric=0.020763, ratio=2.3773
    length 6: observed=0.028647, geometric=0.020291, ratio=1.4118
    length 7: observed=0.032173, geometric=0.019829, ratio=1.6225
    length 8: observed=0.019833, geometric=0.019377, ratio=1.0235
    length 9: observed=0.017629, geometric=0.018936, ratio=0.9310
    length 10: observed=0.012781, geometric=0.018505, ratio=0.6907
  KS test vs geometric: stat=0.524111, p-value=0.000000e+00

--- Sign n-gram analysis ---

  3-grams:
    Total possible: 27
    Observed: 27
    Missing: 0
    Most common:
      ---: 51990 (51.991%)
      +++: 43791 (43.792%)
      --+: 738 (0.738%)
      -++: 697 (0.697%)
      +--: 674 (0.674%)
    Least common:
      00-: 5 (0.005%)
      -00: 4 (0.004%)
      00+: 4 (0.004%)
      +00: 4 (0.004%)
      000: 1 (0.001%)

  4-grams:
    Total possible: 81
    Observed: 74
    Missing: 7
      Missing: -000
      Missing: 00-+
      Missing: 0000
      Missing: 000+
      Missing: 0+00
      Missing: +-00
      Missing: +000
    Most common:
      ----: 51319 (51.321%)
      ++++: 43164 (43.165%)
      ---+: 596 (0.596%)
      -+++: 558 (0.558%)
      +---: 538 (0.538%)
    Least common:
      0-00: 1 (0.001%)
      00+-: 1 (0.001%)
      -0+0: 1 (0.001%)
      0+-0: 1 (0.001%)
      00+0: 1 (0.001%)

--- Delta value n-gram forbidden patterns (4 bins: Q1/Q2/Q3/Q4) ---
  3-grams: 23/64 observed (expected ~1562.5 per pattern)
    Chi-squared uniformity: chi2=1238966.25, dof=63, p=0.000000e+00
  4-grams: 48/256 observed (expected ~390.6 per pattern)
    Chi-squared uniformity: chi2=4928973.50, dof=255, p=0.000000e+00

========================================================================
  SUMMARY
========================================================================

See individual sections above for detailed findings.
