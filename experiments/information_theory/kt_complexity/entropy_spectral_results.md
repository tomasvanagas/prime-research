======================================================================
1. BLOCK ENTROPY SCALING
======================================================================

   k |        H_k |  h_k=H_k/k |   #symbols |  h_k(shuf) | #sym(shuf)
----------------------------------------------------------------------
   1 |     9.4521 |     9.4521 |       1236 |     9.4521 |       1236
   2 |    13.8226 |     6.9113 |      19777 |     7.7297 |      46365
   3 |    14.9167 |     4.9722 |      31578 |     5.0082 |      33331
   4 |    14.6056 |     3.6514 |      24949 |     3.6524 |      25000
   5 |    14.2875 |     2.8575 |      19998 |     2.8575 |      20000
   8 |    13.6096 |     1.7012 |      12500 |     1.7012 |      12500
  10 |    13.2877 |     1.3288 |      10000 |     1.3288 |      10000
  16 |    12.6096 |     0.7881 |       6250 |     0.7881 |       6250
  20 |    12.2877 |     0.6144 |       5000 |     0.6144 |       5000
  32 |    11.6096 |     0.3628 |       3125 |     0.3628 |       3125
  50 |    10.9658 |     0.2193 |       2000 |     0.2193 |       2000

Entropy rate fit: h_k ~ -0.1275 + 14.9554/k
Estimated asymptotic entropy rate h_inf = -0.1275 bits/symbol
Shuffled: h_k ~ -0.1277 + 14.9585/k
Shuffled asymptotic h_inf = -0.1277 bits/symbol

Single-symbol entropy H_1 = 9.4521 bits
Entropy reduction h_1 -> h_50: 9.4521 -> 0.2193 (97.7% reduction)

======================================================================
2. POWER SPECTRAL DENSITY
======================================================================

Spectral slope (full range): beta = 1.6872 (PSD ~ f^-1.6872)
Shuffled spectral slope:     beta = 0.0049 (PSD ~ f^-0.0049)
Low-freq slope (f<0.01):     beta = 1.5783
Mid-freq slope (0.01<f<0.1): beta = 1.8711
High-freq slope (f>0.1):     beta = 1.3817

Mean PSD: 30983.92
Median PSD: 50.10
PSD ratio (low 1% freqs / high 1% freqs): 104223.31
Log-binned PSD coefficient of variation: 2.8323
VERDICT: Spectrum shows strong 1/f^1.69 coloring

======================================================================
3. WAVELET ANALYSIS
======================================================================

Wavelet: db4, decomposition levels: 12
 Level |  #coeffs |        Std |   Kurtosis |   Skewness |  Normal?
----------------------------------------------------------------------
   A12 |       31 |  6937.2298 |     0.3488 |     0.4767 | p=3.02e-01
   D12 |       31 |  5318.0084 |     0.8288 |    -0.5503 | p=1.37e-01
   D11 |       55 |  2614.1136 |     0.2686 |     0.2801 | p=4.63e-01
   D10 |      104 |  1456.7350 |     0.0128 |     0.5107 | p=9.28e-02
    D9 |      202 |   815.9058 |     0.4686 |    -0.1138 | p=2.90e-01
    D8 |      397 |   448.9893 |     0.1896 |    -0.1151 | p=4.24e-01
    D7 |      788 |   263.8135 |     0.4244 |     0.0112 | p=9.19e-02
    D6 |     1569 |   129.5610 |     0.2493 |     0.2179 | p=3.58e-04
    D5 |     3131 |    69.7542 |     0.1946 |     0.1028 | p=6.99e-03
    D4 |     6256 |    35.6438 |     0.3977 |     0.2008 | p=4.01e-16
    D3 |    12506 |    19.1732 |     0.6921 |     0.2575 | p=1.01e-60
    D2 |    25005 |    10.1364 |     1.0245 |     0.3318 | p=3.73e-204
    D1 |    50003 |     6.3059 |     2.2192 |     0.8044 | p=0.00e+00

Energy distribution across scales:
  A12: 44.82%
  D12: 26.90%
  D11: 11.79%
  D10: 6.63%
  D9: 4.05%
  D8: 2.40%
  D7: 1.65%
  D6: 0.79%
  D5: 0.46%
  D4: 0.24%
  D3: 0.14%
  D2: 0.08%
  D1: 0.06%

======================================================================
4. MUTUAL INFORMATION SCALING
======================================================================

   Lag k |   MI(bits) |    MI_shuf
----------------------------------------
       1 |   3.846658 |   0.040528
       2 |   3.530296 |   0.041409
       5 |   2.969284 |   0.041279
      10 |   2.542354 |   0.041190
      20 |   2.117566 |   0.041009
      50 |   1.597132 |   0.041839
     100 |   1.252125 |   0.040990
     200 |   0.956254 |   0.041684
     500 |   0.603003 |   0.042525
    1000 |   0.354266 |   0.040857

MI baseline (shuffled mean): 0.041331 bits
Power law fit (bias-corrected): MI ~ k^(-0.343)
Exponential fit (bias-corrected): MI ~ exp(-k/436.1)
R^2 power law:   0.9419
R^2 exponential: 0.8150
VERDICT: MI decays as POWER LAW ~ k^(-0.34)

======================================================================
5. SUMMARY AND INTERPRETATION
======================================================================

Delta(n) = p(n) - round(R^{-1}(n)), N = 100000
Range: [-664, 654], Mean: -9.6985, Std: 176.0207

Block Entropy:
  h_1 = 9.4521 bits/symbol (1236 distinct values)
  h_2 = 6.9113 bits/symbol vs shuffled 7.7297 => 10.6% reduction
  NOTE: For k >= 4, N/k <= alphabet^k so all blocks are unique and
  H_k ~ log2(N/k) for both real and shuffled. The entropy "decrease"
  is a sampling artifact, not structure. The meaningful signal is at
  k=2 where there is 10.6% entropy reduction vs shuffled, confirming
  short-range correlations in delta.

Power Spectrum:
  PSD ~ f^(-1.69) -- strong 1/f coloring, NOT white noise
  Shuffled PSD is flat (beta ~ 0.005) as expected
  Low-freq slope: 1.58, Mid-freq: 1.87, High-freq: 1.38
  This is close to Brownian noise (beta=2) but slightly less steep.
  Interpretation: delta(n) has long-range power-law correlations,
  consistent with the oscillatory sum over zeta zeros.

Wavelet Analysis:
  Energy heavily concentrated in coarse scales: A12+D12+D11 = 83.5%
  Fine-scale coefficients (D1-D4) are non-Gaussian (excess kurtosis
  up to 2.2, p << 0.001 for normality), indicating heavy tails at
  small scales. Coarse scales (D7-D12) are consistent with Gaussian.
  This is the signature of a process whose large-scale behavior is
  smooth (Central Limit Theorem regime) but whose local fluctuations
  have intermittent, heavy-tailed structure.

Mutual Information:
  MI(lag=1) = 3.85 bits -- extremely high (vs 0.04 bits shuffled)
  MI decays as POWER LAW: MI ~ k^(-0.34) with R^2 = 0.94
  Still significant at lag 1000: MI = 0.35 bits
  This slow power-law decay (not exponential) confirms long-range
  dependence, consistent with the 1/f^1.69 spectrum.

Key Conclusions:
  1. delta(n) is NOT white noise / iid -- it has strong structure
  2. PSD ~ f^(-1.69): long-range power-law correlations
  3. MI ~ k^(-0.34): very slow information decay
  4. Fine-scale non-Gaussianity: intermittent, heavy-tailed local behavior
  5. These properties are consistent with delta encoding contributions
     from ~10^48 zeta zeros with correlated (GUE) phases -- the signal
     is smooth at large scales but complex and non-Gaussian locally
  6. The long-range dependence means compression IS possible in principle
     (unlike iid data), but the power-law decay means no finite-memory
     model (e.g., ARMA) will capture the full structure
