================================================================================
NONLINEAR PREDICTION OF delta(n) FROM n
================================================================================

Data: N=100000, mean=-9.70, std=176.02
Range: [-664, 654]
Train: n=1..80000, Test: n=80001..100000
Train delta: mean=-11.75, std=146.07
Test delta:  mean=-1.49, std=263.61

Baselines:
  Naive (predict 0):    RMSE = 263.61
  Mean baseline:        RMSE = 263.81
  Test std:             263.61

================================================================================
EXPERIMENT 1: Direct prediction of delta(n) from features of n
================================================================================
Feature matrix shape: (80000, 26)

--- 1a. Linear regression on n-features ---
  RMSE = 264.7896, R² = -0.008978

--- 1b. Polynomial regression degree 2 ---
  RMSE = 270.8534, R² = -0.055720

--- 1b. Polynomial regression degree 3 ---
  RMSE = 387.6639, R² = -1.162672

--- 1b. Polynomial regression degree 4 ---
  RMSE = 822.9222, R² = -8.745352

--- 1b. Polynomial regression degree 5 ---
  RMSE = 1691.3431, R² = -40.166454

--- 1b. Polynomial regression degree 6 ---
  RMSE = 3084.1511, R² = -135.883564

--- 1c. Random Forest (n_estimators=200, max_depth=20) ---
  RMSE = 285.0373, R² = -0.169186 (time: 10.4s)
  Top 10 features:
    n*log(log(n))       : 0.1256
    n                   : 0.1206
    n^(1/3)             : 0.1175
    log(n)              : 0.1142
    sqrt(n)             : 0.1141
    1/log(n)            : 0.1049
    sin(2pi*n/10000)    : 0.0656
    cos(2pi*n/10000)    : 0.0592
    sin(2pi*n/5000)     : 0.0553
    cos(2pi*n/2000)     : 0.0318

--- 1d. Gradient Boosting (n_estimators=500, max_depth=5) ---
  RMSE = 280.9565, R² = -0.135947 (time: 258.3s)

--- Experiment 1 Summary ---
Method                          RMSE         R²   vs naive
------------------------------------------------------------
Linear                      264.7896  -0.008978     1.0045
Poly-2                      270.8534  -0.055720     1.0275
GradientBoosting            280.9565  -0.135947     1.0658
RandomForest                285.0373  -0.169186     1.0813
Poly-3                      387.6639  -1.162672     1.4706
Poly-4                      822.9222  -8.745352     3.1217
Poly-5                     1691.3431 -40.166454     6.4160
Poly-6                     3084.1511 -135.883564    11.6995
Naive (predict 0)           263.6130        ---     1.0000
Mean baseline               263.8087        ---     1.0007

================================================================================
EXPERIMENT 2: Local prediction -- AR(k) vs AR(k) + n-features
================================================================================

--- AR lag K=1 ---
  Pure AR(1):              RMSE=11.5971, R²=0.998065
  AR(1)+n (Ridge):         RMSE=11.5976, R²=0.998064
  AR(1)+n (RandomForest):  RMSE=29.1781, R²=0.987748
  Improvement from n:        -0.00% (linear), -151.60% (RF)

--- AR lag K=5 ---
  Pure AR(5):              RMSE=11.5691, R²=0.998074
  AR(5)+n (Ridge):         RMSE=11.5691, R²=0.998074
  AR(5)+n (RandomForest):  RMSE=29.2523, R²=0.987686
  Improvement from n:        0.00% (linear), -152.85% (RF)

--- AR lag K=10 ---
  Pure AR(10):              RMSE=11.5653, R²=0.998075
  AR(10)+n (Ridge):         RMSE=11.5651, R²=0.998075
  AR(10)+n (RandomForest):  RMSE=28.9458, R²=0.987943
  Improvement from n:        0.00% (linear), -150.28% (RF)

--- AR lag K=20 ---
  Pure AR(20):              RMSE=11.5669, R²=0.998075
  AR(20)+n (Ridge):         RMSE=11.5665, R²=0.998075
  AR(20)+n (RandomForest):  RMSE=28.6629, R²=0.988177
  Improvement from n:        0.00% (linear), -147.80% (RF)

================================================================================
EXPERIMENT 3: Sign prediction -- can we predict sign(delta(n)) from n?
================================================================================
Train: 45.1% positive
Test:  48.3% positive

Random baseline accuracy: 50.0%
Majority class baseline: 51.7%

Logistic Regression:  acc=55.86%, AUC=0.6133
Random Forest:        acc=49.14%, AUC=0.5385
Gradient Boosting:    acc=48.19%, AUC=0.5322

--- Sign run lengths (context) ---
Number of sign runs: 2267
Mean run length: 44.1
Median run length: 3.0
Max run length: 4524
Min run length: 1

================================================================================
EXPERIMENT 4: Predict delta(n) mod m from n
================================================================================

--- mod 2 ---
  Train distribution: [0.5005875 0.4994125]
  Test distribution:  [0.50875 0.49125]
  Random baseline:    50.0%
  Majority baseline:  50.9%
  RF accuracy:        50.02%
  Lift over random:   1.0003x

--- mod 3 ---
  Train distribution: [0.3315375 0.3337    0.3347625]
  Test distribution:  [0.32935 0.34065 0.33   ]
  Random baseline:    33.3%
  Majority baseline:  33.0%
  RF accuracy:        33.86%
  Lift over random:   1.0158x

--- mod 4 ---
  Train distribution: [0.2511625 0.2486625 0.249425  0.25075  ]
  Test distribution:  [0.25555 0.24425 0.2532  0.247  ]
  Random baseline:    25.0%
  Majority baseline:  25.6%
  RF accuracy:        24.80%
  Lift over random:   0.9918x

--- mod 5 ---
  Train distribution: [0.199675  0.199525  0.2005125 0.1993875 0.2009   ]
  Test distribution:  [0.1987  0.2006  0.19625 0.19765 0.2068 ]
  Random baseline:    20.0%
  Majority baseline:  20.7%
  RF accuracy:        19.96%
  Lift over random:   0.9980x

--- mod 6 ---
  Train distribution: [0.1651125 0.166375  0.16815   0.166425  0.167325  0.1666125]
  Test distribution:  [0.16735 0.1668  0.16755 0.162   0.17385 0.16245]
  Random baseline:    16.7%
  Majority baseline:  16.8%
  RF accuracy:        16.79%
  Lift over random:   1.0077x

================================================================================
EXPERIMENT 5: Residual analysis
================================================================================

Best predictor: GradientBoosting
Residual stats: mean=-78.0915, std=269.8857
Residual range: [-801.11, 672.12]
Original delta std: 263.61
Residual std / delta std: 1.023811
Variance explained: -0.048189

--- Autocorrelation of residuals ---
  lag    1: autocorr = 0.999075
  lag    2: autocorr = 0.998238
  lag    3: autocorr = 0.997424
  lag    5: autocorr = 0.995866
  lag   10: autocorr = 0.992139
  lag   20: autocorr = 0.984328
  lag   50: autocorr = 0.961576
  lag  100: autocorr = 0.934068

--- Autocorrelation of original delta (test set) ---
  lag    1: autocorr = 0.999083
  lag    2: autocorr = 0.998275
  lag    3: autocorr = 0.997511
  lag    5: autocorr = 0.996091
  lag   10: autocorr = 0.992897
  lag   20: autocorr = 0.986549
  lag   50: autocorr = 0.969886
  lag  100: autocorr = 0.948730

--- Compressibility (gzip) ---
  Raw size:          40000 bytes
  Original delta:    23631 bytes (0.5908 ratio)
  Residuals:         25745 bytes (0.6436 ratio)
  Random baseline:   31486 bytes (0.7872 ratio)
  Residuals closer to random? YES

--- Spectral energy concentration ---
  Top   10 freq: delta=0.9223, residual=0.9193
  Top   50 freq: delta=0.9775, residual=0.9784
  Top  100 freq: delta=0.9881, residual=0.9890
  Top  500 freq: delta=0.9973, residual=0.9976

================================================================================
CONCLUSIONS
================================================================================

1. DIRECT PREDICTION (n -> delta(n)):
   Best method: Linear, RMSE=264.7896, R²=-0.008978
   Naive RMSE: 263.6130
   Improvement over naive: -0.45%

   VERDICT: CANNOT predict delta(n) from n

2. LOCAL PREDICTION: Adding n-features to AR models provides
   negligible improvement over pure AR.

3. SIGN PREDICTION: Best accuracy = 55.9% vs 50% random.
   Some structure in sign vs n

4. MOD PREDICTION: At or near random baseline for all moduli.

5. RESIDUALS: After removing best predictor's contribution,
   residuals have std=269.89 vs original std=263.61.
   Compression ratio: 0.6436 (random baseline: 0.7872)

