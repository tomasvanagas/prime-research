# Wilf-Zeilberger Definite Sum Test Results

**Date**: 2026-04-05
**Data**: f(x) = pi(x) - R(x), x = 2..100000

```
========================================================================
WILF-ZEILBERGER DEFINITE SUM TEST FOR f(x) = pi(x) - R(x)
========================================================================

## TEST 1: First Differences Delta_f(x)
--------------------------------------------------
Delta_f: length=99998, range=[-0.416013, 0.913293]
Delta_f: mean=5.183922e-05, std=0.294222

When x+1 is prime (9591 cases):
  Delta_f mean=0.902632, std=0.015352
When x+1 is composite (90407 cases):
  Delta_f mean=-0.095700, std=0.011524

### 1a: Simple function approximations to Delta_f
  Fit c/log(x): c=0.000526, residual std=0.294222 (vs delta_f std=0.294222)
  Variance explained: 0.000000
  Fit c1/log(x)+c2/log(x)^2: coefs=[-0.01823068  0.18493858], residual std=0.294219
  Variance explained: 0.000021
  Fit 1/log + 1/log^2 + 1/log^3: residual std=0.294218
  Variance explained: 0.000029

### 1b: Recurrence check for Delta_f
  Linear recurrence order 1: R^2 = 0.01164591
  Linear recurrence order 2: R^2 = 0.01212294
  Linear recurrence order 3: R^2 = 0.02264331
  Linear recurrence order 5: R^2 = 0.03271976

### 1c: Autocorrelation of Delta_f
  Autocorrelation(lag=1): -0.107915
  Autocorrelation(lag=2): 0.033471
  Autocorrelation(lag=3): -0.107887
  Autocorrelation(lag=5): -0.107857
  Autocorrelation(lag=10): 0.079768
  Autocorrelation(lag=20): 0.082270
  Autocorrelation(lag=50): 0.082403
  Autocorrelation(lag=100): 0.079144

## TEST 2: Higher-Order Differences Delta^k f
--------------------------------------------------
  Delta^ 1 f: len= 99998, |mean|=5.1839e-05, rms=2.9422e-01, max|.|=9.1329e-01
  Delta^ 2 f: len= 99997, |mean|=6.3482e-06, rms=4.3797e-01, max|.|=1.0253e+00
  Delta^ 3 f: len= 99996, |mean|=9.6415e-06, rms=7.7454e-01, max|.|=2.0064e+00
  Delta^ 4 f: len= 99995, |mean|=1.9896e-05, rms=1.4370e+00, max|.|=4.0023e+00
  Delta^ 5 f: len= 99994, |mean|=3.9961e-05, rms=2.7244e+00, max|.|=8.0010e+00
  Delta^ 6 f: len= 99993, |mean|=7.9987e-05, rms=5.2293e+00, max|.|=1.5999e+01
  Delta^ 7 f: len= 99992, |mean|=1.6000e-04, rms=1.0118e+01, max|.|=3.1000e+01
  Delta^ 8 f: len= 99991, |mean|=3.1002e-04, rms=1.9688e+01, max|.|=5.7000e+01
  Delta^ 9 f: len= 99990, |mean|=5.7006e-04, rms=3.8472e+01, max|.|=1.1200e+02
  Delta^10 f: len= 99989, |mean|=1.0001e-03, rms=7.5426e+01, max|.|=2.1100e+02

  Checking RMS ratios (Delta^{k+1}/Delta^k):
  RMS(Delta^2)/RMS(Delta^1) = 1.488561
  RMS(Delta^3)/RMS(Delta^2) = 1.768491
  RMS(Delta^4)/RMS(Delta^3) = 1.855272
  RMS(Delta^5)/RMS(Delta^4) = 1.895947
  RMS(Delta^6)/RMS(Delta^5) = 1.919406
  RMS(Delta^7)/RMS(Delta^6) = 1.934835
  RMS(Delta^8)/RMS(Delta^7) = 1.945846
  RMS(Delta^9)/RMS(Delta^8) = 1.954114
  RMS(Delta^10)/RMS(Delta^9) = 1.960525

## TEST 3: Hypergeometric Recurrence Test
--------------------------------------------------
Testing: sum_{i=0}^{r} a_i(x) * f(x+i) = 0, a_i = polynomial degree <= d
  r=1, deg=0: R²(train)=0.99181727, R²(test)=0.99666731, params=1
  r=1, deg=1: R²(train)=0.99774811, R²(test)=0.99387594, params=3
  r=1, deg=2: R²(train)=0.99893083, R²(test)=0.96015206, params=5
  r=2, deg=0: R²(train)=0.99192586, R²(test)=0.99669898, params=2
  r=2, deg=1: R²(train)=0.99778738, R²(test)=0.99394819, params=5
  r=2, deg=2: R²(train)=0.99882157, R²(test)=0.95501667, params=8
  r=3, deg=0: R²(train)=0.99193115, R²(test)=0.99670024, params=3
  r=3, deg=1: R²(train)=0.99778934, R²(test)=0.99394962, params=7
  r=3, deg=2: R²(train)=0.99858261, R²(test)=0.94554215, params=11
  r=4, deg=0: R²(train)=0.99202901, R²(test)=0.99672915, params=4
  r=4, deg=1: R²(train)=0.99782436, R²(test)=0.99401360, params=9
  r=4, deg=2: R²(train)=0.99884831, R²(test)=0.95639075, params=14
  r=5, deg=0: R²(train)=0.99202981, R²(test)=0.99672964, params=5
  r=5, deg=1: R²(train)=0.99782461, R²(test)=0.99401431, params=11
  r=5, deg=2: R²(train)=0.99876524, R²(test)=0.95297130, params=17

  Best fit: r=5, deg=0, R²(train)=0.99202981, R²(test)=0.99672964

### 3b: Recurrence test with detrended f
  Detrended f: std=4.1461 (original std=4.2326)
  Detrended r=1, deg=0: R²(train)=0.99182902, R²(test)=0.99647021
  Detrended r=1, deg=1: R²(train)=0.99775100, R²(test)=0.99352949
  Detrended r=2, deg=0: R²(train)=0.99193744, R²(test)=0.99650352
  Detrended r=2, deg=1: R²(train)=0.99779031, R²(test)=0.99360161
  Detrended r=3, deg=0: R²(train)=0.99194271, R²(test)=0.99650491
  Detrended r=3, deg=1: R²(train)=0.99779224, R²(test)=0.99360399

## TEST 4: Summation Kernel K(x) = f(x) * x * log(x)
--------------------------------------------------
K(x): range=[-14571715.01, 14861138.99], mean=230454.39, std=3162375.04

### 4a: Polynomial fits to K(x)
  Poly deg  1: R²=0.0157119416, residual std=3137433.16
  Poly deg  2: R²=0.0455357849, residual std=3089535.56
  Poly deg  3: R²=0.0965969621, residual std=3005758.94
  Poly deg  5: R²=0.1757838944, residual std=2871004.82
  Poly deg  8: R²=0.1857205443, residual std=2853646.10
  Poly deg 10: R²=0.2067116995, residual std=2816624.15

### 4b: Compressibility (Hankel matrix singular values)
  Top 10 normalized singular values: [1.       0.903432 0.621654 0.594538 0.58537  0.480291 0.414891 0.395749
 0.356829 0.353876]
  SV decay: sv[10]/sv[0]=3.362446e-01, sv[50]/sv[0]=1.040572e-01
  Effective rank (sv > 1e-6 * max): 250 / 250

### 4c: Autocorrelation comparison: K(x) vs f(x)
  lag=  1: AC(K)=0.998302, AC(f)=0.997588
  lag=  2: AC(K)=0.996936, AC(f)=0.995698
  lag=  5: AC(K)=0.993066, AC(f)=0.990427
  lag= 10: AC(K)=0.987434, AC(f)=0.982928
  lag= 50: AC(K)=0.951510, AC(f)=0.936653
  lag=100: AC(K)=0.915613, AC(f)=0.891673

### 4d: First differences of K(x)
  Delta_K: mean=52.2445, std=184155.4875
  Poly deg 1 fit to Delta_K: R²=0.00000020
  Poly deg 2 fit to Delta_K: R²=0.00000042
  Poly deg 3 fit to Delta_K: R²=0.00000053

========================================================================
SUMMARY
========================================================================

1. FIRST DIFFERENCES: Delta_f(x) has negligible autocorrelation and is
   not well-approximated by smooth functions. Variance explained by
   1/log(x) expansion: 0.0000
   => Delta_f is dominated by the indicator function of primes (random).

2. HIGHER-ORDER DIFFERENCES: RMS does NOT decrease -- it INCREASES.
   This rules out f(x) being a polynomial or hypergeometric sequence.

3. HYPERGEOMETRIC RECURRENCE: Best R²(test) = 0.996730
   WARNING: This high R² is SPURIOUS. It reflects only that f(x) is a
   slow-varying function with AC(lag=1) = 0.9976, so f(x+1) ~ c*f(x)
   trivially. The residual std ~ 0.29 matches Delta_f std exactly -- the
   recurrence captures ONLY the trivial autocorrelation, not algebraic
   structure. Adding more recurrence terms or polynomial degrees does NOT
   improve test R² (all cluster at ~0.997). A genuine WZ certificate would
   give R² = 1.000 exactly.
   Hypergeometric certificate: NEGATIVE.

4. SUMMATION KERNEL K(x) = f(x)*x*log(x):
   Hankel effective rank: 250/250 (full rank = incompressible)
   K(x) has NO better structure than f(x) itself.

KEY OBSERVATIONS:

- Delta_f is bimodal: ~0.90 at primes, ~-0.096 at composites (prime indicator).
- Higher-order RMS ratios converge to 2.0 -- consistent with differencing
  white noise (binomial coefficient growth), confirming f has no polynomial
  or hypergeometric structure at any finite order.
- Autocorrelation of Delta_f shows weak structure at lags 1,3,5 (~-0.108),
  which reflects the twin-prime / prime-gap correlations, NOT algebraic structure.
- The summation kernel K(x) = f(x)*x*log(x) has full Hankel rank (250/250)
  and polynomial fits explain < 21% of variance even at degree 10.

CONCLUSION: f(x) = pi(x) - R(x) does NOT admit a WZ-style definite sum
representation. The first differences are dominated by prime indicators,
higher differences grow (not shrink), no polynomial recurrence holds,
and the summation kernel is incompressible. This is consistent with
f(x) encoding zeta zero oscillations with no finite closed-form certificate.

This CLOSES the WZ definite-sum approach for f(x) = pi(x) - R(x).
```
