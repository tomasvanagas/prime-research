# Zeta Zeros mod Various Constants: Pattern Search

**Date**: 2026-04-05
**N**: 1000 zeros, range [14.1347, 1419.4225]

## 1. Per-Modulus Uniformity Tests

| Modulus | Value | KS stat | KS p-val | Circ Mean | Circ Var | R_bar | Rayleigh Z | Rayleigh p |
|---------|-------|---------|----------|-----------|----------|-------|------------|------------|
| 1 | 1.000000 | 0.0163 | 0.9497 | 0.0391 | 0.9896 | 0.0104 | 0.1090 | 0.8968 |
| pi | 3.141593 | 0.0113 | 0.9994 | 2.2078 | 0.9977 | 0.0023 | 0.0051 | 0.9949 |
| log(2pi) | 1.837877 | 0.0137 | 0.9904 | 0.4059 | 0.9879 | 0.0121 | 0.1463 | 0.8639 |
| 2pi | 6.283185 | 0.0122 | 0.9980 | 3.3055 | 0.9997 | 0.0003 | 0.0001 | 0.9999 |
| sqrt(2pi) | 2.506628 | 0.0228 | 0.6694 | 2.3452 | 0.9985 | 0.0015 | 0.0023 | 0.9977 |
| e | 2.718282 | 0.0158 | 0.9602 | 0.2697 | 0.9967 | 0.0033 | 0.0111 | 0.9890 |
| log(2) | 0.693147 | 0.0170 | 0.9290 | 0.5948 | 0.9935 | 0.0065 | 0.0424 | 0.9585 |
| log(3) | 1.098612 | 0.0276 | 0.4247 | 0.0815 | 0.9795 | 0.0205 | 0.4194 | 0.6575 |
| log(5) | 1.609438 | 0.0154 | 0.9682 | 1.4079 | 0.9931 | 0.0069 | 0.0477 | 0.9534 |
| log(7) | 1.945910 | 0.0200 | 0.8114 | 1.1578 | 0.9929 | 0.0071 | 0.0507 | 0.9506 |

### Significance Summary

- KS test: NO modulus shows significant deviation from uniformity
- Rayleigh test: NO modulus shows significant circular non-uniformity

## 2. Mod 1 Detailed Analysis

### Weyl Criterion (Equidistribution)

- Checked |S_k| = |(1/N) sum exp(2 pi i k gamma_n)| for k=1..50
- Maximum |S_k| = 0.080124 at k=26
- Mean |S_k| = 0.027585
- Expected for random: 1/sqrt(N) = 0.031623
- Weyl sums above 2/sqrt(N): 1/50
- **Conclusion**: Sequence appears equidistributed mod 1

### Discrepancy

- D_N = 0.016280
- Expected for random uniform: ~ 0.031086
- Ratio D_N / expected: 0.52

### Gap Analysis (mod 1)

- Mean gap: 0.001000 (expected 0.001000)
- Min gap: 0.000001 (ratio to mean: 0.0006)
- Max gap: 0.007120 (ratio to mean: 7.12)
- Normalized gaps vs Exp(1): KS=0.0343, p=0.1845
- Gap distribution consistent with exponential (independent uniform)

## 3. Joint Distribution (mod 1 vs mod pi)

- Chi-squared independence: chi2=62.37, dof=196, p=1.0000
- Pearson r = -0.030078, p = 0.3420
- Spearman rho = -0.030659, p = 0.3328
- **Conclusion**: No evidence of dependence between mod 1 and mod pi

## 4. Overall Conclusions

**No significant patterns found.** The zeta zeros appear uniformly distributed
modulo all tested constants, consistent with theoretical expectations.
This is further evidence that the oscillatory part of pi(x) is pseudorandom.
