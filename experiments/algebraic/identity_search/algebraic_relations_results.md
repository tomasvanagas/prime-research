# Algebraic Relations Search: f(x) = pi(x) - R(x)

**Date:** 2026-04-05  
**Script:** `experiments/algebraic/identity_search/algebraic_relations.py`  
**Data:** `fx_data.npz`, x = 2..100000

---

## 1. Bernoulli Number Relations

**Test:** Does f(x) correlate with sum_{k=1}^{K} B_{2k} * x^{-2k}?

**Result: NO.** The Bernoulli correction terms are astronomically small (order 10^{-7} at x=1000) while f(x) is O(1)-O(10). Pearson r ~ -0.006, p ~ 0.84 for all K tested (1, 2, 3, 5, 10, 20). Even with sqrt(x) scaling, no correlation (r ~ -0.002).

**Interpretation:** The Bernoulli numbers encode contributions from the trivial zeros of zeta (at s = -2, -4, ...). These contribute negligibly to the prime counting function residual in this range. The oscillatory part f(x) is dominated by the nontrivial zeros.

---

## 2. Zeta Value Relations (PSLQ)

**Test:** At sample points, search for integer relations among [f(x), zeta(2..7), log(x)*zeta(2), log(x)*zeta(3), sqrt(x), 1].

**Results:** PSLQ found relations at all sample points, but they are **x-dependent** -- different coefficients at each x. Examples:

- x=100: `sqrt(x) - 10 = 0` (trivially, sqrt(100) = 10; PSLQ ignores f(x))
- x=10000: `sqrt(x) - 100 = 0` (same trivial relation)
- x=500: `-2016*f(x) + 80*zeta(2) + 1436*zeta(3) + ... = 0` (large coefficients, x-specific)
- x=1000: `-1402*f(x) + 4103*zeta(2) - 1592*zeta(3) + ... = 0` (completely different coefficients)

**Interpretation:** The PSLQ relations at perfect-square x values just find sqrt(x) in Z. At other points, the relations are **spurious** -- with 11 basis elements and float64 precision (~15 digits), PSLQ will always find relations with coefficients up to ~10^6. The inconsistency across x values confirms these are numerical artifacts, not algebraic identities.

**Conclusion: NO universal algebraic relation between f(x) and zeta values.**

---

## 3. Dirichlet L-function Values (PSLQ)

**Test:** At sample points, search for integer relations among [f(x), L(1,chi_3), L(1,chi_4), L(1,chi_5), L(1,chi_7), log(x), sqrt(x), 1].

**L-values computed:**
- L(1, chi_3) = 0.6046
- L(1, chi_4) = 0.7854 (= pi/4)
- L(1, chi_5) = 0.4304
- L(1, chi_7) = 1.1874

**Results:** Same pattern as zeta -- x-dependent relations with large coefficients:

- x=100, 10000: trivial sqrt relation only
- x=1000: `311081*f(x) + 568324*L(1,chi_3) + ... = 0` (huge coefficients)
- x=50000: `-416367*f(x) + 149342*L(1,chi_3) + ... = 0` (completely different)

**Conclusion: NO universal algebraic relation between f(x) and Dirichlet L-values.**

---

## 4. Modular Form Coefficients (Ramanujan Tau)

**Test:** Correlate f(n) with tau(n) and sum_{d|n} tau(d) for n=2..100.

**Results:**

| Comparison | Pearson r | p-value | Spearman r |
|---|---|---|---|
| f(n) vs tau(n) | +0.010 | 0.93 | +0.054 |
| f(n) vs sum_{d\|n} tau(d) | +0.008 | 0.94 | +0.053 |
| f(n) vs tau(n)/n^{11/2} | +0.099 | 0.33 | -- |

**Conclusion: NO correlation.** The Ramanujan tau function and f(x) are essentially independent. This is expected: tau(n) encodes modular form coefficients from weight 12, while f(x) encodes zeta zero oscillations. The slight bump in the normalized version (r=0.099) is not significant (p=0.33).

---

## 5. Von Mangoldt / Chebyshev Connection

**Test:** Compute g(x) = psi(x) - x (Chebyshev residual). Is f(x) - g(x)/log(x) simpler than f(x)?

### KEY FINDING: YES, dramatically.

| Statistic | f(x) | g(x)/log(x) | Residual f - g/log |
|---|---|---|---|
| Mean | 0.206 | 0.058 | 0.148 |
| Std dev | 4.233 | 4.256 | **0.387** |
| Mean |.| | 3.238 | -- | **0.337** |

- **Pearson correlation f(x) vs g(x)/log(x): r = 0.9959** (nearly perfect)
- **Standard deviation reduction: 90.87%**
- The residual f(x) - g(x)/log(x) has ~10x smaller amplitude than f(x)

### Scale analysis (improvement grows with x):

| Range | std(f) | std(residual) | Ratio |
|---|---|---|---|
| [2, 100) | 0.40 | 0.20 | 0.500 |
| [100, 1000) | 0.74 | 0.21 | 0.277 |
| [1000, 10000) | 1.80 | 0.26 | 0.146 |
| [10000, 100000) | 4.42 | 0.40 | **0.090** |

At x ~ 100000, the Chebyshev-based approximation captures **91%** of f(x)'s variance, with the residual being only 9% of the original signal.

### Sample values:

| x | f(x) | g(x)/log(x) | residual |
|---|---|---|---|
| 100 | -0.718 | -1.293 | 0.576 |
| 1000 | -0.404 | -0.481 | 0.076 |
| 10000 | 2.032 | 1.455 | 0.577 |
| 100000 | 4.538 | 4.479 | 0.059 |

### Interpretation

This is the **partial summation identity** at work:

    pi(x) = psi(x)/log(x) + integral_2^x psi(t)/(t log^2(t)) dt

So f(x) = pi(x) - R(x) and g(x)/log(x) = (psi(x) - x)/log(x) are linked by:

    f(x) - g(x)/log(x) = [R(x) correction terms] + [integral correction]

The residual f(x) - g(x)/log(x) is O(1) and grows very slowly, confirming that the dominant oscillations in both f(x) and g(x) come from the same source: the nontrivial zeta zeros. The partial summation relates them tightly.

**However:** this does not provide a computational shortcut. Computing g(x) = psi(x) - x requires summing Lambda(n) over all n <= x, which is O(x) -- worse than the O(x^{2/3}) sieve methods.

---

## Overall Conclusions

1. **Bernoulli numbers:** Irrelevant to f(x) in this range (contributions < 10^{-7}).
2. **Zeta values:** No universal algebraic relation. PSLQ finds only x-dependent spurious fits.
3. **Dirichlet L-values:** Same as zeta values -- no universal relation.
4. **Ramanujan tau:** Zero correlation with f(x). Different mathematical objects.
5. **Chebyshev psi(x):** The residual g(x)/log(x) captures ~91-99% of f(x) via partial summation. This is a known identity, not a new shortcut.

**Bottom line:** f(x) does not satisfy any algebraic relation with standard number-theoretic constants (zeta values, L-values, Bernoulli numbers, modular forms). Its oscillatory behavior is entirely governed by the nontrivial zeta zeros, consistent with the explicit formula. The Chebyshev connection is the tightest link, but it is a restatement of the same information (both encode zero oscillations), not a path to faster computation.
