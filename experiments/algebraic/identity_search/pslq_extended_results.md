# Extended PSLQ Identity Search Results

**Date:** 2026-04-05  
**Script:** `experiments/algebraic/identity_search/pslq_extended.py`  
**Data:** `experiments/algebraic/identity_search/fx_data.npz` (x = 2..100000)  
**Precision:** mpmath mp.dps=50  
**Prior work:** `experiments/algebraic/pslq_identity_results.txt` (x up to 10000, no valid relations)

## Summary

**NO UNIVERSAL ALGEBRAIC RELATIONS FOUND.** All 15 PSLQ relations with nonzero f-coefficient and residual < 1e-10 **fail cross-validation** at other x values. They are point-specific numerical coincidences, not identities.

## Test Sections

### Section A: Standard PSLQ at large x

Basis: `[f(x), 1, log(x), sqrt(x), x^{1/3}, x^{1/4}, x^{1/5}, li(x), li(sqrt(x)), li(x^{1/3}), sin(g1*log(x)), cos(g1*log(x)), sin(g2*log(x)), cos(g2*log(x))]`

| x | coeff(f) | Residual | Cross-check | Verdict |
|---|----------|----------|-------------|---------|
| 10000 | 0 | 0.0 | N/A | SPURIOUS (coeff_f=0), found -sqrt(x)+10*x^{1/4}=0 |
| 20000 | 698 | 1.72e-35 | x=21000: res=2911 | POINT-SPECIFIC |
| 50000 | 595 | 3.81e-35 | x=51000: res=1467 | POINT-SPECIFIC |
| 100000 | 0 | 2.14e-50 | N/A | SPURIOUS (coeff_f=0), found 10-x^{1/5}=0 |

At x=10000 and x=100000, PSLQ found trivial relations among the non-f basis elements (sqrt(10000) = 10*10000^{1/4} and 100000^{1/5} = 10), ignoring f entirely. At x=20000 and x=50000, PSLQ found point-specific fits that fail catastrophically at nearby x.

### Section B: Functional relations f(ax) vs f(x)

Basis: `[f(x), f(2x), f(3x), f(4x), 1, log(x)]`

| x | Residual | Cross-check | Verdict |
|---|----------|-------------|---------|
| 1000 | 1.09e-47 | x=1500: res=9800 | POINT-SPECIFIC |
| 2000 | 5.47e-48 | x=2500: res=9268 | POINT-SPECIFIC |
| 5000 | 8.76e-47 | x=5500: res=2234 | POINT-SPECIFIC |
| 10000 | 1.64e-47 | x=10500: res=16078 | POINT-SPECIFIC |

Extended basis (adding sqrt(x), 1/log(x)):

| x | Residual | Cross-check | Verdict |
|---|----------|-------------|---------|
| 5000 | 0.0 | x=6000: res=53172 | POINT-SPECIFIC |
| 10000 | 0.0 | N/A | SPURIOUS (coeff_f=0) |
| 20000 | 0.0 | x=21000: res=33157 | POINT-SPECIFIC |

**No functional relation f(ax) = g(f(x)) exists** in this basis.

### Section C: Shifted relations f(x), f(x+1), ..., f(x+k)

Long window (k=10), basis: `[f(x), f(x+1), ..., f(x+10), 1]`

| x | Residual | Cross-check | Verdict |
|---|----------|-------------|---------|
| 1000 | 0.0 | x=3000: res=13.4 | POINT-SPECIFIC |
| 5000 | 0.0 | x=7000: res=808 | POINT-SPECIFIC |
| 10000 | 1.37e-48 | x=12000: res=123 | POINT-SPECIFIC |

Short window (k=3), basis: `[f(x), f(x+1), f(x+2), f(x+3), 1]`

| x | Residual | Cross-check | Verdict |
|---|----------|-------------|---------|
| 1000 | 0.0 | x=4000: res=10088 | POINT-SPECIFIC |
| 5000 | 2.19e-47 | x=8000: res=52651 | POINT-SPECIFIC |
| 10000 | 1.15e-46 | x=13000: res=391 | POINT-SPECIFIC |
| 50000 | 0.0 | x=53000: res=9405 | POINT-SPECIFIC |

**No shift recurrence** of the form sum_k a_k f(x+k) = c exists.

## Overall Statistics

| Metric | Value |
|--------|-------|
| Total PSLQ tests | 18 |
| Relations with coeff_f != 0 and residual < 1e-10 | 15 |
| Relations surviving cross-validation | **0** |
| Spurious relations (coeff_f = 0) | 3 |

## Interpretation

PSLQ readily finds integer relations at any single point because it is solving a single linear equation with more unknowns than constraints. The 14-element basis (Section A) gives 13 degrees of freedom to fit one value of f(x), making point-specific fits inevitable.

The critical test is **cross-validation**: does the same relation hold at other x values? Every relation found in this experiment fails cross-validation with residuals ranging from 13 to 53,000 -- orders of magnitude larger than the function values themselves (~O(1) to O(10)).

This extends the prior result from x<=10000 to x=100000 and confirms:

1. **f(x) = pi(x) - R(x) satisfies no algebraic identity** in the tested bases
2. **No functional equation** relates f(ax) to f(x) for a=2,3,4
3. **No shift recurrence** relates f(x), f(x+1), ..., f(x+k) for k up to 10
4. The oscillatory correction f(x) behaves as **algebraically independent** of standard number-theoretic functions at each point

This is consistent with f(x) encoding ~10^48 zeta zero phases that vary pseudo-randomly with x, making any fixed-coefficient algebraic identity impossible.

## Path Status

**CLOSED:** Extended PSLQ identity search for f(x) = pi(x) - R(x), x up to 100000, 14-element basis with zeta oscillations, functional relations, and shift recurrences. No valid universal relations exist.
