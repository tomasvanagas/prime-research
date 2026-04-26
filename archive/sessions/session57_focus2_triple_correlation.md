# Session 57 — Focus Task #2 extension: 3-point correlation of zeta zeros

**Date:** 2026-04-26
**Mode:** deep focus, task #2 ("Zeta Zero Structural Patterns")
**Verdict:** Match at order 3. Task 2 remains CLOSED; closure strengthened.

## Context

FOCUS_QUEUE Task #2 was marked COMPLETED in S25 (`SESSION_25_SUMMARY.md`)
with the caveat "structure might exist at scales requiring >1000 zeros to
detect." S45 (`zeta_structure_n2000_results.md`) extended five tests to
N=2000 and resolved the caveat negatively. Both extensions stayed at
correlation order ≤ 2. The residual hypothesis — "could agree with GUE
pair correlation while concealing higher-order cluster structure" — had
not been individually closed.

## Single experiment

`experiments/analytic/zeta_structure/triple_correlation.py`

### Method

1. Load N=2000 non-trivial zeros (height range 14.13 .. 2515.29).
2. Unfold via Riemann-von Mangoldt smooth counting:
   `u_n = (γ_n/2π) log(γ_n/2πe) + 7/8`. Empirical mean spacing 1.0000.
3. Slide a reference zero through u; for each i with full L_max=5
   window ahead, bin all pairs (u_j − u_i, u_k − u_i) on a 2D grid.
   1995 reference zeros usable.
4. Compare to the GUE sine-kernel determinant
   `ρ_3(s_1, s_2) = 1 − K(s_1)² − K(s_2−s_1)² − K(s_2)² + 2 K(s_1) K(s_2−s_1) K(s_2)`
   with `K(t) = sin(πt)/(πt)`.
5. Independent rigidity test: 3rd central moment of zero count in
   disjoint windows of length L ∈ {1,2,4,8,16,32}.

## Results

| Test | Value | Reference | Comment |
|------|-------|-----------|---------|
| Bulk RMS(R_3 − ρ_3), full grid | 0.0875 | 0.0864 (pair, N=2000) | matches noise floor |
| Bulk RMS, s_1, s_2 ≥ 0.5 | 0.0924 | — | away from level-repulsion edge |
| Diagonal RMS, s_2 = 2 s_1 | 0.0972 | — | equally-spaced triples |
| Anti-diagonal RMS, s_2 = s_1+1 | 0.0884 | — | constant 2nd-to-3rd gap |
| c_3 in L=1 disjoint windows | −0.000 | Poisson 1 | GUE rigidity |
| c_3 in L=8 disjoint windows | −0.001 | Poisson 8 | factor 10⁴ separation |
| c_3 in L=32 disjoint windows | +0.000 | Poisson 32 | factor 10⁴ separation |
| Var, L=32 | 0.459 | Poisson 32 / GUE ~0.4 (log L scaling) | matches GUE |

## Interpretation

Three independent third-order tests (bulk grid, two slices, count
cumulants) all return GUE-consistent values. The 3-point correlation
agrees with the determinantal sine-kernel prediction to within the
same RMS noise as the pair correlation at the same sample size,
including on the equally-spaced and constant-second-gap slices.
Cumulant rigidity is the strongest single signal: c_3 stays at ~10⁻³
across L while a Poisson process would give c_3 ranging from 1 to 32,
a factor-10⁴ separation entirely consistent with GUE-class rigidity.

The structural battery on the Riemann zeros now covers:

| Order | Test | Verdict |
|-------|------|---------|
| 1 | mean density, equidistribution mod constants (S25) | GUE / uniform |
| 2 | pair correlation, number variance, DFT (S25 + S45 to N=2000) | GUE |
| 2 | linear (PSLQ), early × late block PSLQ (S25 + S45) | no relations |
| 2 | recurrence in partial sums (S25) | none |
| **3** | **3-point correlation R_3 (S57)** | **GUE** |
| **3** | **count cumulant κ_3(L) (S57)** | **GUE rigidity** |

## Implication

Eliminates the last individually-untested cell of the structural battery.
Task 2 was already CLOSED; this entry strengthens the closure rather
than re-opening anything. No change to FOCUS-3..7 priorities.

The result is also consistent with — but does not test — the broader
Bogomolny-Keating / Conrey-Snaith conjectures for higher-order
correlations of Riemann zeros. The agreement at order 3 with the bare
sine-kernel determinant (as opposed to a sine-kernel-plus-arithmetic
correction) is what would be expected at this scale (γ_n up to 2515,
unfolded scale 1995); arithmetic corrections become visible only at
large lag in the Fourier-conjugate variable, not here.

## Files

```
experiments/analytic/zeta_structure/triple_correlation.py
experiments/analytic/zeta_structure/triple_correlation_results.md
experiments/analytic/zeta_structure/SESSION_25_SUMMARY.md (S57 update block appended)
status/CLOSED_PATHS.md (one entry appended)
status/SESSION_INSIGHTS.md (Session 57 entry appended)
archive/sessions/session57_focus2_triple_correlation.md (this file)
```

## State of project

No breakthrough. Task 2 still CLOSED. The polylog gap remains:
exact O(x^{2/3}) (`algorithms/v10_c_accelerated.py`) vs O(polylog) ~50%
digits (`R^{-1}(n)`). No FOCUS queue items advanced; the project
remains in steady state.
