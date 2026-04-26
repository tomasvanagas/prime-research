# Session 52 — Borel-Padé × Cesàro-Fejér hybrid (FOCUS-4)

**Date:** 2026-04-25
**Mode:** normal (focused experiment from `TODO.md` FOCUS-4)
**Outcome:** strict regression; FOCUS-4 closed.

## Goal

Combine the two convergence-acceleration tricks tested separately in S45
(Borel-Padé regularisation, EDGES §E3.5) and S46 (Cesàro-Fejér window,
EDGES §E3.7). The hypothesis from `TODO.md` FOCUS-4: if the combination
beats Fejér alone (current best 80% recovery at T=300), it is a
publishable engineering result for analytic prime counting.

## Setup

- 8 x-values × 4 T-values × 3 modes (sharp / fejér / borel-fejér).
- 2000 Riemann zeros (`data/zeta_zeros_2000.txt`, max γ ≈ 2515).
- mp.dps = 30; per-zero `mpmath.ei(ρ · log x)` (branch-correct, EDGES §E3.6).
- π(x) ≈ R(x) − Σ_γ w(γ)·2·Re(li(x^ρ)).
- Borel-Padé applied to Fejér-windowed increment sequence with median over
  Padé orders (5,5)..(15,15) and selected diagonals.

## Measurements

### Recovery rate ⌊round(S)⌋ = π(x), x ∈ {1000, 2500, 5000, 10000, 25000, 50000, 75000, 100000}

| mode         | T=50 | T=100 | T=300 | T=1000 |
|--------------|-----:|------:|------:|-------:|
| sharp        | 3/8  | 4/8   | 6/8   | 4/8    |
| fejér        | 4/8  | 4/8   | 5/8   | **6/8**|
| borel-fejér  | 3/8  | 3/8   | 3/8   | 3/8    |

### Mean rounding gap |S − round(S)|

| mode         | T=50  | T=100 | T=300 | T=1000 |
|--------------|------:|------:|------:|-------:|
| sharp        | 0.207 | 0.296 | 0.191 | 0.222  |
| fejér        | 0.309 | 0.270 | 0.246 | 0.249  |
| borel-fejér  | 0.230 | 0.246 | 0.222 | 0.219  |

## Verdict

**Hybrid loses.** Recovery is flat at 3/8 (37.5%) regardless of T, vs
Fejér-alone 6/8 (75%) at T=1000. Hybrid never uniquely beats Fejér by
more than one case at any T; regressions outnumber wins at every T.

## Diagnosis

Borel-Padé locks into a T-independent asymptote. At x=10000:
hybrid returns S ∈ {1229.61, 1229.70, 1229.70, 1229.70} across
T ∈ {50, 100, 300, 1000} — identical to 4 decimals at three of four T.
Padé fits the leading envelope of the increment sequence; the Borel
integral ∫₀^∞ e^{−z} P(z)/Q(z) dz extracts a tail-completion that depends
mostly on low-order coefficients. New zeros added by larger T are
exponentially suppressed by the 1/k! factor in the Borel transform, so
they cannot move Padé's leading-order fit.

This is the same closure mechanism as S45's Borel-Padé alone:
"postprocessing of partial sums; cannot extract information not present
in the K zeros consumed." Combined with Fejér, the windowed sum already
converges to a definite Fejér-mean answer; Borel-Padé's "tail completion"
is then a strict modification of an already-converged sum, with no
theoretical justification.

## Failure-mode classification

**E (Equivalence)** — convergence-acceleration interventions cannot be
stacked. Each re-parametrises the same √x-bounded zero-sum information.

## Actions taken

1. Wrote `experiments/proposals/borel_fejer_hybrid.py` (single script,
   no v2/v3 variants) with `--quick` flag for smoke testing.
2. Wrote `experiments/proposals/borel_fejer_hybrid_results.md` with the
   table above + diagnosis + per-T win/regression breakdown.
3. Appended CLOSED_PATHS.md entry (now 695 lines).
4. Updated `EDGES.md` §E3.7 with the S52 follow-up paragraph.
5. Marked FOCUS-4 as DONE in `TODO.md` with cross-link to results.

## State of FOCUS queue after this session

- FOCUS-1 (Connes operator scaling): unchanged, still highest leverage.
- FOCUS-2 (π(x) mod q for fixed q): unchanged.
- FOCUS-3 (Liouville parity polylog): unchanged.
- **FOCUS-4 (Borel + Fejér hybrid): CLOSED, this session.**
- FOCUS-5 (AKS growing-dim MPOW alternative attacks): unchanged.
- FOCUS-6 (4th encoding of π(x)): unchanged.
- FOCUS-7 (literature watch): unchanged; next pass earliest 2026-05-02.

## Conclusion for the S43–S52 line of analytic-side acceleration

After this session, every documented combination of the standard
convergence-acceleration techniques on the truncated explicit formula
has been measured:

| Technique                              | Source | Best on x ∈ [10², 10⁵]            |
|----------------------------------------|--------|------------------------------------|
| raw partial sum (sharp truncation)     | folklore | 6/8 @ T=300 (this session)        |
| Cesàro mean of partial sums            | S45    | beat Borel only at small x         |
| Borel-Padé alone                       | S45    | mixed; wins at moderate x          |
| Cesàro-Fejér window                    | S46    | 80% @ T=300 (S46 range x ≤ 3000)   |
| Cesàro-Fejér × Borel-Padé hybrid       | S52    | strict regression; 37.5% @ all T   |
| Hermite/Gaussian/Riesz mollification   | S48    | 5–55× WORSE than unmollified       |
| Stein-method sub-Gaussian decay attack | S43    | no sub-Gaussian decay observed     |

**No combination has produced sub-√x asymptotic gain.** Cesàro-Fejér at
T=1000 remains the best-tested constant-factor improvement on the
explicit-formula approach.

## Hygiene checks

- 1 new .py + 1 new _results.md (companion verified).
- No `__pycache__` directories left behind.
- No multi-version scripts (no `_v2`, `_quick`, etc).
- Updated: TODO.md (FOCUS-4 marked DONE), CLOSED_PATHS.md (+1 entry),
  EDGES.md (§E3.7 follow-up paragraph), SESSION_INSIGHTS.md (S52 entry).
