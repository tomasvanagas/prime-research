# Session 490 — §C6 wild-swing: Pfaffian / α-DPP point-process structure of ζ zeros at order n=4

**Date:** 2026-05-01.
**Mode:** WILD SWING (single-attempt full session, permission to fail).
**Run #** 492 (per `.run_state` post-session).
**Self-graded:** **B-grade**.

## Mission

WILD_SWING harness: pick the SINGLE most ambitious open ATTACK_VECTORS
target whose A-grade probability outcomes have not yet been attempted.
Default order from prompt: §C1 → §A1 → §B1 → §A3 → §D4 → §C2.
**All six default-priority items are CLOSED already** (§C1/S71, §A1/S84
partial, §B1/S92, §A3/S79, §D4/S80, §C2/S123). Picked the next-highest-
A-grade open §C-section item: **§C6 — Pfaffian / α-determinantal point
process structure of ζ zeros at order n=4** (PROPOSED in
CROSS_DOMAIN_TECHNIQUES §3, untouched in 488 prior sessions).

## What was produced

Two new experiments:

- `experiments/analytic/zeta_structure/c6_pfaffian_alpha_dpp_n4/c6_pfaffian_alpha_dpp_n4.py`
  (main experiment: empirical R_4 on 8000 ζ zeros + matched GUE/GOE/GSE
  Monte-Carlo pools at 32400/32400/16200 eigenvalues; 96 random non-equally
  -spaced 4-tuples; α-DPP scan with tuple-bootstrap CI; per-tuple χ²
  discrimination)
- `c6_alpha_bias_control.py` (16 independent matched-size 8000-eigenvalue
  GUE-MC trials with same 96-tuple α-DPP fit — critical bias control for
  the α\* = -1.10 finding)

Two new results files:
- `c6_pfaffian_alpha_dpp_n4_results.json` + `_results.md`
- `c6_alpha_bias_control_results.json`

## Headline findings

### 1. Decisive structural rejection of GOE/GSE Pfaffian models on ζ zeros at n=4

Per-tuple χ² discrimination (96 tuples, dof=96, per-tuple SE=0.135 from
GUE-MC batches of 1200 evs):

| Model | χ²/dof | Interpretation |
|-------|--------|----------------|
| ζ vs sine-kernel det theory | **0.96** | Perfect fit |
| ζ vs GUE-MC pool | 1.24 | Within sampling noise of GUE |
| ζ vs GOE-MC pool | **1.99** | Rejects GOE Pfaffian |
| ζ vs GSE-MC pool | **3.11** | Strongly rejects GSE Pfaffian |

This is the C6 stated B-grade success criterion (b): "explicit upper bound
on the Pfaffian-vs-determinantal discrimination at order 4." First project
test of GOE/GSE Pfaffian fits to ζ-zero R_4 at order 4. Refines E7.1 from
"DPP-typical at orders 4-6" to "DPP-typical AND not Pfaffian-typical at
order 4."

### 2. Suggestive (non-decisive) α-DPP shift on ζ zeros: α\* ≈ -1.10

α-determinant scan (Vere-Jones α-permanent generalisation) over
α ∈ [-1.5, -0.5] at Δα = 0.025:

- **α\*(zeta) = -1.10**, L2\* = 0.114 vs L2(α=-1) = 0.118 (ΔL2 = 0.0045)
- Tuple-bootstrap (200 resamples) 95% CI: [-1.150, -1.049] formally excludes -1
- **Critical bias control**: 16 independent matched-size 8000-eigenvalue
  GUE-MC trials gave α\* mean -1.013, std 0.035, **0/16 trials at α\* ≤ -1.10**
- ζ at z ≈ -2.5σ from GUE-MC distribution (parametric one-sided p = 0.007;
  non-parametric Laplace p = 0.059)

The α\* = -1.10 shift is OUTSIDE the matched-size GUE-MC distribution at
~2.5σ, but **does NOT meet the 5σ A-grade bar**. This is an A-grade-shaped
partial-positive that requires Odlyzko γ ≥ 10⁶ to push past 5σ.

## Edges cited / refined

- Refines **E7.1** ("ζ zeros are sine-kernel DPP up to order 6") with two
  new clauses:
  - "...AND not GOE/GSE Pfaffian-typical at order 4 (χ²/dof 1.99 / 3.11
    rejection on 96-tuple test)"
  - "α-DPP best-fit α\* = -1.10 ± ~0.04, z ≈ 2.5σ from matched-size GUE-MC,
    suggestive non-decisive shift, sub-frame C6.b OPEN for promotion to A"
- Refines **E3.13** by adding a 3rd discrimination dimension (Pfaffian-vs-
  DPP) to GUE-typicality.

## Cross-domain ingredient

Pfaffian point processes (Borodin 2009 arXiv:0911.1153 §2.4-2.6, Soshnikov
2000 arXiv:math/0002099 §3) and α-determinantal generalisation (Vere-Jones
1997 *NZ J. Math.* 26 α-permanent). Promoted CROSS_DOMAIN_TECHNIQUES §3 row
"Pfaffian / α-determinantal point processes" PROPOSED → **USED-E**.

## Successor sub-frames (OPEN, added to ATTACK_VECTORS as C6.b/C6.c)

- **C6.b** (Odlyzko γ ≥ 10⁶ extension): re-run α-DPP fit on Odlyzko's
  large-height block. With ~10⁶ zeros, per-tuple SE drops by √125 ≈ 11×, so
  L2 noise floor goes from 0.118 to ~0.011, and the ΔL2 ≈ 0.0045 between
  α\* = -1.10 and α = -1 becomes a 5σ effect IF the shift persists at
  higher height. Either A-grade structural fact or B-grade closure of the
  α-DPP sub-frame.
- **C6.c** (denser tuple sampling at fixed N=8000): increase tuples from 96
  to 10000 with adaptive importance sampling on regions of high R_4
  curvature.

## Bugs / pitfalls noted

- Initial GUE/GOE/GSE Monte-Carlo unfolding had a normalisation bug:
  `H /= sqrt(2*N)` gives Wigner radius √2 not 1, so the standard semicircle
  unfolding code gave mean spacing 1.35 instead of 1.0. Fixed by detecting
  empirical spectral support [a, b] and rescaling to enforce mean spacing 1
  exactly post-unfolding.

## Self-evaluation per CLAUDE.md session-end checklist

1. **What did I produce that was not in the project before this session?**
   - First project measurement of GOE/GSE Pfaffian and Vere-Jones α-DPP
     point-process models on ζ-zero R_4 at order 4.
   - First explicit α-determinant deformation measurement on any object in
     the project.
   - Bias-controlled discrimination: ζ-zero α\* = -1.10 OUTSIDE matched-
     size GUE-MC distribution at ~2.5σ.
   - Two new working experiments (`c6_pfaffian_alpha_dpp_n4.py` +
     `c6_alpha_bias_control.py`) with reproducible results.

2. **What edges did my work compose or cite?**
   - Composed: E7.1 (sine-kernel DPP at order 6) + E3.13 (GUE-typicality)
     + cross-domain Pfaffian / α-DPP machinery.
   - Refined E7.1 with two new clauses (Pfaffian rejection + α-DPP shift).
   - Did NOT cite multiple ATTACK_VECTORS edges simultaneously (single-
     attack wild swing).

3. **Why is this not A-grade?**
   - The α-DPP shift at z ≈ -2.5σ is below the project's stated 5σ A-grade
     bar. Bias control gives parametric p = 0.007 (one-sided) but non-
     parametric Laplace p = 0.059 (above 0.025 threshold). Mixed verdict;
     suggestive but not decisive at this finite-N detection floor.
   - Pfaffian rejection at χ²/dof 1.99 / 3.11 is decisive but is the
     EXPECTED outcome (zeros are GUE, not GOE/GSE) — it's strong B-grade
     refinement of E7.1, not A-grade novelty.

4. **Next-action for next agent (written into ATTACK_VECTORS as C6.b):**
   - Run C6.b: download Odlyzko γ ≥ 10⁶ block (~10⁶ zeros), re-run the
     same 96-tuple α-DPP fit + bias control, decide A-grade
     promotion / B-grade closure of the α\* = -1.10 shift.

## Why this was a good wild swing

- All six default-priority §A-§F items were closed before this session, so
  C6 was a forced choice — but it *was* the next-highest-A-grade-probability
  open frontier item.
- Partial-positive shape (z ≈ -2.5σ shift on α-DPP) is characteristic of an
  ambitious attack producing informative B-grade content that opens an
  A-grade follow-up: the technique imported (Pfaffian / α-DPP) DID perform
  real work, and the failure mode (finite-N detection floor) is structural
  and addressable by C6.b.
- Falsifiable falsifier framework (F1/F2/F3 from results.md) was set
  upfront, allowing honest reporting of "F1 PARTIAL, F2 FAILED, F3 FAILED"
  without inflating the grade.

## Files produced

- `experiments/analytic/zeta_structure/c6_pfaffian_alpha_dpp_n4/c6_pfaffian_alpha_dpp_n4.py`
- `experiments/analytic/zeta_structure/c6_pfaffian_alpha_dpp_n4/c6_pfaffian_alpha_dpp_n4_results.json`
- `experiments/analytic/zeta_structure/c6_pfaffian_alpha_dpp_n4/c6_pfaffian_alpha_dpp_n4_results.md`
- `experiments/analytic/zeta_structure/c6_pfaffian_alpha_dpp_n4/c6_alpha_bias_control.py`
- `experiments/analytic/zeta_structure/c6_pfaffian_alpha_dpp_n4/c6_alpha_bias_control_results.json`
- Updated: `status/CLOSED_PATHS.md` (+1 row, ATTACK_VECTORS §C6 closure)
- Updated: `ATTACK_VECTORS.md` (§C6 marked CLOSED with outcome summary)
- Updated: `CROSS_DOMAIN_TECHNIQUES.md` (Pfaffian/α-DPP row PROPOSED → USED-E)
- Updated: `archive/sessions/session490_c6_pfaffian_alpha_dpp_n4.md` (this file)
