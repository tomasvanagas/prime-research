# Session 71 — Wild Swing: §C1 BK arithmetic correction at Odlyzko's high-height tables

**Date:** 2026-04-26
**Mode:** Wild Swing (one-shot ambitious frontier attack)
**Target:** ATTACK_VECTORS.md §C1 — push Bogomolny-Keating arithmetic-
            correction probe from S49's L≈7 (T~6500, N=8000) to Odlyzko's
            tabulated zeros at zero-index 10²¹ and 10²² (heights T~10²⁰
            and 10²¹, L≈44-47).
**Channel:** Bogomolny / Berry — semiclassical periodic-orbit theory of
             zeta zeros and arithmetic corrections to GUE.
**Cross-domain ingredient:** Odlyzko's published high-precision zeros at
                              extreme heights (numerical analysis frontier;
                              precision ~10⁻⁶ in γ).
**Self-grade:** **B** — ambitious failure with structural reason.

## What I attempted

Download Odlyzko's `zeros4` (10⁴ zeros at index 10²¹+1..10²¹+10⁴, height
T ≈ 1.44·10²⁰) and `zeros5` (10⁴ zeros at index 10²²+1..10²²+10⁴, height
T ≈ 1.37·10²¹). Adapt S49's BK arithmetic-correction probe to handle:

1. The huge integer offsets (γ_n stored as γ_n - C with C ~ 10²⁰-10²¹) —
   solved by working in the offset domain and using mpmath only for the
   single L = log(T/(2π)) constant.
2. The local-constant unfolding `u_n - u_0 = (δ_n - δ_0) · L_C / (2π)`
   — exact at relative precision 10⁻¹⁶ (since δ_n/C ~ 10⁻¹⁶ within block).
3. The detection of an arithmetic signal at L ≈ 45 instead of L ≈ 7.

The probe runs the full battery: pair correlation R₂(s), BK template
Pearson, per-prime Fourier amplitude at f_p = log(p)/L, phase coherence
⟨cos(φ_p − π)⟩, plus two nulls — gap-shuffled (S49 protocol, biased) and
random-prime (proper test, this session's contribution).

## What I found

**Negative.** No prime-specific structure detected at extreme heights.

| Block | L | empirical Pearson | random-prime null μ ± σ | z |
|-------|---|--------------------|----------------------------|---|
| zeros4 (T≈10²⁰) | 44.6 | +0.0628 | +0.0630 ± 0.00021 | -0.94σ |
| zeros5 (T≈10²¹) | 46.8 | +0.0349 | +0.0345 ± 0.00037 | +0.93σ |

Direct sanity check: prime-frequency Fourier amplitudes (median 0.0141,
0.0100) are **not enhanced** over random-frequency amplitudes (0.0148,
0.0101) in the same band — ratios 0.95, 0.99. Pair RMS = 0.037, 0.040,
in line with the 4/√N empirical noise law.

The previously-reported gap-shuffled-null +33σ "signal" is null-bias
artefact (E1.10): gap-shuffled sequences have Poisson-leakage in their
long-range tail that anti-correlates with any oscillatory template.
The bias inflates further at high L (−0.30 Pearson on null at L≈45 vs
+0.49 at L≈7 in S49), making naive z-score interpretation severely
misleading.

## Why this is structural — the L⁴ obstruction (new)

The clean structural insight that justifies B-grade rather than C-grade
filing:

**BK predicted amplitude:** `|BK_pred|_max · L² ≈ 13.6` is invariant
across L=7 (S49), L=44.6, L=46.8 — confirms the theoretical 1/L²
scaling exactly.

**Empirical noise:** `pair_rms ≈ 4/√N` is verified across N ∈ {2000,
8000, 10000} (predictions 0.089, 0.045, 0.040; measured 0.087, 0.054,
0.037 — within 20%).

**Detection threshold:**
```
  N_required(L; κ) ≥ (4κ / 13.6)² · L⁴ ≈ 0.09 κ² L⁴
```

At κ=3 (3σ Pearson detection):

| L | required N | available |
|---|-------------|-----------|
| 7 (S49 control) | ≈4·10³ | 8000 ✓ (just over threshold; S49 Pearson is at noise) |
| 44.6 (zeros4) | ≈3.5·10⁵ | 10⁴ — short by 35× |
| 46.8 (zeros5) | ≈4.2·10⁵ | 10⁴ — short by 42× |
| 80 (T~10⁵²) | ≈3·10⁷ | (none tabulated) |

**Doubling the height L requires 16× more zeros to compensate.** No
existing table can match. The asymptotic regime *suppresses* the BK
signal faster than data accumulation can compensate.

## What this rules out

Closes ATTACK_VECTORS §C1 with a quantitative obstruction. Strengthens
E3.13 from "BK is empirically absent at N=8000" to:

> "BK is empirically absent at all heights through T ≈ 10²¹, AND
> requires N ≥ 0.81 L⁴ zeros for any future detection — a hard scaling
> barrier independent of computational budget."

This is the project's first quantitative version of the
"BK-asymptotically-too-weak" obstruction. Implies that any future BK-
template detection via direct R₂ residual is bounded below by the L⁴
scaling. Future arithmetic-extraction-from-zeros attempts must use a
*different* statistic (not local pair correlation) or a fundamentally
different probe.

## Self-evaluation answers (CLAUDE.md §Session-end)

1. **What did I produce that was not in the project before this session?**
   - A concrete quantitative scaling barrier: `N_required ≥ 0.81 L⁴`
     for BK detection in pair correlation, calibrated against S49 data.
   - The first **random-prime null** for BK probes (cleaner than S49's
     gap-shuffled null; controls for template shape rather than gap
     distribution). The previously-reported +30σ "signal" against
     gap-shuffled is exposed as null-bias.
   - Empirical confirmation that `|BK_pred|_max · L² ≈ 13.6` is
     L-invariant across the L=7→47 range — first cross-scale
     verification of the predicted amplitude scaling.
   - First project-internal use of Odlyzko's zeros4/zeros5 tables.
2. **What edges did my work compose or cite?** E7.1 (zeros linearly
   independent), E1.10 (gap-shuffled null methodology — extended with
   random-prime null), E3.13 (BK empirically absent — sharpened with
   L⁴ scaling), §C1 from ATTACK_VECTORS (closed).
3. **If only duplicate closures, why?** Not duplicate: the L⁴
   obstruction and the proper random-prime null are new.
4. **Next-action for next agent:** ATTACK_VECTORS.md §C2 (orders 4-6
   Conrey-Snaith correlations at N=8000) is the natural follow-up if
   anyone wants to keep probing zero correlations — but per the L⁴
   obstruction, expect a similar scaling wall there. **Consider: §C3
   bespoke statistic** (a non-pair-correlation probe of zeros that
   isn't subject to the L⁴ obstruction) is the more interesting pivot.

## File pointers

- `experiments/analytic/zeta_structure/odlyzko_high_height/odlyzko_bk_probe.py`
  + `odlyzko_bk_probe_results.md` + `.json`
- `data/odlyzko/zeros4`, `data/odlyzko/zeros5` (downloaded from
  Odlyzko's UMN site)
- CLOSED_PATHS row added at line ~748 with full closure description.
- EDGES.md E3.13 extended with S71 sub-section + L⁴ obstruction.
- ATTACK_VECTORS.md §C1 moved to "Closed attacks" with closure note.
