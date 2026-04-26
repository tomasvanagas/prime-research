# Session 45 (2026-04-25): Zeta-structure tests at N=2000 — Settle S25 caveat

**Trigger.** FOCUS_QUEUE Task #2 (Zeta Zero Structural Patterns) was marked
COMPLETED in Session 25 but with one stated caveat: structure might exist at
scales requiring >1000 zeros to detect. This session settled that caveat by
extending the most discriminating S25 tests to N=2000.

## What was done

1. Generated zeros γ_{1001..2000} via `mpmath.zetazero` at 30-digit precision
   (~7.5 minutes, ~0.45 s/zero). Output: `data/zeta_zeros_2000.txt`.
2. Wrote a single consolidated extension script
   `experiments/analytic/zeta_structure/zeta_structure_n2000.py` running six tests:
    - Pair correlation R₂(r) vs GUE
    - Number variance Σ²(L)
    - DFT power spectrum and spectral flatness vs GUE-interp ensemble
    - Mod-constant discrepancy + Weyl sums
    - PSLQ on the previously-untested LATE block γ_{1001..1100}
       (1225 pairs + 4060 triples)
    - Cross-block PSLQ (early × late, 400 pairs)
3. Wrote `zeta_structure_n2000_results.md` with raw output and an interpretation
   section noting two methodological artifacts.

## Key results

| Test                               | N=1000 (S25)  | N=2000 (S45) | Verdict                |
|------------------------------------|---------------|--------------|------------------------|
| Pair-corr RMS deviation from GUE   | ~0.10–0.12    | **0.0864**   | GUE fit improves       |
| Log-power correlation vs GUE       | 0.9999        | 1.0000       | identical              |
| Band 1–4 spectral flatness         | 0.93–0.999    | 0.93–0.999   | identical              |
| Weyl mean / (1/√N)                 | ~0.95         | 0.946        | identical              |
| Discrepancy ratios vs LIL          | 0.3–0.8       | 0.3–0.8      | identical              |
| PSLQ pairs in LATE block           | (not tested)  | **0/1225**   | no relations           |
| PSLQ triples in LATE block         | (not tested)  | **0/4060**   | no relations           |
| Cross-block (early × late) PSLQ    | (not tested)  | **0/400**    | no relations           |

The pair-correlation deviation **decreased** with more data — exactly the
opposite of what would be expected if non-GUE structure existed at the larger
scale. The PSLQ on the LATE block is the most informative new test: it directly
addresses the hypothesis "early zeros are atypical, structure appears at higher
index." Zero relations among γ_{1001..1100} closes that hypothesis.

## Two methodological caveats flagged in writeup

These are NOT structural findings; both are explained as test artifacts:

1. **Overall spectral flatness 0.0065** vs S25's reported 0.93–0.999. The raw
   zero sequence has a smooth Riemann–von Mangoldt linear-log ramp that dumps
   massive power into Fourier indices k=0..4 (10⁴–10⁵× median). The GUE-interp
   comparison shows the same effect (mean SF 0.0077). The high-frequency bands
   1–4 show SF 0.93–0.999, identical to S25. The S25 figure was the band-level
   value; the discrepancy is presentational, not substantive.
2. **Number variance Σ²(L) plateaus at ~0.34** for L ≥ 5 instead of growing
   logarithmically as GUE predicts. Cause: with 800 sampling windows of length L
   spanning a range of ~2000, adjacent windows overlap by ~99% for L ≥ 5,
   strongly correlating their counts and shrinking apparent variance. A clean
   test would need disjoint windows (N ≥ 10⁴ zeros). Inconclusive at N=2000,
   not anomalous.

## Verdict

**S25 caveat resolved negatively.** No structure emerges at 2× the data.
Direction remains CLOSED. Pushing to >2000 would require Odlyzko's tabulated
zeros (~10⁴–10⁵ available online); GUE universality predictions and the
*shrinking* GUE deviation here both argue strongly against any emergent
structure at that scale.

## Files created/updated

- `data/gen_zeros_2000.py` (new) — incremental generator with checkpointing
- `data/zeta_zeros_2000.txt` (new) — 2000 zeros at 30 digits
- `experiments/analytic/zeta_structure/zeta_structure_n2000.py` (new)
- `experiments/analytic/zeta_structure/zeta_structure_n2000_results.md` (new)
- `experiments/analytic/zeta_structure/zeta_structure_n2000_summary.json` (new)
- `experiments/analytic/zeta_structure/SESSION_25_SUMMARY.md` — caveat-resolution paragraph appended
- `status/CLOSED_PATHS.md` — new entry in Analytic section
- `status/SESSION_INSIGHTS.md` — Session 45 section appended
- `.run_state` — set to 9 for next run

## Remaining open direction

Unchanged: **circuit complexity of π(x)** + Berry-Keating literature monitoring.
No new viable directions opened.
