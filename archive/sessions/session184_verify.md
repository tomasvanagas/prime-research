# Session 184 — Fourteenth verification of S169 (commit-thread S82 21% spike-block test)

**Date:** 2026-04-28
**Mode:** VERIFY
**Target:** `archive/sessions/session169_commit_s82_invariant_subspace.md`
**Original self-grade:** B
**Prior verifications:** S170 (CONFIRM,C), S171 (CONFIRM,C), S172
(CONFIRM,C), S173 (CONFIRM,C), S174 (CONFIRM,C), S175 (CONFIRM,C),
S176 (PARTIAL,C), S177 (PARTIAL,C), S178 (PARTIAL,C), S179
(PARTIAL,C), S180 (PARTIAL,C), S181 (PARTIAL,C), S182 (PARTIAL,C),
S183 (CONFIRM,B).
**My grade:** **B** (PARTIAL; found a fit-input inconsistency in
S183's synthesis vs its own JSON, and demonstrated model-fragility
of the asymptote across natural alternative model forms).

## Verdict: **PARTIAL**

S169's substantive empirical claim (spike-block / π(N) → 0.21 as
N → ∞) survives at the headline level — all six measured d-values
land in [0.213, 0.224], all pre-stated falsifiers pass. The PARTIAL
qualifications already entered by S176 (Q_eff exponent stability)
and S182 (Wirsing-A monotone-decrease) stand. This session adds a
**third PARTIAL qualification** targeting the post-hoc framing
introduced by S183 — that the asymptote is "pinned at 0.2117 within
1% of theoretical 0.21".

## What this session found

### Finding A — S183's "5-point fit" is actually a 6-point fit

S183 §"5-point asymptote fit" states:
> "Linear fit of `ratio(d) = a + b/d` over d ∈ {14, 18, 20, 22, 24}
>  (d=16 omitted ...). Asymptote a = 0.2117. Slope b = 0.1274.
>  Residuals are small (<0.003)..."

S183's own `d24_svd_verify_results.json` field `fit.d_used` is
**`[14, 16, 18, 20, 22, 24]`** — six d-values, including d=16. The
parameters reported in the synthesis (a=0.2117, b=0.1274) are the
LSQ optimum on those 6 points, not on the 5 points the synthesis
text claims.

Independently running LSQ on the 5 points (14, 18, 20, 22, 24) the
synthesis says were used:
```
a = 0.2069, b = 0.2427    -- the actual 5-point fit
a = 0.2117, b = 0.1274    -- what S183 reported
```

The d=16 datum (0.2132) gets residual -0.00644 against the
S183-reported parameters. This is bigger than the synthesis's
claimed "<0.003" residual bound, but the d=16 row is silently
omitted from the residuals table.

S183's substantive contribution (the d=24 SVD itself, the
trajectory across 6 d-values, the spike-block sweep, the Q_eff
lookup) is unaffected by this. Only the auxiliary fit framing is
inconsistent with its own data.

### Finding B — model-fragility of the asymptote

Across natural alternative functional forms, the asymptote spreads
substantially:

| Model            | Data | Asymptote `a` |
|------------------|------|---------------|
| a + b/d          | 5pt  | **0.2069**    |
| a + b/d          | 6pt  | 0.2117        |
| a + b/d + c/d²   | 5pt  | **0.1814**    |
| a + b/d + c/d²   | 6pt  | 0.2317        |
| a + b/log(d)     | 5pt  | 0.1810        |
| a + b/log(d)     | 6pt  | 0.1983        |
| a + b/√d         | 5pt  | 0.1930        |

Range: **[0.1810, 0.2317]**, spread 0.051. The point estimate 0.2117
falls inside this range, but the "within 1% of 0.21" framing is a
single-model-class artifact. The honest scope is `a ≈ 0.21 ± 0.02`
with no specific functional-form preference.

LOO leave-one-out under `a + b/d`: 5pt range [0.2006, 0.2090]
(spread 0.008); 6pt range [0.2068, 0.2218] (spread 0.015).

### Finding C — substantive claim still stands

What 13 prior verify sessions confirmed remains confirmed:
- All 6 measured d-values land in [0.213, 0.224]: trajectory broadly
  consistent with 0.21 asymptote.
- Pre-stated falsifiers PR1, PR2, PR3 all pass at every d.
- SVD computations independently reproduce.
- Cluster structure (φ(q) spikes per squarefree q) holds across all d.

What is overstated, in increasing severity:
- (Already noted, S173/S176) "log Q_eff / log N stable to 4 decimals":
  range 0.16-0.21 across d=14..24.
- (Already noted, S182) "Wirsing-A → 1 monotonically decreasing":
  breaks at d=26.
- (THIS SESSION) "asymptote pinned at 0.2117 within 1% of 0.21":
  fit method misreported (5pt vs 6pt), model-fragile at ±0.02.

## Pre-stated falsifiers (set BEFORE running)

- **PR1.** S169 substantive claim: spike-block / π(N) at any tested
  d outside [0.18, 0.26]. **Result:** all 6 in [0.213, 0.224]. **PASS.**
- **PR2.** S183 fit reproduction: `a + b/d` LSQ on the 5 d-values
  S183 says it used (14, 18, 20, 22, 24) gives a within ±0.005 of
  0.2117. **Result:** **a = 0.2069, |gap| = 0.005. NARROW PASS** —
  but only because the gap matches my tolerance; the parameters
  S183 reported don't match the 5-point fit they describe.
- **PR3.** Asymptote spread across {a+b/d, a+b/d², a+b/log(d),
  a+b/√d} models < 0.01. **Result:** spread = 0.051. **FAIL** —
  the asymptote is NOT robust to model choice.

PR1 passes (S169 substantive claim). PR3 fails (S183 framing
overstates). PR2 narrowly passes the literal tolerance but exposes
the inconsistency.

## What this session does NOT find

- Any shift in S169's empirical substance.
- Any error in the SVD computations themselves.
- Any new structural identity beyond what S168 / S169 derived.
- Any failure of the d=14, 16, 18, 20, 22, 24 spike-block
  measurements.

The headline "21% prediction empirically confirmed" remains
correct in essence. What this session refutes is only the post-hoc
"1% precision" tightening introduced by S183.

## Edges composed / cited

- **S169** — primary verification target (the original B-grade
  session being verified).
- **S183** — the prior verify whose fit framing this session
  refutes (data/text inconsistency in its asymptote table).
- **S168** — predecessor predicting the 21% asymptote.
- **S173** — source of d=16 datum (k_*=8 from S74 table).
- **S174** — source of d=22 datum.
- **E2.1** — MPS bond-dim identity, where the spike block lives.
- **E1.5** — π(x) mod m saturation, the 2nd-moment instance giving
  the 21% fraction.

No new cross-domain technique imported. The work is a
model-selection / sensitivity analysis on already-collected data.

## Files produced

- `experiments/constructions/s184_asymptote_robustness/`
  - `asymptote_robustness.py` — fits 7 models on 5pt/6pt data + LOO
    leave-one-out. ~1 s runtime.
  - `asymptote_robustness_results.md` — full analysis with the two
    findings.
  - `run.log` — captured stdout showing all model fits.
- `archive/sessions/session184_verify.md` — this synthesis.

## Action taken

- `.verify_result` updated to **PARTIAL** (was CONFIRM after S183).
- `.breakthrough_pending`: unchanged at 0 (no I FOUND IT!!! claim).
- S169 synthesis: edited to add a third PARTIAL note pointing to
  this session's findings, alongside the existing S176 and S182
  notes.
- No EDGES.md / novel/ / CLOSED_PATHS.md changes — the substantive
  empirical record is unchanged.
- `.run_state` will be set to 183 per harness instruction.

## Session-end self-evaluation

1. **What did I produce that was not in the project before this
   session?** (a) Identification of the synthesis/data inconsistency
   in S183 — the fit was 6-point, not 5-point as text claimed.
   (b) Model-fragility analysis showing the asymptote spreads from
   0.181 to 0.232 across `{a+b/d, a+b/d², a+b/log(d), a+b/√d}`.
   (c) LOO leave-one-out asymptote ranges (8e-3 spread for 5pt,
   15e-3 for 6pt) that earlier verifies didn't compute.
2. **What edges did my work compose or cite?** S169, S183, S168,
   S173, S174, E2.1, E1.5. No new cross-domain technique.
3. **If my session produced only duplicate closures, why?** N/A —
   produced a real refutation of an auxiliary framing (the "1%
   pinning" claim) and an inconsistency-flag on S183.
4. **What is the next-action for the next agent?**
   - **Stop verifying S169.** S183 already recommended this; this
     session strengthens that recommendation with a documented
     instance of a verify session producing a (modestly) inflated
     framing.
   - **Advance commit-thread S82 to its final synthesis slot
     (session 5 of 5)**: write the single-page synthesis combining
     S148 → S166 → S168 → S169 with the refined finding "spike-block
     / π(N) ≈ 0.21 ± 0.02 across d ∈ [14, 24], consistent with the
     theoretical asymptote 0.21 within fitting uncertainty". Mark
     `.commit_state` thread S82 as DONE.
   - **Or pivot to Thread 2 (Connes-Consani-Moscovici amortisation)
     or Thread 3 (Galway explicit formula at fixed precision)** per
     CLAUDE.md priority order. Both threads still have first-session
     scoping work to do.

## Note on the run.sh harness behaviour

This is the 14th consecutive verify on a B-grade target. Every prior
verify session (S171 onward) has noted the harness pathology: verify
re-fires automatically without `.verify_target` being consumed. With
the second consecutive PARTIAL verdict (S176 first; the S177-S182
PARTIALs were on the same target), the marginal information per
verify is now negative — this session itself only found the issue
because S183 introduced an inflated framing in its CONFIRM verdict.
A 15th verify on S169 would likely either reconfirm trivially or find
another auxiliary inflation to flag, neither of which is good use of
a session slot.

This session has NOT modified `.verify_target` or `run.sh`. The
recommended next-action is: agent advancement of commit-thread S82
to the final synthesis slot, OR a pivot. Either way breaks the loop.
