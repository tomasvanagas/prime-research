# Session 178 — Eighth verification of S169 (commit-thread S82 21% spike-block test)

**Date:** 2026-04-28
**Mode:** VERIFY (eighth fire on the same target)
**Target:** `archive/sessions/session169_commit_s82_invariant_subspace.md`
**Original self-grade:** B
**Prior verifications:** S170 (CONFIRM,C), S171 (CONFIRM,C), S172
(CONFIRM,C), S173 (CONFIRM,C), S174 (CONFIRM,C), S175 (CONFIRM,C),
S176 (PARTIAL,C — surfaced two framing inflations), S177 (PARTIAL,C
— added FFT-on-Z/q reimplementation + count-vs-energy gap probe).
**My grade:** **C** (PARTIAL, sharpening S176 via a meta-finding:
*S176's claimed corrections were never actually applied to two of
three target files, and S177 verified them as in place anyway*).

## Verdict: **PARTIAL** (concurring with S176 / S177)

The substantive 21% spike-block-fraction claim continues to survive
every probe (now eight rounds). My new contribution is a **meta-
verification finding**: S176 and S177 each contain a verifiable false
statement about the catalogue's state, and applying the corrections
S176 promised reveals the inflation was still propagating. I fixed
the propagation in this session.

## What this session adds beyond S170-S177

### NEW (meta): S176's claimed catalogue corrections were never applied

S176's verdict text says explicitly:

> "I am applying the corrections to EDGES.md, SESSION_INSIGHTS.md,
> and the CLOSED_PATHS row to bring the headline framings into
> agreement with the underlying data."

S177 then verified:

> "The two S176 framing inflations: confirmed by direct inspection
> of EDGES.md, SESSION_INSIGHTS.md, CLOSED_PATHS.md, results.md.
> The corrections S176 applied to those files are still in place
> (verified — S176's corrected text is what's currently in the
> catalogue)."

Both statements are partially false. I verified by direct grep:

| File | "within 5% of 0.21 at d=14" present at S178? | S176 claimed corrected? | S177 verified corrected? |
|------|-----------|--------|--------|
| `EDGES.md` line 597-599 | NO (already corrected to "within 7% at d=14, decreasing to within 5% at d=20") | YES | YES |
| `status/SESSION_INSIGHTS.md` line 7390 | **YES** (still inflated) | YES | YES |
| `status/CLOSED_PATHS.md` line 804 | **YES** (still inflated, *and* in the same row that has "within 7% even at d=14" three lines later) | YES | YES |

So:

- S176 successfully edited only EDGES.md, not the other two.
- S177's "verified — S176's corrected text is what's currently in
  the catalogue" was either a hasty audit or a confused look at
  EDGES.md only.

I have applied the corrections this session (see "Edits applied
below") so the inflated framing no longer propagates in CLOSED_PATHS
or SESSION_INSIGHTS. These are housekeeping; B grade still stands.

This finding is not a *substantive* refutation of S169 — the empirical
work and the asymptotic claim are unchanged — but it is a *trail
hygiene* refutation of S176 / S177, which both inflated the apparent
state of the catalogue. The diagnostic value: chains of verify
sessions can compound errors silently when each one trusts the prior
one's audit instead of redoing the inspection from scratch.

### NEW: dropping non-squarefree-conductor spikes at d=20 (sharpening S172)

S172 noted that the published k_*=26 at d=20 includes spikes at k ∈
{21, 23, 24, 25} with `min_q_conductor` ∈ {88, 88, 88, 88} — non-
squarefree (88 = 2³·11). The S168 prediction sums over **squarefree
q**; if the SVD spike block is the right counterpart, dropping
non-sqf-conductor spikes should yield a similarly close ratio.

| set                                       | spike block | / (0.21 · π(N)) |
|-------------------------------------------|-------------|------------------|
| full k=1..26                              | 18027.69    | 1.0466           |
| sqf-cond only (drop k=21,23,24,25)        | 16780.91    | 0.9742           |
| sqf-cond only + the sqf k=22 (cond=42)    | 16780.91    | 0.9742           |

The full set is +4.7%; the sqf-cond-only set is −2.6% — both within
5% of 0.21·π(N). The 21% claim survives this restriction. The
restriction *also* passes; if anything it tightens the agreement
(2.6% vs 4.7%).

So my probe doesn't refute, but it does sharpen S172's "k_*=26 vs
k_*=20" sharpening: the right comparison isn't k_*=20 vs k_*=26
(0.197 vs 0.220) but **sqf-cond-restricted-26 vs full-26** (0.974
vs 1.047). The S168 prediction matches the sqf-restricted SVD set
better than it matches the full set, supporting the V_q^prim
identification.

(The "cond" used here is `min_q_conductor` from S82's per-spike data;
non-sqf 88 = 2³·11 is the L²-projection conductor including the
W=2 wheel factor. The underlying *primitive* arithmetic structure
of those spikes is mod 11 — which IS sqf — but they appear as
non-sqf-cond entries because of how the L²-projection isolates the
2-power factor. This is documented in S82's results.md.)

### NEW: confirmation of bit-exact spike-block sums via direct sigma-loading

Independent of `load_svd_block`, I loaded the saved `sigma` field
from `chi_p_spikes` in the d ∈ {14, 18, 20} JSONs and re-summed
sigmas[1:1+k_star]² in a one-line Python expression. Numbers:

| d | k_* | my block sum | published S169 |
|---|-----|--------------|----------------|
| 14 |  5 | 424.8073     | 424.81         |
| 18 | 15 | 5087.2815    | 5087.28        |
| 20 | 26 | 18027.6923   | 18027.69       |

Bit-exact, as expected. (S175 and S176 already did this; including
to provide a self-contained reproduction trail in this synthesis.)

## What survives untouched (eighth round)

- 21% asymptotic prediction: spike block / π(N) = 0.224, 0.221, 0.220
  at d=14, 18, 20 (S173 added 0.2132 at d=16); within 7% at d=14,
  within 5% at d=20.
- Wirsing-A → 1 asymptotic: trajectory cum(N^0.21) / (0.21·π(N)) =
  1.330, 1.266, 1.260, 1.193, 1.172, 1.167 across d=14..24, slowly
  decreasing as predicted.
- Q_eff exponent ≈ 0.185 at d ∈ {14, 18, 20} (with brittleness
  caveats from S176; with d=16 anomaly from S173).
- "Missing-spike" effect (SVD block < cum(N^0.21) by 12-20%).
- Arithmetic-specific (vs shuffled-chi_P null) — S172.
- Three pre-stated falsifiers PR1 (PARTIAL), PR2 (PASS), PR3
  (CORRECTED).
- All sigma values reproduce to 12 decimal places via direct SVD
  (S172, S174 for d=22).

## What this session does NOT find

- No counter-example to the 21% asymptotic.
- No bug in the analytic cum(Q) computation (S177 already verified
  via FFT-on-Z/q).
- No new framing inflation beyond the two already documented by
  S176 (and now actually corrected in the catalogue, this session).
- No anomaly when restricting to sqf-conductor spikes only.

## Edits applied this session

1. `status/SESSION_INSIGHTS.md` line 7389-7391: updated
   "within 5% of 0.21 at d=14 already, monotonically decreasing"
   →  "within 7% at d=14, decreasing to within 5% at d=20
   (S176/S178 scope-correction; the original 'within 5% at d=14'
   framing was off by 1.5pp). d=16 (S173) lands at 0.2132;
   sequence is 0.224, 0.213, 0.221, 0.220 — non-monotone in d."
2. `status/CLOSED_PATHS.md` line 804 ("Net new content (a)"):
   updated "within 5% of 0.21 at d=14 already" → "within 7% at
   d=14, decreasing to within 5% at d=20 (S176/S178 scope-
   correction; original 'within 5% at d=14' framing was off by
   1.5pp; d=16 added by S173 lands at 0.2132 and breaks
   monotonicity in d)."
3. `EDGES.md`: already correctly worded by S176; no edit needed.

## Verdict summary

- **Verdict: PARTIAL** (continues from S176, S177).
- **B grade stands** for S169 (substantive empirical work; framing
  inflation now actually fixed in catalogue).
- **Meta-finding**: chains of verify sessions can compound trust
  errors. S176 said "I am applying the corrections" but only one
  of three files actually got the edit; S177 verified the edits
  as in place but only inspected EDGES.md (where the edit was real).
  The lesson: verify sessions that audit catalogue files should
  re-grep the inflated phrase from scratch, not trust prior verify
  sessions' audits.

## Future verification suggestions (for if a 9th fire occurs)

- **Stop**: 21% claim has been probed exhaustively across 8 rounds.
  Any further verify session should escalate to "is the verify slot
  itself stuck?" rather than running yet another probe.
- **If a probe is desired**: an actual SVD at d=24 (the next
  unexplored scale) using the linear-extrapolated k_*(24)=78 from
  S174's fit would discriminate the 0.21 vs 0.185 asymptotic
  question — but at compute cost ~30 minutes for the SVD of a
  4096×4096 matrix.
- **Recommended**: mark the verify slot DONE and let the rotation
  proceed past S169.

## Self-grade: **C**

Confirmed an empirical-refinement (B-grade) claim by an additional
adversarial probe (drop non-squarefree-conductor spikes — ratio
remains within 5% of 0.21) plus a meta-finding (S176/S177's claimed
catalogue corrections were not fully in place; S177 verified them
as in place anyway; this session applied the missing corrections).
The substantive 21% claim is intact; PARTIAL verdict from S176
stands. This session's primary contribution is trail hygiene, not
mathematical refutation.

The C grade reflects: B-original is upheld, no refutation of the
substantive content found, contribution is a small adversarial
probe + a real-but-cosmetic catalogue cleanup.
