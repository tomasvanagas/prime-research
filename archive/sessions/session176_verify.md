# Session 176 — Sixth verification of S169 (commit-thread S82 21% spike-block test)

**Date:** 2026-04-28
**Mode:** VERIFY (sixth fire on the same target)
**Target:** `archive/sessions/session169_commit_s82_invariant_subspace.md`
**Original self-grade:** B
**Prior verifications:** S170 (CONFIRM,C), S171 (CONFIRM,C), S172
(CONFIRM,C), S173 (CONFIRM,C), S174 (CONFIRM,C), S175 (CONFIRM,C)
**My grade for this re-verification:** **C** (confirmed via independent
reproduction + a new adversarial probe — Q_eff sensitivity to ±1 in
k_* at the *same* d — that prior verifies did not run).

## Verdict: **PARTIAL**

S169's substantive empirical claims reproduce exactly under independent
re-implementation. However, the synthesis prose and EDGES.md text
contain **two quantitative inflations** that materially overstate the
underlying numbers:

- **(F1) "within 5% of 0.21 at d=14 already"** — false. The actual
  deviation at d=14 is **+6.47%** (spike block 0.2236 vs 0.21). The
  correct framing (used in PR2 verdict and `results.md`) is
  *"within 7% at d=14"*. The "within 5%" framing is reached only at
  d=20 (+4.66%). The synthesis (line 39), EDGES.md (line 597–598),
  SESSION_INSIGHTS.md (line 7390), and the CLOSED_PATHS row each
  repeat the inflated "within 5% at d=14" framing — **even though
  the same CLOSED_PATHS row also says "within 7% even at d=14"
  three lines later under the PR2 verdict**. Internal contradiction.

- **(F2) "stable to 4 decimals"** for log Q_eff / log N — false.
  The exact values are 0.184640, 0.184552, 0.185022 (spread
  4.7e-4); the rounded-to-4-decimals presentation (0.1846, 0.1846,
  0.1850) **disagrees in the 4th decimal**. The right framing is
  "agrees to 3 decimals after rounding (all → 0.185)", or "spread
  of 5e-4 across three d values". S173 already noted this
  framing degrades at d=16 (drops to 0.1616); my probe (below)
  shows it is also brittle to ±1 in k_* at the same d.

The B grade stands — the work is substantive empirical confirmation
plus a finite-N exponent measurement. But the **scope of the
"stable" framing should be tightened** in EDGES.md / SESSION_INSIGHTS.

## What this session adds beyond prior verifications

### NEW: Q_eff brittleness to ±1 in k_* (at the same d)

Sessions 170-175 verified the d=14, 18, 20 numbers as published, and
S173 added d=16 to show non-monotonicity in d. None of them probed how
sensitive log Q_eff / log N is to small perturbations in the k_*
choice itself (which is ultimately set by the empirical S74
free-cumulant cutoff, not from first principles).

I recomputed the Q_eff exponent at fixed d for k_* ∈ {k_*−1, k_*,
k_*+1, k_*+2} where k_* = {5, 15, 26} are the values used in S169:

| d | k_*−1 | k_* | k_*+1 | k_*+2 |
|---|-------|------|-------|-------|
| 14 | 0.1659 (Q=5) | **0.1846** (Q=6) | 0.1846 (Q=6) | 0.2005 (Q=7) |
| 18 | 0.1846 (Q=10) | **0.1846** (Q=10) | 0.1846 (Q=10) | 0.1922 (Q=11) |
| 20 | 0.1730 (Q=11) | **0.1850** (Q=13) | 0.1904 (Q=14) | 0.1904 (Q=14) |

So the *full range* of log Q_eff / log N across reasonable
neighbouring k_* choices is **[0.1659, 0.2005]** — a 0.034 spread,
NOT 0.0005. The "stable to 4 decimals" finding only holds in the
narrow sense "the chosen k_* triple lands on the same side of an
integer-Q boundary at three d values".

This does not refute the S169 PR3 falsifier (|exp − 0.21| < 0.05) —
all six perturbation values still satisfy it. But it tightens the
correct interpretation: the matching exponent at finite N is
**0.185 ± 0.02 at the k_* perturbation level**, not "stable to four
decimals at 0.185". S168's asymptotic 0.21 is consistent with the
upper end of this band.

### NEW: bit-exact reproduction of the artefact

I re-ran `spike_block_21pct_test.py` end-to-end (~9 minutes) and diffed
the regenerated `spike_block_21pct_test_results.json` against the
saved S169 version. **Diff is empty** — the script is bit-exactly
deterministic and the published numbers are reproducible.

### NEW: framing-inflation audit

The "within 5% at d=14" claim appears 4 places in S169 outputs and
contradicts both PR2's own verdict and `results.md`:

| Location | "within 5%" or "within 7%"? |
|---|---|
| `session169_*.md` line 39 | **5%** (inflated) |
| `session169_*.md` line 65 (PR2 verdict) | 7% (correct) |
| `EDGES.md` S169 paragraph | **5%** (inflated) |
| `SESSION_INSIGHTS.md` line 7390 | **5%** (inflated) |
| `CLOSED_PATHS.md` S169 row "Net new content (a)" | **5%** (inflated) |
| `CLOSED_PATHS.md` S169 row "PR2" verdict | 7% (correct) |
| `results.md` PR2 verdict | 7% (correct) |
| Empirical reality | **6.47%** (so "within 7%" is correct, "within 5%" is not) |

This is a small inflation (1.5 percentage points) but the pattern is
diagnostic: the headline restatement is more optimistic than the
underlying analytical statement. Future synthesis-writing should
prefer the conservative framing from results.md.

## Reproduced from prior verify sessions

- Spike block sums at d=14, 18, 20: reproduced exactly from saved S82
  JSONs (sigmas[1:1+k_*] squared and summed).
- Analytic cum(Q*) at d=14, 16, 18, 20, 22, 24: reproduced exactly via
  bit-identical script rerun.
- Q_eff = 6, 10, 13 at d=14, 18, 20: reproduced.
- Spike-block / π(N) = 0.2236, 0.2212, 0.2198 at d=14, 18, 20:
  reproduced.
- Analytic ratio cum(Q*)/(0.21·π(N)) trajectory 1.330, 1.266, 1.260,
  1.193, 1.172, 1.167: reproduced.

## What survives untouched

- The 21% spike-block-fraction asymptotic prediction is empirically
  confirmed (within 7% at d=14, decreasing to 4.66% at d=20).
- The Q_eff exponent at finite N is meaningfully near 0.185 (with
  the brittleness caveat above).
- The "missing-spike" / negative-leakage observation is real
  (spike block < cum(Q=N^0.21) by 12-20%).
- All three pre-stated falsifiers (PR1 PARTIAL, PR2 PASS, PR3
  CORRECTED) survive at every probed d, including the new
  k_*-perturbation probe here.

## What this session does NOT find

- No counter-example to the substantive 21% claim.
- No d-extrapolation that breaks the asymptotic-to-1 trajectory.
- No discrepancy between the Fourier sieve in the script and an
  independent re-implementation.

## Action items (small)

The synthesis grade B stands. Two prose corrections would tighten
the artefact:

1. In `EDGES.md`'s S169 paragraph, replace "within 5% of 0.21 already
   at d=14" with "within 7% at d=14, decreasing to 5% at d=20".
2. In `SESSION_INSIGHTS.md` and the CLOSED_PATHS row, replace
   "stable across d=14, 18, 20 (values 0.1846, 0.1846, 0.1850)"
   with "spread 5e-4 across d=14, 18, 20 (all rounding to 0.185)";
   note brittleness to ±1 in k_* (range 0.166-0.201) per S176.

These are housekeeping, not substantive. I will apply them inline
under "Edits applied below" so the catalogue is self-consistent.

## Edits applied below

I am applying the corrections to EDGES.md, SESSION_INSIGHTS.md,
and the CLOSED_PATHS row to bring the headline framings into
agreement with the underlying data. I am NOT editing
session169's synthesis itself — it is a historical record of what
the original session claimed and the verification trail (this file
+ session170-175) is the canonical correction.

## Verdict summary

- **Verdict: PARTIAL.** The substantive empirical claims (21% fraction,
  Q_eff ≈ 0.185, "missing-spike" effect) reproduce exactly. The
  synthesis's two headline framings ("within 5% at d=14",
  "stable to 4 decimals") are quantitatively overstated; the correct
  framings ("within 7% at d=14", "spread 5e-4 / agrees to 3 decimals
  after rounding") are present in the same artefacts but were
  inflated in headline form.
- **B grade stands** — substance is real, framing inflation is
  cosmetic.
- **Future verifications should test** asymptotic convergence by
  extending d ≥ 26 if compute permits, and probe whether
  `k_*` admits a from-first-principles definition (currently
  empirical from S74).

## Self-grade: **C**

Confirmed an empirical-refinement (B-grade) claim by independent
reproduction (bit-identical) plus one new adversarial probe (Q_eff
sensitivity to ±1 in k_*) that prior verifies did not run. No
substantive refutation; the framing inflations are real but minor.
The verdict is PARTIAL rather than CONFIRM because the verification
*did* surface specific quantitative claims in the synthesis that
fail when read literally.
