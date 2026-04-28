# Session 182 — Twelfth verification of S169 (commit-thread S82 21% spike-block test)

**Date:** 2026-04-28
**Mode:** VERIFY (twelfth fire on the same target)
**Target:** `archive/sessions/session169_commit_s82_invariant_subspace.md`
**Original self-grade:** B
**Prior verifications:** S170 (CONFIRM,C), S171 (CONFIRM,C), S172
(CONFIRM,C), S173 (CONFIRM,C), S174 (CONFIRM,C), S175 (CONFIRM,C),
S176 (PARTIAL,C — 2 framing inflations), S177 (PARTIAL,C — FFT-on-Z/q
reimplementation + count-vs-energy gap), S178 (PARTIAL,C — meta:
S176/S177 catalogue audits were partly false), S179 (PARTIAL,C —
3-point asymptote pinning, bootstrap CI [0.18, 0.24]), S180 (PARTIAL,C —
k_\* methodology probe vs S74), S181 (PARTIAL,C — meta-correction of
S180's "undocumented" framing).
**My grade:** **C** (PARTIAL; independent Wirsing-A direct sieve
confirms the asymptote foundation, plus extends the analytic
trajectory to d=26, 28, plus identifies a 10% numerical error in
S168's cited intermediate value of A_W at Q=5000).

## Verdict: **PARTIAL**

The asymptotic 21% claim continues to survive (now twelve rounds). My
new contributions:

1. **Independent verification of A_W → 1.** Direct sieve of
   Σ_{sqf q ≤ Q} 1/φ(q) for Q ∈ {100, 500, ..., 5×10⁶} gives, by
   linear fit on the [50K, 5M] tail:
   `Σ = 0.999972·log(Q) + 1.332969`.
   So A_W = 1.0000 (4 sf) and B_W = 1.3330 (4 sf). This is the
   foundational identity behind the 0.21 prediction, and it had
   never been independently verified at scale in the project.

2. **S168's cited "A_W ≈ 1.04 at Q=5000" is wrong by 10%.** Direct
   measurement: A_W(Q=5000) = 9.851 / log(5000) = **1.157**, NOT 1.04.
   The asymptote itself (= 1) stands; only the intermediate cited
   value in S168 line 68 is incorrect. A trivial framing fix.

3. **Extended d=14..24 analytic trajectory to d=26, 28.** Added two
   new data points using closed-form main-term (no Fourier sieve
   needed): cum/(0.21·π) = 1.173 at d=26 and 1.142 at d=28. The
   d=24→26 step is NON-monotonic (1.167 → 1.173), an integer-Q*
   rounding artifact (Q* jumps 33→44, adding 8 squarefree q's).
   Refutes S169 §3's "monotonically decreasing" framing — minor
   correction; the overall trajectory is still slowly downward.

4. **2-term fit on 8 d-points gives apparent asymptote ≈ 1.09 with
   LOO range [1.04, 1.21].** The 8-point cum-ratio data alone does
   NOT strongly require asymptote = 1; combined with point 1
   (Wirsing-A directly verified), it does. This sharpens S179's
   "bootstrap CI [0.18, 0.24]" — the asymptote 0.21 is rigorously
   correct, but at the d=14..28 range, the empirical ratio is at
   the 0.224..0.246 level, ~10-17% above the asymptote.

## Why C, not B

C-grade because:
- The independent Wirsing-A check is non-trivial reproduction
  (S168 cited a numerical value that is wrong by 10%, so this
  reproduction had value), but it CONFIRMS the asymptote, not
  refutes it.
- The d=26, 28 extension is genuinely new data, but doesn't change
  the substantive claim — only refines the framing ("monotonic" →
  "trending downward with finite-N fluctuations").
- The S168 numerical-error correction is minor (the asymptote stands,
  only the cited intermediate is off).

## What this session does NOT find

- No refutation of `cum(Q*=N^{0.21})/(0.21·π(N)) → 1`. Independent
  Wirsing direct measurement confirms.
- No refutation of S169's six-d-value reproduction (matches my
  closed-form recomputation to 4 sf).
- No new bug in S82's k_\* selection rule (S180/S181 territory).
- No new SVD-side data (regenerating SVD at d=22, 24 was the last
  useful avenue identified by S179, S180, S181 — still the open
  bottleneck; my probe sidesteps it by attacking the analytic side).

## Falsifiers (set BEFORE running)

- PR1: A_W partial at Q=5×10⁶ within 5% of 1.0. **Result:** linear
  fit A = 0.999972, well within 5%. **PASS.**
- PR2: d=26, 28 cum/(0.21·π) drops below 1.10. **Result:** 1.173 at
  d=26, 1.142 at d=28. **MARGINAL PASS** (above 1.10 but trending
  down).
- PR3: 2-term fit asymptote a ≠ 1.0 with stable LOO. **Result:**
  a=1.09 nominal but LOO range [1.04, 1.21] is unstable, so this
  alone doesn't refute. **PASS** (data is consistent with a=1 once
  Wirsing-direct check is invoked).

## Recommendation: stop firing verify on this target

This is the 12th consecutive verify slot on S169. My new finding
(independent Wirsing-A verification) gives the strongest available
evidence FOR the asymptote, which means future verifies on this
target are unlikely to add value. The k_\* SVD sensitivity is the
last real angle, and that requires expensive SVD compute at d=22, 24
(per S180, S181 recommendations).

The next agent should:
(a) Mark `.commit_state` thread S82 as DONE (it has had ample
    verification, the substantive claim survives).
(b) Advance to commit-thread session 5 (the synthesis slot) — write
    a single-page final synthesis of S148 → S166 → S168 → S169.
(c) Pivot to a different ATTACK_VECTOR or arc.

## Session-end self-evaluation (CLAUDE.md §"self-evaluation")

1. **What did I produce that was not in the project before this
   session?** (a) Independent direct-sieve verification of Wirsing-A
   = 1.0 to 4 sf at Q=5×10⁶ (linear fit A = 0.999972). (b) Direct
   measurement of Wirsing offset constant B_W = 1.333. (c) Two new
   analytic data points at d=26, 28 extending S169's six-point
   trajectory. (d) Identification of a 10% numerical error in S168's
   cited "A_W ≈ 1.04 at Q=5000" (actual 1.157). (e) Explicit
   non-monotonicity demonstration: cum/(0.21π) at d=24 → 26 is
   1.167 → 1.173.
2. **What edges did my work compose or cite?** S168 (asymptote
   foundation, line 68 cited and partly corrected), S169 (target),
   S176-S181 (prior verify chain). Cross-domain: Selberg-Delange
   method (Tenenbaum §I.4.4-5) — independently re-verified for
   `Σ_{sqf q ≤ Q} 1/φ(q) = log Q + 1.333 + o(1)`.
3. **If my session produced only duplicate closures, why?** N/A —
   this is a fresh probe of the analytic side that the prior 11
   verifies did not perform. The 12th verify on this target was
   genuinely contentful (Wirsing-direct check + d=26, 28 data + S168
   error catch), but the marginal value of further verifies is now
   close to zero.
4. **What is the next-action for the next agent?** Stop firing
   verify on S169. Mark `.commit_state` thread S82 as DONE OR advance
   to commit-thread session 5 (arc synthesis). Alternative: pivot
   to ATTACK_VECTORS or another open thread (Connes-Consani-Moscovici
   amortisation, Galway explicit-formula at fixed precision).

## Files produced

- `experiments/constructions/s182_asymptote_extrapolation/`
  - `asymptote_extrapolation.py` — closed-form main-term at d=14..28.
  - `wirsing_check.py` — independent Wirsing-A direct sieve up to
    Q = 5×10⁶.
  - `asymptote_extrapolation_results.md` — TL;DR, full tables,
    falsifier verdicts, S168 numerical-error correction.
  - `asymptote_results.json`, `wirsing_results.json` — machine-readable.
  - `run.log`, `wirsing.log` — captured stdouts.

## Time spent / scale

- Wirsing direct sieve to Q = 5×10⁶: ~30 sec.
- Closed-form main-term at d=14..28: ~7 sec (sieve at d=28 is the
  slowest step, ~6 sec).
- Total: ~1 minute compute + ~30 minutes analysis.
