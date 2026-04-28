# Session 198 — Re-verify E1.5 (h_2(pi(X)/X) entropy bound)

**Date:** 2026-04-28
**Mode:** re-verify-closure (adversarial)
**Target:** E1.5 — pi(x) mod m saturates at h_2(pi(X)/X) per step. The
only A-graded edge in the project.
**Self-grade:** **C** (closure stands; small wording refinement; no
algorithmic opening found).

## Mission

The re-verify-closure prompt asks: was a previously-closed edge
conservatively framed? Find the missed algorithmic angle. E1.5 was
the candidate I picked because:

1. It is the project's only A-graded edge.
2. The closure language ("CRT reconstruction cannot win") is broad and
   could have hidden a non-incremental shortcut.
3. The bound h_2(pi(X)/X) → 0 asymptotically *might* have constructive
   content (vanishing entropy rate) that wasn't fully extracted.

## Three candidate missed angles, evaluated

I considered three directions before settling on the joint k-moduli
verification:

**(A) Conditional-on-side-information.** Could H(pi(x) mod m | s) for
some side info s be much smaller than h_2(pi(X)/X), giving a per-step
shortcut? Direct computation: conditioning on x mod 6 ∈ {1, 5} doubles
the prime density to 2 pi(X)/X but the binary entropy h_2(2p) does
not vanish faster than h_2(p) — both scale as p log(1/p) for small p.
**Constant factor only; closure unchanged.**

**(B) Larger-step bound.** Could H(pi(x + Δ) mod m | pi(x) mod m) be
smaller than Δ · h_2 for some Δ? No: for Δ such that Δ · pi(X)/X >> 1,
the increment is approximately uniform mod m and the conditional
entropy approaches log_2(min(m, Δ)). Larger steps extract more info
*per step* but at higher per-step cost. **Trade-off, not a shortcut.**

**(C) Direct (non-incremental) computation of pi(X) mod m.** This is
NOT closed by E1.5. The bound is a chain-level entropy statement; it
constrains per-step incremental algorithms but does not bound the
algorithmic cost C(X, m) of computing pi(X) mod m for a fixed (X, m).
**This IS a real distinction the closure language conflated.** See
finding below.

## Empirical contribution: joint k-moduli conditional entropy

To make the framing distinction precise, I measured a quantity that
S12 / S69 did not: the **joint** k-moduli conditional entropy
H(joint(x) | joint(x-1)) where joint(x) = (pi(x) mod m_1, ..., pi(x)
mod m_k). Two competing hypotheses:

  - **(H1) Independence:** joint cond. entropy = sum of marginals
    = k · h_2(pi(X)/X) — would imply k-fold gain per step.
  - **(H2) Perfect correlation:** joint cond. entropy = h_2(pi(X)/X) —
    the k coordinates share a single-bit randomness source.

Result at X = 10^6 (joint state space << pi(X) for k ≤ 4):

| k | moduli | h_2 | sum marginal | JOINT | J/h_2 |
|---|--------|-----|--------------|-------|-------|
| 1 | {2} | 0.3969 | 0.3969 | 0.3969 | 1.0000 |
| 2 | {2,3} | 0.3969 | 0.7937 | 0.3969 | 1.0000 |
| 3 | {2,3,5} | 0.3969 | 1.1906 | 0.3969 | 1.0000 |
| 4 | {2,3,5,7} | 0.3969 | 1.5875 | 0.3968 | 0.9998 |

J/h_2 = 1.0000 to 10^{-4}. **(H1) decisively rejected; (H2) confirmed.**
Higher k (5, 6) shows a slight DROP below h_2 (J/h_2 ≈ 0.96 at k = 6),
which is the same finite-state coverage effect documented in S69 §3
(joint state space approaches pi(X)) — NOT (H1) signal.

**Mechanism:** trivial extension of E1.5. Since
`pi(x) mod m_i = (pi(x-1) mod m_i + 1[x prime]) mod m_i` holds
simultaneously for all i, the joint transition is determined by the
*single* bit b = 1[x prime]. Joint conditional entropy = h_2(P[x prime])
= h_2(pi(X)/X), regardless of k.

## Sharpened closure of E1.5

> **S198 sharpening of E1.5:** The per-step joint conditional entropy
> across k moduli is bounded by `h_2(pi(X)/X)` — **constant in k**, not
> linear in k. The k coordinates of the joint state share the single-bit
> prime indicator as their unique randomness source. Combining k coprime
> moduli incrementally gives **zero** new bits per step over a single
> modulus.

The previous edge wording ("scales linearly (not multiplicatively)") is
loose. The precise statement is "constant in k". Verified empirically
in 10/21 cells at the 10^{-4} precision floor; deviations in the
remaining cells are finite-state coverage artefacts (joint state space
≳ pi(X)/100), already documented in S69.

## What the adversarial probe DID NOT find

A polylog algorithmic opening. **Closure (i) "no incremental CRT
shortcut" stands and is sharpened.** Closure (ii) "no polylog pi(X)
mod m via direct computation" is *not* what E1.5 closes; it remains
the project's central open question, treated correctly by CLOSED_PATHS
row 244 ("CRT is circular: computing p(n) mod q requires pi(x;q,a)
exactly, as hard as pi(x)").

The adversarial probe FAILED to find a polylog opening because:

- The chain-level entropy bound is information-theoretically tight.
  No conditioning, larger-step strategy, or multi-modulus combination
  improves the per-step rate beyond h_2(pi(X)/X).
- The closure does NOT bound non-incremental algorithms, which would
  require a different argument (a circuit-complexity or
  pseudorandomness lower bound).
- E1.5 is correctly framed as a *negative-shape* edge ("incremental
  approaches can't win"), not as a polylog closure.

## Edges composed / cited

- **E1.5** (entropy bound): refined wording. "Linear" → "constant" in k.
  Joint k-moduli statement added.
- **CLOSED_PATHS row 243** (pi(x) mod m shortcut): wording sharpened
  but verdict unchanged.
- **CLOSED_PATHS row 244** (CRT reconstruction): correctly framed as
  the "polylog pi mod q is open" boundary; not affected.
- **E1.6** (parity bisection): same single-bit randomness source
  picture; consistent.

## What this session produced that was not in the project before

1. **Direct empirical measurement of joint k-moduli conditional entropy**
   at k = 1..6 across X = 10^4..10^6. The S69 work measured per-modulus
   marginals. This session confirms the natural extension to joint
   states: ratio J/h_2 = 1.0000 to 10^{-4} in the regime
   prod(m_i) << pi(X). Reproducible script: `joint_kmoduli_entropy.py`.

2. **Sharpened wording of E1.5:** "linear in k" replaced by
   "constant in k". Edge entry edited inline; CLOSED_PATHS row 243
   unchanged (verdict correct).

3. **Explicit two-closure framing.** Distinguishes incremental-CRT
   closure (E1.5 sufficient) from non-incremental polylog pi(X) mod m
   (open). This was implicit in row 244's wording but never made
   explicit in the edge entry.

4. **Three candidate missed angles documented and ruled out.**
   Future agents reading E1.5 know which directions don't work and
   why.

## Why C-grade and not B or A

- **Not A:** no new theorem, algorithm, or polylog opening.
- **Not B:** the joint k-moduli measurement is a triviality given
  the existing E1.5 mechanism (the same single-bit increment argument
  applies). The empirical verification is a sanity check, not a
  refinement that extends the edge's scope.
- **C:** closure stands. Wording sharpened. Three candidate missed
  angles checked and documented. This is exactly what the
  re-verify-closure prompt called the C-grade outcome:
  *"closure stands, no missed angle. Document why the adversarial
  probe failed."*

## Falsifiability statement

The session output is testable:

- The joint conditional entropy measurement at any new X with
  prod(m_i) << pi(X)/100 should give J/h_2 = 1.0000 ± 10^{-4}. Tested
  in 10 cells; all PASS.
- The candidate missed angles (A, B, C above) are all falsifiable.
  (A) and (B) are falsified by closed-form computation; (C) is
  falsifiable in principle by exhibiting a polylog pi(X) mod m
  algorithm.

## Files modified

- `experiments/information_theory/joint_kmoduli_entropy/joint_kmoduli_entropy.py` — created.
- `experiments/information_theory/joint_kmoduli_entropy/joint_kmoduli_entropy_results.md` — created.
- `EDGES.md` E1.5 entry — wording sharpened ("linear" → "constant in k").
- `status/CLOSED_PATHS.md` — appended S198 row.
- `status/SESSION_INSIGHTS.md` — S198 entry appended.
- `archive/sessions/session198_reverify_e1_5.md` — this synthesis.
- `.run_state` — set to 198.

## Recommended next-action

For algorithmic progress on the OPEN question (closure (ii) above —
polylog pi(X) mod m via direct computation), pursue ATTACK_VECTORS.md
frontier targets. E1.5 itself is now reliably *not* the obstacle for
that question; it concerns only incremental algorithms, which are not
the right shape for polylog.

Specifically: a polylog subroutine for pi(X) mod 2 alone would, by
CRT with k = O(log X) calls, give polylog pi(X). E1.5 places no
constraint on the existence of such a subroutine.

## Session-end self-evaluation (CLAUDE.md 4-question)

1. **What did I produce that was not in the project before this session?**
   (a) Direct empirical measurement of joint k-moduli conditional
   entropy (vs. per-modulus marginals measured in S69); (b) sharpened
   wording of E1.5 ("constant in k" rather than "linear"); (c) explicit
   two-closure framing (incremental-CRT vs. non-incremental polylog
   pi mod m); (d) three candidate missed angles documented.

2. **What edges did my work compose or cite?**
   E1.5 (refined); CLOSED_PATHS rows 243 (sharpened wording) and 244
   (consistent); E1.6 (compatible). No cross-domain technique imported.

3. **If my session produced only duplicate closures, why?**
   This session is honestly a closure CONFIRMATION rather than
   duplicate. The adversarial probe did not find an opening. The
   joint k-moduli measurement is new but a triviality given E1.5's
   mechanism. C-grade reflects this honestly.

4. **What is the next-action for the next agent?**
   For E1.5: no further work needed; closure is now sharpened and the
   open frontier (polylog pi(X) mod m via direct computation) is
   correctly distinguished. For project-level next-action: pursue
   ATTACK_VECTORS.md frontier targets, which is where polylog
   progress would actually originate.
