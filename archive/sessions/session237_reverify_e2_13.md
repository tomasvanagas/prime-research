# Session 237 — Re-verify E2.13 closure (Gowers `U^k`(χ_P) = HL singular series)

**Mode:** `re-verify-closure` — adversarial probe of an existing closure.
**Target edge:** E2.13 (closed S85, refined S205).
**Date:** 2026-04-30.
**Self-grade:** **C** — closure stands; the adversarial finding is a
quantitative sharpening, not a re-opening.

## Closure under attack

E2.13 (Gowers `U^k`(χ_P) matches Hardy-Littlewood `{0,1}^k`-cube
singular series) was closed at S85 with the argument:

> Why this is **not** an algorithmic opening: extracting `Q^k(chi_P)`
> empirically requires summing over (x, h) ∈ (Z/N)^{k+1}, cost
> Theta(N² log N) at U². The structural identity Q^k(chi_P) = S_k
> provides only the same information that HL already gives, in
> compressed form.

S205 added a pointwise refinement: the spike approximator
`T_Q(n) := (π(N)/N) · Σ_{q sqf ≤ Q} μ(q)/φ(q) · c_q(n)` is
polylog-evaluable and reproduces the truncated HL twin-prime singular
series within 0.6 % at every (Q, h) cell.

**The adversarial question:** the closure dismisses algorithmic content
because *empirical extraction* of Q^k is expensive. But T_Q is polylog.
Does T_Q give pointwise primality discrimination BETTER than the
trivial wheel sieve W_Q at the same Q?

## Probe (this session)

`experiments/constructions/spike_primality_discriminator/`:

1. Compute T_Q(n) via the corrected Hardy-Wright Hölder identity
   `μ(q)/φ(q) · c_q(n) = μ(g)/φ(q/g)` where `g = gcd(n, q)`.
2. Compute the wheel sieve W_Q(n) = 1[gcd(n, prod_{p ≤ Q} p) = 1].
3. Score AUC of T_Q(n) and W_Q(n) for primality across n ∈ [2, 30000]
   at Q ∈ {6, 30, 210, 2310}.
4. Restrict separately to the window n > Q (the asymptotic sieving
   regime) and the full range.

## Empirical headline

```
                  Window n ∈ (Q, N]            Full range n ∈ [2, N]
  Q       AUC(T_Q)    AUC(W_Q)   diff       AUC(T_Q)    AUC(W_Q)   diff
  6       0.9111      0.9111      0.000     0.9108      0.9106     +0.0001
  30      0.9719      0.9719      0.000     0.9715      0.9703     +0.0011
  210     0.99996     1.0000     -0.0000    0.9994      0.9929     +0.0065
  2310    1.0000      1.0000      0.000     0.9998      0.9471     +0.0526
```

* In the n > Q window: \|AUC(T_Q) − AUC(W_Q)\| ≤ 4×10⁻⁵ across all four Q's.
* In the full range [2, N]: AUC gap grows with Q to 0.053 at Q = 2310.

The full-range gap is dominated by **the small primes p ≤ Q**: T_Q(p)
correctly evaluates to a high score for prime p, while W_Q(p) = 0 (the
wheel sieves out its own primes as multiples-of-themselves). The number
of small primes is O(Q / log Q) — a polylog band — so the "extra"
discrimination is irrelevant asymptotically.

In the asymptotic sieving regime n > Q, **T_Q and W_Q have identical
rank ordering**.

## Verdict — closure stands

Falsification criteria pre-stated:
* (F-A) AUC gap > 0.05 in window n > Q → A-grade reopen. **NOT TRIGGERED.**
* (F-B) \|AUC gap\| < 0.005 in window n > Q → closure confirmed.
  **PASSES** with max gap 4×10⁻⁵.
* (F-C) T_Q strictly better at large Q only → partial reopen.
  **NOT TRIGGERED.**

The closure of E2.13 holds. The polylog spike approximator T_Q's
pointwise content is **rank-AUC equivalent to the small-prime wheel
sieve W_Q** in the asymptotic regime. There is no algorithmic
exploitation of the HL singular series structure beyond what classical
sieve theory already provides.

## Sharpening of E2.13 (NOT a re-opening)

The closure language is sharpened by this probe:

> *Original closure (S85):* "no information beyond HL".
> *S205 refinement:* T_Q autocorrelation reproduces HL.
> *S237 sharpening (this session):* T_Q's POINTWISE content is rank-AUC
> equivalent to the Mertens-product wheel sieve in the n > Q window —
> AUC equality to within 4×10⁻⁵ at every tested Q ∈ {6, 30, 210, 2310}.

This is a meaningful refinement: it rules out a class of "smoothed
wheel + extra structure" exploitations that the prior closure
language did not formally exclude.

## Self-evaluation (CLAUDE.md 4 questions)

1. **Produced that was not in the project before this session:**
   - Empirical AUC-equivalence T_Q ↔ W_Q in n > Q window (4 conductors,
     N = 30000). Not previously measured.
   - Identification that the full-range AUC gap is an edge-artefact at
     n ≤ Q, not a structural advantage.
   - Sharpening of E2.13's closure language to AUC-quantification.

2. **Edges composed or cited:** E2.13 (target), C9 / S191 (T_Q
   construction), S205 (autocorrelation refinement).

3. **Was this only a duplicate closure?** No — but the result is a
   re-confirmation rather than a re-opening, so the grade is C. The
   probe was substantive: it tested a specific algorithmic claim that
   had not been falsified before.

4. **Next-action for next agent:** Test T_Q^k for k ≥ 2 at the same Q's
   — does any moment of T_Q discriminate primes better than the wheel
   in the n > Q window? If T_Q^k = T_Q (because T_Q is essentially
   binary on the wheel partition with bounded smooth corrections),
   no new edge. Single session.

## Updates this session makes

- `EDGES.md` — append S237 annotation to E2.13 with the AUC-equivalence
  result.
- `status/CLOSED_PATHS.md` — annotate row at session 85 with S237
  re-verify result.
- `status/SESSION_INSIGHTS.md` — Session 237 entry.
- `experiments/constructions/spike_primality_discriminator/` — new
  experiment with code, results, and structural finding.

## Files

* `experiments/constructions/spike_primality_discriminator/spike_primality_discriminator.py`
* `experiments/constructions/spike_primality_discriminator/spike_primality_discriminator_results.md`
* `experiments/constructions/spike_primality_discriminator/results.json`
* `experiments/constructions/spike_primality_discriminator/run.log`
* `archive/sessions/session237_reverify_e2_13.md` (this file)
