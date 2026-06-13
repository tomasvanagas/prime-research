# Session 418 — commit Thread 8 (P2) slot 1: batched π_h(x) baseline

**Date:** 2026-04-30
**Mode:** commit (Thread 8 / OPEN_POSITIVE_TARGETS §P2 — prime gap function
batched on h, slot 1 of 4)
**Prior thread:** Thread 7 (P3) closed PARTIAL_POSITIVE_CONDITIONAL at S244.
`.commit_state` at session start: `sessions_used:5_final` on prev thread,
`recommended_next_action: pick P2 (Thread 8) as next commit thread`. Per
commit-mode rules, this session advances to Thread 8 slot 1.
**Self-grade:** **B** — empirical baseline at three anchors and five
H-values; **structural classification of P2 into a clean two-regime
dichotomy** (EXACT = P1-shape negative, APPROX = P3-shape positive);
identifies precise slot-2 next-action. Quantitative claims falsifiable
at single-script reproducibility.

## Mission (slot-1, from `.commit_state` recommended_next_action)

> "Thread 8 = OPEN_POSITIVE_TARGETS.md §P2 — prime gap function pi_h(x) =
> #{p <= x : p+h prime} batched on h in {2, 4, ..., H}. Hypothesis (slot-1
> entry): EXACT pi_h(x) is at least as hard as single-x pi(x), so per-h-
> amortised-polylog-EXACT is impossible without resolving Thread 1-4
> (closed). PARTIAL-POSITIVE candidate: Hardy-Littlewood approximation
> HL_h(x) = S_h * li_2(x) is per-h polylog (S_h via shared small-prime
> sieve table) with empirical error O(sqrt(x)) — Thread-7-shape
> transposed to the gap-family. Slot 1 plan: build batched-h sieve
> baseline AND HL evaluator, measure (a) parallel-h sieve cost vs M*per-h
> cost, (b) HL approximation accuracy vs exact pi_h at x in {10^5, 10^6,
> 10^7}, (c) S_h shared-table polylog cost. Slot-1 self-grade target:
> B-shape (empirical baseline + structural classification of which
> positive-direction subclaims survive)."

All three sub-tasks completed.

## What was built

`experiments/analytic/batched_pi_h/slot1_baseline.py` — single-file
script implementing three estimators of π_h(x) for h ∈ {2, 4, ..., H}:

- **M1 BATCHED SIEVE** — sieve `[1, x+H]` once, then for each prime p
  ≤ x test all H/2 candidates p+h.
- **M2 PER-h SIEVE** — for each h, independent sieve of `[1, x+h]`
  (measured at single h=2, multiplied by M for prediction).
- **M3 HARDY-LITTLEWOOD APPROXIMATION** — `HL_h(x) = S_h · li_2(x)`,
  with `S_h = 2 C_2 Π_{p|h, p odd} (p−1)/(p−2)`, `li_2` via
  scipy.integrate.quad.

Driver runs three anchors (x = 10⁵, 10⁶, 10⁷ at H = ⌊log²x⌋ rounded
to even, capped at 200) plus an H-sweep at x=10⁶ over H ∈ {20, 50, 100,
200, 400}. Outputs `slot1_raw.tsv` and `slot1_samples.tsv`.

`experiments/analytic/batched_pi_h/slot1_baseline_results.md` — slot-1
results doc with the dichotomy table, falsification statement, edge
citations, and slot-2 proposal.

## Headline numbers

| metric                                  | x=10⁵ (M=66) | x=10⁶ (M=95) | x=10⁷ (M=100) |
|-----------------------------------------|--------------|--------------|---------------|
| π(x)                                    | 9 592        | 78 498       | 664 579       |
| T_batched (M1)                          | 0.041 s      | 0.438 s      | 3.893 s       |
| per-h amortised (M1/M)                  | 0.62 ms      | 4.61 ms      | 38.9 ms       |
| M1/M2 speedup (vs predicted M·T_single) | 6.71×        | 8.85×        | (skipped)     |
| per-h HL eval (M3)                      | 0.62 µs      | 1.23 µs      | 0.68 µs       |
| mean \|π_h − HL_h\| / √x                | 0.098        | 0.058        | 0.046         |
| max \|π_h − HL_h\| / √x                 | 0.242        | 0.175        | 0.185         |

H-sweep at x=10⁶ confirms per-h amortised reaches a floor at ~4 ms for
M ≥ 50; the M1/M2 speedup is constant in M, not growing.

## The structural dichotomy (the slot's main result)

```
                  per-h cost          empirical error
   EXACT          Θ(x / log x)        0
   APPROX (HL)    Θ(polylog x)        ≤ 0.25 √x  (worst, slot 1)
```

Cross-h amortisation in the EXACT regime gives only a one-time
sieve-sharing factor (~7-9× M-independent, exactly the P1-shape
pattern of Thread 6). The polylog-per-h positive lives entirely in
the APPROXIMATE regime, which is the **Thread-7 (P3) shape transposed
to the h-axis**. P2 is *not* a Thread-5 (cross-x correlation) shape —
the h-axis is not a "correlation-of-x" axis in the relevant sense.

## Three concrete things this session produced that weren't in the project

1. The `batched_pi_h` baseline script (parameterised over (x, H)),
   reusable by slots 2–4 of this thread.
2. The empirical batched/per-h cost-separation table at three decades,
   with the H-sweep at x=10⁶ showing the M-independence of the
   per-h amortised cost floor.
3. The structural classification of P2 → P3-shape (not P5-shape),
   redirecting the thread budget toward HL-approximation refinement
   (slot 2: KS-test characterisation; slot 3: Q-truncation tradeoff;
   slot 4: theoretical wrap as conditional named-exponent corollary
   for π_h analogous to Thread 7 Corollary B).

## Edges composed / cited

- **E1.5** (information density of π) — sample-complexity argument
  for the EXACT regime lower bound.
- **E6.x** (sieve hierarchy) — Eratosthenes batched primitive.
- **S224 Correlation Dichotomy** — the dichotomy pattern probed for
  the h-axis (NOT replicated; P2 is shape-different from P5).
- **S231 Thread 6 P1 final wrap** — the negative-shape pattern P2
  EXACT matches.
- **S240–S244 Thread 7 P3** — the positive-shape pattern P2 APPROX
  matches.

## Self-evaluation (per CLAUDE.md session-end protocol)

1. **What did I produce that wasn't in the project before?**
   The batched_pi_h script + slot-1 results.md + the structural
   dichotomy table for P2 with three-decade empirical evidence.
2. **What edges did my work compose or cite?**
   E1.5, E6.x, S224, S231, S240-S244 (cited above).
3. **If I produced only duplicate closures, why?** Not the case —
   slot 1 is the first empirical baseline of P2 in the project.
4. **Next-action for next agent:** slot 2 of Thread 8. Run multi-anchor
   N=30 paired sample of |π_h(x) − HL_h(x)|/√x at x ∈ {10⁶, 10⁷, 10⁸}
   for K representative h-values. KS-test against half-Gaussian.
   Identify whether F_GUE-like factor exists for the gap family.
   Script: `experiments/analytic/batched_pi_h/slot2_multisample.py`.
   See OPEN_POSITIVE_TARGETS.md §P2 "Slot 2 candidates" for full plan.

## Self-extension (per CLAUDE.md autonomy invariants)

This session built (not closed) a NOVELTY_CHALLENGES target. The
Thread-8 schedule already plans 3 follow-on slots; no separate
challenges proposed. New cross-domain technique (`Q-truncated
singular series` from Hardy-Littlewood-Vinogradov) is in scope; it
is already in the project's vocabulary (see
`experiments/constructions/spike_pointwise_HL_correlation/sanity_singular_series.py`)
so no `CROSS_DOMAIN_TECHNIQUES.md` update needed.

## Files touched

- created `experiments/analytic/batched_pi_h/slot1_baseline.py`
- created `experiments/analytic/batched_pi_h/slot1_baseline_results.md`
- created `experiments/analytic/batched_pi_h/slot1_raw.tsv`
- created `experiments/analytic/batched_pi_h/slot1_samples.tsv`
- updated `OPEN_POSITIVE_TARGETS.md` §P2 (marked Thread 8 ACTIVE,
  recorded slot 1 outcome)
- updated `.commit_state` (advanced from Thread 7 DONE → Thread 8
  ACTIVE, sessions_used:1, status:ACTIVE)
- updated `status/SESSION_INSIGHTS.md` (slot-1 entry)
- updated `RESEARCH_AGENDA.md` (Arc 8 entry, Thread 8 active)
- updated `status/CLOSED_PATHS.md` (P2 EXACT regime no-amortisation
  entry, citing batched-sieve as P1-shape)
- this synthesis
