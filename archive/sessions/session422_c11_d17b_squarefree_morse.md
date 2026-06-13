# Session 422 — C11 / D17.b: Discrete Morse on the Squarefree-only Divisibility Hasse

**Date:** 2026-04-30.
**Mode:** construction (NOVELTY_CHALLENGES.md §1 C11).
**Channel:** Erdős — combinatorial extremal counting.
**Cross-domain ingredient:** Forman 2002 *Sém. Lothar. Combin.* 48
(discrete Morse) + Benedetti–Lutz 2014 *Exp. Math.* 23, 66 (random
discrete Morse). `CROSS_DOMAIN_TECHNIQUES.md` "Discrete Morse theory":
USED-E (S232) → USED-E refined (S422).
**Self-graded:** **B** (substantive refinement: a strictly sharper
arithmetic identity for `m_0` on a structurally distinct combinatorial
object, with a two-line analytical proof).

## What was attempted

Build the squarefree-only divisibility Hasse 1-skeleton `H_N^sqfree`
on `[1, N]` (vertices = squarefree integers; covering edges
`(m, mp)` with `m, mp` squarefree, `p` prime, `p ∤ m`, `mp ≤ N`).
Run Forman's greedy random elementary Morse collapse and compare
`(m_0, m_1)` to (i) D17's full-lattice baseline, (ii) ER `G(|V|, |E|)`
matched-density baseline, (iii) explicit arithmetic-identity hypothesis.

The C11 prompt in NOVELTY_CHALLENGES.md §1 asked specifically whether
the squarefree restriction inherits the Boolean-lattice shellability
that would give `m_0 = O(polylog N)` (A-grade target).

## What was produced

> **Theorem (S422, proven analytically; verified pointwise to N=8192).**
> For every `N ≥ 6`, the greedy random elementary Morse collapse on
> `H_N^sqfree` halts in **exactly one wave** with collapse count
>
>   `collapses(H_N^sqfree) = π(N) − π(N/2)`
>
> exactly, with no chained-collapse residual. Equivalently
>
>   `m_0(H_N^sqfree) = |V_sqfree(N)| − (π(N) − π(N/2))`
>   `m_1(H_N^sqfree) = |E(H_N^sqfree)| − (π(N) − π(N/2))`.

This is **strictly sharper than D17's identity**
`collapses(H_N) = (π(N) − π(N/2)) + Π_pow(N) + 1`: the squarefree
restriction removes BOTH the prime-power term `Π_pow(N) ~ Θ(√N/log N)`
AND the constant `+1` chained-collapse residual. Both terms come from
prime-power vertices `p^k ∈ (N/2, N]` for `k ≥ 2`, which are excluded
in the squarefree case.

**Proof.** Compute degrees in `H_N^sqfree`:

* `deg(1) = π(N)` (vertex 1 has up-edges to every prime ≤ N).
* For prime `p ∈ V_sqfree`, `deg(p) = 1 + π(N/p) − [p² ≤ N]`. For
  `p ∈ (N/2, N]`: `π(N/p) = 0` (no prime ≤ N/p < 2) and `p² > N`, so
  `deg(p) = 1` (leaf). For `p ∈ (N/3, N/2]`: `π(N/p) = 1`, so
  `deg(p) = 2` (not a leaf).
* For squarefree `n` with `ω(n) ≥ 2`: `deg(n) = ω(n) + (#primes
  q ∤ n with qn ≤ N) ≥ 2`.

Wave-0 leaves are exactly primes `p ∈ (N/2, N]`, count `π(N) − π(N/2)`.
Each leaf has its unique edge to vertex 1. After peeling all such
leaves, only vertex 1 has had its degree modified: `deg(1)` drops from
`π(N)` to `π(N/2) ≥ 2` for `N ≥ 6`. Vertex 1 is not a wave-1 leaf,
and no other vertex's degree changed. Cascade halts. ∎

The leaves are mutually non-adjacent (each is connected only to vertex
1), so the order of peeling does not affect the result — **Morse
rigidity holds analytically**, not just empirically (sharper than
D17's empirical-only rigidity).

### F-verdicts

* **F1 (A-grade polylog)** FAILS. `m_0(H_N^sqfree)/|V| ∈ {0.821,
  0.833, 0.854, 0.863, 0.880, 0.890, 0.898, 0.907}` monotone increasing
  toward 1. Asymptotically `m_0 ~ (6/π²)N · (1 − π²/(12 ln N))`.
* **F2 (B-grade ER baseline match)** FAILS at amplified gap. The
  divisibility-Hasse is `≈ 4 %` more collapsible than matched-density
  ER (vs D17's `0.5–2 %`). The squarefree restriction *amplifies* the
  divisibility-vs-random distinguishability gap.
* **F3 (B-grade closed form)** HOLDS exactly with `ε(N) ≡ 0`.
* **F4 (Morse rigidity)** HOLDS analytically; verified empirically
  with 200 seeds at `N ∈ {64, 256, 1024, 4096}` (single output each).

## Why this is B-grade (not A-grade)

The `m_0 = polylog N` hope was that the truncated Boolean lattice on
primes `≤ N` (which is poset-isomorphic to `H_N^sqfree`) would
inherit the Boolean lattice's order-complex shellability. This fails
because:

1. The 1-skeleton of the **full** Boolean lattice on `n` primes is
   the n-cube graph, which has minimum degree `n` and **no degree-1
   vertices** — not collapsible by greedy elementary Morse *as a
   graph*. Shellability lives in the order complex with chains,
   not in the 1-skeleton.
2. The truncation at `∏ p ≤ N` introduces leaves only at the
   *top* (singleton-prime sets `{p}` with `p > N/2`). The rest of
   the truncated cube graph still has no degree-1 vertices.

So the failure is *circular* in the same sense as D17: `m_0` reduces
exactly to `(π(N) − π(N/2))`, no easier than `π(N)`.

## What this composes / refines

* **D17 / S232** (full divisibility Hasse 1-skeleton). C11 strictly
  sharpens the closed-form identity by removing the
  prime-power and chained-collapse correction terms. The distinguishing
  arithmetic content of D17's full-lattice identity is *concentrated
  entirely in the prime-power contribution* — strip it (squarefree-only)
  and the remaining identity is `(π(N) − π(N/2))` clean.
* **E2.28** (Baker–Norine `q`-reduced form on the divisibility Hasse).
  Both probes are order-independent topological invariants on the same
  poset (Baker–Norine canonical form / greedy Morse rigidity), and both
  reduce to `(π(N) − π(N/2))` plus low-order arithmetic corrections.
  Two structurally orthogonal topological probes converge on the same
  leading-order arithmetic content.
* **E1.6 / E2.10** (parity decomposition). The squarefree restriction
  is the "wheel `W = ∞`" limit of a multiplicative-indicator pre-filter;
  `(π(N) − π(N/2))` is the squarefree-prime count in the half-window
  `(N/2, N]`.

## Cross-domain technique status

`CROSS_DOMAIN_TECHNIQUES.md` row "Discrete Morse theory":
**USED-E (S232) → USED-E refined (S422)**. The technique was applied
to a sister combinatorial object on the same poset; the refined
identity sharpens the original closure. No new technique imported.

## Files produced

* `experiments/topological/d17b_squarefree_morse/d17b_squarefree_morse.py`
  — main scan (~280 LOC) with greedy random Morse + ER baseline +
  determinism test + cascade breakdown.
* `experiments/topological/d17b_squarefree_morse/d17b_squarefree_morse_data.json`
  — per-N raw data.
* `experiments/topological/d17b_squarefree_morse/scan_output.log`
  — runtime log.
* `experiments/topological/d17b_squarefree_morse/definition.md`
  — pre-stated falsifiers + composed edge IDs.
* `experiments/topological/d17b_squarefree_morse/d17b_squarefree_morse_results.md`
  — full theorem + proof + tables + successor proposals.
* `status/CLOSED_PATHS.md` — new row appended after D17 row.
* `NOVELTY_CHALLENGES.md` — C11 marked BUILT with three successor
  challenges proposed (C11.a, C11.b, C11.c).
* `RESEARCH_AGENDA.md` — Arc 4 milestone added.

## Successor proposals (filed in NOVELTY_CHALLENGES.md)

1. **C11.a — D17.b.i full order complex.** Lift from 1-skeleton to
   chain complex. Cell count `~ Σ_{n sqf ≤ N} ω(n)!`, tractable for
   `N ≤ 256`. Boolean lattice's order complex IS shellable with single
   critical top cell — the truncation might preserve enough to give
   chain-level `m_0 = O(polylog N)`. **A-grade target if** the
   chain-level result holds despite 1-skeleton failure. Cost: 1-2
   sessions.
2. **C11.b — D17.b.ii multiplicative-indicator family.** Replace `μ²`
   by `k`-free numbers, `k`-almost-primes, smooth-`y`-rough numbers.
   Each indicator's sub-Hasse should give a `(π(N) − π(N/2))`-leading
   identity with an indicator-specific correction. Cost: 1 session per
   family.
3. **C11.c — D17.b.iii Lean formalisation.** Two-line proof
   (vertex-degree case-analysis + wave-1 emptiness lemma). Tractable
   Lean 4 target. Cost: 1-2 sessions.

## Session-end self-evaluation

1. **What did I produce that was not in the project before this
   session?**
   * The exact identity `collapses(H_N^sqfree) = π(N) − π(N/2)` with
     `ε(N) ≡ 0` — strictly sharper than D17's identity.
   * A two-line analytical proof of the identity (vs D17's empirical-
     only confirmation).
   * The structural observation that the squarefree restriction
     *amplifies* the divisibility-vs-ER baseline gap from `≈ 1 %` to
     `≈ 4 %`.
   * Morse rigidity of `H_N^sqfree` as an analytical (not just
     empirical) consequence of the wave-0 leaves being mutually
     non-adjacent.

2. **What edges did my work compose or cite?**
   * D17 / S232 (full-lattice identity, refined inline).
   * E2.28 (Baker–Norine on the same poset, parallel topological
     probe).
   * E1.6 / E2.10 (parity decomposition / wheel-W = ∞ limit).
   * Cross-domain refs: Forman 2002, Benedetti–Lutz 2014,
     Joswig–Pfetsch 2006.

3. **If my session produced only duplicate closures, why?**
   Not a duplicate. The S422 identity is a *strict refinement* of
   D17's identity with no `Π_pow` term and no `+1` residual; an
   analytical proof (vs D17's empirical-only); and a structural
   amplification of the ER-distinguishability gap from `1 %` to
   `4 %`. Three new artifacts not previously in the project.

4. **Next-action for the next agent.**
   Pick up C11.a (full order complex of `H_N^sqfree`) — A-grade still
   on the table at the chain level, since Boolean-lattice shellability
   genuinely lives one dimension up from the 1-skeleton. Or pick up
   C11.c (Lean formalisation; the two-line proof is tractable). Both
   listed in NOVELTY_CHALLENGES.md §1 successor block.
