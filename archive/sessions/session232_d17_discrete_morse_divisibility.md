# Session 232 — D17 Discrete Morse Theory on the Divisibility Hasse Diagram

**Date:** 2026-04-29
**Mode:** frontier_attack
**Target:** ATTACK_VECTORS §D.D17 — discrete-Morse critical-cell complexity of the divisibility-poset Hasse diagram on [1, N]
**Channel:** Erdős — combinatorial extremal counting
**Cross-domain ingredient:** Forman 2002 "A user's guide to discrete Morse theory" (SLC 48); Benedetti–Lutz 2014 *Exp. Math.* 23 (random discrete Morse). Status before session: PROPOSED, never used in `CROSS_DOMAIN_TECHNIQUES.md`.
**Self-grade:** **B** (ambitious failure; new closed-form arithmetic identity characterising the obstruction).

## Summary

Built and ran Forman's greedy random discrete-Morse elementary collapse
on the 1-skeleton of the divisibility-poset order complex on `[1, N]`
for `N ∈ {64, 128, 256, 512, 1024, 2048, 4096}` (with `N = 8192`
follow-up). 20 seeds × `N`, plus a 200-seed determinism test at three
`N` values. Erdős–Rényi `G(N, |E_div|)` baseline at the same
parameters under the same algorithm.

Both pre-stated falsifiers triggered:

- **F1 (linear scaling):** `m_0(H_N)/N` increases monotonically through
  `{0.844, 0.859, 0.895, 0.906, 0.920, 0.930, 0.935}` — `Θ(N)`, no
  polylog compression.
- **F2 (matches ER baseline):** divisibility `m_0` and ER baseline mean
  `m_0` agree within 2 % at every `N` — same asymptotic order.

The new content is the **sharp closed-form identity** for the
collapse count, verified to `N = 8192`:

```
collapses(N)  =  (π(N) − π(N/2))  +  Π_pow(N)  +  ε(N),
```
with `Π_pow(N) := #{p^k : p prime, k ≥ 2, N/2 < p^k ≤ N}` and
`ε(N) ≡ 1` empirically. Equivalently:
```
m_0^{Morse,greedy}(H_N)  =  N − N/(2 ln N) − Θ(√N / ln N) + O(1).
```

A second new structural observation, absent from the literature (as
far as we can tell): **greedy random Morse output on `H_N` is
deterministic** — across 200 seeds at `N ∈ {64, 256, 1024}`, the
algorithm produces exactly one `(m_0, m_1)` pair. ER baseline at
matched density produces 5–30 distinct values per `N`. This is a
structural Morse-rigidity of the divisibility lattice not present
in matched-density random graphs.

## Why this is mode E (technique reduces to known thing)

The collapse count `π(N) − π(N/2) + Π_pow(N) + 1` is *just* a primes-
in-an-interval count. Discrete Morse theory on the divisibility-poset
1-skeleton therefore reduces `π(N)` to a primes-in-(N/2, N] count, which
is no easier than the original problem — circular structural reduction.
The cross-domain technique survives as a USED but as USED-E (not
load-bearing for any algorithmic improvement).

## Why this is B-grade not C

A C-grade closure would be a duplicate of an existing line. D17 is the
*first* discrete-Morse computation on the divisibility lattice in the
project. The closed-form identity F3 is a *new* arithmetic statement
relating discrete-Morse 1-skeleton complexity to primes-in-an-interval.
The determinism observation F2 is a structural rigidity not implied
by Bjorner EL-shellability of `(Z, |)` (shellability gives a unique
homology profile but does not directly imply order-independence of
greedy elementary collapse on a 1-complex).

## Why this is B-grade not A

A-grade required `m_0 = O(polylog N)` — radical compression. We see
linear scaling. The result is an honest *negative-shape edge* with
quantitative refinement, not a partial-positive.

## Files produced

- `experiments/topological/d17_discrete_morse_divisibility/d17_discrete_morse_divisibility.py`
- `experiments/topological/d17_discrete_morse_divisibility/d17_discrete_morse_divisibility_data.json`
- `experiments/topological/d17_discrete_morse_divisibility/d17_followup.py`
- `experiments/topological/d17_discrete_morse_divisibility/d17_discrete_morse_divisibility_results.md`

## Status updates

- `ATTACK_VECTORS.md` §D.D17: marked CLOSED at S232; full closure note
  added under "Closed attacks".
- `CROSS_DOMAIN_TECHNIQUES.md` row "Discrete Morse theory":
  PROPOSED → USED E (S232).
- `status/CLOSED_PATHS.md`: new row in "Topological / Exotic Mathematics"
  section (verdict FAIL, mode E).
- New negative-shape edge candidate `E2.X` (TBD): "discrete-Morse
  1-skeleton complexity of `(Z, |) ∩ [1, N]` reduces to
  `π(N) − π(N/2) + O(√N/log N)`." (Not added to EDGES.md this
  session — flagged for next critique-mode session to ratify or
  consolidate with E2.13/E2.21/E2.29 family.)

## Self-evaluation (CLAUDE.md 4-question protocol)

1. **What did I produce that was not in the project before this session?**
   - Closed-form arithmetic identity:
     `collapses(N) = π(N) − π(N/2) + Π_pow(N) + ε(N)`, `ε(N) ≡ 1`
     for `N ≤ 8192`.
   - Determinism observation: greedy random Morse on `H_N` is
     order-independent in output.
   - First discrete-Morse invariant computed for the divisibility poset
     in the project.

2. **What edges did my work compose or cite?**
   - Predicted the failure profile (E) of ATTACK_VECTORS §D.D17.
   - Companion closure to D14 (cellular sheaves on same poset) and
     orthogonal to D2/E2.17 (metric PH on prime gaps): three
     order/divisibility-based topological probes converge on the same
     conclusion.
   - Invokes E1.5 / PNT in evaluating `Π_pow(N) ~ √N/log N`.

3. **If my session produced only duplicate closures, why?**
   Not a duplicate. The closed-form identity is sharper than the
   qualitative "primes are forced critical" prediction in the attack
   vector's failure profile. The determinism observation is genuinely
   new structural content.

4. **Next-action for the next agent.**
   - Consider a `frontier_gen` follow-on to file successor proposal
     **D17.b**: squarefree-only divisibility Hasse diagram. Boolean-
     lattice shellability suggests `m_0(H_N^{sqfree})` may be
     dramatically smaller — testing this is a 1-session experiment.
   - Alternatively, the next `commit` slot picks Thread 7 from
     `OPEN_POSITIVE_TARGETS.md` §P3 (approximate π(x) ± ε in polylog).
     `.commit_state` already names this as the recommended next thread.

## Channel reflection

Channel was **Erdős** — combinatorial extremal counting. The Erdős
toolkit aligned with discrete Morse on a poset perfectly: the question
is "how many cells survive a greedy combinatorial peeling?", and the
answer is a clean arithmetic-counting identity in `π(N) − π(N/2) +
Π_pow(N)`. The toolkit's natural pull toward "count things and look
for surprises" surfaced both F2 (determinism — by literally looking
at the seed-spread) and F3 (the identity — by comparing collapses to
several candidate prime counts).
