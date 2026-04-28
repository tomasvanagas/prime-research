# Session 161 — D45 wild swing: Baker-Norine graph Riemann-Roch on the prime divisor

**Date:** 2026-04-28
**Mode:** WILD SWING (frontier attack, ATTACK_VECTORS.md §D45)
**Target:** chip-firing rank `r_G(D_P^N)` of the prime divisor on the
divisibility graph `Γ_N` and the Hasse cover graph `H_N`.
**Self-grade:** **B**

## Why this attack

ATTACK_VECTORS.md preferred order from session prompt: §C1, §A1, §B1,
§A3, §D4, §C2 — all closed in their main forms. Of the OPEN frontier
attacks, **§D45** (Baker-Norine RR / chip-firing, added in S160
frontier_gen) has the highest concreteness × tractability product:
a 1-2 session budget, a non-linear combinatorial invariant orthogonal
to all 41+ prior pseudorandomness measures, and an A-grade outcome
("Brill-Noether specialness witness") that would be a publishable first
on any arithmetic graph.

Channelled mathematician: **Baker** (graph divisor theory).

## Pick rationale

1. CLAUDE.md preamble: prefer A-grade attempts; ambitious failure
   rewarded as B-grade if structural.
2. D45 is fresh (added 2026-04-28 by S160 frontier_gen); no prior
   session has computed any chip-firing invariant on any arithmetic
   graph.
3. Cross-domain ingredient (Baker-Norine 2007 RR + Dhar 1990 burning)
   has well-defined algorithmic content; tractability not an issue.
4. Negative result with explicit structural mechanism would directly
   close D45 and add a negative-shape edge.

## What was built

`experiments/algebraic/baker_norine_chi_p/`:
- `baker_norine_chi_p.py` — divisibility / coprimality graph
  constructors; Dhar's burning algorithm; q-reduction (with chip-
  lending Phase 1 via BFS-path firing + superstabilization Phase 2);
  Baker-Norine Algorithm 4.1 rank computation by exhaustive multiset
  enumeration.
- `baker_norine_hasse.py` — Hasse cover graph constructor; signature
  extraction; matched-density random divisor controls.
- `baker_norine_signatures.py` — q-reduced form signature extraction
  (D'(q), support, max chip, winnable-subtraction set W(D)).
- `run_full.py` — driver producing `full_results.json`.
- `baker_norine_chi_p_results.md` — writeup with falsifiers, theorems,
  proofs.

## What this session produced (what was NOT in the project before)

1. **Theorem (E2.28-1, divisibility graph, NEW):** On
   `Γ_N := ([1, N], divisibility edges)` with sink q=1:
   `D'_P^N = (π(N) − π(N/2)) δ_1 + Σ_{p prime, p ≤ N/2} δ_p`.
   Empirically verified for N ∈ {16, 32, 64, 128, 256, 512}.
   Proof via Dhar burn analysis (one paragraph).

2. **Theorem (E2.28-2, Hasse cover graph, NEW):** On
   `H_N := ([1, N], cover edges)` with sink q=1:
   `D'_P^N = π(N) δ_1`.
   Empirically verified for N ∈ {16, 32, 64, 128, 256, 512}.
   Proof: in Hasse, neighbors of 1 are exactly primes; one round of
   firing V \\ {1} collects all prime chips at q.

3. **Generalised Hasse identity (B-grade structural fact, NEW):** for
   D supported on `{1} ∪ {primes}` only: `D'(1) = D(1) + Σ_{p prime} D(p)`.
   Holds for D_P, D_sqfree, D_μ_pos. Fails for D_λ_pos, D_Omega2 (chips
   on composites cascade through Dhar burning).

4. **Z-score asymptotic (NEW prediction):** D'(q) for D_P beats matched-
   density random divisors at z-score scaling like `√π(N) ~ √(N/log N)`
   on `H_N`. Empirical:

   | Graph        | N=32  | N=64  | N=128 |
   |--------------|-------|-------|-------|
   | Γ_N (div)    | +3.4σ | +3.2σ | +5.8σ |
   | Hasse(N)     | +8.6σ | +11.1σ| +14.6σ|

5. **Rank verdict:** `r(D_P^N) = max(0, deg(D_P) − g + 1)` for both
   graphs at N ≥ 32 — generic Riemann-Roch value, no Brill-Noether
   specialness. Failure of `r ≥ 1` occurs at the SAME bottleneck
   vertex (typically 4) for D_P and random matched divisors → graph
   topology, not arithmetic.

6. **Closure CLOSED_PATHS row + ATTACK_VECTORS closed attacks entry +
   new edge E2.28 in EDGES.md** documenting the negative-shape result:
   chip-firing graph divisor rank is determined by `(deg(D), graph G)`
   modulo Riemann-Roch + gonality and is BLIND to deeper prime
   arithmetic.

## What this session did NOT produce

- No A-grade chip-firing-derived polylog π(x) algorithm.
- No Brill-Noether-special D_P (rank does not exceed RR bound).
- No identification of D_P with a "geometric-special" Brill-Noether
  class on either graph.

## Edges cited / composed / refined

- **E2.28** (NEW edge added in this session) — chip-firing q-reduced
  form identity for the prime divisor on arithmetic graphs.
- **E2.13 / E2.14 / E2.15 / E2.16 / E2.17 / E2.20 / E2.21 / E2.24**
  (pseudorandomness measures) — the Baker-Norine rank invariant joins
  this family but contributes a NEGATIVE-SHAPE edge: rank is determined
  by Riemann-Roch + gonality, blind to arithmetic.
- **E1.10 / E3.13** (GUE / random pseudorandomness measures all match)
  — REINFORCED at the chip-firing rank level; structurally negative
  for chip-firing as pseudorandomness detector.
- **D45** (this attack vector) — **CLOSED** mode I (identification with
  graph-topology / Riemann-Roch).

## Why B-grade and not A or C

**Why not A:** the session did not produce a Brill-Noether specialness
witness or a chip-firing-derived polylog-promising algorithm. The
closed-form identities (Theorems 1, 2, 3) are clean but reduce to
graph-NEIGHBORHOOD structure (chips on N(q) flow to q), not deeper
arithmetic. The signal grows but the growth is structurally explained
by the obvious "primes occupy specific lattice positions" fact.

**Why not C:** the session committed to a D45 wild swing as instructed
(no mid-session pivot), produced 4+ artefacts not previously in the
project (closed-form identities for both graphs, generalisation,
z-score asymptotic, rank verdict), wrote 5 source files + results
markdown, added a new edge (E2.28), updated 4 status files
(CLOSED_PATHS, ATTACK_VECTORS, EDGES, CROSS_DOMAIN_TECHNIQUES) with
non-trivial content. The negative result has a CLEAN STRUCTURAL
MECHANISM (Riemann-Roch symmetry + graph-neighborhood reduction) that
constrains future chip-firing attacks.

**Why B (substantive structural finding + ambitious failure):** the
session attempted an A-grade target from ATTACK_VECTORS.md and failed
informatively — the failure mode is structural (Baker-Norine rank is
graph-topology-determined), not "I ran out of time" or "the technique
was a bad fit". The chip-firing technique IS applicable; it just does
not surface deeper arithmetic. Two clean closed-form identities are
contributed as B-grade content (Theorems 1, 2 above).

## Self-evaluation against CLAUDE.md session-end checklist

1. **What did I produce that was not in the project before?** Two
   closed-form chip-firing q-reduced identities for D_P on arithmetic
   graphs (E2.28); generalisation to other arithmetic divisors;
   z-score growth asymptotic √π(N) ~ √(N/log N); empirical rank
   verdict; closure of D45.
2. **What edges did my work compose or cite?** E2.28 (new); E2.13/14/
   15/16/17/20/21/24 (pseudorandomness measures, joined by chip-firing
   negative-shape); E1.10 / E3.13 (GUE/random reinforced); E5.3 /
   E7.10 (TC⁰ PRIMES, untouched).
3. **If my session produced only duplicate closures, why?** N/A —
   produced new edge E2.28 with two clean closed forms not previously
   in the project.
4. **Next-action for the next agent:** D45.a (chip-firing on more
   non-trivially-arithmetic graphs like prime-gap multigraph or
   prime-Cayley), D45.b (Brill-Noether locus W^r_d(G) explicit study),
   D45.c (gonality computation) — proposed in ATTACK_VECTORS.md
   "Closed attacks" section under §D.D45.

## Reproducibility

```bash
cd experiments/algebraic/baker_norine_chi_p
python3 run_full.py    # produces full_results.json + full_run.log
```

Total runtime: ~3 seconds on standard machine. Fully deterministic
(seed=42 for random divisor sampling).

## Status file updates this session

- `EDGES.md` — added new edge **E2.28** (chip-firing q-reduced form
  identity for D_P on `Γ_N` and `H_N`, EVS M shape).
- `CLOSED_PATHS.md` — added §D45 closure row (mode I, B-grade case (i)).
- `ATTACK_VECTORS.md` — marked §D45 closed; added detailed closure note
  in "Closed attacks" section.
- `CROSS_DOMAIN_TECHNIQUES.md` — promoted Baker-Norine row from
  PROPOSED to USED-I with edge link to E2.28.
- `archive/sessions/session161_d45_baker_norine_chi_p.md` — this file.
