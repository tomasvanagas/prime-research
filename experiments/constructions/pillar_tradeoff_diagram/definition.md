# C6 — Three-Pillars × HKM Time-Space Tradeoff Diagram

**Target:** `NOVELTY_CHALLENGES.md` §1.C6.
**Edges composed:** E7.7 (three-pillars informational closure), E6.7
(HKM Pareto point at `T = N^{8/15}`, `S = N^{1/3}`).
**Adjacent edges:** E5.1 (BPSW-TC⁰), E6.5 (HKM elementary), E6.6
(Aggarwal binary search), E6.8 (Dusart bracket), E6.9 (primecount
hierarchy), E7.6 (Meissel-Lehmer pebbling lower bound).
**Session:** S81 (B-grade target).
**Type:** structural classification of the algorithmic landscape.

---

## Background

E7.7 says: across 15+ tested intermediate-quantity families, every
algorithm for π(x) routes back to one of three **informationally-complete
encodings**:

```
   prime positions   |   zeta zeros   |   floor values  {⌊x/k⌋ : 1 ≤ k ≤ √x}
```

E6.7 says: the Helfgott-Kessler-Mendlovic (HKM) algorithm achieves
`(time = N^{8/15}, space = N^{1/3})`, the **Pareto-best simultaneously
sub-x^{2/3} time and sub-x^{1/2} space** algorithm known.

What is missing: a **structural placement** of HKM on the three-pillar
diagram. Concretely, three sub-questions:

1. Which pillar does HKM live on?
2. Is HKM on its pillar's Pareto frontier?
3. Is there ANY pillar dominating HKM at HKM's `(T, S)` point — i.e.,
   could a different pillar achieve `T ≤ N^{8/15}` AND `S ≤ N^{1/3}`
   simultaneously?

The third question matters because if the answer is no, then HKM's
specific (T, S) point is achievable on the floor pillar **only**, and
the algorithmic landscape decomposes cleanly: each pillar has a
distinct Pareto-optimal region.

---

## Construction

A `(T, S)`-classified catalog of the project's tracked π(x) algorithms,
plus the per-pillar Pareto frontier, plus a cross-pillar dominance check.

### Algorithm catalog (input data)

For each algorithm we record:

- `name` — descriptive name + citation source.
- `pillar` — one of {prime_positions, zeta_zeros, floor_values}, by
  the algorithm's primary information source.
- `time_exp` — `α` such that `T = O(N^α · polylog(N))`.
- `space_exp` — `β` such that `S = O(N^β · polylog(N))`.
- `assumption` — unconditional / GRH / BPSW / etc.
- `source_edge` — the EDGES.md edge ID where the algorithm is documented.

### Pareto frontier per pillar

A point `(α, β)` is on pillar `P`'s Pareto frontier iff there is an
algorithm with that `(α, β)` and no other algorithm on the same pillar
has `(α', β')` with `α' ≤ α`, `β' ≤ β`, and at least one strict.

### Cross-pillar dominance

A point `(α₀, β₀)` is **dominated cross-pillar** iff there exists an
algorithm on a *different* pillar with `(α', β')` strictly better
elementwise.

### HKM-specific structural placement

We compute:

- `pillar(HKM)` — by algorithm definition (HKM uses NTT-accelerated
  Dirichlet convolution = floor-values).
- Whether HKM is on `pillar(HKM)`'s Pareto frontier.
- Whether HKM is dominated cross-pillar.
- The minimum-time and minimum-space algorithm per pillar.

---

## Pre-stated falsification criteria

The construction is **falsified** if any of:

- **F1** — Some algorithm in the catalog has `(α, β) <` HKM elementwise
  with at least one strict inequality. ⇒ HKM is not Pareto-optimal.
- **F2** — Some algorithm on a non-floor pillar has `(α, β) ≤` HKM
  elementwise. ⇒ Floor-pillar placement is artifactual; HKM is
  achievable cross-pillar at the same point.
- **F3** — All three pillars have IDENTICAL Pareto frontiers in the
  region `α ≤ 1` AND `β ≤ 1`. ⇒ The pillar classification is
  algorithmically vacuous (E7.7 is informationally non-trivial but
  algorithmically meaningless on the (T,S) plane).
- **F4** — Some algorithm achieves `T·S < N^{5/6}/log N`. ⇒
  Contradicts E7.6 (Lucy-DAG pebbling lower bound) at the floor
  pillar; the construction's input data is corrupt.

The construction **passes** iff none of F1-F4 fires.

---

## What this construction adds beyond E7.7 + E6.7 individually

E7.7 is purely *informational*: it says the three pillars exhaust
information sources, but says nothing about which (T, S) regions are
achievable on each pillar.

E6.7 places ONE algorithm (HKM) at ONE point (`(8/15, 1/3)`), but does
not classify it.

The composition adds:

- **A pillar-by-pillar (T, S) Pareto map**, computed from the project's
  catalog of tracked algorithms.
- **A cross-pillar dominance check** at HKM's specific point.
- **A structural claim** about which Pareto-frontier region each pillar
  controls.

The *novel content* is the structural placement, not the individual
algorithms (which are all in the literature).

---

## Why this is B-grade, not A-grade

This is a **structural classification of existing algorithms** — exactly
the kind of synthesis CLAUDE.md flags as B-grade. An A-grade result
would require constructing a *new* algorithm whose (T, S) point falls
outside all current pillar frontiers, OR proving a pillar-restricted
lower bound that strictly improves E7.6.

This construction's value is in (a) verifying that the floor pillar is
asymptotically tight via E7.6, (b) showing the zero pillar wins on time
alone but loses on T·S, and (c) producing a permanent reference diagram
the project can update as new algorithms are added.

---

## File layout

- `definition.md` — this file (precise statement, falsification
  criterion).
- `pillar_tradeoff_diagram.py` — algorithm catalog + Pareto + dominance
  computation; emits a JSON table and a tabular ASCII diagram.
- `pillar_tradeoff_diagram_results.md` — outcomes + closure outcome.
- `pillar_tradeoff_diagram_results.json` — raw numerical outputs.
