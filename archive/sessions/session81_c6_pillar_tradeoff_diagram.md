# Session 81 — C6 Three-Pillars × HKM Tradeoff Diagram

**Mode:** CONSTRUCTION (composition target).
**Target:** `NOVELTY_CHALLENGES.md` §1.C6.
**Edges composed:** E7.7 (three-pillars informational closure) + E6.7
(HKM at `(N^{8/15}, N^{1/3})`).
**Adjacent edges cited:** E5.1, E6.5, E6.6, E6.8, E6.9, E7.6, E7.9.
**Channeled mathematician:** Tao (Pareto-frontier classification across
distinct algorithmic frames).

---

## What the session built

`experiments/constructions/pillar_tradeoff_diagram/` containing:

- `definition.md` — precise statement of the composition C6 with
  pre-stated F1-F4 falsification criteria.
- `pillar_tradeoff_diagram.py` — algorithm catalog (14 entries across
  the three pillars), per-pillar Pareto frontier computation,
  cross-pillar dominance check at HKM's point, E7.6 saturation table,
  ASCII tradeoff diagram, JSON export.
- `pillar_tradeoff_diagram_results.md` — results, falsification verdict,
  subtleties, closure outcome.
- `pillar_tradeoff_diagram_results.json` — raw numerical outputs.
- `run.log` — full stdout from the script.

The catalog encodes each algorithm as `(alpha, beta)` with
`T = O~(N^alpha)`, `S = O~(N^beta)` (polylog absorbed), tagged by
its primary information-source pillar:

```
prime_positions  4 entries (naive sieve, segmented sieve, wheel sieve, Aggarwal)
zeta_zeros       4 entries (direct expl-formula, Galway, Odlyzko-Schoenhage, Connes)
floor_values     6 entries (Lucy DP, Lucy+Fenwick, Meissel-Lehmer, Gourdon, primecount, HKM)
```

---

## Three structural findings

**(1) HKM dominates the floor pillar.** All five other floor-pillar
entries are dominated by HKM elementwise. Meissel-Lehmer at
`(2/3, 1/3)` is dominated in time only; primecount at `(2/3, 1/2)` and
Gourdon at `(2/3, 1/2)` are dominated in both axes. The floor-pillar
Pareto frontier collapses to the single HKM point.

**(2) HKM uniqueness (the headline composition output).** No
zero-pillar or prime-pillar algorithm in the catalog achieves both
`T <= N^{8/15}` AND `S <= N^{1/3}` simultaneously. HKM's
`(8/15, 1/3)` point is *uniquely accessible on the floor pillar*. This
is the structural placement claim that did not previously exist in
the project; it sharpens E6.7 from "HKM is Pareto-best" (true on the
combined diagram) to "HKM is uniquely on the floor pillar at this
point" (which is a stronger structural assertion).

**(3) Non-overlapping pillar dominance regions.** The three pillars
control different axes:

- **Time-only minimum** is shared by prime+zeta pillars at `alpha=1/2`
  (Aggarwal-via-zero-pi, Galway, direct explicit formula). Floor pillar
  is `0.034` worse on time alone.
- **Space-only minimum** is unique to floor at `beta = 1/3` (HKM,
  Meissel-Lehmer). No zero/prime entry has `beta < 1/2`.
- **T*S minimum** is unique to floor at `13/15 = 0.867` (HKM). Zero
  and prime pillars are flat at `T*S = 1`.

Floor pillar saturates the E7.6 lower bound (`T*S >= N^{5/6}` =
`N^{0.833}`) to within `N^{0.034}`.

---

## Aggarwal as a meta-algorithm (S81 sub-observation)

E6.6 (Aggarwal binary search) makes `O(log n)` calls to a `pi(x)`
sub-routine. Its time exponent inherits from the chosen sub-routine.
Coupling Aggarwal to HKM (instead of zero-pillar pi) gives
Aggarwal-via-HKM at `(8/15, 1/3)` — Aggarwal effectively migrates
into the floor pillar as a meta-algorithm.

This is a small structural observation worth noting: the prime pillar
in our catalog is dominated by what is, in effect, a wrapper around
the other pillars. The "true" prime pillar (sieve-of-Eratosthenes
class) is at `(1, 1)` or `(1, sqrt)` and is far from competitive.

---

## Falsification verdict

| F# | criterion | result |
|----|-----------|--------|
| F1 | HKM dominated by some catalog entry | **PASS** (no violator) |
| F2 | HKM cross-pillar dominated | **PASS** (no violator) |
| F3 | All three pillars have identical Pareto frontiers | **PASS** (frontiers differ) |
| F4 | Floor-pillar entry achieves T*S < 5/6 | **PASS** (no violator) |

**Construction passes.**

---

## Self-evaluation (CLAUDE.md required questions)

**Q1. What did this session produce that was not in the project before?**

- An explicit `(alpha, beta)` catalog of 14 pi(x) algorithms classified
  by the three E7.7 pillars.
- A per-pillar Pareto frontier computation showing HKM dominates the
  entire floor pillar elementwise.
- The **HKM uniqueness lemma**: HKM's `(8/15, 1/3)` point has no
  cross-pillar dominator in the catalog.
- A structural "non-overlapping pillar dominance" claim: each pillar
  controls a distinct corner of the (T, S) plane.
- A reproducible Python catalog `(CATALOG)` that future agents can
  extend by appending entries.

**Q2. What edges did this work compose or cite?**

- **Composed:** E7.7 + E6.7 (the C6 target).
- **Cited (catalog source / lower bound):** E5.1, E6.5, E6.6, E6.8,
  E6.9, E7.6, E7.9.
- **Annotated:** E6.7 and E7.7 in `EDGES.md` with back-pointers to the
  construction.

**Q3. If only duplicate closures, why?**

Not applicable — produced a real construction with a structural
placement claim (HKM uniqueness). However the result is
**REFINEMENT** mode: it does not introduce a new edge. The placement
claim follows from existing edges; what's new is the explicit
classification, not a new mathematical fact about pi(x).

**Q4. Next action for next agent.**

Pick **C3** (Brandt × per-bit; remaining open challenge with novelty
potential) or **C4** (Aggarwal × Dusart × BPSW unified library;
engineering work, belongs in algorithms/). C3 likely closes within
~30 minutes as a structural duplicate of E5.8: the per-bit version
of Brandt TRAVERSE cannot escape obstructions O1-O4 because the
parametric J ∈ {0..log N} provides only polylog space, not the
Chaitin-Ω continuum that O1 requires. But that structural reason itself
would be a clean refinement of E5.8 worth filing.

After C3/C4 are closed, all six §1 composition targets (C1-C6) will
be built. At that point, **the natural next step is the §1 meta-synthesis**
("which edge pairs yielded structure, which collapsed?") promised in
RESEARCH_AGENDA.md Arc 4 milestones, or **a new ATTACK_VECTORS.md
entry** (frontier targets are running thin per the post-S77 critique).

---

## Self-grade

**B-grade (refinement, not novelty).** This was a structural
classification of existing algorithms — exactly what CLAUDE.md flags
as "B-grade at best". There is no new mathematical fact beyond what
E6.7 and E7.7 already imply individually; what is new is the explicit
placement and the cross-pillar dominance check.

For an A-grade outcome the construction would have needed to either:

- Discover an algorithm with `(alpha, beta)` strictly outside all three
  current pillar frontiers (would have suggested a fourth pillar);
- Prove a pillar-restricted lower bound strictly better than E7.6 on
  the zero or prime pillars (would have established structural
  pillar separation);
- Find a cross-pillar transformation that maps HKM to a polylog
  algorithm (would have been the breakthrough).

None happened. The construction passes its falsifiers and produces a
permanent reference diagram, but does not move the polylog frontier.

Honestly self-graded B, not A.

---

## Files written

```
experiments/constructions/pillar_tradeoff_diagram/
├── definition.md
├── pillar_tradeoff_diagram.py
├── pillar_tradeoff_diagram_results.md
├── pillar_tradeoff_diagram_results.json
└── run.log

status/CLOSED_PATHS.md         (S81 row appended)
EDGES.md                       (E6.7 and E7.7 annotated)
NOVELTY_CHALLENGES.md          (C6 marked BUILT)
RESEARCH_AGENDA.md             (Arc 4 C6 marked done; next-action updated)
archive/sessions/session81_c6_pillar_tradeoff_diagram.md  (this file)
```
