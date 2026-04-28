# D45 — Baker-Norine graph Riemann-Roch / chip-firing rank of the prime divisor

**Session:** 161 (wild swing on `ATTACK_VECTORS.md §D45`)
**Date:** 2026-04-28
**Mode:** wild swing (one-attempt, full-session, frontier attack)
**Self-grade:** **B** (substantive structural finding + clean closed-form identity + Riemann-Roch ceiling on rank invariant)

## Goal

Test whether the divisor rank `r_G(D_P^N)` of the prime divisor on an
arithmetic graph `G` exhibits non-trivial structure (Brill-Noether
specialness, deviation from random matched-degree controls, or higher-rank
cohomological content). Cross-domain: Baker-Norine 2007 *Adv. Math.* 215,
766 = arXiv:math/0608360 (graph Riemann-Roch); Dhar 1990 *Phys. Rev. Lett.*
64, 1613 (chip-firing burning algorithm).

## What was built

`baker_norine_chi_p.py` — Dhar's burning algorithm for q-reduction on a
chip-firing graph; Baker-Norine Algorithm 4.1-style rank computation by
exhaustive multiset enumeration; matched-density random-divisor controls.

`baker_norine_hasse.py` — same on the Hasse cover graph (a — b iff b/a is
prime).

`run_full.py` — driver that produces `full_results.json`.

Graphs tested:
- **`Γ_N`** (divisibility): `V = [1, N]`, edge (a, b) iff a | b, a ≠ b.
- **`H_N`** (Hasse cover): `V = [1, N]`, edge (a, b) iff b = pa for some
  prime p.

Sink chosen `q = 1` for both (highest-degree vertex; both graphs are
connected because 1 divides everything, hence is universal).

## Findings

### F1 — Closed-form q-reduced identity on Γ_N (divisibility graph)

**Theorem (verified for N ∈ {16, 32, 64, 128, 256, 512}):** On the
divisibility graph `Γ_N`, with sink `q = 1`, the q-reduced form of the
prime divisor `D_P^N` is

  **D'_P^N = (π(N) − π(N/2)) · δ_1 + Σ_{p prime, p ≤ N/2} δ_p**

i.e. `D'(1) = π(N) − π(N/2)`, `D'(p) = 1` for primes `p ≤ N/2`, and
`D'(v) = 0` for all other v.

**Proof (sketch):** Dhar's burn from q=1 first burns all composites c
(adjacent to 1 because 1 | c, with `D(c) = 0` < count = 1) and also all
primes p ≤ N/2 (after composites burn, p has burnt-count ≥ 2 from its
multiples 2p, 3p, ... > 1 = D(p)). Primes p > N/2 are degree-1 vertices
(only neighbor is 1) with `D(p) = 1`, burnt-count = 1, fail to burn.
Firing them moves their chips to q. ∎

### F2 — Closed-form q-reduced identity on Hasse(N)

**Theorem (verified for N ∈ {16, 32, 64, 128, 256, 512}):** On the Hasse
cover graph `H_N`, with sink `q = 1`, the q-reduced form of the prime
divisor `D_P^N` is

  **D'_P^N = π(N) · δ_1.**

**Proof (sketch):** In `H_N`, the neighbors of 1 are exactly the primes
≤ N (since (1, k) is a cover edge iff k = p prime). Initial burn from
q=1: every prime p with `D(p) = 1` has burnt-count 1, fails 1 > 1. Every
composite c has burnt-count 0 (1 ∉ adj(c) in Hasse). So `unburnt = V \\
{1}`. Firing it once: each prime contributes 1 chip to q (the only
edge from p to a non-fired vertex is (1, p)), losing its chip in the
process; each composite has all its edges in `V \\ {1}`, so its
contributions cancel and D(c) is unchanged at 0. Result: `D'(1) =
π(N)`, all others 0. ∎

**This expresses π(N) as the q-component of the q-reduced divisor of
D_P on the Hasse graph — a chip-firing-theoretic representation of the
prime-counting function.**

### F3 — Generalisation to arithmetic divisors on Hasse(N)

**Generalised theorem:** On `H_N`, q=1, for any divisor D supported on
{primes only} ∪ {1} with D(v) = 0 for non-prime v ≠ 1:

  **D' = (D(1) + Σ_{p prime} D(p)) · δ_1.**

Verified on `D_P` (chips on primes) and `D_sqfree restricted to primes`
(same as `D_P` since primes are squarefree). Both reduce to π(N) δ_1.

For divisors with chips on COMPOSITES, the identity does not extend
cleanly: chips on composites with cover-divisors (= prime factors) that
have D = 0 burn IMMEDIATELY in Dhar, then cascade. For example:
- `D_λ_pos` (chips at v with λ(v) = +1, excluding 1) — chips concentrated
  on composites with even Ω; Dhar burning cascades chips to/from primes
  during reduction, yielding `D' ≠ D|composites`.
- `D_Omega2` (chips at semiprimes / prime squares) — same cascade.

For `D_μ_pos` (chips at squarefree v with μ(v) = +1, excluding 1):
empirically `D'(1) = 0` and `D'(c) = D(c)` for all composites c, so the
generalised identity DOES hold. Reason: μ-positive composites have all
prime factors with `D = 0` initially (since the chip is at the composite,
not its factors), so after cascading, no chip flow reaches q.

### F4 — Z-scores vs matched-density random divisors

For each N and graph, compared `D_P^N` against 20 random divisors of
the same degree (= π(N)), placed by uniformly choosing π(N) vertices to
hold one chip each.

| Graph        | N=32  | N=64  | N=128 |
|--------------|-------|-------|-------|
| Γ_N (div)    | +3.39 | +3.22 | +5.80 |
| Hasse(N)     | +8.62 | +11.11| +14.61|

**Z-score grows with N for Hasse — z(N=128) = +14.6σ.** The signal is
SUPER-LINEAR in N because:
- For D_P: D'(q) = π(N) ~ N / log N exactly.
- For random matched-degree: E[D'(q)] = π(N) · (deg(D))/(|V|) ≈ (π(N))² /
  N ~ N / (log N)².
- Gap = π(N) − E[random] ≈ π(N) (1 − π(N)/N) ~ π(N).
- Std ~ √(π(N) · (1 − π(N)/N)) ~ √π(N).
- Z-score ~ π(N) / √π(N) = √π(N) ~ √(N/log N) — **grows as √N / √log N**.

This z-score asymptotic is itself a structural prediction.

### F5 — Rank verdicts

The Baker-Norine **rank** `r_G(D_P^N)` is **0** for all N ∈ {32, 64, 128}
on both `Γ_N` and `H_N`, except `r_{H_16}(D_P^16) = 1` (small-N regime
where deg ≥ g, so Riemann-Roch forces r ≥ 1).

For each (N, graph), the failure of `r ≥ 1` was observed at a SPECIFIC
vertex: subtracting one chip from vertex 4 (or, for coprime, vertex 4 / 6)
gives a non-winnable divisor. The bottleneck is the same for `D_P` and
random matched divisors — it's a graph-topology property, not arithmetic.

**Riemann-Roch implication:** `r(K - D_P) = r(D_P) - π(N) - 1 + g = g - π(N) - 1`,
which is identical for D_P and any other effective degree-π(N) divisor by
Riemann-Roch's symmetry. **Specialness in the Baker-Norine sense is
trivially the same for D_P and random matched-degree divisors.**

## What this rules out (closure mode E / B-grade)

The Baker-Norine **rank** invariant on graphs `Γ_N` and `H_N` does NOT
distinguish `D_P^N` from a random matched-degree divisor: both have
`r = max(0, deg − g + 1)` modulo small-N exceptions.

The q-reduced **FORM** does distinguish them, but the distinguishing
features (D'(q), |support|) reduce to graph-neighborhood combinatorics:
- Γ_N: D'(q) = (# chips on degree-1 leaves) = π(N) − π(N/2) for D_P.
- H_N: D'(q) = (# chips on N(q)) = π(N) for D_P.

In both cases the "signal" is a chip-firing reformulation of a known
structural fact about the arithmetic graph (primes > N/2 are leaves /
primes are at distance 1 from 1 in the Hasse diagram).

## What would falsify this closure

- (F1) A divisor `D'` of degree d ≤ π(N) on `Γ_N` or `H_N` with
  `r(D') ≥ 1` whose chip-firing structure depends on a deep arithmetic
  invariant (e.g., the singular series, Hardy-Littlewood density). NOT
  observed.
- (F2) An ARITHMETIC graph `G` where the prime divisor's rank deviates
  STRUCTURALLY from random matched-degree divisors. The two graphs
  tested don't show this, but other graphs (e.g., the prime-gap
  multigraph, or Cayley with non-trivial generators) might.

## Edges composed / cited / refined

- **E2.x (NEW edge proposed)** — chip-firing q-reduced identity for the
  prime divisor on arithmetic graphs. Two closed forms, both derivable
  via Dhar's burn from elementary lattice structure.
- **E2.17 / E2.20 / E2.21 / E2.24** (pseudorandomness measures) —
  ADDS chip-firing q-reduced D'(q) as a *non-trivially nonzero* signal
  but the signal reduces to graph-neighborhood counts — does NOT break
  the broader pseudorandomness thesis.
- **E1.10 / E3.13** (GUE / random pseudorandomness measures all match) —
  REINFORCED at chip-firing **rank** level; structurally negative result
  for chip-firing-as-pseudorandomness-detector.
- **D45** (this attack vector) — **CLOSED** mode I (identification with
  graph-topology / Riemann-Roch): graph divisor rank is determined by
  (deg(D), graph G); the q-reduced FORM gives clean closed forms but
  encodes only the lattice-structural fact "primes occupy specific
  graph positions".

## Distinguishing features and what's left open

- The Hasse-graph identity `D'_P^N = π(N) δ_1` is the cleanest
  expression. It is NOT algorithmically useful for computing π(N)
  (building D_P^N requires knowing primes), but it's a clean structural
  representation.
- The signal **grows with N** at rate ~√(π(N)) on Hasse (super-linear in
  z-score), so the structural identity is asymptotically increasingly
  visible. This is itself a quantitative structural fact about the
  Hasse cover graph + prime divisor.
- The rank `r(D_P) = 0` is the generic value implied by Riemann-Roch
  for degree π(N) << genus; D_P is **not Brill-Noether-special** in
  the sense `r > deg - g`. This is the negative content.
- An **A-grade** outcome would have required: (a) finding a graph where
  `r(D_P)` exceeds Brill-Noether by ≥ 1 *with* a structural arithmetic
  reason (e.g., Hardy-Littlewood density), or (b) constructing a
  chip-firing-based polylog π(x) algorithm. Neither was achieved.

## Reproducibility

```bash
cd experiments/algebraic/baker_norine_chi_p
python3 run_full.py   # produces full_results.json
```

All scripts are deterministic (seed=42 for random divisor sampling).
Times: total run ~ 3 seconds for N up to 128 on the standard machine.

## Files

- `baker_norine_chi_p.py` — Dhar burn, q-reduce, rank, divisibility graph.
- `baker_norine_hasse.py` — Hasse cover graph + ratio graph variants.
- `baker_norine_signatures.py` — q-reduced signature extraction battery.
- `run_full.py` — full-experiment driver.
- `full_results.json` — all numbers in machine-readable form.
- `full_run.log` — captured stdout from `run_full.py`.

## Cross-domain references

- Baker, Norine 2007 *Adv. Math.* 215, 766 = arXiv:math/0608360 — graph
  Riemann-Roch (canonical reference).
- Dhar 1990 *Phys. Rev. Lett.* 64, 1613 — chip-firing / sandpile burning.
- Bjorner, Lovász, Shor 1991 *Eur. J. Combin.* 12, 283 — chip-firing
  abelian property.
- Backman, Baker, Yuen 2019 — algorithmic q-reduction.
