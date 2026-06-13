# C11 / D17.b — Discrete Morse on the Squarefree Divisibility Hasse: Results

**Session:** 422 (construction mode, follow-up to D17 / S232).
**Channel:** Erdős — combinatorial extremal counting (same as D17).
**Cross-domain ingredient:** Forman 2002 *Sém. Lothar. Combin.* 48
(discrete Morse theory) + Benedetti–Lutz 2014 *Exp. Math.* 23, 66
(random discrete Morse). Status: **USED-E (S232)** at D17;
**USED-E refinement (S422)** with a sharpened identity.
**Self-graded:** **B** (substantive refinement: an exact arithmetic
identity for `m_0` on a structurally distinct combinatorial object).

## Setup

Built `H_N^sqfree` (squarefree-only divisibility Hasse 1-skeleton) at
`N ∈ {64, 128, 256, 512, 1024, 2048, 4096, 8192}`. Ran greedy random
elementary Morse collapse (Forman 2002) for `R = 20` seeds per `N`,
plus Erdős–Rényi `G(|V|, |E|)` baselines (20 seeds each). Determinism
test: 200 random seeds at `N ∈ {64, 256, 1024, 4096}`. Cascade
breakdown: deterministic-order collapse with wave assignment per
vertex.

| N    | `|V_sqf|` | `|E|`  | `χ`    | `π(N)` | `π(N) − π(N/2)` |
|------|-----------|--------|--------|--------|------------------|
| 64   | 39        | 60     | −21    | 18     | 7                |
| 128  | 78        | 132    | −54    | 31     | 13               |
| 256  | 157       | 285    | −128   | 54     | 23               |
| 512  | 314       | 601    | −287   | 97     | 43               |
| 1024 | 624       | 1248   | −624   | 172    | 75               |
| 2048 | 1245      | 2589   | −1344  | 309    | 137              |
| 4096 | 2491      | 5373   | −2882  | 564    | 255              |
| 8192 | 4982      | 11096  | −6114  | 1028   | 464              |

`|V_sqf|/N`: 0.609, 0.609, 0.613, 0.613, 0.609, 0.608, 0.608, 0.608 —
matches `6/π² ≈ 0.6079`.

## Empirical findings

### F1 (A-grade) — FAILS: linear scaling, no polylog compression

| N    | `m_0(div)` | `m_0(div)/|V|` |
|------|------------|------------------|
| 64   | 32         | 0.8205           |
| 128  | 65         | 0.8333           |
| 256  | 134        | 0.8535           |
| 512  | 271        | 0.8631           |
| 1024 | 549        | 0.8798           |
| 2048 | 1108       | 0.8900           |
| 4096 | 2236       | 0.8976           |
| 8192 | 4518       | 0.9069           |

Ratio `m_0/|V|` monotone-increases toward 1; no polylog compression.

### F3 (B-grade) — HOLDS exactly: a sharp two-term identity

> **Theorem (S422, proven analytically; verified pointwise to N=8192).**
> For every `N ≥ 6`, the greedy random elementary Morse collapse on
> `H_N^sqfree` halts in **exactly one wave** with collapse count
>
>   `collapses(H_N^sqfree) = π(N) − π(N/2)`
>
> exactly, with no chained-collapse residual. Equivalently:
>
>   `m_0(H_N^sqfree) = |V_sqfree(N)| − (π(N) − π(N/2))`
>   `m_1(H_N^sqfree) = |E(H_N^sqfree)| − (π(N) − π(N/2))`.

Pointwise verification:

| N    | collapses | `π(N) − π(N/2)` | `ε(N)` |
|------|-----------|------------------|--------|
| 64   | 7         | 7                | 0      |
| 128  | 13        | 13               | 0      |
| 256  | 23        | 23               | 0      |
| 512  | 43        | 43               | 0      |
| 1024 | 75        | 75               | 0      |
| 2048 | 137       | 137              | 0      |
| 4096 | 255       | 255              | 0      |
| 8192 | 464       | 464              | 0      |

`ε(N) ≡ 0` across every measured `N` — strictly sharper than D17's
`ε(N) ≡ 1` chained-collapse residual, and the `Π_pow(N) ~ Θ(√N/log N)`
prime-power term is also absent (no prime powers among squarefree
vertices).

**Proof sketch.** Compute degrees in `H_N^sqfree`:

* `deg(1) = π(N)` (1 has up-edges to every prime ≤ N).
* For prime `p ∈ V_sqfree(N)` with `p ≤ N`, `deg(p) = 1 + |{q prime
  : q ≠ p, qp ≤ N}| = 1 + π(N/p) − [p² ≤ N]`.
* For squarefree `n ∈ V_sqfree(N)` with `ω(n) ≥ 2`, `deg(n) = ω(n)
  + |{q prime : q ∤ n, qn ≤ N}| ≥ ω(n) ≥ 2`.

Wave 0 leaves are vertices with `deg = 1`. Only candidates:

* `p ∈ (N/2, N]`: `π(N/p) = 0` and `p² > N`, so `deg(p) = 1`. PASS.
* `p ∈ (N/3, N/2]`: `π(N/p) = 1`, so `deg(p) = 2`. FAIL.
* All other vertices: degree ≥ 2.

So the initial leaves are exactly the primes in `(N/2, N]`, count
`= π(N) − π(N/2)`.

Each leaf `p > N/2` has a single incident edge `(1, p)`. After peeling
all `π(N) − π(N/2)` leaves, the only modified vertex is `1`, whose
degree drops to `π(N) − (π(N) − π(N/2)) = π(N/2) ≥ 2` for `N ≥ 6`.

Vertex `1` is not a wave-1 leaf, and no other vertex has had its degree
modified. Therefore wave 1 is empty, and the cascade halts.

The collapse count is order-independent (the leaves are mutually
non-adjacent), so the greedy random Morse output is deterministic. ∎

### F2 (B-grade) — FAILS: ER baseline systematically larger

| N    | `m_0(div)` | `m_0(ER mean)` | `Δ = div − ER` | `|Δ|/|V|` |
|------|------------|------------------|------------------|--------------|
| 64   | 32         | 32.9             | −0.85            | 0.022        |
| 128  | 65         | 68.2             | −3.25            | 0.042        |
| 256  | 134        | 141.2            | −7.15            | 0.046        |
| 512  | 271        | 283.9            | −12.85           | 0.041        |
| 1024 | 549        | 572.8            | −23.80           | 0.038        |
| 2048 | 1108       | 1159.9           | −51.90           | 0.042        |
| 4096 | 2236       | 2333.8           | −97.85           | 0.039        |
| 8192 | 4518       | 4707.8           | −189.75          | 0.038        |

The squarefree-Hasse is **MORE collapsible** than the matched-density
ER baseline by a stable `≈ 4 % of |V|`. This is *opposite* to D17's
near-matched outcome (`≈ 0.5 %–2 %`) and shows the squarefree
restriction *amplifies* the structural distinguishability of
divisibility-poset Morse complexity from random graphs of the same
density.

The `4 %` gap is `(π(N) − π(N/2)) − E[ER initial leaves]` to leading
order: ER on `(|V|, |E|)` with average degree `2|E|/|V| ≈ 4.5` for
this regime has expected initial-leaf count `|V| · (2|E|/|V|)
exp(−2|E|/|V|) ≈ 0.045 |V|`, slightly less than `(π(N) − π(N/2))/|V|
≈ 0.085`. (More precisely, ER also has some chained collapse, which
narrows but does not close the gap.)

### F4 — Morse-rigidity HOLDS: deterministic across 200 seeds

| N    | distinct `(m_0, m_1)` outputs / 200 seeds |
|------|------------------------------------------|
| 64   | 1                                        |
| 256  | 1                                        |
| 1024 | 1                                        |
| 4096 | 1                                        |

This *follows from the proof above* — the wave-0 leaves are mutually
non-adjacent so the collapse is order-independent. The squarefree
restriction inherits the divisibility lattice's Morse rigidity exactly.

## What this refines

* **D17 / S232** (full divisibility Hasse): identity
  `collapses(H_N) = (π(N) − π(N/2)) + Π_pow(N) + 1`. The squarefree
  restriction kills the `Π_pow(N)` term (no prime-power vertices) AND
  the constant `+1` chained-collapse residual (no chain triggered by a
  level-2 prime power). The distinguishing arithmetic content of the
  full-lattice identity is *concentrated entirely in the prime-power
  contribution* — the squarefree-only baseline is "primes-in-(N/2, N]"
  cleanly, no residual.
* **E2.28** (Baker–Norine `q`-reduced form on `H_N`): provides exact
  closed forms for the *q*-reduced divisor of `D_P^N` on the same
  combinatorial objects. Both probes find that the order-based
  topology of (`Z`, `|`) reduces to `π(N) − π(N/2)` plus low-order
  arithmetic corrections. The Morse path here is order-independent
  (greedy Morse rigidity); E2.28's reduced form is also order-
  independent (Baker–Norine canonical-form theorem). Two independent
  topological invariants on the same poset converge on
  `(π(N) − π(N/2))` as the leading combinatorial term.
* **E1.5 / E1.6** wheel-bisection: at "wheel `W = ∞`" (no wheel
  filtering), the squarefree restriction here can be read as a
  multiplicative indicator pre-filter analogous to the `μ²` indicator;
  the identity's `(π(N) − π(N/2))` term is the "squarefree-coprime"
  (in fact squarefree-prime) count in the half-window `(N/2, N]`.

## Why this is B-grade (and not A-grade)

* **F1 (polylog) FAILS** — `m_0 = Θ(N)` linear, no algorithmic opening.
* **F3 (closed form) HOLDS exactly** with no residual term — a
  *cleaner* arithmetic identity than D17's, but it is again
  *circular*: `m_0(H_N^sqfree)` is `|V_sqfree(N)|` (a Mertens-density
  count, polylog-evaluable) MINUS `π(N) − π(N/2)` (a primes-in-interval
  count, no easier than `π(N)`). The Morse-complexity-of-squarefree-
  Hasse problem reduces *exactly* to `π(N) − π(N/2)`, which is the
  same hard subroutine as D17's.
* **F4 (Morse rigidity) HOLDS analytically** by the proof above —
  not just empirically. This is sharper than D17's empirical-only
  rigidity claim.

The combinatorial restriction to squarefree integers is a *clean way
to isolate the prime-counting content* of the divisibility-poset
Morse complexity from the prime-power and chained-collapse noise.
It does not, however, reduce the core complexity barrier.

## What would have been A-grade (and why it's blocked)

The A-grade hope was that the Boolean-lattice shellability *of the
order complex* would propagate to the *1-skeleton* of the truncated
Boolean lattice, giving `m_0 = O(polylog N)`.

This fails for two structural reasons:

1. The 1-skeleton of the **full** Boolean lattice on `n` primes is
   the `n`-cube graph, which has minimum degree `n` and **no
   degree-1 vertices** — not collapsible by greedy elementary Morse
   *as a graph*. Shellability lives one dimension up (in the order
   complex with chains, not the 1-skeleton).
2. The truncation to `∏ p ≤ N` introduces leaves only at the *top*
   (singleton-prime sets `{p}` with `p > N/2`). The rest of the
   truncated cube graph still has no degree-1 vertices.

Any improvement on the greedy upper bound would require *non-greedy*
Morse matchings (Joswig–Pfetsch 2006: NP-hard in general) and would
have to encode `μ²(n)` in the dual cohomology basis — equivalent to
the squarefree-counting problem itself, which is `Θ(√N)` in the
Pillai sieve. So the failure is *circular* in the same sense as D17.

## Composes

* **D17 / S232** (full-lattice identity)
* **E2.28** (Baker–Norine on `H_N` and `Γ_N`)
* **E1.6 / E2.10** (parity-bisection wheel decomposition)
* Cross-domain: Forman 2002, Benedetti–Lutz 2014, Joswig–Pfetsch 2006

## Files

* `d17b_squarefree_morse.py` — main scan (Hasse + ER + determinism +
  cascade breakdown).
* `d17b_squarefree_morse_data.json` — raw data per N.
* `definition.md` — pre-stated definition + falsifiers.
* `scan_output.log` — runtime log.

## Successor proposals

**D17.b.i — Δ(H_N^sqfree) full order complex.** Lift from 1-skeleton
to the full chain complex of the squarefree-poset order complex.
Cell count is `~ Σ_{n sqf ≤ N} ω(n)!` (chains correspond to permutations
of `n`'s prime factorisation), still tractable for `N ≤ 256`. The
*order complex* of the unrestricted Boolean lattice IS shellable with
single critical top cell; the truncation at `∏ p ≤ N` may genuinely
reduce `m_0` below `Θ(N)`. **A-grade target if** the chain-level
`m_0 = O(polylog N)` even though the 1-skeleton failed. Cost: 1-2
sessions; save under
`experiments/topological/d17b_squarefree_morse/order_complex/`.

**D17.b.ii — Multiplicative-indicator generalisation.** Replace
`μ²(n)` by other multiplicative indicators (k-free numbers,
k-almost-primes, smooth-y-rough numbers). Each defines a sub-Hasse
of `(Z, |) ∩ [1, N]`. Expected: the same `π(N) − π(N/2)`-leading-
term identity holds with a different second-order correction
matching the indicator's level-2 structure. Cost: 1 session per
indicator family; save under
`experiments/topological/d17b_squarefree_morse/multiplicative_indicators/`.

**D17.b.iii — Lean formalisation.** The proof is short (vertex
degree case-analysis + wave-1 emptiness lemma). Tractable Lean 4
target; add to L1 / L6 queue. Cost: 1-2 sessions; save under
`experiments/formalisations/L_d17b_morse_identity/`.

## Session-end self-evaluation

1. **What did I produce that was not in the project before this
   session?**
   * The closed-form identity `collapses(H_N^sqfree) = π(N) − π(N/2)`
     with `ε(N) ≡ 0` (sharpening D17's identity that had a
     `Π_pow + 1` correction).
   * A two-line analytical proof of the identity via vertex
     degree case-analysis and wave-1 emptiness (D17 was empirically
     established only).
   * The structural observation that the squarefree restriction
     *amplifies* the gap between divisibility-Hasse and ER-baseline
     Morse complexity (`≈ 4 %` vs D17's `0.5–2 %`).
   * Morse rigidity of `H_N^sqfree` as an analytical (not just
     empirical) consequence of the wave-0 leaves being mutually
     non-adjacent.

2. **What edges did my work compose or cite?**
   * D17 / S232 (full-lattice identity, refined).
   * E2.28 (Baker–Norine on the same poset).
   * E1.6 / E2.10 (parity decomposition, "W = ∞" wheel limit).
   * Cross-domain: Forman 2002, Benedetti–Lutz 2014, Joswig–Pfetsch
     2006.

3. **If my session produced only duplicate closures, why?**
   Not a duplicate. The identity is a *strict refinement* of D17 with
   no `Π_pow` and no `+1` residual, an analytical proof, and a
   structural amplification of the ER-distinguishability gap.

4. **Next-action for the next agent.**
   Pick up D17.b.i (full order complex) — A-grade still on the table
   for the chain-level `m_0`. Or D17.b.iii (Lean formalisation).
   Listed in NOVELTY_CHALLENGES.md §1 successor block.
