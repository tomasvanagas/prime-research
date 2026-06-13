# C11 / D17.b — Squarefree-Only Divisibility Hasse Discrete Morse: Definition

## Object

`H_N^sqfree` — the 1-skeleton of the order complex of the divisibility
poset restricted to squarefree integers in `[1, N]`.

* Vertices: `V_sqfree(N) := { n ∈ [1, N] : μ(n) ≠ 0 }`. By Mertens,
  `|V_sqfree(N)| = (6/π²) N + O(√N)`.
* Edges (covering relations): `E := { (m, mp) : m, mp ∈ V_sqfree(N),
  p prime, p ∤ m, mp ≤ N }`. Equivalently, `m | mp` with `mp/m` prime
  and the result still squarefree.

`H_N^sqfree` is isomorphic to the truncation of the Boolean lattice on
the primes `≤ N` to subsets `S ⊆ {primes ≤ N}` with `∏_{p ∈ S} p ≤ N`.

## Discrete Morse functional

Forman's greedy random elementary collapse (Forman 2002 *Sém. Lothar.
Combin.* 48; Benedetti–Lutz 2014 *Exp. Math.* 23, 66) on the 1-skeleton:

1. While `∃ v` of degree 1 in the surviving subgraph:
2. pick uniformly at random a degree-1 vertex `v` and its unique
   incident edge `e = (v, u)`,
3. discard the pair `(v, e)` (an *elementary collapse*),
4. update degrees.

Output: `(m_0, m_1) = (#surviving vertices, #surviving edges)` — the
Morse critical-cell counts of a (random) discrete Morse function on the
graph. Forman's identity `m_0 − m_1 = χ(H_N^sqfree) = |V| − |E|` holds
at every step.

## Edge IDs composed

* **E2.28** — Baker–Norine `q`-reduced divisor form on the divisibility
  Hasse / cover graph: another order-based topological probe on the
  same poset, with explicit prime-supported closed forms. C11 is the
  Morse-theoretic counterpart for the squarefree restriction.
* **D17 closure (S232)** — discrete-Morse 1-skeleton complexity of the
  full divisibility Hasse `H_N`. Identity:
  `collapses(H_N) = (π(N) − π(N/2)) + Π_pow(N) + 1`. C11 sharpens this
  by removing both the `Π_pow(N) ~ Θ(√N/log N)` prime-power term and
  the `+1` chained-collapse residual.
* **E1.6 / E2.10** — primality-vs-Liouville parity decomposition feeds
  the wheel-coprimality lattice C10 / S219 explores; the squarefree
  restriction here is the wheel-W-free (i.e., `W = ∞`) limit of a
  multiplicative-indicator Hasse.

## Intended relationship to π(x)

If `m_0(H_N^sqfree) = polylog N`, the *dual cohomology basis* of the
collapsed complex would encode `μ²(n)` (the squarefree primality
witness) in `polylog N` bits — a partial-positive
Correlation-Dichotomy-shaped result on the *adjacent* problem of
squarefree counting. The construction does NOT a-priori reduce `π(x)`
itself, since χ_P ⊊ μ² in support; but a polylog topological
compression of any natural multiplicative indicator would be
unprecedented.

## Falsifiers (pre-stated, BEFORE running the script)

* **F1 — A-grade target.** `m_0(H_N^sqfree) = O(polylog N)`. Expected
  shape if the truncated Boolean lattice inherits enough of the
  full lattice's shellability to cascade past initial leaves.
* **F2 — B-grade ER baseline match.** `|m_0(H_N^sqfree) − m_0(ER(|V|, |E|))|
  / |V| ≤ 0.02` uniformly in `N`. Expected if the Morse complexity is
  not distinguishable from matched-density random graphs (analogue of
  D17's outcome on the full lattice).
* **F3 — B-grade closed-form arithmetic identity.** `m_0(H_N^sqfree)`
  admits an exact decomposition into terms involving `π`, `π/2`, and
  small ω-strata counts (analogue of D17's
  `(π(N) − π(N/2)) + Π_pow(N) + 1` identity).
* **F4 — Morse-rigidity.** Greedy random Morse output deterministic
  across ≥ 200 random seeds at every measured `N` (analogue of D17 F2).

## Outcome direction

* **F1 FAILS** → expected (linear scaling), no algorithmic opening.
* **F2 PASSES, F3 fails** → 39th pseudorandomness measure, B-grade
  negative-shape closure orthogonal to D17.
* **F2 fails, F3 PASSES** → B-grade refinement of D17 with a sharper
  identity; locates the squarefree restriction's distinguishing
  arithmetic content.
* **F2 and F3 BOTH pass** → strongest B-grade outcome; refinement of
  D17 + new ER-baseline-distinguishability fact.
* **F4 — independent property** measured alongside the above.

See `d17b_squarefree_morse_results.md` for the experimental verdict.
