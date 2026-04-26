# N1 — Tensor-Network Compression Family Closure

**Target:** `NOVELTY_CHALLENGES.md` §4.N1.
**Edges composed:** E2.1 (MPS bond-dim identity), E1.9 (φ(x,a) 2D rank),
E6.3 (DCT/wavelet of δ).
**Session:** S77 (B-grade target).
**Type:** family-level negative-shape conjecture.

---

## Background

E2.1 establishes for the prime indicator `χ_P : [1, W^d] → {0,1}`,
reshaped as a `(W, W, ..., W)` tensor (d copies):

```
   rank M^(j)  =  min ( W^j , φ(W) · W^(d-j-1) + 1 )
```

at every cut `1 ≤ j < d` and every `W ≥ 2`. This is a single-decomposition
statement (MPS / tensor-train).

E1.9 and E6.3 cover individual additional cases: 2D `φ(x, a)` rank,
DCT / wavelet sparsity. Each was closed in isolation. There is no
family-level closure that says *the entire class of polynomial-bond-dim
tensor-network ansätze fails* on χ_P.

The leverage: most "what if we use ANSATZ-X" proposals reduce to a single
unfolding-rank check. A clean meta-theorem reducing the entire family to
E2.1 closes dozens of single proposals at once.

---

## Conjecture N1 (precise)

Let `T : [W]^d → ℚ` be the base-`W` reshape of `χ_P`, with `N := W^d`.
Let `χ` be the bond / multilinear / canonical / hierarchical-rank
parameter of one of the following tensor-network ansätze:

1. **MPS / Tensor-Train (TT)** with bond dimension `χ`.
2. **Tucker decomposition** with multilinear rank tuple `(r_1, ..., r_d)`,
   `χ := max_j r_j`.
3. **Hierarchical Tucker (HT)** along a balanced binary tree with bond
   dimension `χ` at every internal node.
4. **CP / canonical polyadic** with rank `χ`.
5. **Tensor Ring (TR)** with bond dimension `χ`.
6. **MERA** (1D, binary) with bond dimension `χ` at every layer.
7. **PEPS** (any 2D reshape into `(W^a, W^a)` for `a = d/2`) with bond
   dimension `χ`.

For each ansatz to **exactly represent** `T`, the bond / rank / multilinear
parameter `χ` must satisfy:

```
   χ  ≥  φ(W) · W^(d-j-1) + 1     for some cut j with 1 ≤ j < d,
```

equivalently `χ ≥ Ω(N^(1 - log_W W/φ(W)))`, which is `Ω(N^c)` with
`c = c(W) > 0` for every fixed `W ≥ 2`.

In particular: there is no `(W, ansatz)` pair from the seven above with
bond dimension `polylog(N)` representing `χ_P` exactly.

The **mechanism** is uniform across all seven cases: each ansatz admits
an *unfolding contraction* whose matricisation-rank is `≥ χ` (or
`≥ χ^{O(1)}` for MERA / PEPS), and this matricisation-rank is bounded
below by E2.1.

---

## Reduction sketch (one-line proofs per ansatz)

| Ansatz | Where the unfolding bound enters | Lower bound on `χ` |
|--------|----------------------------------|--------------------|
| MPS / TT | bond `j` between sites `j` and `j+1` IS the rank of the `(W^j × W^(d-j))` unfolding (definitional) | `≥ φ(W)·W^(d-j-1)+1` for every `j` |
| Tucker | mode-`j` multilinear rank `r_j` IS the rank of the mode-`j` matricisation (definitional) | `r_j ≥ φ(W)·W^(d-j-1)+1` for every `j` (note: mode-`j` matricisation is `(W × W^(d-1))`, rank ≥ rank of the balanced `j` cut after contraction over the other modes; standard bound `r_j ≥ rank(M^(j))` does not apply — see "Subtlety" below) |
| HT | every tree-cut induces a matricisation; the bond dim at that internal node IS the matricisation rank | `≥ φ(W)·W^(d-j-1)+1` for tree cut at depth `j` |
| CP | Kruskal's bound: `CP-rank(T) ≥ max_j rank(M^(j))` | `≥ max_j (φ(W)·W^(d-j-1)+1)` |
| TR | a single bond cut induces a `(W^j × W^(d-j))` matricisation up to a transposition through the loop; rank is the same | `≥ φ(W)·W^(d-j-1)+1` |
| MERA | Renyi-2 entanglement entropy across any cut is `log(rank(M^(j)))`; MERA bond dim `χ` bounds this by `O(log χ)` per layer, with `O(log N)` layers, so `log(rank) ≤ O(polylog) · log(χ)` ⇒ `χ ≥ rank^{1/polylog}` | `≥ N^Ω(1/polylog)` |
| PEPS | any 2D reshape `(W^a × W^a)` is a single matricisation; rank is `≥ φ(W)·W^(d-a-1)+1` for `a = d/2`; PEPS bond dim across the corresponding boundary cut is the matricisation rank | `≥ Ω(N^{1/2})` |

### Subtlety: mode-j Tucker rank vs. cut-j unfolding rank

Tucker's mode-`j` matricisation has shape `(W × W^(d-1))`, and its rank
is `≤ W` trivially. So Tucker's *multilinear rank tuple* per definition
is bounded by `(W, W, ..., W)` and is NOT directly compressible information.
What IS bounded by E2.1 is the **block-mode** matricisation along a
contiguous prefix `{1, ..., j}` vs suffix `{j+1, ..., d}`. Tucker's
multilinear rank is uninformative for χ_P's unfolding-rank lower bound.

**Therefore Tucker is NOT closed by the unfolding bound alone.** Tucker
has total parameter count `prod r_j + sum r_j · W = O(W^d) = O(N)` in the
worst case — it never compresses `χ_P` non-trivially because `χ_P` has
support of cardinality ≈ `N / log N` (rough density of primes), so the
core tensor `G ∈ [W]^d` has at least that many non-zero entries.

So the Tucker case is closed by a different (trivial) argument: a tensor
with `Θ(N / log N)` non-zero entries cannot be Tucker-compressed to fewer
parameters than its support cardinality.

---

## Falsification criterion (pre-stated)

The conjecture is falsified if any of the following empirical or
theoretical observations hold at the tested scales:

- **F1** (MPS / TT): for some `(W, d, j)` with `2 ≤ W ≤ 5`, `4 ≤ d ≤ 12`,
  `1 ≤ j < d`, the actual MPS bond dim `< φ(W)·W^(d-j-1) + 1`.
- **F2** (HT): for some balanced binary tree on `d ≤ 12` modes and
  some tree cut, the matricisation rank is strictly less than the
  corresponding unfolding rank.
- **F3** (CP): for any `(W, d)` tested, the numerical CP-rank estimate
  (via ALS or Kruskal) is `< max_j rank(M^(j))`.
- **F4** (Tensor Ring): the TR bond dim at any cut is `<` MPS bond dim
  at the same cut.
- **F5** (MERA): for some `(W, d)`, a polylog-bond-dim MERA represents
  `χ_P` exactly. (We test this as: the rank of the unfolding is
  achievable by a MERA with `χ = polylog`; checked indirectly via
  parameter count.)
- **F6** (PEPS): the 2D reshape rank is asymptotically smaller than
  the 1D unfolding rank.

The conjecture **passes** iff none of F1–F6 fires at the tested
`(W, d)` grid.

---

## Why this is novel

E2.1 is a single-ansatz statement (MPS only). Prior CLOSED_PATHS rows
close MPS, automata-theoretic, J-fraction, persistent-homology, and
free-probability decompositions one at a time — each row is a single
sentence "this also fails" without a common mechanism. **The N1 closure
extracts a single mechanism (unfolding-rank lower bound from E2.1)
that subsumes all polynomial-spatial-locality tensor ansätze**.

This is the type of family-level closure that EDGES.md sections E7.10
(AKS family) and E5.8 (Brandt family) achieve for their respective
domains. A successful N1 produces:

- A precise theorem statement covering 7 ansätze.
- A unifying reduction in 7 lines.
- An empirical verification at `(W, d)` up to `(5, 10)`.
- A new CLOSED_PATHS row covering the family.
- An annotation to E2.1 promoting it from MPS-only to all-spatial-ansatz.

---

## Why this is B-grade, not A-grade

This is a **unification of existing closures under a shared mechanism** —
exactly the pattern CLAUDE.md flags as "B-grade at best". An A-grade
result would be a *positive* construction: a polylog-bond-dim
representation that DOES work for some new ansatz, or a reduction
showing that a positive complexity bound on `χ_P` representation
implies an algorithmic improvement.

The N1 closure cannot deliver a polylog algorithm — by construction it
forecloses one. Its value is in narrowing the search space.

---

## File layout

- `definition.md` — this file (precise statement, falsification
  criterion).
- `tensor_compression_family.py` — verification script.
- `tensor_compression_family_results.md` — outcomes + closure outcome.
- `tensor_compression_family_results.json` — raw numerical outputs.
