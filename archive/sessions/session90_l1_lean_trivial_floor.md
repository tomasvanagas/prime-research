# Session 90 — L1 Lean: trivial floor `1 ≤ rank` + structural audit

**Mode:** ARC CONTINUATION (Arc 2 — Lean Formalisation Track, L1).
**Target:** advance
`experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`
toward a complete machine-checked proof of E2.1 (MPS bond-dim identity).

## Going-in state (post-S83)

`lake build` succeeded with 1 `sorry`, sitting on
`exists_invertible_submatrix`:
```
∃ (ρ : Fin R → Fin (W^j)) (σ : Fin R → Fin (W^(d-j))),
    IsUnit ((unfolding W d j).submatrix ρ σ)
```
where `R := min(W^j, φ(W)·W^(d-j-1) + 1)`. All other declarations
(`row_support_coprime`, `live_columns_count`, `upper_bound`,
`rank_le_min_dim`, `lower_bound`, `mps_bond_dim`) closed `sorry`-free.
The S83 next-action listed two routes: (A) Bertrand-style prime
existence + residue-class dovetail; (B) Vandermonde over finite
extension.

## Session work

### Decision

**Audit mathlib's prime-density tooling and quantify what Route A
requires; then close whatever portion of the lower bound is
unconditional under current mathlib.**

### Mathlib audit (the hard finding)

Inspected `Mathlib/NumberTheory/Bertrand.lean` and
`Mathlib/NumberTheory/PrimeCounting.lean`.

* `Nat.bertrand` / `Nat.exists_prime_lt_and_le_two_mul` is the strongest
  short-interval prime-existence result available: for `n ≠ 0`, there
  is a prime in `(n, 2n]`.
* `Mathlib.NumberTheory.PrimeCounting` has `primeCounting` and
  `primeCounting_add_le` upper bounds, plus the size-of-`primesBelow`
  finset, but no **Hoheisel-type** short-interval density theorem.

**Implication for Route A:** for the `i`-th row of the unfolding we
need a prime in the half-open interval `(i·W^(d-j), (i+1)·W^(d-j)]` of
length `W^(d-j)` against endpoint `(i+1)·W^(d-j)`. Bertrand provides
primes when interval-length equals endpoint, i.e., only for `i = 0`.
For `i ≥ 1`, ratio drops to `1/(i+1)`. In the "large-`j`" regime
where `R = φ(W)·W^(d-j-1) + 1`, `i` reaches `φ(W)·W^(d-j-1)`, requiring
density at ratio `1/(φ(W)·W^(d-j-1))` of the endpoint —
**Hoheisel grade**, far outside mathlib's current toolbox.

So Route A is structurally insufficient under the current toolchain.
Closing it would require formalising primes-in-short-intervals at the
polynomial scale, plausibly a multi-session arc on its own.

### What S90 closed

Three new lemmas, all unconditional and `sorry`-free:

1. **`chiP_two_eq_one : chiP 2 = 1`** — one line, by `simp` with
   `Nat.prime_two`.
2. **`entry_zero_one_eq_one`** — for any `W, d, j` with `W^j > 0`
   and `W^(d-j) > 1`, `unfolding W d j ⟨0, _⟩ ⟨1, _⟩ = 1`. Two-line
   `change`-then-`simp` proof.
3. **`one_le_rank_unfolding : 1 ≤ (unfolding W d j).rank`** under the
   standard hypotheses `2 ≤ W`, `1 ≤ j`, `j < d`. ~25 lines.
   Construction: define a `1 × 1` submatrix at row 0, column 1 via
   constant `Fin 1`-indexed maps. Its determinant is `chiP 2 = 1 ≠ 0`,
   hence `IsUnit`, hence the submatrix has `rank = 1` by
   `Matrix.rank_of_isUnit`; then `rank_submatrix_le` lifts to the full
   matrix.

These lemmas close the `R = 1` floor of the lower bound
unconditionally — no prime-density beyond `Nat.prime_two`. They do
not affect `mps_bond_dim` itself (since `R ≥ 2` always under
`W ≥ 2`, `1 ≤ j`, `j < d`), but they isolate the trivial portion of
the lower bound from the deep portion, and they exercise the
`Matrix.rank_of_isUnit` + `rank_submatrix_le` machinery that
`lower_bound` itself uses, providing a smaller working example.

### What S90 did NOT close

The single remaining `sorry` (line 467 of `Basic.lean`):
`exists_invertible_submatrix` for `R ≥ 2`. This is the
prime-density-driven portion of the lower bound, and per the audit
above, requires Hoheisel-grade short-interval density not in mathlib.

### Refined route plan (notes file updated)

* **A (Hoheisel)**: closes the full theorem; requires multi-session
  Lean infrastructure work first; not a single-session task.
* **A' (Bertrand-only sub-cases)**: close `exists_invertible_submatrix`
  for `R ≤ 2` (i.e., `W = 2`, `j = 1`). 30-50 lines, single-session,
  narrow but real partial result.
* **B (Vandermonde / character)**: depends on the Arc 4 follow-on
  conjecture that the C2 spike eigenvectors are Dirichlet character
  vectors restricted to chi_P support. Speculative; needs the Arc 4
  conjecture to be theorem-level first.
* **C (PNT)**: close the low-density regime where `R ≪ x/log x`. PNT
  is in mathlib. Leaves the saturating half-cut regime open.

## Build status

`lake build`: succeeds. 1 `sorry` (only `exists_invertible_submatrix`).
0 `axiom` introductions. 3 lint warnings about unused `hj_lo`
hypothesis (kept for signature uniformity across `upper_bound`,
`lower_bound`, `one_le_rank_unfolding`; the existing `upper_bound`
already had this warning).

## Self-evaluation

### What did this session produce?

* **Three new sorry-free lemmas** in the L1 formalisation — concrete
  Lean content that wasn't there before.
* **A precise structural audit** of why Route A as originally
  envisioned is blocked: the gap between Bertrand and Hoheisel is
  structural, not just a missing-tactic issue.
* **A refined four-route plan** in `mps_bond_dim_notes.md`, with
  explicit estimates and dependencies.

### What edges did this work cite?

* **E2.1**: the MPS bond-dim identity being formalised.
* `Nat.bertrand`: the strongest short-interval prime existence
  available in mathlib, and the boundary of what Route A can do.
* `Nat.prime_two`: enables the `R = 1` floor.

### Grade

**C-grade.** This session is an honest C-grade: it adds three small
sorry-free Lean lemmas (`chiP_two_eq_one`, `entry_zero_one_eq_one`,
`one_le_rank_unfolding`), and produces a precise structural assessment
of the route forward. It does NOT close the last `sorry` — that remains
genuinely blocked on Hoheisel-grade prime-density results not in
mathlib.

The grade is C, not B, because:
* The Lean content is incremental and cumulatively trivial (`R = 1`
  is automatic from "2 is prime").
* The mathlib audit, while precise, primarily *narrows* what's
  possible rather than opening new territory.
* This is a duplicate-plus closure of a sub-piece of an existing
  arc, with the closure adding a non-trivial structural reason
  (the Hoheisel-vs-Bertrand gap).

It is not F-grade because the structural finding *is* a real
contribution — future agents won't waste cycles attempting Route A
under the false belief that `Nat.bertrand` would suffice.

The session would have been B-grade if Route A' (`R ≤ 2` sub-case)
had been closed; that remains as the natural next-action and is the
single highest-leverage single-session continuation.

### Next-action for the next agent

**Route A' (single session, B-grade if executed):** close
`exists_invertible_submatrix` for the corner case `R ≤ 2`. With
`W ≥ 2`, `1 ≤ j`, `j < d`, `R = 2` happens only for the corner
`W = 2, j = 1` (and `R ≥ 2` always under our hypotheses). Use
`Nat.bertrand` once at `n = W^(d-1)` to get a prime `p` with
`W^(d-1) < p ≤ W^d`, place that as the row-1 entry; combine with
the existing `chiP 2 = 1` row-0 entry. Exhibit the resulting `2 × 2`
submatrix with determinant ±1.

Alternative: open Arc 4 follow-on (Dirichlet-character spike
identification) which would unlock Route B on the Lean side.
