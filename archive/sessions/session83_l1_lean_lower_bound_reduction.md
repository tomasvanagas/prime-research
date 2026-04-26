# Session 83 — L1 Lean: closing `lower_bound` via prime-exhibit reduction

**Mode:** Lean formalisation (Arc 2, L1).
**Target:** advance `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`
toward a complete machine-checked proof of E2.1 (MPS bond-dim identity).

## Going-in state (post-S76)

`lake build` succeeded with 1 `sorry`, sitting on the `lower_bound`
lemma:
```
theorem lower_bound (...) :
    min (W^j) (φ(W)·W^(d-j-1) + 1) ≤ (unfolding W d j).rank := by
  sorry
```
All other declarations (`row_support_coprime`, `live_columns_count`,
`upper_bound`, `rank_le_min_dim`, `mps_bond_dim`) closed `sorry`-free.
The informal proof in `novel/mps_bond_dimension.md` hand-waves the
lower bound: "exhibit `R` rows whose restriction to the live columns is
linearly independent, via prime-counting density of base-`W` blocks."
The rank-bound *logic* (how an exhibit yields the inequality) and the
prime-density *content* were not separated.

## Session work

### Decision

Close the rank-bound logic completely. Isolate the prime-density
content into a single new declaration. Two sub-decisions drove this:

* **Why not attempt the full proof?** The prime-density argument
  requires either Bertrand-type prime existence in shrinking intervals
  (`[i·W^(d-j)+1, (i+1)·W^(d-j)]` for every `0 ≤ i < R`) plus a
  residue-class dovetail, or a generic Vandermonde exhibit over a
  finite extension. Both are 100-300-line Lean efforts, not
  single-session work.
* **Why split?** The structural reduction is the cleanest possible
  factoring: it converts a rank inequality (a high-level statement
  about linear algebra) into a pure existence statement about an
  invertible submatrix. The next agent inherits a much sharper
  obligation than "prove `lower_bound`."

### Decomposition

Introduce a new theorem:
```
theorem exists_invertible_submatrix (W d j : ℕ)
    (hW : 2 ≤ W) (hj_lo : 1 ≤ j) (hj_hi : j < d) :
    ∃ (ρ : Fin (min (W^j) (φ(W)·W^(d-j-1)+1)) → Fin (W^j))
      (σ : Fin (min (W^j) (φ(W)·W^(d-j-1)+1)) → Fin (W^(d-j))),
      IsUnit ((unfolding W d j).submatrix ρ σ) := by
  sorry
```

Use it to close `lower_bound`:
```
theorem lower_bound (...) :
    min (W^j) (φ(W)·W^(d-j-1) + 1) ≤ (unfolding W d j).rank := by
  classical
  obtain ⟨ρ, σ, hUnit⟩ := exists_invertible_submatrix W d j hW hj_lo hj_hi
  have h_eq :
      ((unfolding W d j).submatrix ρ σ).rank =
        min (W^j) (φ(W)·W^(d-j-1) + 1) := by
    have h := Matrix.rank_of_isUnit ((unfolding W d j).submatrix ρ σ) hUnit
    rw [Fintype.card_fin] at h
    exact h
  calc min (W^j) (φ(W)·W^(d-j-1) + 1)
      = ((unfolding W d j).submatrix ρ σ).rank := h_eq.symm
    _ ≤ (unfolding W d j).rank := Matrix.rank_submatrix_le _ _ _
```

Two mathlib lemmas do the heavy lifting:

* `Matrix.rank_of_isUnit : IsUnit A → A.rank = Fintype.card n` (square `n × n`).
* `Matrix.rank_submatrix_le : (A.submatrix r c).rank ≤ A.rank`.

The proof is six lines and contains no `sorry`. Build:
```
$ lake build
⚠ [8313/8315] Built MPSBondDim.Basic (10s)
warning: MPSBondDim/Basic.lean:398:8: declaration uses `sorry`
✔ [8314/8315] Built MPSBondDim (5.0s)
Build completed successfully (8315 jobs).
```
The single `sorry` warning is on `exists_invertible_submatrix` —
exactly as designed.

### Documentation

Updated:
- `experiments/formalisations/E2_1_mps_bond_dim/mps_bond_dim_notes.md`
  — table of declarations now lists 8 entries; status marker for
  `lower_bound` flipped from `sorry` to **done**; new "What the next
  session needs to do" describes the Bertrand / Vandermonde routes for
  closing `exists_invertible_submatrix`.
- `RESEARCH_AGENDA.md` Arc 2 — milestone for `lower_bound` checked off
  with the structural reduction note; new milestone for
  `exists_invertible_submatrix` added with the two routes.
- `NOVELTY_CHALLENGES.md` §3 L1 — S83 progress block describes the
  structural reduction. §0 backup-leverage line points at
  `exists_invertible_submatrix` as the new closure target.
- `EDGES.md` E2.1 — annotation block recording the Lean formalisation
  status: rank-bound logic closed, prime-existence exhibit remains.

## Output artefact

* `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`
  — type-checks; 1 `sorry` (was 1, but the `sorry` *moved* from a rank
  inequality to a pure prime-density existential).
* No new `axiom`.
* No new mathlib dependencies beyond what S76 already pulled in.

## Self-evaluation

### Did this session produce something not in the project before?

Yes:

* The structural reduction `lower_bound ← exists_invertible_submatrix`
  is now closed in Lean. Six lines of provably-correct rank-bound
  logic that did not exist before.
* The clean factoring "rank-bound logic / prime-density content" did
  not exist before either; both `novel/mps_bond_dimension.md` and the
  earlier Lean skeleton intermixed them. Future sessions inherit a
  much sharper exhibit-existential obligation.
* Two specific routes (A: Bertrand-type; B: Vandermonde) are now
  documented as concrete next-actions.

### Edges composed or cited

E2.1 (MPS bond-dim identity) — formalisation track, target of L1.

### Why not more?

The prime-density existential is genuinely hard (100-300 Lean lines for
either route). Single-session attempts on it would have produced an
incomplete proof or required pulling in mathlib's APs / Bertrand
machinery without time to debug. The structural reduction was the
high-leverage piece available in this session: it converts a sub-task
("prove `lower_bound`") into a strictly cleaner sub-task ("prove
`exists_invertible_submatrix`") without mixing concerns.

### Next-action for the next agent

Prove `exists_invertible_submatrix` in
`experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`.
Two routes outlined in the file's docstring and in
`mps_bond_dim_notes.md`:
* **Route A (Bertrand-style).** For each row `0 ≤ i < R` use
  `Nat.bertrand` to find a prime in `[i·W^(d-j)+1, 2·i·W^(d-j)]` and
  intersect with the row's window. Combine with Dirichlet's theorem on
  primes in arithmetic progressions to dovetail the chosen primes
  across distinct residue classes mod `W^(d-j)`.
* **Route B (generic Vandermonde).** Replace the prime-density appeal
  by a Vandermonde determinant exhibit over a finite extension of ℚ;
  this bypasses arithmetic and reduces to a pure linear-algebra fact.

Once `exists_invertible_submatrix` is closed, `lower_bound` and
`mps_bond_dim` close automatically without further edits.

## Self-graded letter

**B-grade.**

This is substantive refinement of an in-progress arc — the cleanest
possible factoring of the remaining work, with the rank-bound *logic*
fully closed in Lean and the prime-density *content* isolated for
future formalisation. It is not A-grade because the cumulative L1 Lean
track still has a `sorry` (so the "≥ 50 lines, no sorry, no axiom"
A-grade rule for E2.1 is not yet met). It is not C-grade because the
work introduced a concrete piece of provably-correct Lean
(`lower_bound`'s body) that did not exist before, plus a sharply-stated
new theorem signature that is the cleanest possible obligation for the
remaining gap.

The grade reflects honest mid-arc progress: the project's Lean track
has moved one structural step closer to A-grade closure of E2.1, but
that closure is not yet realised.
