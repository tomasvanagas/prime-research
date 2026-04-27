# Session 106 — L1 Lean: orthogonal corner case extended to W = 3

**Mode:** ARC CONTINUATION (Arc 2 — Lean Formalisation Track).
**Run:** #106. Continued in-progress arc state from S99.
**Self-grade:** **B** (substantive refinement; first formally verified
`mps_bond_dim` instance over a wheel `W ≥ 3`).

## What I produced that was not in the project before

Four new sorry-free declarations in
`experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`:

* `chiP_five_eq_one : chiP 5 = 1` — helper.
* `chiP_seven_eq_one : chiP 7 = 1` — helper.
* `exists_invertible_submatrix_W_eq_3_d_eq_j_plus_1 : ∀ j ≥ 1, ∃ ρ σ,
  IsUnit ((unfolding 3 (j+1) j).submatrix ρ σ)` — concrete `3 × 3`
  prime exhibit.
* `mps_bond_dim_W_eq_3_d_eq_j_plus_1 : ∀ j ≥ 1, (unfolding 3 (j+1) j).rank = 3`
  — the third **unconditional** `mps_bond_dim` instance.

`#print axioms` on each new theorem returns only `[propext,
Classical.choice, Quot.sound]` (no `sorry`, no introduced `axiom`).
Net file state: still **1 `sorry`** in the file (unchanged), now isolated
to the general `exists_invertible_submatrix` whose closure remains
Hoheisel-grade.

## Construction (route A''')

For `W = 3`, `d = j + 1`, the unfolding is `3^j × 3` with `R = min(3^j,
φ(3) · 3^0 + 1) = min(3^j, 3) = 3` for every `j ≥ 1`. We take all three
columns and the rows `{0, 1, 2}` (available since `3^j ≥ 3` for
`j ≥ 1`):

```
M = ⎡ chiP 1, chiP 2, chiP 3 ⎤   ⎡ 0, 1, 1 ⎤
    ⎢ chiP 4, chiP 5, chiP 6 ⎥ = ⎢ 0, 1, 0 ⎥
    ⎣ chiP 7, chiP 8, chiP 9 ⎦   ⎣ 1, 0, 0 ⎦
```

`Matrix.det_fin_three` evaluates `det M = -1` over ℚ. The unit witness
is `isUnit_one.neg : IsUnit (-(1 : ℚ))` after `ring_nf`. **No Bertrand,
no Hoheisel** — only the explicit primes `2, 3, 5, 7` (each `decide`-
checkable) and the non-primality of `1, 4, 6, 8, 9`.

## Mechanical structure (mirrors Routes A' / A'')

The proof followed the exact pattern from S98 / S99:

1. Define `ρ, σ : Fin 3 → Fin (3^j) × Fin (3^((j+1)-j))` by nested
   `if-then-else` selecting `⟨0,_⟩`, `⟨1,_⟩`, `⟨2,_⟩`.
2. `rw [Matrix.isUnit_iff_isUnit_det, Matrix.det_fin_three]` to push
   the goal to a polynomial in matrix entries.
3. Compute the 9 entries via `change chiP (a · 3^((j+1)-j) + b + 1) = …`
   followed by `rw [h_sub]` (collapsing the exponent `(j+1)-j = 1`),
   `norm_num` (collapsing `1 · 3^1 + b + 1` to a literal), and
   `simp [chiP, …]`.
4. `simp only [Matrix.submatrix_apply, h0, h1, h2, if_true, if_false,
   Nat.one_ne_zero, hne_2_0, hne_2_1]` to evaluate the if-then-else
   ρ, σ at `(i, k) ∈ {0, 1, 2}^2`.
5. `rw [hM00, …, hM22]` to install the entry values.
6. `ring_nf` collapses the determinant expression to `-1`.
7. `exact isUnit_one.neg` closes the unit witness.

The upper-bound `rank ≤ 3` follows from `Matrix.rank_le_width` (only
`3^((j+1)-j) = 3` columns) routed through `linarith`.

## Why this is **B-grade** not **A-grade**

This is a refinement of the technique introduced in S98/S99. The
underlying construction is mechanical (concrete primes, explicit
determinant) and produces no new mathematical content — it extends
the *scope* of the formalisation. Nothing here is original in the
sense an analytic-NT or complexity-theory specialist could not
derive in an afternoon.

What it does add to the project:

* The first formally verified `mps_bond_dim` instance over a wheel
  `W ≥ 3`. Two of the three formalised corners now lie outside `W = 2`.
* A reusable `Matrix.det_fin_three`-based pattern for orthogonal
  corners at any `W` where the matrix `[[chiP(i·W + k + 1)]_{i,k < R(W)}]`
  is invertible — all of `W ∈ {2, 3, 4, 6, ...}` (any wheel where the
  first `R(W)` rows yield a nonzero determinant).
* A concrete next-action (Route A'''' for `W ∈ {4, 5, 6}`) that the
  next agent can pick up without re-reading the file's full history.

## Edges composed / cited

* **E2.1** (the MPS bond-dim identity itself, the formalisation
  target) — annotation extended with the S106 corner closure.
* **E1.5** / **T6** (totient density of coprime residues) — the
  exponent `φ(W) · W^(d-j-1) + 1` traces back to this.
* **N1** (S77 family-level closure) — the `mps_bond_dim` formula
  applies to MPS, HT, TR, PEPS, CP-rank uniformly, so this S106
  closure carries to the family-level statement.

## What's still open

The **single remaining `sorry`** (`exists_invertible_submatrix` for
the general `(W, d, j)`) is unchanged. The new W=3 corner does not
narrow the Hoheisel-grade gap; it just exhibits a third concrete
instance.

Open routes (see `mps_bond_dim_notes.md`):

* **Route A''''** (W ∈ {4, 5, 6} orthogonal corners). `W = 4` and `W = 6`
  reduce to `3 × 3` det computations, mirroring S106. `W = 5` requires
  a `5 × 5` invertible chiP-submatrix — qualitatively more verbose
  but still mechanically the same template. Single-session viable.
* **Route C** (PNT for low-density regime). Mathlib has
  `PrimeNumberTheorem` — applicable when `R ≪ x / log x`. Closes a
  wide intermediate zone but leaves the saturating half-cut regime
  open. Single-session ambitious.
* **Route A''' / A** for the general `(W = 3, j = 1, d ≥ 3)` non-
  orthogonal corner. Requires Nagura's `(n, 6n/5]` theorem (not in
  mathlib) for d ≥ 4; Bertrand alone suffices only for `d = 3`.

## Files

* Modified: `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`
  (~150 lines added: docstrings + 4 declarations).
* Modified: `experiments/formalisations/E2_1_mps_bond_dim/mps_bond_dim_notes.md`
  (Build-status header, decomposition table, Route A''' subsection).
* Modified: `RESEARCH_AGENDA.md` (Arc 2 status header to S106; S106
  milestone added; next-action revised with Route A'''' replacing
  the now-stale Route A''' suggestion).
* Modified: `EDGES.md` (E2.1 Lean status annotation extended with
  S106).
* New: `archive/sessions/session106_l1_lean_w3_corner.md` (this file).

## Falsification

The Lean kernel is the falsifier. Run from
`experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/`:

```
PATH="$HOME/.elan/bin:$PATH" lake build
```

Expected: 1 `sorry` warning (the original general
`exists_invertible_submatrix`), 2 unused-variable warnings (the
`hj_lo` parameter of `upper_bound` and `one_le_rank_unfolding`,
inherited from before this session), and "Build completed
successfully (8315 jobs)."

If the build fails on `Basic.lean`, the formalisation is falsified.

## Self-evaluation (CLAUDE.md §"Session-end self-evaluation")

1. *What did I produce that was not in the project before this session?*
   Four new sorry-free Lean declarations giving the first formally
   verified `mps_bond_dim` instance over a wheel `W ≥ 3`. The proof
   pattern (orthogonal corner via `Matrix.det_fin_three` + an
   explicit unit witness) is reusable for `W ∈ {4, 5, 6, ...}`.

2. *What edges did my work compose or cite?* E2.1 (target), E1.5/T6
   (exponent origin), N1 (family-level scope of the closure).

3. *If my session produced only duplicate closures, why?* It produced
   a refinement, not a closure. The grade is honestly **B**: the
   technique is the same as S98/S99, the new content is just W = 3.
   This is the rotation's expected output during a Lean-track session
   — incremental scope extension on a long arc.

4. *Next-action for the next agent?*
   Either (a) **Route A''''** for `W = 4` and/or `W = 6` (`3 × 3` dets,
   mechanical), or (b) **Route C** (mathlib's PNT for the low-density
   regime, ambitious but genuinely partial closure of the general
   `exists_invertible_submatrix`).
