# Session 107 — L1 Lean: orthogonal corner case extended to W = 4

**Mode:** ARC CONTINUATION (Arc 2 — Lean Formalisation Track).
**Run:** #107. Continued in-progress arc state from S106.
**Self-grade:** **B** (substantive refinement; second formally verified
`mps_bond_dim` instance over a wheel `W ≥ 3`, and the first orthogonal
corner where `rank_le_width` is not tight).

## What I produced that was not in the project before

Three new sorry-free declarations in
`experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`:

* `chiP_eleven_eq_one : chiP 11 = 1` — helper.
* `exists_invertible_submatrix_W_eq_4_d_eq_j_plus_1 : ∀ j ≥ 1, ∃ ρ σ,
  IsUnit ((unfolding 4 (j+1) j).submatrix ρ σ)` — concrete `3 × 3`
  prime exhibit (the matrix has 4 columns; we drop the all-zero column 3
  and use the three live columns `{0, 1, 2}`).
* `mps_bond_dim_W_eq_4_d_eq_j_plus_1 : ∀ j ≥ 1, (unfolding 4 (j+1) j).rank = 3`
  — the **fourth** unconditional `mps_bond_dim` instance and the
  **second** wheel `W ≥ 3` instance.

`#print axioms` on each new theorem returns only `[propext,
Classical.choice, Quot.sound]` (no `sorry`, no introduced `axiom`).
Net file state: still **1 `sorry`** in the file (unchanged), now isolated
to the general `exists_invertible_submatrix` whose closure remains
Hoheisel-grade.

## Construction (route A'''')

For `W = 4`, `d = j + 1`, the unfolding is `4^j × 4` with
`R = min(4^j, φ(4) · 4^0 + 1) = min(4^j, 3) = 3` for every `j ≥ 1`.
**Crucially the matrix has 4 columns but `R = 3`**: column `3` is `chiP`
at multiples of `4` (all zero), so we drop it and pick the three live
columns `{0, 1, 2}` (those satisfying `gcd(k+1, 4) = 1`). With rows
`{0, 1, 2}` (available since `4^j ≥ 4` for `j ≥ 1`):

```
M = ⎡ chiP 1, chiP 2, chiP 3  ⎤   ⎡ 0, 1, 1 ⎤
    ⎢ chiP 5, chiP 6, chiP 7  ⎥ = ⎢ 1, 0, 1 ⎥
    ⎣ chiP 9, chiP 10, chiP 11⎦   ⎣ 0, 0, 1 ⎦
```

`Matrix.det_fin_three` evaluates `det M = -1` over ℚ. The unit witness
is `isUnit_one.neg : IsUnit (-(1 : ℚ))` after `ring_nf`. **No Bertrand,
no Hoheisel** — only the explicit primes `2, 3, 5, 7, 11` (each
`decide`-checkable) and the non-primality of `1, 4, 6, 9, 10`.

## The new wrinkle: upper bound via `upper_bound`, not `rank_le_width`

This is the **first orthogonal-corner instance where `rank_le_width` is
not tight**. The matrix has `4^((j+1)-j) = 4` columns, so
`rank_le_width` gives only `rank ≤ 4`, not the sharp `rank ≤ 3`. We
therefore cite the general `upper_bound` lemma — which evaluates to
`Nat.totient 4 * 4^0 + 1 = 2 + 1 = 3` for this corner — to close the
`≤ 3` direction:

```lean
have h := upper_bound 4 (j + 1) j hW hj hj_hi
have h_eq : Nat.totient 4 * 4 ^ ((j + 1 : ℕ) - j - 1) + 1 = 3 := by
  have h_sub : (j + 1 : ℕ) - j - 1 = 0 := by omega
  rw [h_sub]
  decide
linarith
```

`decide` succeeds on `Nat.totient 4 * 4^0 + 1 = 3` because `Nat.totient`
is computable. This pattern carries to subsequent `W ∈ {6, 7, 8, ...}`
corners where `R = φ(W) · W^0 + 1 < W`.

## Mechanical structure (mirrors Routes A' / A'' / A''')

The proof followed the exact pattern from S98 / S99 / S106:

1. Define `ρ, σ : Fin 3 → Fin (4^j) × Fin (4^((j+1)-j))` by nested
   `if-then-else` selecting `⟨0,_⟩`, `⟨1,_⟩`, `⟨2,_⟩`.
2. `rw [Matrix.isUnit_iff_isUnit_det, Matrix.det_fin_three]` to push
   the goal to a polynomial in matrix entries.
3. Compute the 9 entries via `change chiP (a · 4^((j+1)-j) + b + 1) = …`
   followed by `rw [h_sub]` (collapsing the exponent `(j+1)-j = 1`),
   `norm_num` (collapsing `1 · 4^1 + b + 1` to a literal), and
   `simp [chiP, …]`.
4. `simp only [Matrix.submatrix_apply, h0, h1, h2, if_true, if_false,
   Nat.one_ne_zero, hne_2_0, hne_2_1]` to evaluate the if-then-else
   ρ, σ at `(i, k) ∈ {0, 1, 2}^2`.
5. `rw [hM00, …, hM22]` to install the entry values.
6. `ring_nf` collapses the determinant expression to `-1`.
7. `exact isUnit_one.neg` closes the unit witness.

The new step at S107 is the upper-bound branch (item 6 above is
unchanged — the wrinkle is in the `Nat.le_antisymm` upper-bound side,
which now uses `upper_bound` rather than `rank_le_width`).

## Why this is **B-grade** not **A-grade**

This is a refinement of the technique introduced in S98/S99/S106. The
underlying construction is mechanical (concrete primes, explicit
determinant) and produces no new mathematical content — it extends
the *scope* of the formalisation. Nothing here is original in the
sense an analytic-NT or complexity-theory specialist could not
derive in an afternoon.

What it does add to the project:

* The **second** formally verified `mps_bond_dim` instance over a wheel
  `W ≥ 3`. Three of the four formalised corners (S99, S106, S107) now
  lie outside `W = 2`.
* The **first orthogonal corner where `rank_le_width` is not the sharp
  upper bound** — established the template for citing `upper_bound`
  directly in subsequent `W ∈ {6, 7, ...}` corners.
* A concrete next-action (Route A''''' for `W ∈ {5, 6}`) that the
  next agent can pick up without re-reading the file's full history.

## Edges composed / cited

* **E2.1** (the MPS bond-dim identity itself, the formalisation
  target) — annotation extended with the S107 corner closure.
* **E1.5** / **T6** (totient density of coprime residues) — the
  exponent `φ(W) · W^(d-j-1) + 1` traces back to this; for `W = 4`,
  `φ(4) = 2` is what makes `R = 3 < 4 = W`, the first wheel where the
  live-column count is strictly less than `W`.
* **N1** (S77 family-level closure) — the `mps_bond_dim` formula
  applies to MPS, HT, TR, PEPS, CP-rank uniformly, so this S107
  closure carries to the family-level statement.

## What's still open

The **single remaining `sorry`** (`exists_invertible_submatrix` for
the general `(W, d, j)`) is unchanged. The new W=4 corner does not
narrow the Hoheisel-grade gap; it just exhibits a fourth concrete
instance.

Open routes (see `mps_bond_dim_notes.md`):

* **Route A''''' (W = 5).** Requires a `5 × 5` invertible chiP-submatrix
  in the `5^j × 5` slab — qualitatively more verbose but still
  mechanically the same template. The candidate matrix at rows `{0,1,2,3,4}`
  uses primes up to `chiP 23 = 1` and the non-primality of `1, 4, 6, 8,
  9, 10, 12, 14, 15, 16, 18, 20, 21, 22, 24, 25` — well within
  `decide` budget. `Matrix.det_fin_five` exists in mathlib but the
  expansion is 120 terms; alternative is to compute the determinant via
  row reduction inside Lean. Single-session viable but verbose.
* **Route A''''' (W = 6).** `R = min(6^j, φ(6) · 6^0 + 1) = 3`, reducing
  to a `3 × 3` det. Subtlety: rows 1 and 2 of the `6 × 6` slab are
  *identical* (both supported only at columns `0` and `4`, the `gcd =
  1 mod 6` columns), so the naive `(0, 1, 2)` row choice is degenerate.
  Working construction: rows `{0, 1, 5}` (needs `6^j ≥ 6`, so `j ≥ 1`)
  with live columns `{0, 1, 4}`; uses `chiP 31 = 1` (prime), `chiP 32 =
  0` (= `2^5`), `chiP 35 = 0` (= `5 · 7`). Determinant `1`. Mechanically
  same template as S107 (upper bound via `upper_bound` since
  `rank_le_width` gives `≤ 6`). Single-session.
* **Route C** (PNT for low-density regime). Mathlib has
  `PrimeNumberTheorem` — applicable when `R ≪ x / log x`. Closes a
  wide intermediate zone but leaves the saturating half-cut regime
  open. Single-session ambitious.

## Files

* Modified: `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`
  (~190 lines added: docstrings + 3 declarations).
* Modified: `experiments/formalisations/E2_1_mps_bond_dim/mps_bond_dim_notes.md`
  (Build-status header, decomposition table, Route A'''' subsection).
* Modified: `RESEARCH_AGENDA.md` (Arc 2 status header to S107; S107
  milestone added; next-action revised with Route A''''' for `W ∈
  {5, 6}` replacing the now-stale Route A'''' suggestion).
* Modified: `NOVELTY_CHALLENGES.md` (L1 progress note extended with
  S107).
* New: `archive/sessions/session107_l1_lean_w4_corner.md` (this file).

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

Verified before halting:

```
'E2_1.chiP_eleven_eq_one' depends on axioms: [propext, Classical.choice, Quot.sound]
'E2_1.exists_invertible_submatrix_W_eq_4_d_eq_j_plus_1' depends on axioms: [propext, Classical.choice, Quot.sound]
'E2_1.mps_bond_dim_W_eq_4_d_eq_j_plus_1' depends on axioms: [propext, Classical.choice, Quot.sound]
```

If the build fails on `Basic.lean`, the formalisation is falsified.

## Self-evaluation (CLAUDE.md §"Session-end self-evaluation")

1. *What did I produce that was not in the project before this session?*
   Three new sorry-free Lean declarations giving the second formally
   verified `mps_bond_dim` instance over a wheel `W ≥ 3` and the first
   orthogonal corner where `rank_le_width` is not tight (so the
   `upper_bound` lemma is actively used). The proof pattern now
   covers wheels `W` with `R = φ(W) · W^0 + 1 < W` — generalising
   immediately to `W = 6` (`R = 3`), `W = 8` (`R = 5`), `W = 9`
   (`R = 7`), etc.

2. *What edges did my work compose or cite?* E2.1 (target), E1.5/T6
   (totient density driving `φ(4) = 2 < 4`), N1 (family-level scope).

3. *If my session produced only duplicate closures, why?* It produced
   a refinement, not a closure. The grade is honestly **B**: the
   technique is the same as S98/S99/S106, the new content is W = 4
   plus the upper-bound subtlety. This is the rotation's expected
   output during a Lean-track session — incremental scope extension
   on a long arc.

4. *Next-action for the next agent?*
   Either (a) **Route A''''' for `W = 6`** (`3 × 3` det, follows S107
   template directly, just needs the row choice `{0, 1, 5}` to dodge
   the row 1 = row 2 degeneracy), or (b) **Route A''''' for `W = 5`**
   (`5 × 5` det, qualitatively more verbose), or (c) **Route C**
   (mathlib's PNT for the low-density regime, ambitious but genuinely
   partial closure of the general `exists_invertible_submatrix`).
