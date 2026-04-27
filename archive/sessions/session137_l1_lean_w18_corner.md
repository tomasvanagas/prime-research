# Session 137 — L1 Lean Formalisation, W=18 corner

**Mode:** Lean formalisation (Arc 2 / L1).
**Target:** Extend the orthogonal-corner closure of `mps_bond_dim`
(Route A^{(9)}) to a new wheel `W ≥ 13`.

## What this session produced

A new sorry-free Lean theorem
`mps_bond_dim_W_eq_18_d_eq_j_plus_1 : ∀ j ≥ 1, (unfolding 18 (j+1) j).rank = 7`,
plus its supporting prime exhibit
`exists_invertible_submatrix_W_eq_18_d_eq_j_plus_1` and five new
`chiP_X_eq_one` helpers for `X ∈ {29, 43, 179, 211, 293}`. All type-check;
`#print axioms` returns only `[propext, Classical.choice, Quot.sound]`.

This is the **ninth unconditional `mps_bond_dim` instance** (`W ∈ {2, 3,
4, 5, 6, 8, 12, 18}` for the orthogonal corner `d = j + 1`, plus W=2's
non-orthogonal corner `j = 1`), the **seventh instance over a wheel
`W ≥ 3`**, the **fourth instance using `det_of_upperTriangular`** (after
W=5, W=8, W=12), and **the first instance with `R = 7`**.

## Cross-domain content

None. The argument uses only mathlib's matrix/linear algebra +
elementary number theory (primality of 7 explicit primes, decidability
of non-primality of 21 explicit composites). The novel mathematical
content is:

1. **The triangulation.** Permutation `ρ ↦ (0, 2, 9, 1, 11, 6, 16)` and
   `σ ↦ (1, 6, 16, 10, 12, 0, 4)` upper-triangularises the relevant
   `7 × 7` submatrix of the `18^j × 18` unfolding. Diagonal primes are
   `{2, 43, 179, 29, 211, 109, 293}`; below-diagonal entries are 21
   composites. This was found by Python pre-search and verified line-by-
   line in Lean.

2. **The W=14 obstruction (negative shape edge).** Pre-search confirmed
   that W=14 (also `R = 7`) admits no leading-row triangulation with
   `ρ < 14`: rows 2 and 5 of the W=14 j=1 slab have identical support
   pattern `(1, 1, 0, 1, 0, 1)` at the seven chosen columns. W=14
   joins `W ∈ {7, 9, 10, 11}` in the "needs `det_of_blockTriangular`"
   set. This was a non-trivial observation since W=14 looked at first
   like the natural single-session next target.

## Self-evaluation against the four questions

### 1. What did I produce that was not in the project before this session?

- A new Lean theorem `mps_bond_dim_W_eq_18_d_eq_j_plus_1` (sorry-free,
  no new axioms).
- A new `exists_invertible_submatrix_W_eq_18_d_eq_j_plus_1`.
- Five new `chiP_X_eq_one` helpers for new primes `{29, 43, 179, 211,
  293}`.
- The first use of `Fin.prod_univ_seven` in this file.
- The first use of `norm_num` (rather than `decide`) for primality of
  large numbers (179, 211, 293) in this file — a recursion-depth lesson
  for future agents working with primes ≥ 150.
- A new negative-shape observation: W=14 is structurally obstructed
  for leading-row triangulation (now the fifth W in the
  needs-`det_of_blockTriangular` set).

### 2. What edges did my work compose or cite?

- **E2.1** (the MPS bond-dim identity itself) — the formalisation
  closes the W=18 corner of E2.1 unconditionally. The rest of E2.1's
  parameter space (general `(W, d, j)`) remains gated by the prime-
  density existential.
- The proof composes the file-internal lemma `upper_bound` (cited for
  `rank ≤ 7` since `rank_le_width` only gives `rank ≤ 18`).

### 3. If my session produced only duplicate closures, why?

It did not. The W=18 corner is genuinely new mathematical content (no
prior session formalised any instance of `mps_bond_dim` for a wheel
with `R = 7`), and the W=14 obstruction is a new structural negative
shape that future agents should not waste a session attempting.

### 4. Next-action for the next agent

In `RESEARCH_AGENDA.md` Arc 2 and `NOVELTY_CHALLENGES.md` §3 L1:

* (a) **Route A^{(10)}**: develop `Matrix.det_of_blockTriangular`
  technique and use it to close W=9 (the cleanest of the
  obstructed set). Multi-session.
* (b) **Route A^{(11)}**: extend the leading-row triangulation regime
  to W=20 via `Fin.prod_univ_nine`. Pre-search at S137 found a clean
  triangulation (`ρ ↦ (0, 2, 9, 14, 12, 1, 7, 16, 10)`,
  `σ ↦ (1, 6, 18, 12, 0, 2, 8, 16, 10)`, diagonal primes `{2, 47, 199,
  293, 241, 23, 149, 337, 211}`). **Caveat**: mathlib provides
  `Fin.prod_univ_X` only up to `X = 8`; W=20 needs an extra manual
  `Fin.prod_univ_succ` step. Single-session viable but ~700 Lean
  lines.
* (c) **Route C (PNT for low-density)**: closes the easy regime but not
  the saturating one. Single-session, ambitious.

The honest `Hoheisel`-grade Route A is unchanged; it remains beyond a
single session and likely requires a separate sub-arc of its own.

## Self-grade

**B** (substantive refinement). 

This is *not* an A-grade session: the W=18 corner is a continuation of
the pattern S98–S129 already established. It exhibits a ninth instance
of an existing technique. The only marginally novel content is:

(i) the discovery that W=14 — the smallest wheel with `R = 7` — is
    structurally obstructed (a small negative-shape edge), and
(ii) the recursion-depth lesson for primality at scale ≥ 150 (using
     `norm_num` instead of `decide`).

Per CLAUDE.md, "Lean translation of an already-informally-proven
argument, with the translation type-checking but introducing no new
mathematical content" is C-grade. This session's W=18 closure has new
*Lean* content (the explicit triangulation, `Fin.prod_univ_seven`
usage, and the `norm_num` workaround), but the mathematical content
(rank of MPS unfolding for W=18) is corollary of E2.1 already verified
empirically. So **B** rather than A — the session advances the formal
catalogue but does not produce frontier results.

The honest path forward via either Route A^{(10)} (block-triangular,
multi-session) or Route A^{(11)} (W=20 via manual `prod_univ_succ`) is
recorded in the next-action.

## Files touched

- `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`
  — added ~520 lines (W=18 block).
- `experiments/formalisations/E2_1_mps_bond_dim/mps_bond_dim_notes.md`
  — appended W=18 progress; added Route A^{(9)} entry.
- `RESEARCH_AGENDA.md` — Arc 2 milestone added, run/session header
  updated.
- `NOVELTY_CHALLENGES.md` — L1 entry updated with S137 progress + next
  action.

## Build verification

```
$ cd experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim
$ lake build
✔ [8314/8315] Built MPSBondDim (6.2s)
Build completed successfully (8315 jobs).

$ # axiom check (run via lake env lean):
'E2_1.mps_bond_dim_W_eq_18_d_eq_j_plus_1' depends on axioms:
  [propext, Classical.choice, Quot.sound]
'E2_1.exists_invertible_submatrix_W_eq_18_d_eq_j_plus_1' depends on axioms:
  [propext, Classical.choice, Quot.sound]
```

The single `sorry` (in `exists_invertible_submatrix`, the general case)
is unchanged — it remains the prime-density obligation gated by
Hoheisel-type density results not yet in mathlib.
