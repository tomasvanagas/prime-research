# Session 75 — L1 Lean Formalisation: `live_columns_count` closed

## Mission

Continue Arc 2 (Lean Formalisation Track). The next-action from S72 was
to prove `live_columns_count` in
`experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`,
the CRT periodicity count
`#{k ∈ [0, W^(d-j)) : gcd(k+1, W) = 1} = φ(W) · W^(d-j-1)`.

## Outcome

- **`live_columns_count` is fully proved.** `lake build` succeeds with
  no `sorry` warning on this declaration. 0 new `axiom` introductions.
- File state: 3 `sorry`s remain (`mps_bond_dim`, `upper_bound`,
  `lower_bound`), down from 4 before this session.
- Repaired a previously-broken proof attempt of `live_columns_count`
  that was committed but did not type-check (4 distinct compile errors:
  `Fin.sum_univ_eq_sum_range` rewrite mismatch, `Nat.count_add` type
  mismatch, two `omega` failures inside an `Ico` bijection where the
  un-beta-reduced `Finset.card_bij'` lambda was hiding the goal).
- The `live_columns_count` proof is ~110 lines, structured in three
  steps:
  1. Stage A: `Finset.card_bij` with the value-projection
     `(k : Fin (W^(d-j))) ↦ k.val` shows the count over `Fin (W^(d-j))`
     equals the count over `Finset.range (W^(d-j))`.
  2. Stage B: induction on `M`, the multi-block count
     `|{n ∈ range(W·M) : gcd(n+1,W)=1}| = M · φ(W)`. The successor
     step splits the range into `range(W·M) ∪ Ico(W·M)(W·M+W)`
     (disjoint), invokes the IH on the first piece, and reduces the
     second piece to `Nat.filter_coprime_Ico_eq_totient W 1` via the
     shift bijection `n ↔ n + 1 − W·M` (`Finset.card_bij'`).
  3. Combine: `W^(d-j) = W · W^(d-j-1)` via single-occurrence
     `conv_lhs` rewrite (avoiding subterm clash with `d-j-1`),
     instantiate Stage B at `M := W^(d-j-1)`.
- Notes file (`mps_bond_dim_notes.md`) and Arc 2 milestones updated.
  Next-action explicitly handed off to `upper_bound`.

## Self-evaluation

**1. What did I produce that was not in the project before this session?**

A fully verified Lean 4 proof of the periodicity count

> `#{k ∈ Fin (W^(d-j)) : gcd(k.val + 1, W) = 1} = φ(W) · W^(d-j-1)`

(for `W ≥ 2` and `1 ≤ d - j`). This is the third Lean lemma in the L1
proof tree to type-check without `sorry`, and the crux of the upcoming
`upper_bound` argument. Before this session this proof did not exist
(the file held a commit-time attempt that did not compile).

The cleaned-up proof exposes a small reusable pattern: the bijection
`n ↔ n + 1 − W·M` shift between `Ico(W·M)(W·M+W)` and `Ico 1 (1+W)`,
combined with `Nat.gcd_add_mul_right_left` to relate gcds across
shifts. This same identity drives `row_support_coprime`, so two of the
three completed lemmas in this file now share their key arithmetic
ingredient.

**2. What edges did my work compose or cite?**

- **E2.1** (target of formalisation).
- Mathlib citations: `Finset.card_bij`, `Finset.card_bij'`,
  `Finset.filter_union`, `Finset.card_union_of_disjoint`,
  `Finset.disjoint_left`, `Nat.filter_coprime_Ico_eq_totient`,
  `Nat.gcd_add_mul_right_left`, `Nat.gcd_comm`, `Nat.sub_add_cancel`,
  `pow_succ`. All from the standard mathlib totient + Finset
  infrastructure — no new dependencies.

**3. If my session produced only duplicate closures, why?**

N/A — the session produced a new machine-verified proof.

That said, this is not an A-grade session by the CLAUDE.md bar
("Lean 4 proof of a non-trivial theorem (≥ 50 lines of Lean content,
no `sorry`)"). The lemma `live_columns_count` *individually* meets the
50-line threshold (it is ~110 lines), but it is a sub-lemma of a
larger result still gated by `sorry`s. The honest grade is **B**:
substantive refinement (closing one of three remaining placeholders in
L1), with the technical contribution being the structural pattern
shared with `row_support_coprime`. It is not "ambitious failure" —
it is a tractable next-step that landed.

A failure-mode worth recording: the original attempt at this proof
(committed pre-S75) silently did not type-check. The four errors
clustered around (i) a tactic chain that pretended to do
"Fin → Nat.count" via `card_eq_sum_ones` + `Fin.sum_univ_eq_sum_range`
rewrites, which do not unify cleanly because of the lambda-argument
form; (ii) `Nat.count_add` being applied with implicit `p` that the
elaborator could not infer; and (iii) `omega` failing under
un-beta-reduced `Finset.card_bij'` lambdas. The repair was structural
(switch from `Nat.count` to direct `Finset.card_bij` + manual block
splitting; force beta-reduction with `change` before each membership
omega goal) rather than incremental.

**4. What is the next-action for the next agent?**

Prove `upper_bound` in `MPSBondDim/Basic.lean`. With
`row_support_coprime` and `live_columns_count` both closed, this
reduces to a single linear-algebra step: rows `i ≥ 1` of the
`unfolding W d j` matrix lie in the
`φ(W)·W^(d-j-1)`-dimensional subspace
`V := {v : Fin(W^(d-j)) → ℚ | ∀ k, gcd(k+1,W) ≠ 1 → v k = 0}`
(by `row_support_coprime`); `dim V = φ(W) · W^(d-j-1)` (by
`live_columns_count`); and adding row `0` increases the rank by at
most one. Mathlib lemmas to consider: `Matrix.rank_le`,
`Submodule.finrank_le_finrank_of_le`, `Matrix.rank_le_card_height`.

After `upper_bound`, the harder remaining piece is `lower_bound` —
the prime-counting-density argument that may require a different
proof strategy in Lean (probably via `PrimeCounting` infrastructure
or an explicit construction of `min(W^j, φ(W)·W^(d-j-1) + 1)` linearly
independent rows).

## Key decisions

- **Switch the proof strategy from `Nat.count` to direct `Finset.card_bij`.**
  The pre-S75 attempt used `Nat.count` and `Function.Periodic` infrastructure.
  This was conceptually clean but had three implementation pitfalls:
  unifying lambda-argument forms across `count_eq_card_filter_range`
  and `Fin.sum_univ_eq_sum_range`, applying `Nat.count_add` with
  implicit predicate, and the `omega` failure on un-beta-reduced
  bijection lambdas. Switching to direct Finset bijections via
  `Finset.card_bij` and `Finset.card_bij'` flattens the proof: each
  step has an explicit Finset on each side and the bijection is
  given as a function on natural numbers.

- **Use `change` (not `show`) to force beta-reduction in `card_bij'`
  membership goals.** The mathlib style linter warns against using
  `show` to alter the goal (it should only be used for readability of
  unchanged goals). The four `change` rewrites in this proof are the
  beta-reductions needed because `Finset.card_bij'` ambient lambdas
  block `omega` from solving membership inequalities.

- **Use `conv_lhs` for the `W^(d-j) = W · W^(d-j-1)` rewrite.** The
  naive `rw [← Nat.sub_add_cancel hdj]` rewrites *every* occurrence of
  `d - j`, including inside `d - j - 1` (which parses as `(d - j) - 1`,
  so `d - j` is a subterm). The result is a goal of the form
  `W^((d-j-1)+1) = W * W^((d-j-1)+1 - 1)` which subsequent rewrites
  do not resolve. `conv_lhs` restricts the rewrite to the LHS only.

## Falsification

The lemma is verified by the Lean 4 kernel. `lake build` from
`experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/` succeeds
with 8315 jobs and emits warnings only on the three remaining `sorry`s
(`mps_bond_dim`, `upper_bound`, `lower_bound`). Any future change that
breaks `live_columns_count` will fail the build.

## Grade self-assessment

**B-grade.** Substantive refinement of an existing arc — closed one of
three remaining `sorry`s in L1 and, in passing, established a shared
arithmetic ingredient (`Nat.gcd_add_mul_right_left`) across two of the
three closed lemmas. Not A-grade because the closure is a sub-lemma of
a result still gated by `sorry`s; not C-grade because the artifact is
a new machine-verified ~110-line theorem that did not exist in the
project before (the prior attempt did not compile, so it provided no
verification).
