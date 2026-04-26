/-
# E2.1: MPS bond dimension of the prime indicator (statement-only)

Formalisation of the theorem in `novel/mps_bond_dimension.md`:

  For chi_P : [1, W^d] -> {0,1} the prime indicator, viewed as a tensor of
  shape (W, W, ..., W) (d copies) by base-W reshape, for every cut
  1 <= j < d and every integer W >= 2:

      rank M^(j) = min ( W^j ,  phi(W) * W^(d-j-1) + 1 )

  where phi is Euler's totient and M^(j) is the (W^j x W^(d-j)) unfolding.

Stage 1: theorem statement only. The proof is a `sorry`. The next stage
splits the proof into upper and lower bounds, each via separate lemmas.
-/

import Mathlib

open Matrix Nat

namespace E2_1

/--
The prime indicator, embedded in `ℚ`. We use `ℚ` (a field) because matrix
`rank` is most informative over a field.
-/
noncomputable def chiP (n : ℕ) : ℚ := if Nat.Prime n then 1 else 0

/--
The `j`-th unfolding of the base-`W` reshape of the prime indicator.
Rows are indexed by `Fin (W^j)`, columns by `Fin (W^(d-j))`.
The `(i, k)` entry is `chiP (i * W^(d-j) + k + 1)`.

The `+ 1` mirrors the informal convention that the prime indicator domain
is `[1, W^d]`. Internally the matrix entries cover `n = 1, 2, ..., W^d`.
-/
noncomputable def unfolding (W d j : ℕ) :
    Matrix (Fin (W ^ j)) (Fin (W ^ (d - j))) ℚ :=
  fun i k => chiP (i.val * W ^ (d - j) + k.val + 1)

/-!
## Decomposition of the proof

The main theorem `mps_bond_dim` (stated and proved at the bottom of this
file) splits into two inequalities:

* **Upper bound:** `rank ≤ φ(W)·W^(d-j-1) + 1` (`upper_bound`).
* **Trivial ceiling:** `rank ≤ W^j` (`rank_le_min_dim`), giving
  `rank ≤ min(W^j, φ(W)·W^(d-j-1) + 1)`.
* **Lower bound:** `min(W^j, φ(W)·W^(d-j-1) + 1) ≤ rank` (`lower_bound`).

Putting them together with `Nat.le_antisymm` yields the main equality.
The lemmas come first; the main theorem closes the file.
-/

/-!
### Lemmas feeding the upper bound

The upper bound is proved by two independent observations:

* For every row `i ≥ 1`, the only columns `k` where the entry can be
  nonzero are those with `gcd(k + 1, W) = 1` (`row_support_coprime`).
* The number of such columns in `[0, W^(d-j))` is exactly
  `φ(W) * W^(d-j-1)` (`live_columns_count`).

Combining: rows `i ≥ 1` span a subspace of dimension at most
`φ(W) * W^(d-j-1)`. Row `i = 0` adds at most one further dimension.
-/

/--
Row-support lemma. For every `i ≥ 1`, every `k` such that the matrix
entry is nonzero (i.e. `chiP (i * W^(d-j) + k + 1) = 1`), the column
index satisfies `gcd(k + 1, W) = 1`.

The proof: nonzero entry ⇒ `n = i * W^(d-j) + k + 1` is prime; for
`i ≥ 1` and `j < d` we have `n > W^(d-j) ≥ W`, so any prime factor of `W`
cannot equal `n` and thus cannot divide `n`; hence `gcd(n, W) = 1`, and
since `n ≡ k + 1 (mod W)` (because `W | i * W^(d-j)` for `j < d`),
`gcd(k + 1, W) = 1`.
-/
theorem row_support_coprime
    (W d j : ℕ) (hW : 2 ≤ W) (hj_hi : j < d)
    (i : Fin (W ^ j)) (hi : 1 ≤ i.val)
    (k : Fin (W ^ (d - j))) (hentry : unfolding W d j i k ≠ 0) :
    Nat.gcd (k.val + 1) W = 1 := by
  -- We need d - j ≥ 1.
  have hdj : 1 ≤ d - j := Nat.sub_pos_of_lt hj_hi
  have hWpos : 1 ≤ W := le_trans (by norm_num) hW
  -- The matrix entry is `chiP n` where n = i * W^(d-j) + k + 1.
  -- `hentry` says this is nonzero, hence n is prime.
  have hprime : Nat.Prime (i.val * W ^ (d - j) + k.val + 1) := by
    by_contra h
    apply hentry
    change chiP (i.val * W ^ (d - j) + k.val + 1) = 0
    simp [chiP, h]
  -- W ≤ W^(d-j) and so W ≤ i * W^(d-j) for i ≥ 1, giving n > W.
  have hWdj : W ≤ W ^ (d - j) := by
    have := Nat.pow_le_pow_right hWpos hdj
    simpa using this
  have hn_gt : W < i.val * W ^ (d - j) + k.val + 1 := by
    have h1 : W ^ (d - j) ≤ i.val * W ^ (d - j) := by
      have : 1 * W ^ (d - j) ≤ i.val * W ^ (d - j) :=
        Nat.mul_le_mul_right _ hi
      simpa using this
    omega
  -- gcd(n, W) = 1: a prime n > W can share no divisor > 1 with W.
  have hgcd_n :
      Nat.Coprime (i.val * W ^ (d - j) + k.val + 1) W := by
    rcases Nat.coprime_or_dvd_of_prime hprime W with hco | hdvd
    · exact hco
    · exfalso
      have := Nat.le_of_dvd hWpos hdvd
      omega
  -- Now reduce mod W: rewrite the expression so the multiple of W is on
  -- the right, then use `gcd_add_mul_right_left` to peel it off.
  have hpow_split : W ^ (d - j) = W ^ (d - j - 1) * W := by
    have hsucc : (d - j - 1) + 1 = d - j := Nat.sub_add_cancel hdj
    calc W ^ (d - j)
        = W ^ ((d - j - 1) + 1) := by rw [hsucc]
      _ = W ^ (d - j - 1) * W := pow_succ W _
  have hreshape :
      i.val * W ^ (d - j) + k.val + 1
        = (k.val + 1) + (i.val * W ^ (d - j - 1)) * W := by
    -- Direct Nat identity, using `hpow_split : W^(d-j) = W^(d-j-1) * W`.
    have := hpow_split
    nlinarith [hpow_split, i.val, k.val]
  -- Apply the gcd identity.
  have : Nat.gcd (i.val * W ^ (d - j) + k.val + 1) W
            = Nat.gcd (k.val + 1) W := by
    rw [hreshape, Nat.gcd_add_mul_right_left]
  -- And `Coprime` is by definition `gcd … = 1`.
  exact this ▸ hgcd_n

/--
CRT-based count. The number of integers `k ∈ [0, W^(d-j))` with
`gcd(k + 1, W) = 1` is exactly `φ(W) * W^(d-j-1)`.

Proof idea: by periodicity in blocks of `W` consecutive `k`, each block
contributes exactly `φ(W)` admissible values. There are `W^(d-j-1)`
such blocks (using `1 ≤ d - j`).
-/
theorem live_columns_count
    (W d j : ℕ) (hW : 2 ≤ W) (hj_hi : j < d) :
    ((Finset.univ : Finset (Fin (W ^ (d - j)))).filter
        (fun k => Nat.gcd (k.val + 1) W = 1)).card =
      Nat.totient W * W ^ (d - j - 1) := by
  have hdj : 1 ≤ d - j := Nat.sub_pos_of_lt hj_hi
  have hWpos : 0 < W := by omega
  -- Step 1: convert the Fin-indexed filter to a `Finset.range`-indexed filter
  -- via the value-injection bijection `k ↦ k.val`.
  have step_fin_to_range :
      ((Finset.univ : Finset (Fin (W ^ (d - j)))).filter
          (fun k => Nat.gcd (k.val + 1) W = 1)).card =
        ((Finset.range (W ^ (d - j))).filter
          (fun n => Nat.gcd (n + 1) W = 1)).card := by
    refine Finset.card_bij (fun (k : Fin _) _ => k.val) ?_ ?_ ?_
    · intros k hk
      rw [Finset.mem_filter] at hk
      rw [Finset.mem_filter, Finset.mem_range]
      exact ⟨k.isLt, hk.2⟩
    · intros k1 _ k2 _ h
      exact Fin.ext h
    · intros n hn
      rw [Finset.mem_filter, Finset.mem_range] at hn
      refine ⟨⟨n, hn.1⟩, ?_, rfl⟩
      rw [Finset.mem_filter]
      exact ⟨Finset.mem_univ _, hn.2⟩
  rw [step_fin_to_range]
  -- Step 2: prove the multi-block count formula by induction on the number
  -- of W-blocks. For every M:
  --   |{n ∈ range(W * M) : gcd(n+1, W) = 1}| = M * φ(W).
  have multi_block : ∀ M : ℕ,
      ((Finset.range (W * M)).filter
          (fun n => Nat.gcd (n + 1) W = 1)).card = M * Nat.totient W := by
    intro M
    induction M with
    | zero => simp
    | succ M ih =>
      -- Split range(W*(M+1)) = range(W*M) ∪ Ico(W*M)(W*M+W).
      have hsplit_range : Finset.range (W * (M + 1)) =
          Finset.range (W * M) ∪ Finset.Ico (W * M) (W * M + W) := by
        ext n
        simp only [Finset.mem_union, Finset.mem_range, Finset.mem_Ico]
        have : W * (M + 1) = W * M + W := by ring
        omega
      have hdisj : Disjoint (Finset.range (W * M))
          (Finset.Ico (W * M) (W * M + W)) := by
        rw [Finset.disjoint_left]
        intros n hn1 hn2
        rw [Finset.mem_range] at hn1
        rw [Finset.mem_Ico] at hn2
        omega
      rw [hsplit_range, Finset.filter_union,
          Finset.card_union_of_disjoint
            (hdisj.mono (Finset.filter_subset _ _) (Finset.filter_subset _ _)),
          ih]
      -- block_count: count over Ico(W*M)(W*M+W) of gcd(n+1, W) = 1 equals φ(W).
      -- Bijection to ((Ico 1 (1+W)).filter (Nat.Coprime W)) via n ↔ n + 1 - W*M.
      have block_count : ((Finset.Ico (W * M) (W * M + W)).filter
          (fun n => Nat.gcd (n + 1) W = 1)).card = Nat.totient W := by
        rw [← Nat.filter_coprime_Ico_eq_totient W 1]
        refine Finset.card_bij' (fun (n : ℕ) _ => n + 1 - W * M)
            (fun (m : ℕ) _ => m + W * M - 1) ?_ ?_ ?_ ?_
        · -- forward maps into target filter
          intros n hn
          rw [Finset.mem_filter, Finset.mem_Ico] at hn
          change n + 1 - W * M ∈
            (Finset.Ico 1 (1 + W)).filter (fun x => Nat.Coprime W x)
          rw [Finset.mem_filter, Finset.mem_Ico]
          have hWM_le_n : W * M ≤ n := hn.1.1
          have hn_lt : n < W * M + W := hn.1.2
          refine ⟨⟨by omega, by omega⟩, ?_⟩
          -- Goal: Nat.Coprime W (n + 1 - W * M).
          show Nat.Coprime W (n + 1 - W * M)
          have h_eq : n + 1 = (n + 1 - W * M) + M * W := by
            rw [Nat.mul_comm M W]; omega
          have h_gcd : Nat.gcd (n + 1) W = Nat.gcd (n + 1 - W * M) W := by
            conv_lhs => rw [h_eq]
            rw [Nat.gcd_add_mul_right_left]
          unfold Nat.Coprime
          rw [Nat.gcd_comm, ← h_gcd]
          exact hn.2
        · -- inverse maps into source filter
          intros m hm
          rw [Finset.mem_filter, Finset.mem_Ico] at hm
          change m + W * M - 1 ∈
            (Finset.Ico (W * M) (W * M + W)).filter (fun x => Nat.gcd (x + 1) W = 1)
          rw [Finset.mem_filter, Finset.mem_Ico]
          have hm_lo : 1 ≤ m := hm.1.1
          have hm_hi : m < 1 + W := hm.1.2
          refine ⟨⟨by omega, by omega⟩, ?_⟩
          -- Goal: Nat.gcd ((m + W * M - 1) + 1) W = 1.
          show Nat.gcd ((m + W * M - 1) + 1) W = 1
          have h1 : m + W * M - 1 + 1 = m + W * M := by omega
          rw [h1]
          have h2 : m + W * M = m + M * W := by ring
          rw [h2, Nat.gcd_add_mul_right_left]
          have hCop : Nat.Coprime W m := hm.2
          rw [Nat.gcd_comm]
          exact hCop
        · -- left_inv: ((n + 1) - W*M) + W*M - 1 = n
          intros n hn
          rw [Finset.mem_filter, Finset.mem_Ico] at hn
          change n + 1 - W * M + W * M - 1 = n
          have : W * M ≤ n := hn.1.1
          omega
        · -- right_inv: (m + W*M - 1) + 1 - W*M = m
          intros m hm
          rw [Finset.mem_filter, Finset.mem_Ico] at hm
          change m + W * M - 1 + 1 - W * M = m
          have : 1 ≤ m := hm.1.1
          omega
      rw [block_count]
      ring
  -- Step 3: apply multi_block with M = W^(d-j-1).
  have hpow_split : W ^ (d - j) = W * W ^ (d - j - 1) := by
    conv_lhs => rw [show (d - j) = (d - j - 1) + 1 from
                      (Nat.sub_add_cancel hdj).symm]
    rw [pow_succ, mul_comm]
  rw [hpow_split, multi_block]
  ring

/--
Upper bound. Every row `i ≥ 1` is supported on the `φ(W) * W^(d-j-1)`
columns whose index `k` satisfies `gcd(k + 1, W) = 1`, because for such
rows `n = i * W^(d-j) + k + 1 > W` and primality forces `gcd(n, W) = 1`.
Row `i = 0` adds at most one further linearly independent direction.
-/
theorem upper_bound
    (W d j : ℕ) (hW : 2 ≤ W) (hj_lo : 1 ≤ j) (hj_hi : j < d) :
    (unfolding W d j).rank ≤ Nat.totient W * W ^ (d - j - 1) + 1 := by
  classical
  have hWpos : 0 < W := by omega
  have hWj_pos : 0 < W ^ j := Nat.pow_pos hWpos
  -- Index 0 in `Fin (W^j)`.
  let i₀ : Fin (W ^ j) := ⟨0, hWj_pos⟩
  -- Row-0 unit vector in `Fin (W^j) → ℚ`.
  let e0 : Fin (W ^ j) → ℚ := Pi.single i₀ (1 : ℚ)
  -- The "live" / good column index set.
  let GoodCols : Finset (Fin (W ^ (d - j))) :=
    (Finset.univ : Finset _).filter (fun k => Nat.gcd (k.val + 1) W = 1)
  have hGoodCard : GoodCols.card = Nat.totient W * W ^ (d - j - 1) :=
    live_columns_count W d j hW hj_hi
  -- Generating set: row-0 unit vector together with the good-column images.
  let S : Finset (Fin (W ^ j) → ℚ) :=
    insert e0 (GoodCols.image (unfolding W d j).col)
  have hS_card : S.card ≤ Nat.totient W * W ^ (d - j - 1) + 1 := by
    have h1 : S.card ≤ (GoodCols.image (unfolding W d j).col).card + 1 :=
      Finset.card_insert_le _ _
    have h2 : (GoodCols.image (unfolding W d j).col).card ≤ GoodCols.card :=
      Finset.card_image_le
    omega
  -- For every "bad" column k (`gcd(k+1, W) ≠ 1`) and every row index `i ≠ i₀`,
  -- the matrix entry vanishes — by `row_support_coprime`.
  have hcol_zero_outside :
      ∀ (k : Fin (W ^ (d - j))) (_ : Nat.gcd (k.val + 1) W ≠ 1)
        (i : Fin (W ^ j)) (_ : i ≠ i₀),
        (unfolding W d j).col k i = 0 := by
    intro k hk i hi
    have hi_pos : 1 ≤ i.val := by
      rcases Nat.eq_zero_or_pos i.val with hi0 | hipos
      · exact (hi (Fin.ext hi0)).elim
      · exact hipos
    change (unfolding W d j) i k = 0
    by_contra hne
    exact hk (row_support_coprime W d j hW hj_hi i hi_pos k hne)
  -- A bad column is a scalar multiple of `e0` (only the i = i₀ entry survives).
  have hbad_col_form :
      ∀ (k : Fin (W ^ (d - j))) (_ : Nat.gcd (k.val + 1) W ≠ 1),
        (unfolding W d j).col k = ((unfolding W d j) i₀ k) • e0 := by
    intro k hk
    funext i
    by_cases hi : i = i₀
    · subst hi
      change (unfolding W d j) i₀ k = ((unfolding W d j) i₀ k) * e0 i₀
      simp [e0]
    · have h_lhs : (unfolding W d j).col k i = 0 := hcol_zero_outside k hk i hi
      have h_rhs : (((unfolding W d j) i₀ k) • e0) i = 0 := by
        change ((unfolding W d j) i₀ k) * e0 i = 0
        simp [e0, hi]
      rw [h_lhs, h_rhs]
  -- Every column lies in `Submodule.span ℚ (S : Set _)`.
  have hcol_in_span : ∀ k : Fin (W ^ (d - j)),
      (unfolding W d j).col k ∈ Submodule.span ℚ (S : Set (Fin (W ^ j) → ℚ)) := by
    intro k
    by_cases hk : Nat.gcd (k.val + 1) W = 1
    · -- Good column: in S directly.
      apply Submodule.subset_span
      change (unfolding W d j).col k ∈ (S : Set _)
      have hmem : (unfolding W d j).col k ∈ S := by
        apply Finset.mem_insert_of_mem
        refine Finset.mem_image.mpr ⟨k, ?_, rfl⟩
        exact Finset.mem_filter.mpr ⟨Finset.mem_univ _, hk⟩
      exact_mod_cast hmem
    · -- Bad column: scalar multiple of e0 ∈ S.
      rw [hbad_col_form k hk]
      apply Submodule.smul_mem
      apply Submodule.subset_span
      change e0 ∈ (S : Set _)
      exact_mod_cast (Finset.mem_insert_self e0 _)
  -- Conclude: rank ≤ |S| ≤ φ(W) · W^(d-j-1) + 1.
  rw [Matrix.rank_eq_finrank_span_cols]
  have h_le : Submodule.span ℚ (Set.range (unfolding W d j).col)
            ≤ Submodule.span ℚ (S : Set (Fin (W ^ j) → ℚ)) := by
    rw [Submodule.span_le]
    rintro v ⟨k, rfl⟩
    exact hcol_in_span k
  calc Module.finrank ℚ (Submodule.span ℚ (Set.range (unfolding W d j).col))
      ≤ Module.finrank ℚ (Submodule.span ℚ (S : Set (Fin (W ^ j) → ℚ))) :=
        Submodule.finrank_mono h_le
    _ ≤ S.card := finrank_span_finset_le_card S
    _ ≤ Nat.totient W * W ^ (d - j - 1) + 1 := hS_card

/-!
### The lower bound, reduced to a prime exhibit

The lower bound is the harder direction: every measured `(W, d, j)`
saturates the upper bound, and the informal argument
(`novel/mps_bond_dimension.md`) hand-waves over a prime-counting density
to exhibit `R := min(W^j, φ(W)·W^(d-j-1) + 1)` linearly independent rows.

We isolate that content into a single existential lemma —
`exists_invertible_submatrix` — which asserts the existence of an
`R × R` invertible submatrix of the unfolding. The structural reduction
from this exhibit to `lower_bound` is mechanical (mathlib's
`Matrix.rank_of_isUnit` plus `Matrix.rank_submatrix_le`) and is closed
unconditionally below. The exhibit itself is the only remaining `sorry`
in the formalisation; it captures the prime-density content of the
informal proof.
-/

/--
**Prime-exhibit existence.** There exist indexings `ρ : Fin R → Fin (W^j)`
and `σ : Fin R → Fin (W^(d-j))` (where `R = min(W^j, φ(W)·W^(d-j-1)+1)`)
such that the `R × R` submatrix of `unfolding W d j` along `(ρ, σ)` is a
unit in `Matrix (Fin R) (Fin R) ℚ` (equivalently has nonzero determinant
over `ℚ`).

Informally, one chooses ρ as the natural inclusion (the first `R` rows)
and σ to pick out columns each of which contains a prime in row 0 only —
guaranteed by the density of primes in residue classes coprime to `W`
modulo `W^(d-j)`. The IsUnit witness then comes from the resulting
"diagonal-dominated" `R × R` exhibit. The full prime-counting argument is
the *only* outstanding piece in the formal proof of E2.1; once
formalised, both `lower_bound` and `mps_bond_dim` close immediately.

Two proof routes are open for a future session:

* (A) Use Bertrand-type prime existence in
  `[i·W^(d-j) + 1, (i+1)·W^(d-j)]` for every `0 ≤ i < R`, plus a
  dovetailing of residue classes mod `W^(d-j)`.
* (B) Replace the prime-density appeal by a generic Vandermonde-style
  exhibit over a finite extension of ℚ, sidestepping arithmetic
  entirely. This is the lighter-weight path discussed in
  `mps_bond_dim_notes.md`.
-/
theorem exists_invertible_submatrix
    (W d j : ℕ) (hW : 2 ≤ W) (hj_lo : 1 ≤ j) (hj_hi : j < d) :
    ∃ (ρ : Fin (min (W ^ j) (Nat.totient W * W ^ (d - j - 1) + 1)) → Fin (W ^ j))
      (σ : Fin (min (W ^ j) (Nat.totient W * W ^ (d - j - 1) + 1)) → Fin (W ^ (d - j))),
      IsUnit ((unfolding W d j).submatrix ρ σ) := by
  sorry

/--
Lower bound. Reduction-only: given the prime exhibit
(`exists_invertible_submatrix`), the lower bound is closed mechanically
by `Matrix.rank_of_isUnit` (an `R × R` unit matrix has rank `R`) followed
by `Matrix.rank_submatrix_le` (rank only decreases when restricting to
a submatrix). This proof contains no `sorry` of its own; its outstanding
content lives entirely inside `exists_invertible_submatrix`.
-/
theorem lower_bound
    (W d j : ℕ) (hW : 2 ≤ W) (hj_lo : 1 ≤ j) (hj_hi : j < d) :
    min (W ^ j) (Nat.totient W * W ^ (d - j - 1) + 1) ≤
      (unfolding W d j).rank := by
  classical
  -- Pull the prime exhibit.
  obtain ⟨ρ, σ, hUnit⟩ := exists_invertible_submatrix W d j hW hj_lo hj_hi
  -- An `R × R` unit matrix over `ℚ` has rank `Fintype.card (Fin R) = R`.
  have h_eq :
      ((unfolding W d j).submatrix ρ σ).rank =
        min (W ^ j) (Nat.totient W * W ^ (d - j - 1) + 1) := by
    have h := Matrix.rank_of_isUnit ((unfolding W d j).submatrix ρ σ) hUnit
    rw [Fintype.card_fin] at h
    exact h
  -- Restricting to a submatrix can only decrease the rank.
  calc min (W ^ j) (Nat.totient W * W ^ (d - j - 1) + 1)
      = ((unfolding W d j).submatrix ρ σ).rank := h_eq.symm
    _ ≤ (unfolding W d j).rank := Matrix.rank_submatrix_le _ _ _

/--
Trivial rank ceiling: any `m × n` matrix over a field has rank at most
`min m n`. Stated as a lemma to be cited in `upper_bound` so the
`min (W^j) _` form of the answer drops out cleanly.

This is the only one of the four stage-1 placeholders that is fully
proved (it's a direct citation of `Matrix.rank_le_height`).
-/
theorem rank_le_min_dim
    (W d j : ℕ) :
    (unfolding W d j).rank ≤ W ^ j :=
  Matrix.rank_le_height (unfolding W d j)

/--
**Main theorem (E2.1).** For every `W ≥ 2` and every cut `1 ≤ j < d`,
the rank of the `j`-th unfolding of the base-`W` prime indicator equals
`min (W^j) (φ(W) * W^(d-j-1) + 1)`.

This is the closed form proved (informally) in
`novel/mps_bond_dimension.md` and saturated empirically for
`W ∈ {2, 6, 30, 210}` and `d` up to 20.

The proof is a direct `Nat.le_antisymm` of `upper_bound` (combined with
`rank_le_min_dim` via `Nat.le_min`) and `lower_bound`. Once `lower_bound`
loses its `sorry`, this theorem is automatically closed without
modification.
-/
theorem mps_bond_dim
    (W d j : ℕ) (hW : 2 ≤ W) (hj_lo : 1 ≤ j) (hj_hi : j < d) :
    (unfolding W d j).rank =
      min (W ^ j) (Nat.totient W * W ^ (d - j - 1) + 1) :=
  Nat.le_antisymm
    (Nat.le_min.mpr ⟨rank_le_min_dim W d j, upper_bound W d j hW hj_lo hj_hi⟩)
    (lower_bound W d j hW hj_lo hj_hi)

end E2_1
