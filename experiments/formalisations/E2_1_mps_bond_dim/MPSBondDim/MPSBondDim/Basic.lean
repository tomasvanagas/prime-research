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

The cheap part of the lower bound — `1 ≤ rank` — is closed
unconditionally using just `Nat.prime_two` and the row-0 entry at column
`k = 1` (which is `chiP 2 = 1`). The hard part is establishing
`R ≤ rank` for `R ≥ 2`, which requires Hoheisel-type prime-density
results (primes in shrinking-density intervals) that mathlib does not
currently have.
-/

/--
Trivial fact: `chiP 2 = 1` since `2` is prime. This is the only piece
of the lower bound that is unconditional — no prime-density content.
-/
theorem chiP_two_eq_one : chiP 2 = 1 := by
  simp [chiP, Nat.prime_two]

/--
The matrix entry `unfolding W d j ⟨0, _⟩ ⟨1, _⟩ = 1`. Concretely, the
`(0, 1)` position of the unfolding equals `chiP(0 · W^(d-j) + 1 + 1) =
chiP 2 = 1`. The hypotheses `h0` and `h1` ensure the indices are valid;
the values `W, d, j` themselves are otherwise unconstrained at this
intermediate step.

This entry will witness `1 ≤ rank` in `one_le_rank_unfolding` below.
-/
theorem entry_zero_one_eq_one
    (W d j : ℕ) (h0 : 0 < W ^ j) (h1 : 1 < W ^ (d - j)) :
    unfolding W d j ⟨0, h0⟩ ⟨1, h1⟩ = 1 := by
  change chiP (0 * W ^ (d - j) + 1 + 1) = 1
  simp [chiP_two_eq_one]

/--
**Trivial lower bound:** `1 ≤ (unfolding W d j).rank`. This uses only
`Nat.prime_two` (no prime-density beyond Bertrand). The witness is the
`1 × 1` submatrix at row `0`, column `1`, whose entry is
`chiP 2 = 1` ≠ 0, hence `IsUnit` over `ℚ`.

This closes the `R = 1` portion of the lower bound. The remaining
content of `exists_invertible_submatrix` (= `R ≤ rank` for `R ≥ 2`)
requires primes in shrinking-density intervals (Hoheisel-type),
which mathlib does not yet have.
-/
theorem one_le_rank_unfolding
    (W d j : ℕ) (hW : 2 ≤ W) (hj_lo : 1 ≤ j) (hj_hi : j < d) :
    1 ≤ (unfolding W d j).rank := by
  have hWpos : 0 < W := by omega
  have hWj_pos : 0 < W ^ j := Nat.pow_pos hWpos
  have hdj : 1 ≤ d - j := Nat.sub_pos_of_lt hj_hi
  -- `1 < W^(d-j)`: since `W ≥ 2` and `d - j ≥ 1`, `W^(d-j) ≥ W ≥ 2`.
  have hWdj_gt_one : 1 < W ^ (d - j) := by
    have : W ≤ W ^ (d - j) := by
      have := Nat.pow_le_pow_right hWpos hdj
      simpa using this
    omega
  -- Build the 1×1 submatrix: row index 0, column index 1.
  let ρ : Fin 1 → Fin (W ^ j) := fun _ => ⟨0, hWj_pos⟩
  let σ : Fin 1 → Fin (W ^ (d - j)) := fun _ => ⟨1, hWdj_gt_one⟩
  -- Its single entry is `chiP 2 = 1`, so it is a unit in `Matrix (Fin 1) (Fin 1) ℚ`.
  have hUnit : IsUnit ((unfolding W d j).submatrix ρ σ) := by
    rw [Matrix.isUnit_iff_isUnit_det]
    rw [Matrix.det_fin_one]
    change IsUnit (unfolding W d j ⟨0, hWj_pos⟩ ⟨1, hWdj_gt_one⟩)
    rw [entry_zero_one_eq_one W d j hWj_pos hWdj_gt_one]
    exact isUnit_one
  -- `rank` of a 1×1 unit matrix is 1; restricting can only decrease rank.
  have h_eq : ((unfolding W d j).submatrix ρ σ).rank = 1 := by
    have h := Matrix.rank_of_isUnit ((unfolding W d j).submatrix ρ σ) hUnit
    simpa using h
  calc 1 = ((unfolding W d j).submatrix ρ σ).rank := h_eq.symm
    _ ≤ (unfolding W d j).rank := Matrix.rank_submatrix_le _ _ _

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

/-!
### Corner case: `W = 2`, `j = 1` — closed unconditionally via Bertrand

When `W = 2` and `j = 1`, the formula gives
`R = min(2, φ(2) · 2^(d-2) + 1) = min(2, 2^(d-2) + 1) = 2` for every `d ≥ 2`.
This corner case can be closed without any Hoheisel-grade prime-density
result: Bertrand's postulate alone exhibits a `2 × 2` invertible
submatrix, hence `lower_bound` (and therefore `mps_bond_dim`) holds
unconditionally for `(W, j) = (2, 1)`.

This is Route A' from `mps_bond_dim_notes.md`: a narrow but real closure
using only mathlib's `Nat.exists_prime_lt_and_le_two_mul`.
-/

/--
**Corner-case prime exhibit (W = 2, j = 1).**

Specialisation of `exists_invertible_submatrix` to the corner case
`W = 2`, `j = 1`, where `R = min(2, 2^(d-2) + 1) = 2`. Closed
unconditionally using only Bertrand's postulate
(`Nat.exists_prime_lt_and_le_two_mul`).

Construction: pick a prime `p ∈ (2^(d-1), 2 · 2^(d-1)]` (Bertrand at
`n = 2^(d-1)`). Set
* `ρ = id` on `Fin 2`,
* `σ 0 = ⟨1, _⟩` (column with `chiP 2 = 1` in row 0),
* `σ 1 = ⟨p − 2^(d-1) − 1, _⟩` (column with `chiP p = 1` in row 1).

The 2×2 submatrix is upper-triangular:
```
   ⎡ chiP 2,                  ?                  ⎤   ⎡ 1, ? ⎤
   ⎣ chiP (2^(d-1) + 2),      chiP p             ⎦ = ⎣ 0, 1 ⎦
```
since `2^(d-1) + 2` is even and `≥ 4`, hence not prime. The determinant
is `1`, so the submatrix is a unit over `ℚ`.
-/
theorem exists_invertible_submatrix_W_eq_2_j_eq_1
    (d : ℕ) (hd : 2 ≤ d) :
    ∃ (ρ : Fin 2 → Fin (2 ^ 1))
      (σ : Fin 2 → Fin (2 ^ (d - 1))),
      IsUnit ((unfolding 2 d 1).submatrix ρ σ) := by
  -- Setup.
  have h_d_minus_one : 1 ≤ d - 1 := by omega
  have h_pow_pos : 0 < 2 ^ (d - 1) := Nat.pow_pos (by norm_num)
  have h_pow_ge_2 : 2 ≤ 2 ^ (d - 1) := by
    calc (2 : ℕ) = 2 ^ 1 := by norm_num
      _ ≤ 2 ^ (d - 1) := Nat.pow_le_pow_right (by norm_num) h_d_minus_one
  -- Bertrand: prime `p ∈ (2^(d-1), 2 · 2^(d-1)]`.
  obtain ⟨p, hp_prime, hp_lo, hp_hi⟩ :=
    Nat.exists_prime_lt_and_le_two_mul (2 ^ (d - 1)) h_pow_pos.ne'
  -- Index bounds.
  have h_one_lt_pow : (1 : ℕ) < 2 ^ (d - 1) := by omega
  have h_p_minus_lt : p - 2 ^ (d - 1) - 1 < 2 ^ (d - 1) := by omega
  have h_zero_lt_pow_one : (0 : ℕ) < 2 ^ 1 := by norm_num
  have h_one_lt_pow_one : (1 : ℕ) < 2 ^ 1 := by norm_num
  -- Key non-primality fact: `2^(d-1) + 2` is even and `> 2`.
  have h_not_prime : ¬ Nat.Prime (2 ^ (d - 1) + 2) := by
    intro hp
    have h_dvd : 2 ∣ 2 ^ (d - 1) + 2 := by
      have h1 : 2 ∣ 2 ^ (d - 1) :=
        dvd_pow_self 2 (Nat.one_le_iff_ne_zero.mp h_d_minus_one)
      exact dvd_add h1 (dvd_refl 2)
    rcases hp.eq_one_or_self_of_dvd 2 h_dvd with h | h
    · omega
    · omega
  -- Build `ρ` and `σ`. We use ordinary `if-then-else` rather than `![_, _]`
  -- to keep beta-reduction transparent for the four entry computations.
  refine ⟨fun i => if i.val = 0 then ⟨0, h_zero_lt_pow_one⟩ else ⟨1, h_one_lt_pow_one⟩,
          fun i => if i.val = 0 then ⟨1, h_one_lt_pow⟩ else
                                        ⟨p - 2 ^ (d - 1) - 1, h_p_minus_lt⟩, ?_⟩
  -- Reduce `IsUnit` of the matrix to `IsUnit` of its determinant, then compute.
  rw [Matrix.isUnit_iff_isUnit_det, Matrix.det_fin_two]
  -- Compute the four entries.
  have hM00 : unfolding 2 d 1 ⟨0, h_zero_lt_pow_one⟩ ⟨1, h_one_lt_pow⟩ = 1 := by
    change chiP (0 * 2 ^ (d - 1) + 1 + 1) = 1
    simp [chiP, Nat.prime_two]
  have hM10 :
      unfolding 2 d 1 ⟨1, h_one_lt_pow_one⟩ ⟨1, h_one_lt_pow⟩ = 0 := by
    change chiP (1 * 2 ^ (d - 1) + 1 + 1) = 0
    have h_eq : 1 * 2 ^ (d - 1) + 1 + 1 = 2 ^ (d - 1) + 2 := by ring
    rw [h_eq]
    simp [chiP, h_not_prime]
  have hM11 :
      unfolding 2 d 1 ⟨1, h_one_lt_pow_one⟩
        ⟨p - 2 ^ (d - 1) - 1, h_p_minus_lt⟩ = 1 := by
    change chiP (1 * 2 ^ (d - 1) + (p - 2 ^ (d - 1) - 1) + 1) = 1
    have h_eq : 1 * 2 ^ (d - 1) + (p - 2 ^ (d - 1) - 1) + 1 = p := by omega
    rw [h_eq]
    simp [chiP, hp_prime]
  -- The (0, 1) entry is irrelevant; the determinant is `1·1 - x·0 = 1`.
  -- Reduce the four `submatrix` accesses; the `if-then-else` with `Fin` literals
  -- collapses, leaving the four explicit `unfolding` entries.
  have h0 : (0 : Fin 2).val = 0 := rfl
  have h1 : (1 : Fin 2).val = 1 := rfl
  simp only [Matrix.submatrix_apply, h0, h1, if_true, Nat.one_ne_zero, if_false]
  rw [hM00, hM10, hM11]
  -- Goal: `IsUnit (1 * 1 - x * 0)` where `x` is the row-0 column-1 entry.
  ring_nf
  exact isUnit_one

/--
**Corner-case main theorem (W = 2, j = 1).** The unfolding rank is exactly
`2` for every `d ≥ 2`. This is the only case of `mps_bond_dim` that is
currently formalisable in mathlib: it follows from `upper_bound` +
`rank_le_min_dim` for the `≤` direction and from
`exists_invertible_submatrix_W_eq_2_j_eq_1` for the `≥` direction.

The general `mps_bond_dim` requires `exists_invertible_submatrix` whose
proof is the only remaining `sorry` in this file (Hoheisel-grade
prime-density content).
-/
theorem mps_bond_dim_W_eq_2_j_eq_1
    (d : ℕ) (hd : 2 ≤ d) :
    (unfolding 2 d 1).rank = 2 := by
  apply Nat.le_antisymm
  · -- `≤ 2`: from `rank_le_min_dim` (i.e. `rank ≤ W^j = 2^1 = 2`).
    have h := rank_le_min_dim 2 d 1
    simpa using h
  · -- `2 ≤ rank`: from the corner-case prime exhibit.
    obtain ⟨ρ, σ, hUnit⟩ := exists_invertible_submatrix_W_eq_2_j_eq_1 d hd
    have h_eq : ((unfolding 2 d 1).submatrix ρ σ).rank = 2 := by
      have h := Matrix.rank_of_isUnit ((unfolding 2 d 1).submatrix ρ σ) hUnit
      simpa using h
    calc (2 : ℕ) = ((unfolding 2 d 1).submatrix ρ σ).rank := h_eq.symm
      _ ≤ (unfolding 2 d 1).rank := Matrix.rank_submatrix_le _ _ _

/-!
### Corner case: `W = 2`, `d = j + 1` — closed unconditionally without Bertrand

When `W = 2` and `d = j + 1` (so `d - j = 1`), the formula gives
`R = min(2^j, φ(2) · 2^0 + 1) = min(2^j, 2) = 2` for every `j ≥ 1`.
**Even simpler than the `(W = 2, j = 1)` case**: the matrix has only two
columns, so we take both, and rows `{0, 1}` give the `2 × 2` submatrix
```
   ⎡ chiP 1, chiP 2 ⎤   ⎡ 0, 1 ⎤
   ⎣ chiP 3, chiP 4 ⎦ = ⎣ 1, 0 ⎦
```
of determinant `−1`. Only `Nat.prime_two` and `Nat.prime_three` are used
— no Bertrand, no prime-density beyond the explicit small primes.

This case overlaps with the previous corner at `(j, d) = (1, 2)` (which
is also `j = 1`, `d = j + 1 = 2`); the new content is `j ≥ 2`.
-/

/--
`chiP 3 = 1` since `3` is prime.
-/
theorem chiP_three_eq_one : chiP 3 = 1 := by
  have h_prime_3 : Nat.Prime 3 := by decide
  simp [chiP, h_prime_3]

/--
**Corner-case prime exhibit (W = 2, d = j + 1).**

Specialisation of `exists_invertible_submatrix` to the corner case
`W = 2`, `d = j + 1`, where `R = min(2^j, 2) = 2` for every `j ≥ 1`.
Closed unconditionally using only `Nat.prime_two` and `Nat.prime_three`
(no Bertrand, no Hoheisel).

Construction: with the column-swap σ
* `ρ = id` on `Fin 2`,
* `σ 0 = ⟨1, _⟩`, `σ 1 = ⟨0, _⟩`,

the `2 × 2` submatrix is
```
   ⎡ unfolding (0, 1),  unfolding (0, 0) ⎤   ⎡ chiP 2, chiP 1 ⎤   ⎡ 1, 0 ⎤
   ⎣ unfolding (1, 1),  unfolding (1, 0) ⎦ = ⎣ chiP 4, chiP 3 ⎦ = ⎣ 0, 1 ⎦
```
of determinant `1`, hence `IsUnit` over `ℚ`.
-/
theorem exists_invertible_submatrix_W_eq_2_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    ∃ (ρ : Fin 2 → Fin (2 ^ j))
      (σ : Fin 2 → Fin (2 ^ ((j + 1) - j))),
      IsUnit ((unfolding 2 (j + 1) j).submatrix ρ σ) := by
  -- The exponent simplifies: `(j + 1) - j = 1`, hence `2 ^ ((j + 1) - j) = 2`.
  have h_sub : (j + 1 : ℕ) - j = 1 := by omega
  have h_pow_j_pos : 0 < 2 ^ j := Nat.pow_pos (by norm_num)
  have h_one_lt_pow_j : 1 < 2 ^ j := by
    calc (1 : ℕ) < 2 := by norm_num
      _ = 2 ^ 1 := by norm_num
      _ ≤ 2 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_pow_dj_pos : 0 < 2 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_one_lt_pow_dj : 1 < 2 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  -- Build `ρ` as identity-on-`Fin 2`, `σ` as the column swap (0 ↔ 1).
  refine ⟨fun i => if i.val = 0 then ⟨0, h_pow_j_pos⟩ else ⟨1, h_one_lt_pow_j⟩,
          fun i => if i.val = 0 then ⟨1, h_one_lt_pow_dj⟩ else ⟨0, h_pow_dj_pos⟩, ?_⟩
  -- Reduce `IsUnit` of the `2 × 2` submatrix to `IsUnit` of its determinant.
  rw [Matrix.isUnit_iff_isUnit_det, Matrix.det_fin_two]
  -- Compute the four entries.
  -- Entry (0, 0) of submatrix = unfolding(0, 1) = chiP(0·2^((j+1)-j) + 1 + 1) = chiP 2.
  have hM00 : unfolding 2 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨1, h_one_lt_pow_dj⟩ = 1 := by
    change chiP (0 * 2 ^ ((j + 1) - j) + 1 + 1) = 1
    simp [chiP, Nat.prime_two]
  -- Entry (0, 1) of submatrix = unfolding(0, 0) = chiP(0·… + 0 + 1) = chiP 1 = 0.
  have hM01 : unfolding 2 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨0, h_pow_dj_pos⟩ = 0 := by
    change chiP (0 * 2 ^ ((j + 1) - j) + 0 + 1) = 0
    simp [chiP, Nat.not_prime_one]
  -- Entry (1, 0) of submatrix = unfolding(1, 1) = chiP(1·2 + 1 + 1) = chiP 4 = 0.
  have hM10 :
      unfolding 2 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (1 * 2 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_not_prime_4 : ¬ Nat.Prime 4 := by decide
    simp [chiP, h_not_prime_4]
  -- Entry (1, 1) of submatrix = unfolding(1, 0) = chiP(1·2 + 0 + 1) = chiP 3 = 1.
  have hM11 :
      unfolding 2 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨0, h_pow_dj_pos⟩ = 1 := by
    change chiP (1 * 2 ^ ((j + 1) - j) + 0 + 1) = 1
    rw [h_sub]
    change chiP (1 * 2 + 0 + 1) = 1
    have h_eq3 : (1 * 2 + 0 + 1 : ℕ) = 3 := by norm_num
    rw [h_eq3]
    exact chiP_three_eq_one
  -- Reduce the `if-then-else` ρ, σ at the four `(0, 0), (0, 1), (1, 0), (1, 1)` cells.
  have h0 : (0 : Fin 2).val = 0 := rfl
  have h1 : (1 : Fin 2).val = 1 := rfl
  simp only [Matrix.submatrix_apply, h0, h1, if_true, Nat.one_ne_zero, if_false]
  rw [hM00, hM01, hM10, hM11]
  -- Goal: `IsUnit (1 * 1 - 0 * 0)` over `ℚ`.
  ring_nf
  exact isUnit_one

/--
**Corner-case main theorem (W = 2, d = j + 1).** The unfolding rank is
exactly `2` for every `j ≥ 1`. Mirrors `mps_bond_dim_W_eq_2_j_eq_1`
and is fully formalised in Lean (no `sorry`, no `axiom`).

This case overlaps the `(W = 2, j = 1)` corner at `j = 1, d = 2`; the
genuinely new content is the family `j ≥ 2`. Together with the previous
corner, we now have unconditional Lean proofs of `mps_bond_dim` whenever
either `j = 1` or `d - j = 1` — i.e. on the entire boundary of the
`(j, d - j)` grid.
-/
theorem mps_bond_dim_W_eq_2_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    (unfolding 2 (j + 1) j).rank = 2 := by
  apply Nat.le_antisymm
  · -- `≤ 2`: `unfolding 2 (j+1) j` has only `2^((j+1)-j) = 2^1 = 2` columns,
    -- so its rank is at most `2` by `Matrix.rank_le_width`.
    have h := (unfolding 2 (j + 1) j).rank_le_width
    have h_eq : 2 ^ ((j + 1 : ℕ) - j) = 2 := by
      have h_sub : (j + 1 : ℕ) - j = 1 := by omega
      rw [h_sub]; norm_num
    linarith
  · -- `2 ≤ rank`: from the corner-case prime exhibit.
    obtain ⟨ρ, σ, hUnit⟩ := exists_invertible_submatrix_W_eq_2_d_eq_j_plus_1 j hj
    have h_eq : ((unfolding 2 (j + 1) j).submatrix ρ σ).rank = 2 := by
      have h := Matrix.rank_of_isUnit ((unfolding 2 (j + 1) j).submatrix ρ σ) hUnit
      simpa using h
    calc (2 : ℕ) = ((unfolding 2 (j + 1) j).submatrix ρ σ).rank := h_eq.symm
      _ ≤ (unfolding 2 (j + 1) j).rank := Matrix.rank_submatrix_le _ _ _

/-!
### Corner case: `W = 3`, `d = j + 1` — closed unconditionally without Bertrand

When `W = 3` and `d = j + 1` (so `d - j = 1`), the formula gives
`R = min(3^j, φ(3) · 3^0 + 1) = min(3^j, 3) = 3` for every `j ≥ 1`.
Like the `W = 2` orthogonal corner, the matrix has only `3^(d-j) = 3`
columns, so we take all of them; rows `{0, 1, 2}` (available since
`3^j ≥ 3` for `j ≥ 1`) give the `3 × 3` submatrix
```
   ⎡ chiP 1, chiP 2, chiP 3 ⎤   ⎡ 0, 1, 1 ⎤
   ⎢ chiP 4, chiP 5, chiP 6 ⎥ = ⎢ 0, 1, 0 ⎥
   ⎣ chiP 7, chiP 8, chiP 9 ⎦   ⎣ 1, 0, 0 ⎦
```
of determinant `−1`, hence `IsUnit` over `ℚ`. Only the explicit primes
`2, 3, 5, 7` and the non-primality of `1, 4, 6, 8, 9` are required —
no Bertrand, no Hoheisel-grade prime density.

Together with the `W = 2` corners, this extends the unconditional
formalisation to a second base wheel.
-/

/--
`chiP 5 = 1` since `5` is prime.
-/
theorem chiP_five_eq_one : chiP 5 = 1 := by
  have h_prime_5 : Nat.Prime 5 := by decide
  simp [chiP, h_prime_5]

/--
`chiP 7 = 1` since `7` is prime.
-/
theorem chiP_seven_eq_one : chiP 7 = 1 := by
  have h_prime_7 : Nat.Prime 7 := by decide
  simp [chiP, h_prime_7]

/--
**Corner-case prime exhibit (W = 3, d = j + 1).**

Specialisation of `exists_invertible_submatrix` to the corner case
`W = 3`, `d = j + 1`, where `R = min(3^j, 3) = 3` for every `j ≥ 1`.
Closed unconditionally using only `Nat.prime_two`, `Nat.prime_three`,
`chiP_five_eq_one`, `chiP_seven_eq_one`, and decidability of the
non-primality of `1, 4, 6, 8, 9`.

Construction:
* `ρ = id` on `Fin 3` (rows `0, 1, 2`),
* `σ = id` on `Fin 3` (columns `0, 1, 2`),

giving the `3 × 3` submatrix
```
   ⎡ chiP 1, chiP 2, chiP 3 ⎤   ⎡ 0, 1, 1 ⎤
   ⎢ chiP 4, chiP 5, chiP 6 ⎥ = ⎢ 0, 1, 0 ⎥
   ⎣ chiP 7, chiP 8, chiP 9 ⎦   ⎣ 1, 0, 0 ⎦
```
of determinant `−1`, hence `IsUnit` over `ℚ`.
-/
theorem exists_invertible_submatrix_W_eq_3_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    ∃ (ρ : Fin 3 → Fin (3 ^ j))
      (σ : Fin 3 → Fin (3 ^ ((j + 1) - j))),
      IsUnit ((unfolding 3 (j + 1) j).submatrix ρ σ) := by
  -- The exponent simplifies: `(j + 1) - j = 1`, hence `3 ^ ((j + 1) - j) = 3`.
  have h_sub : (j + 1 : ℕ) - j = 1 := by omega
  have h_pow_j_pos : 0 < 3 ^ j := Nat.pow_pos (by norm_num)
  have h_one_lt_pow_j : 1 < 3 ^ j := by
    calc (1 : ℕ) < 3 := by norm_num
      _ = 3 ^ 1 := by norm_num
      _ ≤ 3 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_two_lt_pow_j : 2 < 3 ^ j := by
    calc (2 : ℕ) < 3 := by norm_num
      _ = 3 ^ 1 := by norm_num
      _ ≤ 3 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_zero_lt_pow_dj : (0 : ℕ) < 3 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_one_lt_pow_dj : (1 : ℕ) < 3 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_two_lt_pow_dj : (2 : ℕ) < 3 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  -- Build `ρ` and `σ` as identity-on-`Fin 3` mappings.
  refine ⟨fun i => if i.val = 0 then ⟨0, h_pow_j_pos⟩
                   else if i.val = 1 then ⟨1, h_one_lt_pow_j⟩
                   else ⟨2, h_two_lt_pow_j⟩,
          fun i => if i.val = 0 then ⟨0, h_zero_lt_pow_dj⟩
                   else if i.val = 1 then ⟨1, h_one_lt_pow_dj⟩
                   else ⟨2, h_two_lt_pow_dj⟩, ?_⟩
  -- Reduce `IsUnit` of the `3 × 3` submatrix to `IsUnit` of its determinant.
  rw [Matrix.isUnit_iff_isUnit_det, Matrix.det_fin_three]
  -- Non-primality lemmas for the zero entries.
  have h_not_prime_4 : ¬ Nat.Prime 4 := by decide
  have h_not_prime_6 : ¬ Nat.Prime 6 := by decide
  have h_not_prime_8 : ¬ Nat.Prime 8 := by decide
  have h_not_prime_9 : ¬ Nat.Prime 9 := by decide
  -- Compute the nine entries of the submatrix.
  -- Entry (0, 0): chiP(0 · 3^((j+1)-j) + 0 + 1) = chiP 1 = 0.
  have hM00 : unfolding 3 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (0 * 3 ^ ((j + 1) - j) + 0 + 1) = 0
    simp [chiP, Nat.not_prime_one]
  -- Entry (0, 1): chiP(0 · 3^((j+1)-j) + 1 + 1) = chiP 2 = 1.
  have hM01 : unfolding 3 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨1, h_one_lt_pow_dj⟩ = 1 := by
    change chiP (0 * 3 ^ ((j + 1) - j) + 1 + 1) = 1
    simp [chiP, Nat.prime_two]
  -- Entry (0, 2): chiP(0 · 3^((j+1)-j) + 2 + 1) = chiP 3 = 1.
  have hM02 : unfolding 3 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨2, h_two_lt_pow_dj⟩ = 1 := by
    change chiP (0 * 3 ^ ((j + 1) - j) + 2 + 1) = 1
    have h_eq : (0 * 3 ^ ((j + 1) - j) + 2 + 1 : ℕ) = 3 := by simp
    rw [h_eq]
    exact chiP_three_eq_one
  -- Entry (1, 0): chiP(1 · 3^((j+1)-j) + 0 + 1) = chiP 4 = 0.
  have hM10 : unfolding 3 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (1 * 3 ^ ((j + 1) - j) + 0 + 1) = 0
    rw [h_sub]
    have h_eq : (1 * 3 ^ 1 + 0 + 1 : ℕ) = 4 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_4]
  -- Entry (1, 1): chiP(1 · 3^((j+1)-j) + 1 + 1) = chiP 5 = 1.
  have hM11 : unfolding 3 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 1 := by
    change chiP (1 * 3 ^ ((j + 1) - j) + 1 + 1) = 1
    rw [h_sub]
    have h_eq : (1 * 3 ^ 1 + 1 + 1 : ℕ) = 5 := by norm_num
    rw [h_eq]
    exact chiP_five_eq_one
  -- Entry (1, 2): chiP(1 · 3^((j+1)-j) + 2 + 1) = chiP 6 = 0.
  have hM12 : unfolding 3 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨2, h_two_lt_pow_dj⟩ = 0 := by
    change chiP (1 * 3 ^ ((j + 1) - j) + 2 + 1) = 0
    rw [h_sub]
    have h_eq : (1 * 3 ^ 1 + 2 + 1 : ℕ) = 6 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_6]
  -- Entry (2, 0): chiP(2 · 3^((j+1)-j) + 0 + 1) = chiP 7 = 1.
  have hM20 : unfolding 3 (j + 1) j ⟨2, h_two_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 1 := by
    change chiP (2 * 3 ^ ((j + 1) - j) + 0 + 1) = 1
    rw [h_sub]
    have h_eq : (2 * 3 ^ 1 + 0 + 1 : ℕ) = 7 := by norm_num
    rw [h_eq]
    exact chiP_seven_eq_one
  -- Entry (2, 1): chiP(2 · 3^((j+1)-j) + 1 + 1) = chiP 8 = 0.
  have hM21 : unfolding 3 (j + 1) j ⟨2, h_two_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (2 * 3 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (2 * 3 ^ 1 + 1 + 1 : ℕ) = 8 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_8]
  -- Entry (2, 2): chiP(2 · 3^((j+1)-j) + 2 + 1) = chiP 9 = 0.
  have hM22 : unfolding 3 (j + 1) j ⟨2, h_two_lt_pow_j⟩ ⟨2, h_two_lt_pow_dj⟩ = 0 := by
    change chiP (2 * 3 ^ ((j + 1) - j) + 2 + 1) = 0
    rw [h_sub]
    have h_eq : (2 * 3 ^ 1 + 2 + 1 : ℕ) = 9 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_9]
  -- Reduce the `if-then-else` ρ, σ at the nine cells via `Fin.val` of `0, 1, 2`.
  have h0 : (0 : Fin 3).val = 0 := rfl
  have h1 : (1 : Fin 3).val = 1 := rfl
  have h2 : (2 : Fin 3).val = 2 := rfl
  have hne_2_0 : (2 : ℕ) ≠ 0 := by decide
  have hne_2_1 : (2 : ℕ) ≠ 1 := by decide
  simp only [Matrix.submatrix_apply, h0, h1, h2, if_true, if_false,
             Nat.one_ne_zero, hne_2_0, hne_2_1]
  rw [hM00, hM01, hM02, hM10, hM11, hM12, hM20, hM21, hM22]
  -- Goal: `IsUnit (0·1·0 - 0·0·0 - 1·0·0 + 1·0·1 + 1·0·0 - 1·1·1)` over `ℚ`.
  --     = `IsUnit (-1 : ℚ)`.
  ring_nf
  exact isUnit_one.neg

/--
**Corner-case main theorem (W = 3, d = j + 1).** The unfolding rank is
exactly `3` for every `j ≥ 1`. Mirrors `mps_bond_dim_W_eq_2_d_eq_j_plus_1`
and is fully formalised in Lean (no `sorry`, no `axiom`).

Together with the two `W = 2` corners, this is the third unconditional
instance of `mps_bond_dim` and the first instance over the wheel `W = 3`.
The general `mps_bond_dim` still requires `exists_invertible_submatrix`
whose proof is the only remaining `sorry` in this file.
-/
theorem mps_bond_dim_W_eq_3_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    (unfolding 3 (j + 1) j).rank = 3 := by
  apply Nat.le_antisymm
  · -- `≤ 3`: `unfolding 3 (j+1) j` has only `3^((j+1)-j) = 3^1 = 3` columns.
    have h := (unfolding 3 (j + 1) j).rank_le_width
    have h_eq : 3 ^ ((j + 1 : ℕ) - j) = 3 := by
      have h_sub : (j + 1 : ℕ) - j = 1 := by omega
      rw [h_sub]; norm_num
    linarith
  · -- `3 ≤ rank`: from the corner-case prime exhibit.
    obtain ⟨ρ, σ, hUnit⟩ := exists_invertible_submatrix_W_eq_3_d_eq_j_plus_1 j hj
    have h_eq : ((unfolding 3 (j + 1) j).submatrix ρ σ).rank = 3 := by
      have h := Matrix.rank_of_isUnit ((unfolding 3 (j + 1) j).submatrix ρ σ) hUnit
      simpa using h
    calc (3 : ℕ) = ((unfolding 3 (j + 1) j).submatrix ρ σ).rank := h_eq.symm
      _ ≤ (unfolding 3 (j + 1) j).rank := Matrix.rank_submatrix_le _ _ _

/-!
### Corner case: `W = 4`, `d = j + 1` — closed unconditionally without Bertrand

When `W = 4` and `d = j + 1` (so `d - j = 1`), the formula gives
`R = min(4^j, φ(4) · 4^0 + 1) = min(4^j, 3) = 3` for every `j ≥ 1`. The
matrix has `4^(d-j) = 4` columns; column `3` is `chiP` at multiples of
`4` (all zeros), so we pick the three live columns `{0, 1, 2}` (those
with `gcd(k+1, 4) = 1`, i.e. `k + 1 ∈ {1, 2, 3}`) and rows `{0, 1, 2}`
(available since `4^j ≥ 4` for `j ≥ 1`) to get the `3 × 3` submatrix
```
   ⎡ chiP 1, chiP 2, chiP 3  ⎤   ⎡ 0, 1, 1 ⎤
   ⎢ chiP 5, chiP 6, chiP 7  ⎥ = ⎢ 1, 0, 1 ⎥
   ⎣ chiP 9, chiP 10, chiP 11⎦   ⎣ 0, 0, 1 ⎦
```
of determinant `−1`, hence `IsUnit` over `ℚ`. Only the explicit primes
`2, 3, 5, 7, 11` and the non-primality of `1, 4, 6, 9, 10` are required
— no Bertrand, no Hoheisel-grade prime density.

**Upper-bound subtlety:** `rank_le_width` gives only `rank ≤ 4`, not the
sharp `rank ≤ 3`. We therefore cite the general `upper_bound` lemma —
which evaluates to `φ(4) · 4^0 + 1 = 2 + 1 = 3` for this corner — to
close the `≤` direction. This is the first orthogonal-corner instance
where `rank_le_width` is not tight; subsequent `W ∈ {6, 7, …}` corners
will follow the same pattern.
-/

/--
`chiP 11 = 1` since `11` is prime.
-/
theorem chiP_eleven_eq_one : chiP 11 = 1 := by
  have h_prime_11 : Nat.Prime 11 := by decide
  simp [chiP, h_prime_11]

/--
**Corner-case prime exhibit (W = 4, d = j + 1).**

Specialisation of `exists_invertible_submatrix` to the corner case
`W = 4`, `d = j + 1`, where `R = min(4^j, 3) = 3` for every `j ≥ 1`.
Closed unconditionally using only `Nat.prime_two`, `Nat.prime_three`,
`chiP_five_eq_one`, `chiP_seven_eq_one`, `chiP_eleven_eq_one`, and
decidability of the non-primality of `1, 4, 6, 9, 10`.

Construction:
* `ρ = id` on `Fin 3` (rows `0, 1, 2`),
* `σ : Fin 3 → Fin (4^((j+1)-j))` picks the three live columns
  `0, 1, 2`. The dropped column `3` corresponds to `chiP` at multiples
  of `4` and is identically zero (hence not needed for the rank-3
  exhibit).

The `3 × 3` submatrix is
```
   ⎡ chiP 1, chiP 2, chiP 3  ⎤   ⎡ 0, 1, 1 ⎤
   ⎢ chiP 5, chiP 6, chiP 7  ⎥ = ⎢ 1, 0, 1 ⎥
   ⎣ chiP 9, chiP 10, chiP 11⎦   ⎣ 0, 0, 1 ⎦
```
with determinant `-1`, hence `IsUnit` over `ℚ`.
-/
theorem exists_invertible_submatrix_W_eq_4_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    ∃ (ρ : Fin 3 → Fin (4 ^ j))
      (σ : Fin 3 → Fin (4 ^ ((j + 1) - j))),
      IsUnit ((unfolding 4 (j + 1) j).submatrix ρ σ) := by
  -- The exponent simplifies: `(j + 1) - j = 1`, hence `4 ^ ((j + 1) - j) = 4`.
  have h_sub : (j + 1 : ℕ) - j = 1 := by omega
  have h_pow_j_pos : 0 < 4 ^ j := Nat.pow_pos (by norm_num)
  have h_one_lt_pow_j : 1 < 4 ^ j := by
    calc (1 : ℕ) < 4 := by norm_num
      _ = 4 ^ 1 := by norm_num
      _ ≤ 4 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_two_lt_pow_j : 2 < 4 ^ j := by
    calc (2 : ℕ) < 4 := by norm_num
      _ = 4 ^ 1 := by norm_num
      _ ≤ 4 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_zero_lt_pow_dj : (0 : ℕ) < 4 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_one_lt_pow_dj : (1 : ℕ) < 4 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_two_lt_pow_dj : (2 : ℕ) < 4 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  -- Build `ρ` and `σ` as identity-on-`Fin 3` mappings.
  refine ⟨fun i => if i.val = 0 then ⟨0, h_pow_j_pos⟩
                   else if i.val = 1 then ⟨1, h_one_lt_pow_j⟩
                   else ⟨2, h_two_lt_pow_j⟩,
          fun i => if i.val = 0 then ⟨0, h_zero_lt_pow_dj⟩
                   else if i.val = 1 then ⟨1, h_one_lt_pow_dj⟩
                   else ⟨2, h_two_lt_pow_dj⟩, ?_⟩
  -- Reduce `IsUnit` of the `3 × 3` submatrix to `IsUnit` of its determinant.
  rw [Matrix.isUnit_iff_isUnit_det, Matrix.det_fin_three]
  -- Non-primality lemmas for the zero entries.
  have h_not_prime_4 : ¬ Nat.Prime 4 := by decide
  have h_not_prime_6 : ¬ Nat.Prime 6 := by decide
  have h_not_prime_9 : ¬ Nat.Prime 9 := by decide
  have h_not_prime_10 : ¬ Nat.Prime 10 := by decide
  -- Compute the nine entries of the submatrix.
  -- Entry (0, 0): chiP(0 · 4^((j+1)-j) + 0 + 1) = chiP 1 = 0.
  have hM00 : unfolding 4 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (0 * 4 ^ ((j + 1) - j) + 0 + 1) = 0
    simp [chiP, Nat.not_prime_one]
  -- Entry (0, 1): chiP(0 · 4^((j+1)-j) + 1 + 1) = chiP 2 = 1.
  have hM01 : unfolding 4 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨1, h_one_lt_pow_dj⟩ = 1 := by
    change chiP (0 * 4 ^ ((j + 1) - j) + 1 + 1) = 1
    simp [chiP, Nat.prime_two]
  -- Entry (0, 2): chiP(0 · 4^((j+1)-j) + 2 + 1) = chiP 3 = 1.
  have hM02 : unfolding 4 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨2, h_two_lt_pow_dj⟩ = 1 := by
    change chiP (0 * 4 ^ ((j + 1) - j) + 2 + 1) = 1
    have h_eq : (0 * 4 ^ ((j + 1) - j) + 2 + 1 : ℕ) = 3 := by simp
    rw [h_eq]
    exact chiP_three_eq_one
  -- Entry (1, 0): chiP(1 · 4^((j+1)-j) + 0 + 1) = chiP 5 = 1.
  have hM10 : unfolding 4 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 1 := by
    change chiP (1 * 4 ^ ((j + 1) - j) + 0 + 1) = 1
    rw [h_sub]
    have h_eq : (1 * 4 ^ 1 + 0 + 1 : ℕ) = 5 := by norm_num
    rw [h_eq]
    exact chiP_five_eq_one
  -- Entry (1, 1): chiP(1 · 4^((j+1)-j) + 1 + 1) = chiP 6 = 0.
  have hM11 : unfolding 4 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (1 * 4 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (1 * 4 ^ 1 + 1 + 1 : ℕ) = 6 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_6]
  -- Entry (1, 2): chiP(1 · 4^((j+1)-j) + 2 + 1) = chiP 7 = 1.
  have hM12 : unfolding 4 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨2, h_two_lt_pow_dj⟩ = 1 := by
    change chiP (1 * 4 ^ ((j + 1) - j) + 2 + 1) = 1
    rw [h_sub]
    have h_eq : (1 * 4 ^ 1 + 2 + 1 : ℕ) = 7 := by norm_num
    rw [h_eq]
    exact chiP_seven_eq_one
  -- Entry (2, 0): chiP(2 · 4^((j+1)-j) + 0 + 1) = chiP 9 = 0.
  have hM20 : unfolding 4 (j + 1) j ⟨2, h_two_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (2 * 4 ^ ((j + 1) - j) + 0 + 1) = 0
    rw [h_sub]
    have h_eq : (2 * 4 ^ 1 + 0 + 1 : ℕ) = 9 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_9]
  -- Entry (2, 1): chiP(2 · 4^((j+1)-j) + 1 + 1) = chiP 10 = 0.
  have hM21 : unfolding 4 (j + 1) j ⟨2, h_two_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (2 * 4 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (2 * 4 ^ 1 + 1 + 1 : ℕ) = 10 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_10]
  -- Entry (2, 2): chiP(2 · 4^((j+1)-j) + 2 + 1) = chiP 11 = 1.
  have hM22 : unfolding 4 (j + 1) j ⟨2, h_two_lt_pow_j⟩ ⟨2, h_two_lt_pow_dj⟩ = 1 := by
    change chiP (2 * 4 ^ ((j + 1) - j) + 2 + 1) = 1
    rw [h_sub]
    have h_eq : (2 * 4 ^ 1 + 2 + 1 : ℕ) = 11 := by norm_num
    rw [h_eq]
    exact chiP_eleven_eq_one
  -- Reduce the `if-then-else` ρ, σ at the nine cells via `Fin.val` of `0, 1, 2`.
  have h0 : (0 : Fin 3).val = 0 := rfl
  have h1 : (1 : Fin 3).val = 1 := rfl
  have h2 : (2 : Fin 3).val = 2 := rfl
  have hne_2_0 : (2 : ℕ) ≠ 0 := by decide
  have hne_2_1 : (2 : ℕ) ≠ 1 := by decide
  simp only [Matrix.submatrix_apply, h0, h1, h2, if_true, if_false,
             Nat.one_ne_zero, hne_2_0, hne_2_1]
  rw [hM00, hM01, hM02, hM10, hM11, hM12, hM20, hM21, hM22]
  -- Goal: `IsUnit (0·0·1 - 0·1·0 - 1·1·1 + 1·1·0 + 1·1·0 - 1·0·0)` over `ℚ`.
  --     = `IsUnit (-1 : ℚ)`.
  ring_nf
  exact isUnit_one.neg

/--
**Corner-case main theorem (W = 4, d = j + 1).** The unfolding rank is
exactly `3` for every `j ≥ 1`. Mirrors `mps_bond_dim_W_eq_3_d_eq_j_plus_1`
and is fully formalised in Lean (no `sorry`, no `axiom`).

**Upper-bound subtlety:** the matrix has `4^((j+1)-j) = 4` columns, so
`rank_le_width` only gives `rank ≤ 4`. We instead cite the general
`upper_bound` lemma, which evaluates to `φ(4) · 4^0 + 1 = 2 · 1 + 1 = 3`
in this corner — the first instance where the live-column count strictly
beats the column count.

Together with the `W ∈ {2, 3}` corners, this is the fourth unconditional
instance of `mps_bond_dim` and the second instance over a wheel `W ≥ 3`.
The general `mps_bond_dim` still requires `exists_invertible_submatrix`
whose proof is the only remaining `sorry` in this file.
-/
theorem mps_bond_dim_W_eq_4_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    (unfolding 4 (j + 1) j).rank = 3 := by
  apply Nat.le_antisymm
  · -- `≤ 3`: from `upper_bound`. (`rank_le_width` only gives `≤ 4`.)
    have hW : (2 : ℕ) ≤ 4 := by norm_num
    have hj_hi : j < j + 1 := Nat.lt_succ_self j
    have h := upper_bound 4 (j + 1) j hW hj hj_hi
    have h_eq : Nat.totient 4 * 4 ^ ((j + 1 : ℕ) - j - 1) + 1 = 3 := by
      have h_sub : (j + 1 : ℕ) - j - 1 = 0 := by omega
      rw [h_sub]
      decide
    linarith
  · -- `3 ≤ rank`: from the corner-case prime exhibit.
    obtain ⟨ρ, σ, hUnit⟩ := exists_invertible_submatrix_W_eq_4_d_eq_j_plus_1 j hj
    have h_eq : ((unfolding 4 (j + 1) j).submatrix ρ σ).rank = 3 := by
      have h := Matrix.rank_of_isUnit ((unfolding 4 (j + 1) j).submatrix ρ σ) hUnit
      simpa using h
    calc (3 : ℕ) = ((unfolding 4 (j + 1) j).submatrix ρ σ).rank := h_eq.symm
      _ ≤ (unfolding 4 (j + 1) j).rank := Matrix.rank_submatrix_le _ _ _

/-!
### Corner case: `W = 5`, `d = j + 1` — closed unconditionally without Bertrand

When `W = 5` and `d = j + 1` (so `d - j = 1`), the formula gives
`R = min(5^j, φ(5) · 5^0 + 1) = min(5^j, 5) = 5` for every `j ≥ 1`. The
matrix has `5^(d-j) = 5` columns; column `4` (orig col index `4`,
i.e. `chiP` at multiples of `5`) is dead beyond row `0` (only `chiP 5 = 1`
is nonzero in that column). We pick *all five columns* — the dead column
contributes the diagonal entry that bumps the rank from `4 = φ(5)` (live
columns alone) to `5`. Rows `{0, 1, 2, 3, 4}` are available since
`5^j ≥ 5` for `j ≥ 1`.

**Strategy:** mathlib has `Matrix.det_fin_three` but no `det_fin_four` or
`det_fin_five`. Direct evaluation of a `5 × 5` determinant is not
available. We therefore choose `ρ` and `σ` as *permutations* such that
the `5 × 5` submatrix becomes upper triangular with `1` on the diagonal:
```
   ⎡ chiP  5, chiP  4, chiP  1, chiP  2, chiP  3 ⎤   ⎡ 1, 0, 0, 1, 1 ⎤
   ⎢ chiP 20, chiP 19, chiP 16, chiP 17, chiP 18 ⎥   ⎢ 0, 1, 0, 1, 0 ⎥
   ⎢ chiP 15, chiP 14, chiP 11, chiP 12, chiP 13 ⎥ = ⎢ 0, 0, 1, 0, 1 ⎥.
   ⎢ chiP 10, chiP  9, chiP  6, chiP  7, chiP  8 ⎥   ⎢ 0, 0, 0, 1, 0 ⎥
   ⎣ chiP 25, chiP 24, chiP 21, chiP 22, chiP 23 ⎦   ⎣ 0, 0, 0, 0, 1 ⎦
```
Concretely `ρ` permutes original rows to `(0, 3, 2, 1, 4)` and `σ`
permutes original columns to `(4, 3, 0, 1, 2)`. The diagonal hits the
five primes `{5, 19, 11, 7, 23}` (all in `[1, 25]`) and the lower
triangle hits the ten composites `{20, 15, 14, 10, 9, 6, 25, 24, 21, 22}`
(all `≤ 25`). The determinant is then the product of the diagonal,
namely `1`, by `Matrix.det_of_upperTriangular`. Used helpers:
`Nat.prime_two`, `Nat.prime_three`, `chiP_five_eq_one`, `chiP_seven_eq_one`,
`chiP_eleven_eq_one`, `chiP_nineteen_eq_one` (new at S117),
`chiP_twenty_three_eq_one` (new at S117), and decidability of
non-primality for the ten composites listed above. No Bertrand, no
Hoheisel-grade prime density.

This is the **first orthogonal-corner instance with `R = W` (so all `W`
columns are needed)**, and the **first instance using
`Matrix.det_of_upperTriangular` rather than `Matrix.det_fin_three`** —
the proof technique scales to every `W` for which we can exhibit a
permutation of `Fin W` that triangularises the slab `chiP 1 .. chiP W^2`.
-/

/--
`chiP 19 = 1` since `19` is prime.
-/
theorem chiP_nineteen_eq_one : chiP 19 = 1 := by
  have h_prime_19 : Nat.Prime 19 := by decide
  simp [chiP, h_prime_19]

/--
`chiP 23 = 1` since `23` is prime.
-/
theorem chiP_twenty_three_eq_one : chiP 23 = 1 := by
  have h_prime_23 : Nat.Prime 23 := by decide
  simp [chiP, h_prime_23]

/--
**Corner-case prime exhibit (W = 5, d = j + 1).**

Specialisation of `exists_invertible_submatrix` to the corner case
`W = 5`, `d = j + 1`, where `R = min(5^j, 5) = 5` for every `j ≥ 1`.
Closed unconditionally using only `Nat.prime_two`, `Nat.prime_three`,
`chiP_five_eq_one`, `chiP_seven_eq_one`, `chiP_eleven_eq_one`,
`chiP_nineteen_eq_one`, `chiP_twenty_three_eq_one`, and decidability of
the non-primality of `1, 4, 6, 8, 9, 10, 14, 15, 20, 21, 22, 24, 25`.

The construction uses `ρ : Fin 5 → Fin (5^j)` mapping
`0 ↦ 0, 1 ↦ 3, 2 ↦ 2, 3 ↦ 1, 4 ↦ 4` and `σ : Fin 5 → Fin (5^((j+1)-j))`
mapping `0 ↦ 4, 1 ↦ 3, 2 ↦ 0, 3 ↦ 1, 4 ↦ 2`. The `5 × 5` submatrix is
upper triangular with `1` on the diagonal, hence `IsUnit` over `ℚ`.
-/
theorem exists_invertible_submatrix_W_eq_5_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    ∃ (ρ : Fin 5 → Fin (5 ^ j))
      (σ : Fin 5 → Fin (5 ^ ((j + 1) - j))),
      IsUnit ((unfolding 5 (j + 1) j).submatrix ρ σ) := by
  -- The exponent simplifies: `(j + 1) - j = 1`, hence `5 ^ ((j + 1) - j) = 5`.
  have h_sub : (j + 1 : ℕ) - j = 1 := by omega
  have h_pow_j_pos : 0 < 5 ^ j := Nat.pow_pos (by norm_num)
  have h_one_lt_pow_j : 1 < 5 ^ j := by
    calc (1 : ℕ) < 5 := by norm_num
      _ = 5 ^ 1 := by norm_num
      _ ≤ 5 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_two_lt_pow_j : 2 < 5 ^ j := by
    calc (2 : ℕ) < 5 := by norm_num
      _ = 5 ^ 1 := by norm_num
      _ ≤ 5 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_three_lt_pow_j : 3 < 5 ^ j := by
    calc (3 : ℕ) < 5 := by norm_num
      _ = 5 ^ 1 := by norm_num
      _ ≤ 5 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_four_lt_pow_j : 4 < 5 ^ j := by
    calc (4 : ℕ) < 5 := by norm_num
      _ = 5 ^ 1 := by norm_num
      _ ≤ 5 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_zero_lt_pow_dj : (0 : ℕ) < 5 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_one_lt_pow_dj : (1 : ℕ) < 5 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_two_lt_pow_dj : (2 : ℕ) < 5 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_three_lt_pow_dj : (3 : ℕ) < 5 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_four_lt_pow_dj : (4 : ℕ) < 5 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  -- ρ permutes original rows in the order `(0, 3, 2, 1, 4)`.
  -- σ permutes original columns in the order `(4, 3, 0, 1, 2)`.
  refine ⟨fun i => if i.val = 0 then ⟨0, h_pow_j_pos⟩
                   else if i.val = 1 then ⟨3, h_three_lt_pow_j⟩
                   else if i.val = 2 then ⟨2, h_two_lt_pow_j⟩
                   else if i.val = 3 then ⟨1, h_one_lt_pow_j⟩
                   else ⟨4, h_four_lt_pow_j⟩,
          fun i => if i.val = 0 then ⟨4, h_four_lt_pow_dj⟩
                   else if i.val = 1 then ⟨3, h_three_lt_pow_dj⟩
                   else if i.val = 2 then ⟨0, h_zero_lt_pow_dj⟩
                   else if i.val = 3 then ⟨1, h_one_lt_pow_dj⟩
                   else ⟨2, h_two_lt_pow_dj⟩, ?_⟩
  -- Reduce `IsUnit` to `IsUnit (det)` and compute the det as the
  -- product of the diagonal via `Matrix.det_of_upperTriangular`.
  rw [Matrix.isUnit_iff_isUnit_det]
  -- Non-primality lemmas for the ten below-diagonal zeros.
  have h_not_prime_6  : ¬ Nat.Prime 6  := by decide
  have h_not_prime_9  : ¬ Nat.Prime 9  := by decide
  have h_not_prime_10 : ¬ Nat.Prime 10 := by decide
  have h_not_prime_14 : ¬ Nat.Prime 14 := by decide
  have h_not_prime_15 : ¬ Nat.Prime 15 := by decide
  have h_not_prime_20 : ¬ Nat.Prime 20 := by decide
  have h_not_prime_21 : ¬ Nat.Prime 21 := by decide
  have h_not_prime_22 : ¬ Nat.Prime 22 := by decide
  have h_not_prime_24 : ¬ Nat.Prime 24 := by decide
  have h_not_prime_25 : ¬ Nat.Prime 25 := by decide
  -- Compute the 5 diagonal entries (all = 1).
  -- (0,0): unfolding(0, 4) = chiP 5 = 1
  have hD00 : unfolding 5 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨4, h_four_lt_pow_dj⟩ = 1 := by
    change chiP (0 * 5 ^ ((j + 1) - j) + 4 + 1) = 1
    rw [h_sub]
    have h_eq : (0 * 5 ^ 1 + 4 + 1 : ℕ) = 5 := by norm_num
    rw [h_eq]
    exact chiP_five_eq_one
  -- (1,1): unfolding(3, 3) = chiP 19 = 1
  have hD11 : unfolding 5 (j + 1) j ⟨3, h_three_lt_pow_j⟩ ⟨3, h_three_lt_pow_dj⟩ = 1 := by
    change chiP (3 * 5 ^ ((j + 1) - j) + 3 + 1) = 1
    rw [h_sub]
    have h_eq : (3 * 5 ^ 1 + 3 + 1 : ℕ) = 19 := by norm_num
    rw [h_eq]
    exact chiP_nineteen_eq_one
  -- (2,2): unfolding(2, 0) = chiP 11 = 1
  have hD22 : unfolding 5 (j + 1) j ⟨2, h_two_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 1 := by
    change chiP (2 * 5 ^ ((j + 1) - j) + 0 + 1) = 1
    rw [h_sub]
    have h_eq : (2 * 5 ^ 1 + 0 + 1 : ℕ) = 11 := by norm_num
    rw [h_eq]
    exact chiP_eleven_eq_one
  -- (3,3): unfolding(1, 1) = chiP 7 = 1
  have hD33 : unfolding 5 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 1 := by
    change chiP (1 * 5 ^ ((j + 1) - j) + 1 + 1) = 1
    rw [h_sub]
    have h_eq : (1 * 5 ^ 1 + 1 + 1 : ℕ) = 7 := by norm_num
    rw [h_eq]
    exact chiP_seven_eq_one
  -- (4,4): unfolding(4, 2) = chiP 23 = 1
  have hD44 : unfolding 5 (j + 1) j ⟨4, h_four_lt_pow_j⟩ ⟨2, h_two_lt_pow_dj⟩ = 1 := by
    change chiP (4 * 5 ^ ((j + 1) - j) + 2 + 1) = 1
    rw [h_sub]
    have h_eq : (4 * 5 ^ 1 + 2 + 1 : ℕ) = 23 := by norm_num
    rw [h_eq]
    exact chiP_twenty_three_eq_one
  -- Compute the 10 below-diagonal entries (all = 0).
  -- (1,0): unfolding(3, 4) = chiP 20 = 0
  have hL10 : unfolding 5 (j + 1) j ⟨3, h_three_lt_pow_j⟩ ⟨4, h_four_lt_pow_dj⟩ = 0 := by
    change chiP (3 * 5 ^ ((j + 1) - j) + 4 + 1) = 0
    rw [h_sub]
    have h_eq : (3 * 5 ^ 1 + 4 + 1 : ℕ) = 20 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_20]
  -- (2,0): unfolding(2, 4) = chiP 15 = 0
  have hL20 : unfolding 5 (j + 1) j ⟨2, h_two_lt_pow_j⟩ ⟨4, h_four_lt_pow_dj⟩ = 0 := by
    change chiP (2 * 5 ^ ((j + 1) - j) + 4 + 1) = 0
    rw [h_sub]
    have h_eq : (2 * 5 ^ 1 + 4 + 1 : ℕ) = 15 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_15]
  -- (2,1): unfolding(2, 3) = chiP 14 = 0
  have hL21 : unfolding 5 (j + 1) j ⟨2, h_two_lt_pow_j⟩ ⟨3, h_three_lt_pow_dj⟩ = 0 := by
    change chiP (2 * 5 ^ ((j + 1) - j) + 3 + 1) = 0
    rw [h_sub]
    have h_eq : (2 * 5 ^ 1 + 3 + 1 : ℕ) = 14 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_14]
  -- (3,0): unfolding(1, 4) = chiP 10 = 0
  have hL30 : unfolding 5 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨4, h_four_lt_pow_dj⟩ = 0 := by
    change chiP (1 * 5 ^ ((j + 1) - j) + 4 + 1) = 0
    rw [h_sub]
    have h_eq : (1 * 5 ^ 1 + 4 + 1 : ℕ) = 10 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_10]
  -- (3,1): unfolding(1, 3) = chiP 9 = 0
  have hL31 : unfolding 5 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨3, h_three_lt_pow_dj⟩ = 0 := by
    change chiP (1 * 5 ^ ((j + 1) - j) + 3 + 1) = 0
    rw [h_sub]
    have h_eq : (1 * 5 ^ 1 + 3 + 1 : ℕ) = 9 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_9]
  -- (3,2): unfolding(1, 0) = chiP 6 = 0
  have hL32 : unfolding 5 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (1 * 5 ^ ((j + 1) - j) + 0 + 1) = 0
    rw [h_sub]
    have h_eq : (1 * 5 ^ 1 + 0 + 1 : ℕ) = 6 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_6]
  -- (4,0): unfolding(4, 4) = chiP 25 = 0
  have hL40 : unfolding 5 (j + 1) j ⟨4, h_four_lt_pow_j⟩ ⟨4, h_four_lt_pow_dj⟩ = 0 := by
    change chiP (4 * 5 ^ ((j + 1) - j) + 4 + 1) = 0
    rw [h_sub]
    have h_eq : (4 * 5 ^ 1 + 4 + 1 : ℕ) = 25 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_25]
  -- (4,1): unfolding(4, 3) = chiP 24 = 0
  have hL41 : unfolding 5 (j + 1) j ⟨4, h_four_lt_pow_j⟩ ⟨3, h_three_lt_pow_dj⟩ = 0 := by
    change chiP (4 * 5 ^ ((j + 1) - j) + 3 + 1) = 0
    rw [h_sub]
    have h_eq : (4 * 5 ^ 1 + 3 + 1 : ℕ) = 24 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_24]
  -- (4,2): unfolding(4, 0) = chiP 21 = 0
  have hL42 : unfolding 5 (j + 1) j ⟨4, h_four_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (4 * 5 ^ ((j + 1) - j) + 0 + 1) = 0
    rw [h_sub]
    have h_eq : (4 * 5 ^ 1 + 0 + 1 : ℕ) = 21 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_21]
  -- (4,3): unfolding(4, 1) = chiP 22 = 0
  have hL43 : unfolding 5 (j + 1) j ⟨4, h_four_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (4 * 5 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (4 * 5 ^ 1 + 1 + 1 : ℕ) = 22 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_22]
  -- Helpers for `Fin 5` `.val` and the `if-then-else` ρ, σ.
  have h0 : (0 : Fin 5).val = 0 := rfl
  have h1 : (1 : Fin 5).val = 1 := rfl
  have h2 : (2 : Fin 5).val = 2 := rfl
  have h3 : (3 : Fin 5).val = 3 := rfl
  have h4 : (4 : Fin 5).val = 4 := rfl
  have hne_2_0 : (2 : ℕ) ≠ 0 := by decide
  have hne_2_1 : (2 : ℕ) ≠ 1 := by decide
  have hne_3_0 : (3 : ℕ) ≠ 0 := by decide
  have hne_3_1 : (3 : ℕ) ≠ 1 := by decide
  have hne_3_2 : (3 : ℕ) ≠ 2 := by decide
  have hne_4_0 : (4 : ℕ) ≠ 0 := by decide
  have hne_4_1 : (4 : ℕ) ≠ 1 := by decide
  have hne_4_2 : (4 : ℕ) ≠ 2 := by decide
  have hne_4_3 : (4 : ℕ) ≠ 3 := by decide
  -- Establish the upper-triangular property of the submatrix.
  set Mρ : Fin 5 → Fin (5 ^ j) :=
    fun i => if i.val = 0 then ⟨0, h_pow_j_pos⟩
             else if i.val = 1 then ⟨3, h_three_lt_pow_j⟩
             else if i.val = 2 then ⟨2, h_two_lt_pow_j⟩
             else if i.val = 3 then ⟨1, h_one_lt_pow_j⟩
             else ⟨4, h_four_lt_pow_j⟩ with hMρ_def
  set Mσ : Fin 5 → Fin (5 ^ ((j + 1) - j)) :=
    fun i => if i.val = 0 then ⟨4, h_four_lt_pow_dj⟩
             else if i.val = 1 then ⟨3, h_three_lt_pow_dj⟩
             else if i.val = 2 then ⟨0, h_zero_lt_pow_dj⟩
             else if i.val = 3 then ⟨1, h_one_lt_pow_dj⟩
             else ⟨2, h_two_lt_pow_dj⟩ with hMσ_def
  set M : Matrix (Fin 5) (Fin 5) ℚ :=
    (unfolding 5 (j + 1) j).submatrix Mρ Mσ with hM_def
  have h_blocktri : M.BlockTriangular id := by
    intro i k h_lt
    -- `id k < id i` reduces to `k.val < i.val`.
    simp only [id_eq, Fin.lt_def] at h_lt
    fin_cases i <;> fin_cases k <;>
      first
        | (exact absurd h_lt (by decide))
        | (simp only [hM_def, Matrix.submatrix_apply, hMρ_def, hMσ_def,
                       if_true, if_false,
                       Nat.one_ne_zero, hne_2_0, hne_2_1, hne_3_0, hne_3_1,
                       hne_3_2, hne_4_0, hne_4_1, hne_4_2, hne_4_3]
           first | exact hL10 | exact hL20 | exact hL21
                 | exact hL30 | exact hL31 | exact hL32
                 | exact hL40 | exact hL41 | exact hL42 | exact hL43)
  rw [Matrix.det_of_upperTriangular h_blocktri]
  -- Now expand the diagonal product over `Fin 5`.
  rw [Fin.prod_univ_five]
  -- Reduce `M k k` for `k ∈ {0, 1, 2, 3, 4}` to the precomputed diagonal entries.
  simp only [hM_def, Matrix.submatrix_apply, hMρ_def, hMσ_def,
             h0, h1, h2, h3, h4, if_true, if_false,
             Nat.one_ne_zero, hne_2_0, hne_2_1, hne_3_0, hne_3_1, hne_3_2,
             hne_4_0, hne_4_1, hne_4_2, hne_4_3]
  rw [hD00, hD11, hD22, hD33, hD44]
  -- Goal: `IsUnit (1 * 1 * 1 * 1 * 1 : ℚ)`.
  norm_num

/--
**Corner-case main theorem (W = 5, d = j + 1).** The unfolding rank is
exactly `5` for every `j ≥ 1`. Mirrors `mps_bond_dim_W_eq_4_d_eq_j_plus_1`
and is fully formalised in Lean (no `sorry`, no `axiom`).

**Upper-bound subtlety:** the matrix has `5^((j+1)-j) = 5` columns, so
`rank_le_width` gives `rank ≤ 5`, which is already sharp here (since
`R = min(5^j, 5) = 5`). We use `rank_le_width` directly without needing
the general `upper_bound` lemma.

Together with the `W ∈ {2, 3, 4}` corners, this is the fifth
unconditional instance of `mps_bond_dim` and the third instance over a
wheel `W ≥ 3`. The general `mps_bond_dim` still requires
`exists_invertible_submatrix` whose proof is the only remaining `sorry`
in this file.
-/
theorem mps_bond_dim_W_eq_5_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    (unfolding 5 (j + 1) j).rank = 5 := by
  apply Nat.le_antisymm
  · -- `≤ 5`: `unfolding 5 (j+1) j` has only `5^((j+1)-j) = 5^1 = 5` columns.
    have h := (unfolding 5 (j + 1) j).rank_le_width
    have h_eq : 5 ^ ((j + 1 : ℕ) - j) = 5 := by
      have h_sub : (j + 1 : ℕ) - j = 1 := by omega
      rw [h_sub]; norm_num
    linarith
  · -- `5 ≤ rank`: from the corner-case prime exhibit.
    obtain ⟨ρ, σ, hUnit⟩ := exists_invertible_submatrix_W_eq_5_d_eq_j_plus_1 j hj
    have h_eq : ((unfolding 5 (j + 1) j).submatrix ρ σ).rank = 5 := by
      have h := Matrix.rank_of_isUnit ((unfolding 5 (j + 1) j).submatrix ρ σ) hUnit
      simpa using h
    calc (5 : ℕ) = ((unfolding 5 (j + 1) j).submatrix ρ σ).rank := h_eq.symm
      _ ≤ (unfolding 5 (j + 1) j).rank := Matrix.rank_submatrix_le _ _ _

/-!
### Corner case: `W = 6`, `d = j + 1` — closed unconditionally without Bertrand

When `W = 6` and `d = j + 1` (so `d - j = 1`), the formula gives
`R = min(6^j, φ(6) · 6^0 + 1) = min(6^j, 3) = 3` for every `j ≥ 1`. The
matrix has `6^(d-j) = 6` columns; live columns are those with
`gcd(k+1, 6) = 1`, i.e. `k + 1 ∈ {1, 5}`, so `k ∈ {0, 4}` — exactly
`φ(6) = 2` of them. We pick the two live columns plus one dead column —
column `1` (dead since `gcd(2, 6) = 2`) — which contributes `chiP 2 = 1`
at row `0` and zero elsewhere on the chosen row set.

**Row choice subtlety (the key novelty over W ∈ {2, 3, 4, 5}).** Unlike
the smaller wheels, the first three rows of the W=6 slab are *not*
linearly independent: rows `1, 2, 3` of `chiP 7..24` all have the same
support pattern `(1, 0, 0, 0, 1, 0)` because each of `7..12`, `13..18`,
`19..24` contains exactly two primes, both at residues `1, 5 (mod 6)`.
We therefore pick rows `{0, 1, 5}` — the row-`5` window `chiP 31..36`
has support pattern `(1, 0, 0, 0, 0, 0)` because only `31` is prime
among `31..36`, and that prime sits at residue `1 (mod 6)` (column `0`),
not `5 (mod 6)` (column `4`). This is the **first orthogonal-corner
instance where the working row set is not `{0, 1, …, R-1}`**.

Concretely `ρ : Fin 3 → Fin (6^j)` maps `(0, 1, 2) ↦ (0, 1, 5)` and
`σ : Fin 3 → Fin (6^((j+1)-j))` maps `(0, 1, 2) ↦ (0, 1, 4)`. The
`3 × 3` submatrix is
```
   ⎡ chiP  1, chiP  2, chiP  5 ⎤   ⎡ 0, 1, 1 ⎤
   ⎢ chiP  7, chiP  8, chiP 11 ⎥ = ⎢ 1, 0, 1 ⎥
   ⎣ chiP 31, chiP 32, chiP 35 ⎦   ⎣ 1, 0, 0 ⎦
```
of determinant `+1` (computed via `Matrix.det_fin_three`). The unit
witness is `isUnit_one`.

**Upper-bound subtlety.** As with W=4, `rank_le_width` only gives
`rank ≤ 6`, not the sharp `rank ≤ 3 = φ(6) + 1`. We cite the general
`upper_bound` lemma, which evaluates to `φ(6) · 6^0 + 1 = 2 · 1 + 1 = 3`
in this corner.

Used helpers: `chiP_two_eq_one`, `chiP_five_eq_one`, `chiP_seven_eq_one`,
`chiP_eleven_eq_one`, `chiP_thirty_one_eq_one` (new at S121), and
decidability of non-primality for `1, 8, 32, 35`. The prime `31` is the
smallest "row-5 prime" needed in the W=6 slab. **Sixth unconditional
`mps_bond_dim` instance; fourth instance over a wheel `W ≥ 3`.** The
general `mps_bond_dim` still requires `exists_invertible_submatrix`
whose proof is the only remaining `sorry` in this file.
-/

/--
`chiP 31 = 1` since `31` is prime.
-/
theorem chiP_thirty_one_eq_one : chiP 31 = 1 := by
  have h_prime_31 : Nat.Prime 31 := by decide
  simp [chiP, h_prime_31]

/--
**Corner-case prime exhibit (W = 6, d = j + 1).**

Specialisation of `exists_invertible_submatrix` to the corner case
`W = 6`, `d = j + 1`, where `R = min(6^j, 3) = 3` for every `j ≥ 1`.
Closed unconditionally using the chiP-helpers `chiP_two_eq_one`,
`chiP_five_eq_one`, `chiP_seven_eq_one`, `chiP_eleven_eq_one`, and
`chiP_thirty_one_eq_one` (all `decide`-checkable), plus decidability of
non-primality for `1, 8, 32, 35`.

The construction picks rows `(0, 1, 5)` and live columns `(0, 1, 4)`
(i.e. residues `1, 2, 5 (mod 6)`); the row choice skips rows `2, 3, 4`
because their `chiP`-pattern is identical to row `1`'s. The `3 × 3`
submatrix has determinant `+1`, hence `IsUnit` over `ℚ`.
-/
theorem exists_invertible_submatrix_W_eq_6_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    ∃ (ρ : Fin 3 → Fin (6 ^ j))
      (σ : Fin 3 → Fin (6 ^ ((j + 1) - j))),
      IsUnit ((unfolding 6 (j + 1) j).submatrix ρ σ) := by
  -- The exponent simplifies: `(j + 1) - j = 1`, hence `6 ^ ((j + 1) - j) = 6`.
  have h_sub : (j + 1 : ℕ) - j = 1 := by omega
  have h_pow_j_pos : 0 < 6 ^ j := Nat.pow_pos (by norm_num)
  have h_one_lt_pow_j : 1 < 6 ^ j := by
    calc (1 : ℕ) < 6 := by norm_num
      _ = 6 ^ 1 := by norm_num
      _ ≤ 6 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_five_lt_pow_j : 5 < 6 ^ j := by
    calc (5 : ℕ) < 6 := by norm_num
      _ = 6 ^ 1 := by norm_num
      _ ≤ 6 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_zero_lt_pow_dj : (0 : ℕ) < 6 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_one_lt_pow_dj : (1 : ℕ) < 6 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_four_lt_pow_dj : (4 : ℕ) < 6 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  -- Build `ρ` mapping `(0, 1, 2) ↦ (0, 1, 5)` and `σ` mapping
  -- `(0, 1, 2) ↦ (0, 1, 4)`.
  refine ⟨fun i => if i.val = 0 then ⟨0, h_pow_j_pos⟩
                   else if i.val = 1 then ⟨1, h_one_lt_pow_j⟩
                   else ⟨5, h_five_lt_pow_j⟩,
          fun i => if i.val = 0 then ⟨0, h_zero_lt_pow_dj⟩
                   else if i.val = 1 then ⟨1, h_one_lt_pow_dj⟩
                   else ⟨4, h_four_lt_pow_dj⟩, ?_⟩
  -- Reduce `IsUnit` of the `3 × 3` submatrix to `IsUnit` of its determinant.
  rw [Matrix.isUnit_iff_isUnit_det, Matrix.det_fin_three]
  -- Non-primality lemmas for the zero entries.
  have h_not_prime_8 : ¬ Nat.Prime 8 := by decide
  have h_not_prime_32 : ¬ Nat.Prime 32 := by decide
  have h_not_prime_35 : ¬ Nat.Prime 35 := by decide
  -- Compute the nine entries of the submatrix.
  -- Entry (0, 0): chiP(0 · 6^((j+1)-j) + 0 + 1) = chiP 1 = 0.
  have hM00 : unfolding 6 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (0 * 6 ^ ((j + 1) - j) + 0 + 1) = 0
    simp [chiP, Nat.not_prime_one]
  -- Entry (0, 1): chiP(0 · 6^((j+1)-j) + 1 + 1) = chiP 2 = 1.
  have hM01 : unfolding 6 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨1, h_one_lt_pow_dj⟩ = 1 := by
    change chiP (0 * 6 ^ ((j + 1) - j) + 1 + 1) = 1
    simp [chiP, Nat.prime_two]
  -- Entry (0, 2): chiP(0 · 6^((j+1)-j) + 4 + 1) = chiP 5 = 1.
  have hM02 : unfolding 6 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨4, h_four_lt_pow_dj⟩ = 1 := by
    change chiP (0 * 6 ^ ((j + 1) - j) + 4 + 1) = 1
    have h_eq : (0 * 6 ^ ((j + 1) - j) + 4 + 1 : ℕ) = 5 := by simp
    rw [h_eq]
    exact chiP_five_eq_one
  -- Entry (1, 0): chiP(1 · 6^((j+1)-j) + 0 + 1) = chiP 7 = 1.
  have hM10 : unfolding 6 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 1 := by
    change chiP (1 * 6 ^ ((j + 1) - j) + 0 + 1) = 1
    rw [h_sub]
    have h_eq : (1 * 6 ^ 1 + 0 + 1 : ℕ) = 7 := by norm_num
    rw [h_eq]
    exact chiP_seven_eq_one
  -- Entry (1, 1): chiP(1 · 6^((j+1)-j) + 1 + 1) = chiP 8 = 0.
  have hM11 : unfolding 6 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (1 * 6 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (1 * 6 ^ 1 + 1 + 1 : ℕ) = 8 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_8]
  -- Entry (1, 2): chiP(1 · 6^((j+1)-j) + 4 + 1) = chiP 11 = 1.
  have hM12 : unfolding 6 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨4, h_four_lt_pow_dj⟩ = 1 := by
    change chiP (1 * 6 ^ ((j + 1) - j) + 4 + 1) = 1
    rw [h_sub]
    have h_eq : (1 * 6 ^ 1 + 4 + 1 : ℕ) = 11 := by norm_num
    rw [h_eq]
    exact chiP_eleven_eq_one
  -- Entry (2, 0): chiP(5 · 6^((j+1)-j) + 0 + 1) = chiP 31 = 1.
  have hM20 : unfolding 6 (j + 1) j ⟨5, h_five_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 1 := by
    change chiP (5 * 6 ^ ((j + 1) - j) + 0 + 1) = 1
    rw [h_sub]
    have h_eq : (5 * 6 ^ 1 + 0 + 1 : ℕ) = 31 := by norm_num
    rw [h_eq]
    exact chiP_thirty_one_eq_one
  -- Entry (2, 1): chiP(5 · 6^((j+1)-j) + 1 + 1) = chiP 32 = 0.
  have hM21 : unfolding 6 (j + 1) j ⟨5, h_five_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (5 * 6 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (5 * 6 ^ 1 + 1 + 1 : ℕ) = 32 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_32]
  -- Entry (2, 2): chiP(5 · 6^((j+1)-j) + 4 + 1) = chiP 35 = 0.
  have hM22 : unfolding 6 (j + 1) j ⟨5, h_five_lt_pow_j⟩ ⟨4, h_four_lt_pow_dj⟩ = 0 := by
    change chiP (5 * 6 ^ ((j + 1) - j) + 4 + 1) = 0
    rw [h_sub]
    have h_eq : (5 * 6 ^ 1 + 4 + 1 : ℕ) = 35 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_35]
  -- Reduce the `if-then-else` ρ, σ at the nine cells via `Fin.val` of `0, 1, 2`.
  have h0 : (0 : Fin 3).val = 0 := rfl
  have h1 : (1 : Fin 3).val = 1 := rfl
  have h2 : (2 : Fin 3).val = 2 := rfl
  have hne_2_0 : (2 : ℕ) ≠ 0 := by decide
  have hne_2_1 : (2 : ℕ) ≠ 1 := by decide
  simp only [Matrix.submatrix_apply, h0, h1, h2, if_true, if_false,
             Nat.one_ne_zero, hne_2_0, hne_2_1]
  rw [hM00, hM01, hM02, hM10, hM11, hM12, hM20, hM21, hM22]
  -- Goal: `IsUnit (0·0·0 - 0·1·0 - 1·1·0 + 1·1·1 + 1·1·0 - 1·0·1)` over `ℚ`.
  --     = `IsUnit (1 : ℚ)`.
  ring_nf
  exact isUnit_one

/--
**Corner-case main theorem (W = 6, d = j + 1).** The unfolding rank is
exactly `3` for every `j ≥ 1`. Mirrors `mps_bond_dim_W_eq_4_d_eq_j_plus_1`
(both have `R = 3`) and is fully formalised in Lean (no `sorry`,
no `axiom`).

**Upper-bound subtlety:** the matrix has `6^((j+1)-j) = 6` columns, so
`rank_le_width` only gives `rank ≤ 6`. We instead cite the general
`upper_bound` lemma, which evaluates to `φ(6) · 6^0 + 1 = 2 · 1 + 1 = 3`
in this corner.

Together with the `W ∈ {2, 3, 4, 5}` corners, this is the sixth
unconditional instance of `mps_bond_dim` and the fourth instance over a
wheel `W ≥ 3`. The general `mps_bond_dim` still requires
`exists_invertible_submatrix` whose proof is the only remaining `sorry`
in this file.
-/
theorem mps_bond_dim_W_eq_6_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    (unfolding 6 (j + 1) j).rank = 3 := by
  apply Nat.le_antisymm
  · -- `≤ 3`: from `upper_bound`. (`rank_le_width` only gives `≤ 6`.)
    have hW : (2 : ℕ) ≤ 6 := by norm_num
    have hj_hi : j < j + 1 := Nat.lt_succ_self j
    have h := upper_bound 6 (j + 1) j hW hj hj_hi
    have h_eq : Nat.totient 6 * 6 ^ ((j + 1 : ℕ) - j - 1) + 1 = 3 := by
      have h_sub : (j + 1 : ℕ) - j - 1 = 0 := by omega
      rw [h_sub]
      decide
    linarith
  · -- `3 ≤ rank`: from the corner-case prime exhibit.
    obtain ⟨ρ, σ, hUnit⟩ := exists_invertible_submatrix_W_eq_6_d_eq_j_plus_1 j hj
    have h_eq : ((unfolding 6 (j + 1) j).submatrix ρ σ).rank = 3 := by
      have h := Matrix.rank_of_isUnit ((unfolding 6 (j + 1) j).submatrix ρ σ) hUnit
      simpa using h
    calc (3 : ℕ) = ((unfolding 6 (j + 1) j).submatrix ρ σ).rank := h_eq.symm
      _ ≤ (unfolding 6 (j + 1) j).rank := Matrix.rank_submatrix_le _ _ _

/-!
### Corner case: `W = 8`, `d = j + 1` — closed unconditionally without Bertrand

When `W = 8` and `d = j + 1` (so `d - j = 1`), the formula gives
`R = min(8^j, φ(8) · 8^0 + 1) = min(8^j, 5) = 5` for every `j ≥ 1`. The
matrix has `8^(d-j) = 8` columns; live columns are those with
`gcd(k+1, 8) = 1`, i.e. `k + 1 ∈ {1, 3, 5, 7}` (`k + 1` odd), so
`k ∈ {0, 2, 4, 6}` — exactly `φ(8) = 4` of them. We pick the four live
columns plus one dead column — column `1`, which is the unique dead
column whose `chiP` value at row `0` is `1` (since `chiP 2 = 1` and
the other dead columns `{3, 5, 7}` yield `chiP 4, 6, 8` — all zero).

**Triangulation strategy (BlockTriangular route, à la W=5).** The full
`5 × 8` row-{0..4} restricted matrix at columns `{0, 1, 2, 4, 6}` admits
a permutation `(ρ, σ)` that triangularises it:
```
  Original 5×5 submatrix at rows {0,1,2,3,4}, cols {0,1,2,4,6}:
     rows  0  1  2  3  4
     col 0   0, 0, 1, 0, 0
     col 1   1, 0, 0, 0, 0
     col 2   1, 1, 1, 0, 0
     col 4   1, 1, 0, 1, 1
     col 6   1, 0, 1, 1, 0

  After ρ : Fin 5 → Fin (8^j) mapping `(0, 1, 2, 3, 4) ↦ (2, 0, 1, 3, 4)`
  and σ : Fin 5 → Fin (8^((j+1)-j)) mapping `(0, 1, 2, 3, 4) ↦ (0, 1, 2, 6, 4)`:

     ⎡ chiP 17, chiP 18, chiP 19, chiP 23, chiP 21 ⎤   ⎡ 1, 0, 1, 1, 0 ⎤
     ⎢ chiP  1, chiP  2, chiP  3, chiP  7, chiP  5 ⎥   ⎢ 0, 1, 1, 1, 1 ⎥
     ⎢ chiP  9, chiP 10, chiP 11, chiP 15, chiP 13 ⎥ = ⎢ 0, 0, 1, 0, 1 ⎥
     ⎢ chiP 25, chiP 26, chiP 27, chiP 31, chiP 29 ⎥   ⎢ 0, 0, 0, 1, 1 ⎥
     ⎣ chiP 33, chiP 34, chiP 35, chiP 39, chiP 37 ⎦   ⎣ 0, 0, 0, 0, 1 ⎦
```
i.e. upper-triangular with `1` on the diagonal. The five diagonal
primes are `{17, 2, 11, 31, 37}`; the ten below-diagonal entries are
`chiP` of `{1, 9, 10, 25, 26, 27, 33, 34, 35, 39}`, all composite.

**Determinant via `det_of_upperTriangular`** (mathlib has no `det_fin_five`,
same situation as W=5). The proof mirrors `W = 5` exactly: we precompute
the 5 diagonal entries (each `= 1`) and the 10 below-diagonal entries
(each `= 0`), establish `BlockTriangular id` by `fin_cases i <;> fin_cases k`,
apply `det_of_upperTriangular`, expand the diagonal product via
`Fin.prod_univ_five`, and finish with `norm_num`.

**Upper-bound subtlety.** As with W=4 and W=6, `rank_le_width` gives only
`rank ≤ 8`, not the sharp `rank ≤ 5 = φ(8) + 1`. We cite the general
`upper_bound`, which evaluates to `φ(8) · 8^0 + 1 = 4 · 1 + 1 = 5`.

Used helpers: `chiP_two_eq_one`, `chiP_eleven_eq_one`,
`chiP_thirty_one_eq_one`, `chiP_seventeen_eq_one` (new at S128),
`chiP_thirty_seven_eq_one` (new at S128), and decidability of
non-primality for `1, 9, 10, 25, 26, 27, 33, 34, 35, 39`.
-/

/--
`chiP 17 = 1` since `17` is prime.
-/
theorem chiP_seventeen_eq_one : chiP 17 = 1 := by
  have h_prime_17 : Nat.Prime 17 := by decide
  simp [chiP, h_prime_17]

/--
`chiP 37 = 1` since `37` is prime.
-/
theorem chiP_thirty_seven_eq_one : chiP 37 = 1 := by
  have h_prime_37 : Nat.Prime 37 := by decide
  simp [chiP, h_prime_37]

/--
**Corner-case prime exhibit (W = 8, d = j + 1).**

Specialisation of `exists_invertible_submatrix` to the corner case
`W = 8`, `d = j + 1`, where `R = min(8^j, 5) = 5` for every `j ≥ 1`.
Closed unconditionally using `chiP_two_eq_one`, `chiP_eleven_eq_one`,
`chiP_thirty_one_eq_one`, `chiP_seventeen_eq_one`, and
`chiP_thirty_seven_eq_one`.

The construction uses `ρ : Fin 5 → Fin (8^j)` mapping
`(0, 1, 2, 3, 4) ↦ (2, 0, 1, 3, 4)` and `σ : Fin 5 → Fin (8^((j+1)-j))`
mapping `(0, 1, 2, 3, 4) ↦ (0, 1, 2, 6, 4)`. The `5 × 5` submatrix is
upper triangular with `1` on the diagonal, hence `IsUnit` over `ℚ`.
-/
theorem exists_invertible_submatrix_W_eq_8_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    ∃ (ρ : Fin 5 → Fin (8 ^ j))
      (σ : Fin 5 → Fin (8 ^ ((j + 1) - j))),
      IsUnit ((unfolding 8 (j + 1) j).submatrix ρ σ) := by
  -- The exponent simplifies: `(j + 1) - j = 1`, hence `8 ^ ((j + 1) - j) = 8`.
  have h_sub : (j + 1 : ℕ) - j = 1 := by omega
  have h_pow_j_pos : 0 < 8 ^ j := Nat.pow_pos (by norm_num)
  have h_one_lt_pow_j : 1 < 8 ^ j := by
    calc (1 : ℕ) < 8 := by norm_num
      _ = 8 ^ 1 := by norm_num
      _ ≤ 8 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_two_lt_pow_j : 2 < 8 ^ j := by
    calc (2 : ℕ) < 8 := by norm_num
      _ = 8 ^ 1 := by norm_num
      _ ≤ 8 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_three_lt_pow_j : 3 < 8 ^ j := by
    calc (3 : ℕ) < 8 := by norm_num
      _ = 8 ^ 1 := by norm_num
      _ ≤ 8 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_four_lt_pow_j : 4 < 8 ^ j := by
    calc (4 : ℕ) < 8 := by norm_num
      _ = 8 ^ 1 := by norm_num
      _ ≤ 8 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_zero_lt_pow_dj : (0 : ℕ) < 8 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_one_lt_pow_dj : (1 : ℕ) < 8 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_two_lt_pow_dj : (2 : ℕ) < 8 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_four_lt_pow_dj : (4 : ℕ) < 8 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_six_lt_pow_dj : (6 : ℕ) < 8 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  -- ρ permutes original rows in the order `(2, 0, 1, 3, 4)`.
  -- σ permutes original columns in the order `(0, 1, 2, 6, 4)`.
  refine ⟨fun i => if i.val = 0 then ⟨2, h_two_lt_pow_j⟩
                   else if i.val = 1 then ⟨0, h_pow_j_pos⟩
                   else if i.val = 2 then ⟨1, h_one_lt_pow_j⟩
                   else if i.val = 3 then ⟨3, h_three_lt_pow_j⟩
                   else ⟨4, h_four_lt_pow_j⟩,
          fun i => if i.val = 0 then ⟨0, h_zero_lt_pow_dj⟩
                   else if i.val = 1 then ⟨1, h_one_lt_pow_dj⟩
                   else if i.val = 2 then ⟨2, h_two_lt_pow_dj⟩
                   else if i.val = 3 then ⟨6, h_six_lt_pow_dj⟩
                   else ⟨4, h_four_lt_pow_dj⟩, ?_⟩
  -- Reduce `IsUnit` to `IsUnit (det)` and compute the det as the
  -- product of the diagonal via `Matrix.det_of_upperTriangular`.
  rw [Matrix.isUnit_iff_isUnit_det]
  -- Non-primality lemmas for the ten below-diagonal zeros.
  have h_not_prime_9  : ¬ Nat.Prime 9  := by decide
  have h_not_prime_10 : ¬ Nat.Prime 10 := by decide
  have h_not_prime_25 : ¬ Nat.Prime 25 := by decide
  have h_not_prime_26 : ¬ Nat.Prime 26 := by decide
  have h_not_prime_27 : ¬ Nat.Prime 27 := by decide
  have h_not_prime_33 : ¬ Nat.Prime 33 := by decide
  have h_not_prime_34 : ¬ Nat.Prime 34 := by decide
  have h_not_prime_35 : ¬ Nat.Prime 35 := by decide
  have h_not_prime_39 : ¬ Nat.Prime 39 := by decide
  -- Compute the 5 diagonal entries (all = 1).
  -- (0,0): unfolding(2, 0) = chiP 17 = 1
  have hD00 : unfolding 8 (j + 1) j ⟨2, h_two_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 1 := by
    change chiP (2 * 8 ^ ((j + 1) - j) + 0 + 1) = 1
    rw [h_sub]
    have h_eq : (2 * 8 ^ 1 + 0 + 1 : ℕ) = 17 := by norm_num
    rw [h_eq]
    exact chiP_seventeen_eq_one
  -- (1,1): unfolding(0, 1) = chiP 2 = 1
  have hD11 : unfolding 8 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨1, h_one_lt_pow_dj⟩ = 1 := by
    change chiP (0 * 8 ^ ((j + 1) - j) + 1 + 1) = 1
    simp [chiP, Nat.prime_two]
  -- (2,2): unfolding(1, 2) = chiP 11 = 1
  have hD22 : unfolding 8 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨2, h_two_lt_pow_dj⟩ = 1 := by
    change chiP (1 * 8 ^ ((j + 1) - j) + 2 + 1) = 1
    rw [h_sub]
    have h_eq : (1 * 8 ^ 1 + 2 + 1 : ℕ) = 11 := by norm_num
    rw [h_eq]
    exact chiP_eleven_eq_one
  -- (3,3): unfolding(3, 6) = chiP 31 = 1
  have hD33 : unfolding 8 (j + 1) j ⟨3, h_three_lt_pow_j⟩ ⟨6, h_six_lt_pow_dj⟩ = 1 := by
    change chiP (3 * 8 ^ ((j + 1) - j) + 6 + 1) = 1
    rw [h_sub]
    have h_eq : (3 * 8 ^ 1 + 6 + 1 : ℕ) = 31 := by norm_num
    rw [h_eq]
    exact chiP_thirty_one_eq_one
  -- (4,4): unfolding(4, 4) = chiP 37 = 1
  have hD44 : unfolding 8 (j + 1) j ⟨4, h_four_lt_pow_j⟩ ⟨4, h_four_lt_pow_dj⟩ = 1 := by
    change chiP (4 * 8 ^ ((j + 1) - j) + 4 + 1) = 1
    rw [h_sub]
    have h_eq : (4 * 8 ^ 1 + 4 + 1 : ℕ) = 37 := by norm_num
    rw [h_eq]
    exact chiP_thirty_seven_eq_one
  -- Compute the 10 below-diagonal entries (all = 0).
  -- (1,0): unfolding(0, 0) = chiP 1 = 0
  have hL10 : unfolding 8 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (0 * 8 ^ ((j + 1) - j) + 0 + 1) = 0
    simp [chiP, Nat.not_prime_one]
  -- (2,0): unfolding(1, 0) = chiP 9 = 0
  have hL20 : unfolding 8 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (1 * 8 ^ ((j + 1) - j) + 0 + 1) = 0
    rw [h_sub]
    have h_eq : (1 * 8 ^ 1 + 0 + 1 : ℕ) = 9 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_9]
  -- (2,1): unfolding(1, 1) = chiP 10 = 0
  have hL21 : unfolding 8 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (1 * 8 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (1 * 8 ^ 1 + 1 + 1 : ℕ) = 10 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_10]
  -- (3,0): unfolding(3, 0) = chiP 25 = 0
  have hL30 : unfolding 8 (j + 1) j ⟨3, h_three_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (3 * 8 ^ ((j + 1) - j) + 0 + 1) = 0
    rw [h_sub]
    have h_eq : (3 * 8 ^ 1 + 0 + 1 : ℕ) = 25 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_25]
  -- (3,1): unfolding(3, 1) = chiP 26 = 0
  have hL31 : unfolding 8 (j + 1) j ⟨3, h_three_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (3 * 8 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (3 * 8 ^ 1 + 1 + 1 : ℕ) = 26 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_26]
  -- (3,2): unfolding(3, 2) = chiP 27 = 0
  have hL32 : unfolding 8 (j + 1) j ⟨3, h_three_lt_pow_j⟩ ⟨2, h_two_lt_pow_dj⟩ = 0 := by
    change chiP (3 * 8 ^ ((j + 1) - j) + 2 + 1) = 0
    rw [h_sub]
    have h_eq : (3 * 8 ^ 1 + 2 + 1 : ℕ) = 27 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_27]
  -- (4,0): unfolding(4, 0) = chiP 33 = 0
  have hL40 : unfolding 8 (j + 1) j ⟨4, h_four_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (4 * 8 ^ ((j + 1) - j) + 0 + 1) = 0
    rw [h_sub]
    have h_eq : (4 * 8 ^ 1 + 0 + 1 : ℕ) = 33 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_33]
  -- (4,1): unfolding(4, 1) = chiP 34 = 0
  have hL41 : unfolding 8 (j + 1) j ⟨4, h_four_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (4 * 8 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (4 * 8 ^ 1 + 1 + 1 : ℕ) = 34 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_34]
  -- (4,2): unfolding(4, 2) = chiP 35 = 0
  have hL42 : unfolding 8 (j + 1) j ⟨4, h_four_lt_pow_j⟩ ⟨2, h_two_lt_pow_dj⟩ = 0 := by
    change chiP (4 * 8 ^ ((j + 1) - j) + 2 + 1) = 0
    rw [h_sub]
    have h_eq : (4 * 8 ^ 1 + 2 + 1 : ℕ) = 35 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_35]
  -- (4,3): unfolding(4, 6) = chiP 39 = 0
  have hL43 : unfolding 8 (j + 1) j ⟨4, h_four_lt_pow_j⟩ ⟨6, h_six_lt_pow_dj⟩ = 0 := by
    change chiP (4 * 8 ^ ((j + 1) - j) + 6 + 1) = 0
    rw [h_sub]
    have h_eq : (4 * 8 ^ 1 + 6 + 1 : ℕ) = 39 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_39]
  -- Helpers for `Fin 5` `.val` and the `if-then-else` ρ, σ.
  have h0 : (0 : Fin 5).val = 0 := rfl
  have h1 : (1 : Fin 5).val = 1 := rfl
  have h2 : (2 : Fin 5).val = 2 := rfl
  have h3 : (3 : Fin 5).val = 3 := rfl
  have h4 : (4 : Fin 5).val = 4 := rfl
  have hne_2_0 : (2 : ℕ) ≠ 0 := by decide
  have hne_2_1 : (2 : ℕ) ≠ 1 := by decide
  have hne_3_0 : (3 : ℕ) ≠ 0 := by decide
  have hne_3_1 : (3 : ℕ) ≠ 1 := by decide
  have hne_3_2 : (3 : ℕ) ≠ 2 := by decide
  have hne_4_0 : (4 : ℕ) ≠ 0 := by decide
  have hne_4_1 : (4 : ℕ) ≠ 1 := by decide
  have hne_4_2 : (4 : ℕ) ≠ 2 := by decide
  have hne_4_3 : (4 : ℕ) ≠ 3 := by decide
  -- Establish the upper-triangular property of the submatrix.
  set Mρ : Fin 5 → Fin (8 ^ j) :=
    fun i => if i.val = 0 then ⟨2, h_two_lt_pow_j⟩
             else if i.val = 1 then ⟨0, h_pow_j_pos⟩
             else if i.val = 2 then ⟨1, h_one_lt_pow_j⟩
             else if i.val = 3 then ⟨3, h_three_lt_pow_j⟩
             else ⟨4, h_four_lt_pow_j⟩ with hMρ_def
  set Mσ : Fin 5 → Fin (8 ^ ((j + 1) - j)) :=
    fun i => if i.val = 0 then ⟨0, h_zero_lt_pow_dj⟩
             else if i.val = 1 then ⟨1, h_one_lt_pow_dj⟩
             else if i.val = 2 then ⟨2, h_two_lt_pow_dj⟩
             else if i.val = 3 then ⟨6, h_six_lt_pow_dj⟩
             else ⟨4, h_four_lt_pow_dj⟩ with hMσ_def
  set M : Matrix (Fin 5) (Fin 5) ℚ :=
    (unfolding 8 (j + 1) j).submatrix Mρ Mσ with hM_def
  have h_blocktri : M.BlockTriangular id := by
    intro i k h_lt
    -- `id k < id i` reduces to `k.val < i.val`.
    simp only [id_eq, Fin.lt_def] at h_lt
    fin_cases i <;> fin_cases k <;>
      first
        | (exact absurd h_lt (by decide))
        | (simp only [hM_def, Matrix.submatrix_apply, hMρ_def, hMσ_def,
                       if_true, if_false,
                       Nat.one_ne_zero, hne_2_0, hne_2_1, hne_3_0, hne_3_1,
                       hne_3_2, hne_4_0, hne_4_1, hne_4_2, hne_4_3]
           first | exact hL10 | exact hL20 | exact hL21
                 | exact hL30 | exact hL31 | exact hL32
                 | exact hL40 | exact hL41 | exact hL42 | exact hL43)
  rw [Matrix.det_of_upperTriangular h_blocktri]
  -- Now expand the diagonal product over `Fin 5`.
  rw [Fin.prod_univ_five]
  -- Reduce `M k k` for `k ∈ {0, 1, 2, 3, 4}` to the precomputed diagonal entries.
  simp only [hM_def, Matrix.submatrix_apply, hMρ_def, hMσ_def,
             h0, h1, h2, h3, h4, if_true, if_false,
             Nat.one_ne_zero, hne_2_0, hne_2_1, hne_3_0, hne_3_1, hne_3_2,
             hne_4_0, hne_4_1, hne_4_2, hne_4_3]
  rw [hD00, hD11, hD22, hD33, hD44]
  -- Goal: `IsUnit (1 * 1 * 1 * 1 * 1 : ℚ)`.
  norm_num

/--
**Corner-case main theorem (W = 8, d = j + 1).** The unfolding rank is
exactly `5` for every `j ≥ 1`. Mirrors `mps_bond_dim_W_eq_4_d_eq_j_plus_1`
and `mps_bond_dim_W_eq_6_d_eq_j_plus_1` (both citing the general
`upper_bound`) and is fully formalised in Lean (no `sorry`, no `axiom`).

**Upper-bound subtlety:** the matrix has `8^((j+1)-j) = 8` columns, so
`rank_le_width` only gives `rank ≤ 8`. We instead cite the general
`upper_bound` lemma, which evaluates to `φ(8) · 8^0 + 1 = 4 · 1 + 1 = 5`
in this corner.

Together with the `W ∈ {2, 3, 4, 5, 6}` corners, this is the seventh
unconditional instance of `mps_bond_dim` and the fifth instance over a
wheel `W ≥ 3`. The general `mps_bond_dim` still requires
`exists_invertible_submatrix` whose proof is the only remaining `sorry`
in this file.
-/
theorem mps_bond_dim_W_eq_8_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    (unfolding 8 (j + 1) j).rank = 5 := by
  apply Nat.le_antisymm
  · -- `≤ 5`: from `upper_bound`. (`rank_le_width` only gives `≤ 8`.)
    have hW : (2 : ℕ) ≤ 8 := by norm_num
    have hj_hi : j < j + 1 := Nat.lt_succ_self j
    have h := upper_bound 8 (j + 1) j hW hj hj_hi
    have h_eq : Nat.totient 8 * 8 ^ ((j + 1 : ℕ) - j - 1) + 1 = 5 := by
      have h_sub : (j + 1 : ℕ) - j - 1 = 0 := by omega
      rw [h_sub]
      decide
    linarith
  · -- `5 ≤ rank`: from the corner-case prime exhibit.
    obtain ⟨ρ, σ, hUnit⟩ := exists_invertible_submatrix_W_eq_8_d_eq_j_plus_1 j hj
    have h_eq : ((unfolding 8 (j + 1) j).submatrix ρ σ).rank = 5 := by
      have h := Matrix.rank_of_isUnit ((unfolding 8 (j + 1) j).submatrix ρ σ) hUnit
      simpa using h
    calc (5 : ℕ) = ((unfolding 8 (j + 1) j).submatrix ρ σ).rank := h_eq.symm
      _ ≤ (unfolding 8 (j + 1) j).rank := Matrix.rank_submatrix_le _ _ _

/-!
### Corner case: `W = 12`, `d = j + 1` — closed unconditionally

When `W = 12` and `d = j + 1`, the formula gives
`R = min(12^j, φ(12) * 12^0 + 1) = min(12^j, 5) = 5` for every `j ≥ 1`
(since `12^j ≥ 12 ≥ 5` and `φ(12) = 4`).

**Construction.** Live columns at `W = 12` are
`k ∈ {0, 4, 6, 10}` (residues `1, 5, 7, 11 (mod 12)`, where `gcd(k+1, 12) = 1`),
giving `φ(12) = 4` live columns. We need `R = φ(12) + 1 = 5` columns total,
so we add one **dead column** with `chiP` non-zero at row `0`. Two dead
candidates are available (`k = 1`: `chiP 2 = 1`, and `k = 2`: `chiP 3 = 1`);
we pick `k = 1` (the dead column with `chiP 2 = 1`) — paralleling the
W=8 choice of `k = 1` with `chiP 2 = 1`.

The five chosen columns `{1, 0, 6, 10, 4}` (in this order) and five
chosen rows `{0, 9, 10, 4, 7}` (also in this order) form the `5 × 5`
upper-triangular submatrix
```
     col   k+1=2 k+1=1 k+1=7 k+1=11 k+1=5
   row 0 ⎡ chiP 2,  chiP 1,  chiP 7,  chiP 11, chiP 5  ⎤   ⎡ 1, 0, 1, 1, 1 ⎤
   row 9 ⎢ chiP 110,chiP 109,chiP 115,chiP 119,chiP 113⎥ = ⎢ 0, 1, 0, 0, 1 ⎥
   row 10⎢ chiP 122,chiP 121,chiP 127,chiP 131,chiP 125⎥   ⎢ 0, 0, 1, 1, 0 ⎥
   row 4 ⎢ chiP 50, chiP 49, chiP 55, chiP 59, chiP 53 ⎥   ⎢ 0, 0, 0, 1, 1 ⎥
   row 7 ⎣ chiP 86, chiP 85, chiP 91, chiP 95, chiP 89 ⎦   ⎣ 0, 0, 0, 0, 1 ⎦
```
i.e. upper-triangular with `1` on the diagonal. The five diagonal
primes are `{2, 109, 127, 59, 89}`; the ten below-diagonal entries are
`chiP` of `{110, 122, 121, 50, 49, 55, 86, 85, 91, 95}`, all composite.

**Determinant via `det_of_upperTriangular`** (mathlib has no `det_fin_five`,
same situation as W=5 and W=8). The proof mirrors `W = 8` exactly: we
precompute the 5 diagonal entries (each `= 1`) and the 10 below-
diagonal entries (each `= 0`), establish `BlockTriangular id` by
`fin_cases i <;> fin_cases k`, apply `det_of_upperTriangular`, expand
the diagonal product via `Fin.prod_univ_five`, and finish with `norm_num`.

**Upper-bound subtlety.** As with W=4, W=6, W=8, `rank_le_width` gives only
`rank ≤ 12`, not the sharp `rank ≤ 5 = φ(12) + 1`. We cite the general
`upper_bound`, which evaluates to `φ(12) · 12^0 + 1 = 4 · 1 + 1 = 5`.

Used helpers: `chiP_two_eq_one`, `chiP_fifty_nine_eq_one` (new),
`chiP_eighty_nine_eq_one` (new), `chiP_one_hundred_nine_eq_one` (new),
`chiP_one_hundred_twenty_seven_eq_one` (new), and decidability of
non-primality for `1, 49, 50, 55, 85, 86, 91, 95, 110, 121, 122`.

**Why W=12 uses a different row set than W=8.**
W=8 uses leading rows `{0, 1, 2, 3, 4}` (after row permutation
`(2, 0, 1, 3, 4)`). W=12 uses non-leading rows `{0, 4, 7, 9, 10}` —
the first leading-row permutations all fail because the first 5 rows
have multiple identical support patterns at the chosen 5 columns (e.g.,
rows 1 and 3 both have support `(0, 1, 0, 1, 0)` at cols
`(k+1) = (2, 1, 7, 11, 5)`). Non-leading rows `{4, 7, 9, 10}` each
contribute a distinct prime (`59, 89, 109, 127`) at a distinct column,
breaking the linear-dependence and admitting a clean triangulation. This
mirrors the W=6 row-`{0, 1, 5}` trick from S122 but with **four** non-
leading rows instead of one.
-/

/--
`chiP 59 = 1` since `59` is prime.
-/
theorem chiP_fifty_nine_eq_one : chiP 59 = 1 := by
  have h_prime_59 : Nat.Prime 59 := by decide
  simp [chiP, h_prime_59]

/--
`chiP 89 = 1` since `89` is prime.
-/
theorem chiP_eighty_nine_eq_one : chiP 89 = 1 := by
  have h_prime_89 : Nat.Prime 89 := by decide
  simp [chiP, h_prime_89]

/--
`chiP 109 = 1` since `109` is prime.
-/
theorem chiP_one_hundred_nine_eq_one : chiP 109 = 1 := by
  have h_prime_109 : Nat.Prime 109 := by decide
  simp [chiP, h_prime_109]

/--
`chiP 127 = 1` since `127` is prime.
-/
theorem chiP_one_hundred_twenty_seven_eq_one : chiP 127 = 1 := by
  have h_prime_127 : Nat.Prime 127 := by decide
  simp [chiP, h_prime_127]

/--
**Corner-case prime exhibit (W = 12, d = j + 1).**

Specialisation of `exists_invertible_submatrix` to the corner case
`W = 12`, `d = j + 1`, where `R = min(12^j, 5) = 5` for every `j ≥ 1`.
Closed unconditionally using `chiP_two_eq_one`, `chiP_fifty_nine_eq_one`,
`chiP_eighty_nine_eq_one`, `chiP_one_hundred_nine_eq_one`, and
`chiP_one_hundred_twenty_seven_eq_one`.

The construction uses `ρ : Fin 5 → Fin (12^j)` mapping
`(0, 1, 2, 3, 4) ↦ (0, 9, 10, 4, 7)` and
`σ : Fin 5 → Fin (12^((j+1)-j))` mapping `(0, 1, 2, 3, 4) ↦ (1, 0, 6, 10, 4)`.
The `5 × 5` submatrix is upper triangular with `1` on the diagonal,
hence `IsUnit` over `ℚ`.
-/
theorem exists_invertible_submatrix_W_eq_12_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    ∃ (ρ : Fin 5 → Fin (12 ^ j))
      (σ : Fin 5 → Fin (12 ^ ((j + 1) - j))),
      IsUnit ((unfolding 12 (j + 1) j).submatrix ρ σ) := by
  -- The exponent simplifies: `(j + 1) - j = 1`, hence `12 ^ ((j + 1) - j) = 12`.
  have h_sub : (j + 1 : ℕ) - j = 1 := by omega
  have h_pow_j_pos : 0 < 12 ^ j := Nat.pow_pos (by norm_num)
  have h_four_lt_pow_j : 4 < 12 ^ j := by
    calc (4 : ℕ) < 12 := by norm_num
      _ = 12 ^ 1 := by norm_num
      _ ≤ 12 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_seven_lt_pow_j : 7 < 12 ^ j := by
    calc (7 : ℕ) < 12 := by norm_num
      _ = 12 ^ 1 := by norm_num
      _ ≤ 12 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_nine_lt_pow_j : 9 < 12 ^ j := by
    calc (9 : ℕ) < 12 := by norm_num
      _ = 12 ^ 1 := by norm_num
      _ ≤ 12 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_ten_lt_pow_j : 10 < 12 ^ j := by
    calc (10 : ℕ) < 12 := by norm_num
      _ = 12 ^ 1 := by norm_num
      _ ≤ 12 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_zero_lt_pow_dj : (0 : ℕ) < 12 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_one_lt_pow_dj : (1 : ℕ) < 12 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_four_lt_pow_dj : (4 : ℕ) < 12 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_six_lt_pow_dj : (6 : ℕ) < 12 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_ten_lt_pow_dj : (10 : ℕ) < 12 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  -- ρ permutes original rows in the order `(0, 9, 10, 4, 7)`.
  -- σ permutes original columns in the order `(1, 0, 6, 10, 4)`.
  refine ⟨fun i => if i.val = 0 then ⟨0, h_pow_j_pos⟩
                   else if i.val = 1 then ⟨9, h_nine_lt_pow_j⟩
                   else if i.val = 2 then ⟨10, h_ten_lt_pow_j⟩
                   else if i.val = 3 then ⟨4, h_four_lt_pow_j⟩
                   else ⟨7, h_seven_lt_pow_j⟩,
          fun i => if i.val = 0 then ⟨1, h_one_lt_pow_dj⟩
                   else if i.val = 1 then ⟨0, h_zero_lt_pow_dj⟩
                   else if i.val = 2 then ⟨6, h_six_lt_pow_dj⟩
                   else if i.val = 3 then ⟨10, h_ten_lt_pow_dj⟩
                   else ⟨4, h_four_lt_pow_dj⟩, ?_⟩
  -- Reduce `IsUnit` to `IsUnit (det)` and compute the det as the
  -- product of the diagonal via `Matrix.det_of_upperTriangular`.
  rw [Matrix.isUnit_iff_isUnit_det]
  -- Non-primality lemmas for the ten below-diagonal zeros.
  have h_not_prime_49  : ¬ Nat.Prime 49  := by decide
  have h_not_prime_50  : ¬ Nat.Prime 50  := by decide
  have h_not_prime_55  : ¬ Nat.Prime 55  := by decide
  have h_not_prime_85  : ¬ Nat.Prime 85  := by decide
  have h_not_prime_86  : ¬ Nat.Prime 86  := by decide
  have h_not_prime_91  : ¬ Nat.Prime 91  := by decide
  have h_not_prime_95  : ¬ Nat.Prime 95  := by decide
  have h_not_prime_110 : ¬ Nat.Prime 110 := by decide
  have h_not_prime_121 : ¬ Nat.Prime 121 := by decide
  have h_not_prime_122 : ¬ Nat.Prime 122 := by decide
  -- Compute the 5 diagonal entries (all = 1).
  -- (0,0): unfolding(0, 1) = chiP 2 = 1
  have hD00 : unfolding 12 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨1, h_one_lt_pow_dj⟩ = 1 := by
    change chiP (0 * 12 ^ ((j + 1) - j) + 1 + 1) = 1
    simp [chiP, Nat.prime_two]
  -- (1,1): unfolding(9, 0) = chiP 109 = 1
  have hD11 : unfolding 12 (j + 1) j ⟨9, h_nine_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 1 := by
    change chiP (9 * 12 ^ ((j + 1) - j) + 0 + 1) = 1
    rw [h_sub]
    have h_eq : (9 * 12 ^ 1 + 0 + 1 : ℕ) = 109 := by norm_num
    rw [h_eq]
    exact chiP_one_hundred_nine_eq_one
  -- (2,2): unfolding(10, 6) = chiP 127 = 1
  have hD22 : unfolding 12 (j + 1) j ⟨10, h_ten_lt_pow_j⟩ ⟨6, h_six_lt_pow_dj⟩ = 1 := by
    change chiP (10 * 12 ^ ((j + 1) - j) + 6 + 1) = 1
    rw [h_sub]
    have h_eq : (10 * 12 ^ 1 + 6 + 1 : ℕ) = 127 := by norm_num
    rw [h_eq]
    exact chiP_one_hundred_twenty_seven_eq_one
  -- (3,3): unfolding(4, 10) = chiP 59 = 1
  have hD33 : unfolding 12 (j + 1) j ⟨4, h_four_lt_pow_j⟩ ⟨10, h_ten_lt_pow_dj⟩ = 1 := by
    change chiP (4 * 12 ^ ((j + 1) - j) + 10 + 1) = 1
    rw [h_sub]
    have h_eq : (4 * 12 ^ 1 + 10 + 1 : ℕ) = 59 := by norm_num
    rw [h_eq]
    exact chiP_fifty_nine_eq_one
  -- (4,4): unfolding(7, 4) = chiP 89 = 1
  have hD44 : unfolding 12 (j + 1) j ⟨7, h_seven_lt_pow_j⟩ ⟨4, h_four_lt_pow_dj⟩ = 1 := by
    change chiP (7 * 12 ^ ((j + 1) - j) + 4 + 1) = 1
    rw [h_sub]
    have h_eq : (7 * 12 ^ 1 + 4 + 1 : ℕ) = 89 := by norm_num
    rw [h_eq]
    exact chiP_eighty_nine_eq_one
  -- Compute the 10 below-diagonal entries (all = 0).
  -- (1,0): unfolding(9, 1) = chiP 110 = 0
  have hL10 : unfolding 12 (j + 1) j ⟨9, h_nine_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (9 * 12 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (9 * 12 ^ 1 + 1 + 1 : ℕ) = 110 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_110]
  -- (2,0): unfolding(10, 1) = chiP 122 = 0
  have hL20 : unfolding 12 (j + 1) j ⟨10, h_ten_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (10 * 12 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (10 * 12 ^ 1 + 1 + 1 : ℕ) = 122 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_122]
  -- (2,1): unfolding(10, 0) = chiP 121 = 0
  have hL21 : unfolding 12 (j + 1) j ⟨10, h_ten_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (10 * 12 ^ ((j + 1) - j) + 0 + 1) = 0
    rw [h_sub]
    have h_eq : (10 * 12 ^ 1 + 0 + 1 : ℕ) = 121 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_121]
  -- (3,0): unfolding(4, 1) = chiP 50 = 0
  have hL30 : unfolding 12 (j + 1) j ⟨4, h_four_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (4 * 12 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (4 * 12 ^ 1 + 1 + 1 : ℕ) = 50 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_50]
  -- (3,1): unfolding(4, 0) = chiP 49 = 0
  have hL31 : unfolding 12 (j + 1) j ⟨4, h_four_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (4 * 12 ^ ((j + 1) - j) + 0 + 1) = 0
    rw [h_sub]
    have h_eq : (4 * 12 ^ 1 + 0 + 1 : ℕ) = 49 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_49]
  -- (3,2): unfolding(4, 6) = chiP 55 = 0
  have hL32 : unfolding 12 (j + 1) j ⟨4, h_four_lt_pow_j⟩ ⟨6, h_six_lt_pow_dj⟩ = 0 := by
    change chiP (4 * 12 ^ ((j + 1) - j) + 6 + 1) = 0
    rw [h_sub]
    have h_eq : (4 * 12 ^ 1 + 6 + 1 : ℕ) = 55 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_55]
  -- (4,0): unfolding(7, 1) = chiP 86 = 0
  have hL40 : unfolding 12 (j + 1) j ⟨7, h_seven_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (7 * 12 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (7 * 12 ^ 1 + 1 + 1 : ℕ) = 86 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_86]
  -- (4,1): unfolding(7, 0) = chiP 85 = 0
  have hL41 : unfolding 12 (j + 1) j ⟨7, h_seven_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (7 * 12 ^ ((j + 1) - j) + 0 + 1) = 0
    rw [h_sub]
    have h_eq : (7 * 12 ^ 1 + 0 + 1 : ℕ) = 85 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_85]
  -- (4,2): unfolding(7, 6) = chiP 91 = 0
  have hL42 : unfolding 12 (j + 1) j ⟨7, h_seven_lt_pow_j⟩ ⟨6, h_six_lt_pow_dj⟩ = 0 := by
    change chiP (7 * 12 ^ ((j + 1) - j) + 6 + 1) = 0
    rw [h_sub]
    have h_eq : (7 * 12 ^ 1 + 6 + 1 : ℕ) = 91 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_91]
  -- (4,3): unfolding(7, 10) = chiP 95 = 0
  have hL43 : unfolding 12 (j + 1) j ⟨7, h_seven_lt_pow_j⟩ ⟨10, h_ten_lt_pow_dj⟩ = 0 := by
    change chiP (7 * 12 ^ ((j + 1) - j) + 10 + 1) = 0
    rw [h_sub]
    have h_eq : (7 * 12 ^ 1 + 10 + 1 : ℕ) = 95 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_95]
  -- Helpers for `Fin 5` `.val` and the `if-then-else` ρ, σ.
  have h0 : (0 : Fin 5).val = 0 := rfl
  have h1 : (1 : Fin 5).val = 1 := rfl
  have h2 : (2 : Fin 5).val = 2 := rfl
  have h3 : (3 : Fin 5).val = 3 := rfl
  have h4 : (4 : Fin 5).val = 4 := rfl
  have hne_2_0 : (2 : ℕ) ≠ 0 := by decide
  have hne_2_1 : (2 : ℕ) ≠ 1 := by decide
  have hne_3_0 : (3 : ℕ) ≠ 0 := by decide
  have hne_3_1 : (3 : ℕ) ≠ 1 := by decide
  have hne_3_2 : (3 : ℕ) ≠ 2 := by decide
  have hne_4_0 : (4 : ℕ) ≠ 0 := by decide
  have hne_4_1 : (4 : ℕ) ≠ 1 := by decide
  have hne_4_2 : (4 : ℕ) ≠ 2 := by decide
  have hne_4_3 : (4 : ℕ) ≠ 3 := by decide
  -- Establish the upper-triangular property of the submatrix.
  set Mρ : Fin 5 → Fin (12 ^ j) :=
    fun i => if i.val = 0 then ⟨0, h_pow_j_pos⟩
             else if i.val = 1 then ⟨9, h_nine_lt_pow_j⟩
             else if i.val = 2 then ⟨10, h_ten_lt_pow_j⟩
             else if i.val = 3 then ⟨4, h_four_lt_pow_j⟩
             else ⟨7, h_seven_lt_pow_j⟩ with hMρ_def
  set Mσ : Fin 5 → Fin (12 ^ ((j + 1) - j)) :=
    fun i => if i.val = 0 then ⟨1, h_one_lt_pow_dj⟩
             else if i.val = 1 then ⟨0, h_zero_lt_pow_dj⟩
             else if i.val = 2 then ⟨6, h_six_lt_pow_dj⟩
             else if i.val = 3 then ⟨10, h_ten_lt_pow_dj⟩
             else ⟨4, h_four_lt_pow_dj⟩ with hMσ_def
  set M : Matrix (Fin 5) (Fin 5) ℚ :=
    (unfolding 12 (j + 1) j).submatrix Mρ Mσ with hM_def
  have h_blocktri : M.BlockTriangular id := by
    intro i k h_lt
    -- `id k < id i` reduces to `k.val < i.val`.
    simp only [id_eq, Fin.lt_def] at h_lt
    fin_cases i <;> fin_cases k <;>
      first
        | (exact absurd h_lt (by decide))
        | (simp only [hM_def, Matrix.submatrix_apply, hMρ_def, hMσ_def,
                       if_true, if_false,
                       Nat.one_ne_zero, hne_2_0, hne_2_1, hne_3_0, hne_3_1,
                       hne_3_2, hne_4_0, hne_4_1, hne_4_2, hne_4_3]
           first | exact hL10 | exact hL20 | exact hL21
                 | exact hL30 | exact hL31 | exact hL32
                 | exact hL40 | exact hL41 | exact hL42 | exact hL43)
  rw [Matrix.det_of_upperTriangular h_blocktri]
  -- Now expand the diagonal product over `Fin 5`.
  rw [Fin.prod_univ_five]
  -- Reduce `M k k` for `k ∈ {0, 1, 2, 3, 4}` to the precomputed diagonal entries.
  simp only [hM_def, Matrix.submatrix_apply, hMρ_def, hMσ_def,
             h0, h1, h2, h3, h4, if_true, if_false,
             Nat.one_ne_zero, hne_2_0, hne_2_1, hne_3_0, hne_3_1, hne_3_2,
             hne_4_0, hne_4_1, hne_4_2, hne_4_3]
  rw [hD00, hD11, hD22, hD33, hD44]
  -- Goal: `IsUnit (1 * 1 * 1 * 1 * 1 : ℚ)`.
  norm_num

/--
**Corner-case main theorem (W = 12, d = j + 1).** The unfolding rank is
exactly `5` for every `j ≥ 1`. Mirrors `mps_bond_dim_W_eq_8_d_eq_j_plus_1`
(both citing the general `upper_bound`) and is fully formalised in Lean
(no `sorry`, no `axiom`).

**Upper-bound subtlety:** the matrix has `12^((j+1)-j) = 12` columns, so
`rank_le_width` only gives `rank ≤ 12`. We instead cite the general
`upper_bound` lemma, which evaluates to `φ(12) · 12^0 + 1 = 4 · 1 + 1 = 5`
in this corner.

Together with the `W ∈ {2, 3, 4, 5, 6, 8}` corners, this is the eighth
unconditional instance of `mps_bond_dim` and the sixth instance over a
wheel `W ≥ 3`. The general `mps_bond_dim` still requires
`exists_invertible_submatrix` whose proof is the only remaining `sorry`
in this file.
-/
theorem mps_bond_dim_W_eq_12_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    (unfolding 12 (j + 1) j).rank = 5 := by
  apply Nat.le_antisymm
  · -- `≤ 5`: from `upper_bound`. (`rank_le_width` only gives `≤ 12`.)
    have hW : (2 : ℕ) ≤ 12 := by norm_num
    have hj_hi : j < j + 1 := Nat.lt_succ_self j
    have h := upper_bound 12 (j + 1) j hW hj hj_hi
    have h_eq : Nat.totient 12 * 12 ^ ((j + 1 : ℕ) - j - 1) + 1 = 5 := by
      have h_sub : (j + 1 : ℕ) - j - 1 = 0 := by omega
      rw [h_sub]
      decide
    linarith
  · -- `5 ≤ rank`: from the corner-case prime exhibit.
    obtain ⟨ρ, σ, hUnit⟩ := exists_invertible_submatrix_W_eq_12_d_eq_j_plus_1 j hj
    have h_eq : ((unfolding 12 (j + 1) j).submatrix ρ σ).rank = 5 := by
      have h := Matrix.rank_of_isUnit ((unfolding 12 (j + 1) j).submatrix ρ σ) hUnit
      simpa using h
    calc (5 : ℕ) = ((unfolding 12 (j + 1) j).submatrix ρ σ).rank := h_eq.symm
      _ ≤ (unfolding 12 (j + 1) j).rank := Matrix.rank_submatrix_le _ _ _

/-!
### Corner case `(W = 18, d = j + 1)` (S137)

Closes `exists_invertible_submatrix` for `W = 18, d = j + 1` by exhibiting
an invertible `7 × 7` submatrix of the `18^j × 18` unfolding.

`R = min(18^j, φ(18) · 18^0 + 1) = min(18^j, 7) = 7` for every `j ≥ 1`,
since `18^j ≥ 18 > 7`.

The construction permutes rows and columns to upper-triangularise the
slab. The row permutation `ρ : Fin 7 → Fin (18^j)` maps
`(0, 1, 2, 3, 4, 5, 6) ↦ (0, 2, 9, 1, 11, 6, 16)` and the column
permutation `σ : Fin 7 → Fin 18` maps
`(0, 1, 2, 3, 4, 5, 6) ↦ (1, 6, 16, 10, 12, 0, 4)`.

The resulting `7 × 7` submatrix is upper-triangular with `1`s on the
diagonal:

```
   ⎡ chiP   2, chiP   7, chiP  17, chiP  11, chiP  13, chiP   1, chiP   5 ⎤   ⎡ 1, 1, 1, 1, 1, 0, 1 ⎤
   ⎢ chiP  38, chiP  43, chiP  53, chiP  47, chiP  49, chiP  37, chiP  41 ⎥   ⎢ 0, 1, 1, 1, 0, 1, 1 ⎥
   ⎢ chiP 164, chiP 169, chiP 179, chiP 173, chiP 175, chiP 163, chiP 167 ⎥   ⎢ 0, 0, 1, 1, 0, 1, 1 ⎥
   ⎢ chiP  20, chiP  25, chiP  35, chiP  29, chiP  31, chiP  19, chiP  23 ⎥ = ⎢ 0, 0, 0, 1, 1, 1, 1 ⎥
   ⎢ chiP 200, chiP 205, chiP 215, chiP 209, chiP 211, chiP 199, chiP 203 ⎥   ⎢ 0, 0, 0, 0, 1, 1, 0 ⎥
   ⎢ chiP 110, chiP 115, chiP 125, chiP 119, chiP 121, chiP 109, chiP 113 ⎥   ⎢ 0, 0, 0, 0, 0, 1, 1 ⎥
   ⎣ chiP 290, chiP 295, chiP 305, chiP 299, chiP 301, chiP 289, chiP 293 ⎦   ⎣ 0, 0, 0, 0, 0, 0, 1 ⎦
```

with diagonal primes `{2, 43, 179, 29, 211, 109, 293}` and
below-diagonal composites
`{20, 25, 35, 38, 110, 115, 119, 121, 125, 164, 169, 200, 205, 209, 215,
 289, 290, 295, 299, 301, 305}`.

Determinant via `Matrix.det_of_upperTriangular` and `Fin.prod_univ_seven`
(which mathlib provides). The unit witness comes from `det = 1`.

Used helpers: `chiP_two_eq_one`, `chiP_one_hundred_nine_eq_one` (existing
from S129) and `chiP_twenty_nine_eq_one`, `chiP_forty_three_eq_one`,
`chiP_one_hundred_seventy_nine_eq_one`, `chiP_two_hundred_eleven_eq_one`,
`chiP_two_hundred_ninety_three_eq_one` (new at S137), and decidability of
non-primality for the 21 below-diagonal composites.

**Why W=18 closes a new corner over the leading-row regime.**
At W=14 (also `R = 7`), the 14 rows of the j=1 slab admit no
upper-triangulation: rows `2` and `5` of the W=14 slab have identical
support pattern at the 7 chosen columns, and no permutation (with all
rho values < 14) achieves an upper-triangular form. W=18 is the
smallest `W` with `R = 7` and `j = 1` row-count `W^j = W ≥ 17` (since
the search at W=18 found valid triangulations using `ρ ∈ {0,1,2,6,9,11,16}`,
all `< 18`). This makes W=18 the first `R = 7` corner closure;
W ∈ {7, 14} remain structurally obstructed at the leading-row stage and
W = 9 needs `Matrix.det_of_blockTriangular`.

**Ninth unconditional `mps_bond_dim` instance; seventh instance over a
wheel `W ≥ 3`; fourth instance using `det_of_upperTriangular`; first
instance with `R = 7`.**
-/

/--
`chiP 29 = 1` since `29` is prime.
-/
theorem chiP_twenty_nine_eq_one : chiP 29 = 1 := by
  have h_prime_29 : Nat.Prime 29 := by norm_num
  simp [chiP, h_prime_29]

/--
`chiP 43 = 1` since `43` is prime.
-/
theorem chiP_forty_three_eq_one : chiP 43 = 1 := by
  have h_prime_43 : Nat.Prime 43 := by norm_num
  simp [chiP, h_prime_43]

/--
`chiP 179 = 1` since `179` is prime.
-/
theorem chiP_one_hundred_seventy_nine_eq_one : chiP 179 = 1 := by
  have h_prime_179 : Nat.Prime 179 := by norm_num
  simp [chiP, h_prime_179]

/--
`chiP 211 = 1` since `211` is prime.
-/
theorem chiP_two_hundred_eleven_eq_one : chiP 211 = 1 := by
  have h_prime_211 : Nat.Prime 211 := by norm_num
  simp [chiP, h_prime_211]

/--
`chiP 293 = 1` since `293` is prime.
-/
theorem chiP_two_hundred_ninety_three_eq_one : chiP 293 = 1 := by
  have h_prime_293 : Nat.Prime 293 := by norm_num
  simp [chiP, h_prime_293]

/--
**Corner-case prime exhibit (W = 18, d = j + 1).**

Specialisation of `exists_invertible_submatrix` to the corner case
`W = 18`, `d = j + 1`, where `R = min(18^j, 7) = 7` for every `j ≥ 1`.
Closed unconditionally using `chiP_two_eq_one`, `chiP_twenty_nine_eq_one`,
`chiP_forty_three_eq_one`, `chiP_one_hundred_nine_eq_one`,
`chiP_one_hundred_seventy_nine_eq_one`, `chiP_two_hundred_eleven_eq_one`,
and `chiP_two_hundred_ninety_three_eq_one`.

The construction uses `ρ : Fin 7 → Fin (18^j)` mapping
`(0, 1, 2, 3, 4, 5, 6) ↦ (0, 2, 9, 1, 11, 6, 16)` and
`σ : Fin 7 → Fin (18^((j+1)-j))` mapping
`(0, 1, 2, 3, 4, 5, 6) ↦ (1, 6, 16, 10, 12, 0, 4)`.
The `7 × 7` submatrix is upper triangular with `1` on the diagonal,
hence `IsUnit` over `ℚ`.
-/
theorem exists_invertible_submatrix_W_eq_18_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    ∃ (ρ : Fin 7 → Fin (18 ^ j))
      (σ : Fin 7 → Fin (18 ^ ((j + 1) - j))),
      IsUnit ((unfolding 18 (j + 1) j).submatrix ρ σ) := by
  -- The exponent simplifies: `(j + 1) - j = 1`, hence `18 ^ ((j + 1) - j) = 18`.
  have h_sub : (j + 1 : ℕ) - j = 1 := by omega
  have h_pow_j_pos : 0 < 18 ^ j := Nat.pow_pos (by norm_num)
  have h_one_lt_pow_j : 1 < 18 ^ j := by
    calc (1 : ℕ) < 18 := by norm_num
      _ = 18 ^ 1 := by norm_num
      _ ≤ 18 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_two_lt_pow_j : 2 < 18 ^ j := by
    calc (2 : ℕ) < 18 := by norm_num
      _ = 18 ^ 1 := by norm_num
      _ ≤ 18 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_six_lt_pow_j : 6 < 18 ^ j := by
    calc (6 : ℕ) < 18 := by norm_num
      _ = 18 ^ 1 := by norm_num
      _ ≤ 18 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_nine_lt_pow_j : 9 < 18 ^ j := by
    calc (9 : ℕ) < 18 := by norm_num
      _ = 18 ^ 1 := by norm_num
      _ ≤ 18 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_eleven_lt_pow_j : 11 < 18 ^ j := by
    calc (11 : ℕ) < 18 := by norm_num
      _ = 18 ^ 1 := by norm_num
      _ ≤ 18 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_sixteen_lt_pow_j : 16 < 18 ^ j := by
    calc (16 : ℕ) < 18 := by norm_num
      _ = 18 ^ 1 := by norm_num
      _ ≤ 18 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_zero_lt_pow_dj : (0 : ℕ) < 18 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_one_lt_pow_dj : (1 : ℕ) < 18 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_four_lt_pow_dj : (4 : ℕ) < 18 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_six_lt_pow_dj : (6 : ℕ) < 18 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_ten_lt_pow_dj : (10 : ℕ) < 18 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_twelve_lt_pow_dj : (12 : ℕ) < 18 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_sixteen_lt_pow_dj : (16 : ℕ) < 18 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  -- ρ permutes original rows in the order `(0, 2, 9, 1, 11, 6, 16)`.
  -- σ permutes original columns in the order `(1, 6, 16, 10, 12, 0, 4)`.
  refine ⟨fun i => if i.val = 0 then ⟨0, h_pow_j_pos⟩
                   else if i.val = 1 then ⟨2, h_two_lt_pow_j⟩
                   else if i.val = 2 then ⟨9, h_nine_lt_pow_j⟩
                   else if i.val = 3 then ⟨1, h_one_lt_pow_j⟩
                   else if i.val = 4 then ⟨11, h_eleven_lt_pow_j⟩
                   else if i.val = 5 then ⟨6, h_six_lt_pow_j⟩
                   else ⟨16, h_sixteen_lt_pow_j⟩,
          fun i => if i.val = 0 then ⟨1, h_one_lt_pow_dj⟩
                   else if i.val = 1 then ⟨6, h_six_lt_pow_dj⟩
                   else if i.val = 2 then ⟨16, h_sixteen_lt_pow_dj⟩
                   else if i.val = 3 then ⟨10, h_ten_lt_pow_dj⟩
                   else if i.val = 4 then ⟨12, h_twelve_lt_pow_dj⟩
                   else if i.val = 5 then ⟨0, h_zero_lt_pow_dj⟩
                   else ⟨4, h_four_lt_pow_dj⟩, ?_⟩
  -- Reduce `IsUnit` to `IsUnit (det)` and compute the det as the
  -- product of the diagonal via `Matrix.det_of_upperTriangular`.
  rw [Matrix.isUnit_iff_isUnit_det]
  -- Non-primality lemmas for the 21 below-diagonal zeros.
  have h_not_prime_20  : ¬ Nat.Prime 20  := by decide
  have h_not_prime_25  : ¬ Nat.Prime 25  := by decide
  have h_not_prime_35  : ¬ Nat.Prime 35  := by decide
  have h_not_prime_38  : ¬ Nat.Prime 38  := by decide
  have h_not_prime_110 : ¬ Nat.Prime 110 := by decide
  have h_not_prime_115 : ¬ Nat.Prime 115 := by decide
  have h_not_prime_119 : ¬ Nat.Prime 119 := by decide
  have h_not_prime_121 : ¬ Nat.Prime 121 := by decide
  have h_not_prime_125 : ¬ Nat.Prime 125 := by decide
  have h_not_prime_164 : ¬ Nat.Prime 164 := by norm_num
  have h_not_prime_169 : ¬ Nat.Prime 169 := by norm_num
  have h_not_prime_200 : ¬ Nat.Prime 200 := by norm_num
  have h_not_prime_205 : ¬ Nat.Prime 205 := by norm_num
  have h_not_prime_209 : ¬ Nat.Prime 209 := by norm_num
  have h_not_prime_215 : ¬ Nat.Prime 215 := by norm_num
  have h_not_prime_289 : ¬ Nat.Prime 289 := by norm_num
  have h_not_prime_290 : ¬ Nat.Prime 290 := by norm_num
  have h_not_prime_295 : ¬ Nat.Prime 295 := by norm_num
  have h_not_prime_299 : ¬ Nat.Prime 299 := by norm_num
  have h_not_prime_301 : ¬ Nat.Prime 301 := by norm_num
  have h_not_prime_305 : ¬ Nat.Prime 305 := by norm_num
  -- Compute the 7 diagonal entries (all = 1).
  -- (0,0): unfolding(0, 1) = chiP 2 = 1
  have hD00 : unfolding 18 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨1, h_one_lt_pow_dj⟩ = 1 := by
    change chiP (0 * 18 ^ ((j + 1) - j) + 1 + 1) = 1
    simp [chiP, Nat.prime_two]
  -- (1,1): unfolding(2, 6) = chiP 43 = 1
  have hD11 : unfolding 18 (j + 1) j ⟨2, h_two_lt_pow_j⟩ ⟨6, h_six_lt_pow_dj⟩ = 1 := by
    change chiP (2 * 18 ^ ((j + 1) - j) + 6 + 1) = 1
    rw [h_sub]
    have h_eq : (2 * 18 ^ 1 + 6 + 1 : ℕ) = 43 := by norm_num
    rw [h_eq]
    exact chiP_forty_three_eq_one
  -- (2,2): unfolding(9, 16) = chiP 179 = 1
  have hD22 : unfolding 18 (j + 1) j ⟨9, h_nine_lt_pow_j⟩ ⟨16, h_sixteen_lt_pow_dj⟩ = 1 := by
    change chiP (9 * 18 ^ ((j + 1) - j) + 16 + 1) = 1
    rw [h_sub]
    have h_eq : (9 * 18 ^ 1 + 16 + 1 : ℕ) = 179 := by norm_num
    rw [h_eq]
    exact chiP_one_hundred_seventy_nine_eq_one
  -- (3,3): unfolding(1, 10) = chiP 29 = 1
  have hD33 : unfolding 18 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨10, h_ten_lt_pow_dj⟩ = 1 := by
    change chiP (1 * 18 ^ ((j + 1) - j) + 10 + 1) = 1
    rw [h_sub]
    have h_eq : (1 * 18 ^ 1 + 10 + 1 : ℕ) = 29 := by norm_num
    rw [h_eq]
    exact chiP_twenty_nine_eq_one
  -- (4,4): unfolding(11, 12) = chiP 211 = 1
  have hD44 : unfolding 18 (j + 1) j ⟨11, h_eleven_lt_pow_j⟩ ⟨12, h_twelve_lt_pow_dj⟩ = 1 := by
    change chiP (11 * 18 ^ ((j + 1) - j) + 12 + 1) = 1
    rw [h_sub]
    have h_eq : (11 * 18 ^ 1 + 12 + 1 : ℕ) = 211 := by norm_num
    rw [h_eq]
    exact chiP_two_hundred_eleven_eq_one
  -- (5,5): unfolding(6, 0) = chiP 109 = 1
  have hD55 : unfolding 18 (j + 1) j ⟨6, h_six_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 1 := by
    change chiP (6 * 18 ^ ((j + 1) - j) + 0 + 1) = 1
    rw [h_sub]
    have h_eq : (6 * 18 ^ 1 + 0 + 1 : ℕ) = 109 := by norm_num
    rw [h_eq]
    exact chiP_one_hundred_nine_eq_one
  -- (6,6): unfolding(16, 4) = chiP 293 = 1
  have hD66 : unfolding 18 (j + 1) j ⟨16, h_sixteen_lt_pow_j⟩ ⟨4, h_four_lt_pow_dj⟩ = 1 := by
    change chiP (16 * 18 ^ ((j + 1) - j) + 4 + 1) = 1
    rw [h_sub]
    have h_eq : (16 * 18 ^ 1 + 4 + 1 : ℕ) = 293 := by norm_num
    rw [h_eq]
    exact chiP_two_hundred_ninety_three_eq_one
  -- Compute the 21 below-diagonal entries (all = 0).
  -- (1,0): unfolding(2, 1) = chiP 38 = 0
  have hL10 : unfolding 18 (j + 1) j ⟨2, h_two_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (2 * 18 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (2 * 18 ^ 1 + 1 + 1 : ℕ) = 38 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_38]
  -- (2,0): unfolding(9, 1) = chiP 164 = 0
  have hL20 : unfolding 18 (j + 1) j ⟨9, h_nine_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (9 * 18 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (9 * 18 ^ 1 + 1 + 1 : ℕ) = 164 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_164]
  -- (2,1): unfolding(9, 6) = chiP 169 = 0
  have hL21 : unfolding 18 (j + 1) j ⟨9, h_nine_lt_pow_j⟩ ⟨6, h_six_lt_pow_dj⟩ = 0 := by
    change chiP (9 * 18 ^ ((j + 1) - j) + 6 + 1) = 0
    rw [h_sub]
    have h_eq : (9 * 18 ^ 1 + 6 + 1 : ℕ) = 169 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_169]
  -- (3,0): unfolding(1, 1) = chiP 20 = 0
  have hL30 : unfolding 18 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (1 * 18 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (1 * 18 ^ 1 + 1 + 1 : ℕ) = 20 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_20]
  -- (3,1): unfolding(1, 6) = chiP 25 = 0
  have hL31 : unfolding 18 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨6, h_six_lt_pow_dj⟩ = 0 := by
    change chiP (1 * 18 ^ ((j + 1) - j) + 6 + 1) = 0
    rw [h_sub]
    have h_eq : (1 * 18 ^ 1 + 6 + 1 : ℕ) = 25 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_25]
  -- (3,2): unfolding(1, 16) = chiP 35 = 0
  have hL32 : unfolding 18 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨16, h_sixteen_lt_pow_dj⟩ = 0 := by
    change chiP (1 * 18 ^ ((j + 1) - j) + 16 + 1) = 0
    rw [h_sub]
    have h_eq : (1 * 18 ^ 1 + 16 + 1 : ℕ) = 35 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_35]
  -- (4,0): unfolding(11, 1) = chiP 200 = 0
  have hL40 : unfolding 18 (j + 1) j ⟨11, h_eleven_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (11 * 18 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (11 * 18 ^ 1 + 1 + 1 : ℕ) = 200 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_200]
  -- (4,1): unfolding(11, 6) = chiP 205 = 0
  have hL41 : unfolding 18 (j + 1) j ⟨11, h_eleven_lt_pow_j⟩ ⟨6, h_six_lt_pow_dj⟩ = 0 := by
    change chiP (11 * 18 ^ ((j + 1) - j) + 6 + 1) = 0
    rw [h_sub]
    have h_eq : (11 * 18 ^ 1 + 6 + 1 : ℕ) = 205 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_205]
  -- (4,2): unfolding(11, 16) = chiP 215 = 0
  have hL42 : unfolding 18 (j + 1) j ⟨11, h_eleven_lt_pow_j⟩ ⟨16, h_sixteen_lt_pow_dj⟩ = 0 := by
    change chiP (11 * 18 ^ ((j + 1) - j) + 16 + 1) = 0
    rw [h_sub]
    have h_eq : (11 * 18 ^ 1 + 16 + 1 : ℕ) = 215 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_215]
  -- (4,3): unfolding(11, 10) = chiP 209 = 0
  have hL43 : unfolding 18 (j + 1) j ⟨11, h_eleven_lt_pow_j⟩ ⟨10, h_ten_lt_pow_dj⟩ = 0 := by
    change chiP (11 * 18 ^ ((j + 1) - j) + 10 + 1) = 0
    rw [h_sub]
    have h_eq : (11 * 18 ^ 1 + 10 + 1 : ℕ) = 209 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_209]
  -- (5,0): unfolding(6, 1) = chiP 110 = 0
  have hL50 : unfolding 18 (j + 1) j ⟨6, h_six_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (6 * 18 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (6 * 18 ^ 1 + 1 + 1 : ℕ) = 110 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_110]
  -- (5,1): unfolding(6, 6) = chiP 115 = 0
  have hL51 : unfolding 18 (j + 1) j ⟨6, h_six_lt_pow_j⟩ ⟨6, h_six_lt_pow_dj⟩ = 0 := by
    change chiP (6 * 18 ^ ((j + 1) - j) + 6 + 1) = 0
    rw [h_sub]
    have h_eq : (6 * 18 ^ 1 + 6 + 1 : ℕ) = 115 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_115]
  -- (5,2): unfolding(6, 16) = chiP 125 = 0
  have hL52 : unfolding 18 (j + 1) j ⟨6, h_six_lt_pow_j⟩ ⟨16, h_sixteen_lt_pow_dj⟩ = 0 := by
    change chiP (6 * 18 ^ ((j + 1) - j) + 16 + 1) = 0
    rw [h_sub]
    have h_eq : (6 * 18 ^ 1 + 16 + 1 : ℕ) = 125 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_125]
  -- (5,3): unfolding(6, 10) = chiP 119 = 0
  have hL53 : unfolding 18 (j + 1) j ⟨6, h_six_lt_pow_j⟩ ⟨10, h_ten_lt_pow_dj⟩ = 0 := by
    change chiP (6 * 18 ^ ((j + 1) - j) + 10 + 1) = 0
    rw [h_sub]
    have h_eq : (6 * 18 ^ 1 + 10 + 1 : ℕ) = 119 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_119]
  -- (5,4): unfolding(6, 12) = chiP 121 = 0
  have hL54 : unfolding 18 (j + 1) j ⟨6, h_six_lt_pow_j⟩ ⟨12, h_twelve_lt_pow_dj⟩ = 0 := by
    change chiP (6 * 18 ^ ((j + 1) - j) + 12 + 1) = 0
    rw [h_sub]
    have h_eq : (6 * 18 ^ 1 + 12 + 1 : ℕ) = 121 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_121]
  -- (6,0): unfolding(16, 1) = chiP 290 = 0
  have hL60 : unfolding 18 (j + 1) j ⟨16, h_sixteen_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (16 * 18 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (16 * 18 ^ 1 + 1 + 1 : ℕ) = 290 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_290]
  -- (6,1): unfolding(16, 6) = chiP 295 = 0
  have hL61 : unfolding 18 (j + 1) j ⟨16, h_sixteen_lt_pow_j⟩ ⟨6, h_six_lt_pow_dj⟩ = 0 := by
    change chiP (16 * 18 ^ ((j + 1) - j) + 6 + 1) = 0
    rw [h_sub]
    have h_eq : (16 * 18 ^ 1 + 6 + 1 : ℕ) = 295 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_295]
  -- (6,2): unfolding(16, 16) = chiP 305 = 0
  have hL62 : unfolding 18 (j + 1) j ⟨16, h_sixteen_lt_pow_j⟩ ⟨16, h_sixteen_lt_pow_dj⟩ = 0 := by
    change chiP (16 * 18 ^ ((j + 1) - j) + 16 + 1) = 0
    rw [h_sub]
    have h_eq : (16 * 18 ^ 1 + 16 + 1 : ℕ) = 305 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_305]
  -- (6,3): unfolding(16, 10) = chiP 299 = 0
  have hL63 : unfolding 18 (j + 1) j ⟨16, h_sixteen_lt_pow_j⟩ ⟨10, h_ten_lt_pow_dj⟩ = 0 := by
    change chiP (16 * 18 ^ ((j + 1) - j) + 10 + 1) = 0
    rw [h_sub]
    have h_eq : (16 * 18 ^ 1 + 10 + 1 : ℕ) = 299 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_299]
  -- (6,4): unfolding(16, 12) = chiP 301 = 0
  have hL64 : unfolding 18 (j + 1) j ⟨16, h_sixteen_lt_pow_j⟩ ⟨12, h_twelve_lt_pow_dj⟩ = 0 := by
    change chiP (16 * 18 ^ ((j + 1) - j) + 12 + 1) = 0
    rw [h_sub]
    have h_eq : (16 * 18 ^ 1 + 12 + 1 : ℕ) = 301 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_301]
  -- (6,5): unfolding(16, 0) = chiP 289 = 0
  have hL65 : unfolding 18 (j + 1) j ⟨16, h_sixteen_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (16 * 18 ^ ((j + 1) - j) + 0 + 1) = 0
    rw [h_sub]
    have h_eq : (16 * 18 ^ 1 + 0 + 1 : ℕ) = 289 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_289]
  -- Helpers for `Fin 7` `.val` and the `if-then-else` ρ, σ.
  have h0 : (0 : Fin 7).val = 0 := rfl
  have h1 : (1 : Fin 7).val = 1 := rfl
  have h2 : (2 : Fin 7).val = 2 := rfl
  have h3 : (3 : Fin 7).val = 3 := rfl
  have h4 : (4 : Fin 7).val = 4 := rfl
  have h5 : (5 : Fin 7).val = 5 := rfl
  have h6 : (6 : Fin 7).val = 6 := rfl
  have hne_2_0 : (2 : ℕ) ≠ 0 := by decide
  have hne_2_1 : (2 : ℕ) ≠ 1 := by decide
  have hne_3_0 : (3 : ℕ) ≠ 0 := by decide
  have hne_3_1 : (3 : ℕ) ≠ 1 := by decide
  have hne_3_2 : (3 : ℕ) ≠ 2 := by decide
  have hne_4_0 : (4 : ℕ) ≠ 0 := by decide
  have hne_4_1 : (4 : ℕ) ≠ 1 := by decide
  have hne_4_2 : (4 : ℕ) ≠ 2 := by decide
  have hne_4_3 : (4 : ℕ) ≠ 3 := by decide
  have hne_5_0 : (5 : ℕ) ≠ 0 := by decide
  have hne_5_1 : (5 : ℕ) ≠ 1 := by decide
  have hne_5_2 : (5 : ℕ) ≠ 2 := by decide
  have hne_5_3 : (5 : ℕ) ≠ 3 := by decide
  have hne_5_4 : (5 : ℕ) ≠ 4 := by decide
  have hne_6_0 : (6 : ℕ) ≠ 0 := by decide
  have hne_6_1 : (6 : ℕ) ≠ 1 := by decide
  have hne_6_2 : (6 : ℕ) ≠ 2 := by decide
  have hne_6_3 : (6 : ℕ) ≠ 3 := by decide
  have hne_6_4 : (6 : ℕ) ≠ 4 := by decide
  have hne_6_5 : (6 : ℕ) ≠ 5 := by decide
  -- Establish the upper-triangular property of the submatrix.
  set Mρ : Fin 7 → Fin (18 ^ j) :=
    fun i => if i.val = 0 then ⟨0, h_pow_j_pos⟩
             else if i.val = 1 then ⟨2, h_two_lt_pow_j⟩
             else if i.val = 2 then ⟨9, h_nine_lt_pow_j⟩
             else if i.val = 3 then ⟨1, h_one_lt_pow_j⟩
             else if i.val = 4 then ⟨11, h_eleven_lt_pow_j⟩
             else if i.val = 5 then ⟨6, h_six_lt_pow_j⟩
             else ⟨16, h_sixteen_lt_pow_j⟩ with hMρ_def
  set Mσ : Fin 7 → Fin (18 ^ ((j + 1) - j)) :=
    fun i => if i.val = 0 then ⟨1, h_one_lt_pow_dj⟩
             else if i.val = 1 then ⟨6, h_six_lt_pow_dj⟩
             else if i.val = 2 then ⟨16, h_sixteen_lt_pow_dj⟩
             else if i.val = 3 then ⟨10, h_ten_lt_pow_dj⟩
             else if i.val = 4 then ⟨12, h_twelve_lt_pow_dj⟩
             else if i.val = 5 then ⟨0, h_zero_lt_pow_dj⟩
             else ⟨4, h_four_lt_pow_dj⟩ with hMσ_def
  set M : Matrix (Fin 7) (Fin 7) ℚ :=
    (unfolding 18 (j + 1) j).submatrix Mρ Mσ with hM_def
  have h_blocktri : M.BlockTriangular id := by
    intro i k h_lt
    -- `id k < id i` reduces to `k.val < i.val`.
    simp only [id_eq, Fin.lt_def] at h_lt
    fin_cases i <;> fin_cases k <;>
      first
        | (exact absurd h_lt (by decide))
        | (simp only [hM_def, Matrix.submatrix_apply, hMρ_def, hMσ_def,
                       if_true, if_false,
                       Nat.one_ne_zero, hne_2_0, hne_2_1, hne_3_0, hne_3_1,
                       hne_3_2, hne_4_0, hne_4_1, hne_4_2, hne_4_3,
                       hne_5_0, hne_5_1, hne_5_2, hne_5_3, hne_5_4,
                       hne_6_0, hne_6_1, hne_6_2, hne_6_3, hne_6_4, hne_6_5]
           first | exact hL10 | exact hL20 | exact hL21
                 | exact hL30 | exact hL31 | exact hL32
                 | exact hL40 | exact hL41 | exact hL42 | exact hL43
                 | exact hL50 | exact hL51 | exact hL52 | exact hL53 | exact hL54
                 | exact hL60 | exact hL61 | exact hL62 | exact hL63 | exact hL64 | exact hL65)
  rw [Matrix.det_of_upperTriangular h_blocktri]
  -- Now expand the diagonal product over `Fin 7`.
  rw [Fin.prod_univ_seven]
  -- Reduce `M k k` for `k ∈ {0, 1, 2, 3, 4, 5, 6}` to the precomputed diagonal entries.
  simp only [hM_def, Matrix.submatrix_apply, hMρ_def, hMσ_def,
             h0, h1, h2, h3, h4, h5, h6, if_true, if_false,
             Nat.one_ne_zero, hne_2_0, hne_2_1, hne_3_0, hne_3_1, hne_3_2,
             hne_4_0, hne_4_1, hne_4_2, hne_4_3,
             hne_5_0, hne_5_1, hne_5_2, hne_5_3, hne_5_4,
             hne_6_0, hne_6_1, hne_6_2, hne_6_3, hne_6_4, hne_6_5]
  rw [hD00, hD11, hD22, hD33, hD44, hD55, hD66]
  -- Goal: `IsUnit (1 * 1 * 1 * 1 * 1 * 1 * 1 : ℚ)`.
  norm_num

/--
**Corner-case main theorem (W = 18, d = j + 1).** The unfolding rank is
exactly `7` for every `j ≥ 1`. Mirrors `mps_bond_dim_W_eq_12_d_eq_j_plus_1`
(both citing the general `upper_bound`) and is fully formalised in Lean
(no `sorry`, no `axiom`).

**Upper-bound subtlety:** the matrix has `18^((j+1)-j) = 18` columns, so
`rank_le_width` only gives `rank ≤ 18`. We instead cite the general
`upper_bound` lemma, which evaluates to `φ(18) · 18^0 + 1 = 6 · 1 + 1 = 7`
in this corner.

Together with the `W ∈ {2, 3, 4, 5, 6, 8, 12}` corners, this is the
ninth unconditional instance of `mps_bond_dim` and the seventh instance
over a wheel `W ≥ 3`. **First instance with `R = 7`.** The general
`mps_bond_dim` still requires `exists_invertible_submatrix` whose proof
is the only remaining `sorry` in this file.
-/
theorem mps_bond_dim_W_eq_18_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    (unfolding 18 (j + 1) j).rank = 7 := by
  apply Nat.le_antisymm
  · -- `≤ 7`: from `upper_bound`. (`rank_le_width` only gives `≤ 18`.)
    have hW : (2 : ℕ) ≤ 18 := by norm_num
    have hj_hi : j < j + 1 := Nat.lt_succ_self j
    have h := upper_bound 18 (j + 1) j hW hj hj_hi
    have h_eq : Nat.totient 18 * 18 ^ ((j + 1 : ℕ) - j - 1) + 1 = 7 := by
      have h_sub : (j + 1 : ℕ) - j - 1 = 0 := by omega
      rw [h_sub]
      decide
    linarith
  · -- `7 ≤ rank`: from the corner-case prime exhibit.
    obtain ⟨ρ, σ, hUnit⟩ := exists_invertible_submatrix_W_eq_18_d_eq_j_plus_1 j hj
    have h_eq : ((unfolding 18 (j + 1) j).submatrix ρ σ).rank = 7 := by
      have h := Matrix.rank_of_isUnit ((unfolding 18 (j + 1) j).submatrix ρ σ) hUnit
      simpa using h
    calc (7 : ℕ) = ((unfolding 18 (j + 1) j).submatrix ρ σ).rank := h_eq.symm
      _ ≤ (unfolding 18 (j + 1) j).rank := Matrix.rank_submatrix_le _ _ _

/-! ## W = 20 corner section (S143)

The `(W = 20, d = j + 1)` corner case of `mps_bond_dim`.
-/

/-- `chiP 47 = 1` since `47` is prime. -/
theorem chiP_forty_seven_eq_one : chiP 47 = 1 := by
  have h_prime_47 : Nat.Prime 47 := by norm_num
  simp [chiP, h_prime_47]

/-- `chiP 149 = 1` since `149` is prime. -/
theorem chiP_one_hundred_forty_nine_eq_one : chiP 149 = 1 := by
  have h_prime_149 : Nat.Prime 149 := by norm_num
  simp [chiP, h_prime_149]

/-- `chiP 199 = 1` since `199` is prime. -/
theorem chiP_one_hundred_ninety_nine_eq_one : chiP 199 = 1 := by
  have h_prime_199 : Nat.Prime 199 := by norm_num
  simp [chiP, h_prime_199]

/-- `chiP 241 = 1` since `241` is prime. -/
theorem chiP_two_hundred_forty_one_eq_one : chiP 241 = 1 := by
  have h_prime_241 : Nat.Prime 241 := by norm_num
  simp [chiP, h_prime_241]

/-- `chiP 337 = 1` since `337` is prime. -/
theorem chiP_three_hundred_thirty_seven_eq_one : chiP 337 = 1 := by
  have h_prime_337 : Nat.Prime 337 := by norm_num
  simp [chiP, h_prime_337]

/--
Local `prod_univ_nine` lemma. Mathlib provides `Fin.prod_univ_eight` but
not `Fin.prod_univ_nine`. We follow mathlib's pattern verbatim.
-/
private theorem prod_univ_nine' {M : Type*} [CommMonoid M] (f : Fin 9 → M) :
    ∏ i, f i = f 0 * f 1 * f 2 * f 3 * f 4 * f 5 * f 6 * f 7 * f 8 := by
  rw [Fin.prod_univ_castSucc, Fin.prod_univ_eight]
  rfl

set_option maxHeartbeats 2000000 in
-- The W=20 corner has R = 9, producing 81 fin_cases × fin_cases subgoals in the
-- BlockTriangular check (vs 49 at W=18, R=7). The simp blow-up scales as R² and
-- pushes simp past the default 200000 heartbeats. See S143 synthesis.
/--
**Corner-case prime exhibit (W = 20, d = j + 1).**

Specialisation of `exists_invertible_submatrix` to the corner case
`W = 20`, `d = j + 1`, where `R = min(20^j, φ(20)·1 + 1) = 9` for every `j ≥ 1`.
Closed unconditionally using `chiP_two_eq_one`, `chiP_twenty_three_eq_one`,
`chiP_forty_seven_eq_one`, `chiP_one_hundred_forty_nine_eq_one`,
`chiP_one_hundred_ninety_nine_eq_one`, `chiP_two_hundred_eleven_eq_one`,
`chiP_two_hundred_forty_one_eq_one`, `chiP_two_hundred_ninety_three_eq_one`,
and `chiP_three_hundred_thirty_seven_eq_one`.

The construction uses `ρ : Fin 9 → Fin (20^j)` mapping
`(0..8) ↦ (0, 2, 9, 14, 1, 7, 12, 16, 10)` and
`σ : Fin 9 → Fin (20^((j+1)-j))` mapping
`(0..8) ↦ (1, 6, 18, 12, 2, 8, 0, 16, 10)`.
The `9 × 9` submatrix is upper triangular with `1` on the diagonal
(diagonal primes `{2, 47, 199, 293, 23, 149, 241, 337, 211}`),
hence `IsUnit` over `ℚ`.

`maxHeartbeats` is bumped because the `R = 9` corner produces 81
`fin_cases` × `fin_cases` subgoals in the BlockTriangular check and
9-deep `if-then-else` chains in the diagonal expansion; the default
200000 is insufficient.
-/
theorem exists_invertible_submatrix_W_eq_20_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    ∃ (ρ : Fin 9 → Fin (20 ^ j))
      (σ : Fin 9 → Fin (20 ^ ((j + 1) - j))),
      IsUnit ((unfolding 20 (j + 1) j).submatrix ρ σ) := by
  -- The exponent simplifies: `(j + 1) - j = 1`, hence `20 ^ ((j + 1) - j) = 20`.
  have h_sub : (j + 1 : ℕ) - j = 1 := by omega
  have h_pow_j_pos : 0 < 20 ^ j := Nat.pow_pos (by norm_num)
  have h_one_lt_pow_j : 1 < 20 ^ j := by
    calc (1 : ℕ) < 20 := by norm_num
      _ = 20 ^ 1 := by norm_num
      _ ≤ 20 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_two_lt_pow_j : 2 < 20 ^ j := by
    calc (2 : ℕ) < 20 := by norm_num
      _ = 20 ^ 1 := by norm_num
      _ ≤ 20 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_seven_lt_pow_j : 7 < 20 ^ j := by
    calc (7 : ℕ) < 20 := by norm_num
      _ = 20 ^ 1 := by norm_num
      _ ≤ 20 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_nine_lt_pow_j : 9 < 20 ^ j := by
    calc (9 : ℕ) < 20 := by norm_num
      _ = 20 ^ 1 := by norm_num
      _ ≤ 20 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_ten_lt_pow_j : 10 < 20 ^ j := by
    calc (10 : ℕ) < 20 := by norm_num
      _ = 20 ^ 1 := by norm_num
      _ ≤ 20 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_twelve_lt_pow_j : 12 < 20 ^ j := by
    calc (12 : ℕ) < 20 := by norm_num
      _ = 20 ^ 1 := by norm_num
      _ ≤ 20 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_fourteen_lt_pow_j : 14 < 20 ^ j := by
    calc (14 : ℕ) < 20 := by norm_num
      _ = 20 ^ 1 := by norm_num
      _ ≤ 20 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_sixteen_lt_pow_j : 16 < 20 ^ j := by
    calc (16 : ℕ) < 20 := by norm_num
      _ = 20 ^ 1 := by norm_num
      _ ≤ 20 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_zero_lt_pow_dj : (0 : ℕ) < 20 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_one_lt_pow_dj : (1 : ℕ) < 20 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_two_lt_pow_dj : (2 : ℕ) < 20 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_six_lt_pow_dj : (6 : ℕ) < 20 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_eight_lt_pow_dj : (8 : ℕ) < 20 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_ten_lt_pow_dj : (10 : ℕ) < 20 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_twelve_lt_pow_dj : (12 : ℕ) < 20 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_sixteen_lt_pow_dj : (16 : ℕ) < 20 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_eighteen_lt_pow_dj : (18 : ℕ) < 20 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  -- ρ permutes original rows in the order `(0, 2, 9, 14, 1, 7, 12, 16, 10)`.
  -- σ permutes original columns in the order `(1, 6, 18, 12, 2, 8, 0, 16, 10)`.
  refine ⟨fun i =>
    if i.val = 0 then ⟨0, h_pow_j_pos⟩
                   else if i.val = 1 then ⟨2, h_two_lt_pow_j⟩
                   else if i.val = 2 then ⟨9, h_nine_lt_pow_j⟩
                   else if i.val = 3 then ⟨14, h_fourteen_lt_pow_j⟩
                   else if i.val = 4 then ⟨1, h_one_lt_pow_j⟩
                   else if i.val = 5 then ⟨7, h_seven_lt_pow_j⟩
                   else if i.val = 6 then ⟨12, h_twelve_lt_pow_j⟩
                   else if i.val = 7 then ⟨16, h_sixteen_lt_pow_j⟩
                   else ⟨10, h_ten_lt_pow_j⟩,
          fun i =>
    if i.val = 0 then ⟨1, h_one_lt_pow_dj⟩
                   else if i.val = 1 then ⟨6, h_six_lt_pow_dj⟩
                   else if i.val = 2 then ⟨18, h_eighteen_lt_pow_dj⟩
                   else if i.val = 3 then ⟨12, h_twelve_lt_pow_dj⟩
                   else if i.val = 4 then ⟨2, h_two_lt_pow_dj⟩
                   else if i.val = 5 then ⟨8, h_eight_lt_pow_dj⟩
                   else if i.val = 6 then ⟨0, h_zero_lt_pow_dj⟩
                   else if i.val = 7 then ⟨16, h_sixteen_lt_pow_dj⟩
                   else ⟨10, h_ten_lt_pow_dj⟩, ?_⟩
  -- Reduce `IsUnit` to `IsUnit (det)` and compute via `det_of_upperTriangular`.
  rw [Matrix.isUnit_iff_isUnit_det]
  -- Non-primality lemmas for the 36 below-diagonal zeros.
  have h_not_prime_22 : ¬ Nat.Prime 22 := by decide
  have h_not_prime_27 : ¬ Nat.Prime 27 := by decide
  have h_not_prime_33 : ¬ Nat.Prime 33 := by decide
  have h_not_prime_39 : ¬ Nat.Prime 39 := by decide
  have h_not_prime_42 : ¬ Nat.Prime 42 := by decide
  have h_not_prime_142 : ¬ Nat.Prime 142 := by norm_num
  have h_not_prime_143 : ¬ Nat.Prime 143 := by norm_num
  have h_not_prime_147 : ¬ Nat.Prime 147 := by norm_num
  have h_not_prime_153 : ¬ Nat.Prime 153 := by norm_num
  have h_not_prime_159 : ¬ Nat.Prime 159 := by norm_num
  have h_not_prime_182 : ¬ Nat.Prime 182 := by norm_num
  have h_not_prime_187 : ¬ Nat.Prime 187 := by norm_num
  have h_not_prime_201 : ¬ Nat.Prime 201 := by norm_num
  have h_not_prime_202 : ¬ Nat.Prime 202 := by norm_num
  have h_not_prime_203 : ¬ Nat.Prime 203 := by norm_num
  have h_not_prime_207 : ¬ Nat.Prime 207 := by norm_num
  have h_not_prime_209 : ¬ Nat.Prime 209 := by norm_num
  have h_not_prime_213 : ¬ Nat.Prime 213 := by norm_num
  have h_not_prime_217 : ¬ Nat.Prime 217 := by norm_num
  have h_not_prime_219 : ¬ Nat.Prime 219 := by norm_num
  have h_not_prime_242 : ¬ Nat.Prime 242 := by norm_num
  have h_not_prime_243 : ¬ Nat.Prime 243 := by norm_num
  have h_not_prime_247 : ¬ Nat.Prime 247 := by norm_num
  have h_not_prime_249 : ¬ Nat.Prime 249 := by norm_num
  have h_not_prime_253 : ¬ Nat.Prime 253 := by norm_num
  have h_not_prime_259 : ¬ Nat.Prime 259 := by norm_num
  have h_not_prime_282 : ¬ Nat.Prime 282 := by norm_num
  have h_not_prime_287 : ¬ Nat.Prime 287 := by norm_num
  have h_not_prime_299 : ¬ Nat.Prime 299 := by norm_num
  have h_not_prime_321 : ¬ Nat.Prime 321 := by norm_num
  have h_not_prime_322 : ¬ Nat.Prime 322 := by norm_num
  have h_not_prime_323 : ¬ Nat.Prime 323 := by norm_num
  have h_not_prime_327 : ¬ Nat.Prime 327 := by norm_num
  have h_not_prime_329 : ¬ Nat.Prime 329 := by norm_num
  have h_not_prime_333 : ¬ Nat.Prime 333 := by norm_num
  have h_not_prime_339 : ¬ Nat.Prime 339 := by norm_num
  -- Compute the 9 diagonal entries (all = 1).
  -- (0,0): unfolding(0, 1) = chiP 2 = 1
  have hD00 : unfolding 20 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨1, h_one_lt_pow_dj⟩ = 1 := by
    change chiP (0 * 20 ^ ((j + 1) - j) + 1 + 1) = 1
    simp [chiP, Nat.prime_two]
  -- (1,1): unfolding(2, 6) = chiP 47 = 1
  have hD11 : unfolding 20 (j + 1) j ⟨2, h_two_lt_pow_j⟩ ⟨6, h_six_lt_pow_dj⟩ = 1 := by
    change chiP (2 * 20 ^ ((j + 1) - j) + 6 + 1) = 1
    rw [h_sub]
    have h_eq : (2 * 20 ^ 1 + 6 + 1 : ℕ) = 47 := by norm_num
    rw [h_eq]
    exact chiP_forty_seven_eq_one
  -- (2,2): unfolding(9, 18) = chiP 199 = 1
  have hD22 : unfolding 20 (j + 1) j ⟨9, h_nine_lt_pow_j⟩ ⟨18, h_eighteen_lt_pow_dj⟩ = 1 := by
    change chiP (9 * 20 ^ ((j + 1) - j) + 18 + 1) = 1
    rw [h_sub]
    have h_eq : (9 * 20 ^ 1 + 18 + 1 : ℕ) = 199 := by norm_num
    rw [h_eq]
    exact chiP_one_hundred_ninety_nine_eq_one
  -- (3,3): unfolding(14, 12) = chiP 293 = 1
  have hD33 : unfolding 20 (j + 1) j ⟨14, h_fourteen_lt_pow_j⟩ ⟨12, h_twelve_lt_pow_dj⟩ = 1 := by
    change chiP (14 * 20 ^ ((j + 1) - j) + 12 + 1) = 1
    rw [h_sub]
    have h_eq : (14 * 20 ^ 1 + 12 + 1 : ℕ) = 293 := by norm_num
    rw [h_eq]
    exact chiP_two_hundred_ninety_three_eq_one
  -- (4,4): unfolding(1, 2) = chiP 23 = 1
  have hD44 : unfolding 20 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨2, h_two_lt_pow_dj⟩ = 1 := by
    change chiP (1 * 20 ^ ((j + 1) - j) + 2 + 1) = 1
    rw [h_sub]
    have h_eq : (1 * 20 ^ 1 + 2 + 1 : ℕ) = 23 := by norm_num
    rw [h_eq]
    exact chiP_twenty_three_eq_one
  -- (5,5): unfolding(7, 8) = chiP 149 = 1
  have hD55 : unfolding 20 (j + 1) j ⟨7, h_seven_lt_pow_j⟩ ⟨8, h_eight_lt_pow_dj⟩ = 1 := by
    change chiP (7 * 20 ^ ((j + 1) - j) + 8 + 1) = 1
    rw [h_sub]
    have h_eq : (7 * 20 ^ 1 + 8 + 1 : ℕ) = 149 := by norm_num
    rw [h_eq]
    exact chiP_one_hundred_forty_nine_eq_one
  -- (6,6): unfolding(12, 0) = chiP 241 = 1
  have hD66 : unfolding 20 (j + 1) j ⟨12, h_twelve_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 1 := by
    change chiP (12 * 20 ^ ((j + 1) - j) + 0 + 1) = 1
    rw [h_sub]
    have h_eq : (12 * 20 ^ 1 + 0 + 1 : ℕ) = 241 := by norm_num
    rw [h_eq]
    exact chiP_two_hundred_forty_one_eq_one
  -- (7,7): unfolding(16, 16) = chiP 337 = 1
  have hD77 : unfolding 20 (j + 1) j ⟨16, h_sixteen_lt_pow_j⟩ ⟨16, h_sixteen_lt_pow_dj⟩ = 1 := by
    change chiP (16 * 20 ^ ((j + 1) - j) + 16 + 1) = 1
    rw [h_sub]
    have h_eq : (16 * 20 ^ 1 + 16 + 1 : ℕ) = 337 := by norm_num
    rw [h_eq]
    exact chiP_three_hundred_thirty_seven_eq_one
  -- (8,8): unfolding(10, 10) = chiP 211 = 1
  have hD88 : unfolding 20 (j + 1) j ⟨10, h_ten_lt_pow_j⟩ ⟨10, h_ten_lt_pow_dj⟩ = 1 := by
    change chiP (10 * 20 ^ ((j + 1) - j) + 10 + 1) = 1
    rw [h_sub]
    have h_eq : (10 * 20 ^ 1 + 10 + 1 : ℕ) = 211 := by norm_num
    rw [h_eq]
    exact chiP_two_hundred_eleven_eq_one
  -- Compute the 36 below-diagonal entries (all = 0).
  -- (1,0): unfolding(2, 1) = chiP 42 = 0
  have hL10 : unfolding 20 (j + 1) j ⟨2, h_two_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (2 * 20 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (2 * 20 ^ 1 + 1 + 1 : ℕ) = 42 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_42]
  -- (2,0): unfolding(9, 1) = chiP 182 = 0
  have hL20 : unfolding 20 (j + 1) j ⟨9, h_nine_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (9 * 20 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (9 * 20 ^ 1 + 1 + 1 : ℕ) = 182 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_182]
  -- (2,1): unfolding(9, 6) = chiP 187 = 0
  have hL21 : unfolding 20 (j + 1) j ⟨9, h_nine_lt_pow_j⟩ ⟨6, h_six_lt_pow_dj⟩ = 0 := by
    change chiP (9 * 20 ^ ((j + 1) - j) + 6 + 1) = 0
    rw [h_sub]
    have h_eq : (9 * 20 ^ 1 + 6 + 1 : ℕ) = 187 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_187]
  -- (3,0): unfolding(14, 1) = chiP 282 = 0
  have hL30 : unfolding 20 (j + 1) j ⟨14, h_fourteen_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (14 * 20 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (14 * 20 ^ 1 + 1 + 1 : ℕ) = 282 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_282]
  -- (3,1): unfolding(14, 6) = chiP 287 = 0
  have hL31 : unfolding 20 (j + 1) j ⟨14, h_fourteen_lt_pow_j⟩ ⟨6, h_six_lt_pow_dj⟩ = 0 := by
    change chiP (14 * 20 ^ ((j + 1) - j) + 6 + 1) = 0
    rw [h_sub]
    have h_eq : (14 * 20 ^ 1 + 6 + 1 : ℕ) = 287 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_287]
  -- (3,2): unfolding(14, 18) = chiP 299 = 0
  have hL32 : unfolding 20 (j + 1) j ⟨14, h_fourteen_lt_pow_j⟩ ⟨18, h_eighteen_lt_pow_dj⟩ = 0 := by
    change chiP (14 * 20 ^ ((j + 1) - j) + 18 + 1) = 0
    rw [h_sub]
    have h_eq : (14 * 20 ^ 1 + 18 + 1 : ℕ) = 299 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_299]
  -- (4,0): unfolding(1, 1) = chiP 22 = 0
  have hL40 : unfolding 20 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (1 * 20 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (1 * 20 ^ 1 + 1 + 1 : ℕ) = 22 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_22]
  -- (4,1): unfolding(1, 6) = chiP 27 = 0
  have hL41 : unfolding 20 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨6, h_six_lt_pow_dj⟩ = 0 := by
    change chiP (1 * 20 ^ ((j + 1) - j) + 6 + 1) = 0
    rw [h_sub]
    have h_eq : (1 * 20 ^ 1 + 6 + 1 : ℕ) = 27 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_27]
  -- (4,2): unfolding(1, 18) = chiP 39 = 0
  have hL42 : unfolding 20 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨18, h_eighteen_lt_pow_dj⟩ = 0 := by
    change chiP (1 * 20 ^ ((j + 1) - j) + 18 + 1) = 0
    rw [h_sub]
    have h_eq : (1 * 20 ^ 1 + 18 + 1 : ℕ) = 39 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_39]
  -- (4,3): unfolding(1, 12) = chiP 33 = 0
  have hL43 : unfolding 20 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨12, h_twelve_lt_pow_dj⟩ = 0 := by
    change chiP (1 * 20 ^ ((j + 1) - j) + 12 + 1) = 0
    rw [h_sub]
    have h_eq : (1 * 20 ^ 1 + 12 + 1 : ℕ) = 33 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_33]
  -- (5,0): unfolding(7, 1) = chiP 142 = 0
  have hL50 : unfolding 20 (j + 1) j ⟨7, h_seven_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (7 * 20 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (7 * 20 ^ 1 + 1 + 1 : ℕ) = 142 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_142]
  -- (5,1): unfolding(7, 6) = chiP 147 = 0
  have hL51 : unfolding 20 (j + 1) j ⟨7, h_seven_lt_pow_j⟩ ⟨6, h_six_lt_pow_dj⟩ = 0 := by
    change chiP (7 * 20 ^ ((j + 1) - j) + 6 + 1) = 0
    rw [h_sub]
    have h_eq : (7 * 20 ^ 1 + 6 + 1 : ℕ) = 147 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_147]
  -- (5,2): unfolding(7, 18) = chiP 159 = 0
  have hL52 : unfolding 20 (j + 1) j ⟨7, h_seven_lt_pow_j⟩ ⟨18, h_eighteen_lt_pow_dj⟩ = 0 := by
    change chiP (7 * 20 ^ ((j + 1) - j) + 18 + 1) = 0
    rw [h_sub]
    have h_eq : (7 * 20 ^ 1 + 18 + 1 : ℕ) = 159 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_159]
  -- (5,3): unfolding(7, 12) = chiP 153 = 0
  have hL53 : unfolding 20 (j + 1) j ⟨7, h_seven_lt_pow_j⟩ ⟨12, h_twelve_lt_pow_dj⟩ = 0 := by
    change chiP (7 * 20 ^ ((j + 1) - j) + 12 + 1) = 0
    rw [h_sub]
    have h_eq : (7 * 20 ^ 1 + 12 + 1 : ℕ) = 153 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_153]
  -- (5,4): unfolding(7, 2) = chiP 143 = 0
  have hL54 : unfolding 20 (j + 1) j ⟨7, h_seven_lt_pow_j⟩ ⟨2, h_two_lt_pow_dj⟩ = 0 := by
    change chiP (7 * 20 ^ ((j + 1) - j) + 2 + 1) = 0
    rw [h_sub]
    have h_eq : (7 * 20 ^ 1 + 2 + 1 : ℕ) = 143 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_143]
  -- (6,0): unfolding(12, 1) = chiP 242 = 0
  have hL60 : unfolding 20 (j + 1) j ⟨12, h_twelve_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (12 * 20 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (12 * 20 ^ 1 + 1 + 1 : ℕ) = 242 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_242]
  -- (6,1): unfolding(12, 6) = chiP 247 = 0
  have hL61 : unfolding 20 (j + 1) j ⟨12, h_twelve_lt_pow_j⟩ ⟨6, h_six_lt_pow_dj⟩ = 0 := by
    change chiP (12 * 20 ^ ((j + 1) - j) + 6 + 1) = 0
    rw [h_sub]
    have h_eq : (12 * 20 ^ 1 + 6 + 1 : ℕ) = 247 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_247]
  -- (6,2): unfolding(12, 18) = chiP 259 = 0
  have hL62 : unfolding 20 (j + 1) j ⟨12, h_twelve_lt_pow_j⟩ ⟨18, h_eighteen_lt_pow_dj⟩ = 0 := by
    change chiP (12 * 20 ^ ((j + 1) - j) + 18 + 1) = 0
    rw [h_sub]
    have h_eq : (12 * 20 ^ 1 + 18 + 1 : ℕ) = 259 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_259]
  -- (6,3): unfolding(12, 12) = chiP 253 = 0
  have hL63 : unfolding 20 (j + 1) j ⟨12, h_twelve_lt_pow_j⟩ ⟨12, h_twelve_lt_pow_dj⟩ = 0 := by
    change chiP (12 * 20 ^ ((j + 1) - j) + 12 + 1) = 0
    rw [h_sub]
    have h_eq : (12 * 20 ^ 1 + 12 + 1 : ℕ) = 253 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_253]
  -- (6,4): unfolding(12, 2) = chiP 243 = 0
  have hL64 : unfolding 20 (j + 1) j ⟨12, h_twelve_lt_pow_j⟩ ⟨2, h_two_lt_pow_dj⟩ = 0 := by
    change chiP (12 * 20 ^ ((j + 1) - j) + 2 + 1) = 0
    rw [h_sub]
    have h_eq : (12 * 20 ^ 1 + 2 + 1 : ℕ) = 243 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_243]
  -- (6,5): unfolding(12, 8) = chiP 249 = 0
  have hL65 : unfolding 20 (j + 1) j ⟨12, h_twelve_lt_pow_j⟩ ⟨8, h_eight_lt_pow_dj⟩ = 0 := by
    change chiP (12 * 20 ^ ((j + 1) - j) + 8 + 1) = 0
    rw [h_sub]
    have h_eq : (12 * 20 ^ 1 + 8 + 1 : ℕ) = 249 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_249]
  -- (7,0): unfolding(16, 1) = chiP 322 = 0
  have hL70 : unfolding 20 (j + 1) j ⟨16, h_sixteen_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (16 * 20 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (16 * 20 ^ 1 + 1 + 1 : ℕ) = 322 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_322]
  -- (7,1): unfolding(16, 6) = chiP 327 = 0
  have hL71 : unfolding 20 (j + 1) j ⟨16, h_sixteen_lt_pow_j⟩ ⟨6, h_six_lt_pow_dj⟩ = 0 := by
    change chiP (16 * 20 ^ ((j + 1) - j) + 6 + 1) = 0
    rw [h_sub]
    have h_eq : (16 * 20 ^ 1 + 6 + 1 : ℕ) = 327 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_327]
  -- (7,2): unfolding(16, 18) = chiP 339 = 0
  have hL72 : unfolding 20 (j + 1) j ⟨16, h_sixteen_lt_pow_j⟩ ⟨18, h_eighteen_lt_pow_dj⟩ = 0 := by
    change chiP (16 * 20 ^ ((j + 1) - j) + 18 + 1) = 0
    rw [h_sub]
    have h_eq : (16 * 20 ^ 1 + 18 + 1 : ℕ) = 339 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_339]
  -- (7,3): unfolding(16, 12) = chiP 333 = 0
  have hL73 : unfolding 20 (j + 1) j ⟨16, h_sixteen_lt_pow_j⟩ ⟨12, h_twelve_lt_pow_dj⟩ = 0 := by
    change chiP (16 * 20 ^ ((j + 1) - j) + 12 + 1) = 0
    rw [h_sub]
    have h_eq : (16 * 20 ^ 1 + 12 + 1 : ℕ) = 333 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_333]
  -- (7,4): unfolding(16, 2) = chiP 323 = 0
  have hL74 : unfolding 20 (j + 1) j ⟨16, h_sixteen_lt_pow_j⟩ ⟨2, h_two_lt_pow_dj⟩ = 0 := by
    change chiP (16 * 20 ^ ((j + 1) - j) + 2 + 1) = 0
    rw [h_sub]
    have h_eq : (16 * 20 ^ 1 + 2 + 1 : ℕ) = 323 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_323]
  -- (7,5): unfolding(16, 8) = chiP 329 = 0
  have hL75 : unfolding 20 (j + 1) j ⟨16, h_sixteen_lt_pow_j⟩ ⟨8, h_eight_lt_pow_dj⟩ = 0 := by
    change chiP (16 * 20 ^ ((j + 1) - j) + 8 + 1) = 0
    rw [h_sub]
    have h_eq : (16 * 20 ^ 1 + 8 + 1 : ℕ) = 329 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_329]
  -- (7,6): unfolding(16, 0) = chiP 321 = 0
  have hL76 : unfolding 20 (j + 1) j ⟨16, h_sixteen_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (16 * 20 ^ ((j + 1) - j) + 0 + 1) = 0
    rw [h_sub]
    have h_eq : (16 * 20 ^ 1 + 0 + 1 : ℕ) = 321 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_321]
  -- (8,0): unfolding(10, 1) = chiP 202 = 0
  have hL80 : unfolding 20 (j + 1) j ⟨10, h_ten_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (10 * 20 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (10 * 20 ^ 1 + 1 + 1 : ℕ) = 202 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_202]
  -- (8,1): unfolding(10, 6) = chiP 207 = 0
  have hL81 : unfolding 20 (j + 1) j ⟨10, h_ten_lt_pow_j⟩ ⟨6, h_six_lt_pow_dj⟩ = 0 := by
    change chiP (10 * 20 ^ ((j + 1) - j) + 6 + 1) = 0
    rw [h_sub]
    have h_eq : (10 * 20 ^ 1 + 6 + 1 : ℕ) = 207 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_207]
  -- (8,2): unfolding(10, 18) = chiP 219 = 0
  have hL82 : unfolding 20 (j + 1) j ⟨10, h_ten_lt_pow_j⟩ ⟨18, h_eighteen_lt_pow_dj⟩ = 0 := by
    change chiP (10 * 20 ^ ((j + 1) - j) + 18 + 1) = 0
    rw [h_sub]
    have h_eq : (10 * 20 ^ 1 + 18 + 1 : ℕ) = 219 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_219]
  -- (8,3): unfolding(10, 12) = chiP 213 = 0
  have hL83 : unfolding 20 (j + 1) j ⟨10, h_ten_lt_pow_j⟩ ⟨12, h_twelve_lt_pow_dj⟩ = 0 := by
    change chiP (10 * 20 ^ ((j + 1) - j) + 12 + 1) = 0
    rw [h_sub]
    have h_eq : (10 * 20 ^ 1 + 12 + 1 : ℕ) = 213 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_213]
  -- (8,4): unfolding(10, 2) = chiP 203 = 0
  have hL84 : unfolding 20 (j + 1) j ⟨10, h_ten_lt_pow_j⟩ ⟨2, h_two_lt_pow_dj⟩ = 0 := by
    change chiP (10 * 20 ^ ((j + 1) - j) + 2 + 1) = 0
    rw [h_sub]
    have h_eq : (10 * 20 ^ 1 + 2 + 1 : ℕ) = 203 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_203]
  -- (8,5): unfolding(10, 8) = chiP 209 = 0
  have hL85 : unfolding 20 (j + 1) j ⟨10, h_ten_lt_pow_j⟩ ⟨8, h_eight_lt_pow_dj⟩ = 0 := by
    change chiP (10 * 20 ^ ((j + 1) - j) + 8 + 1) = 0
    rw [h_sub]
    have h_eq : (10 * 20 ^ 1 + 8 + 1 : ℕ) = 209 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_209]
  -- (8,6): unfolding(10, 0) = chiP 201 = 0
  have hL86 : unfolding 20 (j + 1) j ⟨10, h_ten_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (10 * 20 ^ ((j + 1) - j) + 0 + 1) = 0
    rw [h_sub]
    have h_eq : (10 * 20 ^ 1 + 0 + 1 : ℕ) = 201 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_201]
  -- (8,7): unfolding(10, 16) = chiP 217 = 0
  have hL87 : unfolding 20 (j + 1) j ⟨10, h_ten_lt_pow_j⟩ ⟨16, h_sixteen_lt_pow_dj⟩ = 0 := by
    change chiP (10 * 20 ^ ((j + 1) - j) + 16 + 1) = 0
    rw [h_sub]
    have h_eq : (10 * 20 ^ 1 + 16 + 1 : ℕ) = 217 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_217]
  -- Helpers for `Fin 9` `.val` and the `if-then-else` ρ, σ.
  have h0 : (0 : Fin 9).val = 0 := rfl
  have h1 : (1 : Fin 9).val = 1 := rfl
  have h2 : (2 : Fin 9).val = 2 := rfl
  have h3 : (3 : Fin 9).val = 3 := rfl
  have h4 : (4 : Fin 9).val = 4 := rfl
  have h5 : (5 : Fin 9).val = 5 := rfl
  have h6 : (6 : Fin 9).val = 6 := rfl
  have h7 : (7 : Fin 9).val = 7 := rfl
  have h8 : (8 : Fin 9).val = 8 := rfl
  have hne_2_0 : (2 : ℕ) ≠ 0 := by decide
  have hne_2_1 : (2 : ℕ) ≠ 1 := by decide
  have hne_3_0 : (3 : ℕ) ≠ 0 := by decide
  have hne_3_1 : (3 : ℕ) ≠ 1 := by decide
  have hne_3_2 : (3 : ℕ) ≠ 2 := by decide
  have hne_4_0 : (4 : ℕ) ≠ 0 := by decide
  have hne_4_1 : (4 : ℕ) ≠ 1 := by decide
  have hne_4_2 : (4 : ℕ) ≠ 2 := by decide
  have hne_4_3 : (4 : ℕ) ≠ 3 := by decide
  have hne_5_0 : (5 : ℕ) ≠ 0 := by decide
  have hne_5_1 : (5 : ℕ) ≠ 1 := by decide
  have hne_5_2 : (5 : ℕ) ≠ 2 := by decide
  have hne_5_3 : (5 : ℕ) ≠ 3 := by decide
  have hne_5_4 : (5 : ℕ) ≠ 4 := by decide
  have hne_6_0 : (6 : ℕ) ≠ 0 := by decide
  have hne_6_1 : (6 : ℕ) ≠ 1 := by decide
  have hne_6_2 : (6 : ℕ) ≠ 2 := by decide
  have hne_6_3 : (6 : ℕ) ≠ 3 := by decide
  have hne_6_4 : (6 : ℕ) ≠ 4 := by decide
  have hne_6_5 : (6 : ℕ) ≠ 5 := by decide
  have hne_7_0 : (7 : ℕ) ≠ 0 := by decide
  have hne_7_1 : (7 : ℕ) ≠ 1 := by decide
  have hne_7_2 : (7 : ℕ) ≠ 2 := by decide
  have hne_7_3 : (7 : ℕ) ≠ 3 := by decide
  have hne_7_4 : (7 : ℕ) ≠ 4 := by decide
  have hne_7_5 : (7 : ℕ) ≠ 5 := by decide
  have hne_7_6 : (7 : ℕ) ≠ 6 := by decide
  have hne_8_0 : (8 : ℕ) ≠ 0 := by decide
  have hne_8_1 : (8 : ℕ) ≠ 1 := by decide
  have hne_8_2 : (8 : ℕ) ≠ 2 := by decide
  have hne_8_3 : (8 : ℕ) ≠ 3 := by decide
  have hne_8_4 : (8 : ℕ) ≠ 4 := by decide
  have hne_8_5 : (8 : ℕ) ≠ 5 := by decide
  have hne_8_6 : (8 : ℕ) ≠ 6 := by decide
  have hne_8_7 : (8 : ℕ) ≠ 7 := by decide
  -- Establish the upper-triangular property of the submatrix.
  set Mρ : Fin 9 → Fin (20 ^ j) :=
    fun i => if i.val = 0 then ⟨0, h_pow_j_pos⟩
             else if i.val = 1 then ⟨2, h_two_lt_pow_j⟩
             else if i.val = 2 then ⟨9, h_nine_lt_pow_j⟩
             else if i.val = 3 then ⟨14, h_fourteen_lt_pow_j⟩
             else if i.val = 4 then ⟨1, h_one_lt_pow_j⟩
             else if i.val = 5 then ⟨7, h_seven_lt_pow_j⟩
             else if i.val = 6 then ⟨12, h_twelve_lt_pow_j⟩
             else if i.val = 7 then ⟨16, h_sixteen_lt_pow_j⟩
             else ⟨10, h_ten_lt_pow_j⟩ with hMρ_def
  set Mσ : Fin 9 → Fin (20 ^ ((j + 1) - j)) :=
    fun i => if i.val = 0 then ⟨1, h_one_lt_pow_dj⟩
             else if i.val = 1 then ⟨6, h_six_lt_pow_dj⟩
             else if i.val = 2 then ⟨18, h_eighteen_lt_pow_dj⟩
             else if i.val = 3 then ⟨12, h_twelve_lt_pow_dj⟩
             else if i.val = 4 then ⟨2, h_two_lt_pow_dj⟩
             else if i.val = 5 then ⟨8, h_eight_lt_pow_dj⟩
             else if i.val = 6 then ⟨0, h_zero_lt_pow_dj⟩
             else if i.val = 7 then ⟨16, h_sixteen_lt_pow_dj⟩
             else ⟨10, h_ten_lt_pow_dj⟩ with hMσ_def
  set M : Matrix (Fin 9) (Fin 9) ℚ :=
    (unfolding 20 (j + 1) j).submatrix Mρ Mσ with hM_def
  have h_blocktri : M.BlockTriangular id := by
    intro i k h_lt
    simp only [id_eq, Fin.lt_def] at h_lt
    fin_cases i <;> fin_cases k <;>
      first
        | (exact absurd h_lt (by decide))
        | (simp only [hM_def, Matrix.submatrix_apply, hMρ_def, hMσ_def,
                       if_true, if_false,
                       Nat.one_ne_zero, hne_2_0, hne_2_1, hne_3_0, hne_3_1, hne_3_2, hne_4_0, hne_4_1, hne_4_2, hne_4_3, hne_5_0, hne_5_1, hne_5_2, hne_5_3, hne_5_4, hne_6_0, hne_6_1, hne_6_2, hne_6_3, hne_6_4, hne_6_5, hne_7_0, hne_7_1, hne_7_2, hne_7_3, hne_7_4, hne_7_5, hne_7_6, hne_8_0, hne_8_1, hne_8_2, hne_8_3, hne_8_4, hne_8_5, hne_8_6, hne_8_7]
           first | exact hL10 | exact hL20 | exact hL21 | exact hL30 | exact hL31 | exact hL32 | exact hL40 | exact hL41 | exact hL42 | exact hL43 | exact hL50 | exact hL51 | exact hL52 | exact hL53 | exact hL54 | exact hL60 | exact hL61 | exact hL62 | exact hL63 | exact hL64 | exact hL65 | exact hL70 | exact hL71 | exact hL72 | exact hL73 | exact hL74 | exact hL75 | exact hL76 | exact hL80 | exact hL81 | exact hL82 | exact hL83 | exact hL84 | exact hL85 | exact hL86 | exact hL87)
  rw [Matrix.det_of_upperTriangular h_blocktri]
  rw [prod_univ_nine']
  simp only [hM_def, Matrix.submatrix_apply, hMρ_def, hMσ_def,
             h0, h1, h2, h3, h4, h5, h6, h7, h8, if_true, if_false,
             Nat.one_ne_zero, hne_2_0, hne_2_1, hne_3_0, hne_3_1, hne_3_2, hne_4_0, hne_4_1, hne_4_2, hne_4_3, hne_5_0, hne_5_1, hne_5_2, hne_5_3, hne_5_4, hne_6_0, hne_6_1, hne_6_2, hne_6_3, hne_6_4, hne_6_5, hne_7_0, hne_7_1, hne_7_2, hne_7_3, hne_7_4, hne_7_5, hne_7_6, hne_8_0, hne_8_1, hne_8_2, hne_8_3, hne_8_4, hne_8_5, hne_8_6, hne_8_7]
  rw [hD00, hD11, hD22, hD33, hD44, hD55, hD66, hD77, hD88]
  norm_num

/--
**Corner-case main theorem (W = 20, d = j + 1).** The unfolding rank is
exactly `9` for every `j ≥ 1`. Mirrors `mps_bond_dim_W_eq_18_d_eq_j_plus_1`
(both citing the general `upper_bound`) and is fully formalised in Lean
(no `sorry`, no `axiom`).

**Upper-bound subtlety:** the matrix has `20^((j+1)-j) = 20` columns, so
`rank_le_width` only gives `rank ≤ 20`. We instead cite the general
`upper_bound` lemma, which evaluates to `φ(20) · 20^0 + 1 = 8 · 1 + 1 = 9`
in this corner.

Together with the `W ∈ {2, 3, 4, 5, 6, 8, 12, 18}` corners, this is the
**tenth** unconditional instance of `mps_bond_dim` and the **eighth**
instance over a wheel `W ≥ 3`. **First instance with `R = 9`.** The general
`mps_bond_dim` still requires `exists_invertible_submatrix` whose proof
is the only remaining `sorry` in this file.
-/
theorem mps_bond_dim_W_eq_20_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    (unfolding 20 (j + 1) j).rank = 9 := by
  apply Nat.le_antisymm
  · -- `≤ 9`: from `upper_bound`. (`rank_le_width` only gives `≤ 20`.)
    have hW : (2 : ℕ) ≤ 20 := by norm_num
    have hj_hi : j < j + 1 := Nat.lt_succ_self j
    have h := upper_bound 20 (j + 1) j hW hj hj_hi
    have h_eq : Nat.totient 20 * 20 ^ ((j + 1 : ℕ) - j - 1) + 1 = 9 := by
      have h_sub : (j + 1 : ℕ) - j - 1 = 0 := by omega
      rw [h_sub]
      decide
    linarith
  · -- `9 ≤ rank`: from the corner-case prime exhibit.
    obtain ⟨ρ, σ, hUnit⟩ := exists_invertible_submatrix_W_eq_20_d_eq_j_plus_1 j hj
    have h_eq : ((unfolding 20 (j + 1) j).submatrix ρ σ).rank = 9 := by
      have h := Matrix.rank_of_isUnit ((unfolding 20 (j + 1) j).submatrix ρ σ) hUnit
      simpa using h
    calc (9 : ℕ) = ((unfolding 20 (j + 1) j).submatrix ρ σ).rank := h_eq.symm
      _ ≤ (unfolding 20 (j + 1) j).rank := Matrix.rank_submatrix_le _ _ _

/-!
### Corner case: `W = 10`, `d = j + 1` — closed unconditionally (S144)

When `W = 10` and `d = j + 1`, the formula gives
`R = min(10^j, φ(10) · 10^0 + 1) = min(10^j, 5) = 5` for every `j ≥ 1`.

**Construction.** Live columns at `W = 10` are `k ∈ {0, 2, 6, 8}`
(residues `1, 3, 7, 9 (mod 10)`), giving `φ(10) = 4` live columns. We
add the dead column `k = 1` (with `chiP 2 = 1` at row `0`), giving
`R = φ(10) + 1 = 5`.

The permutation `ρ ↦ (1, 0, 4, 3, 9)` and `σ ↦ (8, 1, 2, 0, 6)`
triangularises the `5 × 5` submatrix:
```
   row 1 ⎡ chiP 19, chiP 12, chiP 13, chiP 11, chiP 17 ⎤   ⎡ 1, 0, 1, 1, 1 ⎤
   row 0 ⎢ chiP  9, chiP  2, chiP  3, chiP  1, chiP  7 ⎥   ⎢ 0, 1, 1, 0, 1 ⎥
   row 4 ⎢ chiP 49, chiP 42, chiP 43, chiP 41, chiP 47 ⎥ = ⎢ 0, 0, 1, 1, 1 ⎥.
   row 3 ⎢ chiP 39, chiP 32, chiP 33, chiP 31, chiP 37 ⎥   ⎢ 0, 0, 0, 1, 1 ⎥
   row 9 ⎣ chiP 99, chiP 92, chiP 93, chiP 91, chiP 97 ⎦   ⎣ 0, 0, 0, 0, 1 ⎦
```
Diagonal primes are `{19, 2, 43, 31, 97}`; below-diagonal composites
are `{9, 32, 33, 39, 42, 49, 91, 92, 93, 99}`.

**S128/S129 retraction.** Earlier session notes listed W=10 as
"structurally obstructed" alongside `W ∈ {7, 9, 11}`, citing a
"multiplicity-2 residue pattern" argument. The S144 DP-based search
refuted this: a triangulation exists when row `9` is included (replacing
some leading row), with dead column `k = 1`. The earlier search was
inadvertently restricted to row prefixes `{0, 1, 2, 3, 4}`, which are
indeed insufficient at W=10; row `9`'s window `chiP 91..100` provides
`chiP 97` at the unique position breaking the rank-2 obstruction.

**Determinant via `det_of_upperTriangular`** (mathlib has no
`det_fin_five`, same situation as W=5/W=8/W=12). The proof mirrors
`W = 8` exactly.

**Upper-bound subtlety.** As with W=4/6/8/12, `rank_le_width` gives only
`rank ≤ 10`, not the sharp `rank ≤ 5 = φ(10) + 1`. We cite the general
`upper_bound`, which evaluates to `φ(10) · 10^0 + 1 = 4 · 1 + 1 = 5`.
-/

/--
`chiP 97 = 1` since `97` is prime.
-/
theorem chiP_ninety_seven_eq_one : chiP 97 = 1 := by
  have h_prime_97 : Nat.Prime 97 := by decide
  simp [chiP, h_prime_97]

/--
**Corner-case prime exhibit (W = 10, d = j + 1).**

Specialisation of `exists_invertible_submatrix` to the corner case
`W = 10`, `d = j + 1`, where `R = min(10^j, 5) = 5` for every `j ≥ 1`.
Closed unconditionally using `chiP_two_eq_one`, `chiP_nineteen_eq_one`
(S117), `chiP_thirty_one_eq_one` (S122), `chiP_forty_three_eq_one` (S137),
and `chiP_ninety_seven_eq_one` (new at S144).

The construction uses `ρ : Fin 5 → Fin (10^j)` mapping
`(0, 1, 2, 3, 4) ↦ (1, 0, 4, 3, 9)` and `σ : Fin 5 → Fin (10^((j+1)-j))`
mapping `(0, 1, 2, 3, 4) ↦ (8, 1, 2, 0, 6)`. The `5 × 5` submatrix is
upper triangular with `1` on the diagonal, hence `IsUnit` over `ℚ`.
-/
theorem exists_invertible_submatrix_W_eq_10_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    ∃ (ρ : Fin 5 → Fin (10 ^ j))
      (σ : Fin 5 → Fin (10 ^ ((j + 1) - j))),
      IsUnit ((unfolding 10 (j + 1) j).submatrix ρ σ) := by
  have h_sub : (j + 1 : ℕ) - j = 1 := by omega
  have h_pow_j_pos : 0 < 10 ^ j := Nat.pow_pos (by norm_num)
  have h_one_lt_pow_j : 1 < 10 ^ j := by
    calc (1 : ℕ) < 10 := by norm_num
      _ = 10 ^ 1 := by norm_num
      _ ≤ 10 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_three_lt_pow_j : 3 < 10 ^ j := by
    calc (3 : ℕ) < 10 := by norm_num
      _ = 10 ^ 1 := by norm_num
      _ ≤ 10 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_four_lt_pow_j : 4 < 10 ^ j := by
    calc (4 : ℕ) < 10 := by norm_num
      _ = 10 ^ 1 := by norm_num
      _ ≤ 10 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_nine_lt_pow_j : 9 < 10 ^ j := by
    calc (9 : ℕ) < 10 := by norm_num
      _ = 10 ^ 1 := by norm_num
      _ ≤ 10 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_zero_lt_pow_dj : (0 : ℕ) < 10 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_one_lt_pow_dj : (1 : ℕ) < 10 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_two_lt_pow_dj : (2 : ℕ) < 10 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_six_lt_pow_dj : (6 : ℕ) < 10 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_eight_lt_pow_dj : (8 : ℕ) < 10 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  -- ρ permutes original rows in the order `(1, 0, 4, 3, 9)`.
  -- σ permutes original columns in the order `(8, 1, 2, 0, 6)`.
  refine ⟨fun i => if i.val = 0 then ⟨1, h_one_lt_pow_j⟩
                   else if i.val = 1 then ⟨0, h_pow_j_pos⟩
                   else if i.val = 2 then ⟨4, h_four_lt_pow_j⟩
                   else if i.val = 3 then ⟨3, h_three_lt_pow_j⟩
                   else ⟨9, h_nine_lt_pow_j⟩,
          fun i => if i.val = 0 then ⟨8, h_eight_lt_pow_dj⟩
                   else if i.val = 1 then ⟨1, h_one_lt_pow_dj⟩
                   else if i.val = 2 then ⟨2, h_two_lt_pow_dj⟩
                   else if i.val = 3 then ⟨0, h_zero_lt_pow_dj⟩
                   else ⟨6, h_six_lt_pow_dj⟩, ?_⟩
  rw [Matrix.isUnit_iff_isUnit_det]
  -- Non-primality lemmas for the ten below-diagonal zeros.
  have h_not_prime_9  : ¬ Nat.Prime 9  := by decide
  have h_not_prime_32 : ¬ Nat.Prime 32 := by decide
  have h_not_prime_33 : ¬ Nat.Prime 33 := by decide
  have h_not_prime_39 : ¬ Nat.Prime 39 := by decide
  have h_not_prime_42 : ¬ Nat.Prime 42 := by decide
  have h_not_prime_49 : ¬ Nat.Prime 49 := by decide
  have h_not_prime_91 : ¬ Nat.Prime 91 := by decide
  have h_not_prime_92 : ¬ Nat.Prime 92 := by decide
  have h_not_prime_93 : ¬ Nat.Prime 93 := by decide
  have h_not_prime_99 : ¬ Nat.Prime 99 := by decide
  -- Compute the 5 diagonal entries (all = 1).
  -- (0,0): unfolding(1, 8) = chiP 19 = 1
  have hD00 : unfolding 10 (j + 1) j ⟨1, h_one_lt_pow_j⟩ ⟨8, h_eight_lt_pow_dj⟩ = 1 := by
    change chiP (1 * 10 ^ ((j + 1) - j) + 8 + 1) = 1
    rw [h_sub]
    have h_eq : (1 * 10 ^ 1 + 8 + 1 : ℕ) = 19 := by norm_num
    rw [h_eq]
    exact chiP_nineteen_eq_one
  -- (1,1): unfolding(0, 1) = chiP 2 = 1
  have hD11 : unfolding 10 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨1, h_one_lt_pow_dj⟩ = 1 := by
    change chiP (0 * 10 ^ ((j + 1) - j) + 1 + 1) = 1
    simp [chiP, Nat.prime_two]
  -- (2,2): unfolding(4, 2) = chiP 43 = 1
  have hD22 : unfolding 10 (j + 1) j ⟨4, h_four_lt_pow_j⟩ ⟨2, h_two_lt_pow_dj⟩ = 1 := by
    change chiP (4 * 10 ^ ((j + 1) - j) + 2 + 1) = 1
    rw [h_sub]
    have h_eq : (4 * 10 ^ 1 + 2 + 1 : ℕ) = 43 := by norm_num
    rw [h_eq]
    exact chiP_forty_three_eq_one
  -- (3,3): unfolding(3, 0) = chiP 31 = 1
  have hD33 : unfolding 10 (j + 1) j ⟨3, h_three_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 1 := by
    change chiP (3 * 10 ^ ((j + 1) - j) + 0 + 1) = 1
    rw [h_sub]
    have h_eq : (3 * 10 ^ 1 + 0 + 1 : ℕ) = 31 := by norm_num
    rw [h_eq]
    exact chiP_thirty_one_eq_one
  -- (4,4): unfolding(9, 6) = chiP 97 = 1
  have hD44 : unfolding 10 (j + 1) j ⟨9, h_nine_lt_pow_j⟩ ⟨6, h_six_lt_pow_dj⟩ = 1 := by
    change chiP (9 * 10 ^ ((j + 1) - j) + 6 + 1) = 1
    rw [h_sub]
    have h_eq : (9 * 10 ^ 1 + 6 + 1 : ℕ) = 97 := by norm_num
    rw [h_eq]
    exact chiP_ninety_seven_eq_one
  -- Compute the 10 below-diagonal entries (all = 0).
  -- (1,0): unfolding(0, 8) = chiP 9 = 0
  have hL10 : unfolding 10 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨8, h_eight_lt_pow_dj⟩ = 0 := by
    change chiP (0 * 10 ^ ((j + 1) - j) + 8 + 1) = 0
    rw [h_sub]
    have h_eq : (0 * 10 ^ 1 + 8 + 1 : ℕ) = 9 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_9]
  -- (2,0): unfolding(4, 8) = chiP 49 = 0
  have hL20 : unfolding 10 (j + 1) j ⟨4, h_four_lt_pow_j⟩ ⟨8, h_eight_lt_pow_dj⟩ = 0 := by
    change chiP (4 * 10 ^ ((j + 1) - j) + 8 + 1) = 0
    rw [h_sub]
    have h_eq : (4 * 10 ^ 1 + 8 + 1 : ℕ) = 49 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_49]
  -- (2,1): unfolding(4, 1) = chiP 42 = 0
  have hL21 : unfolding 10 (j + 1) j ⟨4, h_four_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (4 * 10 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (4 * 10 ^ 1 + 1 + 1 : ℕ) = 42 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_42]
  -- (3,0): unfolding(3, 8) = chiP 39 = 0
  have hL30 : unfolding 10 (j + 1) j ⟨3, h_three_lt_pow_j⟩ ⟨8, h_eight_lt_pow_dj⟩ = 0 := by
    change chiP (3 * 10 ^ ((j + 1) - j) + 8 + 1) = 0
    rw [h_sub]
    have h_eq : (3 * 10 ^ 1 + 8 + 1 : ℕ) = 39 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_39]
  -- (3,1): unfolding(3, 1) = chiP 32 = 0
  have hL31 : unfolding 10 (j + 1) j ⟨3, h_three_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (3 * 10 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (3 * 10 ^ 1 + 1 + 1 : ℕ) = 32 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_32]
  -- (3,2): unfolding(3, 2) = chiP 33 = 0
  have hL32 : unfolding 10 (j + 1) j ⟨3, h_three_lt_pow_j⟩ ⟨2, h_two_lt_pow_dj⟩ = 0 := by
    change chiP (3 * 10 ^ ((j + 1) - j) + 2 + 1) = 0
    rw [h_sub]
    have h_eq : (3 * 10 ^ 1 + 2 + 1 : ℕ) = 33 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_33]
  -- (4,0): unfolding(9, 8) = chiP 99 = 0
  have hL40 : unfolding 10 (j + 1) j ⟨9, h_nine_lt_pow_j⟩ ⟨8, h_eight_lt_pow_dj⟩ = 0 := by
    change chiP (9 * 10 ^ ((j + 1) - j) + 8 + 1) = 0
    rw [h_sub]
    have h_eq : (9 * 10 ^ 1 + 8 + 1 : ℕ) = 99 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_99]
  -- (4,1): unfolding(9, 1) = chiP 92 = 0
  have hL41 : unfolding 10 (j + 1) j ⟨9, h_nine_lt_pow_j⟩ ⟨1, h_one_lt_pow_dj⟩ = 0 := by
    change chiP (9 * 10 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (9 * 10 ^ 1 + 1 + 1 : ℕ) = 92 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_92]
  -- (4,2): unfolding(9, 2) = chiP 93 = 0
  have hL42 : unfolding 10 (j + 1) j ⟨9, h_nine_lt_pow_j⟩ ⟨2, h_two_lt_pow_dj⟩ = 0 := by
    change chiP (9 * 10 ^ ((j + 1) - j) + 2 + 1) = 0
    rw [h_sub]
    have h_eq : (9 * 10 ^ 1 + 2 + 1 : ℕ) = 93 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_93]
  -- (4,3): unfolding(9, 0) = chiP 91 = 0
  have hL43 : unfolding 10 (j + 1) j ⟨9, h_nine_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (9 * 10 ^ ((j + 1) - j) + 0 + 1) = 0
    rw [h_sub]
    have h_eq : (9 * 10 ^ 1 + 0 + 1 : ℕ) = 91 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_91]
  -- Helpers for `Fin 5` `.val` and the `if-then-else` ρ, σ.
  have h0 : (0 : Fin 5).val = 0 := rfl
  have h1 : (1 : Fin 5).val = 1 := rfl
  have h2 : (2 : Fin 5).val = 2 := rfl
  have h3 : (3 : Fin 5).val = 3 := rfl
  have h4 : (4 : Fin 5).val = 4 := rfl
  have hne_2_0 : (2 : ℕ) ≠ 0 := by decide
  have hne_2_1 : (2 : ℕ) ≠ 1 := by decide
  have hne_3_0 : (3 : ℕ) ≠ 0 := by decide
  have hne_3_1 : (3 : ℕ) ≠ 1 := by decide
  have hne_3_2 : (3 : ℕ) ≠ 2 := by decide
  have hne_4_0 : (4 : ℕ) ≠ 0 := by decide
  have hne_4_1 : (4 : ℕ) ≠ 1 := by decide
  have hne_4_2 : (4 : ℕ) ≠ 2 := by decide
  have hne_4_3 : (4 : ℕ) ≠ 3 := by decide
  -- Establish the upper-triangular property of the submatrix.
  set Mρ : Fin 5 → Fin (10 ^ j) :=
    fun i => if i.val = 0 then ⟨1, h_one_lt_pow_j⟩
             else if i.val = 1 then ⟨0, h_pow_j_pos⟩
             else if i.val = 2 then ⟨4, h_four_lt_pow_j⟩
             else if i.val = 3 then ⟨3, h_three_lt_pow_j⟩
             else ⟨9, h_nine_lt_pow_j⟩ with hMρ_def
  set Mσ : Fin 5 → Fin (10 ^ ((j + 1) - j)) :=
    fun i => if i.val = 0 then ⟨8, h_eight_lt_pow_dj⟩
             else if i.val = 1 then ⟨1, h_one_lt_pow_dj⟩
             else if i.val = 2 then ⟨2, h_two_lt_pow_dj⟩
             else if i.val = 3 then ⟨0, h_zero_lt_pow_dj⟩
             else ⟨6, h_six_lt_pow_dj⟩ with hMσ_def
  set M : Matrix (Fin 5) (Fin 5) ℚ :=
    (unfolding 10 (j + 1) j).submatrix Mρ Mσ with hM_def
  have h_blocktri : M.BlockTriangular id := by
    intro i k h_lt
    simp only [id_eq, Fin.lt_def] at h_lt
    fin_cases i <;> fin_cases k <;>
      first
        | (exact absurd h_lt (by decide))
        | (simp only [hM_def, Matrix.submatrix_apply, hMρ_def, hMσ_def,
                       if_true, if_false,
                       Nat.one_ne_zero, hne_2_0, hne_2_1, hne_3_0, hne_3_1,
                       hne_3_2, hne_4_0, hne_4_1, hne_4_2, hne_4_3]
           first | exact hL10 | exact hL20 | exact hL21
                 | exact hL30 | exact hL31 | exact hL32
                 | exact hL40 | exact hL41 | exact hL42 | exact hL43)
  rw [Matrix.det_of_upperTriangular h_blocktri]
  rw [Fin.prod_univ_five]
  simp only [hM_def, Matrix.submatrix_apply, hMρ_def, hMσ_def,
             h0, h1, h2, h3, h4, if_true, if_false,
             Nat.one_ne_zero, hne_2_0, hne_2_1, hne_3_0, hne_3_1, hne_3_2,
             hne_4_0, hne_4_1, hne_4_2, hne_4_3]
  rw [hD00, hD11, hD22, hD33, hD44]
  norm_num

/--
**Corner-case main theorem (W = 10, d = j + 1).** The unfolding rank
is exactly `5` for every `j ≥ 1`. Mirrors the W=8/W=12 corners (both
citing the general `upper_bound`) and is fully formalised in Lean
(no `sorry`, no `axiom`).

**Upper-bound subtlety:** the matrix has `10^((j+1)-j) = 10` columns,
so `rank_le_width` only gives `rank ≤ 10`. We instead cite the general
`upper_bound`, which evaluates to `φ(10) · 10^0 + 1 = 4 · 1 + 1 = 5`.

Together with `W ∈ {2, 3, 4, 5, 6, 8, 12, 18, 20}` corners, this is the
**eleventh** unconditional instance of `mps_bond_dim` and the **ninth**
instance over a wheel `W ≥ 3`. **First instance refuting an entry on the
S128/S129 "structurally obstructed" list** — W=10 turned out to be
algorithmically reachable via DP search (S144).
-/
theorem mps_bond_dim_W_eq_10_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    (unfolding 10 (j + 1) j).rank = 5 := by
  apply Nat.le_antisymm
  · -- `≤ 5`: from `upper_bound`. (`rank_le_width` only gives `≤ 10`.)
    have hW : (2 : ℕ) ≤ 10 := by norm_num
    have hj_hi : j < j + 1 := Nat.lt_succ_self j
    have h := upper_bound 10 (j + 1) j hW hj hj_hi
    have h_eq : Nat.totient 10 * 10 ^ ((j + 1 : ℕ) - j - 1) + 1 = 5 := by
      have h_sub : (j + 1 : ℕ) - j - 1 = 0 := by omega
      rw [h_sub]
      decide
    linarith
  · -- `5 ≤ rank`: from the corner-case prime exhibit.
    obtain ⟨ρ, σ, hUnit⟩ := exists_invertible_submatrix_W_eq_10_d_eq_j_plus_1 j hj
    have h_eq : ((unfolding 10 (j + 1) j).submatrix ρ σ).rank = 5 := by
      have h := Matrix.rank_of_isUnit ((unfolding 10 (j + 1) j).submatrix ρ σ) hUnit
      simpa using h
    calc (5 : ℕ) = ((unfolding 10 (j + 1) j).submatrix ρ σ).rank := h_eq.symm
      _ ≤ (unfolding 10 (j + 1) j).rank := Matrix.rank_submatrix_le _ _ _

/-!
### W=9 corner support (S152)

The W=9 corner case `(W = 9, d = j + 1)` requires four new prime witnesses:
`13, 41, 53, 61`. Together with existing helpers, these support the
`(1 + 3 + 3)` block-DIAGONAL decomposition discovered in the S151
pre-search. The corner closure itself is the file's first use of
`Matrix.det_fromBlocks_zero_21` (a different determinant technique
from the previous nine corner closures, all of which used
`Matrix.det_of_upperTriangular`).
-/

/-- `chiP 13 = 1` since `13` is prime. -/
theorem chiP_thirteen_eq_one : chiP 13 = 1 := by
  have h : Nat.Prime 13 := by decide
  simp [chiP, h]

/-- `chiP 41 = 1` since `41` is prime. -/
theorem chiP_forty_one_eq_one : chiP 41 = 1 := by
  have h : Nat.Prime 41 := by decide
  simp [chiP, h]

/-- `chiP 53 = 1` since `53` is prime. -/
theorem chiP_fifty_three_eq_one : chiP 53 = 1 := by
  have h : Nat.Prime 53 := by decide
  simp [chiP, h]

/-- `chiP 61 = 1` since `61` is prime. -/
theorem chiP_sixty_one_eq_one : chiP 61 = 1 := by
  have h : Nat.Prime 61 := by decide
  simp [chiP, h]



/-
**Corner-case prime exhibit (W = 9, d = j + 1).**

Closed at S152 via the BlockTriangular `(1 + 3 + 3)` decomposition
identified in S151's pre-search. The construction uses
`ρ : Fin 7 → Fin (9^j)` mapping `(0..6) ↦ (0, 1, 3, 5, 2, 4, 6)` and
`σ : Fin 7 → Fin (9^((j+1)-j))` mapping `(0..6) ↦ (2, 1, 3, 7, 0, 4, 6)`.
The `7 × 7` submatrix decomposes (via `finSumFinEquiv : Fin 4 ⊕ Fin 3 ≃ Fin 7`)
into `Matrix.fromBlocks A B 0 D` with `A : Matrix (Fin 4) (Fin 4) ℚ` (det -1)
and `D : Matrix (Fin 3) (Fin 3) ℚ` (det -1), giving total det = `(-1) · (-1) = 1`,
hence `IsUnit` over `ℚ`. Uses `chiP_thirteen_eq_one`, `chiP_forty_one_eq_one`,
`chiP_fifty_three_eq_one`, `chiP_sixty_one_eq_one` (new at S152) plus existing
helpers for primes `{2, 3, 5, 7, 11, 17, 19, 23, 29, 31, 37, 43, 47, 59}`.
First instance using `Matrix.det_fromBlocks_zero_21` (a determinant technique
orthogonal to the previous nine corner closures).
-/
set_option maxHeartbeats 800000 in
theorem exists_invertible_submatrix_W_eq_9_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    ∃ (ρ : Fin 7 → Fin (9 ^ j))
      (σ : Fin 7 → Fin (9 ^ ((j + 1) - j))),
      IsUnit ((unfolding 9 (j + 1) j).submatrix ρ σ) := by
  have h_sub : (j + 1 : ℕ) - j = 1 := by omega
  have h_pow_j_pos : 0 < 9 ^ j := Nat.pow_pos (by norm_num)
  have h_1_lt_pow_j : 1 < 9 ^ j := by
    calc (1 : ℕ) < 9 := by norm_num
      _ = 9 ^ 1 := by norm_num
      _ ≤ 9 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_2_lt_pow_j : 2 < 9 ^ j := by
    calc (2 : ℕ) < 9 := by norm_num
      _ = 9 ^ 1 := by norm_num
      _ ≤ 9 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_3_lt_pow_j : 3 < 9 ^ j := by
    calc (3 : ℕ) < 9 := by norm_num
      _ = 9 ^ 1 := by norm_num
      _ ≤ 9 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_4_lt_pow_j : 4 < 9 ^ j := by
    calc (4 : ℕ) < 9 := by norm_num
      _ = 9 ^ 1 := by norm_num
      _ ≤ 9 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_5_lt_pow_j : 5 < 9 ^ j := by
    calc (5 : ℕ) < 9 := by norm_num
      _ = 9 ^ 1 := by norm_num
      _ ≤ 9 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_6_lt_pow_j : 6 < 9 ^ j := by
    calc (6 : ℕ) < 9 := by norm_num
      _ = 9 ^ 1 := by norm_num
      _ ≤ 9 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_zero_lt_pow_dj : (0 : ℕ) < 9 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_1_lt_pow_dj : (1 : ℕ) < 9 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_2_lt_pow_dj : (2 : ℕ) < 9 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_3_lt_pow_dj : (3 : ℕ) < 9 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_4_lt_pow_dj : (4 : ℕ) < 9 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_6_lt_pow_dj : (6 : ℕ) < 9 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_7_lt_pow_dj : (7 : ℕ) < 9 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  -- ρ ↦ (0, 1, 3, 5, 2, 4, 6)
  -- σ ↦ (2, 1, 3, 7, 0, 4, 6)
  refine ⟨fun i =>
           if i.val = 0 then ⟨0, h_pow_j_pos⟩
           else if i.val = 1 then ⟨1, h_1_lt_pow_j⟩
           else if i.val = 2 then ⟨3, h_3_lt_pow_j⟩
           else if i.val = 3 then ⟨5, h_5_lt_pow_j⟩
           else if i.val = 4 then ⟨2, h_2_lt_pow_j⟩
           else if i.val = 5 then ⟨4, h_4_lt_pow_j⟩
           else ⟨6, h_6_lt_pow_j⟩,
         fun i =>
           if i.val = 0 then ⟨2, h_2_lt_pow_dj⟩
           else if i.val = 1 then ⟨1, h_1_lt_pow_dj⟩
           else if i.val = 2 then ⟨3, h_3_lt_pow_dj⟩
           else if i.val = 3 then ⟨7, h_7_lt_pow_dj⟩
           else if i.val = 4 then ⟨0, h_zero_lt_pow_dj⟩
           else if i.val = 5 then ⟨4, h_4_lt_pow_dj⟩
           else ⟨6, h_6_lt_pow_dj⟩, ?_⟩
  rw [Matrix.isUnit_iff_isUnit_det]
  -- Non-primality for the composites used in the entries.
  have h_not_prime_1 : ¬ Nat.Prime 1 := Nat.not_prime_one
  have h_not_prime_4 : ¬ Nat.Prime 4 := by decide
  have h_not_prime_8 : ¬ Nat.Prime 8 := by decide
  have h_not_prime_10 : ¬ Nat.Prime 10 := by decide
  have h_not_prime_12 : ¬ Nat.Prime 12 := by decide
  have h_not_prime_14 : ¬ Nat.Prime 14 := by decide
  have h_not_prime_16 : ¬ Nat.Prime 16 := by decide
  have h_not_prime_20 : ¬ Nat.Prime 20 := by decide
  have h_not_prime_21 : ¬ Nat.Prime 21 := by decide
  have h_not_prime_22 : ¬ Nat.Prime 22 := by decide
  have h_not_prime_25 : ¬ Nat.Prime 25 := by decide
  have h_not_prime_26 : ¬ Nat.Prime 26 := by decide
  have h_not_prime_28 : ¬ Nat.Prime 28 := by decide
  have h_not_prime_30 : ¬ Nat.Prime 30 := by decide
  have h_not_prime_32 : ¬ Nat.Prime 32 := by decide
  have h_not_prime_34 : ¬ Nat.Prime 34 := by decide
  have h_not_prime_35 : ¬ Nat.Prime 35 := by decide
  have h_not_prime_38 : ¬ Nat.Prime 38 := by decide
  have h_not_prime_39 : ¬ Nat.Prime 39 := by decide
  have h_not_prime_40 : ¬ Nat.Prime 40 := by decide
  have h_not_prime_44 : ¬ Nat.Prime 44 := by decide
  have h_not_prime_46 : ¬ Nat.Prime 46 := by decide
  have h_not_prime_48 : ¬ Nat.Prime 48 := by decide
  have h_not_prime_49 : ¬ Nat.Prime 49 := by decide
  have h_not_prime_50 : ¬ Nat.Prime 50 := by decide
  have h_not_prime_52 : ¬ Nat.Prime 52 := by decide
  have h_not_prime_55 : ¬ Nat.Prime 55 := by decide
  have h_not_prime_56 : ¬ Nat.Prime 56 := by decide
  have h_not_prime_57 : ¬ Nat.Prime 57 := by decide
  have h_not_prime_58 : ¬ Nat.Prime 58 := by decide
  have h_not_prime_62 : ¬ Nat.Prime 62 := by decide
  -- Compute the 49 entries of the 7×7 submatrix.
  have hE00 : unfolding 9 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨2, h_2_lt_pow_dj⟩ = 1 := by
    change chiP (0 * 9 ^ ((j + 1) - j) + 2 + 1) = 1
    rw [h_sub]
    have h_eq : (0 * 9 ^ 1 + 2 + 1 : ℕ) = 3 := by norm_num
    rw [h_eq]
    exact chiP_three_eq_one
  have hE01 : unfolding 9 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨1, h_1_lt_pow_dj⟩ = 1 := by
    change chiP (0 * 9 ^ ((j + 1) - j) + 1 + 1) = 1
    rw [h_sub]
    have h_eq : (0 * 9 ^ 1 + 1 + 1 : ℕ) = 2 := by norm_num
    rw [h_eq]
    exact chiP_two_eq_one
  have hE02 : unfolding 9 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨3, h_3_lt_pow_dj⟩ = 0 := by
    change chiP (0 * 9 ^ ((j + 1) - j) + 3 + 1) = 0
    rw [h_sub]
    have h_eq : (0 * 9 ^ 1 + 3 + 1 : ℕ) = 4 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_4]
  have hE03 : unfolding 9 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨7, h_7_lt_pow_dj⟩ = 0 := by
    change chiP (0 * 9 ^ ((j + 1) - j) + 7 + 1) = 0
    rw [h_sub]
    have h_eq : (0 * 9 ^ 1 + 7 + 1 : ℕ) = 8 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_8]
  have hE04 : unfolding 9 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (0 * 9 ^ ((j + 1) - j) + 0 + 1) = 0
    rw [h_sub]
    have h_eq : (0 * 9 ^ 1 + 0 + 1 : ℕ) = 1 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_1]
  have hE05 : unfolding 9 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨4, h_4_lt_pow_dj⟩ = 1 := by
    change chiP (0 * 9 ^ ((j + 1) - j) + 4 + 1) = 1
    rw [h_sub]
    have h_eq : (0 * 9 ^ 1 + 4 + 1 : ℕ) = 5 := by norm_num
    rw [h_eq]
    exact chiP_five_eq_one
  have hE06 : unfolding 9 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨6, h_6_lt_pow_dj⟩ = 1 := by
    change chiP (0 * 9 ^ ((j + 1) - j) + 6 + 1) = 1
    rw [h_sub]
    have h_eq : (0 * 9 ^ 1 + 6 + 1 : ℕ) = 7 := by norm_num
    rw [h_eq]
    exact chiP_seven_eq_one
  have hE10 : unfolding 9 (j + 1) j ⟨1, h_1_lt_pow_j⟩ ⟨2, h_2_lt_pow_dj⟩ = 0 := by
    change chiP (1 * 9 ^ ((j + 1) - j) + 2 + 1) = 0
    rw [h_sub]
    have h_eq : (1 * 9 ^ 1 + 2 + 1 : ℕ) = 12 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_12]
  have hE11 : unfolding 9 (j + 1) j ⟨1, h_1_lt_pow_j⟩ ⟨1, h_1_lt_pow_dj⟩ = 1 := by
    change chiP (1 * 9 ^ ((j + 1) - j) + 1 + 1) = 1
    rw [h_sub]
    have h_eq : (1 * 9 ^ 1 + 1 + 1 : ℕ) = 11 := by norm_num
    rw [h_eq]
    exact chiP_eleven_eq_one
  have hE12 : unfolding 9 (j + 1) j ⟨1, h_1_lt_pow_j⟩ ⟨3, h_3_lt_pow_dj⟩ = 1 := by
    change chiP (1 * 9 ^ ((j + 1) - j) + 3 + 1) = 1
    rw [h_sub]
    have h_eq : (1 * 9 ^ 1 + 3 + 1 : ℕ) = 13 := by norm_num
    rw [h_eq]
    exact chiP_thirteen_eq_one
  have hE13 : unfolding 9 (j + 1) j ⟨1, h_1_lt_pow_j⟩ ⟨7, h_7_lt_pow_dj⟩ = 1 := by
    change chiP (1 * 9 ^ ((j + 1) - j) + 7 + 1) = 1
    rw [h_sub]
    have h_eq : (1 * 9 ^ 1 + 7 + 1 : ℕ) = 17 := by norm_num
    rw [h_eq]
    exact chiP_seventeen_eq_one
  have hE14 : unfolding 9 (j + 1) j ⟨1, h_1_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (1 * 9 ^ ((j + 1) - j) + 0 + 1) = 0
    rw [h_sub]
    have h_eq : (1 * 9 ^ 1 + 0 + 1 : ℕ) = 10 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_10]
  have hE15 : unfolding 9 (j + 1) j ⟨1, h_1_lt_pow_j⟩ ⟨4, h_4_lt_pow_dj⟩ = 0 := by
    change chiP (1 * 9 ^ ((j + 1) - j) + 4 + 1) = 0
    rw [h_sub]
    have h_eq : (1 * 9 ^ 1 + 4 + 1 : ℕ) = 14 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_14]
  have hE16 : unfolding 9 (j + 1) j ⟨1, h_1_lt_pow_j⟩ ⟨6, h_6_lt_pow_dj⟩ = 0 := by
    change chiP (1 * 9 ^ ((j + 1) - j) + 6 + 1) = 0
    rw [h_sub]
    have h_eq : (1 * 9 ^ 1 + 6 + 1 : ℕ) = 16 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_16]
  have hE20 : unfolding 9 (j + 1) j ⟨3, h_3_lt_pow_j⟩ ⟨2, h_2_lt_pow_dj⟩ = 0 := by
    change chiP (3 * 9 ^ ((j + 1) - j) + 2 + 1) = 0
    rw [h_sub]
    have h_eq : (3 * 9 ^ 1 + 2 + 1 : ℕ) = 30 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_30]
  have hE21 : unfolding 9 (j + 1) j ⟨3, h_3_lt_pow_j⟩ ⟨1, h_1_lt_pow_dj⟩ = 1 := by
    change chiP (3 * 9 ^ ((j + 1) - j) + 1 + 1) = 1
    rw [h_sub]
    have h_eq : (3 * 9 ^ 1 + 1 + 1 : ℕ) = 29 := by norm_num
    rw [h_eq]
    exact chiP_twenty_nine_eq_one
  have hE22 : unfolding 9 (j + 1) j ⟨3, h_3_lt_pow_j⟩ ⟨3, h_3_lt_pow_dj⟩ = 1 := by
    change chiP (3 * 9 ^ ((j + 1) - j) + 3 + 1) = 1
    rw [h_sub]
    have h_eq : (3 * 9 ^ 1 + 3 + 1 : ℕ) = 31 := by norm_num
    rw [h_eq]
    exact chiP_thirty_one_eq_one
  have hE23 : unfolding 9 (j + 1) j ⟨3, h_3_lt_pow_j⟩ ⟨7, h_7_lt_pow_dj⟩ = 0 := by
    change chiP (3 * 9 ^ ((j + 1) - j) + 7 + 1) = 0
    rw [h_sub]
    have h_eq : (3 * 9 ^ 1 + 7 + 1 : ℕ) = 35 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_35]
  have hE24 : unfolding 9 (j + 1) j ⟨3, h_3_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (3 * 9 ^ ((j + 1) - j) + 0 + 1) = 0
    rw [h_sub]
    have h_eq : (3 * 9 ^ 1 + 0 + 1 : ℕ) = 28 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_28]
  have hE25 : unfolding 9 (j + 1) j ⟨3, h_3_lt_pow_j⟩ ⟨4, h_4_lt_pow_dj⟩ = 0 := by
    change chiP (3 * 9 ^ ((j + 1) - j) + 4 + 1) = 0
    rw [h_sub]
    have h_eq : (3 * 9 ^ 1 + 4 + 1 : ℕ) = 32 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_32]
  have hE26 : unfolding 9 (j + 1) j ⟨3, h_3_lt_pow_j⟩ ⟨6, h_6_lt_pow_dj⟩ = 0 := by
    change chiP (3 * 9 ^ ((j + 1) - j) + 6 + 1) = 0
    rw [h_sub]
    have h_eq : (3 * 9 ^ 1 + 6 + 1 : ℕ) = 34 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_34]
  have hE30 : unfolding 9 (j + 1) j ⟨5, h_5_lt_pow_j⟩ ⟨2, h_2_lt_pow_dj⟩ = 0 := by
    change chiP (5 * 9 ^ ((j + 1) - j) + 2 + 1) = 0
    rw [h_sub]
    have h_eq : (5 * 9 ^ 1 + 2 + 1 : ℕ) = 48 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_48]
  have hE31 : unfolding 9 (j + 1) j ⟨5, h_5_lt_pow_j⟩ ⟨1, h_1_lt_pow_dj⟩ = 1 := by
    change chiP (5 * 9 ^ ((j + 1) - j) + 1 + 1) = 1
    rw [h_sub]
    have h_eq : (5 * 9 ^ 1 + 1 + 1 : ℕ) = 47 := by norm_num
    rw [h_eq]
    exact chiP_forty_seven_eq_one
  have hE32 : unfolding 9 (j + 1) j ⟨5, h_5_lt_pow_j⟩ ⟨3, h_3_lt_pow_dj⟩ = 0 := by
    change chiP (5 * 9 ^ ((j + 1) - j) + 3 + 1) = 0
    rw [h_sub]
    have h_eq : (5 * 9 ^ 1 + 3 + 1 : ℕ) = 49 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_49]
  have hE33 : unfolding 9 (j + 1) j ⟨5, h_5_lt_pow_j⟩ ⟨7, h_7_lt_pow_dj⟩ = 1 := by
    change chiP (5 * 9 ^ ((j + 1) - j) + 7 + 1) = 1
    rw [h_sub]
    have h_eq : (5 * 9 ^ 1 + 7 + 1 : ℕ) = 53 := by norm_num
    rw [h_eq]
    exact chiP_fifty_three_eq_one
  have hE34 : unfolding 9 (j + 1) j ⟨5, h_5_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (5 * 9 ^ ((j + 1) - j) + 0 + 1) = 0
    rw [h_sub]
    have h_eq : (5 * 9 ^ 1 + 0 + 1 : ℕ) = 46 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_46]
  have hE35 : unfolding 9 (j + 1) j ⟨5, h_5_lt_pow_j⟩ ⟨4, h_4_lt_pow_dj⟩ = 0 := by
    change chiP (5 * 9 ^ ((j + 1) - j) + 4 + 1) = 0
    rw [h_sub]
    have h_eq : (5 * 9 ^ 1 + 4 + 1 : ℕ) = 50 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_50]
  have hE36 : unfolding 9 (j + 1) j ⟨5, h_5_lt_pow_j⟩ ⟨6, h_6_lt_pow_dj⟩ = 0 := by
    change chiP (5 * 9 ^ ((j + 1) - j) + 6 + 1) = 0
    rw [h_sub]
    have h_eq : (5 * 9 ^ 1 + 6 + 1 : ℕ) = 52 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_52]
  have hE40 : unfolding 9 (j + 1) j ⟨2, h_2_lt_pow_j⟩ ⟨2, h_2_lt_pow_dj⟩ = 0 := by
    change chiP (2 * 9 ^ ((j + 1) - j) + 2 + 1) = 0
    rw [h_sub]
    have h_eq : (2 * 9 ^ 1 + 2 + 1 : ℕ) = 21 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_21]
  have hE41 : unfolding 9 (j + 1) j ⟨2, h_2_lt_pow_j⟩ ⟨1, h_1_lt_pow_dj⟩ = 0 := by
    change chiP (2 * 9 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (2 * 9 ^ 1 + 1 + 1 : ℕ) = 20 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_20]
  have hE42 : unfolding 9 (j + 1) j ⟨2, h_2_lt_pow_j⟩ ⟨3, h_3_lt_pow_dj⟩ = 0 := by
    change chiP (2 * 9 ^ ((j + 1) - j) + 3 + 1) = 0
    rw [h_sub]
    have h_eq : (2 * 9 ^ 1 + 3 + 1 : ℕ) = 22 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_22]
  have hE43 : unfolding 9 (j + 1) j ⟨2, h_2_lt_pow_j⟩ ⟨7, h_7_lt_pow_dj⟩ = 0 := by
    change chiP (2 * 9 ^ ((j + 1) - j) + 7 + 1) = 0
    rw [h_sub]
    have h_eq : (2 * 9 ^ 1 + 7 + 1 : ℕ) = 26 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_26]
  have hE44 : unfolding 9 (j + 1) j ⟨2, h_2_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 1 := by
    change chiP (2 * 9 ^ ((j + 1) - j) + 0 + 1) = 1
    rw [h_sub]
    have h_eq : (2 * 9 ^ 1 + 0 + 1 : ℕ) = 19 := by norm_num
    rw [h_eq]
    exact chiP_nineteen_eq_one
  have hE45 : unfolding 9 (j + 1) j ⟨2, h_2_lt_pow_j⟩ ⟨4, h_4_lt_pow_dj⟩ = 1 := by
    change chiP (2 * 9 ^ ((j + 1) - j) + 4 + 1) = 1
    rw [h_sub]
    have h_eq : (2 * 9 ^ 1 + 4 + 1 : ℕ) = 23 := by norm_num
    rw [h_eq]
    exact chiP_twenty_three_eq_one
  have hE46 : unfolding 9 (j + 1) j ⟨2, h_2_lt_pow_j⟩ ⟨6, h_6_lt_pow_dj⟩ = 0 := by
    change chiP (2 * 9 ^ ((j + 1) - j) + 6 + 1) = 0
    rw [h_sub]
    have h_eq : (2 * 9 ^ 1 + 6 + 1 : ℕ) = 25 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_25]
  have hE50 : unfolding 9 (j + 1) j ⟨4, h_4_lt_pow_j⟩ ⟨2, h_2_lt_pow_dj⟩ = 0 := by
    change chiP (4 * 9 ^ ((j + 1) - j) + 2 + 1) = 0
    rw [h_sub]
    have h_eq : (4 * 9 ^ 1 + 2 + 1 : ℕ) = 39 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_39]
  have hE51 : unfolding 9 (j + 1) j ⟨4, h_4_lt_pow_j⟩ ⟨1, h_1_lt_pow_dj⟩ = 0 := by
    change chiP (4 * 9 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (4 * 9 ^ 1 + 1 + 1 : ℕ) = 38 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_38]
  have hE52 : unfolding 9 (j + 1) j ⟨4, h_4_lt_pow_j⟩ ⟨3, h_3_lt_pow_dj⟩ = 0 := by
    change chiP (4 * 9 ^ ((j + 1) - j) + 3 + 1) = 0
    rw [h_sub]
    have h_eq : (4 * 9 ^ 1 + 3 + 1 : ℕ) = 40 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_40]
  have hE53 : unfolding 9 (j + 1) j ⟨4, h_4_lt_pow_j⟩ ⟨7, h_7_lt_pow_dj⟩ = 0 := by
    change chiP (4 * 9 ^ ((j + 1) - j) + 7 + 1) = 0
    rw [h_sub]
    have h_eq : (4 * 9 ^ 1 + 7 + 1 : ℕ) = 44 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_44]
  have hE54 : unfolding 9 (j + 1) j ⟨4, h_4_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 1 := by
    change chiP (4 * 9 ^ ((j + 1) - j) + 0 + 1) = 1
    rw [h_sub]
    have h_eq : (4 * 9 ^ 1 + 0 + 1 : ℕ) = 37 := by norm_num
    rw [h_eq]
    exact chiP_thirty_seven_eq_one
  have hE55 : unfolding 9 (j + 1) j ⟨4, h_4_lt_pow_j⟩ ⟨4, h_4_lt_pow_dj⟩ = 1 := by
    change chiP (4 * 9 ^ ((j + 1) - j) + 4 + 1) = 1
    rw [h_sub]
    have h_eq : (4 * 9 ^ 1 + 4 + 1 : ℕ) = 41 := by norm_num
    rw [h_eq]
    exact chiP_forty_one_eq_one
  have hE56 : unfolding 9 (j + 1) j ⟨4, h_4_lt_pow_j⟩ ⟨6, h_6_lt_pow_dj⟩ = 1 := by
    change chiP (4 * 9 ^ ((j + 1) - j) + 6 + 1) = 1
    rw [h_sub]
    have h_eq : (4 * 9 ^ 1 + 6 + 1 : ℕ) = 43 := by norm_num
    rw [h_eq]
    exact chiP_forty_three_eq_one
  have hE60 : unfolding 9 (j + 1) j ⟨6, h_6_lt_pow_j⟩ ⟨2, h_2_lt_pow_dj⟩ = 0 := by
    change chiP (6 * 9 ^ ((j + 1) - j) + 2 + 1) = 0
    rw [h_sub]
    have h_eq : (6 * 9 ^ 1 + 2 + 1 : ℕ) = 57 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_57]
  have hE61 : unfolding 9 (j + 1) j ⟨6, h_6_lt_pow_j⟩ ⟨1, h_1_lt_pow_dj⟩ = 0 := by
    change chiP (6 * 9 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (6 * 9 ^ 1 + 1 + 1 : ℕ) = 56 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_56]
  have hE62 : unfolding 9 (j + 1) j ⟨6, h_6_lt_pow_j⟩ ⟨3, h_3_lt_pow_dj⟩ = 0 := by
    change chiP (6 * 9 ^ ((j + 1) - j) + 3 + 1) = 0
    rw [h_sub]
    have h_eq : (6 * 9 ^ 1 + 3 + 1 : ℕ) = 58 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_58]
  have hE63 : unfolding 9 (j + 1) j ⟨6, h_6_lt_pow_j⟩ ⟨7, h_7_lt_pow_dj⟩ = 0 := by
    change chiP (6 * 9 ^ ((j + 1) - j) + 7 + 1) = 0
    rw [h_sub]
    have h_eq : (6 * 9 ^ 1 + 7 + 1 : ℕ) = 62 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_62]
  have hE64 : unfolding 9 (j + 1) j ⟨6, h_6_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (6 * 9 ^ ((j + 1) - j) + 0 + 1) = 0
    rw [h_sub]
    have h_eq : (6 * 9 ^ 1 + 0 + 1 : ℕ) = 55 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_55]
  have hE65 : unfolding 9 (j + 1) j ⟨6, h_6_lt_pow_j⟩ ⟨4, h_4_lt_pow_dj⟩ = 1 := by
    change chiP (6 * 9 ^ ((j + 1) - j) + 4 + 1) = 1
    rw [h_sub]
    have h_eq : (6 * 9 ^ 1 + 4 + 1 : ℕ) = 59 := by norm_num
    rw [h_eq]
    exact chiP_fifty_nine_eq_one
  have hE66 : unfolding 9 (j + 1) j ⟨6, h_6_lt_pow_j⟩ ⟨6, h_6_lt_pow_dj⟩ = 1 := by
    change chiP (6 * 9 ^ ((j + 1) - j) + 6 + 1) = 1
    rw [h_sub]
    have h_eq : (6 * 9 ^ 1 + 6 + 1 : ℕ) = 61 := by norm_num
    rw [h_eq]
    exact chiP_sixty_one_eq_one
  -- Fin 7 .val and ≠ helpers for the if-then-else ρ, σ.
  have hv0 : (0 : Fin 7).val = 0 := rfl
  have hv1 : (1 : Fin 7).val = 1 := rfl
  have hv2 : (2 : Fin 7).val = 2 := rfl
  have hv3 : (3 : Fin 7).val = 3 := rfl
  have hv4 : (4 : Fin 7).val = 4 := rfl
  have hv5 : (5 : Fin 7).val = 5 := rfl
  have hv6 : (6 : Fin 7).val = 6 := rfl
  have hne_1_0 : (1 : ℕ) ≠ 0 := by decide
  have hne_2_0 : (2 : ℕ) ≠ 0 := by decide
  have hne_2_1 : (2 : ℕ) ≠ 1 := by decide
  have hne_3_0 : (3 : ℕ) ≠ 0 := by decide
  have hne_3_1 : (3 : ℕ) ≠ 1 := by decide
  have hne_3_2 : (3 : ℕ) ≠ 2 := by decide
  have hne_4_0 : (4 : ℕ) ≠ 0 := by decide
  have hne_4_1 : (4 : ℕ) ≠ 1 := by decide
  have hne_4_2 : (4 : ℕ) ≠ 2 := by decide
  have hne_4_3 : (4 : ℕ) ≠ 3 := by decide
  have hne_5_0 : (5 : ℕ) ≠ 0 := by decide
  have hne_5_1 : (5 : ℕ) ≠ 1 := by decide
  have hne_5_2 : (5 : ℕ) ≠ 2 := by decide
  have hne_5_3 : (5 : ℕ) ≠ 3 := by decide
  have hne_5_4 : (5 : ℕ) ≠ 4 := by decide
  have hne_6_0 : (6 : ℕ) ≠ 0 := by decide
  have hne_6_1 : (6 : ℕ) ≠ 1 := by decide
  have hne_6_2 : (6 : ℕ) ≠ 2 := by decide
  have hne_6_3 : (6 : ℕ) ≠ 3 := by decide
  have hne_6_4 : (6 : ℕ) ≠ 4 := by decide
  have hne_6_5 : (6 : ℕ) ≠ 5 := by decide
  -- Set up the named ρ, σ, M for the rest of the proof.
  set Mρ : Fin 7 → Fin (9 ^ j) :=
    fun i => if i.val = 0 then ⟨0, h_pow_j_pos⟩
             else if i.val = 1 then ⟨1, h_1_lt_pow_j⟩
             else if i.val = 2 then ⟨3, h_3_lt_pow_j⟩
             else if i.val = 3 then ⟨5, h_5_lt_pow_j⟩
             else if i.val = 4 then ⟨2, h_2_lt_pow_j⟩
             else if i.val = 5 then ⟨4, h_4_lt_pow_j⟩
             else ⟨6, h_6_lt_pow_j⟩ with hMρ_def
  set Mσ : Fin 7 → Fin (9 ^ ((j + 1) - j)) :=
    fun i => if i.val = 0 then ⟨2, h_2_lt_pow_dj⟩
             else if i.val = 1 then ⟨1, h_1_lt_pow_dj⟩
             else if i.val = 2 then ⟨3, h_3_lt_pow_dj⟩
             else if i.val = 3 then ⟨7, h_7_lt_pow_dj⟩
             else if i.val = 4 then ⟨0, h_zero_lt_pow_dj⟩
             else if i.val = 5 then ⟨4, h_4_lt_pow_dj⟩
             else ⟨6, h_6_lt_pow_dj⟩ with hMσ_def
  set M : Matrix (Fin 7) (Fin 7) ℚ :=
    (unfolding 9 (j + 1) j).submatrix Mρ Mσ with hM_def
  -- Explicit 7x7 matrix matching the entries.
  let Mexp : Matrix (Fin 7) (Fin 7) ℚ :=
    !![1, 1, 0, 0, 0, 1, 1;
       0, 1, 1, 1, 0, 0, 0;
       0, 1, 1, 0, 0, 0, 0;
       0, 1, 0, 1, 0, 0, 0;
       0, 0, 0, 0, 1, 1, 0;
       0, 0, 0, 0, 1, 1, 1;
       0, 0, 0, 0, 0, 1, 1]
  -- Prove M = Mexp.
  have hMeq : M = Mexp := by
    ext i k
    simp only [hM_def, Matrix.submatrix_apply, hMρ_def, hMσ_def]
    fin_cases i <;> fin_cases k <;>
      simp only [hv0, hv1, hv2, hv3, hv4, hv5, hv6,
                 if_true, if_false, Nat.one_ne_zero,
                 hne_1_0, hne_2_0, hne_2_1, hne_3_0, hne_3_1, hne_3_2,
                 hne_4_0, hne_4_1, hne_4_2, hne_4_3,
                 hne_5_0, hne_5_1, hne_5_2, hne_5_3, hne_5_4,
                 hne_6_0, hne_6_1, hne_6_2, hne_6_3, hne_6_4, hne_6_5]
    all_goals first
      | exact hE00
      | exact hE01
      | exact hE02
      | exact hE03
      | exact hE04
      | exact hE05
      | exact hE06
      | exact hE10
      | exact hE11
      | exact hE12
      | exact hE13
      | exact hE14
      | exact hE15
      | exact hE16
      | exact hE20
      | exact hE21
      | exact hE22
      | exact hE23
      | exact hE24
      | exact hE25
      | exact hE26
      | exact hE30
      | exact hE31
      | exact hE32
      | exact hE33
      | exact hE34
      | exact hE35
      | exact hE36
      | exact hE40
      | exact hE41
      | exact hE42
      | exact hE43
      | exact hE44
      | exact hE45
      | exact hE46
      | exact hE50
      | exact hE51
      | exact hE52
      | exact hE53
      | exact hE54
      | exact hE55
      | exact hE56
      | exact hE60
      | exact hE61
      | exact hE62
      | exact hE63
      | exact hE64
      | exact hE65
      | exact hE66
      | (simp [Mexp, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons, Matrix.cons_val_fin_one])
  rw [hMeq]
  -- Now show Mexp.det = 1 via two nested fromBlocks decompositions:
  --   Layer 1: Mexp = fromBlocks A B 0 D (1 + 6 split, lower-left zero).
  --   Layer 2: D   = fromBlocks D1 0 0 D2 (3 + 3 split, both off-diagonals zero).
  -- The 1×1 block A has det 1 (`det_fin_one`); the two 3×3 blocks D1, D2
  -- each have det -1 (`det_fin_three`); total det = 1 * (-1 * -1) = 1.
  -- Layer-1 blocks: 1+6 split.
  let A : Matrix (Fin 1) (Fin 1) ℚ := !![1]
  let B : Matrix (Fin 1) (Fin 6) ℚ := !![1, 0, 0, 0, 1, 1]
  let D : Matrix (Fin 6) (Fin 6) ℚ :=
    !![1, 1, 1, 0, 0, 0;
       1, 1, 0, 0, 0, 0;
       1, 0, 1, 0, 0, 0;
       0, 0, 0, 1, 1, 0;
       0, 0, 0, 1, 1, 1;
       0, 0, 0, 0, 1, 1]
  -- Layer-2 blocks: 3+3 split of D.
  let D1 : Matrix (Fin 3) (Fin 3) ℚ := !![1, 1, 1; 1, 1, 0; 1, 0, 1]
  let D2 : Matrix (Fin 3) (Fin 3) ℚ := !![1, 1, 0; 1, 1, 1; 0, 1, 1]
  -- Layer 1 equality: under finSumFinEquiv : Fin 1 ⊕ Fin 6 ≃ Fin 7.
  have h_fromBlocks :
      Mexp.submatrix finSumFinEquiv finSumFinEquiv = Matrix.fromBlocks A B 0 D := by
    ext i j
    rcases i with i | i <;> rcases j with j | j
    all_goals fin_cases i <;> fin_cases j <;> rfl
  -- Layer 2 equality: under finSumFinEquiv : Fin 3 ⊕ Fin 3 ≃ Fin 6.
  have h_fromBlocks_D :
      D.submatrix finSumFinEquiv finSumFinEquiv =
        Matrix.fromBlocks D1 (0 : Matrix (Fin 3) (Fin 3) ℚ) 0 D2 := by
    ext i j
    rcases i with i | i <;> rcases j with j | j
    all_goals fin_cases i <;> fin_cases j <;> rfl
  have h_Mexp_det :
      Mexp.det = (Matrix.fromBlocks A B (0 : Matrix (Fin 6) (Fin 1) ℚ) D).det := by
    rw [← h_fromBlocks]
    simp [Matrix.det_submatrix_equiv_self]
  rw [h_Mexp_det, Matrix.det_fromBlocks_zero₂₁]
  -- A.det = 1 (1×1 case).
  have hA_det : A.det = 1 := by
    show (!![(1 : ℚ)] : Matrix (Fin 1) (Fin 1) ℚ).det = 1
    simp [Matrix.det_fin_one]
  -- D.det = 1 via the layer-2 split: D1.det * D2.det = (-1) * (-1) = 1.
  have hD_det : D.det = 1 := by
    have h_D_det_eq :
        D.det =
          (Matrix.fromBlocks D1 (0 : Matrix (Fin 3) (Fin 3) ℚ) 0 D2).det := by
      rw [← h_fromBlocks_D]
      simp [Matrix.det_submatrix_equiv_self]
    rw [h_D_det_eq, Matrix.det_fromBlocks_zero₂₁]
    have hD1 : D1.det = -1 := by
      show (!![(1 : ℚ), 1, 1; 1, 1, 0; 1, 0, 1] :
              Matrix (Fin 3) (Fin 3) ℚ).det = -1
      simp [Matrix.det_fin_three]
    have hD2 : D2.det = -1 := by
      show (!![(1 : ℚ), 1, 0; 1, 1, 1; 0, 1, 1] :
              Matrix (Fin 3) (Fin 3) ℚ).det = -1
      simp [Matrix.det_fin_three]
    rw [hD1, hD2]
    norm_num
  rw [hA_det, hD_det]
  -- 1 * 1 = 1, IsUnit (1 : ℚ).
  norm_num

/--
**Corner-case main theorem (W = 9, d = j + 1).** The unfolding rank is
exactly `7` for every `j ≥ 1`. Mirrors `mps_bond_dim_W_eq_10_d_eq_j_plus_1`
(both citing the general `upper_bound`, since `rank_le_width` only gives
`rank ≤ 9`, not the sharp `rank ≤ 7 = φ(9) + 1`). Fully formalised in
Lean (no `sorry`, no `axiom`).

**Upper-bound subtlety:** the matrix has `9^((j+1)-j) = 9` columns, so
`rank_le_width` gives only `rank ≤ 9`. We instead cite the general
`upper_bound`, which evaluates to `φ(9) · 9^0 + 1 = 6 · 1 + 1 = 7`.

**Proof technique novelty.** This is the **first** corner closure using
`Matrix.det_fromBlocks_zero₂₁` (i.e., a non-triangular block-diagonal
decomposition). Previous nine corner closures (W ∈ {2..6, 8, 10, 12, 18,
20}) all used `Matrix.det_of_upperTriangular`, which is structurally
unavailable here: S144's exhaustive DP search confirmed the W=9 slab
admits NO leading-row + dead-col upper-triangulation in rows `[0, 9)`.
The S151 pre-search identified the `(1 + 3 + 3)` block-DIAGONAL pattern
at row set `{0, 1, 3, 5, 2, 4, 6}`, col set `{2, 1, 3, 7, 0, 4, 6}`.

Together with prior corners, this is the **eleventh** unconditional
instance of `mps_bond_dim` and the **tenth** instance over a wheel
`W ≥ 3`.
-/
theorem mps_bond_dim_W_eq_9_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    (unfolding 9 (j + 1) j).rank = 7 := by
  apply Nat.le_antisymm
  · -- `≤ 7`: from `upper_bound`. (`rank_le_width` only gives `≤ 9`.)
    have hW : (2 : ℕ) ≤ 9 := by norm_num
    have hj_hi : j < j + 1 := Nat.lt_succ_self j
    have h := upper_bound 9 (j + 1) j hW hj hj_hi
    have h_eq : Nat.totient 9 * 9 ^ ((j + 1 : ℕ) - j - 1) + 1 = 7 := by
      have h_sub : (j + 1 : ℕ) - j - 1 = 0 := by omega
      rw [h_sub]
      decide
    linarith
  · -- `7 ≤ rank`: from the corner-case prime exhibit.
    obtain ⟨ρ, σ, hUnit⟩ := exists_invertible_submatrix_W_eq_9_d_eq_j_plus_1 j hj
    have h_eq : ((unfolding 9 (j + 1) j).submatrix ρ σ).rank = 7 := by
      have h := Matrix.rank_of_isUnit ((unfolding 9 (j + 1) j).submatrix ρ σ) hUnit
      simpa using h
    calc (7 : ℕ) = ((unfolding 9 (j + 1) j).submatrix ρ σ).rank := h_eq.symm
      _ ≤ (unfolding 9 (j + 1) j).rank := Matrix.rank_submatrix_le _ _ _

/-!
### W=7 corner support (S159)

The W=7 corner case `(W = 7, d = j + 1)` reuses the `(1 + 3 + 3)`
block-DIAGONAL Lean assembly developed at S152 (W=9). W=7 is
**structurally a small prime corner with NO new prime helpers**: the
S159 pre-search (`w7_blocktriangular_search.py`) shows that the
permutation `ρ ↦ (0, 1, 3, 5, 2, 4, 6)` and `σ ↦ (6, 1, 3, 5, 0, 2, 4)`
gives a `1 + 3 + 3` block-DIAGONAL structure using only primes
`{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47}` — every one
of which is already declared in the file.

Determinants: `D₁.det = -1`, `D₂.det = -2`, `A.det = 1`, total
`det = 1 · (-1) · (-2) = 2 ≠ 0`, hence `IsUnit` over `ℚ`. (W=9 had
total det `1`; W=7's nonzero-but-not-`±1` total is the first instance
where the closing `IsUnit` requires `Ne.isUnit` rather than `IsUnit 1`.)
-/

/-
**Corner-case prime exhibit (W = 7, d = j + 1).**

Closed at S159 via the BlockTriangular `(1 + 3 + 3)` decomposition
identified in the S159 pre-search. The construction uses
`ρ : Fin 7 → Fin (7^j)` mapping `(0..6) ↦ (0, 1, 3, 5, 2, 4, 6)` and
`σ : Fin 7 → Fin (7^((j+1)-j))` mapping `(0..6) ↦ (6, 1, 3, 5, 0, 2, 4)`.
The `7 × 7` submatrix decomposes (via `finSumFinEquiv : Fin 1 ⊕ Fin 6 ≃ Fin 7`)
into `Matrix.fromBlocks A B 0 D` with `A : Matrix (Fin 1) (Fin 1) ℚ` (det 1)
and `D : Matrix (Fin 6) (Fin 6) ℚ` (det 2 via inner 3+3 split: D₁.det = -1,
D₂.det = -2), giving total det = `1 · 2 = 2`, hence `IsUnit` over `ℚ`.
NO new prime helpers — uses existing primes
`{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47}`.
-/
set_option maxHeartbeats 800000 in
theorem exists_invertible_submatrix_W_eq_7_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    ∃ (ρ : Fin 7 → Fin (7 ^ j))
      (σ : Fin 7 → Fin (7 ^ ((j + 1) - j))),
      IsUnit ((unfolding 7 (j + 1) j).submatrix ρ σ) := by
  have h_sub : (j + 1 : ℕ) - j = 1 := by omega
  have h_pow_j_pos : 0 < 7 ^ j := Nat.pow_pos (by norm_num)
  have h_1_lt_pow_j : 1 < 7 ^ j := by
    calc (1 : ℕ) < 7 := by norm_num
      _ = 7 ^ 1 := by norm_num
      _ ≤ 7 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_2_lt_pow_j : 2 < 7 ^ j := by
    calc (2 : ℕ) < 7 := by norm_num
      _ = 7 ^ 1 := by norm_num
      _ ≤ 7 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_3_lt_pow_j : 3 < 7 ^ j := by
    calc (3 : ℕ) < 7 := by norm_num
      _ = 7 ^ 1 := by norm_num
      _ ≤ 7 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_4_lt_pow_j : 4 < 7 ^ j := by
    calc (4 : ℕ) < 7 := by norm_num
      _ = 7 ^ 1 := by norm_num
      _ ≤ 7 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_5_lt_pow_j : 5 < 7 ^ j := by
    calc (5 : ℕ) < 7 := by norm_num
      _ = 7 ^ 1 := by norm_num
      _ ≤ 7 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_6_lt_pow_j : 6 < 7 ^ j := by
    calc (6 : ℕ) < 7 := by norm_num
      _ = 7 ^ 1 := by norm_num
      _ ≤ 7 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h_zero_lt_pow_dj : (0 : ℕ) < 7 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_1_lt_pow_dj : (1 : ℕ) < 7 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_2_lt_pow_dj : (2 : ℕ) < 7 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_3_lt_pow_dj : (3 : ℕ) < 7 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_4_lt_pow_dj : (4 : ℕ) < 7 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_5_lt_pow_dj : (5 : ℕ) < 7 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  have h_6_lt_pow_dj : (6 : ℕ) < 7 ^ ((j + 1) - j) := by rw [h_sub]; norm_num
  -- ρ ↦ (0, 1, 3, 5, 2, 4, 6)
  -- σ ↦ (6, 1, 3, 5, 0, 2, 4)
  refine ⟨fun i =>
           if i.val = 0 then ⟨0, h_pow_j_pos⟩
           else if i.val = 1 then ⟨1, h_1_lt_pow_j⟩
           else if i.val = 2 then ⟨3, h_3_lt_pow_j⟩
           else if i.val = 3 then ⟨5, h_5_lt_pow_j⟩
           else if i.val = 4 then ⟨2, h_2_lt_pow_j⟩
           else if i.val = 5 then ⟨4, h_4_lt_pow_j⟩
           else ⟨6, h_6_lt_pow_j⟩,
         fun i =>
           if i.val = 0 then ⟨6, h_6_lt_pow_dj⟩
           else if i.val = 1 then ⟨1, h_1_lt_pow_dj⟩
           else if i.val = 2 then ⟨3, h_3_lt_pow_dj⟩
           else if i.val = 3 then ⟨5, h_5_lt_pow_dj⟩
           else if i.val = 4 then ⟨0, h_zero_lt_pow_dj⟩
           else if i.val = 5 then ⟨2, h_2_lt_pow_dj⟩
           else ⟨4, h_4_lt_pow_dj⟩, ?_⟩
  rw [Matrix.isUnit_iff_isUnit_det]
  -- Non-primality for the composites used in the entries.
  have h_not_prime_1 : ¬ Nat.Prime 1 := Nat.not_prime_one
  have h_not_prime_4 : ¬ Nat.Prime 4 := by decide
  have h_not_prime_6 : ¬ Nat.Prime 6 := by decide
  have h_not_prime_8 : ¬ Nat.Prime 8 := by decide
  have h_not_prime_9 : ¬ Nat.Prime 9 := by decide
  have h_not_prime_10 : ¬ Nat.Prime 10 := by decide
  have h_not_prime_12 : ¬ Nat.Prime 12 := by decide
  have h_not_prime_14 : ¬ Nat.Prime 14 := by decide
  have h_not_prime_15 : ¬ Nat.Prime 15 := by decide
  have h_not_prime_16 : ¬ Nat.Prime 16 := by decide
  have h_not_prime_18 : ¬ Nat.Prime 18 := by decide
  have h_not_prime_20 : ¬ Nat.Prime 20 := by decide
  have h_not_prime_21 : ¬ Nat.Prime 21 := by decide
  have h_not_prime_22 : ¬ Nat.Prime 22 := by decide
  have h_not_prime_24 : ¬ Nat.Prime 24 := by decide
  have h_not_prime_25 : ¬ Nat.Prime 25 := by decide
  have h_not_prime_26 : ¬ Nat.Prime 26 := by decide
  have h_not_prime_27 : ¬ Nat.Prime 27 := by decide
  have h_not_prime_28 : ¬ Nat.Prime 28 := by decide
  have h_not_prime_30 : ¬ Nat.Prime 30 := by decide
  have h_not_prime_32 : ¬ Nat.Prime 32 := by decide
  have h_not_prime_33 : ¬ Nat.Prime 33 := by decide
  have h_not_prime_34 : ¬ Nat.Prime 34 := by decide
  have h_not_prime_35 : ¬ Nat.Prime 35 := by decide
  have h_not_prime_36 : ¬ Nat.Prime 36 := by decide
  have h_not_prime_38 : ¬ Nat.Prime 38 := by decide
  have h_not_prime_39 : ¬ Nat.Prime 39 := by decide
  have h_not_prime_40 : ¬ Nat.Prime 40 := by decide
  have h_not_prime_42 : ¬ Nat.Prime 42 := by decide
  have h_not_prime_44 : ¬ Nat.Prime 44 := by decide
  have h_not_prime_45 : ¬ Nat.Prime 45 := by decide
  have h_not_prime_46 : ¬ Nat.Prime 46 := by decide
  have h_not_prime_48 : ¬ Nat.Prime 48 := by decide
  have h_not_prime_49 : ¬ Nat.Prime 49 := by decide
  -- Compute the 49 entries of the 7×7 submatrix.
  -- Row 0 (rho=0): n = 0*7 + sigma[k] + 1.
  have hE00 : unfolding 7 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨6, h_6_lt_pow_dj⟩ = 1 := by
    change chiP (0 * 7 ^ ((j + 1) - j) + 6 + 1) = 1
    rw [h_sub]
    have h_eq : (0 * 7 ^ 1 + 6 + 1 : ℕ) = 7 := by norm_num
    rw [h_eq]
    exact chiP_seven_eq_one
  have hE01 : unfolding 7 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨1, h_1_lt_pow_dj⟩ = 1 := by
    change chiP (0 * 7 ^ ((j + 1) - j) + 1 + 1) = 1
    rw [h_sub]
    have h_eq : (0 * 7 ^ 1 + 1 + 1 : ℕ) = 2 := by norm_num
    rw [h_eq]
    exact chiP_two_eq_one
  have hE02 : unfolding 7 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨3, h_3_lt_pow_dj⟩ = 0 := by
    change chiP (0 * 7 ^ ((j + 1) - j) + 3 + 1) = 0
    rw [h_sub]
    have h_eq : (0 * 7 ^ 1 + 3 + 1 : ℕ) = 4 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_4]
  have hE03 : unfolding 7 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨5, h_5_lt_pow_dj⟩ = 0 := by
    change chiP (0 * 7 ^ ((j + 1) - j) + 5 + 1) = 0
    rw [h_sub]
    have h_eq : (0 * 7 ^ 1 + 5 + 1 : ℕ) = 6 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_6]
  have hE04 : unfolding 7 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (0 * 7 ^ ((j + 1) - j) + 0 + 1) = 0
    rw [h_sub]
    have h_eq : (0 * 7 ^ 1 + 0 + 1 : ℕ) = 1 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_1]
  have hE05 : unfolding 7 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨2, h_2_lt_pow_dj⟩ = 1 := by
    change chiP (0 * 7 ^ ((j + 1) - j) + 2 + 1) = 1
    rw [h_sub]
    have h_eq : (0 * 7 ^ 1 + 2 + 1 : ℕ) = 3 := by norm_num
    rw [h_eq]
    exact chiP_three_eq_one
  have hE06 : unfolding 7 (j + 1) j ⟨0, h_pow_j_pos⟩ ⟨4, h_4_lt_pow_dj⟩ = 1 := by
    change chiP (0 * 7 ^ ((j + 1) - j) + 4 + 1) = 1
    rw [h_sub]
    have h_eq : (0 * 7 ^ 1 + 4 + 1 : ℕ) = 5 := by norm_num
    rw [h_eq]
    exact chiP_five_eq_one
  -- Row 1 (rho=1): n = 1*7 + sigma[k] + 1 = sigma[k] + 8.
  have hE10 : unfolding 7 (j + 1) j ⟨1, h_1_lt_pow_j⟩ ⟨6, h_6_lt_pow_dj⟩ = 0 := by
    change chiP (1 * 7 ^ ((j + 1) - j) + 6 + 1) = 0
    rw [h_sub]
    have h_eq : (1 * 7 ^ 1 + 6 + 1 : ℕ) = 14 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_14]
  have hE11 : unfolding 7 (j + 1) j ⟨1, h_1_lt_pow_j⟩ ⟨1, h_1_lt_pow_dj⟩ = 0 := by
    change chiP (1 * 7 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (1 * 7 ^ 1 + 1 + 1 : ℕ) = 9 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_9]
  have hE12 : unfolding 7 (j + 1) j ⟨1, h_1_lt_pow_j⟩ ⟨3, h_3_lt_pow_dj⟩ = 1 := by
    change chiP (1 * 7 ^ ((j + 1) - j) + 3 + 1) = 1
    rw [h_sub]
    have h_eq : (1 * 7 ^ 1 + 3 + 1 : ℕ) = 11 := by norm_num
    rw [h_eq]
    exact chiP_eleven_eq_one
  have hE13 : unfolding 7 (j + 1) j ⟨1, h_1_lt_pow_j⟩ ⟨5, h_5_lt_pow_dj⟩ = 1 := by
    change chiP (1 * 7 ^ ((j + 1) - j) + 5 + 1) = 1
    rw [h_sub]
    have h_eq : (1 * 7 ^ 1 + 5 + 1 : ℕ) = 13 := by norm_num
    rw [h_eq]
    exact chiP_thirteen_eq_one
  have hE14 : unfolding 7 (j + 1) j ⟨1, h_1_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (1 * 7 ^ ((j + 1) - j) + 0 + 1) = 0
    rw [h_sub]
    have h_eq : (1 * 7 ^ 1 + 0 + 1 : ℕ) = 8 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_8]
  have hE15 : unfolding 7 (j + 1) j ⟨1, h_1_lt_pow_j⟩ ⟨2, h_2_lt_pow_dj⟩ = 0 := by
    change chiP (1 * 7 ^ ((j + 1) - j) + 2 + 1) = 0
    rw [h_sub]
    have h_eq : (1 * 7 ^ 1 + 2 + 1 : ℕ) = 10 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_10]
  have hE16 : unfolding 7 (j + 1) j ⟨1, h_1_lt_pow_j⟩ ⟨4, h_4_lt_pow_dj⟩ = 0 := by
    change chiP (1 * 7 ^ ((j + 1) - j) + 4 + 1) = 0
    rw [h_sub]
    have h_eq : (1 * 7 ^ 1 + 4 + 1 : ℕ) = 12 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_12]
  -- Row 2 (rho=3): n = 3*7 + sigma[k] + 1 = sigma[k] + 22.
  have hE20 : unfolding 7 (j + 1) j ⟨3, h_3_lt_pow_j⟩ ⟨6, h_6_lt_pow_dj⟩ = 0 := by
    change chiP (3 * 7 ^ ((j + 1) - j) + 6 + 1) = 0
    rw [h_sub]
    have h_eq : (3 * 7 ^ 1 + 6 + 1 : ℕ) = 28 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_28]
  have hE21 : unfolding 7 (j + 1) j ⟨3, h_3_lt_pow_j⟩ ⟨1, h_1_lt_pow_dj⟩ = 1 := by
    change chiP (3 * 7 ^ ((j + 1) - j) + 1 + 1) = 1
    rw [h_sub]
    have h_eq : (3 * 7 ^ 1 + 1 + 1 : ℕ) = 23 := by norm_num
    rw [h_eq]
    exact chiP_twenty_three_eq_one
  have hE22 : unfolding 7 (j + 1) j ⟨3, h_3_lt_pow_j⟩ ⟨3, h_3_lt_pow_dj⟩ = 0 := by
    change chiP (3 * 7 ^ ((j + 1) - j) + 3 + 1) = 0
    rw [h_sub]
    have h_eq : (3 * 7 ^ 1 + 3 + 1 : ℕ) = 25 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_25]
  have hE23 : unfolding 7 (j + 1) j ⟨3, h_3_lt_pow_j⟩ ⟨5, h_5_lt_pow_dj⟩ = 0 := by
    change chiP (3 * 7 ^ ((j + 1) - j) + 5 + 1) = 0
    rw [h_sub]
    have h_eq : (3 * 7 ^ 1 + 5 + 1 : ℕ) = 27 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_27]
  have hE24 : unfolding 7 (j + 1) j ⟨3, h_3_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (3 * 7 ^ ((j + 1) - j) + 0 + 1) = 0
    rw [h_sub]
    have h_eq : (3 * 7 ^ 1 + 0 + 1 : ℕ) = 22 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_22]
  have hE25 : unfolding 7 (j + 1) j ⟨3, h_3_lt_pow_j⟩ ⟨2, h_2_lt_pow_dj⟩ = 0 := by
    change chiP (3 * 7 ^ ((j + 1) - j) + 2 + 1) = 0
    rw [h_sub]
    have h_eq : (3 * 7 ^ 1 + 2 + 1 : ℕ) = 24 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_24]
  have hE26 : unfolding 7 (j + 1) j ⟨3, h_3_lt_pow_j⟩ ⟨4, h_4_lt_pow_dj⟩ = 0 := by
    change chiP (3 * 7 ^ ((j + 1) - j) + 4 + 1) = 0
    rw [h_sub]
    have h_eq : (3 * 7 ^ 1 + 4 + 1 : ℕ) = 26 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_26]
  -- Row 3 (rho=5): n = 5*7 + sigma[k] + 1 = sigma[k] + 36.
  have hE30 : unfolding 7 (j + 1) j ⟨5, h_5_lt_pow_j⟩ ⟨6, h_6_lt_pow_dj⟩ = 0 := by
    change chiP (5 * 7 ^ ((j + 1) - j) + 6 + 1) = 0
    rw [h_sub]
    have h_eq : (5 * 7 ^ 1 + 6 + 1 : ℕ) = 42 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_42]
  have hE31 : unfolding 7 (j + 1) j ⟨5, h_5_lt_pow_j⟩ ⟨1, h_1_lt_pow_dj⟩ = 1 := by
    change chiP (5 * 7 ^ ((j + 1) - j) + 1 + 1) = 1
    rw [h_sub]
    have h_eq : (5 * 7 ^ 1 + 1 + 1 : ℕ) = 37 := by norm_num
    rw [h_eq]
    exact chiP_thirty_seven_eq_one
  have hE32 : unfolding 7 (j + 1) j ⟨5, h_5_lt_pow_j⟩ ⟨3, h_3_lt_pow_dj⟩ = 0 := by
    change chiP (5 * 7 ^ ((j + 1) - j) + 3 + 1) = 0
    rw [h_sub]
    have h_eq : (5 * 7 ^ 1 + 3 + 1 : ℕ) = 39 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_39]
  have hE33 : unfolding 7 (j + 1) j ⟨5, h_5_lt_pow_j⟩ ⟨5, h_5_lt_pow_dj⟩ = 1 := by
    change chiP (5 * 7 ^ ((j + 1) - j) + 5 + 1) = 1
    rw [h_sub]
    have h_eq : (5 * 7 ^ 1 + 5 + 1 : ℕ) = 41 := by norm_num
    rw [h_eq]
    exact chiP_forty_one_eq_one
  have hE34 : unfolding 7 (j + 1) j ⟨5, h_5_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (5 * 7 ^ ((j + 1) - j) + 0 + 1) = 0
    rw [h_sub]
    have h_eq : (5 * 7 ^ 1 + 0 + 1 : ℕ) = 36 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_36]
  have hE35 : unfolding 7 (j + 1) j ⟨5, h_5_lt_pow_j⟩ ⟨2, h_2_lt_pow_dj⟩ = 0 := by
    change chiP (5 * 7 ^ ((j + 1) - j) + 2 + 1) = 0
    rw [h_sub]
    have h_eq : (5 * 7 ^ 1 + 2 + 1 : ℕ) = 38 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_38]
  have hE36 : unfolding 7 (j + 1) j ⟨5, h_5_lt_pow_j⟩ ⟨4, h_4_lt_pow_dj⟩ = 0 := by
    change chiP (5 * 7 ^ ((j + 1) - j) + 4 + 1) = 0
    rw [h_sub]
    have h_eq : (5 * 7 ^ 1 + 4 + 1 : ℕ) = 40 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_40]
  -- Row 4 (rho=2): n = 2*7 + sigma[k] + 1 = sigma[k] + 15.
  have hE40 : unfolding 7 (j + 1) j ⟨2, h_2_lt_pow_j⟩ ⟨6, h_6_lt_pow_dj⟩ = 0 := by
    change chiP (2 * 7 ^ ((j + 1) - j) + 6 + 1) = 0
    rw [h_sub]
    have h_eq : (2 * 7 ^ 1 + 6 + 1 : ℕ) = 21 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_21]
  have hE41 : unfolding 7 (j + 1) j ⟨2, h_2_lt_pow_j⟩ ⟨1, h_1_lt_pow_dj⟩ = 0 := by
    change chiP (2 * 7 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (2 * 7 ^ 1 + 1 + 1 : ℕ) = 16 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_16]
  have hE42 : unfolding 7 (j + 1) j ⟨2, h_2_lt_pow_j⟩ ⟨3, h_3_lt_pow_dj⟩ = 0 := by
    change chiP (2 * 7 ^ ((j + 1) - j) + 3 + 1) = 0
    rw [h_sub]
    have h_eq : (2 * 7 ^ 1 + 3 + 1 : ℕ) = 18 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_18]
  have hE43 : unfolding 7 (j + 1) j ⟨2, h_2_lt_pow_j⟩ ⟨5, h_5_lt_pow_dj⟩ = 0 := by
    change chiP (2 * 7 ^ ((j + 1) - j) + 5 + 1) = 0
    rw [h_sub]
    have h_eq : (2 * 7 ^ 1 + 5 + 1 : ℕ) = 20 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_20]
  have hE44 : unfolding 7 (j + 1) j ⟨2, h_2_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 0 := by
    change chiP (2 * 7 ^ ((j + 1) - j) + 0 + 1) = 0
    rw [h_sub]
    have h_eq : (2 * 7 ^ 1 + 0 + 1 : ℕ) = 15 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_15]
  have hE45 : unfolding 7 (j + 1) j ⟨2, h_2_lt_pow_j⟩ ⟨2, h_2_lt_pow_dj⟩ = 1 := by
    change chiP (2 * 7 ^ ((j + 1) - j) + 2 + 1) = 1
    rw [h_sub]
    have h_eq : (2 * 7 ^ 1 + 2 + 1 : ℕ) = 17 := by norm_num
    rw [h_eq]
    exact chiP_seventeen_eq_one
  have hE46 : unfolding 7 (j + 1) j ⟨2, h_2_lt_pow_j⟩ ⟨4, h_4_lt_pow_dj⟩ = 1 := by
    change chiP (2 * 7 ^ ((j + 1) - j) + 4 + 1) = 1
    rw [h_sub]
    have h_eq : (2 * 7 ^ 1 + 4 + 1 : ℕ) = 19 := by norm_num
    rw [h_eq]
    exact chiP_nineteen_eq_one
  -- Row 5 (rho=4): n = 4*7 + sigma[k] + 1 = sigma[k] + 29.
  have hE50 : unfolding 7 (j + 1) j ⟨4, h_4_lt_pow_j⟩ ⟨6, h_6_lt_pow_dj⟩ = 0 := by
    change chiP (4 * 7 ^ ((j + 1) - j) + 6 + 1) = 0
    rw [h_sub]
    have h_eq : (4 * 7 ^ 1 + 6 + 1 : ℕ) = 35 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_35]
  have hE51 : unfolding 7 (j + 1) j ⟨4, h_4_lt_pow_j⟩ ⟨1, h_1_lt_pow_dj⟩ = 0 := by
    change chiP (4 * 7 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (4 * 7 ^ 1 + 1 + 1 : ℕ) = 30 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_30]
  have hE52 : unfolding 7 (j + 1) j ⟨4, h_4_lt_pow_j⟩ ⟨3, h_3_lt_pow_dj⟩ = 0 := by
    change chiP (4 * 7 ^ ((j + 1) - j) + 3 + 1) = 0
    rw [h_sub]
    have h_eq : (4 * 7 ^ 1 + 3 + 1 : ℕ) = 32 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_32]
  have hE53 : unfolding 7 (j + 1) j ⟨4, h_4_lt_pow_j⟩ ⟨5, h_5_lt_pow_dj⟩ = 0 := by
    change chiP (4 * 7 ^ ((j + 1) - j) + 5 + 1) = 0
    rw [h_sub]
    have h_eq : (4 * 7 ^ 1 + 5 + 1 : ℕ) = 34 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_34]
  have hE54 : unfolding 7 (j + 1) j ⟨4, h_4_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 1 := by
    change chiP (4 * 7 ^ ((j + 1) - j) + 0 + 1) = 1
    rw [h_sub]
    have h_eq : (4 * 7 ^ 1 + 0 + 1 : ℕ) = 29 := by norm_num
    rw [h_eq]
    exact chiP_twenty_nine_eq_one
  have hE55 : unfolding 7 (j + 1) j ⟨4, h_4_lt_pow_j⟩ ⟨2, h_2_lt_pow_dj⟩ = 1 := by
    change chiP (4 * 7 ^ ((j + 1) - j) + 2 + 1) = 1
    rw [h_sub]
    have h_eq : (4 * 7 ^ 1 + 2 + 1 : ℕ) = 31 := by norm_num
    rw [h_eq]
    exact chiP_thirty_one_eq_one
  have hE56 : unfolding 7 (j + 1) j ⟨4, h_4_lt_pow_j⟩ ⟨4, h_4_lt_pow_dj⟩ = 0 := by
    change chiP (4 * 7 ^ ((j + 1) - j) + 4 + 1) = 0
    rw [h_sub]
    have h_eq : (4 * 7 ^ 1 + 4 + 1 : ℕ) = 33 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_33]
  -- Row 6 (rho=6): n = 6*7 + sigma[k] + 1 = sigma[k] + 43.
  have hE60 : unfolding 7 (j + 1) j ⟨6, h_6_lt_pow_j⟩ ⟨6, h_6_lt_pow_dj⟩ = 0 := by
    change chiP (6 * 7 ^ ((j + 1) - j) + 6 + 1) = 0
    rw [h_sub]
    have h_eq : (6 * 7 ^ 1 + 6 + 1 : ℕ) = 49 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_49]
  have hE61 : unfolding 7 (j + 1) j ⟨6, h_6_lt_pow_j⟩ ⟨1, h_1_lt_pow_dj⟩ = 0 := by
    change chiP (6 * 7 ^ ((j + 1) - j) + 1 + 1) = 0
    rw [h_sub]
    have h_eq : (6 * 7 ^ 1 + 1 + 1 : ℕ) = 44 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_44]
  have hE62 : unfolding 7 (j + 1) j ⟨6, h_6_lt_pow_j⟩ ⟨3, h_3_lt_pow_dj⟩ = 0 := by
    change chiP (6 * 7 ^ ((j + 1) - j) + 3 + 1) = 0
    rw [h_sub]
    have h_eq : (6 * 7 ^ 1 + 3 + 1 : ℕ) = 46 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_46]
  have hE63 : unfolding 7 (j + 1) j ⟨6, h_6_lt_pow_j⟩ ⟨5, h_5_lt_pow_dj⟩ = 0 := by
    change chiP (6 * 7 ^ ((j + 1) - j) + 5 + 1) = 0
    rw [h_sub]
    have h_eq : (6 * 7 ^ 1 + 5 + 1 : ℕ) = 48 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_48]
  have hE64 : unfolding 7 (j + 1) j ⟨6, h_6_lt_pow_j⟩ ⟨0, h_zero_lt_pow_dj⟩ = 1 := by
    change chiP (6 * 7 ^ ((j + 1) - j) + 0 + 1) = 1
    rw [h_sub]
    have h_eq : (6 * 7 ^ 1 + 0 + 1 : ℕ) = 43 := by norm_num
    rw [h_eq]
    exact chiP_forty_three_eq_one
  have hE65 : unfolding 7 (j + 1) j ⟨6, h_6_lt_pow_j⟩ ⟨2, h_2_lt_pow_dj⟩ = 0 := by
    change chiP (6 * 7 ^ ((j + 1) - j) + 2 + 1) = 0
    rw [h_sub]
    have h_eq : (6 * 7 ^ 1 + 2 + 1 : ℕ) = 45 := by norm_num
    rw [h_eq]
    simp [chiP, h_not_prime_45]
  have hE66 : unfolding 7 (j + 1) j ⟨6, h_6_lt_pow_j⟩ ⟨4, h_4_lt_pow_dj⟩ = 1 := by
    change chiP (6 * 7 ^ ((j + 1) - j) + 4 + 1) = 1
    rw [h_sub]
    have h_eq : (6 * 7 ^ 1 + 4 + 1 : ℕ) = 47 := by norm_num
    rw [h_eq]
    exact chiP_forty_seven_eq_one
  -- Fin 7 .val and ≠ helpers for the if-then-else ρ, σ.
  have hv0 : (0 : Fin 7).val = 0 := rfl
  have hv1 : (1 : Fin 7).val = 1 := rfl
  have hv2 : (2 : Fin 7).val = 2 := rfl
  have hv3 : (3 : Fin 7).val = 3 := rfl
  have hv4 : (4 : Fin 7).val = 4 := rfl
  have hv5 : (5 : Fin 7).val = 5 := rfl
  have hv6 : (6 : Fin 7).val = 6 := rfl
  have hne_1_0 : (1 : ℕ) ≠ 0 := by decide
  have hne_2_0 : (2 : ℕ) ≠ 0 := by decide
  have hne_2_1 : (2 : ℕ) ≠ 1 := by decide
  have hne_3_0 : (3 : ℕ) ≠ 0 := by decide
  have hne_3_1 : (3 : ℕ) ≠ 1 := by decide
  have hne_3_2 : (3 : ℕ) ≠ 2 := by decide
  have hne_4_0 : (4 : ℕ) ≠ 0 := by decide
  have hne_4_1 : (4 : ℕ) ≠ 1 := by decide
  have hne_4_2 : (4 : ℕ) ≠ 2 := by decide
  have hne_4_3 : (4 : ℕ) ≠ 3 := by decide
  have hne_5_0 : (5 : ℕ) ≠ 0 := by decide
  have hne_5_1 : (5 : ℕ) ≠ 1 := by decide
  have hne_5_2 : (5 : ℕ) ≠ 2 := by decide
  have hne_5_3 : (5 : ℕ) ≠ 3 := by decide
  have hne_5_4 : (5 : ℕ) ≠ 4 := by decide
  have hne_6_0 : (6 : ℕ) ≠ 0 := by decide
  have hne_6_1 : (6 : ℕ) ≠ 1 := by decide
  have hne_6_2 : (6 : ℕ) ≠ 2 := by decide
  have hne_6_3 : (6 : ℕ) ≠ 3 := by decide
  have hne_6_4 : (6 : ℕ) ≠ 4 := by decide
  have hne_6_5 : (6 : ℕ) ≠ 5 := by decide
  -- Set up the named ρ, σ, M for the rest of the proof.
  set Mρ : Fin 7 → Fin (7 ^ j) :=
    fun i => if i.val = 0 then ⟨0, h_pow_j_pos⟩
             else if i.val = 1 then ⟨1, h_1_lt_pow_j⟩
             else if i.val = 2 then ⟨3, h_3_lt_pow_j⟩
             else if i.val = 3 then ⟨5, h_5_lt_pow_j⟩
             else if i.val = 4 then ⟨2, h_2_lt_pow_j⟩
             else if i.val = 5 then ⟨4, h_4_lt_pow_j⟩
             else ⟨6, h_6_lt_pow_j⟩ with hMρ_def
  set Mσ : Fin 7 → Fin (7 ^ ((j + 1) - j)) :=
    fun i => if i.val = 0 then ⟨6, h_6_lt_pow_dj⟩
             else if i.val = 1 then ⟨1, h_1_lt_pow_dj⟩
             else if i.val = 2 then ⟨3, h_3_lt_pow_dj⟩
             else if i.val = 3 then ⟨5, h_5_lt_pow_dj⟩
             else if i.val = 4 then ⟨0, h_zero_lt_pow_dj⟩
             else if i.val = 5 then ⟨2, h_2_lt_pow_dj⟩
             else ⟨4, h_4_lt_pow_dj⟩ with hMσ_def
  set M : Matrix (Fin 7) (Fin 7) ℚ :=
    (unfolding 7 (j + 1) j).submatrix Mρ Mσ with hM_def
  -- Explicit 7x7 matrix matching the entries.
  let Mexp : Matrix (Fin 7) (Fin 7) ℚ :=
    !![1, 1, 0, 0, 0, 1, 1;
       0, 0, 1, 1, 0, 0, 0;
       0, 1, 0, 0, 0, 0, 0;
       0, 1, 0, 1, 0, 0, 0;
       0, 0, 0, 0, 0, 1, 1;
       0, 0, 0, 0, 1, 1, 0;
       0, 0, 0, 0, 1, 0, 1]
  -- Prove M = Mexp.
  have hMeq : M = Mexp := by
    ext i k
    simp only [hM_def, Matrix.submatrix_apply, hMρ_def, hMσ_def]
    fin_cases i <;> fin_cases k <;>
      simp only [hv0, hv1, hv2, hv3, hv4, hv5, hv6,
                 if_true, if_false, Nat.one_ne_zero,
                 hne_1_0, hne_2_0, hne_2_1, hne_3_0, hne_3_1, hne_3_2,
                 hne_4_0, hne_4_1, hne_4_2, hne_4_3,
                 hne_5_0, hne_5_1, hne_5_2, hne_5_3, hne_5_4,
                 hne_6_0, hne_6_1, hne_6_2, hne_6_3, hne_6_4, hne_6_5]
    all_goals first
      | exact hE00
      | exact hE01
      | exact hE02
      | exact hE03
      | exact hE04
      | exact hE05
      | exact hE06
      | exact hE10
      | exact hE11
      | exact hE12
      | exact hE13
      | exact hE14
      | exact hE15
      | exact hE16
      | exact hE20
      | exact hE21
      | exact hE22
      | exact hE23
      | exact hE24
      | exact hE25
      | exact hE26
      | exact hE30
      | exact hE31
      | exact hE32
      | exact hE33
      | exact hE34
      | exact hE35
      | exact hE36
      | exact hE40
      | exact hE41
      | exact hE42
      | exact hE43
      | exact hE44
      | exact hE45
      | exact hE46
      | exact hE50
      | exact hE51
      | exact hE52
      | exact hE53
      | exact hE54
      | exact hE55
      | exact hE56
      | exact hE60
      | exact hE61
      | exact hE62
      | exact hE63
      | exact hE64
      | exact hE65
      | exact hE66
      | (simp [Mexp, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons, Matrix.cons_val_fin_one])
  rw [hMeq]
  -- Now show Mexp.det ≠ 0 via two nested fromBlocks decompositions:
  --   Layer 1: Mexp = fromBlocks A B 0 D (1 + 6 split, lower-left zero).
  --   Layer 2: D   = fromBlocks D1 0 0 D2 (3 + 3 split, both off-diagonals zero).
  -- The 1×1 block A has det 1 (`det_fin_one`); the two 3×3 blocks D1, D2
  -- have det -1 and -2 (`det_fin_three`); total det = 1 * ((-1) * (-2)) = 2.
  -- Layer-1 blocks: 1+6 split.
  let A : Matrix (Fin 1) (Fin 1) ℚ := !![1]
  let B : Matrix (Fin 1) (Fin 6) ℚ := !![1, 0, 0, 0, 1, 1]
  let D : Matrix (Fin 6) (Fin 6) ℚ :=
    !![0, 1, 1, 0, 0, 0;
       1, 0, 0, 0, 0, 0;
       1, 0, 1, 0, 0, 0;
       0, 0, 0, 0, 1, 1;
       0, 0, 0, 1, 1, 0;
       0, 0, 0, 1, 0, 1]
  -- Layer-2 blocks: 3+3 split of D.
  let D1 : Matrix (Fin 3) (Fin 3) ℚ := !![0, 1, 1; 1, 0, 0; 1, 0, 1]
  let D2 : Matrix (Fin 3) (Fin 3) ℚ := !![0, 1, 1; 1, 1, 0; 1, 0, 1]
  -- Layer 1 equality: under finSumFinEquiv : Fin 1 ⊕ Fin 6 ≃ Fin 7.
  have h_fromBlocks :
      Mexp.submatrix finSumFinEquiv finSumFinEquiv = Matrix.fromBlocks A B 0 D := by
    ext i j
    rcases i with i | i <;> rcases j with j | j
    all_goals fin_cases i <;> fin_cases j <;> rfl
  -- Layer 2 equality: under finSumFinEquiv : Fin 3 ⊕ Fin 3 ≃ Fin 6.
  have h_fromBlocks_D :
      D.submatrix finSumFinEquiv finSumFinEquiv =
        Matrix.fromBlocks D1 (0 : Matrix (Fin 3) (Fin 3) ℚ) 0 D2 := by
    ext i j
    rcases i with i | i <;> rcases j with j | j
    all_goals fin_cases i <;> fin_cases j <;> rfl
  have h_Mexp_det :
      Mexp.det = (Matrix.fromBlocks A B (0 : Matrix (Fin 6) (Fin 1) ℚ) D).det := by
    rw [← h_fromBlocks]
    simp [Matrix.det_submatrix_equiv_self]
  rw [h_Mexp_det, Matrix.det_fromBlocks_zero₂₁]
  -- A.det = 1 (1×1 case).
  have hA_det : A.det = 1 := by
    show (!![(1 : ℚ)] : Matrix (Fin 1) (Fin 1) ℚ).det = 1
    simp [Matrix.det_fin_one]
  -- D.det = 2 via the layer-2 split: D1.det * D2.det = (-1) * (-2) = 2.
  have hD_det : D.det = 2 := by
    have h_D_det_eq :
        D.det =
          (Matrix.fromBlocks D1 (0 : Matrix (Fin 3) (Fin 3) ℚ) 0 D2).det := by
      rw [← h_fromBlocks_D]
      simp [Matrix.det_submatrix_equiv_self]
    rw [h_D_det_eq, Matrix.det_fromBlocks_zero₂₁]
    have hD1 : D1.det = -1 := by
      show (!![(0 : ℚ), 1, 1; 1, 0, 0; 1, 0, 1] :
              Matrix (Fin 3) (Fin 3) ℚ).det = -1
      simp [Matrix.det_fin_three]
    have hD2 : D2.det = -2 := by
      show (!![(0 : ℚ), 1, 1; 1, 1, 0; 1, 0, 1] :
              Matrix (Fin 3) (Fin 3) ℚ).det = -2
      simp [Matrix.det_fin_three]
      norm_num
    rw [hD1, hD2]
    norm_num
  rw [hA_det, hD_det]
  -- 1 * 2 = 2, IsUnit (2 : ℚ) since (2 : ℚ) ≠ 0.
  have h_ne_zero : (1 * 2 : ℚ) ≠ 0 := by norm_num
  exact h_ne_zero.isUnit

/--
**Corner-case main theorem (W = 7, d = j + 1).** The unfolding rank is
exactly `7` for every `j ≥ 1`. Mirrors `mps_bond_dim_W_eq_9_d_eq_j_plus_1`
(both citing the general `upper_bound`, since `rank_le_width` only gives
`rank ≤ 7`, but `R = φ(7) + 1 = 7` and the slab has 7 cols, so they
coincide here). Fully formalised in Lean (no `sorry`, no `axiom`).

**Upper-bound subtlety:** the matrix has `7^((j+1)-j) = 7` columns, so
`rank_le_width` gives `rank ≤ 7`, which equals the target. We still cite
the general `upper_bound` for uniformity with the other corner cases.

**Proof technique novelty.** Twelfth corner closure; second use of
`Matrix.det_fromBlocks_zero₂₁` (after S152's W=9). W=7 was on the S128
"block-triangular-required" list; S159's pre-search found the same
1+(3+3) decomposition shape as W=9 but with `D₂.det = -2` rather than
`-1`. This is the **first instance with total `det = ±2`** (W=9 had
total `det = 1`), demonstrating that the BlockTriangular template is
robust to nontrivial determinant magnitudes — `IsUnit` over `ℚ` follows
from `Ne.isUnit` since `(2 : ℚ) ≠ 0`.

**Helper reuse.** S159 added zero new `chiP_X_eq_one` helpers — every
prime in the 7×7 submatrix (`{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
37, 41, 43, 47}`) was already declared. Ten new `h_not_prime_X` lemmas
internal to the proof (for composites `{6, 9, 15, 18, 24, 27, 33, 36,
42, 45}` not seen in the W=9 closure).

Together with prior corners, this is the **twelfth** unconditional
instance of `mps_bond_dim` and the **eleventh** instance over a wheel
`W ≥ 3`.
-/
theorem mps_bond_dim_W_eq_7_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    (unfolding 7 (j + 1) j).rank = 7 := by
  apply Nat.le_antisymm
  · -- `≤ 7`: from `upper_bound`. (`rank_le_width` also gives `≤ 7` here,
    -- but we route through `upper_bound` for uniformity.)
    have hW : (2 : ℕ) ≤ 7 := by norm_num
    have hj_hi : j < j + 1 := Nat.lt_succ_self j
    have h := upper_bound 7 (j + 1) j hW hj hj_hi
    have h_eq : Nat.totient 7 * 7 ^ ((j + 1 : ℕ) - j - 1) + 1 = 7 := by
      have h_sub : (j + 1 : ℕ) - j - 1 = 0 := by omega
      rw [h_sub]
      decide
    linarith
  · -- `7 ≤ rank`: from the corner-case prime exhibit.
    obtain ⟨ρ, σ, hUnit⟩ := exists_invertible_submatrix_W_eq_7_d_eq_j_plus_1 j hj
    have h_eq : ((unfolding 7 (j + 1) j).submatrix ρ σ).rank = 7 := by
      have h := Matrix.rank_of_isUnit ((unfolding 7 (j + 1) j).submatrix ρ σ) hUnit
      simpa using h
    calc (7 : ℕ) = ((unfolding 7 (j + 1) j).submatrix ρ σ).rank := h_eq.symm
      _ ≤ (unfolding 7 (j + 1) j).rank := Matrix.rank_submatrix_le _ _ _

end E2_1
