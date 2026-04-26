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

/--
**Main theorem (E2.1).** For every `W ≥ 2` and every cut `1 ≤ j < d`,
the rank of the `j`-th unfolding of the base-`W` prime indicator equals
`min (W^j) (φ(W) * W^(d-j-1) + 1)`.

This is the closed form proved (informally) in
`novel/mps_bond_dimension.md` and saturated empirically for
`W ∈ {2, 6, 30, 210}` and `d` up to 20.
-/
theorem mps_bond_dim
    (W d j : ℕ) (hW : 2 ≤ W) (hj_lo : 1 ≤ j) (hj_hi : j < d) :
    (unfolding W d j).rank =
      min (W ^ j) (Nat.totient W * W ^ (d - j - 1) + 1) := by
  sorry

/-!
## Decomposition of the proof

The informal proof in `novel/mps_bond_dimension.md` splits into two
inequalities, separated as `upper_bound` and `lower_bound` below. The
main theorem then follows by `Nat.le_antisymm` plus a small case split on
which of the two arguments to `min` is the active one.
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
  -- Step 1: the predicate `n ↦ gcd(n+1, W) = 1` is W-periodic.
  have hP_periodic :
      Function.Periodic (fun n => Nat.gcd (n + 1) W = 1) W := fun x => by
    show (Nat.gcd (x + W + 1) W = 1) = (Nat.gcd (x + 1) W = 1)
    rw [show x + W + 1 = (x + 1) + W from by ring, Nat.gcd_add_self_left]
  -- Step 2: convert the Fin-indexed cardinality to `Nat.count`.
  have step_fin_to_count :
      ((Finset.univ : Finset (Fin (W ^ (d - j)))).filter
          (fun k => Nat.gcd (k.val + 1) W = 1)).card =
        Nat.count (fun n => Nat.gcd (n + 1) W = 1) (W ^ (d - j)) := by
    rw [Nat.count_eq_card_filter_range, Finset.card_eq_sum_ones,
        Finset.card_eq_sum_ones, Finset.sum_filter, Finset.sum_filter,
        Fin.sum_univ_eq_sum_range]
  rw [step_fin_to_count]
  -- Step 3: split W^(d-j) = W * W^(d-j-1).
  have W_pow : W ^ (d - j) = W * W ^ (d - j - 1) := by
    have h : (d - j - 1) + 1 = d - j := Nat.sub_add_cancel hdj
    calc W ^ (d - j)
        = W ^ ((d - j - 1) + 1) := by rw [h]
      _ = W * W ^ (d - j - 1) := by rw [pow_succ, mul_comm]
  rw [W_pow]
  -- Step 4: multi-block periodicity. count P (W * m) = m * count P W.
  have block_count : ∀ m,
      Nat.count (fun n => Nat.gcd (n + 1) W = 1) (W * m) =
        m * Nat.count (fun n => Nat.gcd (n + 1) W = 1) W := by
    intro m
    induction m with
    | zero => simp
    | succ k ih =>
      have step_a : W * (k + 1) = W * k + W := by ring
      rw [step_a, Nat.count_add (W * k) W]
      -- The shifted predicate counts the same number of survivors as P over W
      -- by W-periodicity (W * k is a multiple of W).
      have shift :
          Nat.count (fun n => Nat.gcd (W * k + n + 1) W = 1) W =
            Nat.count (fun n => Nat.gcd (n + 1) W = 1) W := by
        rw [Nat.count_eq_card_filter_range, Nat.count_eq_card_filter_range]
        congr 1
        apply Finset.filter_congr
        intros n _
        rw [show W * k + n + 1 = n + k * W + 1 from by ring]
        exact Iff.of_eq (hP_periodic.nat_mul k n)
      rw [shift, ih]
      ring
  rw [block_count]
  -- Step 5: count over a single W-block equals totient W.
  have base_count :
      Nat.count (fun n => Nat.gcd (n + 1) W = 1) W = Nat.totient W := by
    rw [Nat.count_eq_card_filter_range, ← Nat.filter_coprime_Ico_eq_totient W 1]
    refine Finset.card_bij' (fun n _ => n + 1) (fun m _ => m - 1) ?_ ?_ ?_ ?_
    · intro n hn
      rw [Finset.mem_filter, Finset.mem_range] at hn
      rw [Finset.mem_filter, Finset.mem_Ico]
      refine ⟨⟨Nat.succ_le_succ (Nat.zero_le _), by omega⟩, ?_⟩
      show Nat.Coprime W (n + 1)
      rw [Nat.Coprime, Nat.gcd_comm]; exact hn.2
    · intro m hm
      rw [Finset.mem_filter, Finset.mem_Ico] at hm
      rw [Finset.mem_filter, Finset.mem_range]
      have hm_pos : 1 ≤ m := hm.1.1
      refine ⟨by omega, ?_⟩
      have m_eq : m - 1 + 1 = m := Nat.sub_add_cancel hm_pos
      rw [m_eq]
      show Nat.gcd m W = 1
      rw [Nat.gcd_comm]; exact hm.2
    · intros a _; omega
    · intro m hm
      rw [Finset.mem_filter, Finset.mem_Ico] at hm
      have : 1 ≤ m := hm.1.1
      omega
  rw [base_count]
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
  sorry

/--
Lower bound. Empirically saturated at every measured `(W, d, j)`. The
informal argument identifies, among the `φ(W) * W^(d-j-1)` "live" columns,
enough rows whose restriction is linearly independent over `ℚ` —
exhibited via the prime-counting density of base-`W` blocks. The full
formal argument is open.
-/
theorem lower_bound
    (W d j : ℕ) (hW : 2 ≤ W) (hj_lo : 1 ≤ j) (hj_hi : j < d) :
    min (W ^ j) (Nat.totient W * W ^ (d - j - 1) + 1) ≤
      (unfolding W d j).rank := by
  sorry

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

end E2_1
