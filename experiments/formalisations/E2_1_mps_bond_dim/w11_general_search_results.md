# w11_general_search — results

**Result (S206):** Exhaustive search across 23 ordered partitions of 11
with block sizes ≤ 5 over rows [0, 11) × all cols [0, 11).

* **Partitions with all blocks ≤ 3** (no det_fin_four/det_fin_five
  needed): **0 candidates** for every partition tried (3+3+3+2,
  2+3+3+3, 3+3+3+1+1, 1+3+3+3+1, 1+1+3+3+3, 3+3+1+1+3, 3+3+2+2+1,
  1+3+2+2+3). The W=11 11×11 admits NO block-triangular decomposition
  with all blocks ≤ 3 over rows [0, 11).
* **Partitions including a size-4 block** (4+3+3+1, 1+3+3+4, 3+4+3+1,
  1+4+3+3, 4+4+3, 4+3+4, 3+4+4, 4+4+2+1, 4+3+2+2): **0 candidates**.
* **Partitions including a size-5 block**: candidates exist for
  (5, 5, 1), (1, 5, 5), (5, 1, 5), (3, 5, 3), (3, 3, 5), (4, 5, 2),
  (2, 5, 4). **In every successful candidate the 5×5 block is the
  same**: rows `(1, 3, 5, 7, 9)` × cols `(1, 3, 5, 7, 9)` — the odd-row
  × odd-col parity block.

The triangulability check confirmed:
* **Top 5×5 of (5, 5, 1)**: rows `(0, 2, 4, 6, 8)` × cols
  `(0, 4, 6, 8, 10)` IS leading-row triangulable
  (`ρ ↦ (0, 6, 2, 8, 4)`, `σ ↦ (10, 4, 6, 0, 8)`, diagonal primes
  `{11, 71, 29, 89, 53}`).
* **Bottom 5×5 (odd parity block)**: NOT triangulable. Verified
  separately atomically irreducible by `w11_odd_block_atomicity.py`.

**Falsification statement:** if the W=11 11×11 admitted a block-
triangular decomposition with all blocks ≤ 3 over rows [0, 11), the
single-session arc-anticipated closure would proceed via composed
`det_fin_three`/`det_fin_two`/`det_fin_one` calls with no need for any
mathlib extension. The general search confirms the decomposition does
not exist.

**Implication:** Single-session closure with rows [0, 11) requires
either a `det_fin_five` capability or pivoting to a different W
(e.g., W = 14, composite, no parity-atomicity issue).

See `w11_blocktriangular_search_results.md` for context.
