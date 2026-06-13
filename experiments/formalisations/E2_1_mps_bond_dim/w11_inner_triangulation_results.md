# w11_inner_triangulation — results

**Result (S206):** Leading-row upper-triangulation of the inner 10×10
(rows ≥ 1, live cols of W=11) is OBSTRUCTED in rows [1, 11) but FOUND
in rows [1, 22) with `ρ ↦ (3, 2, 1, 6, 5, 7, 16, 10, 19, 18)`,
`σ ↦ (9, 8, 7, 6, 3, 5, 4, 2, 1, 0)` (max row 19 → forces j ≥ 2).

Diagonal primes `{19, 31, 43, 59, 73, 83, 113, 181, 199, 211}`. Six new
prime helpers needed: `{67, 71, 73, 79, 83, 113, 181}`.

**Falsification statement:** if the matrix `M[r][c] = chiP(11r + c + 1)`
admitted a leading-row triangulation over `r ∈ [1, 11)` and live cols,
the proposed Path A multi-session split would not be needed. Pre-search
confirms the obstruction.

See `w11_blocktriangular_search_results.md` for the consolidated S206
findings (parity-block decomposition, atomicity of odd 5×5, three
multi-session paths forward).
