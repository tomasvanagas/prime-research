# extract_annihilator.py — Results

This is an auxiliary script supporting the main S92 experiment. It
extracts the explicit smallest-degree annihilator polynomial g of
chi_P over F_2 via row-reduction with nullspace recovery, and verifies
g * chi_P == 0 (or g * (1 + chi_P) == 0) over the full input space.

**Key result:** the annihilator is the SAME polynomial
`g(x) = (1 + x_0)(1 + x_1)` for every N from 5 through 13, encoding
"n ≡ 0 mod 4". Verified by direct evaluation across N=5..13.

See `algebraic_immunity_chi_p_results.md` for the full session writeup
and `extract_output.log` for raw extracted annihilators at each N.
