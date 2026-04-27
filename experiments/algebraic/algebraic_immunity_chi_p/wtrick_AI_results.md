# wtrick_AI.py — Results

This is an auxiliary script supporting the main S92 experiment. It
computes algebraic immunity AI(chi_P_{W,b}) of the W-tricked prime
indicator `chi_P_{W,b}(n) := chi_P(W*n + b)` for `gcd(b, W) = 1`,
and the matched random matched-density Bernoulli baseline.

**Key result:** the AI deviation `AI(chi_P) = 2 << AI(random)`
(documented in `algebraic_immunity_chi_p_results.md`) is **fully
removed** by the W-trick at W >= 6. AI(chi_P_{W,b}) tracks
AI(random matched-density) within ±1 step, often exact match
(zero std) at multiple W ∈ {6, 30}, b ∈ {1, 5, 7, 11}, N=8..11.

The W-trick correction provides the closure: chi_P's polynomial-method
deviation from random IS the polynomial-method encoding of small-modulus
residue-class structure (the mod-4 sieve fact and its mod-q extensions).
Same picture as E2.13 (Gowers W-trick) and E2.14 (Anderson W-cascade).

See `algebraic_immunity_chi_p_results.md` for the full writeup and
`wtrick_AI_data.json` for raw numerical data.
