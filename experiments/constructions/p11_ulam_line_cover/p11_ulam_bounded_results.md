# Results: p11_ulam_bounded.py

Bounded-direction numpy variant of the Ulam-spiral line cover.
Restricts the line search to canonical directions (a, b) with
gcd(|a|,|b|) = 1 and max(|a|, |b|) ≤ K. Used at large N (10⁵, 10⁶)
where all-pairs O(π(N)²) is infeasible. Gives an UPPER BOUND on the
true greedy line cover.

See `p11_ulam_line_cover_results.md` for the full empirical
findings, scaling table, and analysis.

Sanity check at N=10⁴, K=5: bounded gives L_primes = 95 vs exact
L = 91 (≤ 5% over). At K=20 the bound is tighter still.

Output CSVs:
- `summary_bounded_N100000_K20.csv` — N=10⁵ scaling
- `summary_bounded_N1000000_K20.csv` — N=10⁶ scaling
- `top_lines_bounded_*.csv` — top-K dominant lines at each N

## What would falsify

If the bounded estimate at moderately large K (e.g. K=50) gives a
materially smaller L than at K=20 (say, ≥10% reduction), then K=20
is missing important "knight's move" diagonals and the empirical
scaling needs revisiting. Slot 2/3 will sweep K to verify saturation.
