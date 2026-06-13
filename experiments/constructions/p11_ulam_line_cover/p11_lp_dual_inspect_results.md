# p11_lp_dual_inspect — LP solution and direction-weight extraction

Helper script for slot 3 (S486). Solves the same set-cover LP as
`p11_lp_relaxation.py` and extracts:

- active LP variables `x_l > 1e-7` (lines with positive weight)
- total LP weight per direction `(a, b)` (sums `x_l` over intercept c)
- count of integer-1 lines vs purely-fractional lines

This was used to identify that the prime-Ulam LP solution at N=10⁴
puts 69% of weight on slope-±1 Hardy-Littlewood-rich diagonals
(direction (1, −1) carries 28.03; direction (1, +1) carries 25.74)
versus the random-baseline LP at N=10⁴ which puts 100% weight on
direction (1, 0) (axis-aligned vertical columns) at integer weight 1.0.

This is a **diagnostic tool**, not a primary measurement. Results
appear in `p11_lp_relaxation_results.md` "LP solution structure"
table and in the slot-3 session synthesis (S486).

## Falsifier

If running on a different `N` or `K` produces a qualitatively different
direction breakdown (e.g. weight no longer concentrating on slope-±1
for primes-Ulam), the slot-3 finding "primes have 69% LP weight on
HL diagonals" would be embedding-fragile. Tested at K=5, 10, 20 in
slot 3 with stable conclusion (LP value changes by < 1%, direction
breakdown unchanged).

## Edges cited

- Stanley 1989 (cross-domain matroid LP)
- HL Conjecture F (off-edge, explains the weight concentration on slope-±1)

## CLI examples

```
python p11_lp_dual_inspect.py --embedding ulam --N 10000 --K 5 --top 20
python p11_lp_dual_inspect.py --embedding ulam --N 10000 --K 5 --baseline-seed 42
python p11_lp_dual_inspect.py --embedding residue --q 210 --N 10000 --K 5
```
