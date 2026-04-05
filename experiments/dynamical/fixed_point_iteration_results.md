# Fixed-Point Iteration Approach: Results

**Script:** fixed_point_iteration.py

## What Was Tested
Fixed-point iteration methods to compute p(n): (1) prime walk guided by R(x) as direction indicator, (2) corrected Newton with bias estimation for R(x) vs pi(x), (3) iterative snap-to-prime refinement, (4) ensemble voting with multiple R_inv shifts. Tested on n=10..1000.

## Key Findings
- g5 prime walk: ~70-80% correct; R(x) gives correct direction most of the time but misdirects ~20-30%
- g6 corrected Newton: bias R(p(n)) - n has mean ~0.55, std ~1.3; using R_inv(n + 0.55) gives marginal improvement
- Iterative refinement: after 5 iterations, ~75-85% correct (diminishing returns after iteration 2)
- Ensemble voting (5 shifted R_inv values): ~80-85% correct, +5-10% improvement over single R_inv
- Best single approach: nearest_prime(R_inv(n + c)) with optimal c ~ 0.5 gives ~82% exact
- All methods fail on the ~15-20% of cases where R_inv(n) lands on the wrong side of the prime gap
- Performance degrades with n as gaps grow relative to R_inv precision

## Verdict
**CLOSED**
**Failure Mode:** I (Information Loss -- R_inv gives only the smooth approximation; the ~15-20% failure rate corresponds to cases where delta(n) exceeds half the local gap)

## One-Line Summary
Fixed-point/snap-to-prime iteration achieves ~80-85% exact via ensemble voting but cannot reach 100% without computing pi(x) exactly.
