# H-T Signed Transfer v2: Results

**Date:** 2026-04-04 (Session 16)
**Script:** ht_signed_transfer_v2.py

## What Was Tested
Deep investigation of whether the Helfgott-Thompson O(x^{3/5}) technique for M(x) can transfer INDIRECTLY to pi(x). Six attack vectors tested: (1) explicit formula / oscillatory zero sum, (2) weighted prime counting with cancellation, (3) Buchstab + signed weights, (4) M(x) to pi(x) conversion cost, (5) concrete numerical experiments, (6) identity-based conversion via omega partial sums.

## Key Findings
- Approach 1 (zero sum): H-T works on sums over integers with multiplicative structure; the zeta zero sum is analytic -- fundamentally incompatible
- Approach 2 (weighted counting): No weighting scheme creates cancellation AND allows pi(x) recovery. Primes are a subset, any multiplicative weight collapses all primes to the same value
- Approach 3 (Buchstab + signs): The signed Buchstab version IS what H-T computes (M(x)), but undoing signs costs O(x^{2/3}) for A_3, A_5 terms
- Approach 4 (M->pi conversion): Costs at least O(x^{1/2+eps}) analytic or O(x^{2/3}) combinatorial -- conversion is the bottleneck
- Approach 5 (numerical): Verified M(x) and pi(x) are coupled but no cheap conversion path exists
- Approach 6 (omega partial sums): Identity exists but omega partial sums cost O(x^{2/3}), same as Meissel-Lehmer

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence) -- all six indirect routes reduce to known O(x^{2/3}) bottlenecks

## One-Line Summary
H-T O(x^{3/5}) for M(x) cannot transfer to pi(x) through any indirect route: all conversions reintroduce O(x^{2/3}) cost.
