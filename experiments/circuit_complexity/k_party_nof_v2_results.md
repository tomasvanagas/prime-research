# k-Party NOF Communication Complexity v2: Results

**Script:** `k_party_nof_v2.py`
**Session:** 20
**See also:** `k_party_nof_results.md`

## What Was Tested
Refined analysis of k-party NOF communication complexity of pi(x). Investigated rank over GF(2) and GF(p), the normalized "information per party," and whether rank drops below full for the largest group. Tested N up to 18.

## Key Findings
- Verified the exact formula: max mode-unfolding rank = 2^{ceil(N/k)} for k >= 3
- Rank over GF(2): also full for all tested (N,k) pairs with k >= 3
- Rank over GF(p): full for p = 3, 5, 7 as well -- no algebraic structure at any characteristic
- Normalized information per party = log2(max_rank) / (N/k) = 1.0 exactly (full information)
- For k = N/c parties, NOF complexity = c (constant), but this requires N/c parties -- not useful for circuits

## Verdict
**CLOSED** -- subsumed by `k_party_nof_results.md`.
**Failure Mode:** Information loss (all mode-unfolding ranks are full; no multiparty structure exploitable by ACC circuits)

## One-Line Summary
k-party NOF refined analysis confirms full rank over all fields; no algebraic structure at any characteristic.
