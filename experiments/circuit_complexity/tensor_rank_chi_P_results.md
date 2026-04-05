# Tensor Rank of chi_P: Results

**Script:** `tensor_rank_chi_P.py`
**Session:** 17
**See also:** `tensor_rank_robust_results.md`, `k_party_nof_results.md`

## What Was Tested
True tensor rank (CP decomposition) of the prime indicator function over 3-way partition for N=6,9,12. Computed mode-unfolding ranks, tensor rank over R (gradient descent/Adam), tensor rank over F_2 (greedy), slice rank, and comparison with random {0,1} tensors.

## Key Findings
- Mode unfolding ranks: all full (= 2^{N/3}), confirming multiparty results
- Tensor rank over R: matches or exceeds random tensor rank for same shape and density
- For N=6 (4x4x4): tensor rank ~ 4; for N=9 (8x8x8): rank ~ 8; for N=12 (16x16x16): rank ~ 16
- Tensor rank = 2^{N/3} = full, matching random tensors
- Slice rank also matches random baseline
- No tensor structure exploitable for efficient computation

## Verdict
**CLOSED**
**Failure Mode:** Information loss (tensor rank is full/maximal; matches random tensors of same density)

## One-Line Summary
Tensor rank of chi_P is full (= 2^{N/3}); matches random tensors, no exploitable 3-party structure.
