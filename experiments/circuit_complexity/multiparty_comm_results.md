# Multi-party Communication Complexity of pi(x): Results

**Script:** `multiparty_comm.py`
**Session:** 17
**See also:** `k_party_nof_results.md`, `k_party_nof_v2_results.md`

## What Was Tested
3-party NOF communication complexity of pi(x) via tensor construction and mode-unfolding rank analysis. Built 3D tensors T[a,b,c] = pi(x) for balanced 3-way split of N bits.

## Key Findings
- All mode-unfolding ranks are FULL (= dimension of that mode)
- For N=6,9,12: max rank = 2^{ceil(N/3)} exactly
- Cut rank for every partition is maximal
- The 3-party NOF complexity grows linearly as ceil(N/3)
- No evidence of tensor structure exploitable by ACC circuits

## Verdict
**CLOSED** -- subsumed by `k_party_nof_results.md`.
**Failure Mode:** Information loss (3-party tensor has full mode-unfolding ranks; no ACC structure)

## One-Line Summary
3-party NOF tensor for pi(x) has full mode-unfolding ranks; complexity grows as ceil(N/3).
