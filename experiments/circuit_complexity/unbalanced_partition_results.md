# Unbalanced Communication Complexity Partitions: Results

**Script:** `unbalanced_partition.py`
**Session:** 17
**See also:** `unbalanced_comm_complexity_results.md`

## What Was Tested
Communication matrix rank for unbalanced partitions of pi(x) where Alice has top k bits and Bob has bottom (N-k) bits, for all k=1..N-1.

## Key Findings
- For balanced partition (k = N/2): rank = 2^{N/2-1} + 2 (matches Session 17)
- For all other partitions: rank = min(2^k, 2^{N-k}) (full rank)
- The balanced partition is the ONLY one with slightly sub-maximal rank (by factor ~2 from the smooth part)
- No partition achieves polynomial rank in N
- Communication complexity lower bound: Omega(min(k, N-k)) bits for any partition

## Verdict
**CLOSED** -- confirms `unbalanced_comm_complexity_results.md`.
**Failure Mode:** Information loss (all partition ranks are full or near-full; no efficient protocol)

## One-Line Summary
All partition ranks for pi(x) are full or near-full; balanced partition has rank 2^{N/2-1}+2, all others are exactly full.
