# Unbalanced Communication Complexity (Comprehensive Sweep): Results

**Script:** `unbalanced_comm_complexity.py`
**Session:** 17
**See also:** `unbalanced_partition_results.md`

## What Was Tested
Comprehensive sweep of communication matrix rank for ALL partitions k=1..N-1 of N-bit inputs, where Alice holds top (N-k) bits and Bob holds bottom k bits. Computed for all feasible N.

## Key Findings
- Rank = min(rows, cols) for ALL partitions k = 1..N-1 (FULL RANK everywhere)
- No "good" partition exists where rank grows polynomially in N
- For balanced k = N/2: rank = 2^{N/2-1} + 2 (Session 17 result; slightly below full due to smooth contribution)
- For unbalanced k: rank is exactly min(2^k, 2^{N-k}) -- full rank
- The function is "full rank everywhere" -- no partition helps

## Verdict
**CLOSED**
**Failure Mode:** Information loss (communication rank is full for all partitions; no efficient 2-party protocol exists)

## One-Line Summary
Communication matrix rank is full for ALL bit partitions of pi(x); no unbalanced split reduces complexity.
