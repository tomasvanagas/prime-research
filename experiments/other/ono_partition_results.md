# Ono-Craig-van Ittersum Partition-Based Prime Detection: Results

**Script:** ono_partition.py

## What Was Tested
The Ono-Craig-van Ittersum (PNAS 2024, arXiv:2405.06451) partition characterization: n >= 2 is prime iff (n^2 - n + 1)*sigma_1(n) = sigma_3(n). Tested whether this criterion provides a computational speedup for finding or counting primes.

## Key Findings
- The criterion is mathematically correct: verified for all n up to tested range
- For prime p: sigma_1(p) = p+1, sigma_3(p) = p^3+1, and (p^2-p+1)(p+1) = p^3+1 identically
- Computing sigma_k(n) requires factoring n or enumerating divisors: O(sqrt(n)) per number
- Using this as a primality test costs O(sqrt(n)) per candidate, worse than Miller-Rabin (O(log^2 n))
- For counting primes (pi(x)): must test each n individually, total cost O(x*sqrt(x)) -- much worse than Meissel-Lehmer
- Follow-up papers (MacMahonesque, Botkin) provide approximations but not exact formulas in polylog

## Verdict
**CLOSED** -- Failure Mode: C (Circularity)

## One-Line Summary
Ono's partition criterion is elegant but computationally slower than Miller-Rabin; no speedup for pi(x) or p(n).
