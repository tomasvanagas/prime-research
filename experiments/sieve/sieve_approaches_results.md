# Sieve Approaches: Results

**Date:** 2026-04-04 (Session context)
**Script:** sieve_approaches.py

## What Was Tested
Three sieve-based implementations for computing the nth prime: (1) segmented Sieve of Eratosthenes, (2) Sieve of Atkin, (3) combined tight-upper-bound + optimized odd-only sieve. All validated against known prime values and benchmarked.

## Key Findings
- All three pass validation against known primes (p(1)=2 through p(10^6)=15485863)
- Benchmarks at n=10^6:
  - Segmented Sieve: 1.00s
  - Sieve of Atkin: 3.13s (slowest in pure Python)
  - Combined (bound + sieve): 0.35s (fastest)
- Upper bound tightness: Dusart bound overhead decreases from ~14% at n=10 to <1% at n=10^6
- All are O(N) complexity (sieve the range up to p(n)), unsuitable for large n
- These are reference/baseline implementations, not improvements over O(x^{2/3}) algorithms

## Verdict
**CLOSED**
**Failure Mode:** N/A -- baseline implementations, not novel approaches; O(N) complexity cannot compete with O(x^{2/3}) Meissel-Lehmer

## One-Line Summary
Three sieve implementations benchmarked; combined approach fastest (0.35s at n=10^6); all O(N) -- baseline only.
