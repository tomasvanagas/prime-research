# Wheel Formula: Results

**Date:** 2026-04-04 (Session 8)
**Script:** wheel_formula.py

## What Was Tested
Whether wheel factorization can reduce entropy per prime or enable an analytical formula for pi(x). Tested: (1) prime distribution within wheel residue classes, (2) Legendre sieve complexity, (3) entropy reduction via wheel, (4) analytical sieve continuation, (5) fundamental limitations.

## Key Findings
- Within each residue class mod 30, gap std/mean ratio ~0.6-0.7 and entropy ~4-5 bits -- same randomness as overall prime gaps
- Legendre sieve: 2^{pi(sqrt(x))} subsets needed; for x=10^100 this is 2^{10^48} -- impossible
- Wheel concentrates primes but does NOT reduce the number to check: after mod 2310, bits/prime unchanged at ~5 bits
- Selberg parity barrier: upper/lower bounds only, NOT exact pi(x)
- Information-theoretic: removing multiples of primes up to P eliminates only ~log2(ln(P)) bits per position -- negligible
- Remaining: ~0.5*log2(x) bits per prime needed regardless of wheel size

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence) -- wheel reduces constant factors but not complexity class; O(x^{2/3}) barrier is intrinsic

## One-Line Summary
Wheel factorization reduces constants but not asymptotics; entropy per prime remains ~5 bits regardless of wheel size.
