# H-T Transfer Attempt: Results

**Date:** 2026-04-04 (Session 11)
**Script:** ht_transfer_attempt.py

## What Was Tested
Whether the Helfgott-Thompson O(x^{3/5}) technique for computing M(x) = sum mu(n) can directly transfer to computing pi(x). Analyzed the H-T decomposition (sieve by small primes, count k-almost-primes) and why it works for M(x) but not pi(x).

## Key Findings
- H-T exploits alternating sign (-1)^{omega(n)} in mu(n) -- massive cancellation between even/odd squarefree counts
- For pi(x): only need A_1 (1-almost-primes = primes), which is irreducible -- no cancellation possible
- pi(x) = pi(z) + A_1(x,z) where computing A_1 IS the same problem as computing pi(x)
- For M(x), the alternating sum A_0 - A_1 + A_2 - A_3 + A_4 has cancellation; for pi(x) we need A_1 alone
- Meissel-Lehmer IS the analogous approach for pi(x) and it costs O(x^{2/3}), not O(x^{3/5})

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence) -- pi(x) lacks the signed cancellation that H-T exploits for M(x)

## One-Line Summary
Direct H-T transfer fails: pi(x) counts primes (positive, no sign) while H-T requires alternating-sign cancellation in mu(n).
