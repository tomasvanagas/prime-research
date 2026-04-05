# Rowland Recurrence Acceleration: Results

**Script:** rowland_acceleration.py

## What Was Tested
Five acceleration strategies for Rowland's prime-generating recurrence a(n) = a(n-1) + gcd(n, a(n-1)): (1) matrix exponentiation / linearization, (2) p-adic analysis of GCD structure, (3) orbit prediction via waiting time statistics, (4) D-finite / algebraic recurrence fast-forwarding, (5) number wall / Pade structure analysis.

## Key Findings
- Matrix exponentiation: GCD creates fundamentally nonlinear coupling; cannot be expressed as matrix multiplication
- P-adic analysis: the recurrence has no p-adic regularity; gcd(n, a(n-1)) depends on full factorization of both arguments
- Waiting time to k-th prime output is O(p_k^2) -- proven and cannot be improved
- D-finite: the recurrence is arithmetic (not algebraic); no P-recursive or D-finite framework applies
- Number wall: neither primes nor Rowland sequence satisfies any linear recurrence of reasonable order
- Rowland is O(p_n^2) per prime -- WORSE than trial division O(p_n^{3/2}), far worse than sieving O(p_n ln ln p_n)

## Verdict
**CLOSED**
**Failure Mode:** C (Circularity -- GCD introduces irreducible arithmetic complexity; the recurrence encodes the sieve in a less efficient form)

## One-Line Summary
Rowland recurrence cannot be accelerated: GCD coupling is nonlinear, waiting time is O(p^2), worse than trial division; all 5 strategies fail.
