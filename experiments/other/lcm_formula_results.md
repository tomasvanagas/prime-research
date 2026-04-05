# LCM / Primorial / Korselt Approach: Results

**Script:** lcm_formula.py

## What Was Tested
Whether LCM/primorial identities, Korselt's criterion, Carmichael's function, and the prime race / Chebyshev bias can provide shortcuts for computing p(n). Also tested von Mangoldt function extraction and psi(x) mod m computation.

## Key Findings
- lcm(1,...,n)/lcm(1,...,n-1) = p if n=p^k, else 1; but computing this ratio requires factoring n (circular)
- Carmichael's lambda(n): n is prime iff lambda(n) = n-1, but computing lambda requires factoring n (circular)
- Chebyshev psi(x) = ln(lcm(1,...,x)) requires factoring all n up to x or sieving
- psi(x) mod m has no known shortcut over computing psi(x) directly
- The prime race (Chebyshev bias) gives only ~52% prediction accuracy for residue classes, far from exact
- The final approach of searching near R^{-1}(n) using Miller-Rabin still requires pi(x) to determine the index

## Verdict
**CLOSED** -- Failure Mode: C (Circularity)

## One-Line Summary
LCM, primorial, and Carmichael function approaches all require factorization, which is circular.
