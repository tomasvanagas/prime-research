# Unconventional Computational Models for p(n): Results

**Script:** unconventional_models.py

## What Was Tested
Seven unconventional approaches: precomputed oracle constant (Copeland-Erdos), computable prime-encoding constants (Chaitin analog), Euler-Maclaurin with prime zeta function, Diophantine representation (DPRM / Matiyasevich), Selberg's formula deconvolution, Mertens function shortcut, and additive number theory / Goldbach inversion.

## Key Findings
- Copeland-Erdos constant: computing its digits requires enumerating primes (circular)
- Chaitin-type encoding: a real constant alpha_P encoding all primes exists but is not computable in polylog (Omega-type incomputability)
- Euler-Maclaurin + prime zeta: P(s) can be computed from zeta via Mobius inversion, but extracting pi(x) from P(s) requires the explicit formula
- DPRM / Matiyasevich: a Diophantine equation for primes exists but has exponential-size solutions; no computational speedup
- Selberg's formula: pi(x)*log(x) + sum_{p<=x} log(p)*log(x/p)*pi(x/p) = 2x*log(x) + O(x); deconvolving gives pi(x) but at O(x) cost
- Mertens function: M(x) = sum mu(n) has the same computational complexity as pi(x)
- Goldbach inversion: decomposing 2n into p+q requires knowing which summands are prime (circular)

## Verdict
**CLOSED** -- Failure Mode: C (Circularity) / E (Equivalence)

## One-Line Summary
All seven unconventional models (oracle constants, DPRM, Selberg deconvolution, Goldbach) are circular or reduce to known methods.
