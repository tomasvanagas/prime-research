# Prunescu-Shunia Arithmetic Term Implementation: Results

**Script:** prunescu_shunia_impl.py

## What Was Tested
Implementation and analysis of the Prunescu-Shunia (arXiv:2412.14594v2) arithmetic term for p(n). Implemented building blocks: GCD as arithmetic term, Robinson's binomial coefficient, and analyzed intermediate value sizes at each layer.

## Key Findings
- The formula is theoretically correct: a fixed-length expression using {+, -, *, //, ^} computes p(n) for all n
- GCD arithmetic term: gcd(a,b) involves 2^{ab(ab+a+b)} -- for a=b=3, this is 2^{135}, already enormous
- Robinson's binomial: C(a,b) involves (2^a + 1)^a -- double-exponential growth
- Full Prunescu-Shunia formula: intermediate values for n=1 have ~10^78913 digits (quadruple-exponential tower)
- The formula is cosmically impractical: cannot be evaluated even for n=1 on any real computer
- Confirms that arithmetic term existence is a theoretical result with no computational implications

## Verdict
**CLOSED** -- Failure Mode: I (Information Loss)

## One-Line Summary
The Prunescu-Shunia arithmetic term is correct but has quadruple-exponential intermediate values, making it cosmically impractical.
