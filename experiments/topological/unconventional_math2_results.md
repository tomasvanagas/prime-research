# Unconventional Mathematics Session 6: Results

**Script:** unconventional_math2.py

## What Was Tested
Seven radical approaches: (1) tropical/min-plus algebra for prime sieving, (2) p-adic interpolation of p(n), (3) surreal number structure analysis, (4) non-standard analysis / hyperreal telescoping, (5) matroid-theoretic exchange for prime selection, (6) sheaf-theoretic approach on Spec(Z), (7) physics-inspired oscillatory neural architecture at zeta-zero frequencies. Tested against reference primes up to 200000.

## Key Findings
- Tropical algebra: tropical sieve is equivalent to addition-based exclusion (trivially all n >= 4); no prime detection capability
- P-adic interpolation: p(n) is not p-adically continuous for any fixed p; Mahler expansion coefficients do not converge
- Surreal numbers: no computational advantage; surreal arithmetic is isomorphic to real arithmetic for standard inputs
- Hyperreal telescoping: nonstandard sum telescopes formally but standard part extraction requires the original sum
- Matroid exchange: primes form an independent set in the "divisibility matroid" but exchange axiom gives no O(polylog) access
- Sheaf theory: global sections of structure sheaf on Spec(Z) minus finite set of primes encode those primes, but computing sections is circular
- Oscillatory neural net: linear combination of sin/cos at zeta zero frequencies with ridge regression; train RMSE ~3-5, test RMSE ~8-15; generalization fails due to finite zero truncation

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence -- tropical/p-adic/surreal are reformulations; sheaf theory is circular; oscillatory net = truncated explicit formula)

## One-Line Summary
Seven radical approaches (tropical, p-adic, surreal, hyperreal, matroid, sheaf, oscillatory net) are either reformulations or circular; no path to polylog.
