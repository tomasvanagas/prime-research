# Stern-Brocot / Continued Fraction / Farey Sequence Approach: Results

**Script:** stern_brocot_cf.py

## What Was Tested
Four approaches using continued fractions and related structures: CF expansion of p(n)/n for patterns, prime constant and Copeland-Erdos constant computation, Stern-Brocot tree encoding of primes, and Farey fractions connection to RH (Franel-Landau theorem).

## Key Findings
- CF expansion of p(n)/n: first partial quotient converges to 1 (since p(n)/n ~ ln(n)), but subsequent terms show no exploitable pattern
- Prime constant C_p = sum 2^{-p_k}: computing its digits requires knowing primes (circular)
- Copeland-Erdos constant: known to be normal but computing its digits requires prime enumeration
- Stern-Brocot encoding: path length to p/1 is O(p), just a different representation of the integer
- Farey fractions: the Franel-Landau connection to RH is about the distribution of all Farey fractions, not about individual primes
- No CF-based approach avoids the fundamental information barrier

## Verdict
**CLOSED** -- Failure Mode: C (Circularity) / I (Information Loss)

## One-Line Summary
Continued fraction and Farey approaches require knowing primes to compute; no pattern in CF(p(n)/n) is exploitable.
