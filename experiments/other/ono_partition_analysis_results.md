# Ono-Craig-van Ittersum Computational Analysis: Results

**Script:** ono_partition_analysis.py

## What Was Tested
Detailed computational analysis of the Ono partition characterization using the correct definitions: M_1(n) = sigma_1(n), M_2(n) computed via explicit two-part partition enumeration. Tested whether M_2(n) has any efficiently computable closed form.

## Key Findings
- M_2(n) = sum over s1 < s2 with m1*s1 + m2*s2 = n of m1*m2; this is O(n^2) to compute naively
- The criterion (n^2-3n+2)*M_1(n) - 8*M_2(n) = 0 is verified correct for primes
- Reduces algebraically to (n^2-n+1)*sigma_1(n) = sigma_3(n), which needs only divisor sums
- sigma_k(n) computation requires factoring n: O(sqrt(n)) or O(n^{1/3}) with Pollard's rho
- No way to use this for prime counting without testing each candidate individually
- The partition-theoretic reformulation adds complexity without reducing computational cost

## Verdict
**CLOSED** -- Failure Mode: C (Circularity)

## One-Line Summary
The Ono partition characterization reduces to a divisor-sum identity that requires factorization, providing no speedup.
