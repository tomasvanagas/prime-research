# Matrix Power Encoding of pi(x): Results

**Script:** `matrix_power_encoding.py`
**Session:** 14

## What Was Tested
Whether pi(x) can be expressed as a matrix power entry (A^x)_{ij} or trace Tr(A^x * B) for a fixed-size matrix A over Z. Also tested whether pi(x) mod m is a linear recurrence sequence (LRS) for small m.

## Key Findings
- Skolem-Mahler-Lech theorem: if pi(x) were LRS, then 1_prime(x) = pi(x) - pi(x-1) would also be LRS, but the prime indicator is not eventually periodic in its support -- contradiction
- pi(x) mod m is NOT a linear recurrence mod m for any m tested (m = 2..20), with orders up to 20
- Tr(A^x * B) and det(I - zA^x) extractors cannot yield pi(x) because they also produce LRS/quasi-LRS
- The cumulative sum pi(x) inherits the non-LRS nature of the prime indicator

## Verdict
**CLOSED**
**Failure Mode:** Information loss (pi(x) is provably not a linear recurrence sequence; no fixed-size matrix power encoding exists)

## One-Line Summary
pi(x) is not an LRS (proved via Skolem-Mahler-Lech + non-periodicity of primes); no fixed-size matrix power encoding.
