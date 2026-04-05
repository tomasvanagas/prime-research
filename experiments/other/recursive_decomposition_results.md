# Recursive Decomposition of pi(x): Results

**Script:** recursive_decomposition.py

## What Was Tested
Seven approaches to beat O(x^{1/2}) for computing pi(x) or p(n): hierarchical decomposition, parallel decomposition, algebraic decomposition via number fields, analytic continuation tricks, approximate counting with certified error bounds, binary splitting on the explicit formula, and FFT acceleration of the explicit formula.

## Key Findings
- Hierarchical decomposition: pi(x) = pi(x/2) + pi(x) - pi(x/2), but the difference requires sieving the interval [x/2, x]
- Parallel decomposition: splitting [1,x] into sub-intervals allows parallelism but total work is unchanged
- Algebraic decomposition via number fields: splitting primes by Frobenius element gives partial counts but each requires O(x^{2/3}) work
- Analytic continuation: Riemann explicit formula converges conditionally; acceleration techniques don't reduce the O(sqrt(x)) zero requirement
- Approximate counting: certified error bounds are O(sqrt(x)*log(x)) even under RH, too large for exact result
- Binary splitting: applicable to the explicit formula but each term involves a zeta zero, and O(sqrt(x)) terms are needed
- FFT acceleration: can speed up the explicit formula sum by a log factor but doesn't change the O(sqrt(x)) term count

## Verdict
**CLOSED** -- Failure Mode: E (Equivalence)

## One-Line Summary
All seven decomposition/acceleration strategies leave O(x^{1/2+epsilon}) as the best analytic complexity for exact pi(x).
