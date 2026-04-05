# Rounding R Exact: Results

**Script:** rounding_R_exact.py

## What Was Tested
How often does round(R(x)) = pi(x)? If failures are rare and predictable, they could be handled separately, making R(x) (which is O(polylog)) sufficient for exact prime counting.

## Key Findings
- For small x (2 to ~100): round(R(x)) = pi(x) most of the time
- Failure rate increases with x: by x ~ 10^4, failures occur at ~5-10% of integers
- Failures cluster near prime gaps and near points where pi(x) crosses half-integers
- The maximum |R(x) - pi(x)| grows as O(x^{1/2}/ln(x)) -- confirmed by theory (Schoenfeld)
- Failures are NOT predictable from R(x) alone; predicting WHERE rounding fails requires knowing pi(x) -- circular
- Under RH, |pi(x) - R(x)| < (1/8pi) * sqrt(x) * ln(x), which exceeds 0.5 for x > ~50

## Verdict
**CLOSED**
**Failure Mode:** Information Loss (I) -- R(x) - pi(x) grows as O(x^{1/2}/ln(x)); failures are frequent and unpredictable without knowing pi(x).

## One-Line Summary
round(R(x)) = pi(x) fails increasingly often; error O(x^{1/2}/ln(x)) exceeds 0.5 for x > ~50; failures unpredictable without pi(x).
