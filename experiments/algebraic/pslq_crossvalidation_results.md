# PSLQ Cross-Validation: Results

**Script:** pslq_crossvalidation.py

## What Was Tested

Whether PSLQ relations found at a single point x (e.g., x=200) for f(x) = pi(x) - R(x) hold at other x values (x=100,200,...,10000). Tests both linear and degree-2 polynomial relations across multiple basis functions including sqrt(x), x^{1/3}, log(x), and zeta-zero oscillatory terms.

## Key Findings

- Single-point PSLQ relations are SPURIOUS: with 12 basis functions and 1 point, PSLQ is guaranteed to find a relation (underdetermined system)
- Linear relations found at x=200 fail badly at x=300, 500, 1000, etc. (residuals blow up)
- Degree-2 polynomial relations from x=500 similarly fail at other points
- Multi-point relations need as many coefficients as points, so they are trivially solvable and carry no information
- Only a relation holding across ALL x with FEWER coefficients than test points would be genuine -- none found

## Verdict

**CLOSED** -- Failure Mode: Information Loss (I). PSLQ relations for f(x) are artifacts of underdetermined systems; no genuine cross-validated identity exists.

## One-Line Summary

Cross-validation confirms all PSLQ relations for f(x) = pi(x) - R(x) are spurious -- they fail at points other than where they were found.
