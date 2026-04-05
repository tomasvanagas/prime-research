# Integer Rounding Shortcut: Results

**Script:** integer_rounding_shortcut.py

## What Was Tested
Whether a correction C(x) simpler than the full zero sum can get |li(x) + C(x) - pi(x)| < 0.5 for exact rounding: (A) R(x) instead of li(x), (B) empirical fitted polynomials in 1/log(x), (C) integer snapping, (D) hybrid few-zeros + smooth + rounding.

## Key Findings
- R(x) vs li(x): R(x) is significantly more accurate (typical error ~x^{1/2}/ln(x) vs x^{1/2}) but still far from < 0.5 for large x
- Empirical corrections in 1/log(x): fit well locally but do NOT generalize; the correction is not a smooth function of x
- Integer snapping: works trivially when error < 0.5, but error grows as O(x^{1/2}/ln(x))
- Hybrid approach: using K zeros reduces error, but still need K = O(x^{1/2}/ln(x)) zeros to get error < 0.5
- The rounding approach requires the SAME number of zeros as the standard explicit formula -- no shortcut

## Verdict
**CLOSED**
**Failure Mode:** Information Loss (I) -- the correction needed for < 0.5 error encodes O(x^{1/2}/ln(x)) bits of zeta-zero information; no simpler representation exists.

## One-Line Summary
Integer rounding shortcut (R(x) + corrections): error is O(x^{1/2}/ln(x)) from smooth terms alone; still need O(x^{1/2}/ln(x)) zeros for exact rounding.
