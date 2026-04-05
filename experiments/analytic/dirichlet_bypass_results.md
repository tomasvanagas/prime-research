# Dirichlet Series Bypass: Results

**Script:** dirichlet_bypass.py

## What Was Tested
Whether pi(x) can be computed from evaluations of -zeta'(s)/zeta(s) via Perron's contour integral while bypassing individual zero computation. Tests shorter contours (T << x^{1/2}) combined with combinatorial corrections (hybrid approach).

## Key Findings
- Perron's formula contour at height T gives error O(x/T); for error < 1, need T ~ x^{1/2} (Lagarias-Odlyzko)
- Short contour + combinatorial correction: this is exactly the hybrid method, giving O(x^{3/5+eps}) at best
- zeta(s) can be evaluated in O(|t|^{1/3+eps}) time, so the contour integral at height T costs O(T^{4/3+eps})
- Total cost: T=x^{1/2} gives O(x^{2/3+eps}), matching the best known combinatorial methods
- No shorter contour exists that preserves exactness -- the information in zeros between T and x^{1/2} is lost
- The hybrid approach (short contour + sieve) cannot improve past O(x^{1/2+eps}) without new ideas

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence -- Perron contour integral IS the explicit formula; shorter contours lose information from unsampled zeros)

## One-Line Summary
Perron contour integral bypass fails: short contour needs T~x^{1/2} for exactness, reducing to known O(x^{2/3}) methods.
