# Optimized Linear Combinations of li(x^{1/k}): Results

**Script:** optimized_li_combination.py

## What Was Tested
Whether optimized coefficients c_1,...,c_K in pi(x) = round(sum c_k * li(x^{1/k})) can improve on the standard Mobius coefficients mu(k)/k, potentially achieving exactness for a wider range. Tests whether the barrier is in the basis functions or in the coefficients.

## Key Findings
- Standard R(x) uses c_k = mu(k)/k; least-squares optimized coefficients differ significantly
- Optimized coefficients reduce RMS error ~30% over standard R(x) for x in [100, 10000]
- BUT: optimized coefficients are specific to the training range and do not generalize
- For larger x, the optimal coefficients converge back toward mu(k)/k -- R(x) is already optimal
- The barrier is NOT in the coefficients but in the **basis functions**: li(x^{1/k}) are all smooth
- No linear combination of smooth functions can reproduce the O(sqrt(x)) oscillatory correction
- The information-theoretic gap: K smooth basis functions carry K*O(polylog) bits; the correction needs O(sqrt(x)) bits
- Confirms that the fundamental limit is in the BASIS (smooth functions), not the COEFFICIENTS

## Verdict
**CLOSED**
**Failure Mode:** I (Information Loss -- smooth basis functions li(x^{1/k}) cannot encode oscillatory correction regardless of coefficients)

## One-Line Summary
Optimized li(x^{1/k}) coefficients: ~30% improvement over R(x) on training set, but barrier is in smooth basis functions, not coefficients.
