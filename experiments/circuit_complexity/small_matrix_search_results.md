# Small Matrix Search for det(M) = pi(x): Results

**Script:** `small_matrix_search.py`
**Session:** 14
**See also:** `dc_extended_results.md`

## What Was Tested
For small x, searched for integer matrices M of size k x k (k << x) whose determinant equals pi(x). Tried random search, structured search (companion/triangular matrices), gradient-free optimization, and constructive number-theoretic approaches.

## Key Findings
- Trivial 1x1: M = [[pi(x)]] always works but is not systematic
- 2x2 with integer entries: many solutions exist for each target value pi(x) (factorization of ad-bc = target)
- 3x3+: solutions exist but no systematic construction that scales with x
- Constructive approaches using number-theoretic data fail to generalize beyond specific small cases
- The ad-hoc solutions do not form a pattern that extends to arbitrary x

## Verdict
**CLOSED**
**Failure Mode:** Circularity (constructing the matrix requires knowing pi(x); no systematic scalable construction found)

## One-Line Summary
Ad-hoc small integer matrices with det = pi(x) exist but no systematic construction scales; the GapL question remains open.
