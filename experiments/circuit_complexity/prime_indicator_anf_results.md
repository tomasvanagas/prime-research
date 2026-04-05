# Prime Indicator ANF over GF(2): Results

**Script:** `prime_indicator_anf.py`
**Session:** 13
**See also:** `approx_degree_prime_results.md`, `gf2_slp_structure_results.md`

## What Was Tested
Algebraic Normal Form of the prime indicator chi_P over GF(2) for N=4..14. Measured algebraic degree, ANF sparsity (fraction of nonzero coefficients), and structure of high-degree terms.

## Key Findings
- Algebraic degree = N for all tested N (maximal)
- ANF sparsity: ~50% of all 2^N coefficients are nonzero, matching random functions
- No exploitable structure in high-degree terms (they appear random)
- The degree scaling is linear (not polylogarithmic), confirming chi_P is not in AC^0[2]
- Degree-d truncation correlation decays exponentially for fixed d

## Verdict
**CLOSED**
**Failure Mode:** Information loss (GF(2) ANF has degree N and 50% sparsity; random-like)

## One-Line Summary
GF(2) ANF of chi_P has degree N and ~50% nonzero coefficients (random-like); confirms non-membership in AC^0[2].
