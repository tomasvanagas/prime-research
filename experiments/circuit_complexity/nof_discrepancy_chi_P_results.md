# NOF Discrepancy of chi_P: Results

**Script:** `nof_discrepancy_chi_P.py`
**Session:** 20
**See also:** `k_party_nof_results.md`

## What Was Tested
3-party NOF discrepancy of the prime indicator function (a novel computation). Built tensor T[a][b][c] = (-1)^{chi_P(x)}, computed spectral norm of unfoldings, exact discrepancy (N=6), estimated discrepancy (N=9,12), and compared with random functions.

## Key Findings
- Spectral norm of unfoldings matches random +/-1 tensors of same density
- Exact discrepancy for N=6: ~0.25 (close to random baseline of ~0.22)
- Estimated discrepancy for N=9,12: decays roughly as 2^{-N/6}, consistent with random tensors
- Correlation with low-degree F_2 polynomials: negligible (< 2^{-N/4})
- No evidence of "high discrepancy" that would enable efficient NOF protocols

## Verdict
**CLOSED**
**Failure Mode:** Information loss (discrepancy matches random tensors; no exploitable NOF structure)

## One-Line Summary
NOF discrepancy of chi_P matches random tensors; spectral norm and correlation with low-degree polynomials are negligible.
