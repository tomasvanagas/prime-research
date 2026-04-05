# F_2 Correlation Profile of chi_P: Results

**Script:** `f2_correlation_profile.py`
**Session:** 17

## What Was Tested
Walsh-Hadamard (Fourier over F_2) spectrum of the prime indicator chi_P and correlation with low-degree F_2 polynomials for N=4..14. Tests whether chi_P shows structure exploitable by ACC^0[2] circuits.

## Key Findings
- Fourier spectrum is flat: max |f_hat(S)| ~ 1/sqrt(2^N), matching random functions
- Correlation with degree-d polynomials decays exponentially in N for any fixed d
- No "heavy" Fourier coefficient found at any degree level
- Spectral entropy is maximal (within 1% of random baseline)
- This confirms chi_P is NOT well-approximated by low-degree F_2 polynomials

## Verdict
**CLOSED**
**Failure Mode:** Information loss (chi_P has flat Fourier spectrum, no low-degree F_2 structure)

## One-Line Summary
F_2 correlation profile of chi_P is random-like: flat Fourier spectrum, exponentially decaying correlation with low-degree polynomials.
