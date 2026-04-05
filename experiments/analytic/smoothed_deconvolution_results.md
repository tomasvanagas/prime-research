# Smoothed Explicit Formula + Deconvolution: Results

**Script:** smoothed_deconvolution.py

## What Was Tested
Gaussian-smoothed explicit formula (converges exponentially: zero at gamma contributes exp(-gamma^2*sigma^2/2)) followed by deconvolution to recover the step function pi(x). With sigma=0.1, only ~5 zeros needed for the smoothed version.

## Key Findings
- Smoothed psi_sigma(x) with sigma=0.1: zero at gamma=14.13 contributes ~0.37, gamma=101 contributes ~10^{-22}
- Effectively only 3-5 zeros needed for smoothed version -- huge reduction
- BUT: deconvolution to recover the step function amplifies noise exponentially
- Gaussian smoothing in frequency domain: F[G_sigma](k) = exp(-k^2*sigma^2/2)
- Deconvolution divides by this: amplification factor exp(k^2*sigma^2/2) for frequency k
- The step function pi(x) has Fourier content at ALL frequencies (it's a staircase)
- High-frequency content (which encodes individual prime positions) is amplified exponentially
- The deconvolution undoes the smoothing exactly, recovering the need for ALL zeros
- Staircase fitting: fitting step locations to smoothed data requires knowing approximate prime positions -- circular
- Fundamental: smoothing and deconvolution are inverse operations; no free lunch

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence -- deconvolution exactly undoes the smoothing; net effect is the original explicit formula)

## One-Line Summary
Smoothed explicit formula + deconvolution: smoothing reduces zeros but deconvolution restores the full O(sqrt(x)) requirement (no free lunch).
