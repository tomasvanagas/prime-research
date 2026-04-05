# Critique: Verification Separation -- "Verification is Cheap" Claim -- Results

**Script:** critique_verification_separation.py

## What Was Tested
Quantify the claim that verifying a candidate g = p(n) is cheaper than finding p(n). Measures: |p(n) - R^{-1}(n)| (prediction error), prime gaps at p(n), cost breakdown of primality testing vs counting in the search interval, and verification cost as a function of prediction error D.

## Key Findings
- |p(n) - R^{-1}(n)| grows as ~sqrt(p(n))/log(p(n)) -- the prediction error from the smooth approximation.
- Primality testing a single candidate: O(log(x)^6) via AKS -- genuinely cheap.
- But verifying that a prime g is specifically the NTH prime requires computing pi(g) = n. This costs O(x^{2/3}) via Meissel-Lehmer -- the SAME hard problem.
- Alternative: verify by checking no primes between g and p(n-1). But this requires knowing p(n-1), which is equally hard.
- Local verification (count primes in [g-D, g]): costs O(D) if D is small, but D = O(sqrt(x)/log(x)) from the prediction error.
- The verification-prediction separation is FALSE for this problem: verification (pi(g) = n?) is as hard as prediction (find g with pi(g) = n).
- Contrast with NP problems where verification is cheap: here the "certificate" would need to encode pi(x), which has the same complexity.

## Verdict
**CLOSED**
**Failure Mode:** Circularity (C) -- verification that g = p(n) requires computing pi(g), which is the same hard problem as finding p(n).

## One-Line Summary
Verification separation is false: verifying pi(g)=n costs O(x^{2/3}), same as finding p(n); no cheap certificate for prime index.
