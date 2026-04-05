# Critique: Is delta(n) Autocorrelation Exploitable for Incremental Computation? -- Results

**Script:** critique_incremental_delta.py

## What Was Tested
Test whether the high autocorrelation of delta(n) = p(n) - R^{-1}(n) at lag 1 (~0.79) enables an incremental algorithm: predict delta(n+1) from delta(n) via AR(1) or AR(5) models, and measure conditional entropy H(delta(n+1) | delta(n)).

## Key Findings
- delta(n) computed for n = 1..5000. Autocorrelation at lag 1 is indeed high (~0.79).
- AR(1) prediction: delta_hat(n+1) = a * delta(n) + b. Residual standard deviation is still O(sqrt(log n)), too large for exact prediction.
- AR(5) prediction marginally better but residuals remain O(sqrt(log n)).
- Conditional entropy H(delta(n+1) | delta(n)) ~ 0.5 * log(n) bits -- still O(log n), meaning each step requires O(log n) new bits of information.
- The autocorrelation reflects the smooth variation of R^{-1}(n), which is already captured. The RESIDUAL after removing the smooth trend is independent noise.
- High autocorrelation is misleading: it comes from the smooth part, not from exploitable structure in the random part.
- Prediction error scales as O(sqrt(log n)), which means the predicted delta is off by O(sqrt(log n)) -- too large for exact prime identification.

## Verdict
**CLOSED**
**Failure Mode:** Information Loss (I) -- autocorrelation comes from the smooth component; residual after detrending is independent noise with O(log n) bits of entropy per step.

## One-Line Summary
delta(n) autocorrelation: high r(1) is from smooth part; residual is independent noise, O(log n) bits per step, not exploitable.
