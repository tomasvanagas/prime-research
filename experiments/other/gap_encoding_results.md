# Prime Gap Encoding / Prediction: Results

**Script:** gap_encoding.py

## What Was Tested
Whether prime gaps g(k) = p(k+1) - p(k) can be predicted exactly, enabling p(n) = 2 + sum g(k). Tested gap autocorrelation, local history prediction, and whether g(n) = 2*round(f(n)) for a smooth function f.

## Key Findings
- Raw gap autocorrelation is approximately 0 (gaps are essentially uncorrelated)
- Local history (k previous gaps) predicts next gap with only ~20% accuracy
- Cumulative error from gap prediction grows as sqrt(n), consistent with random walk
- PNT gives f(n) ~ ln(n*ln(n))/2 on average, but rounding errors accumulate
- The fluctuations in gaps around their mean carry the same ~5 bits/prime entropy as delta(n)

## Verdict
**CLOSED** -- Failure Mode: I (Information Loss)

## One-Line Summary
Prime gaps are essentially unpredictable from local history; cumulative prediction error grows as sqrt(n).
