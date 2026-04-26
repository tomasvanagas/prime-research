# P3 — RMT moment-matching for Δ(x) = π(x) − li(x)

**Script:** `session63fresh_rmt_moments.py`

## What was tested

Hypothesis: under GUE-statistics for ζ-zeros, Δ(x) = π(x) − li(x)
behaves as a smooth-on-small-scales process whose local moments
determine its central value. If true, knowing Δ in a window
[X − H, X + H] (excluding X) should predict Δ(X) within ±0.5,
suggesting a polylog algorithm: probe H = polylog(X) values, predict
the center.

## Method

For X ∈ {1000, 5000, 10000, 50000} and H = 100:
1. Compute Δ(x) = π(x) − li(x) for x in window via direct sieve.
2. Predict Δ(X) using estimators that *exclude* the central value:
   arithmetic mean, distance-weighted (1/|x−X|) mean, median, linreg.
3. Compare to true Δ(X).

## Result

| X | true Δ(X) | mean err | weighted err | median err | std in window |
|---|-----------|----------|--------------|------------|----------------|
| 1000 | −9.61 | +0.34 | +0.00 | +0.25 | 0.68 |
| 5000 | −15.28 | −0.59 | −0.16 | −0.34 | 1.48 |
| 10000 | −17.14 | +0.70 | +0.41 | +0.57 | 1.76 |
| 50000 | −33.55 | −0.35 | −0.48 | −0.36 | 0.89 |

**RMSE across estimators:**
- zero-baseline (Δ ≈ 0): 20.89
- arithmetic mean: 0.52
- weighted mean (1/|x−X|): **0.33** ← best
- median: 0.40
- linear regression: 0.52

The weighted mean PASSES the ±0.5 criterion. Δ(x) is locally smooth
enough that a 200-point window predicts the center within half-prime
accuracy.

## The catch — circularity

This is encouraging on its own, but the algorithm built on this
predictor is **circular**: computing Δ(x) for x in [X−H, X+H] requires
π(x) at every x in the window. Two evaluation strategies:

1. **Global sieve up to X** — costs O(X log log X), worse than
   Meissel–Lehmer (already O(X^{2/3})).
2. **Anchor + primality tests** — start from a known π(x_0) (e.g.,
   from Meissel-Lehmer, cost O(X^{2/3})), incrementally update
   through (x_0, X+H] via Miller-Rabin (cost H · log²X = polylog).

Strategy 2 ostensibly seems polylog *given* an anchor. But the anchor
costs O(X^{2/3}). To compute π(N) for a *single* value of N from
scratch, we still pay O(N^{2/3}). So this is *not* a polylog
algorithm for π(X) given only X — it amortises cost only when many
queries fall in the same neighbourhood, which is not the use case.

## Verdict

**CLOSED.** Failure mode: **Circularity (C)** — the predictor works
on local data, but obtaining the local data is the same problem we
are trying to solve. There is no anchor cheaper than O(X^{2/3})
known unconditionally.

**Useful by-product:** the experiment confirms that Δ(x) is locally
smooth in a precise sense — std(Δ) over a 200-window is consistently
< 2 even at X = 50000. This is consistent with predictions from
GUE-statistics literature and strengthens the case that the
"randomness" in π(x) lives at scales ≥ √x, not ~1.

## What would change the verdict

A method to compute Δ at ONE anchor point in polylog (e.g., conditional
on RH + LLR-style faster zero summation) would, combined with this
predictor, give polylog π(X). The predictor is unconditional;
the bottleneck is the anchor.

## One-line summary

A weighted-mean predictor recovers Δ(X) from a 200-point window
within ±0.5 (RMSE 0.33), confirming local smoothness of Δ — but
computing the window itself requires solving the original problem
at one anchor, so the proposal is circular for *single* π(X) queries.
