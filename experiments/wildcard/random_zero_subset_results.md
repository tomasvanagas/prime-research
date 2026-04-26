# Random vs First-K Zero Sampling — Results

**Session:** 23 (fresh-perspective)
**Verdict:** CLOSED (Failure mode: Equivalence — confirms low-pass is optimal)
**Hypothesis tested:** Whether random subsets of zeta zeros yield better psi(x)
recovery than the deterministic "first K zeros" choice (compressed-sensing
analog: random measurements vs Nyquist-style low-pass).

## Setup
- Pool: M=2000 zeros (gamma in [14.13, 2515.29]).
- Subset size K=50 (= 2.5% of pool).
- 50 random subsets, 200 consecutive x values near 1e6.
- Compared: first-K, random-K (50 trials), and full-pool (M=2000).

## Results
```
Method              rms        max|res|
First-K=50          97.7       120.2
Random-K=50         404.1      ~460   (mean over 50 trials)
Full pool M=2000    42.8       64.1
```
**Random-K is ~4x worse than first-K.** This is the opposite of what
naive compressed-sensing intuition would predict.

## Why first-K dominates
The contribution of zero rho_k = 1/2 + i*gamma_k to psi(x) is weighted by
`sqrt(x) / |rho_k|`. Since `|rho_k| ~ gamma_k` and gamma_k grows roughly
linearly in k (`gamma_k ~ 2*pi*k/log k`), the *amplitudes are not uniform*:
the first 50 zeros carry several times more magnitude than zeros 1950-2000.

A random subset misses most low-frequency zeros, so it estimates the
*full pool sum* with high variance. Properly rescaling the random sum
by M/K = 40 would correct the bias but blow up the variance further:
this gives an unbiased Monte-Carlo estimate of the M=2000 partial sum,
with variance proportional to (M/K) * sum over single-zero amplitudes —
strictly worse than just taking first-K.

Compressed-sensing recovery requires the signal to be **K-sparse in some
basis**. Here the relevant signal (psi(x) - x) is *dense* in the zero
basis (every zero contributes), but with magnitudes that decay smoothly.
Compressed sensing doesn't apply.

## Falsification
The hypothesis "random subsets find lucky configurations of zeros that
nearly span the contribution space" is refuted: the **best** random
subset (rms=310) is still 3x worse than first-K (rms=98). No subset
of 50 zeros can match the magnitude budget of the first 50.

## Implication for shortcut search
A future polylog algorithm cannot simply "pick a clever K=polylog subset
of zeros". Any polylog-zero subset has total amplitude budget at most
~`polylog * (1/gamma_max)` over its support, far smaller than the
required ~`sqrt(x)`. The information bound holds.

## Files
- Script: `random_zero_subset.py`
- This results file: `random_zero_subset_results.md`

## One-line takeaway
Low-frequency zeros dominate; first-K is *near-optimal* among K-subsets;
no compressed-sensing-style shortcut is available with current zero data.
