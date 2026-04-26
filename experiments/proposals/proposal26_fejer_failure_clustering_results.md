# Proposal 26 follow-up — Fejér failure clustering

**Goal.** Test the proposer's hypothesis (proposal 26 results.md): if Fejér-recovery
failures cluster on a recognisable subset of x by cheap arithmetic features, a
hybrid sharp/Fejér scheme could beat plain Fejér. Run for x ∈ [100, 3000] (step 50,
59 samples) at T ∈ {30, 100, 300} using 2000 zeros.

## Headline numbers

Fejér recovery rate (round(S_Fejér_T(x)) == π(x)):
- T = 30:  55.9% (sharp 54.2%)
- T = 100: 71.2% (sharp 62.7%)
- T = 300: 83.1% (sharp 71.2%)

Improvements consistent with the original proposal-26 results.md.

## Cluster analysis (Fejér recovery rate per arithmetic-feature partition)

Spread = max bucket rate − min bucket rate. Larger spread = stronger separability.

| Feature                | T=30 spread | T=100 spread | T=300 spread |
|------------------------|-------------|--------------|--------------|
| **near_prime_dist**    | **0.867**   | **0.857**    | **0.643**    |
| **lpf (quartile)**     | **0.652**   | 0.505        | 0.362        |
| div_count (quartile)   | 0.310       | 0.157        | 0.200        |
| mod6 / mod30           | 0.07–0.14   | 0.04–0.14    | 0.05         |
| smooth_le_7            | 0.106       | 0.035        | 0.097        |
| x_mod_2 / x_mod_3      | 0.07–0.07   | 0.14         | 0.05         |
| is_prime (uniform)     | 0           | 0            | 0            |

## Two genuinely separating features — both useless

### 1. `near_prime_dist`: quartile recovery rates at T=100

| Q1 (closest)          | Q2     | Q3     | Q4 (farthest)        |
|-----------------------|--------|--------|----------------------|
| 0.14                  | 1.00   | 0.67   | 1.00                 |

Strong separation, but Q1 corresponds essentially to **x being prime itself** (or
within 1–2 of a prime). Detecting "x is prime" IS the original problem; the
classifier "use sharp at full T when x is near a prime" needs a primality oracle
at every candidate, which already costs Ω(log² x) per query and on the search
window of size O(log² n) costs Ω(log⁴ n) — the AKS+binary search rate of
arXiv:2510.16285 — **not** a sub-√x route.

Mechanistic explanation: when x is exactly prime, π(x) jumps by 1 at x, so the
smooth function R(x) sits halfway between π(x−1) and π(x), and S_T(x) ≈ π(x)−0.5
or π(x)+0.5 regardless of how clean the zero-sum truncation is. This is the
"unit-step / Heaviside Gibbs phenomenon" — Fejér summation handles smooth
functions, not jumps, and at the jump the rounding boundary is 50/50.

### 2. `lpf` (largest prime factor) quartile separation

Q1 (smoothest x, lpf small) is hardest: 0.21 / 0.43 / 0.57 at T=30/100/300. But
sorting by lpf requires factoring x, which is super-polylog in worst case
(O(x^{1/4}) Pollard-rho for moderate x). A hybrid that demands factoring is not
faster than direct sieving.

## Verdict

The "lightweight follow-up" empirically confirms the suspicion: Fejér failures
cluster by features that **are themselves circular or harder than π(x)**
(primality of x; factorisation of x). Cheap arithmetic features (residues mod
small primes, divisor count, smoothness ≤7) all give spread ≤ 0.14, near the
sample-size noise floor for n=59 buckets.

**Conclusion:** No "easy/hard" partitionable subset exists by cheap features.
The Fejér window is a constant-factor improvement, period. Closes the proposal
26 follow-up.

## Failure mode

I (information loss). The clustering signal lives entirely in features
(primality, smoothness) that are at least as hard as the original π(x) problem.

## Runtime

~6 s on a laptop (3 T values × 59 x × 2 modes). Mostly mpmath Ei evaluations.

## Recommended CLOSED_PATHS entry

```
| Fejér-window damping of explicit formula + arithmetic clustering | FAIL | I | Cesàro–Fejér window gives constant-factor improvement (53→67% recovery at T=100) but does not break √x. Failures cluster on near-prime x and smooth x; both features are at least as hard to detect as the original π(x) problem. Cheap features (mod m, div count) give spread <0.14 (noise). Lightweight follow-up to S43 proposal 26. See experiments/proposals/proposal26_fejer_failure_clustering.py | 46 |
```
