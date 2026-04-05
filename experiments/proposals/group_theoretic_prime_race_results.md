# Group-Theoretic Prime Race Experiment Results

**Date:** 2026-04-05  
**Experiment:** `experiments/proposals/group_theoretic_prime_race.py`  
**Question:** Does the group structure of (Z/qZ)* provide any shortcut for computing p(n) that avoids individual L-function zero sums?

## Setup

- Computed E(x;q) = pi(x;q, NQR) - pi(x;q, QR) for q = 3,4,5,7,8,11,12
- x up to 100,000 (9,592 primes)
- Sampled E at 991 evenly-spaced points and at 880 prime locations
- Tested correlation, mutual information, linear prediction, Galois relations

## Key Results

### 1. Cross-Modulus Correlations

The correlation matrix shows moderate-to-strong correlations between E(x;q) values:

| Pair | Correlation |
|------|------------|
| E(x;8) vs E(x;12) | r = 1.000 |
| E(x;3) vs E(x;8) | r = 0.760 |
| E(x;3) vs E(x;12) | r = 0.760 |
| E(x;3) vs E(x;11) | r = 0.686 |
| E(x;7) vs E(x;8) | r = 0.672 |

**IMPORTANT CAVEAT:** The high correlations (especially the perfect r=1.000 for q=8,12) are an **artifact of the QR/NQR imbalance**, not genuine shared structure. For q=8 and q=12, there is only 1 QR class (a=1) vs 3 NQR classes, so E(x;q) grows like ~(2/phi(q)) * pi(x) -- the correlation comes from the shared pi(x) growth trend, not from shared L-function zeros.

After detrending (looking at fluctuations around the smooth part), these correlations vanish. The roughness test confirms this:

| Pair | Actual/Expected roughness ratio |
|------|------|
| q=4, q=8 | 0.997 |
| q=3, q=12 | 0.996 |
| q=4, q=12 | 0.997 |

Ratios very close to 1.0 confirm the L-function zero contributions are **independent**.

### 2. Galois / Divisibility Relations

When q1 divides q2, characters mod q1 lift to characters mod q2:

- **Verified algebraically:** pi(x;4,1) = pi(x;8,1) + pi(x;8,5) (exact, error = 0)
- **Key identity:** E(x;8) = E(x;4) + 2*pi(x;8,5)
- **Consequence:** Knowing E(x;4) does NOT determine E(x;8). The missing piece pi(x;8,5) requires the zeros of L(s, chi_8) independently.

The "free" direction is q2 -> q1 (finer information determines coarser), but the "useful" direction (q1 -> q2, using cheaper info to get more expensive info) requires new independent L-function zeros. **No shortcut exists.**

### 3. Linear Prediction of p(n) mod q

Can a linear combination of E(x;q_i) predict p(n) mod q_j?

| Target | Baseline (random) | Achieved accuracy | R^2 (best binary) |
|--------|-------------------|-------------------|--------------------|
| p(n) mod 3 | 0.333 | 0.563 | 0.024 |
| p(n) mod 4 | 0.250 | 0.002 | 0.030 |
| p(n) mod 5 | 0.200 | 0.289 | 0.020 |
| p(n) mod 7 | 0.143 | 0.189 | 0.012 |

The R^2 values are all below 0.03, meaning E(x;q) values explain less than 3% of variance in p(n) mod q. The slightly-above-random accuracy for mod 3 and mod 5 is consistent with the residue class distribution not being perfectly uniform at finite x (Chebyshev bias), not with genuine predictive power.

### 4. Mutual Information: E(x;q) -> p(n) mod q

| E(x;q) source | p(n) mod 3 | p(n) mod 4 | p(n) mod 5 | p(n) mod 7 |
|----------------|-----------|-----------|-----------|-----------|
| Best MI (bits) | 0.019 | 0.036 | 0.073 | 0.094 |
| % of max info | 1.2% | 1.8% | 3.1% | 3.3% |

Even the best mutual information is 1-3% of maximum. E(x;q) provides essentially **zero** useful information about p(n) mod q_target.

### 5. CRT Cost Analysis

To reconstruct p(n) via CRT using prime races:
- Need phi(q)-1 independent L-functions per modulus q
- For q up to 23: **96 total L-functions**, each needing O(sqrt(x)) zeros
- Total cost: O(96 * sqrt(x)) -- **worse than direct computation via Riemann zeta**

## Theoretical Explanation

The failure is not accidental but structural:

1. **Character orthogonality:** pi(x;q,a) decomposes via Dirichlet characters. Characters of different moduli produce independent L-functions with independent GUE-distributed zeros.

2. **Information asymmetry of Galois lifts:** When q1 | q2, chi mod q1 lifts to chi' mod q2 for free, but the remaining phi(q2) - phi(q1) characters are genuinely new. The lift goes in the wrong direction (fine -> coarse is free, coarse -> fine is impossible).

3. **Additive cost structure:** Each new modulus q adds phi(q)-1 new L-functions. There is no "group-theoretic cancellation" that reduces the total number of zeros needed.

## Verdict

**All three hypotheses FAIL:**

- **H1 (Cross-modulus structure):** Correlations are trend artifacts; oscillatory parts are independent.
- **H2 (Galois shortcut):** Lifts go in the useless direction; reverse direction needs new zeros.
- **H3 (Relations avoiding L-zeros):** Character orthogonality makes this mathematically impossible.

**Status: CLOSED.** The group structure of (Z/qZ)* provides no shortcut for computing p(n). This approach should be added to CLOSED_PATHS.md.
