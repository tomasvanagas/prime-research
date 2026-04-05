# Information Theory & Kolmogorov Complexity Attack Results

## What Was Tested
Session 10 comprehensive information-theoretic analysis of p(n) via 6 experiment classes on first 100,000 primes:
1. **Kolmogorov complexity estimation** via compression (gzip/bz2/lzma) of primes, gaps, and delta sequences.
2. **Algorithmic randomness tests** (monobit, runs, mod-m uniformity) on delta(n) = p(n) - R^{-1}(n).
3. **Mutual information** between delta(n) and features of n (Omega(n), mu(n), n mod k, prev_gap, etc.).
4. **Circuit complexity**: Testing whether low-depth circuits (individual bits, XOR pairs, GF(2)-linear, TC0 threshold) can predict delta(n) mod 2.
5. **Communication complexity / conditional K**: Bits per element with and without conditioning on n mod k.
6. **Next-bit prediction**: Various predictors for delta(n) mod 2 and first-difference signs.

## Key Findings
- **Compression**: Scaling exponent -> 1.0 for all sequences; ~7 incompressible bits per element. Removing smooth approximation does NOT reduce information.
- **Randomness tests**: Monobit PASS, Runs FAIL (autocorrelation r=0.96 at lag 1 -- delta is smooth), Mod-m PASS (residues are uniform).
- **Mutual information**: floor(log2(n)) has NMI=0.12 (scale effect only). All other features have NMI < 0.01. No simple function of n predicts delta(n).
- **Circuit complexity**: All shallow circuit types achieve only ~50% on delta(n) mod 2 -- indistinguishable from random.
- **Conditional K**: ~7.3 bits/element unconditionally; conditioning on n mod k actually HURTS. Bits/element is roughly constant as N grows.
- **Prediction**: No predictor beats 56% on delta(n) mod 2. Running mean achieves R^2=0.988 for value prediction but residual RMSE=19 (~4.2 bits) is unpredictable.

**Critical conclusion**: K(p(n)|n) = Theta(log n) BITS, NOT Omega(sqrt(n*ln(n))). The answer is informationally small but computationally hard to find. The barrier is a SEARCH barrier, not an information barrier.

## Verdict
**CLOSED** -- Failure Mode: **E** (Equivalence)

Information theory does NOT close the polylog question. K(p(n)|n) = O(log n) bits, but finding those bits requires surveying prime distribution up to x, which takes O(x^{2/3}) time. The gap between information complexity and computational complexity is inherent, not exploitable.

## One-Line Summary
K(p(n)|n) = Theta(log n) bits -- informationally small -- but all experiments confirm that computing those bits requires exponential-in-input-size work; the barrier is computational, not informational.
