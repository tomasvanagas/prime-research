# Proposal 11: Compressed Sensing on delta(n) -- Results

**Script:** proposal11_delta_compressed_sensing.py

## What Was Tested
Three approaches to exploit possible structure in delta(n) = p(n) - R^{-1}(n): (1) compressed sensing -- sample delta at strategic points and fit a sparse model; (2) interpolation using piecewise-constant structure of delta between prime gap events; (3) recursive halving via binary search on "is pi(x_mid) >= n?".

## Key Findings
- Compressed sensing requires delta to be sparse in some basis. delta(n) is a sum over ALL zeta zeros -- not sparse in Fourier, wavelet, or polynomial bases.
- Interpolation fails because delta changes at every prime gap, and gaps are themselves unpredictable (encode the same zeta-zero information).
- Recursive halving (binary search on pi(x) = n) takes O(log x) steps, each requiring sign(pi(x) - n). The sign depends on the full oscillatory correction, which costs O(sqrt(x)) to determine exactly.
- Even the COMPARISON query "is pi(x) >= n?" is as hard as computing pi(x) to unit accuracy, because the smooth approximation R(x) has error >> 1.
- delta(n) has O(log n) bits but those bits are information-theoretically incompressible (GUE-random phases).

## Verdict
**CLOSED**
**Failure Mode:** Information Loss (I) -- delta(n) is not compressible in any known basis; comparison oracle is as hard as exact computation.

## One-Line Summary
Compressed sensing / interpolation / binary search on delta(n): not sparse in any basis; comparison query as hard as exact pi(x).
