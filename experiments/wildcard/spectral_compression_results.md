# Spectral Compression: Results

**Script:** spectral_compression.py

## What Was Tested
Whether the zeta zero sum in the explicit formula is compressible: top-k zeros by contribution magnitude, random subsampling, hierarchical grouping, low-rank approximation of the "zero matrix" Z[i,j] = x_i^{rho_j}/rho_j, and scaling analysis of K(x) (minimum zeros for |error| < 0.5).

## Key Findings
- Top-k zeros: contribution magnitudes decay as ~1/gamma_j (slowly); need substantial fraction of zeros below T for accuracy
- Random subsampling with reweighting: introduces O(K^{-1/2}) variance -- need K = O(1/epsilon^2) samples for epsilon accuracy
- Hierarchical grouping: nearby zeros partially cancel, but groups still contribute O(1/log(group_size)) -- not negligible
- Low-rank approximation of zero matrix: rank grows as O(sqrt(x)) for error < 0.5 -- no compression
- K(x) scaling: empirically K(x) = Theta(x^{1/2}/ln(x)), matching the theoretical prediction from the explicit formula truncation error

## Verdict
**CLOSED**
**Failure Mode:** Information Loss (I) -- zero sum contributions decay too slowly (1/gamma_j); K(x) = Theta(x^{1/2}/ln(x)) zeros needed for exact pi(x).

## One-Line Summary
Spectral compression of zero sum: K(x) = Theta(x^{1/2}/ln(x)); contributions decay as 1/gamma; no subsampling/grouping/low-rank shortcut.
