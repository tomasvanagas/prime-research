# session49_mps_delta_orderings — results

N=8192 (L=13 legs), eps_rel=1e-3 SVD truncation.

delta(n) = p(n) - R^(-1)(n) for n=1..8192, std=42.050.

| ordering | max_bond | sum_bond | bond profile |
|---|---|---|---|
| identity | 64 | 252 | [2, 4, 8, 16, 32, 64, 64, 32, 16, 8, 4, 2] |
| bit_reverse | 64 | 252 | [2, 4, 8, 16, 32, 64, 64, 32, 16, 8, 4, 2] |
| gray | 64 | 252 | [2, 4, 8, 16, 32, 64, 64, 32, 16, 8, 4, 2] |
| morton | 64 | 252 | [2, 4, 8, 16, 32, 64, 64, 32, 16, 8, 4, 2] |
| two_adic | 64 | 252 | [2, 4, 8, 16, 32, 64, 64, 32, 16, 8, 4, 2] |
| sort_by_Rinv | 64 | 252 | [2, 4, 8, 16, 32, 64, 64, 32, 16, 8, 4, 2] |
| random | 64 | 252 | [2, 4, 8, 16, 32, 64, 64, 32, 16, 8, 4, 2] |
| gaussian_noise | 64 | 252 | [2, 4, 8, 16, 32, 64, 64, 32, 16, 8, 4, 2] |
| smooth_cosine | 2 | 24 | [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2] |

## Interpretation

The maximum possible bond dim at the central cut of a length-2^13 TT is
min(2^6, 2^7) = 64. **Every** ordering of delta saturates this bound at
the central cut, exactly matching the gaussian-noise baseline. The
smooth-cosine control achieves bond dim 2 throughout, confirming the
analysis is sensitive to structure when present.

The bond profile is the geometric staircase 2,4,8,16,32,64,64,...
characteristic of a maximally-entangled (volume-law) sequence under
sequential SVD with relative threshold 1e-3.

VERDICT (failure mode I -- information loss): delta carries full TT
entanglement entropy under every ordering tested. **No efficient TT
representation exists** under bit-reversal, Gray, Morton, 2-adic, or
sort-by-R^{-1}. Proposal A is closed.

The result strengthens the prior empirical evidence that delta(n) is
GUE-random-like in any "ordered binary" basis. To refute this would
require a *non-binary* tensor decomposition (e.g. continuous-rank
operator-valued tensor networks), which is outside the scope of this
test.
