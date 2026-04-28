# Nisan-Wigderson partial-derivative space dim of f_chi_P^{(n)}

Sub-experiment of S156 GCT orbit-dim attack on A7. **Main results file:**
`gct_chi_p_orbit_results.md` §4.4.

## Summary

For each n ∈ {4, 5, 6, 7} and each k ∈ {0, ..., n}, computed
dim PD_k(f_chi_P^{(n)}) := dim span{∂^k f / ∂x_T : |T| = k}, plus a
matched-support random-coefficient baseline (100 trials per cell).

## Results (raw)

```
n = 4:  PD_k(chi_P) = (1, 4, 4, 1, 0); baseline std = 0.00 across all k
n = 5:  PD_k(chi_P) = (1, 5, 10, 10, 5, 1); baseline std = 0.00
n = 6:  PD_k(chi_P) = (1, 6, 15, 16, 7, 1, 0); baseline std = 0.00
n = 7:  PD_k(chi_P) = (1, 7, ...); baseline std = 0.00
```

## Conclusion

dim PD_k of f_chi_P **exactly matches** matched-baseline mean with std = 0
across 100 random matched-support trials, at every n, every k. The dim is
fully determined by the support hypergraph, NOT by the χ_P-specific
coefficient choice.

This contributes to E2.26 (sub-claim (iii)) and is part of the orbit-dim
sub-frame closure of A7 (S156, mode E, B-grade).

## Falsification status

NW-PDS extension test does NOT produce an A-grade signal for χ_P. The
partial-derivative-space invariant is "matched-baseline-typical" and
the formula-lower-bound machinery of Nisan-Wigderson does NOT
discriminate χ_P from random-coefficient versions of the same support.

## See also

- `gct_chi_p_orbit_results.md` (main results).
- `archive/sessions/session156_a7_gct_chi_p_orbit.md` (synthesis).
- EDGES.md E2.26 (full edge entry).
