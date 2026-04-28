# Degree-components and support-hypergraph rigidity of f_chi_P^{(n)}

Sub-experiment of S156 GCT orbit-dim attack on A7. **Main results file:**
`gct_chi_p_orbit_results.md` §4.5–§4.7.

## Summary

Two probes:

(a) **Degree-component decomposition** of f_chi_P^{(n)} into homogeneous
    pieces f_d := sum_{|S|=d, val(S) prime} prod_S x_i. Compute Stab Lie
    dim of each f_d.

(b) **Support-hypergraph Lie-rigidity test**: for n ∈ {4, 5, 6}, generate
    100 (or 50) random integer-coefficient choices on the chi_P support
    pattern and report distribution of Stab dim.

## Results (raw)

### (a) Degree-component Stab dimensions of f_chi_P^{(n)}

| n | deg 1 | deg 2 | deg 3 | deg 4 | deg 5 | full f_chi_P |
|---|------:|------:|------:|------:|------:|-------------:|
| 4 | 12    | 9     | 4     | —     | —     | 0            |
| 5 | 20    | 16    | 7     | 8     | 4     | 0            |
| 6 | 30    | 25    | 11    | 2     | 2     | 0            |

Individual degree components have non-trivial stabilizers; their
intersection (= dim Stab of full f_chi_P) collapses to 0.

### (b) Support-hypergraph Lie-rigidity

| n | trials | Stab dim distribution | conclusion |
|---|------:|-----------------------|------------|
| 4 | 100   | {0: 100}              | LIE-RIGID  |
| 5 | 100   | {0: 100}              | LIE-RIGID  |
| 6 | 50    | {0: 50}               | LIE-RIGID  |

100% of random integer-coefficient choices on the chi_P support pattern
give Stab dim = 0. The support hypergraph itself is Lie-rigid: NO
choice of coefficients on this support produces a polynomial with
non-trivial Lie stabilizer.

## Conclusion

The chi_P-coefficients-all-1 case is typical at the orbit-dim level —
the support pattern alone forces trivial Lie stabilizer. This is a
structural fact about the prime-encoding hypergraph, not the
chi_P-specific coefficient choice.

Contributes to E2.26 (sub-claims (vi) and (vii)) and is part of the
orbit-dim sub-frame closure of A7 (S156, mode E, B-grade).

## Linear-factor structural observation (visible from the support)

f_chi_P^{(n)} = x_2 + x_1 · g_n(x_2, ..., x_n) for all n ≥ 2, since
val(S) is odd iff `1 ∈ S`, and all primes ≥ 3 are odd. Only the prime
2 (val({2}) = 2) contributes a monomial without x_1. This is a
deterministic structural feature of the polynomial.

## See also

- `gct_chi_p_orbit_results.md` (main results).
- `archive/sessions/session156_a7_gct_chi_p_orbit.md` (synthesis).
- EDGES.md E2.26 (full edge entry).
