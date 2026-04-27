# op_count_scaling.py — results pointer

Power-law fit `mean_ops ∝ N^α` for the mean coprime simplex tuple
count per single-n Maynard weight evaluation, restricted to odd `n`,
across N ∈ {10^4, 10^4.5, 10^5, 10^5.5, 10^6}.

**Headline result:**

```
                        mean ops   p99   max
 θ=0.20, N=10^6           4.12     8     9    α = 0.10
 θ=0.30, N=10^6           6.89    17    20    α = 0.11
 θ=0.40, N=10^6          10.77    32    40    α = 0.12
```

α is well above 0 (polylog requires α=0). Sub-θ but positive — the
simplex constraint `d_1 d_2 d_3 ≤ R` makes mean ops grow more slowly
than the box `d_i ≤ R`, but it still scales polynomially in N.

**Full analysis:** see `maynard_weight_pointwise_results.md` §
"Obstruction 2 — divisor enumeration is sub-poly but not polylog".
Per-θ result files: `op_t{020,030,040}.json`.
