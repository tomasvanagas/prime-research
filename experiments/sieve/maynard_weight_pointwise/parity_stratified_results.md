# parity_stratified.py — results pointer

Stratify pointwise AUC by `n mod 2` to isolate genuine sieve content
from trivial parity detection. With H = {0, 2, 6} all coordinates
share `n`'s parity, so even `n` ⇒ no primes in window.

**Headline result** (selberg_gpy F, N = 10^5, window = 20000):

```
 θ=0.10:  AUC(all n) = 0.876,  AUC(odd n only) = 0.657
 θ=0.15:  AUC(all n) = 0.884,  AUC(odd n only) = 0.679
 θ=0.20:  AUC(all n) = 0.667,  AUC(odd n only) = 0.691
 θ=0.30:  AUC(all n) = 0.432,  AUC(odd n only) = 0.661
```

The 0.876 → 0.657 drop on θ = 0.10 is exactly the parity-detection
component. Restricted to odd n, AUC stays in [0.66, 0.69] across
all θ — well below the 0.95+ required for a primality witness.

**Full analysis:** see `maynard_weight_pointwise_results.md` §
"Obstruction 1 — the weight does NOT pointwise separate primes".
Per-θ result files: `parity_t{010,015,020,030}.json`.
