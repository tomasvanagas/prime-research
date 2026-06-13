# s425_verify_omega1_identity — results

Companion to `s425_inverse_chipfiring_probe.py`. Extends D_Ω₁
verification to N = 256 and prints the full non-zero chip distribution
of the q-reduced form, separating chips at prime powers k ≥ 2 from
chips at other vertices.

**Result table.** N ∈ {16, 32, 64, 128, 256}:

| N   | π(N) | empirical sink | sink == π(N) | off-sink total | off-sink predicted (= #{k≥2 prime powers}) | match |
|-----|-----:|---------------:|:------------:|---------------:|-------------------------------------------:|:-----:|
|  16 |    6 |              6 |   ✓          |              4 |                                          4 |   ✓   |
|  32 |   11 |             11 |   ✓          |              7 |                                          7 |   ✓   |
|  64 |   18 |             18 |   ✓          |              9 |                                          9 |   ✓   |
| 128 |   31 |             31 |   ✓          |             13 |                                         13 |   ✓   |
| 256 |   54 |             54 |   ✓          |             16 |                                         16 |   ✓   |

**Caveat (correction to initial conjecture).** The off-sink chips are
*scattered* — they do NOT all land on prime powers k ≥ 2. E.g., at
N=32 chips land at vertices {2, 3, 5, 6, 18, 20, 28} with total 7;
none of these are prime powers k ≥ 2 individually. The TOTAL count is
preserved (chip conservation: input − sink = #{higher prime powers}),
but the distribution is graph-dynamical.

**Run command.**

```
python3 s425_verify_omega1_identity.py --N 16 32 64 128 256
```

Wall time ~30 s (dominated by N=256 q-reduction).

See `s425_inverse_chipfiring_results.md` for the closure analysis
including the depth-1 / leaf collection mechanism and the decomposition
lemma `D'_Ω₁ = D'_P + D'_higher`.
