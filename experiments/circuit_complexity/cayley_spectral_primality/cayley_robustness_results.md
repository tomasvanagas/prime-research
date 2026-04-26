# cayley_robustness.py — Generator-set robustness check for A3 closure

This script complements `cayley_spectral_primality.py` by repeating the
spectral-primality probe under three different generator sets, to
verify that the closure mechanism (n_int_eigs ≥ 2^ω(n)) is generator-
independent.

See `cayley_spectral_primality_results.md` for the main result write-up,
falsification criterion, theorem, proof sketch, and self-grading.

## What this run produced

For 342 values of n in [11, 1500] coprime to 210 (so {2,3,5,7} are all
units mod n), three generator sets:

| Generator set | Lower bound violations (n_int_eigs ≥ 2^ω(n)) |
|---------------|----------------------------------------------|
| {2,3,5} ∪ inv | 0 / 342 |
| {2,3,7} ∪ inv | 0 / 342 |
| {2,3,5,7} ∪ inv | 0 / 342 |

Mann-Whitney AUC for primes vs prime powers (the hard ω=1 sub-case):

| Generator set | AUC |
|---------------|-----|
| {2,3,5} ∪ inv | 0.509 |
| {2,3,7} ∪ inv | 0.566 |
| {2,3,5,7} ∪ inv | 0.673 |

All AUC values are far from 1.0 (perfect discrimination); 0.509 and
0.566 are statistically indistinguishable from 0.5 (no signal); 0.673
is weak signal but on a small (9-prime-power) sample. **No generator
set tested produces a clean primality discrimination.**

## Conclusion

The closure of A3 is robust to generator-set choice. The 2^ω(n) lower
bound and the cyclic-(Z/p^kZ)* obstruction apply uniformly. Adding
more generators does not help because the obstruction is at the
group-structure (cyclicity) level, not the sampling-resolution level.

## Files

- `cayley_robustness.py` — this script.
- `robustness_features.json` — the full data.
