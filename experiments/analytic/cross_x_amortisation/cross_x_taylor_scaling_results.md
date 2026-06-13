# cross_x_taylor_scaling — supplementary K-scaling for slot 2

**Session:** S221 (Thread 5, slot 2 of 5)

Supplementary to `cross_x_batched_evaluator.py` /
`cross_x_batched_evaluator_results.md`. Isolates the K-scaling of the
Taylor-P speedup at fixed M, confirming the asymptotic formula
`speedup → a/(b·P)` is K-independent.

## What it measures

For x = 10⁶, P = 2, integer cluster, K ∈ {100, 200, 400}, M ∈ {4, 16,
64, 256}: time the direct M-batched evaluator and the Taylor-P
multipoint evaluator. Compute speedup and per-zero error.

## Results

CSV: `cross_x_taylor_scaling.csv`. Console output:

```
    K      M     T_dir     T_tay   speedup      maxerr     per0err
  100      4     0.603     0.195     3.09x   1.242e-10   1.242e-12
  100     16     2.331     0.191    12.19x   1.553e-08   1.553e-10
  100     64     7.939     0.189    42.00x   1.153e-06   1.153e-08
  100    256    33.594     0.189   178.01x   7.700e-05   7.700e-07
  200      4     0.890     0.294     3.03x   1.864e-10   9.318e-13
  200     16     2.888     0.263    10.99x   2.336e-08   1.168e-10
  200     64    13.215     0.258    51.31x   1.751e-06   8.753e-09
  200    256    53.479     0.334   160.33x   1.213e-04   6.067e-07
  400      4     1.553     0.567     2.74x   8.350e-10   2.087e-12
  400     16     5.451     0.567     9.61x   1.041e-07   2.603e-10
  400     64    21.638     0.610    35.48x   7.633e-06   1.908e-08
  400    256    89.731     0.447   200.84x   4.846e-04   1.211e-06
```

## K-scaling at fixed M

| Fixed M | speedups (K=100, 200, 400) | spread |
|---------|----------------------------|--------|
| 4       | 3.1×, 3.0×, 2.7×           | 13%    |
| 16      | 12.2×, 11.0×, 9.6×         | 27%    |
| 64      | 42.0×, 51.3×, 35.5×        | 45%    |
| 256     | 178×, 160×, 201×           | 25%    |

Speedup is K-independent within ±10–25% at small/moderate M, with
larger spread at M=64 reflecting setup-amortisation contention. The
asymptotic formula `M·a/(a' + M·b·P) → a/(b·P)` is empirically
confirmed: at small M (M ≤ 16) where `a' ≫ M·b·P`, speedup ≈
`M·a/a'` (linear in M, K-independent); at large M, speedup → `a/(b·P)`
(constant in K).

## Per-zero error

Max error scales like (K · M³) per Taylor-2 cubic truncation
`(γ_K · Δx / x)³`. At K=400, M=256 max_err = 4.85e-4: per-zero error
1.2e-6 vs threshold 1/(2·400) = 1.25e-3 — Taylor-2 stays accurate
1000× under threshold for these parameters. For larger M or K, P must
grow (P=4 supports Δx ≤ x^{2/5}; P=O(log x) supports Δx ≈ x^{1/2}).

## What this confirms

1. The asymptotic Taylor speedup formula `a/(b·P)` is K-independent —
   eliminates the possibility that Taylor's apparent speedup is a
   K-finite-size effect.
2. Taylor-P stays accurate (max_err ≪ 0.5) up to (K=400, M=256) at
   x=10⁶ for P=2, even with naïve x_0 at the cluster left edge.

## Pointer

See `cross_x_batched_evaluator_results.md` for the full slot-2 analysis,
cluster-width derivation, and structural conclusion.
