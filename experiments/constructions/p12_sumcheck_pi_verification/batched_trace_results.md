# batched_trace.py — results (S502)

**Question (from PROGRAM.md NEXT ACTION / the S501 profile).** `run_chain`'s
wall is `fmul` (56%, 2.3M calls, ~30-elem arrays) — op-count-bound on many small
per-layer cubes. The single fattest cleanly-batchable kernel is the **trace
zero-test** (`verify_constraints`): **53% of wall** in K = π(√x) *independent*
per-layer sum-checks, each over its own nb-cube witness `build_witness(x, pₗ,
nb)`, each paying full Python+numpy dispatch on a ~30-element array. Batch the K
zero-tests into ONE wide sum-check → cut the `fmul` **count** ~K-fold while
**widening** each `fmul` ~K-fold (op-count-bound → width-bound), the regime where
the S500 fast-Mersenne mulmod finally pays off.

## Construction

The K witnesses are stacked along a layer axis of size `L = 2^Lk`,
`Lk = ⌈log₂K⌉`; the joint cube is indexed `(ℓ<<nb)|e`. ONE shared `tau`
(e-cube point), ONE shared `alpha` (intra-layer constraint combiner), ONE `beta`
(inter-layer combiner). A single degree-3 sum-check of claim 0:

```
0 = Σ_{ℓ,e} BETA_EQ[ℓ,e] · F(ℓ,e),   F = Σ_c alpha^c·constraint_c(stacked W)
BETA_EQ[ℓ,e] = (beta^ℓ if ℓ<K else 0) · eq~(tau,e)        (multilinear)
```

— the SAME flat term list as the per-layer test (`build_terms`), with the
eq-selector `EQ` replaced by `BETA_EQ`. Each honest layer's inner sum is
`F_ℓ~(tau) = 0`, so the whole sum is `Σ_{ℓ<K} beta^ℓ·0 = 0`.

Two facts keep one term list serving all layers:
- **MSB-zero-padding.** The width-`Lv_max`/`Lr_max` bit tables embed every
  layer's shorter decomposition (`u_e < 2^{Lv_ℓ} ≤ 2^{Lv_max}` ⇒ high bits 0),
  so one recomposition-weight set is correct for all.
- **Shared X, prover-supplied A.** `X` is the same global x for every layer;
  the wiring `a` is a prover table. So only the table *contents* vary per layer.

**Padding rows** `ℓ ≥ K` carry `BETA_EQ = 0` (contribute 0 to the sum regardless
of their zeroed witness tables); `ONE` is all-ones over the whole cube (its MLE
must be 1 for the −X term's identity at the random point).

**Verifier** recomputes, in O(K), `BETA_EQ~(r) = [Σ_{ℓ<K} beta^ℓ χ_ℓ(r_L)] ·
eq~(tau,r_e)` and the true wiring MLE `av = (Int(r_e)+dstart) · Σ_{ℓ<K} pₗ
χ_ℓ(r_L)` (the wiring `a(ℓ,e)=(e+dstart)pₗ` factors across the disjoint e/ℓ
variable blocks, so its MLE is the product of the per-block MLEs), then checks
`final == BETA_EQ~(r)·constraint_eval(scalars, av, X, Lv_max, Lr_max, alpha)`.
The prover's `A`/`ONE`/`BETA_EQ` scalars are ignored; `av` and the literal `X`
anchor the multiply identity `U·a+R−X=0` to the true wiring, exactly as the
per-layer test does.

## Costs (re-confirmed Õ(√x), per the NEXT-ACTION caveat)

| quantity | per-layer × K (unbatched) | batched (this) |
|---|---|---|
| sum-checks | K | **1** |
| rounds | K·nb | **Lk + nb** |
| communication (field elems) | K·4·nb | **4(Lk+nb)** |
| verifier extra | — | O(K) for `BETA_EQ~`/`av` (within the chain's existing O(K) layer budget) |

At n=16: comm **1728 → 56** field elems. The O(K) verifier term is the same
O(K) the chain already spends on K layer-reductions, so the whole-chain verifier
stays **Õ(√x)**. Soundness error ≈ `(K·#constraints + nb + K + 3(Lk+nb))/q` —
the `beta` random-linear-combination over layers (deg < K) adds `(K−1)/q`, the
shared `alpha` adds `K·#constraints/q`, the sum-check adds `3(Lk+nb)/q`; all
poly(n)/q, ≈ **1.0e-15** at n=16 over BIG_Q (2⁶¹−1).

## Headline measurement (`--bench --field big`)

K per-layer sum-checks (today's chain) vs ONE batched sum-check, on the BIG_Q
object path and the fast-Mersenne uint64 path:

| n | K | nb | Lk | path | unb_calls | bat_calls | c_drop | unb_w | bat_w | w_rise | unb_ms | bat_ms | speedup |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 8 | 6 | 4 | 3 | object | 40 944 | 13 490 | 3.0× | 3.8 | 18.1 | 4.8× | 246 | 94 | 2.62× |
| 8 | 6 | 4 | 3 | fast | 40 944 | 13 490 | 3.0× | 3.8 | 18.1 | 4.8× | 207 | 342 | **0.60×** |
| 10 | 11 | 5 | 4 | object | 113 810 | 21 340 | 5.3× | 6.2 | 56.8 | 9.2× | 653 | 242 | 2.70× |
| 10 | 11 | 5 | 4 | fast | 113 810 | 21 340 | 5.3× | 6.2 | 56.8 | 9.2× | 537 | 456 | 1.18× |
| 12 | 18 | 6 | 5 | object | 261 948 | 30 966 | 8.5× | 10.5 | 186.0 | 17.7× | 1558 | 742 | 2.10× |
| 12 | 18 | 6 | 5 | fast | 261 948 | 30 966 | 8.5× | 10.5 | 186.0 | 17.7× | 1306 | 861 | 1.52× |
| 14 | 31 | 7 | 5 | object | 605 892 | 39 110 | 15.5× | 18.1 | 341.1 | 18.8× | 3144 | 1396 | 2.25× |
| 14 | 31 | 7 | 5 | fast | 605 892 | 39 110 | 15.5× | 18.1 | 341.1 | 18.8× | 2522 | 1048 | 2.41× |
| 16 | 54 | 8 | 6 | object | 1 365 696 | 51 844 | **26.3×** | 31.9 | 1169.9 | **36.7×** | 8526 | 5484 | 1.55× |
| 16 | 54 | 8 | 6 | fast | 1 365 696 | 51 844 | **26.3×** | 31.9 | 1169.9 | 36.7× | 6042 | **1707** | **3.54×** |

**Readings.**
1. **The structural target is met exactly.** `fmul` COUNT drops ~K-fold
   (n=16: 26.3×, K=54) and per-`fmul` WIDTH rises ~K-fold (36.7×). Identical
   call/width numbers on object and fast (same transcript; the field only
   changes per-multiply cost).
2. **The fast-Mersenne sign flip S501 predicted is real.** Unbatched, the fast
   path *loses* on the chain (small arrays). Batched, the width win amortises
   `_mul61`'s ~24 vectorised ops/multiply: fast `speedup` goes
   0.60× → 1.18× → 1.52× → 2.41× → **3.54×** as n grows. At n=8 it still loses
   (width 18 too small) — honest and expected.
3. **Best combination = batched + fast.** At n=16 the trace zero-test goes from
   **8526 ms** (unbatched object — today's chain, ≈ S501's 53%-of-15s kernel) to
   **1707 ms** (batched + fast) — a **5.0×** reduction in the chain's single
   biggest kernel. Batched-object alone is 1.55×; the fast path needs the
   batch's width to pay off, and the batch's width-bound regime needs the fast
   mulmod to maximise — the two are complementary, exactly as S501 argued.

## End-to-end (`compressed_layer.run_chain(..., batch_trace=True)`)

Wired into the real chain: the K per-layer `verify_constraints` calls are
replaced by ONE `verify_constraints_batched` run before the layer loop, then
skipped inside each `large_reduce` (the masked B1/B2 lookup that pins
`g1_trace` to the certified `u_e` is unchanged). The trace witnesses are
deterministic in `(x, pₗ)` and independent of the carried chain claim, so this
is a faithful detach; the chain's soundness vs `delta_pi` / self-consistent liar
lives in phase-A and the base check, not the trace test. Selftest case 18 in
`compressed_layer.py` confirms `batch_trace=True` leaves the verdict UNCHANGED
(honest accepts & matches the sieve; `delta_pi` and the self-consistent liar
still rejected) over Q and BIG_Q, automaton and delegated+structured.

Measured wall, `run_chain` delegate+structured over BIG_Q, **object** dtype:

| n | K | baseline (today) | batch_trace | speedup |
|---|---|---|---|---|
| 12 | 18 | 2.31 s | 1.72 s | 1.34× |
| 14 | 31 | 5.34 s | 3.71 s | 1.44× |
| 16 | 54 | 15.52 s | **13.60 s** | 1.14× |

And all four corner configs at n=16 (single run each):

| config | wall | claimed π | match |
|---|---|---|---|
| baseline (object, today's chain) | 15.97 s | 6542 | ✓ |
| batch_trace (object) | 13.79 s | 6542 | ✓ |
| fast only (no batch) | 52.21 s | 6542 | ✓ |
| batch_trace + fast | 22.42 s | 6542 | ✓ |

**Honest reading — the trace batch is HALF the job.** Batching the trace (one of
the chain's two big kernels) shrinks it ~5× in isolation, but the **end-to-end
win is only ~1.1–1.4×** because the *second* big kernel — the wiring delegation
(30% of wall, 76% of all sum-check calls, all tiny 2^ℓ cubes; S501) — is still K
independent per-layer sum-checks. And **globally enabling the fast-Mersenne path
is a NET LOSS** (22 s vs 16 s baseline) precisely because `_cpmt.FAST_BIG` is a
*global* flag: it speeds the now-wide batched trace but penalises every
still-tiny unbatched wiring cube (`fast only` is catastrophic, 52 s — the chain
is dominated by small cubes). So the fast path's 5× on the trace cannot be
captured end-to-end until the **wiring is batched too**. The trace batch is the
landed first half; the wiring defer-and-batch (collect all K layers'
`(r_v,r_u,claim)` `inner_verify_div` obligations, discharge in one inner
protocol) is the precondition to globally enable the fast path and capture the
combined win — the second target S501 named.

## Selftest (`--selftest`)

1. **Structural** — `Lk = ⌈log₂K⌉`, one sum-check of `Lk+nb` rounds,
   `comm = 4(Lk+nb)`, and `comm < K·4·nb` (the comm drop), n ∈ {8,12,16}.
2. **Honest + soundness** over Q and BIG_Q (object): honest accepts; corrupting
   ANY single layer (first / middle / last) with `u_consistent` (self-consistent
   wrong quotient → multiply identity fails), `u_value` (recomposition fails),
   `r_value`, or `nonbit` is rejected 5/5.
3. **Agreement** — batched accept/reject equals the AND of the K per-layer
   `verify_constraints` it replaces (honest + a corrupted-witness control).
4. **K=1 edge** (x=7: only prime 2 ≤ V=2; Lk=0, no layer bits) accepts honest,
   rejects a cheat.
5. **Fast == object** — the uint64 Mersenne path is bit-for-bit identical to the
   object reference (same accept/reject, honest + every cheat), BIG_Q.
6. **The win itself** — instrumented `fmul`: batched calls < unbatched, batched
   mean width > unbatched, and `unbatched/batched calls > 0.4·K` (~K-fold).

## What would falsify this

- The batched test rejecting an honest chain's witnesses, or accepting when ANY
  single layer's witness is corrupted (selftest 2/3).
- Batched accept/reject disagreeing with the AND of the K per-layer tests.
- The fast uint64 path disagreeing with the object reference (selftest 5).
- The `fmul` COUNT not dropping ~K-fold, or WIDTH not rising ~K-fold (selftest 6,
  bench) — which would refute the op-count→width conversion.
- The fast `speedup` column failing to cross 1.0× and grow with n — which would
  refute the predicted sign flip (the reason batching is the precondition for
  the fast-Mersenne path to be a net win on the chain).
