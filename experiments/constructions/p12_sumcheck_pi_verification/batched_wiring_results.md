# batched_wiring.py — results (S503)

## What this is

The S501 wall-clock profile of the compressed π(x) chain (`run_chain`,
`--delegate`) found **two** dominant kernels. S502 batched the first (the trace
zero-test). The **second** is the **wiring delegation**: 30% of `run_chain` wall
and **76% of all sum-check calls** — `K = π(√x)` *independent* `inner_verify_div`
chains (`lucy_dp_delegated_wiring`), each on a **tiny** `2^{l_ℓ}` state cube
(`l_ℓ = ⌈log₂ p_ℓ⌉`), each paying full Python+numpy per-call dispatch on a ~30-element
array. That is **op-count-bound** — the exact regime (S500) where the fast-Mersenne
mulmod is a **net loss** on the chain.

`batched_wiring.py` is the **defer-and-batch** fix: collect the `K` wiring
obligations `(p_ℓ, r_v_ℓ, r_u_ℓ, claim_ℓ, accept_rem_ℓ)` and discharge **all of
them in ONE batched inner GKR chain**.

## Construction — data-parallel GKR with a *non-uniform, per-layer* wiring

Each `inner_verify_div` is a GKR chain: sum-check #0 (degree-2, over the `2^l`
state cube), then `n` backward-sweep layers (each a `2l`-var structured
sum-check, S496), then a base check. Batching `K` of them is **not** standard
data-parallel GKR, because the per-layer wiring is **heterogeneous**: the
transition `s' = 2s + a − b·p_ℓ` and the junk mask `[σ < p_ℓ]` depend on `p_ℓ`,
and the bit-weights `wv_ℓ[j], wu_ℓ[j]` depend on each layer's own query point
`(r_v_ℓ, r_u_ℓ)`.

The resolution (worked out and verified here):

- Add `Lk = ⌈log₂ K⌉` **layer-index variables**; MSB-zero-pad every layer's
  `2^{l_ℓ}` cube to a common `2^{l_max}` (`l_max = ⌈log₂ √x⌉`; valid states stay
  `< p_ℓ`, junk stays 0). The `K` chains share **one** `(r_L, ρ)` trajectory; a
  random `γ` combines the `K` claims.
- Keep the layer index a **genuine sum-check variable** throughout. Then the
  round-poly final factorizes as a product of stacked-table MLEs, and the
  **non-uniform** wiring is recomputed by the **verifier** in `O(K·l)` via
  `Σ_ℓ χ_ℓ(r_L) · ⟨per-layer wiring MLE⟩` (`lt_const_eval`, `r_tilde_eval` per
  `p_ℓ`, `wv_ℓ`, `wu_ℓ`) — the trace batch's `av` trick generalized to a
  point-dependent operator.
- **No dense `2^{2l}` table is ever built.** Each backward layer rides the layer
  axis *inside* the S496 **structured phase split** (phase 1 binds `σ'`, phase 2
  binds `σ`), so every prover table is `K·2^{l_max}` **wide**, not
  `2^{Lk+2l_max}` dense. The carried claim is a single scalar
  `Σ_ℓ χ_ℓ(r_L) v_j^ℓ~(ρ)`; the per-layer wiring factors of one step are
  **absorbed into the next step's verifier-known weight** (the standard GKR
  eq-tracking, here `eq(r_L1, r_L2)` between the two phase bindings).

`verify_wiring_batched` = sum-check #0 (batched, deg 2) + `n` batched backward
steps (`batch_backward_step`, two deg-4 phase sum-checks each) + a base check
`(Σ_{ℓ<K} χ_ℓ(r_L))·Π_i(1−ρ_i)` (`v_0^ℓ = e_0` for every real layer).

This **exactly** matches the two wirings `run_chain --delegate` emits: the small
division `inner_verify_div(accept_rem=None)` (`compressed_layer.small_reduce`,
line 657) and the large affine image `inner_verify_div(accept_rem=p−1)` (the core
of `verify_waff_value`, line 273) — both at the uniform depth `n = nb`.
`synth_obligations(..., 'div' | 'affine' | 'mixed')` reproduces both.

## Soundness / correctness — verified by AGREEMENT, not assertion

The construction's correctness is pinned by `--selftest` to the protocol it
replaces:

1. **Agreement** with the AND of `K` independent `inner_verify_div` — over the
   demo prime `Q=2³¹−1` **and** `BIG_Q=2⁶¹−1` (object), for `div`/`affine`/`mixed`
   wirings: honest ⇒ both accept; corrupting **any single** obligation (first,
   middle, last) ⇒ batched rejects every trial.
   - `claim` (wrong wiring value, honest chain) and `lie` (self-consistent lying
     chain) are expressible against both ⇒ asserted that the **baseline also
     rejects**.
   - `vbreak` (corrupt a leaf of the prover's submitted chain, keep the honest
     claim) is a **batched-only** soundness cheat — `inner_verify_div`
     regenerates its chain from `(r_v, r_u)` and never consumes the prover's
     vectors, so it cannot see it; the batched chain must (and does) reject it.
2. **`K=1` edge** (`x=7`: only `p=2`, `Lk=0`, no layer bits).
3. **`fast`==`object` bit-identical** accept/reject (honest + every cheat) for
   `BIG_Q`.
4. The **structural win**: one batched chain has strictly fewer fmul calls and a
   strictly larger mean width than the `K`-fold baseline.

Soundness error `~ (K + n·K·deg·l_max + …)/q`, negligible at `BIG_Q` (`~2⁻⁶¹`).
Verifier work `O(n·K·l_max) = Õ(√x)` — unchanged asymptotics (the chain already
pays `O(K)` per layer); batching is for **prover efficiency / fast-path
enablement**, not verifier scaling.

## Measurement — the op-count→width conversion + the fast-path sign flip

`--bench --field big`: the wiring delegation **in isolation**, `K` independent
`inner_verify_div` (today's chain) vs ONE batched chain. `c_drop` = fmul CALL
ratio, `w_rise` = per-fmul WIDTH ratio, `speedup` = wall ratio (within a path).

```
  n    K lmax  Lk    path ||  unb_calls  bat_calls c_drop ||   unb_w   bat_w w_rise ||    unb_ms    bat_ms speedup
  8    6    4   3  object ||      15987       7811   2.0x ||     2.9    17.1   6.0x ||     132.8      76.9   1.73x
  8    6    4   3    fast ||      15987       7811   2.0x ||     2.9    17.1   6.0x ||     436.5     211.2   2.07x
 10   11    5   4  object ||      48828      12772   3.8x ||     4.9    52.2  10.6x ||     404.2     185.2   2.18x
 10   11    5   4    fast ||      48828      12772   3.8x ||     4.9    52.2  10.6x ||    1277.6     347.2   3.68x
 12   18    6   5  object ||     116436      19127   6.1x ||     7.9   166.4  21.2x ||     942.1     545.3   1.73x
 12   18    6   5    fast ||     116436      19127   6.1x ||     7.9   166.4  21.2x ||    2889.4     471.9   6.12x
 14   31    7   5  object ||     281631      25585  11.0x ||    13.5   292.6  21.7x ||    2103.8     971.1   2.17x
 14   31    7   5    fast ||     281631      25585  11.0x ||    13.5   292.6  21.7x ||    6918.6     769.3   8.99x
 16   54    8   6  object ||     657925      36368  18.1x ||    23.4   936.6  40.0x ||    6173.6    4772.7   1.29x
 16   54    8   6    fast ||     657925      36368  18.1x ||    23.4   936.6  40.0x ||   15851.3    1308.6  12.11x
```

Reading it:

- **fmul CALL count drops** `2.0×→18.1×`, growing with `K` (toward the K-fold
  target; it tracks ~`K/3` because each batched call spans `Lk+l_max` vars at
  degree 4–5 = more work *per* call). **Per-fmul WIDTH rises** `6×→40×`. Net
  total fmul-element-work rises only ~2× — the **op-count→width conversion**
  S501 named: trade modestly more arithmetic for far fewer dispatches.
- **The fast-Mersenne sign flip is reproduced** (S500's prediction, S502's trace
  precedent): **unbatched, the fast path LOSES** (n=16: `fast` 15851 ms vs
  `object` 6174 ms — 2.6× *slower*, op-count-bound); **batched, the fast path
  WINS** (n=16: `fast` 1309 ms vs `object` 4773 ms — 3.6× faster). Speedup vs the
  unbatched fast baseline reaches **12.1× at n=16** and grows with `K`.
- Communication: the batched chain's transcript is far shorter than `K`
  independent chains (n=10: 927 vs 5166 field elems; checked `< baseline` in the
  selftest).

## Integration into `run_chain` (S504) — DONE

The defer-and-batch restructure is built. Two entry points were added to this
file and wired into `compressed_layer.run_chain`:

- `chain_obligation(p, n, l_max, r_v, r_u, accept_rem, claim, q, lie)` — builds
  ONE obligation dict from a deferred chain wiring check. Unlike
  `make_obligation`, the **`claim` is the prover-supplied wiring scalar** the
  layer's outer sum-check already pinned (honest, or a forged scalar the batch
  must reject at sum-check #0); `vs` is the honest chain over the shared `l_max`
  cube, or a self-consistent lying chain when `lie=True` (mirroring
  `inner_verify_div(lie=True)` verbatim — the `small_chain` / `waff_chain` cheat).
- `verify_wiring_obligations(raw, n, l_max, rng, stats, q)` — discharges a list
  of `(p, r_v, r_u, accept_rem, claim, lie)` tuples in ONE `verify_wiring_batched`
  call. The single replacement for the `len(raw)` inline `inner_verify_div`.

`compressed_layer.py` now threads a `defer` accumulator through
`small_reduce` / `large_reduce` / `verify_affine_region` / `verify_waff_value`:
with `run_chain(batch_wiring=True)` (delegate only), each layer **appends** its
wiring obligation instead of calling `inner_verify_div` inline, and all `2K−1`
obligations (large affine `accept_rem=p−1`, small division `accept_rem=None`,
small skipped at layer K) discharge in one batched chain after the layer loop.
The `verify_waff_value` `[e≤ecut]` comparator fold (O(2^nb), p-independent) stays
per-layer, as planned. `compressed_layer.py` **selftest §19/§19b** assert the
chain's verdict is unchanged (honest accepts & matches the sieve; `delta_pi` and
the self-consistent liar rejected) over `q` & `BIG_Q`, structured & dense, alone
and composed with `batch_trace`, plus the wiring-specific liars
(`small_forge`/`small_chain`/`waff_chain`/`waff_forge`) rejected **through** the
batched discharge.

**End-to-end re-measure** (`compressed_layer.py --n 16 --bench-combined`, BIG_Q,
delegate+structured): with both big kernels batched the global fast path is
finally a **net end-to-end win** — the S502 sign flip is resolved.

  | config | wall (ms) | comm | vs baseline |
  |---|---|---|---|
  | baseline (no batch, obj) | 17075 | 87226 | 1.00× |
  | batch_trace (obj) | 13998 | 85554 | 1.22× |
  | batch_trace+wiring (obj) | 12784 | **16509** | 1.34× |
  | batch_trace (FAST) | 25452 | 85554 | 0.67× ← **S502 loss reproduced** |
  | **batch_trace+wiring (FAST)** | **8873** | **16509** | **1.92×** |

The wiring batch alone cuts comm **5.3×** (87226→16509: `K` chain transcripts
collapse to one). FAST_BIG is a **loss** with only the trace batched (tiny
per-layer wiring cubes penalise the 24-op Mersenne mulmod), and a **1.92× win**
once both kernels are widened — 8.87 s, well under the 15.5 s target.

- The leaf/`S_0` openings remain `mle_eval` stand-ins (unchanged, orthogonal).

## What would falsify the result

- The batched chain **rejecting an honest** obligation set, or **accepting** when
  any single obligation is corrupted (`claim` / `lie` / `vbreak`) — i.e. the
  batched accept/reject disagreeing with the AND of the `K` independent
  `inner_verify_div`.
- The `fast` uint64 path disagreeing with the `object` reference (any `n`, any
  cheat).
- The fmul **CALL count not dropping** with `K` or the per-fmul **WIDTH not
  rising** with `K` vs the `K`-independent baseline (the op-count→width
  conversion failing).
- The verifier work growing faster than `O(n·K·l_max)` (breaking the `Õ(√x)`
  budget).

## Reproduce

```
python3 batched_wiring.py --selftest          # agreement + cheats + fast==object
python3 batched_wiring.py --n 14 --field big  # one honest run + cheat rejections
python3 batched_wiring.py --bench --field big # the table above
```
