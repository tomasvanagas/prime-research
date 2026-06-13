# Delegated-Wiring Upgrade: the Õ(√x) Verifier for Exact π(x) (S491)

> **S500 note:** `_modmul_arr` now routes `BIG_Q=2⁶¹−1` array×scalar products
> through `compressed_prover_mult_trace.fmul` (uint64 Mersenne under `FAST_BIG`,
> else exact object `%q`); imports add `fmul`/`fast_big`/`BIG_Q`. The other
> arithmetic here is Python-int scalar (`wv[a]·wu[b]`, the carry-chain
> automata), exact at any modulus, so it is untouched. Default-`q` and
> object-`big` transcripts are **bit-identical** (selftest unchanged). See the
> S500 sections of `compressed_prover_mult_trace_results.md` (the win) and
> `compressed_layer_results.md` (the chain is op-count-bound, fast path slower).

## What was built

The base protocol (`lucy_dp_verification.py`) verifies exact π(x)
unconditionally, but its verifier evaluates each layer's division-wiring
MLE W̃_p(r_v, r_u) by a p-state automaton DP — O(n·p) per layer, total
Õ(x/log x) over the wheel, and the width-dichotomy measurement proved
2p+1 is TIGHT: no explicit automaton can do better. This experiment
implements the necessary escape: **delegate the automaton DP itself to
the prover** and verify it by GKR over the chain of n matrix-vector
products, exploiting the fact that the transition matrix is known to the
verifier *symbolically* — its support is the affine relation
s′ = 2s + a − b·p, whose MLE the verifier evaluates in O(log p) by a
**width-4 carry automaton** (`aff_rel_eval`), independent of p.

Agreed layer extension: M̃(σ′,σ) = LT̃(σ′)·LT̃(σ)·R̃(σ′,σ), with
R̃ = Σ_{a,b} wv(a)·wu(b)·AFF̃_{a−bp}. Each chain layer is one
2ℓ-variable sum-check (ℓ = ⌈log₂ p⌉) with the eq̃-selector inserted
(GKR form — the layer identity only extends multilinearly through
Boolean points). Verifier per Lucy layer: O(n·log p). **Total exact-π(x)
verifier: Σ_{p≤√x} O(n·log p) = Õ(√x)** — down from Õ(x/log x).
Soundness unconditional throughout.

## Measured (run_delegated_n16.log, run_wiring_bench.log)

End-to-end at x = 2¹⁶−1: π = 6542 verified exactly (54 layers, all
wiring delegated); self-consistent DP liar rejected 10/10; **false
wiring-claim liar** (lies precisely about the delegated value) rejected
10/10.

Wiring-check verifier cost per evaluation, n = 16:

| p | automaton (ms) | delegated (ms) | delegated comm (elems) |
|---|---|---|---|
| 3 | 0.058 | 1.035 | 390 |
| 13 | 0.165 | 1.775 | 780 |
| 31 | 0.333 | 2.042 | 975 |
| 61 | 0.600 | 2.571 | 1170 |
| 127 | 1.091 | 2.769 | 1365 |
| 251 | 1.970 | 3.309 | 1560 |

Automaton grows 34× over the 84× p-range (linear, as the tightness
theorem requires); delegated grows 3.2× (tracking ℓ: 2 → 8). Crossover
at p ≈ 130 in this implementation; the gap then widens linearly in p —
at wheel scale for x = 10¹² (p up to 10⁶) the per-layer comparison is
~8 s vs ~3 ms.

Honest trade-offs: communication rises (301 KB vs 28 KB at n = 16 —
the chain transcripts; still Õ(√x·polylog) total) and the prover gains
O(n·p²)-class extra work per layer (the 4^ℓ joint-cube tables; still
dominated by its sieve cost). At small x the base protocol is the
better instrument; the delegated one wins asymptotically and at large
wheel primes — which is exactly the regime exact π(x) needs.

## Debugging record (both bugs were boundary cases, now in selftest)

1. The honest chain originally allowed transitions into junk states
   [p, 2^ℓ); junk never flows back (so the chain *total* matched the
   automaton — a trap), but the per-layer MLE identity failed off the
   support. Fix: restrict the chain to valid states, matching the
   agreed LT-killed extension on the whole cube.
2. `aff_rel_eval` ignored constant bits at positions ≥ ℓ — exposed only
   at p = 2 (d = 2 = 2^ℓ); and LT̃(σ, p) with p = 2^ℓ truncated to the
   constantly-true predicate's complement. Both p = 2 boundary cases;
   the selftest now sweeps ℓ ∈ {1..4} with constants beyond the ℓ-bit
   range exhaustively.

## What would falsify

Honest delegated run rejected or accepting a wrong count (none across
selftest + n = 16); either liar class accepted above the soundness
bound (0/20); delegated wiring cost growing linearly in p (refuted:
3.2× over an 84× range).

## S495 generalization — `inner_verify_div(accept_rem=…)`

`inner_verify_W` is now a thin wrapper over `inner_verify_div`, which adds an
`accept_rem` parameter so the SAME carry chain proves not only the division
wiring `[u=⌊v/p⌋]` (`accept_rem=None`: accept summed over all remainders) but
also the **fixed-remainder** relation `[v=p·u+m]` (`accept_rem=m`: accept only
when the final remainder is `m`). The latter is exactly the affine-image
wiring `e'=p·e+(p−1)` (`m=p−1`, `v←e'`, `u←e`) that the compressed large-side
layer needs. Only the accept weighting of sum-check #0 changes (`[σ=m]`, MLE
`eq_point`, instead of `[σ<p]`, MLE `lt_const_eval`); the backward chain sweep
(whose `LT~` masks junk states in the *agreed extension*, not acceptance) is
unchanged. Selftest case 3b checks `D~(v,u)=[v=p·u+(p−1)]` against a brute
full-table MLE and against the chain accept weight `v_n[p−1]`, with truth+1 and
a lying chain rejected. This is the primitive `compressed_layer.py --delegate`
uses to make the compressed-π(x) verifier cleanly Õ(√x).

## S496 — structured chain PROVER: Õ(x^{3/2}) → Õ(x), transcript bit-identical

After S495 made the *verifier* cleanly Õ(√x), the binding cost was the chain
**prover**: the dense backward sweep (`inner_verify_div`) materialized five
`2^{2ℓ}`-size tables per chain layer (`np.repeat`/`np.tile`/`inner_r_table`) and
ran a 2ℓ-variable sum-check over them — **O(n·p²) per delegated wiring
evaluation**, so the whole chain was `Σ_{p≤√x} n·p² ≈ Õ(x^{3/2})`, exactly the
exponent the compressed DP had removed.

`chain_layer_reduce_structured` (opt-in `structured=True`, threaded through
`inner_verify_W` / `inner_verify_div` / `verify_layer_delegated` /
`run_protocol_delegated`) removes it. The summand
`F(σ′,σ)=E(σ′)P(σ′)R_j(σ′,σ)L(σ)v_j(σ)` factors — E,P depend only on σ′; L,v
only on σ; **R_j is the sparse width-4 affine operator** (4 nonzeros per column
σ via σ′=2σ+a−bp) and v_j has support ≤ p. Splitting the 2ℓ-variable sum-check
at the σ′/σ boundary runs BOTH phases over **O(p)-size** tables:
- **phase 1** (bind σ′=y): product `E(y)P(y)G(y)`, where
  `G(y)=Σ_σ R_j(y,σ)L(σ)v_j(σ)` is built in **O(p)** (each σ<p sends weight to
  ≤4 targets y=2σ+a−bp — the multilinear y-table of the folded σ-sum);
- **phase 2** (bind σ=w): product `(κ·R*(w))·L(w)·v_j(w)`, with
  `R*(w)=R_j(ρ_out,w)`, `κ=E(ρ_out)P(ρ_out)`, both from O(p) work over
  `eq_table(ρ_out)`.

**Why the transcript is bit-identical** (so `structured=True` is a drop-in and
every downstream protocol — incl. `compressed_layer --delegate` — is byte-for-
byte unchanged): (i) sum-check challenges are drawn from `rng`, *never* from the
eval values, and both the dense (one 2ℓ-round call) and the split (ℓ+ℓ rounds)
draw exactly 2ℓ challenges in the same order from the same `rng`; (ii) the round
polynomial is the unique multilinear sum-check polynomial of the *same* function
F — G merely regroups the σ-sum, which commutes with the multilinear extension
in y — so the threaded claim, the final value and `scal['V']` coincide. The
verifier checks (`eq_point`, `lt_const_eval`, `r_tilde_eval`) are untouched, and
`accept_rem` (the affine-image path) is orthogonal (only sum-check #0's weight),
so it is a drop-in there too. Cost per chain layer drops O(p²) → **O(p)**;
whole chain `Σ_{p≤√x} n·p ≈ Õ(x)`.

### Measured (`--prover-bench`, n = 16, per delegated wiring evaluation)

| p | ℓ | dense (ms) | structured (ms) | ratio | dense/p² (µs) | struct/p (µs) |
|---|---|---|---|---|---|---|
| 7 | 3 | 26.6 | 21.1 | 1.3 | 543 | 3011 |
| 31 | 5 | 46.7 | 33.4 | 1.4 | 48.5 | 1078 |
| 127 | 7 | 120.5 | 47.6 | 2.5 | 7.47 | 374 |
| 251 | 8 | 315.4 | 58.5 | 5.4 | 5.01 | 233 |
| 509 | 9 | 1010.8 | 74.0 | 13.7 | 3.90 | 145 |
| 1021 | 10 | 3640.3 | 82.9 | **43.9** | 3.49 | 81.2 |

`dense/p²` **converges to a constant (~3.5)** → dense prover is Θ(p²).
`structured` grows ~4× over the 146× p-range (overhead-bound by the O(n·ℓ)
round structure over O(p) tables — the genuine O(p) arithmetic is still small),
so `struct/p` falls monotonically; the **ratio widens 1.3× → 43.9×**, i.e. the
exponent dropped by one. End-to-end at x=2¹⁶−1, structured and dense produce the
**identical** transcript (π=6542, comm=77133 elems, both liars rejected); the
selftest asserts per-round-polynomial equality exhaustively for small p
(cases 5–6b) and equal accept/reject + comm on truth / truth+1 / lying-chain /
corrupt-DP / false-wiring across the division and affine-image paths.

This is the last piece of open item 1: with S495's Õ(√x) verifier and S496's
Õ(x) prover, the compressed delegated π(x) chain is **Õ(√x) verifier, Õ(x)
prover, unconditional** — both wirings delegated, no p-linear term either side.
`structured` is kept opt-in (default dense) so existing artifacts stay verbatim;
flipping the default in `run_protocol_delegated` is a safe one-line change
(transcript-identical) when wiring the production prover.

### What would falsify

A structured run whose per-round polynomials differ from the dense ones (would
break the bit-identical claim — selftest cases 5/6/6b would fail); structured
and dense disagreeing on accept/reject or comm for any honest/cheat input;
`dense/p²` failing to flatten (would refute the Θ(p²) baseline) or `struct`
growing like p² (would refute the exponent drop). None observed.

## S499 — field-parameterised inner wiring (arbitrary prime modulus)

The whole inner protocol now carries a `q=Q` parameter so it runs over an
arbitrary prime field, the field lift the compressed chain needs past `n≈30`
(see `compressed_layer_results.md` Step 6 and `compressed_prover_mult_trace`
S498). What changed:

- The module imports the field-parameterised `lagrange_eval / eq_table /
  eq_point / mle_eval / sumcheck` (and `_dt / _asum`) from
  `compressed_prover_mult_trace`, which default to `q = Q = 2³¹−1`; the 2ⁿ-table
  demo (`run_protocol_delegated`) and every existing caller are therefore
  **byte-for-byte unchanged**.
- `q` threads through `lt_const_eval`, `aff_rel_eval`, `r_tilde_eval`,
  `bit_weights`, `inner_chain_vectors`, `inner_r_table`, `_modmul_arr`,
  `chain_layer_reduce_structured`, `inner_verify_div`, `inner_verify_W`.
- Arrays use `_dt(q)`: `uint64` for `q ≤ 2³¹−1`, exact-int `object` arrays
  above (so `BIG_Q = 2⁶¹−1` is exact). `_modmul_arr` keeps the `np.uint64(s)`
  scalar cast on the fast path (a Python-int scalar would promote uint64 arrays
  to float64) and a plain Python-int multiply on the object path.
- The test-only `_sumcheck_rec` (round-poly instrumentation) stays at the
  default prime; the structured reduction uses the real field-parameterised
  `sumcheck` on its production (`rec=None`) path.

Verified: the full selftest passes verbatim at the default prime, and
`inner_verify_W` / `inner_verify_div` over `BIG_Q` (object dtype, dense AND
structured) match the q-parameterised automaton ground truth, accept truth, and
reject truth+1 / a lying chain (driven from `compressed_layer.py` selftest 16).

## Cross-references

`lucy_dp_verification.py` (base protocol; this experiment imports its
primitives), `../automaton_width_dichotomy/` (the 2p+1 tightness
theorem that makes delegation *necessary*),
novel/succinct_verification_of_pi.md (updated landscape: practical
verifier exponent now √x, matching the layer count — the remaining gap
to the in-principle GKR polylog verifier is the K = π(√x) sequential
layer structure itself), GKR 2008 / Thaler PAZK (the layer-reduction
form with eq̃-selector).
