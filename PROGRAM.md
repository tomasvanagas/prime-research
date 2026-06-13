# PROGRAM.md — the live research program

This is the single state document of the continuous research loop.
Every cycle reads it first and updates it last. Keep it honest, keep
it current, keep it short — history goes to `archive/`, results go to
`experiments/` and `novel/`.

## Goal

An O(polylog) algorithm computing the n-th prime p(n) (equivalently
π(x)) exactly. Target: p(10^100) in < 1 s, 100% accurate.

**Honest position:** every known computation route is blocked — 730+
approaches catalogued in `status/CLOSED_PATHS.md`, four obstruction
theorems, the information-theoretic barrier (~√x zeta-zero phases,
GUE-random, behind the last ~50% of digits). The goal may require
mathematics that does not exist. The deliverable of each cycle is
honest progress: working algorithms/protocols, new structural facts
with measurements, precisely-closed questions with mechanisms,
adjacent-problem wins. Never a faked breakthrough.

## ACTIVE LINE — the verification program (started S491, 2026-06-13)

Computation of π(x) is blocked at √x; *verification* is not. Built and
measured (all unconditionally sound, all with cheating-prover tests):

| protocol | verifies | verifier | status |
|---|---|---|---|
| wheel sum-check | Φ(x, polylog wheel) | polylog | done |
| Lucy base | exact π(x) | Õ(x/log x) | done — π(2²⁰−1) verified |
| Lucy + delegated wiring | exact π(x) | **Õ(√x)** | done — crossover p≈130 measured |
| compressed-prover mult-trace (step 1) | batched `Σ_d eq(ρ,d)S̃(⌊X/(dp)⌋)` | polylog `O(m·Lv)` | done — trace-certified quotient, 6 cheat classes rejected |
| compressed Lucy layer (step 2) | one large-side layer over the √x state | polylog + O(nb·p) wiring | done — affine/trace dispatch, 5 cheat classes, outputs = true MLEs |
| compressed two-sided chain (step 3) | **exact π(x) over the √x state, all K layers** | Õ(√x) (waff O(nb·p) undelegated) | done — π(x)=sieve at n≤16; 9+3 cheat classes incl. self-consistent liar |
| compressed chain + delegated wiring (step 4) | same, **both wirings delegated** | **cleanly Õ(√x)**, O(nb·log p)/layer, no p-linear term | done — π(x)=sieve at n≤16; 13 cheat classes incl. 4 delegation-specific liars |
| structured chain prover (step 5) | same, **prover Õ(x)** not Õ(x^{3/2}) | unchanged Õ(√x) | done — bit-identical transcript; prover/chain-layer O(p²)→O(p), ratio→44× @p=1021 |
| structured prover IN the compressed chain (step 6) | same, on the REAL compressed √x-state chain | unchanged Õ(√x) | done — `--structured` threaded through `compressed_layer.py`; bit-identical (selftest 15); `--prover-bench` ratio→46× @p=1021, both wirings |
| field-lifted compressed chain (step 7) | same, over an **arbitrary prime** (no `u·a<q` shortcut) | unchanged Õ(√x) | done (S499) — `--field {q,big,small}`; bit-identical at default q; chain over `BIG_Q=2⁶¹−1` (object dtype) `claimed π==sieve` n≤12 all modes; chain-config alias rejected by lift, accepted by demo prime |
| fast Mersenne-61 field path (step 8) | same, `BIG_Q=2⁶¹−1` on uint64 not object | unchanged Õ(√x) | done (S500) — numpy `_mul61`/`_sum61` (`2⁶¹≡1` fold), threaded chain-wide via `fmul`; bit-identical fast-vs-object (selftests §8/§17). **Speedup is large-array-only:** certifier wins 1.7–3.8× at n=24–28; the **chain is SLOWER** (0.2–0.3×, op-count-bound on many small cubes) — see correction below |

Canonical doc: `novel/succinct_verification_of_pi.md` (protocol stack,
width-spectrum theory, census law, open program). Code:
`experiments/constructions/p12_sumcheck_pi_verification/`,
`experiments/constructions/automaton_width_dichotomy/`,
`experiments/analytic/local_pattern_census/`.

### Open items (priority order)

1. **Compressed prover** — prover Õ(x) instead of Õ(x^{3/2}); the gate
   to verifying π(10⁹⁺) against real recomputation. Design: small side
   value-addressed (division-by-constant wiring); large side
   d-addressed; the cross map u = ⌊x/(dp)⌋ (variable×variable, the
   small-side crossing when dp>√x) certified by a trace-table sum-check.
   **Step 1 DONE (S492):**
   `experiments/constructions/p12_sumcheck_pi_verification/compressed_prover_mult_trace.py`
   — the batched primitive verifying C = Σ_d eq(ρ,d)·S̃(⌊X/(d·p)⌋) with
   each u_d certified through prover trace tables (multiply+remainder
   identity U·a_d+R−X=0, bit-range recompositions, remainder bound) — a
   degree-3 constraint zero-test + a routing lookup (B1 deg-2, B2
   deg-Lv+1). Verifier polylog O(m·Lv), never divides; prover works the
   D≈√x cube. 6 cheat classes rejected; verifier flat in batch size D
   (16→256× D, ~2× verifier). Field caveat (demo q=2³¹−1 holds because
   products<q at n≤16; production lifts |F|>X or a small-field carry
   trace) and leaf-opening stand-ins documented in the results file.
   **Step 2 DONE (S493):**
   `experiments/constructions/p12_sumcheck_pi_verification/compressed_layer.py`
   — the compressed Lucy DP substrate (`small[v]`, `large[e]`, validated
   π(x)=large[0] vs sieve n≤20) plus ONE large-side layer verifier. Phase-A
   `B(e)=[e≤Bcut]` folds validity+p² gate into one threshold; phase-B splits
   `g1=g1_aff+g1_trace` and dispatches each large key by `(e+1)p≤V` (affine
   map e→p·e+(p−1), verified by a long-division `waff_eval` automaton, NO
   trace) vs `>V` (step-1 band-masked trace lookup into the small side). Two
   load-bearing facts keep sizes at √x: the verified band mask T(e) zeroes the
   affine/padding columns out of the routing, and the lookup uses the LOW
   nb=bits(V) bits of the certified u_e (so omega/S_small are 2^nb≈√x even
   though build_witness certifies u_e with Lv≈bits(x/p) bits). Line-batches
   s1@r_v and s2_aff@r_u; outputs verified == true MLEs of S_{i-1}^{small,large};
   5 cheat classes rejected; comm O(nb²); verifier polylog except waff's O(nb·p)
   (delegatable to O(nb·log p) by lucy_dp_delegated_wiring's carry chain).
   **Step 3 DONE (S494):** same file — `small_reduce` (value→value division
   wiring, gate is a BAND `[p²≤v≤V]` not a bare comparator — the planned
   "comparator suffices" holds only when V=2^nb−1/no-padding; empty-band p²>V
   guarded), `run_two_sided_layer` (width-2 claim vector, same-table small
   claims `s_B`+2 folded by line restriction), and `run_chain` — the full
   K=π(√x)-layer two-sided chain from `S_K^large(e=0)=π(x)` to S₀ (opened both
   sides). `claimed π(x)==sieve π(x)` for n≤16; prover Õ(x) (works the √x cube,
   never the x table — the open-item-1 target vs the demo's Õ(x^{3/2}));
   self-consistent-liar caught at the corruption layer (i₀∈{1,K/2,K}); 9 layer
   + 3 chain cheat classes all rejected. **Core of open item 1 is now built.**
   **Step 4 DONE (S495):** same file, `--delegate` — both `O(nb·p)` wiring
   automata removed. Small side: `w_div_eval` routed through the EXISTING
   `lucy_dp_delegated_wiring.inner_verify_W` carry chain (drop-in, as planned).
   Large side `waff`: generalized that chain to `inner_verify_div(accept_rem=m)`
   (fixed-remainder relation `[v=p·u+m]`, only sum-check #0's accept weighting
   changes); `verify_waff_value` folds the `[e≤ecut]` comparator AS-IS via an
   inner e-cube sum-check (its MLE is the O(nb) `le_eval`) and delegates only
   the pure affine-image relation `[e'=p·e+(p−1)]` (`accept_rem=p−1`). **Verifier
   wiring now O(nb·log p)/layer — `--wiring-bench` confirms the delegated cost
   tracks `l=⌈log₂p⌉` while the automaton tracks `p` (waff crossover p≈251,
   division p≈1021); whole-chain verifier is cleanly Õ(√x)** (`Σ_{p≤√x}
   nb·log p ≈ nb·θ(√x)`, vs the automaton's secret `Σ nb·p ≈ Õ(x)`). π(x)=sieve
   unchanged at n≤16; 13 cheat classes rejected incl. 4 new delegation liars
   (`waff_forge`/`small_forge` = corrupt-table-and-forge-matching-scalar, caught
   by the inner sum-check; `waff_chain`/`small_chain` = lie in the chain, caught
   by its base check). **The verifier target of open item 1 is met.** Honest
   trade-offs: at the demo's tiny x (p≤V≤255, below crossover) the delegated
   verifier is ~2.3× slower wall-clock and comm rises ~5× — the win is the
   *scaling*; AND the delegated chain's dense `2^{2l}` prover tables cost
   `O(nb·p²)`/layer → prover `Õ(x^{3/2})`. So the *new* bottleneck is the chain
   PROVER (reducible to `Õ(x)` via the sparse width-4 affine operator).
   **Step 5 DONE (S496):** `lucy_dp_delegated_wiring.py`,
   `chain_layer_reduce_structured` (opt-in `structured=True`, threaded through
   `inner_verify_W`/`inner_verify_div`/`verify_layer_delegated`/
   `run_protocol_delegated`) — the dense backward sweep's five `2^{2ℓ}`-size
   tables are gone. F(σ′,σ)=E(σ′)P(σ′)R_j(σ′,σ)L(σ)v_j(σ) factors; splitting the
   2ℓ-var sum-check at the σ′/σ boundary runs phase 1 over `E·P·G` with
   `G(y)=Σ_σ R_j(y,σ)L(σ)v_j(σ)` built in O(p) (R_j is sparse: 4 nonzeros/column)
   and phase 2 over `(κR*)·L·v` with `R*(w)=R_j(ρ_out,w)`, κ=E(ρ_out)P(ρ_out),
   both O(p). **Per chain layer O(p²)→O(p); whole chain Σ_{p≤√x} n·p ≈ Õ(x).**
   **Transcript is bit-identical** (challenges come from rng not evals; round
   polys are the same multilinear sum-check of the same F — G just regroups the
   σ-sum), so it is a verbatim drop-in (verifier checks untouched; works the
   `accept_rem` affine path too). `--prover-bench`: dense/p² flattens to ~3.5
   (Θ(p²)); structured stays overhead-bound, ratio widens 1.3×→**43.9×** at
   p=1021. Selftest asserts per-round-poly equality exhaustively (small p) +
   equal accept/reject+comm on truth/truth+1/lying-chain/corrupt-DP/false-wiring.
   **Open item 1 is now complete: verifier Õ(√x) (S495) + prover Õ(x) (S496),
   unconditional, both wirings delegated, no p-linear term either side.**
   Kept opt-in (default dense) so existing artifacts stay verbatim; flipping the
   default is a safe one-line change. Next: thread into `compressed_layer.py`.
   **Step 6 DONE (S497):** propagated the structured Õ(x) wiring prover INTO the
   real compressed-√x-state π(x) chain. `structured` threaded through
   `compressed_layer.py`'s `verify_waff_value → verify_affine_region →
   large_reduce` and `small_reduce`, reaching the S496-ready
   `inner_verify_div(..., structured=...)` for BOTH delegated wirings (small
   division `accept_rem=None`, large affine image `accept_rem=p−1`). The
   `verify_waff_value` comparator fold is already `O(2^nb)` (p-independent),
   untouched. Bit-identical transcript: selftest case 15 asserts identical
   accept/reject + comm for the two primitives, + comm for the two-sided layer
   over 11 cheat classes, + comm + claimed π(x) for the full chain (honest,
   delta_pi, self-consistent liar) at n∈{8,10,12}; default stays dense so all
   prior artifacts are verbatim. New `--prover-bench` (mirrors S496's): the
   compressed chain's wiring prover flattens `p²→p` — `dense/p²`→~3.3 (Θ(p²)),
   structured grows ×4 while dense grows ×149 over p∈[7,1021], **ratio 1.25→46×**
   for BOTH wirings. **Whole compressed delegated π(x) chain prover is now Õ(x)**
   (`Σ_{p≤√x} nb·p² ≈ Õ(x^{3/2}) → Σ nb·p ≈ Õ(x)`), verifier still Õ(√x). Open
   item 1 is complete end-to-end on the real compressed chain.
   **Step 7 DONE (S499):** field lift threaded through the WHOLE chain — `q` param
   on `compressed_layer.run_chain` + `lucy_dp_delegated_wiring` inner wiring + the
   two `lucy_dp_verification` automata, importing the field-parameterised helpers
   from `compressed_prover_mult_trace` (default `q=Q`, bit-identical; `_dt(q)`
   object arrays for `q>2³¹−1`). `--field {q,big,small}`. Chain over `BIG_Q=2⁶¹−1`
   gives `claimed π==sieve` (n≤12, all 3 modes), all cheats rejected; the
   chain-config wrap-around alias the demo prime accepts is rejected by the lift
   (selftest 16). **The demo's `u·a<q` shortcut is removed** — the protocol is now
   sound for `n≲60` over `big`. Remaining gap to a large-n run is pure speed (a
   numpy 2⁶¹−1 mulmod; see NEXT ACTION), not soundness.
2. **Layer batching: CLOSED NEGATIVE (S491-final)** — balanced product
   tree fails on fill-in growth 2^j per row; cost m·K²/2 ≈ x^{3/2}/ln²x,
   same as unbatched. Remaining form of the question: algebraic
   compression of the semigroup generated by v ↦ ⌊v/p⌋. Do not re-run
   the tree analysis.
3. **L_π = {(x,c): π(x)=c} ∈ NP?** — succinct unconditional
   certificates: no construction, no barrier. Complementary formal
   target: is π(x) #P-hard? Either answer is major.
4. **Census entropy law at k = 128** — conjecture A_k = e^{(1+o(1))π(k)}
   (enumeration-verified to k = 64, geometric convergence). Needs an
   analytic count: surjection inclusion-exclusion has catastrophic
   cancellation; DFS is enumeration-bound at A₁₂₈ ~ e²⁵. Transfer-
   matrix / generating-function attack is the well-posed route.
5. **Large-x benchmark** — gated on item 1.

## SECONDARY LINES (valid, not currently active)

- `OPEN_POSITIVE_TARGETS.md`: P7.b (batched dyadic π), P8 (sparse
  precision), P9 (quantum batched), P10 (adaptive queries) remain
  open. **P6 carries an S202-duplication warning — read it before
  spending anything there.**
- Width-spectrum follow-ons: D_k(x) census at k = 32 still converging
  (189 dense constellations unrealized at 2²⁶ — first-occurrence
  probe); sliding-window variant unmeasured.
- The guess-comparison geography (`experiments/analytic/
  guess_comparison_oracle/`) — decision-version facts, complete as is.

## Done this era (S491–S501 cycles, 2026-06-13)

- Chain wall-clock attribution + premise correction (S501):
  `experiments/constructions/p12_sumcheck_pi_verification/chain_profile.py` — a
  read-only two-pass profiler over `run_chain` (cProfile for tottime/cumtime/
  ncalls + caller breakdowns; a count-only monkeypatch for per-primitive mean
  array size). **Corrected the NEXT-ACTION premise:** `sumcheck` never touches
  `stats` and every chain `sumcheck` call is OUTSIDE any t_prover/t_verifier
  region, so the reported `t_prover+t_verifier` is only **6.4% of wall** (n=16,
  BIG_Q; 7.1% uint64) — the per-layer stats earlier cycles printed undercount
  the true wall ~15×. **Attribution (n=16, BIG_Q object, K=54):** `fmul` is
  **55.6% of wall** over **2.32M calls, mean array size ~30 elems** (max 256) —
  the op-count-bound signature S500 predicted; +27% is `sumcheck`'s own Python
  loop (cumtime 93%). Logical: `verify_trace_region`/`verify_constraints` (the
  trace zero-test) is the **fattest single kernel at 53% of wall** in just 54
  calls (deg-3, up to 133 terms over Lv+2Lr bit tables); the wiring delegation
  (`chain_layer_reduce_structured`) is **76% of all sumcheck CALLS** (1712) on
  tiny 2^l-cubes = 30% of wall. **Highest-leverage target: batch the K
  independent per-layer trace zero-tests into ONE random-linear-combination
  sum-check over a (K·2^nb) stacked cube** — cuts fmul count ~K-fold, widens
  each ~K-fold, converting op-count-bound→width-bound, the precondition that
  makes S500's fast-Mersenne path a net win. Selftest (~19s): faithfulness
  (instrumented==reference π/comm/accept over q+BIG_Q, both modes), instrument
  agreement (cProfile ncalls==size-pass counts), the 6–7% gap, fmul-dominant +
  small-array, sumcheck-cumtime-dominant. See `chain_profile_results.md`.
- Fast Mersenne-61 field path + the chain-bottleneck correction (S500): built the
  numpy `mulmod` for `BIG_Q=2⁶¹−1` (`_mul61`/`_sum61`/`_reduce61`, split-limb
  schoolbook folded by `2⁶¹≡1`, all intermediates `<2⁶⁴`) and threaded it
  chain-wide — `compressed_prover_mult_trace`'s kernels route array products
  through `fmul`; `compressed_layer`/`lucy_dp_delegated_wiring` inline products
  through `fmul`/`_modmul_arr`, with object-accumulation fallbacks for the
  `np.add.at` sums. **Bit-identical** fast-vs-object (certifier selftest §8 +
  primitives vs Python-int; whole-chain selftest §17 over 3 modes incl. cheats +
  alias). **Two honest corrections of the S499 premise:** (1) the pure-Python
  Lucy DP is **negligible** (`<0.1%` of `run_chain`: 0.5 ms @n=16 vs 7 s chain) —
  speeding it is NOT the lever; `run_chain` is ~100% field arithmetic. (2) the
  mulmod is a **large-array-only** win: `_mul61` is ~24 vectorised ops/multiply,
  so it amortises only at `D≳2000` — the certifier (one `√x`-cube) wins
  **1.7–3.8×** at n=24–28 (→10×+ at the n=32–34 target), but the **full chain is
  SLOWER** (0.2–0.3× at n≤16) because it is `π(√x)` layers of **many small**
  cubes, i.e. *operation-count-bound*. Kept opt-in (`--fast-big`); BIG_Q chain
  stays on the object reference by default (faster at reachable n). The real
  lever for a fast large-x chain run is reducing the **count** of small per-layer
  field-ops (cross-cube batching), not the per-multiply cost. (`--bench-big` on
  both files; results in the two `_results.md`.)
- Field lift threaded through the compressed chain (S499): added a prime-modulus
  parameter `q` to `compressed_layer.py`'s entire `run_chain` call graph AND
  `lucy_dp_delegated_wiring.py`'s inner wiring (`inner_verify_div`/`_W`,
  `chain_layer_reduce_structured`, the carry-chain automata) + the two
  `lucy_dp_verification` automata (`ge_const_eval`, `w_div_eval`). Both files now
  import the field-parameterised helpers from `compressed_prover_mult_trace`
  (default `q=Q=2³¹−1`, uint64 fast path; `_dt(q)`→`object` arrays for larger
  primes), so the default-q transcript is **bit-identical** (every prior selftest
  passes verbatim, incl. the structured-vs-dense bit-identity). `--field {q,big,
  small}`. **Verified (selftest 16):** the whole chain runs over `BIG_Q=2⁶¹−1`
  (object dtype) at n∈{8,10,12}, all three modes — `claimed π==sieve`, all chain
  cheats rejected — and the lift CLOSES the wrap-around aliasing hole in the
  chain's exact trace config (`build_witness(x,p,nb,dstart=1)`): at n=22 (x>SMALL_Q)
  the demo-too-small prime accepts a forged alias that `BIG_Q` rejects 5/5.
  **Honest scope:** the n>31 *soundness* failure lives in the trace certifier
  (field-parameterised, tested at chain config); a *full* n>31 `run_chain` is
  prover-bound (`compressed_lucy` is an O(x/log x) pure-Python DP + √x-size object
  cubes × π(√x) layers), not field-bound — so whole-chain correctness over the
  lift is shown at small n, the n>31 evidence at the trace certifier (the seat of
  the hole). Speed fix for a large-n object-dtype run = a numpy Mersenne mulmod
  for 2⁶¹−1 (no structural change).
- Field lift + F_{q²} refutation (S498): made `compressed_prover_mult_trace.py`
  run over an arbitrary prime modulus (`--field {q,big,small}`, uint64 fast path
  for q≤2³¹−1, exact Python-int object arrays for larger primes; backward-compatible
  defaults so `compressed_layer.py` is untouched and its selftest still passes).
  **Corrected the S497 NEXT ACTION**: the planned **extension field F_{q²} does NOT
  lift the integer range** — |F_{q²}|=q² counts ELEMENTS but the CHARACTERISTIC stays
  q, so a wrap-around alias (off by k·q) embeds as the SAME element of F_{q²} as the
  truth and is undetectable (`--refute-q2`: zero-test sum = (0+0t) for honest AND
  forged; ≠0 only over a true larger prime). The genuine hole is real and now
  demonstrated: at n=33 (X>2³¹) the demo prime q=2³¹−1 ACCEPTS its own wrap-around
  alias (8/8), which BIG_Q=2⁶¹−1 rejects 8/8 — the concrete reason n≳31 was blocked.
  Correct fix = larger-**characteristic** prime (implemented, `big`=2⁶¹−1, sound to
  n≲60; honest+6 cheats pass over both q and big; SZ error 3.3e−8→3.1e−17 as a
  bonus). New cheat class `u_alias` + `forge_alias`; selftest cases 5–7. CLOSED row
  added; canonical-doc field caveat corrected. Cost note: q>2³² needs a Mersenne
  mulmod to keep numpy speed (object arrays ~10–50× slower; structural skeleton
  unchanged).
- Compressed-prover step 6 (S497): propagated the S496 structured Õ(x) wiring
  prover INTO the real compressed-√x-state π(x) chain — `structured` threaded
  through `compressed_layer.py`'s reduction stack to both delegated wirings
  (`inner_verify_div`, small division + large affine `accept_rem=p−1`).
  Bit-identical transcript (selftest 15: identical accept/reject + comm +
  claimed π(x), structured vs dense, over primitives/layer/chain incl. cheats).
  New `--prover-bench`: wiring prover flattens `p²→p` (`dense/p²`→~3.3, ratio
  1.25→46× over p∈[7,1021], both wirings). **Whole compressed delegated π(x)
  chain prover now Õ(x)**, verifier still Õ(√x) — open item 1 complete on the
  real compressed chain. Default stays dense (artifacts verbatim).
- Compressed-prover step 5 (S496): structured chain PROVER —
  `chain_layer_reduce_structured` replaces the dense backward sweep's
  five `2^{2ℓ}`-size tables with two O(p)-size phase sum-checks (fold the
  σ-sum into G(y) via the sparse 4-nonzero R_j; phase 2 over R*(w)).
  Per chain layer O(p²)→O(p); whole chain Õ(x^{3/2})→**Õ(x)**.
  Bit-identical transcript (rng-drawn challenges + same multilinear
  round polys), so a verbatim opt-in drop-in. `--prover-bench`:
  dense/p²→const (Θ(p²)), ratio 1.3×→43.9× over p∈[7,1021]. **Open
  item 1 complete** (Õ(√x) verifier + Õ(x) prover, unconditional).
- Compressed-prover step 4 (S495): delegated BOTH wirings (small
  `w_div_eval`, large `waff`) from O(nb·p) to O(nb·log p) via the
  carry-chain — verifier of the compressed π(x) chain is now cleanly
  Õ(√x), no p-linear term (`--wiring-bench` shows delegated∼log p,
  automaton∼p; crossovers p≈251/1021). Generalized
  `inner_verify_div(accept_rem)`; `verify_waff_value` folds the
  comparator as-is + delegates the pure affine-image relation. 13 cheat
  classes (incl. forge + bad-chain liars). New bottleneck: the chain
  prover's dense O(nb·p²)/layer tables (Õ(x^{3/2}); reducible to Õ(x)).
- Compressed-prover step 3 (S494): small-side layer (value→value
  division wiring; band gate; empty-band guard) + the full two-sided
  K-layer chain verifying π(x)==sieve over the √x state (prover Õ(x));
  self-consistent-liar + 11 other cheat classes rejected (open item 1
  core complete).
- Compressed-prover step 2 (S493): compressed Lucy DP substrate +
  one large-side layer verifier with affine/trace dispatch; outputs
  equal true MLEs of S_{i-1}; 5 cheat classes rejected (open item 1).
- Compressed-prover step 1 (S492): batched multiplication-trace
  sum-check certifying ⌊X/(dp)⌋ lookups, polylog verifier, Õ(√x)
  prover, 6 cheat classes rejected (see open item 1).
- Verification dichotomy + 3 working protocols (see table above).
- Width-spectrum theory: transfer-DP MLE evaluation (T1); division
  wiring Nerode width = 2p+1 TIGHT; χ_P profile = maximal plateau +
  admissibility cliff; one-shot verification dead (constructibility +
  width); E2.1 middle-cut refinement (width maximal, rank deficiency
  linear-algebraic).
- Census law: D_k(x) = A_k + 1 saturation (k=8,16 closed; k=32
  converging); A_k exact to k=80 via pruned DFS; entropy conjecture
  A_k = e^{(1+o(1))π(k)}.
- Scope corrections: S18/S48 sumcheck closures (computation ≠
  verification); GKR-over-AKS gives in-principle polylog unconditional
  verifier (literature corollary).
- Old multi-stage framework archived: `archive/framework_v1/`
  (run_v1.sh, FOCUS_QUEUE.md, commit_state_final, run/verify state).
  Sessions 1-490 syntheses: `archive/sessions/`.

## NEXT ACTION (single, concrete)

The S501 profile is in: `run_chain`'s wall is **`fmul` (56%, 2.3M calls,
~30-elem arrays) + `sumcheck`'s Python loop (27%)** — op-count-bound, exactly as
S500 predicted. The fattest cleanly-batchable kernel is the **trace zero-test
(`verify_constraints`, 53% of wall in K independent per-layer sum-checks)**.

**NEXT: build the batched trace zero-test** — replace the K per-layer
`verify_constraints` calls (each a deg-3 sum-check over its own nb-cube witness,
fresh tau/alpha) with ONE random-linear-combination sum-check over the stacked
`(K·2^nb)` cube. Each layer's witness `build_witness(x, pₗ, nb)` is independent
of the carried chain claim, so this is a textbook batched sum-check (combine the
K zero-tests by powers of a random β, share one tau across the stacked cube).
Deliverable: a new experiment (e.g. `experiments/.../batched_trace.py`, or a
`--batch-trace` path in `compressed_layer.py`) with `--selftest` proving
bit-equivalent soundness (honest accepts; every per-layer cheat — wrong
quotient, corrupted bit, wrong claimed C — still rejected) and a `--bench`
showing the fmul COUNT drop ~K-fold and the per-fmul array WIDTH rise ~K-fold,
measured on BOTH the object and `--fast-big` uint64 paths (the latter is where
the width win should flip S500's "fast path loses on the chain" to a net gain).
Target metric: end-to-end `run_chain` wall at n=16 BIG_Q down from ~15 s.

Note the verifier-cost caveat: the stacked sum-check has nb+⌈log₂K⌉ rounds and
the layers share challenges, so re-confirm the verifier stays Õ(√x) and
soundness error stays ~(deg·rounds)/q. If the trace batch lands, the second
target is the wiring delegation (30%, 76% of sumcheck calls) via defer-and-batch
(collect all `(r_v,r_u,claim)` obligations, discharge in one inner protocol).

Fallback headline (if batching proves hard): a **slow but sound** item-5 run —
the object chain already verifies `π(2ⁿ)==sieve` over `BIG_Q` end-to-end; n=22–24
is ~7–20 min of pure object arithmetic. That is the headline application of open
item 1 (exact π at x≈10⁷ behind an Õ(√x) verifier over a sound prime), just not
yet fast.

Also still open (multi-cycle, lower priority): a real outer poly-commitment for
the leaf/`S_0` openings (currently `mle_eval` stand-ins).
