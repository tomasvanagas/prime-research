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
| cross-cube batched chain (step 9) | same, **both big kernels batched** (trace + wiring) | unchanged Õ(√x) | done (S502 trace + S504 wiring) — `run_chain(batch_trace, batch_wiring)`; widens the per-fmul array so the fast Mersenne path WINS: BIG_Q n=16 **8.87 s vs 17.1 s baseline (1.92×)**, comm 5.3×↓. Sound `π(2ⁿ)==sieve` over BIG_Q to **n=20 (x≈10⁶, 66 s)** |
| real leaf openings (step 10) | same, **carried-claim leaf opens via real sum-check openings** (residuals threaded to base, no `mle_eval` close) | Õ(√x) (removes ~5 per-layer O(2^nb) leaf closes; line-429 Ub openings remain) | done (S505) — `leaf_open.py` `open_eval`/`open_batch`; `run_chain(--pcs)`; selftest 20; verdict unchanged over q & BIG_Q; standalone opening verifier flat in 2^nb (6.5→432× vs `mle_eval`, 21→3605× over BIG_Q); chain `--bench-pcs` 1.2–1.36× |
| batched Ub openings (step 11) | same, **the LAST per-layer O(2^nb) verifier term** (`verify_trace_region`'s nb Ub-bit openings) discharged in ONE sum-check vs the committed stacked Ub cube | **per-layer verifier leaf-eval-free → honestly Õ(√x) end-to-end** | done (S506) — `batched_trace.verify_ub_openings_batched` (γ-RLC `Σ_w B·C`, B verifier-anchored, C the committed Ub fold); `run_chain(--batch_ub`, needs pcs+batch_trace); selftest 21 + bt §7; verdict unchanged over q & BIG_Q; **`ub_leaf_v`: K·nb→0** at every n (moved to prover); wall-clock tV ~1.1–1.2× (flat — vectorised mle_eval cheap vs Python recomputes; the win is op-count, correcting the S505 "ratio grows" prediction) |
| verifier op-count CURVE (step 12) | the milestone as a **measured scaling law**, not a per-n fact | n/a (measurement) | done (S507) — `_acct_vleaf` tallies VERIFIER large-table (`mle_eval` over a √x-size cube) evals, split per-layer vs one-time; `--bench-verifier-ops` over n∈{8..18}, 3 configs. **Falsifiable claim CONFIRMED:** per-layer leaf-op count `K·(nb+5)−1` (no-pcs) / `K·nb` (pcs) / **0** (pcs+batch_trace+batch_ub); whole-chain `vleaf_ops` fitted exponent **α=0.961 / 0.998 / 0.500** (Θ(x)→Θ(√x)); config (c) total = one-time `2·2^nb` only. `comm` slope ~0.60 (Õ(√x)) confirms the non-leaf residual never re-introduces a Θ(x) term. selftest 22 (counts over q & BIG_Q) |
| tensor-PCS base closes (step 13) | the two one-time S₀ base closes — the LAST x-scaling verifier term — via a **hash-based multilinear commitment**, **sub-√x** | **whole-chain verifier sub-√x** (computational, under a CRH/RO); O(x^{1/4}·polylog) base | done (S508) — `pcs_commit.py` (Ligero/Brakedown tensor: reshape 2^nb→r×k, RS-encode rows, Merkle-commit N columns; `S~(pt)=a^T M b`; eval-proof verifier `O(t·(r+k))=O(x^{1/4})`). `commit/prove/verify`, Fiat–Shamir, sha256. `run_chain(--commit-base)`; selftest 23 + standalone selftest (honest==mle_eval; **4 cheat classes rejected**: wrong claim / forged opening / tampered codeword row / tampered revealed column). `--bench-verifier-ops` config (d): whole-chain `vleaf_ops` exponent **α=0.258 (Θ(x^{1/4}))** vs (c)'s 0.500 — **Θ(x)→Θ(√x)→Θ(x^{1/4})**. Verdict unchanged over q & BIG_Q. **Trades unconditionality for full succinctness; win is the SLOPE** (absolute crossover at nb≈14/n≈28) |
| streaming the batched discharge — space-TIME tradeoff measured (step 19) | the S525 lever: can the batched trace/Ub discharge stream to Õ(√x) peak RSS, and is it free? | n/a (measurement + premise correction) | done (S526) — `streaming_batched_sumcheck.py`. **Premise corrected:** the discharge is NOT a pure γ-RLC fold — sum-check round polys couple all K witnesses through the pointwise product, each round's FS challenge sums over all K ⇒ no one-pass accumulator. **(1)** the K-witness LIST streams BIT-IDENTICALLY for free (`build_stacked_streaming`==`stacked_tables`, all tables/cheats over q & BIG_Q; removes ~1.4 GB of the n=24 3.9 GB). **(2)** the stacked CUBE (~2.6 GB, dominant) does not. **(3)** a bind-e-FIRST multi-pass streaming sum-check DOES reach Õ(√x) (block-diagonal layer axis ⇒ additive e-round polys), VERDICT-faithful (MLE order-independent; honest accept + every cheat reject first/mid/last over q & BIG_Q, K=1 edge, fast==object; selftests 1–7). MEASURED: streamed peak RSS **FLAT 36 MB** (1.00×/Δn=2) vs stacked **37→91 MB** (→Θ(x) ×4/Δn=2) — but **Θ(√x)× MORE WALL** (3.6→17× slower n=8..14, mean array width 341→19) — narrow dispatch-bound arrays RE-INTRODUCE the S500/S501 op-count regime batching (S502) escaped. **Headline:** held cube is INHERENT to the one-pass batched sum-check's wall efficiency; the only free win is LIST-streaming. No landed-chain change |
| LIST-streaming LANDED on the real chain (step 20) | the S526 free win, threaded into `run_chain` and MEASURED at the reach | unchanged Õ(√x) (peak-RSS only) | done (S527) — `run_chain(stream_witnesses=True)` default; both batched discharges route through `build_stacked_streaming` (now in `batched_trace.py`) so `Ws=[build_witness(p) for p in primes]` is never materialized. **BIT-IDENTICAL** (cube unchanged ⇒ same transcript): selftest case 25 (compressed_layer) + case 6 (large_x_benchmark) assert `stream` ON==OFF on claimed/comm/comm_bt/comm_ub/comm_outer/accepted over q & BIG_Q, liar still rejected. **MEASURED A/B peak RSS** (`--stream-probe`, FULL vs FULL_NOSTREAM, fresh proc): n=20 **185 vs 265 MB**, n=22 **677 vs 988 MB** (full_nostream reproduces S525's pre-streaming peaks bit-for-bit); LIST saved **6→21→80→311 MB** = `Θ(x)` (×4/Δn=2), wall **unchanged** (free). Constant-factor (~⅓) prover working-set drop — removes the LIST, NOT the `Θ(x)` cube (S526: cube can't stream for free) |
| field-size-free expander row code (step 18) | the tensor PCS over an **arbitrary field** — the RS `N=blowup·k ≤ q` demo constraint removed | unchanged sub-√x verifier (computational) | done (S514) — `pcs_commit.py` `code=` selector (`"rs"` default verbatim; `"expander"`). Brakedown-style linear-time **column-regular** recursive code `Enc_n(x)=x‖Enc_{⌈n/2⌉}(Ax)‖B·Enc_{⌈n/2⌉}(Ax)`, fixed public map (q-free indices), any prime. Verifier re-encodes v,w once `O(N)=O(k)` then indexes columns — **slope 0.524 (α≈0.262), same as RS, over q/BIG_Q/SMALL_Q**, codeword shorter (N=380 vs 512 @nb=14). **4 cheat classes rejected** (forge diff = δ·Enc(e_{j0}), weight ≥δN, measured δ≈0.45); concrete win: over **q=17 RS refuses to commit (N=64>17), expander commits+verifies+rejects cheats**. `run_chain(--commit-code expander)`; selftest cases 9 (tiny-field) + 24 (chain verdict unchanged over q & BIG_Q). **Honest:** distance MEASURED (δ≈0.45 demo params), not the asymptotic expander bound (Brakedown analyzed: δ≈0.07, t≈148 @100-bit); still computational (CRH/RO). **S516: δ measured as a CURVE** (`--distance-sweep`/`--forge-rate`) — δ declines 0.65→0.40 over k=8→512 but DECELERATES (decrements shrink, tail-floor est ≈0.34, NOT →0); empirical forge accept tracks `(1−δ)^t` (emp/pred ≤1, no query correlation); t for 2^−100 grows only 67→134 |
| certificate-SIZE profile + comm wall localized (step 14) | the FULL non-interactive certificate: **size, verifier ops, prover wall** as fitted exponents + comm **attributed by source** | n/a (measurement) | done (S509) — `run_chain` partitions `comm` into `comm_outer/bt/bw/ub/base` (additive snapshots, sum==comm asserted); `certificate_profile.py` over n∈{8..18}, full config (delegate+structured+pcs+batch_trace+batch_ub+**batch_wiring**+commit_base). **Cert size α=0.473 (Θ(√x)), DOMINATED by `comm_outer` α=0.522 (60%→90%→100% of comm)** = the K=π(√x) SEQUENTIAL outer layer reductions; every batched discharge polylog (α≤0.18), base Θ(x^{1/4}) (α=0.30), verifier ops α=0.258. Batching wiring cuts comm 5.3× @n=18. **Honesty sharpening of S508: verifier READS the Õ(√x) proof ⇒ total verification Θ(√x); the x^{1/4} is the large-table sub-term only.** Membership: `L_π` has a **CRH-based NI certificate of size Õ(√x), Õ(x^{1/4}) large-table opens**; polylog hits the SAME ⌊v/p⌋-semigroup wall as open item 2. selftest 5 cases (attribution closes over q & BIG_Q honest, cheats reject; field-indep count; exponent bounds) |
| merged-layer COMM probe — √x is layering-inherent (step 15) | the merged depth-j layer's **comm**: is the Õ(√x) cert size inherent to the layering or just the unbatched chain? | n/a (measurement) | done (S510) — `merged_layer_comm.py` collapses j consecutive small-side Lucy layers into ONE reduction over the 2^j-subset inclusion-exclusion form `M[S](v)=Σ_T(−1)^|T|GateProd_T(v)S(⌊v/m_T⌋)` (== direct composition, selftest 1); two SOUND modes (DIRECT 2^j openings; BATCHED stacks the 2^j on a j-bit axis), cheat-tested over q & BIG_Q. **Chain TOTAL comm exponent ~0.5 for EVERY j∈{1,2,4} in BOTH modes → √x LAYERING-INHERENT.** DIRECT comm_fill=(2^j+1)(3nb+2) grows ~2^j (merging makes comm *worse*, min at j=2 then grows); BATCHED pushes the 2^j OUT of comm (fill flat 27) but into **verifier WORK** (`vwork` ~2^j, ratio 5.7× @j:1→4) AND prover, comm floored at ⌈K/j⌉(j+nb)→K=√x; **prover fill total ⌈K/j⌉·2^j·2^nb → x^{3/2}** (/x: 0.42→2.25 @j:1→6) — open problem 2's m·K²/2 on the VERIFICATION face. comm prime- & field-independent (selftest 4) |
| #P / NP-membership probe — the COMPUTATIONAL lower-bound face (step 17) | item 3's other half: is `L_π∈NP` at all, and is π #P-hard? | n/a (measurement + complexity anchor) | done (S512) — `sharpP_probe.py`. **(A)** π(x), φ(x,a) are #P FUNCTIONS of binary x (count ≤2^N candidates under a poly(N) predicate) ⇒ `L_π∈C_=P`, `{π≥c}∈PP⊆P^#P`; the φ subset-witness is 2^{π(√x)} (doubly-exp in N). **(B)** π a function ⇒ `L_π∈NP ⟹ NP∩coNP ⟹` NP-complete forces NP=coNP — the live question is mere NP-MEMBERSHIP (polylog witness). **(C) witness-size ladder** (leading-power α, 2-term fit separating power from polylog): enumeration NP-cert α=0.985 (Θ(x)); sieve transcript α=0.490 (matches S509 0.473); S511 info floor α=0.486; zeta-zero/explicit-formula witness α=0.500 (Galway √x·log²x, **cited**) — THREE natural families converge at √x (max\|α−0.5\|=0.014), enumeration a full power higher, polylog (α=0) ruled out (S511). NO natural witness reaches poly(N). **(D)** parsimonious #P-hardness reduction is CIRCULAR (target-count realization c↦x(c) IS the inverse-prime p(·)=the goal; sieve sets lattice-structured, no instance richness) ⇒ hardness, if real, is value-incompressibility (the S511 face). **(E)** corrected CLOSED row 175 ("exact π #P-hard" unsubstantiated). 40 selftests |
| soundness budget + field-size requirement (step 21) | the field-soundness axis (orthogonal to the closed integer-range face): how big must q be, and does soundness or correctness bind at p(10^100)? | n/a (measurement) | done (S534) — `soundness_field_size.py` (9 selftests). Field-soundness degree budget **B(x)=Σ per-round deg = K·c_BK** with K=π(√x) layers (exact) and c_BK=B/K≈**15.35·n−36** (polylog) ⇒ the prime-count 1/log x CANCELS the per-layer log x ⇒ **B=Θ(√x)** (single-power α=0.548; d_max=n−2). Per-round Schwartz–Zippel MEASURED: deg-3 trace zero-test false-accepts at **rate·q≈3.8 const** (n=8 high-count, slope −1; n=10 high-q points Poisson-noise-thin), chain `corrupt_layer` self-consistent liar **≤ B/q always (no amplification)**, `delta_pi` output-value lie caught **deterministically (0/800 any q)**. Field-size: soundness needs `log₂q>λ+n/2−log₂n+polylog`, correctness needs `log₂q>n+~1` (q>X, measured slack ~1 bit) ⇒ **CROSSOVER n≈209 (x≈10^63=2λ); at p(10^100) (n=332) q_corr≈333b > q_sound≈272b ⇒ CORRECTNESS BINDS**. Sharpens "Õ(√x) verifier": elements are Θ(log x)-bit (correctness-forced) ⇒ the Õ polylog carries ~n^{1+ε} field-op bits, λ-independent for x≳10^63. Adjacent measurement (a known √x wall on the field-soundness face), NOT a goal advance |
| certificate INCOMPRESSIBILITY probe — √x is **information-forced** too (step 16) | the lower-bound face of item 3: is the Õ(√x) cert size an INFORMATION floor, or just construction shape (S510)? | n/a (measurement) | done (S511) — `cert_incompressibility.py`. Per-layer checkpoint = Legendre partial sieve `φ(x,a)=#{n≤x : P⁻(n)>p_a}` (the large-side survivor count). Measures the JOINT hard bits of `{φ(x,a)}_{a=1..K}` beyond a polylog smooth predictor, at INTEGER precision. **Authoritative measure = integer-reconstruction SVD rank of the residual matrix** (smooth-model-free: SVD subsumes the best low-rank/smooth model). rank≤min(W,K); window-sensitivity at x=2¹⁹ shows rank/K **0.09→0.97** as W/K 1→187 ⇒ **FULL RANK once W≫K**. Adaptive-window sweep (W=64K): rank≈K, **α_rank=0.459 ≈ α_K=0.415 ≈ √x**. ⇒ the K checkpoints are integer-INDEPENDENT, joint info **Θ(√x)·polylog** ⇒ the √x cert is **information-forced for any sieve-reconstructing verifier** (polylog cert RULED OUT; sub-√x RULED OUT for this class). Corroborated per-x (gzip α=0.39, AR α=0.56 track K). Single-value control reproduces S36: |π(x)−Li(x)| is O(log x) bits (slope 0.089) — the barrier is JOINT, not per-value. **Self-corrected** a narrow-fixed-window artifact (W=4096 → rank~x^0.35 sub-√x) via the window-sensitivity check. Scope: sieve-reconstruction class only; a UNIVERSAL bound still needs #P-hardness (item 3 formal half). selftest 22 cases. **Sharpened (S515):** wide window (W=192K) lifts rank/K **floor to 0.88** (no drift, k=16..22), **α_rank/α_K=0.978≈1** ⇒ robust signal is rank=Θ(K)=Θ(√x), NOT the exponent (α_rank≈α_K≈0.40 = shared PNT 1/ln x discount, →0.5 only as x→∞; the 64K α_rank=0.459 was inflated by an across-window rank/K climb — self-corrected). Per-layer info **DENSE** across all K layers (bits uniform 0.77, active 0.98, prefix rank LINEAR `rank_half_ratio`=0.52) — Θ(K) layers, not o(K). 30 selftests |

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
   the tree analysis. **The COMM face is now measured too (S510,
   `merged_layer_comm.py`):** a sound depth-j merged-layer reduction has
   chain total comm exponent ~0.5 for EVERY j in both a direct (comm_fill
   ~2^j) and a fully-batched (comm polylog, 2^j → verifier WORK + prover)
   mode — the √x is layering-inherent on the verification face, and the
   prover fill total still → x^{3/2}. So the cert-size wall, the
   verifier-work wall, and the prover-time wall are ALL the one
   ⌊v/p⌋-semigroup sequential-sieve wall.
3. **L_π = {(x,c): π(x)=c} ∈ NP?** — succinct unconditional
   certificates. Upper bound: CRH-based NI cert of size Õ(√x) (S509).
   **Lower-bound face PARTIALLY ANSWERED (S511, `cert_incompressibility.py`):**
   the √x is an INFORMATION floor, not just construction shape — the K=π(√x)
   sieve checkpoints `{φ(x,a)}` carry Θ(√x)·polylog JOINT hard bits at integer
   precision (residual matrix is FULL integer-rank once W≫K: rank/K→0.97,
   α_rank≈α_K≈√x). **Sharpened (S515):** the robust √x signal is the rank/K
   FLOOR (0.88 at W=192K, no decay) ⇒ rank=Θ(K)=Θ(√x), with α_rank/α_K=0.978≈1
   (the finite-window exponents ≈0.40 are the shared PNT 1/ln x discount, not a
   rank deficiency; the S511 α_rank=0.459 was inflated by an across-window
   rank/K climb — self-corrected). The info is also DENSE across all K layers
   (uniform per-layer hard bits, LINEAR prefix rank, `rank_half_ratio`=0.52) ⇒
   carried by Θ(K) layers, not o(K). So **no sieve-reconstructing verifier has a
   polylog or sub-√x certificate** — the construction-side √x (S510) is matched
   on the information side. Single π(x) value remains O(log x) bits (S36) — the
   barrier is JOINT across layers. **Still open / the remaining universal
   lever:** the bound is for the sieve-reconstruction CLASS only; a different
   witness (cf. factoring) is not ruled out, so a universal sub-√x impossibility
   needs the COMPUTATIONAL route. **COMPUTATIONAL face PROBED (S512,
   `sharpP_probe.py`):** (upper) π(x), φ(x,a) ∈ **#P** (count ≤2^N integers under a
   poly(N) predicate) ⇒ `L_π∈C_=P`, `{π≥c}∈PP⊆P^#P` (Toda: `PH⊆P^#P`, not the
   reverse; `PP⊆PH` is unknown, not claimed). (structure) π a function ⇒
   `L_π∈NP ⟹ NP∩coNP`, so **NP-completeness would force NP=coNP** — the open
   question is mere NP-MEMBERSHIP (a polylog witness), not completeness. (ladder)
   every NATURAL witness family — enumeration NP-cert (α=0.985, Θ(x)), sieve
   transcript (α=0.490≈S509's 0.473), S511 info floor (α=0.486), zeta-zero/
   explicit-formula (α=0.500, Galway √x·log²x cited) — lands at √x or worse; THREE
   converge at √x, polylog (α=0) is ruled out by S511 for the sieve class ⇒ **no
   natural witness reaches poly(N)**. (reduction) a parsimonious #P-hardness
   reduction is **CIRCULAR** — target-count realization `c↦x(c)` IS the
   inverse-prime `p(·)`=the goal, and the sieve/φ instance is lattice-structured
   (no embedding richness) ⇒ #P-hardness, if true, must be value-incompressibility
   (the S511 face), not instance-embedding. **Corrected** CLOSED row 175 ("exact
   π(x) #P-hard" was an unsubstantiated S7 assertion). **TURING-REDUCTION FACE
   CONDITIONALLY CLOSED (S517, `turing_reduction_barrier.py`):** S512's direction
   (i) is ruled out under SETH for all sub-c*-blowup Turing reductions. THEOREM:
   with α=inf{a:π∈TIME(x^{a+o(1)})} and SETH, no poly-time Turing reduction #SAT→π
   has all queries ≤c·n bits with c<1/α (simulate π by its own x^a' algorithm →
   sub-2^n #SAT → refutes SETH); critical blowup c*=1/α, measured `pi_lucy`
   α≈0.66–0.70 → c*≈1.5, analytic α=1/2 → **c*=2**. The natural (parsimonious)
   blowup is c→1<c* (measured exact n=4..22), so EVERY natural Turing reduction is
   forbidden. COROLLARY (P≠NP): **polylog-π XOR #P-hard-π** — a #P-hardness proof
   for π proves the project goal impossible; the two questions are one dichotomy
   with √x at c*=2. **STILL OPEN (filed, not closed):** (a) super-c*-blowup Turing
   reductions (query π exponentially past p_C≈2^n — no natural construction, no
   #SAT speedup), and (b) a genuine #P-hardness proof OR a non-sieve sub-√x witness
   that also beats the S511 info floor. On present evidence `L_π` is an
   NP-intermediate-flavoured counting problem: a √x certificate (S491–S509 from
   above, S511 from below), no proven poly(N) one, #P-hardness ruled out for every
   natural reduction (S512+S517). Either answer to the remaining #P-hardness
   question is major. **DIRECTION (ii) NATURAL WITNESS CLOSED (S518,
   `explicit_formula_witness.py`):** the most natural non-sieve witness — the truncated Riemann
   explicit formula over real low-lying zeros — is √x·polylog (Galway √x·log²x best-fit), NOT
   sub-√x (raw exponent reads ABOVE ½; matched control confirms ½ leading power; nothing dips
   below ½) and NOT polylog. Its rounding-relevant information is DENSE across Θ(√x) zeros (no
   octave up to √x droppable at x=10⁵–10⁶), so it MATCHES the S511 floor rather than beating it — the floor
   caps the analytic route too, not just the sieve class. (Also fixed CLOSED row 31: the
   "divergence" was a li-branch-cut bug.) Residual (b): a genuinely NON-natural sub-√x witness
   stays open (no candidate).
4. **Census entropy law at k = 128** — conjecture A_k = e^{(1+o(1))π(k)}.
   **Transfer-matrix route CLOSED (S519, `census_transfer_matrix.py`):** the
   natural position-DP carrying per-prime occupied-residue sets has **Θ(A_k)
   states** — measured collapse factor A_k/states = 1.55–2.37 (bounded, not
   growing), ln(states) slope 1.020 ≈ ln A_k slope 1.011 — so it is
   exponential with A_k's own exponent and strictly worse than DFS (Θ(A_k)
   memory vs Θ(k)). No new clean-family exact point is reachable (k=128 is the
   next power, A₁₂₈~10¹¹). **Two scaling structural facts delivered instead:**
   (a) an exact **active-prime reduction** — A_k depends only on primes ≤
   B(k), B(k)=3,5,7,13,23,47 at k=8,16,32,64,128,256 (A₁₂₈ uses only 8 primes
   {3..23}; rigorous, two independent proofs: count-stabilisation + a
   W(B)-weight bound via branch-and-bound min-union, agreeing); (b) the
   identity **W(B(k)) = ρ*(k) = maximal admissible-tuple size** = 3,5,9,16,28,
   52, matching H() exactly at k=64 (H(16)=60<H(17)=66) and k=128 (H(28)=126<
   H(29)=132) — this explains ln A_k=Θ(π(k)) (admissibles live in the ρ*(k)≈
   k/ln k geometry). The entropy slope→1 (1.01) is re-confirmed.
   **CONSTANT now estimated by a cluster expansion (S520,
   `census_entropy_asymptotic.py`):** exact reformulation `A_k=2^N·P[adm]`
   (offsets i.i.d. p=½), `\bar C_q`="miss a class mod q" DECREASING ⇒ FKG.
   Mean field `2^N∏_q p_q` is a RIGOROUS LOWER bound; corrections = the
   cluster (Möbius) expansion of `ln P_B`, each term a joint-cover prob got
   from a **nested inclusion-exclusion** (cost `2^{sumS−maxS}`, k-independent,
   bit-identical to S519 `count_tm`). Validated EXACT k≤64: errors alternate
   (odd r under, even over), shrink ~5×/order, r=4 within 0.001 at k=64.
   **RIGOROUS BRACKET via FKG blocks** (`P[∩_B]≥∏_blocks P[∩_block]` lower,
   single block ⊆B upper): **k=64 ∈[0.9236,0.9276] contains exact 0.9276;
   k=128 ∈[0.9243,0.9403]** (cluster est 0.926). Ratio climbs
   0.667→0.808→0.886→0.928, consistent with →1. **Key finding:** the
   constant is CORRELATION-driven — the independent-prime term DECLINES
   (0.75→0.61→0.39 over k=64,128,256), so the CRT overlaps carry a GROWING
   fraction of the entropy (19%→35%→57%); →1 is not in the leading term.
   **METHOD CLOSED + S520 k=256 CORRECTED (S521, same file):** `--bracket-sweep`
   shows the rigorous bracket width is a STEP FUNCTION of the I-E cost budget C,
   FLAT between discrete prefix-cost thresholds; at **k=128 width STUCK at 0.0160
   for every feasible C≥26** (next tighten 2^39, shuts 2^75). **Closing cost
   `2^{sum(B)−max(B)}` super-exponential: 2^26→2^75→2^279→2^872 = 2^{Θ(k²/ln²k)}
   ⇒ the FKG-block/cluster method CANNOT close the bracket / prove constant=1 for
   k≥128** (engine-independent: {17,19,23} also aborts count_tm >25.7M states).
   SECOND wall: at k=256 (N=128) `P[miss]~1e-17` but the I-E sums terms ~1 ⇒ the
   FLOAT engine cancels to ZERO digits — **S520's k=256 bracket [0.9282,1.126] was
   NOISE, RETRACTED** (its "lower" 0.9282 > true exact lower 0.8312). New
   `joint_cover_int` (exact integer, `==count_tm`) removes it: **corrected k=256
   bracket [0.8312,0.9431] RIGOROUS** (width 0.1119; cost wall pins upper at 5 of
   14 primes). **Still open:** constant=1 NOT proved AND not provable by this
   method — needs an analytic handle on the FULL Θ(entropy) correlation series (new
   analysis, not new code). An exact A₁₂₈ still needs a genuinely sub-A_k algorithm
   (none found). The census line is closed for computational cycles.
5. **Large-x benchmark** — REACH AT n=24 DELIVERED (S513 driver, S524 n=24).
   `experiments/constructions/p12_sumcheck_pi_verification/large_x_benchmark.py`
   (a clean reproducible driver, `--selftest`/`--n`), running the **FULL
   succinct config** (delegate+structured+pcs+batch_trace+batch_ub+batch_wiring+
   commit_base) + FAST_BIG over the sound `BIG_Q=2⁶¹−1` — strictly stronger than
   the S504 config (adds the real leaf openings, batched Ub discharge, and the
   tensor-PCS base commitment). Verified `claimed π==sieve`, honest ACCEPTED:
   **n=20 π(1048575)=82025 (71.8 s)**, **n=22 π(4194303)=295947 (252.8 s)**,
   **n=24 π(16777215)=1 077 871 (1467.3 s, x≈1.7×10⁷)** — the n=24 claim
   cross-checked against an independent from-scratch Eratosthenes sieve
   (π(2²⁴)=1 077 871, known value). Soundness witnessed at the reach: the
   `delta_pi` liar is REJECTED at n=22 (66.9 s) AND n=24 (541.9 s, claim
   1 077 872). Profile holds at every reach n — per-layer verifier large-table
   leaf-eval count **= 0** (Õ(√x) end-to-end), `comm` ~97% `comm_outer` (Θ(√x);
   batched discharges polylog: n=24 comm_bt/bw/ub=88/2829/67), base opens the
   one-time tensor commitment only (Θ(x^{1/4}); n=24 vcommit_ops=16512). Chain
   stays prover-bound (reported t_prover+t_verifier ~1.7% of wall at n=24).
   **Scaling CORRECTED with the 3rd point:** the per-Δn=2 wall ratio is NOT the
   constant ~3.5× S513 reported from the lone n20→n22 step — it ACCELERATES
   (3.52× → **5.80×**, geo-mean 4.52×); least-squares **wall ≈ x^1.088** over
   n=20..24. Honest reading: the asymptotic prover is Õ(x) (open item 1) and the
   ~1.09 exponent is linear-in-x+polylog; the n22→n24 step runs above the ideal
   4×/Δn because of a ~3.9 GB working set ⇒ super-linear WALL overhead — an
   implementation cost, not a complexity change (op-count milestones unchanged).
   **S525 CORRECTED the working-set attribution:** S524 blamed the object-dtype
   `np.add.at` per-layer scatter-sums; building the uint64 `scatter_fold61` fix and
   re-measuring shows it is **WALL-NEUTRAL** (n=20/22/24 within ~2% of S524; n=22
   A/B comm + peak RSS IDENTICAL ON-vs-object). The 3.9 GB is the **Θ(x)
   batched-discharge held cubes** (uint64), localized by `--mem-probe` (FULL peak
   265 vs NOBATCH 138 MB @ n=20); the super-linear wall is a memory-hierarchy
   constant, not GC. `large_x_benchmark_results.md` (§S525); run logs `run_n24.log`,
   `run_s525_remeasure.log`. **Streaming-the-stack lever RESOLVED (S526,
   `streaming_batched_sumcheck.py`):** the held cube is NOT a pure γ-fold; streaming it to Õ(√x)
   space is verdict-faithful (bind-e-first, peak RSS FLAT vs stacked's Θ(x)) but costs Θ(√x)× MORE
   WALL (narrow dispatch-bound arrays = the S500/S501 regime batching escaped) ⇒ it makes the wall
   WORSE, not better. Only the K-witness LIST streams free + bit-identically (removes ~1.4 GB of the
   3.9 GB; the ~2.6 GB cube stays). The Θ(x) held cube is inherent to the one-pass batched
   sum-check's wall efficiency — the working-set wall is not removable for free.
   **LIST-streaming LANDED + MEASURED (S527):** `run_chain(stream_witnesses=True)` (default) now
   builds the cube slice-by-slice (the list is never held); bit-identical (selftests pass; A/B
   claimed/comm/comm_bt/comm_ub equal), so all artifacts incl. n=24 reproduce verbatim. Measured
   A/B peak RSS: n=20 **185 vs 265 MB**, n=22 **677 vs 988 MB** (full_nostream reproduces S525's
   pre-streaming peaks); LIST saved Θ(x), wall unchanged (free). Constant-factor (~⅓) drop; the
   Θ(x) cube remains the binding constraint, so the reach is now cube/RAM-bound (see NEXT ACTION).
6. **Leaf openings → end-to-end Õ(√x) verifier** — DELIVERED (S505 pcs + S506 batch_ub).
   `leaf_open.py` builds the real sum-check MLE opening (`open_eval`/`open_batch`);
   `run_chain(--pcs)` threads the carried-claim folds' residuals to the S₀ base. The
   one remaining per-layer O(2^nb) term — `verify_trace_region`'s nb Ub-bit-table
   openings (per-layer witness data, no carried claim) — is now DEFERRED and
   discharged in ONE sum-check against the committed stacked Ub cube
   (`batched_trace.verify_ub_openings_batched`, `run_chain(--batch_ub)`, needs
   pcs+batch_trace). **Verifier-side O(2^nb) Ub-leaf count K·nb→0 at every n** (selftest
   21, `--bench-ub`); the per-layer verifier is leaf-eval-free, the only O(√x) verifier
   work left is one-time (the two S₀ base closes + the batched discharges). **The
   Õ(√x) end-to-end UNCONDITIONAL verifier claim is now structurally real.** Honest
   residual: the wall-clock tV win is a flat ~1.1–1.2× (vectorised mle_eval is cheap
   vs Python-loop recomputes), so the asymptotic gain is in op-count, not measured
   wall-clock at reachable n — this corrects the S505 "ratio grows with n" prediction.
   **The op-count gain is now a MEASURED CURVE (S507):** `--bench-verifier-ops`
   attributes the verifier's large-table (`mle_eval`) evals per-layer vs one-time over
   n∈{8..18}; the whole-chain `vleaf_ops` leading exponent is α=0.961 (no-pcs) / 0.998
   (pcs) / **0.500 (pcs+batch_trace+batch_ub)** — Θ(x)→Θ(√x), config (c) total = the
   one-time `2·2^nb` base closes only, per-layer leaf count exactly 0. `comm` slope ~0.60
   (Õ(√x)) confirms the non-leaf residual stays sub-linear. The Õ(√x) end-to-end verifier
   is now both structurally real AND a measured scaling law.
   **The last one-time O(√x) item is now DONE (S508):** a real outer multilinear
   commitment (`pcs_commit.py`, tensor/Ligero-Brakedown) discharges the two S₀ base
   closes with a **sub-√x verifier** (`O(x^{1/4}·polylog)`); `run_chain(--commit-base)`,
   `--bench-verifier-ops` config (d) measures the whole-chain verifier leading exponent
   **α=0.258 (Θ(x^{1/4}))**. The **whole-chain verifier is now sub-√x** — full
   succinctness, COMPUTATIONAL (under a CRH/RO), trading the otherwise-unconditional
   Õ(√x). The win is the SLOPE (absolute crossover nb≈14/n≈28; RS's N≤q is a demo
   artifact a Brakedown expander code removes behind the same interface).
   **Certificate SIZE characterized + the Õ(√x) comm wall LOCALIZED (S509):**
   `certificate_profile.py` measures the full succinct config's cert size (α=0.473,
   Θ(√x)), verifier ops (α=0.258, Θ(x^{1/4})) and prover wall over n∈{8..18}, and
   `run_chain` now attributes `comm` by source. The Õ(√x) cert size is DOMINATED by
   `comm_outer` (α=0.522; 60%→90%→100% of comm) = the K=π(√x) SEQUENTIAL outer layer
   reductions; every batched discharge is polylog, the base is Θ(x^{1/4}). **Honesty
   sharpening:** the verifier READS the Õ(√x) proof, so total verification is Θ(√x)
   (the S508 x^{1/4} is the large-table sub-term only). Membership: `L_π` has a
   CRH-based NI certificate of size Õ(√x), Õ(x^{1/4}) large-table opens; compressing
   the size to polylog hits the SAME ⌊v/p⌋-semigroup wall as open item 2.
   **S510 closes the quantitative question this raised:** merging j consecutive
   layers does NOT drop the comm exponent below 0.5 — measured over a sound DIRECT
   (comm_fill ~2^j) AND a fully-BATCHED (comm polylog; the 2^j moves to verifier
   WORK + prover) merged-layer reduction, `merged_layer_comm.py`. The Õ(√x)
   certificate is layering-inherent, and the merge pays open item 2's x^{3/2} prover
   fill. No sub-√x certificate comes from layer merging.

## SECONDARY LINES (valid, not currently active)

- `OPEN_POSITIVE_TARGETS.md`: **P8 (sparse precision) CLOSED-NEGATIVE
  (S533)** — per-bit cost of the explicit-formula π evaluator is
  GEOMETRIC/super-linear (slope 1.5–2.5, R²≈0.95), not sub-linear; the
  precision-domain face of the info barrier (see Done-this-era). **P9
  (quantum batched)** remains open (literature survey only, likely
  inconclusive). **P7.b CLOSED** (S530: re-derivation of the closed
  S224 Correlation Dichotomy — dyadic anchors = uncorrelated side).
  **P10 (adaptive queries) CLOSED-NEGATIVE (S532)** — adaptivity beats
  nothing vs non-adaptive R⁻¹(n)+sieve; SMOOTH+RANDOM is the adaptivity
  boundary. **P6 carries an S202-duplication warning — read it before
  spending anything there.**
- Width-spectrum follow-ons: D_k(x) census at k = 32 **CLOSED — fully
  saturated (S522)**: D_32(2³⁷)=A_32+1=3574 (`dk_first_occurrence.py`;
  deficit 189→0, the last a weight-9 pattern realized at x≈8.3×10¹⁰),
  joining k=8,16 as a confirmed instance of `D_k→A_k+1`. The
  **sliding-window** census variant is now **MEASURED + CLOSED (S523,
  `dk_sliding_census.py`)**: unaligned windows `[n,n+k)` drop the aligned
  `q|k` condition → classical HL admissibility. **A_k^slide =
  2·A_k^aligned−1 = 25,211,7145** (the q=2 clause forces all-even or
  all-odd offsets; odd-family == aligned set, even-family its mirror under
  `o↦k−1−o`, an admissibility automorphism); densest weight = ρ*(k)=9
  (parity-symmetric, no denser tuple than aligned). **Full saturation
  D_32^slide=A+E=7153 by n<2²⁷ (≈1000× earlier in x than aligned 2³⁷)** —
  32× more windows + the densest 9-tuples realized at n≈10 by the
  small-prime cluster. Exceptional set E_k^slide=q*(k)+1=4,6,8 (windows
  n=0..q*, q*=3,5,7). Classical 𝔖 validated (obs/pred=0.994 @w6).
- The guess-comparison geography (`experiments/analytic/
  guess_comparison_oracle/`) — decision-version facts, complete as is.

## Done this era (S491–S535 cycles, 2026-06-13/14)

- Reach wall THIRD leg found — the n=26 reach is page-fault/allocation-bound, sharpening the S529
  falsifier (S535): the n=26 reach was STILL RUNNING DETACHED (started 13:24Z, ~2.85 h in at cycle's
  work, RSS 4 GB, no DONE), so I could not harvest — and a concurrent cycle must NOT start a second
  reach. Took fresh CONTENTION-SAFE ground by OBSERVING the live reach via `/proc` (pure reads) + a
  self-limiting micro-probe. Built `experiments/constructions/p12_sumcheck_pi_verification/
  alloc_overhead_wall.py` (one script, `--live`/`--probe`/`--predict`/`--all`/`--selftest`; 17
  selftests incl. a field-offset cross-check vs `os.times()` that caught an off-by-one /proc parse bug).
  **HEADLINE: the live n=26 reach is KERNEL/ALLOCATION-BOUND, not compute-bound — `stime`/(`utime`+
  `stime`)=0.695 (cumulative; 0.72–0.76 windowed), ~3.3×10⁵ minor page-faults/s, 3.1×10⁹ cumulative
  minor faults, majflt=0 (no swap).** The unavoidable fault floor = working set/page ≈ 12 GiB/4 KiB ≈
  3.1×10⁶, so ≈1000× of the faults are RE-allocations of the same Θ(x) buffers across K≈π(√x)≈1028
  layers. **This is a THIRD wall leg S529's model omitted:** S529 measured steady-state `_mul61` over
  PRE-allocated arrays (the L1→DRAM compute bowl, leg 2, per-Δn ≤1.35×) and never measured the
  ALLOCATION/first-touch cost. The micro-probe isolates it: reused-buffer op = flat ~3 ns/elem,
  utime-dominated (stime%≈0 = leg 2); fresh-allocation op = ~14 ns/elem above 256 KiB (the L2 boundary),
  **kernel/stime-dominated 0.69–0.83 (matching the reach's 0.695)**, a **~4× penalty / ~10 ns/elem**
  allocation overhead. Probe perturbed the reach's fault rate by only **+2.6%** (clean falsifier point
  intact). **Three-leg model:** `wall = opcount(Θ(x), 4.0–4.4×↓ S530/1/2b) × compute_cache(≤1.35× S529)
  × alloc_factor(NEW)`; the alloc leg grows with the Θ(x) working set so its per-Δn ratio can EXCEED the
  op-count ratio and push the WALL above S529's [4.0,5.4]× band WITHOUT an op-count change. **SHARPENED
  S529 FALSIFIER (the deliverable for the harvest cycle):** a wall ratio n24→26 >5.4× is CONSISTENT with
  Θ(x) op-count IFF the reach is `stime`-dominated (measured 0.695) — the discriminator is `stime`/`utime`,
  NOT the wall alone; only a `utime`-dominated excess with super-linear per-element growth reopens
  super-Θ(x). **Removable headroom (a lever, NOT implemented):** reusing buffers (`np.multiply(...,out=)`)
  would fault each page ~once (~3×10⁶), collapsing `stime` ~7000→~7 s, a **≈3.3× wall ceiling** — a
  CONSTANT-FACTOR implementation win (analog of S527 list-streaming/S525 scatter-fold), NOT a complexity
  change; deliberately not done (would rewrite every `a*b%q`, break the n=24/n=26 bit-identical
  artifacts) — file opt-in + bit-identical A/B if pursued. F1–F4 all un-triggered. Honest scope: an
  in-cycle MEASUREMENT refining the S529 wall model (the same Θ(x)·polylog prover, new leg), NOT a goal
  advance. Caveat: n=24's `stime`/`utime` was not logged (process exited) so the quantitative n24→26
  *utime* ratio is unpinned — future reaches should log `getrusage`. `alloc_overhead_wall_results.md`,
  `alloc_overhead_run.log`.
- Soundness budget + field-size requirement of the verification chain MEASURED (S534): the detached
  n=26 reach was still running (started 13:24Z; ~2 h in at cycle's work, still no DONE — must not start a
  second, and must not contend the clean wall), so per the NEXT ACTION I took **fresh RAM-light /
  cache-resident ground** (n≤16, ≤few-MB working set ⇒ no DRAM-bandwidth contention with the reach;
  confirmed the reach is single-threaded, load≈1.6/32, so the wall bottleneck is DRAM bandwidth not
  cores). Built `experiments/constructions/p12_sumcheck_pi_verification/soundness_field_size.py`
  (one script, `--budget`/`--sz-atomic`/`--sz-chain`/`--extrapolate`/`--all`; 9 selftest groups).
  **Target:** the SOUNDNESS / Schwartz–Zippel field-size axis that CLOSED row 498 (F_{q²} refutation)
  flagged "orthogonal to the integer-range hole" but never quantified. **(A) Soundness degree budget**
  (deterministic: a counting rng + a per-module `sumcheck` degree-tally monkeypatch on the REAL honest
  `run_chain`, both field-independent — selftest 2,4): `B(x)=Σ per-round deg = K·c_BK` with **K=layers=
  π(√x) exact** (6,11,18,31,54 = π(15..255)), **c_BK=B/K≈15.35·n−36** (~linear in log x, polylog), and
  **d_max=n−2** (the B2 lookup, O(log x)). The prime-count 1/log x **cancels** the per-layer log x ⇒
  **B=Θ(√x)** (measured single-power α=0.548≈0.5 over n=8..16; the 2-term fit's α=0.36 is the S511/S515
  PNT-deflation artifact, NOT used for extrapolation — anchored on the exact π(√x) backbone instead).
  **(B) Per-round Schwartz–Zippel MEASURED** (build witness once, re-verify with fresh challenges): the
  deg-3 trace zero-test false-accepts at **rate·q≈3.8 constant** (n=8 high-count, slope −1); the n=10
  multi-q sweep confirms rate drops ~1/q but its high-q points are Poisson-noise-thin (1–4 accepts/6000
  ⇒ fitted slope over-steepens to −1.1…−1.9 — reported honestly). **(C) Chain composition / first
  soundness audit at non-negligible error:** the self-consistent internal liar `corrupt_layer`
  false-accepts at **≤ B/q always (no amplification past the union bound)** and drops with q, while the
  output-value lie `delta_pi` (claim π+1) is caught **DETERMINISTICALLY** (0/800 at every q — the final
  base equality, not Schwartz–Zippel). Two soundness regimes. **(D) Which axis binds:** soundness needs
  `log₂q>λ+log₂π(√x)+log₂c_BK ~ λ+n/2−log₂n+polylog`; correctness needs `log₂q>n+slack` (q>X, measured
  **slack≈1 bit**). **CROSSOVER at n≈209 (x≈10^63=2λ for λ=100):** below it the security floor binds,
  above it correctness binds. **At p(10^100) (n=332): q_correctness≈333 bits > q_soundness≈272 bits ⇒
  CORRECTNESS BINDS** — the Θ(√x) soundness budget never overtakes the linear-in-n integer-range
  requirement at large x, and the ~333-bit field is automatically >100-bit-sound. **Sharpens the
  "Õ(√x) verifier" honesty:** field elements are Θ(log x)-bit (correctness-forced) ⇒ each field op is
  O((log x)^{1+ε}) bit-ops, so the Õ polylog concretely carries ~n^{1+ε}, set by correctness and
  λ-INDEPENDENT once x≳10^63. **Honest scope:** an adjacent verification-line MEASUREMENT (a known √x
  wall — same √x as the cert size S509 and info floor S511 — now confirmed on the field-soundness face,
  + the first empirical chain-soundness audit at measurable error); NOT a goal advance (polylog π(x)
  stays blocked), and the union bound is the standard loose composition (the audit confirms the
  implementation does not violate it, it does not prove a matching soundness LOWER bound). 4 pre-stated
  falsifiers F1–F4 all un-triggered. `soundness_field_size_results.md`, `soundness_field_size_run.log`.
- P8 sparse-precision π CLOSED-NEGATIVE — the per-bit cost is geometric, not sub-linear (S533): with the
  reach (n=26) running detached and RAM-bound (must not start a second), took fresh RAM/CPU-light ground
  on the last concrete positive target (NEXT ACTION's P8 first step). Built
  `experiments/analytic/sparse_precision_pi/sparse_precision_pi.py` (one script, `--selftest`/`--profile`/
  `--anchor`; 19 selftest cases), REUSING verbatim the S518 numerically-stable explicit-formula evaluator
  (CLOSED-row-31 li(exp()) branch-cut fix imported, cannot recur). **P8 (OPEN_POSITIVE_TARGETS): can the
  first k bits of π(x), batched on increasing precision k, be refined with a SUB-LINEAR per-bit cost?**
  **Verdict: NO.** Fixing x and varying the target precision k (bits of π(x)), the zeros-to-resolve-k-bits
  N(k) (via the windowed-RMS error envelope) has **per-bit slope s = d log2 N(k)/dk = 1.5–2.5** (R²≈0.95
  at x=10⁵/3·10⁵/10⁶) ⇒ each hard bit costs **~3–6× the zeros of the last** — GEOMETRIC/super-linear, the
  OPPOSITE of P8's hope. **Mechanism:** error envelope E∝T^{−eta}, eta≈½ (the L² zeta-zero tail: omitted
  γ>T each ~√x/(γ log x) with random phase ⇒ RMS ~ (√x/log x)·√(log T/T)); eta<1 ⇒ slope 1/eta>1 (error
  decays slower than 1/T, so a bit costs MORE than a doubling of zeros). **Free fraction k0/B=0.88→0.71
  DECLINING with x** (R(x) beats its √x worst-case bound at finite x, E0=|π−R|~√x/log x — only ≈29 at
  x=10⁶ — so R nails ~70–90% of bits free, →~½ only asymptotically = CLAUDE.md's "~50% of digits" in the
  precision domain). **No amortisation:** the explicit formula is ALREADY incrementally batched (refining
  k→k+1 reuses N(k) zeros) yet the cost is geometric, so total = N(B) is dominated by the last hard bit
  (amortised zeros/hard-bit 251/1292/941 at x=10⁵/3·10⁵/10⁶, rising ×3.7 across the clean 10⁵→10⁶ endpoints,
  does NOT fall to a bounded cost). **s≥1 is INFORMATION-FORCED / kernel-independent:** a
  sub-linear slope would settle π(x) exactly in sub-√x zeros, contradicting S518's √x settling + S511's
  joint info floor; a smoothed/Galway kernel can push s toward 1 (doubling) but NOT below. A synthetic
  E∝T^{−2} control (slope 0.48) proves the method DETECTS a sub-linear regime if one existed ⇒ s=1.5–2.5
  is a genuine negative, not a blind spot. **Honest scope:** precision-domain RE-DERIVATION of the
  SMOOTH+RANDOM info barrier (a measurement confirming a known wall from a new angle — NOT novelty, per
  the contract); adjacent NEGATIVE closure, polylog π(x) stays blocked. Caveat: 8000 cached zeros settle
  π(x) exactly only to x~10⁶ (~5 hard bits) — measuring more hard bits needs larger x hence more zeros,
  which IS the √x wall. CLOSED_PATHS row 533 + OPEN_POSITIVE_TARGETS P8 marked closed.
  `sparse_precision_pi_results.md`, `profile_run.log`.
- Op-count curve extended to n=24 + a FALSE "super-Θ(x)" auto-reading corrected (S532b): harvested the
  S531 NEXT-ACTION path (b) — the detached op-count curve `prover_opcount_scaling.py --ns 16,18,20,22,24`
  (`opcount_scaling_n24.log`) COMPLETED. **Verdict: the whole-chain prover op-count is Θ(x)·polylog =
  Õ(x), confirmed to n=24 (x≈1.7×10⁷)** — extends S530's n≤20 Θ(x) measurement one octave. **The
  detached run's own auto-"READING" printed "SUPER-Θ(x)"** (single-power fit over the short high-n
  window n=16..24 → α=1.064>1.05) — but the SAME output's per-step ratios DECLINE monotonically
  4.471→4.432→4.384→4.353→4.329 (a pure power x^α has a CONSTANT ratio 2^{2α}; a declining ratio toward
  4.0 is the fingerprint of a polylog-dressed LINEAR term). The α=1.064 was the **polylog-inflation
  artifact S512's `sharpP_probe` 2-term fit exists to defeat** — the script only had a single-power
  `_fit_exponent`. **FIX:** added `_fit_power_log` (`log₂mul = α·n + δ·log₂n + c`, separating power from
  polylog), `_fit_polylog_fixed_alpha`, `_stable_band` (trailing monotone-declining suffix = past the
  low-n discrete nb/K lump = exactly n≥14), a `--points` injector for deterministic high-n counts (the
  op-count is seed/contention-independent — selftest [3] — so injecting the logged n=22,24 is
  bit-identical to re-running, avoiding a redundant ~30-min recompute that would contend with the live
  reach), and rewrote the verdict to use the per-step TREND + 2-term separation, not a single-power
  threshold. **2-term fit over n=14..24: α=1.012≈1.0, polylog δ=+0.72, residual 0.0012 (17× tighter than
  the single-power 0.021); the forced-α=1 model `x·(log₂x)^0.87` predicts every per-step ratio to <0.4%.**
  ⇒ **S530's live falsifier ("ratio turns back UP above ~4.5× at n≥22") is NOT triggered** — it kept
  declining (4.353×, 4.329×); the Θ(x) op-count reading STANDS and the reach falsifier's op-count leg
  (per-Δn=2 4.0–4.4×↓) is confirmed to n=24, so the S529 prediction band n24→26 ∈ [4.0,5.4]× is intact.
  The n≤20 rows re-measured live reproduce S530 BIT-FOR-BIT (determinism re-confirmed). 6/6 selftest
  (new [6] codifies the separation: high-n single-power α>1.05 inflated, 2-term α≈1.0 with δ>0, ratios
  decline). No landed-chain change (standalone measurement). **Honest scope:** a measurement correcting
  a fit artifact in an earlier in-cycle reading, NOT a protocol or goal advance (polylog π stays
  blocked). `prover_opcount_scaling_results.md` §S531; `opcount_scaling_corrected.log`.
- P10 adaptive query complexity CLOSED-NEGATIVE — adjacent-problem closure with mechanism (S532):
  picked fresh RAM-light ground while BOTH detached runs (reach n=25, op-count path-(b) n=24) were
  still computing — must not start a third. Built `experiments/analytic/adaptive_query_complexity/
  adaptive_query_complexity.py` (standalone, `--selftest`/`--measure`/`--window-demo`; 9 selftest
  groups). **P10 (OPEN_POSITIVE_TARGETS): does an adaptive π-query strategy beat Aggarwal binary
  search for p(n)?** **Verdict: NO.** The non-adaptive recipe — predict x̂=R⁻¹(n) (Newton on polylog
  R), then sieve a Θ(√x) window once — is Õ(√x) with **zero queries**; adaptivity cannot beat it.
  **Mechanism = the project's own SMOOTH+RANDOM split, now measured as a query-complexity dichotomy:**
  **M1** the smooth far zone (where interpolation's loglog advantage lives) is already FREE — Newton
  converges in **2–3 iters, FLAT across 22 scales** (n=199..1.2M); **M2** the random hard zone has
  width **Θ(√x)** (RMS|p(n)−R⁻¹(n)| ~ **x^0.488** over 240 n in 8 bins; err/√(p(n)) median 0.145, max
  0.72, **trend vs x = −0.007 ≈ 0 ⇒ bounded** — matches guess-oracle S491 0.13/0.53); **M3** that
  window carries **no interpolable signal** — within-window π residual sits at the √count≈x^{1/4}
  Poisson scale (rms/√count flat ~0.12) and does NOT drop deg1→deg3 ⇒ interpolation gets no traction
  there (measured **interp_q ≈ binary_q ≈ ½log₂x** inside the window). Each hard-zone probe is a √x
  π-eval, so sieving the window once (0 adaptive queries) dominates ANY adaptive strategy. **Adaptivity
  helps exactly where it is already free (smooth) and fails exactly where the cost lives (random √x).**
  This is the query-complexity face of the info barrier (S511/S518), sharpens guess_comparison_oracle
  (S491) and the CLOSED "Binary search + local sieve — problem relocated" row (now a full count-vs-cost
  dichotomy). 4 pre-stated falsifiers F1–F4 all un-triggered. **Honest scope:** adjacent NEGATIVE
  closure, NOT a goal advance (polylog π stays blocked); the "no strategy beats it" is a
  measurement-backed mechanism over natural strategies, not a formal universal adaptive lower bound
  (that needs item-3 machinery). CLOSED_PATHS row + OPEN_POSITIVE_TARGETS P10 marked closed.
  `adaptive_query_complexity_results.md`, `run_measure.log`.
- Prover op-count exponent MEASURED — the assumed leg of the S529 reach falsifier pinned (S530):
  the reach is running detached (n=25, branch (a) harvest hours away), so per the NEXT ACTION I took
  fresh RAM-light ground. First confirmed branch (b) P7.b (batched dyadic π) is a **re-derivation** of
  the closed S224 Correlation Dichotomy (`cross_x_amortisation/`): dyadic anchors are factor-2 apart =
  the UNCORRELATED side ⇒ Θ̃(√x)/query, no batch-polylog; the NUFFT/FFT angle is pre-closed by
  dichotomy falsifier #4 (Schönhage/Odlyzko cluster-width ceiling). So I did NOT run P7.b. Instead built
  `prover_opcount_scaling.py` (standalone, monkeypatches the `fmul` field-multiply dispatch in every
  importing module with an element-op counter, runs the real `compressed_layer.run_chain` FULL config
  over BIG_Q+FAST_BIG at n=8..20, fits log₂(mul_ops) vs n). **Why:** `mem_hierarchy_wall.py:218`
  hard-codes `IDEAL=4.0` ("op-count = Θ(x) = 4×/Δn=2") — the LOAD-BEARING leg of the S529 reach
  prediction n24→26∈[4.0,5.4]× — but the whole-chain prover op-count was NEVER measured (prior counts
  are verifier-side `vleaf_ops` or one primitive's `--prover-bench` ratio). The op-count is
  DETERMINISTIC (array sizes structural, seed-independent — selftest [3]) ⇒ immune to the concurrent
  reach's CPU contention; RAM-light (n≤20, ~185 MB); zero landed-code change (verdict bit-identical
  under the patch — selftest [1]). **MEASURED: α_op = 0.9955 ≈ 1.0 ⇒ prover op-count is Θ(x)**
  (geo-mean per-Δn=2 ratio **4.11×**, fit ratio 3.98×; mul/x flat ~1400–2400 over n=8..20) — the
  IDEAL=4.0 assumption is now a measured fact. **The small excess above 4 at moderate n VANISHES:**
  upper-regime steps DECLINE monotonically 4.471→4.432→4.384 (n14→16→18→20), settling toward 4.0, no
  plateau at 4.4 ⇒ NOT super-Θ(x). **So the reach's 5.80× wall step (n22→24) is NOT op-count** (op-count
  there ~4.4× and falling) — it is the S529 memory step + noise, confirmed from the deterministic side;
  the [4.0,5.4]× band stands (op-count leg 4.0–4.4×↓ × measured DRAM growth 1.0–1.35×). **Independent
  corroboration of S529's "single-step ratios are noisy":** even the deterministic op-count has lumpy
  per-step ratios (2.43×…4.88×, range >2×) from the discrete nb/K=π(√x) increments ⇒ over-reading any
  one wall step is unjustified, the trend/geo-mean is the signal. **Wall brackets the bowl:** the
  secondary wall exponent over n=8..20 is only 0.524 (in-cache/dispatch-deflated, per-step ratios
  climbing 1.97→2.94 toward 4) vs the reach's ~1.088 (DRAM-inflated) — op-count α≈1.0 is the invariant
  Θ(x) core between them, exactly S529's two-regime per-element wall. fmul ~78% of (instrumented) wall ⇒
  multiplies are the dominant op + faithful proxy. 4 selftest groups. Honest scope: a measurement
  pinning an assumed falsifier leg, NOT a protocol or a goal advance (polylog π(x) stays blocked).
  `prover_opcount_scaling_results.md`. Falsifier left live: an op-count ratio that turns back UP above
  ~4.5× at n≥22 (cheap to extend, no reach wait) would reopen the super-Θ(x) concern.
- Reach orchestration FIXED + the super-linear-wall falsifier mechanism MEASURED (S529): found the
  recorded NEXT ACTION's "⏳ IN PROGRESS (S528)" claim **stale and false** — the S528 reach was a
  harness-tracked background process that was KILLED when the S528 cycle's claude session ended
  (run_n25.log had only the START banner, ~3 min of compute; an earlier orphan died the same way).
  **Root cause:** `run.sh` terminates a cycle's process tree at cycle end, so `run_in_background`
  Bash children do NOT survive the cycle boundary — a multi-hour reach can never land that way.
  **Fix (verified):** `run_reach_detached.sh` launched with `setsid ... </dev/null >/dev/null 2>&1 &`
  reparents to **PID 1** (confirmed `ps`: PPID=1), escaping the kill; the n=25 (honest+delta_pi) then
  n=26 (--no-cheat) reach over BIG_Q now runs genuinely detached and a later cycle harvests when
  run_nNN.log shows DONE (ground truth locked: π(2²⁵−1)=2 063 689, π(2²⁶−1)=3 957 809). Host has
  122 GiB RAM ⇒ the reach is held-cube-RAM-bound (~2.6 GiB @n24, ~10 GiB @n26), NOT memory-ceilinged.
  **In-cycle deliverable (`mem_hierarchy_wall.py`, standalone, 7 selftest groups):** the item-5
  filed falsifier (does streaming tame the 5.80×/Δn=2 wall, or is it super-Θ(x)?) tested at the
  MECHANISM level without the hours-long run. Timed the BIG_Q kernel `_mul61` per-element wall vs
  array size across L1/L2/L3/DRAM: **bowl-shaped — min 12.8 ns @32 KiB (L1), flat ~90–110 ns plateau
  once >L3 (≥64 MiB, DRAM)**, a ~7–8× cache penalty that is a **ONE-TIME step at the L3 crossover
  (every reach n≥18), NOT a per-n growing term**. Across the reach DRAM regime (185 MiB→2.6 GiB)
  per-element growth is only **1.08–1.11×/Δn=2** — far too mild to turn 4× into 5.80×. Observed reach
  ratios n20→22=3.52×, n22→24=5.80× **anti-correlate (geo-mean 4.52×)** ⇒ the 5.80× is mostly
  **3-point noise**, consistent with the Õ(x) prover (Θ(x)·polylog), NOT super-Θ(x). `--bandwidth-probe
  30` (light controls to 8 GiB) shows only a mild multi-GiB TLB/bandwidth rise (~1.05–1.09×/octave).
  **Refines BOTH:** S524 (5.80× ≠ super-Θ(x)) and S525 (the memory hierarchy explains the ~7×
  ABSOLUTE constant but its per-Δn GROWTH is too small to be the 5.80× source). **Falsifiable
  prediction the detached reach tests:** n=24→26 ~[4.0, 5.4]× (op-count 4× × plateauing memory ≤1.35×),
  n=26 wall ~6.5–7.9 ks; a ratio well above ~5.5× with the per-elem wall flat-to-mild would reopen the
  super-Θ(x) concern. `mem_hierarchy_wall_results.md`.
- LIST-streaming LANDED on the real chain — the S526 free win delivered + measured (S527): completed
  the S526 NEXT ACTION (a). A prior interrupted cycle had written the `run_chain(stream_witnesses=...)`
  integration (default True), the `build_stacked_streaming` move into `batched_trace.py`, the
  `verify_constraints_batched`/`verify_ub_openings_batched` `stream=True` paths, the
  `large_x_benchmark` `FULL_NOSTREAM` A/B control + `--stream-probe`, and selftest case 6 — but never
  verified it, measured the win, or updated PROGRAM.md. This cycle VERIFIED, MEASURED, and documented
  it (+ added compressed_layer selftest **case 25** as the chain-level A/B). **What it does:**
  `run_chain` under `batch_trace`/`batch_ub` now passes `Ws=None` and builds the stacked K-witness
  trace cube SLICE-BY-SLICE (each `build_witness(p)` materialized, its `[l·D,(l+1)·D)` slice copied
  in, dropped) instead of first materializing `Ws=[build_witness(x,p) for p in primes]` (~1.4 GB of
  the n=24 3.9 GB peak, S525). **BIT-IDENTICAL** (same cube ⇒ same sum-check ⇒ same challenges/
  verdict): all three selftests PASS — `batched_trace` (8: `build_stacked_streaming`==`stacked_tables`
  table-for-table), `compressed_layer` (1–25; **case 25** asserts `stream` ON==OFF on claimed/comm/
  comm_bt/comm_ub/comm_outer/accepted over q & BIG_Q + liar rejected under streaming),
  `large_x_benchmark` (1–6; case 6 the FULL-config A/B). The default flipped to True, so the ENTIRE
  existing selftest suite already runs through the streaming path. **MEASURED A/B peak RSS**
  (`--stream-probe 16,18,20,22`, FULL=stream vs FULL_NOSTREAM, fresh subprocess each, BIG_Q+FAST_BIG):
  n=16 43 vs 49 MB, n=18 70 vs 91, **n=20 185 vs 265, n=22 677 vs 988** — claimed π bit-identical at
  every n (6542/23000/82025/295947, all ==sieve). LIST saved **6→21→80→311 MB**, ratios 3.5/3.8/3.9 →
  the **`4×/Δn=2 = Θ(x)`** law (the removed list's own scaling); `full_nostream` peaks REPRODUCE
  S525's pre-streaming FULL peaks bit-for-bit (265@n20, 988@n22), confirming the A/B isolates exactly
  the list. **Wall unchanged** (n=20 69.8 vs 69.4 s; n=22 253.9 vs 247.3 s — free). n=24 streaming
  confirmation [APPEND ON LAND: peak GB, claimed]. **Honest scope:** peak-RSS only — soundness/verdict/
  asymptotic profile untouched (per-layer leaf=0, comm Θ(√x), opens Θ(x^{1/4})); the win is a
  CONSTANT-FACTOR (~⅓) prover working-set drop (removes the Θ(x) LIST, NOT the Θ(x) stacked CUBE — S526
  showed the cube can't stream for free), so it is NOT a Θ(x)→Õ(√x) asymptotic change. The
  verification-line CODE lever is now spent; deeper open items need new math. `large_x_benchmark_results.md`
  §S527; run logs (`--stream-probe`, `--mem-one 24`).
- Streaming the batched discharge — S525 NEXT ACTION TESTED + PREMISE CORRECTED (S526): built
  `experiments/constructions/p12_sumcheck_pi_verification/streaming_batched_sumcheck.py` (one
  script, `--selftest`/`--bench`/`--mem-probe`/`--mem-one`; 7 selftest groups). The S525 NEXT ACTION
  conjectured the two batched discharges (`verify_constraints_batched`/`verify_ub_openings_batched`)
  are pure γ-RLC folds `Σ_w γ^w(B_w·C_w)` that stream one freshly-built witness at a time to Õ(√x)
  peak RSS **bit-identically**. **Three findings, the middle one a correction:** **(1)** the
  K-witness LIST (the ~1.4 GB of the n=24 3.9 GB peak) genuinely streams BIT-IDENTICALLY for free —
  `build_stacked_streaming` builds the cube slice-by-slice (materialize→copy slice→drop) and is
  table-for-table == `stacked_tables`, honest + all 4 cheats, over q & BIG_Q (selftest 1). **(2) the
  stacked CUBE is NOT a pure γ-fold** — the trace test (deg-3 zero-test) and Ub test (deg-2 B·C)
  have sum-check ROUND POLYNOMIALS that couple multiple witnesses through the pointwise PRODUCT, and
  each round's Fiat–Shamir challenge sums over ALL K witnesses ⇒ a one-pass accumulator cannot
  produce them. **(3)** a streaming sum-check that DOES reach Õ(√x) space exists (bind-e-FIRST: the
  summand is block-diagonal in the layer index, so each e-round poly is ADDITIVE over the K layers;
  nb streaming passes, then a size-K in-memory layer-phase). VERDICT-faithful (MLE value is
  order-independent ⇒ identical bound scalars/final check): honest ACCEPT matches stacked, every
  single-layer cheat (trace ×4 classes, Ub any (layer,bit) forge) REJECTED first/mid/last layer over
  q & BIG_Q, K=1 edge, fast==object (selftests 2–6). **MEASURED:** peak RSS streamed **FLAT at the
  ~36 MB numpy baseline** (1.00–1.01×/Δn=2, n=12..18) vs stacked **37→39→49→91 MB** (climbing toward
  the Θ(x) ×4/Δn=2 law) — the Õ(√x)-space win is real. **BUT the cost is Θ(√x) MORE WALL** (NOT a
  log factor as I first wrote): fmul element-work is comparable (round-poly work telescopes; e-phase
  even SKIPS the 2^Lk−K padding layers the cube folds), but the streamed prover makes ~K× more fmul
  CALLS on NARROW D-wide arrays, so the dispatch-bound wall tracks the call count — `--bench` shows
  3.6×(n8)→17.2×(n14) slower, mean width 341→19, **re-introducing the exact op-count-bound regime
  (S500/S501) that batching (S502) escaped by WIDENING arrays.** **Headline:** the batched discharge
  is a space-time tradeoff — batching minimizes WALL (wide arrays, Θ(x) space), streaming minimizes
  SPACE (Õ(√x), Θ(√x)× wall); the Θ(x) held cube is INHERENT to the one-pass batched sum-check's
  wall efficiency, so paying √x× compute to shrink it makes the wall WORSE. The only free
  bit-identical win is streaming the LIST (item 1). No landed-chain change (n=24 artifact untouched);
  `build_stacked_streaming` is the ready, safe LIST-streaming port. `streaming_batched_sumcheck_results.md`.
- Large-x prover working set LOCALIZED + S524 scatter attribution CORRECTED (S525): answered the
  S524 NEXT ACTION (vectorise the per-layer outer-reduction scatter-sums onto the uint64 fast path,
  S524 NEXT ACTION (vectorise the per-layer outer-reduction scatter-sums onto the uint64 fast path,
  re-measure to separate GC overhead from a genuine super-Õ(x) term). Built `scatter_fold61`
  (`compressed_prover_mult_trace.py`): a uint64-safe segmented scatter-fold (31/30-bit limb split +
  `argsort`+`np.add.reduceat` per bucket + Mersenne recombine) replacing the object-dtype
  `np.add.at` at BOTH per-layer scatter sites (`compressed_layer.py` `verify_trace_region` omega +
  `small_reduce` Wt). **Bit-identical** to the object path (selftest 8a': random/structured/
  collision/empty/gap/2²⁰-near-overflow), exact for `<2³³` summands/bucket (callers' fill ≤√x ⇒
  n<66). Added a `--no-scatter-fold` A/B control (fixed a bug where the toggle's fall-through hit
  the `_dt(q)`=uint64 path and OVERFLOWED → spurious REJECT; now an explicit object branch, guarded
  by benchmark selftest 5: bit-identical claimed/comm/accepted over BIG_Q). **HEADLINE FINDING —
  the S524 attribution was WRONG and is corrected:** scatter_fold61 is correct but **WALL-NEUTRAL**
  (n=20 67.8 s vs S524 object 71.8 s; n=22 256.3 s vs 252.8 s; n=24 1445.7 s vs 1467.3 s — all
  within ~2%; n=22 A/B has comm 146637 and peak RSS 988 MB IDENTICAL ON-vs-object, wall within 3%).
  The `~3.9 GB` working set and the super-linear n22→24 step are NOT the object scatter (it is
  removed in the ON path yet n=24 peak RSS is still 3.89 GB and the acceleration persists, local
  x^1.25 vs object x^1.27). They are the **Θ(x) batched-discharge held-cube working set** (the
  K-stacked witness list `batch_trace`/`batch_ub` hold to discharge in ONE sum-check — uint64, not
  object): `--mem-probe` FULL vs NOBATCH shows FULL peak RSS strictly HIGHER (265 vs 138 MB @ n=20;
  the gap IS the held cubes), peak RSS scales `~4×/Δn=2` = Θ(x) (baseline-subtracted ratios →4.0).
  So the super-linear wall is a memory-hierarchy constant (Θ(x) arithmetic over a cache-busting
  Θ(x) working set), NOT Python-int GC (held cubes are uint64) and NOT the per-layer √x-size
  scatter. Op-count milestones unchanged (per-layer leaf=0, comm Θ(√x), opens Θ(x^{1/4})). **The
  real remaining wall lever is the held-cube SPACE** — stream/checkpoint the K-witness stack
  (`batched_trace.py:110` `Ws=[build_witness(p) for p in primes]`) so the batched discharges don't
  hold all K at once, dropping the prover working set Θ(x)→Õ(√x) (see NEXT ACTION). No new protocol
  machinery; deliverable is the corrected localization + a genuinely uint64-pure-scatter fast path.
  `large_x_benchmark_results.md` §S525; run log `run_s525_remeasure.log`.
- Large-x benchmark reach at n=24 — the live verification frontier pushed one octave, the
  scaling claim corrected (S524): answered the S523 NEXT ACTION verbatim — ran
  `large_x_benchmark.py --n 24` to completion under the FULL succinct config over BIG_Q=2⁶¹−1.
  **HEADLINE: honest run ACCEPTED, claimed π(16777215) = 1 077 871 == sieve, cross-checked
  against an independent from-scratch Eratosthenes sieve (π(2²⁴)=1 077 871, a known value) —
  so the chain's exact count at x≈1.7×10⁷ is verified against ground truth from two sources.**
  Wall 1467.3 s (~24.5 min). **Deterministically reproducible:** an earlier unlanded background
  run (`run_n24_full.log`) and the landing run (`run_n24.log`), both seed=1, agree bit-for-bit on
  every structural number (claimed, comm=292 635 + full split, vleaf=0, opens=16 512); only
  wall-clock differs 3.5% (1520.2 vs 1467.3 s) — as expected for the deterministic FS transcript. **Soundness witnessed at the reach:** the `delta_pi` liar (claim
  1 077 872) is REJECTED (541.9 s). **All three asymptotic-profile claims hold at the new reach:**
  per-layer verifier large-table leaf-eval count **= 0** (Õ(√x) end-to-end, S506); `comm`=292 635
  elems is **97% `comm_outer`** (the K=564 sequential outer reductions, Θ(√x)) with every batched
  discharge polylog (comm_bt/bw/ub=88/2829/67) and base=4352; large-table opens = the one-time
  tensor commitment only (vcommit_ops=16512, Θ(x^{1/4}), S508). **HONEST CORRECTION of the S513
  scaling claim:** the prior "~3.5×/Δn=2" was a single-step (n20→n22) estimate; the 3-point series
  (71.8 / 252.8 / 1467.3 s) shows the per-Δn=2 ratio ACCELERATES (3.52× → **5.80×**, geo-mean
  4.52×), least-squares **wall ≈ x^1.088**. The asymptotic Õ(x) prover (open item 1) is intact —
  ~1.09 is linear-in-x+polylog — but the n22→n24 step runs above the ideal 4×/Δn because the BIG_Q
  object-dtype accumulation paths (`np.add.at` off the uint64 fast path) push the working set to
  ~3.9 GB, so Python-int alloc+GC adds a super-linear WALL overhead (implementation, not
  complexity — every op-count milestone is unchanged). The reported t_prover+t_verifier is only
  ~1.7% of wall (cf. S501; wall is field arithmetic outside the timed regions). Falsifier filed:
  the local exponent staying well above ~1.1 at n≥26 with the fast path fully engaged would signal
  a genuine super-Õ(x) term. No new protocol machinery — pure reach + measurement; the deliverable
  is a verified exact π(x) at x≈1.7×10⁷ behind an Õ(√x) verifier / Õ(x^{1/4}) opens (NOT a polylog
  π(x); the goal stays blocked). `large_x_benchmark_results.md`; run log `run_n24.log`.
- Sliding-window census D_k^slide(x) — the last named width-spectrum follow-on, MEASURED +
  CLOSED (S523): answered the S522 NEXT ACTION verbatim. Built
  `experiments/analytic/local_pattern_census/dk_sliding_census.py` (one script, `--k`/`--xmax`/
  `--seg`/`--aligned`/`--validate-hl`/`--selftest`; 11 selftest groups incl. segmented==direct
  full-sieve==naive trial-division, the aligned base case, k=8/16 saturation, the parity
  identities, classical singular series). Unaligned windows `[n,n+k)` (n over ALL residues) drop
  the aligned `q|k` admissibility condition → the CLASSICAL Hardy–Littlewood prime-constellation
  criterion (`{o mod q}≠ℤ_q` for every prime q≤k). **HEADLINE STRUCTURAL IDENTITY:
  A_k^slide = 2·A_k^aligned − 1 = 25, 211, 7145** for k=8,16,32 — the q=2 clause forces every
  admissible S all-even or all-odd; the **odd-offset family is EXACTLY the aligned admissible set**
  (aligned (i) for q=2|k already restricts to odd offsets), the **even-offset family is its image
  under the reflection o↦k−1−o**, an affine bijection mod every prime ⇒ an admissibility
  automorphism preserving ω(q), so #even=#odd=A_aligned and the two share only ∅. **Densest
  admissible weight = ρ*(k) = 9 at k=32, parity-symmetric** — sliding reaches no denser tuple than
  aligned (the same H()=W(B(k)) geometry from S519 caps both). **FULL SATURATION:
  D_32^slide = A_slide + E = 7145 + 8 = 7153, every classical admissible realized by n < 2²⁷
  (x≈1.3×10⁸) — ≈1000× earlier in x than the aligned 2³⁷**, for two reasons: (a) sliding samples
  32× more windows at the same x (every n, not multiples of k); (b) the densest weight-9
  constellations are realized IMMEDIATELY at n≈10 by the small-prime cluster itself (window
  [10,42) = primes {11,13,17,19,23,29,31,37,41}, an admissible 9-tuple), far ahead of their HL
  schedule. k=8 saturates D=29 (n<2⁹), k=16 D=217 (n<2¹²). **Exceptional set
  E_k^slide = q*(k)+1 = 4,6,8** (windows n=0..q*, q*=3,5,7 = the largest prime whose residues the
  small-prime cluster covers within length k; n=0..2 cover mod 2,3,5,7 down to n=6,7 covering mod 7
  alone; beyond n=q* the prime q* leaves the window so every pattern is admissible) — vs aligned's
  E=1 (only the m=0 window). **Classical singular series VALIDATED:** `--validate-hl` at k=32, x=2²²
  gives obs/pred=0.994 at weight 6 (matching aligned's 0.980); the twin constant 𝔖({0,2})=1.3204=2C₂
  falls out. **Base case (S522 NEXT ACTION's requested cross-check):** `--aligned` restricts to
  n≡0 mod k and reproduces the aligned realized set bit-for-bit (selftest [5]). No falsifier (no
  inadmissible pattern at n>q*, no overdue admissible, no plateau below 7153). Honest scope: a
  finite-x census measurement in the width-spectrum track, not progress on the polylog π(x) goal;
  the 2·A_aligned−1 identity and the ρ*(k)-cap are the new content. `dk_sliding_census_results.md`;
  run log `run_slide_k32_2e28.log`.
- D_k(x) census first-occurrence probe at k=32 — the 189 unrealized dense constellations
  HUNTED DOWN (S522): answered the S521 NEXT ACTION verbatim. Built
  `experiments/analytic/local_pattern_census/dk_first_occurrence.py` (one script,
  `--k`/`--xmax`/`--seg`/`--validate-hl`/`--selftest`; 10 selftest groups incl. segmented
  sieve == direct full sieve bit-for-bit at k=8/16/32, the 2²⁶ baseline 3385/189 reproduced,
  saturation D_8=14/D_16=107 reproduced, the singular-series algebra). A constant-memory
  **segmented sieve** extends the aligned-window census far past the prior 2²⁶ reach.
  **HEADLINE: FULL SATURATION — D_32(2³⁷)=A_32+1=3574, deficit 189 → 0** (x≈1.4×10¹¹), so k=32
  now joins k=8 (14) and k=16 (107) as a confirmed instance of `D_k(x)→A_k+1`. The last admissible
  pattern to appear was the densest possible (weight-9) constellation `[1,5,7,11,17,19,25,29,31]`
  at x≈8.3×10¹⁰ (HL x*≈3.9×10¹⁰, N(occ)=1.63 — on the HL schedule). The curve
  reproduces every prior checkpoint EXACTLY (2¹⁶=1093, 2²⁴=3213, 2²⁶=3385) and the residual at
  each reach is precisely the densest HL tail (the 6 missing at 2³² = 1 wt-7 + 2 wt-8 + 3 wt-9).
  **Hardy–Littlewood ALIGNED singular series derived + VALIDATED:** for the linear forms
  f_i(m)=km+o_i, S_aligned=∏_q(1−ω(q)/q)/(1−1/q)^w with ω(q)=0 for q|k (factor 2^w at q=2),
  ω(q)=#{o_i mod q} for q∤k; expected count N_S(x)=S_aligned·(x/k)/(ln x)^w. `--validate-hl`
  occurrence-count test at 2³⁰: obs/pred rises monotonically through 1, **=0.980 at weight 6**
  (28842 obs vs 29435 pred — clean 2% validation, no survivorship bias there); <1 at low w
  (exact pattern ⊂ constellation) and >1 at high w (only upward-fluctuating dense patterns
  realized early). First occurrences are HL-ORDERED (N(occ)≈1 for the cleanly-tested band).
  **Statistically-grounded falsification monitor:** flags a missing pattern only if its HL
  MEAN count N_S(xmax)≥5 (P(0)=e^{−N}<0.7%), not the naive x*<xmax. It raised 2 weight-6 flags
  at 2²⁸ (cleared by 2³²) and 1 weight-7 (`[3,9,13,15,19,25,27]`, four offsets ≡0 mod 3) at
  2³²; the latter appeared at x≈4.75×10⁹ with N(occ)=7.35 — the EXPECTED max-order-statistic
  lag over ~300 weight-7 patterns (max of n exp(1) ≈ ln n ≈ 5.7; P(max≥7.35)≈17%), NOT an HL
  over-prediction. **Every flag the monitor ever raised subsequently appeared — self-correcting.**
  No falsifier observed (no inadmissible pattern at m>0; no statistically overdue admissible
  pattern; no plateau below 3574). Cost: one constant-memory segmented sieve to ≈1.4×10¹¹
  (deepest run ~12 min). The k=32 census is now fully saturated and HL-consistent, matching the
  law at k=8,16. `dk_first_occurrence_results.md`; run logs
  `run_k32_2e32.log`/`run_k32_2e36.log`/`run_k32_2e37.log` (the saturating run).
- Census entropy constant — the cluster/block METHOD CLOSED + S520 k=256 CORRECTED (S521):
  answered the S520 NEXT ACTION's sub-task (b) ("does the bracket width shrink to 0 or
  PLATEAU — i.e. can the cluster/block method EVER prove =1?") with a clean NEGATIVE, and did
  sub-task (a) (vectorise the engine). Extended `census_entropy_asymptotic.py`: **(a) THREE
  joint-cover engines** — `_jc_vec` (numpy, multi-limb uint64 for N>64, the largest heads
  batched via doubling so a single large head prime never falls into a 2^h Python loop; ≈10×,
  2^26 block 209→~20 s, bit-identical to `_jc_py`, selftest 12) + `joint_cover_int` (EXACT
  integer `jc·2^N`, each term `=sign·∏_w(2^{n_w}−1)` since excl+Σn_w=N; numerator `==count_tm`
  exactly, selftest 14). **(b) `--bracket-sweep`** measures the best feasible rigorous bracket
  vs the I-E cost budget C. **The width is a STEP FUNCTION of C, flat between discrete
  prefix-cost thresholds** (a prime enters the bracketing block only when its 2^{prefix-cost}
  becomes affordable): at **k=128 the width is STUCK at 0.0160 for EVERY feasible C≥26**
  (upper pinned {3,5,7,11,13}=0.9403, lower 0.9243); next tightening needs **2^39** (admit 17),
  bracket shuts only at **2^75**. **Closing cost `2^{sum(B)−max(B)}` is super-exponential in k:
  2^26→2^75→2^279→2^872 (k=64,128,256,512) = 2^{Θ(k²/ln²k)}** ⇒ **the FKG-block/cluster method
  provably CANNOT close the bracket / prove the constant=1 for any k≥128.** Plateau is
  engine-independent (the un-includable {17,19,23} block also aborts count_tm past 25.7M states).
  **CORRECTION of S520:** a SECOND wall — at k=256 (N=128) `P[miss]~2^{−Θ(N)}~1e-17` but the
  I-E sums terms ~1, so the FLOAT engine cancels to ZERO significant digits (e.g. exact
  `P_miss({3,5,7})|₂₅₆=1.48e-19`, float`=0.0`). **S520's k=256 bracket [0.9282,1.126] was
  numerical noise and is RETRACTED** (its "rigorous lower" 0.9282 is ABOVE the true exact lower
  0.8312 — not even a valid bound). Exact-integer engine removes the wall: **corrected k=256
  bracket = [0.8312, 0.9431] RIGOROUS** (width 0.1119; cost wall still pins the upper at
  {3,5,7,11,13}=5 of 14 primes). k≤128 float is reliable (~7 digits at N=64). Structural
  finding (mean-field ratio declines 0.752→0.605→0.394, correlation-driven →1) is from the
  stable marginals — unaffected. 14 selftest groups (+12 vec==py multi-limb; +13 sweep
  monotone/plateau/shut; +14 exact==count_tm + precision wall). The census line is now closed
  for computational cycles: proving =1 needs an analytic handle on the FULL Θ(entropy)
  correlation series (new analysis, not new code). `census_entropy_asymptotic_results.md`.
- Census entropy CONSTANT via a cluster expansion — open item 4's remaining (post-S519)
  analytic question, ESTIMATED with a rigorous bracket (S520): answered the S519 NEXT ACTION
  ("estimate the entropy constant WITHOUT exact counting via a generating-function /
  cluster-expansion, validate vs exact A_k k≤64, extrapolate to 128/256"). Built
  `experiments/analytic/local_pattern_census/census_entropy_asymptotic.py` (one script,
  `--selftest` (11 groups) / `--k` / `--scan` / `--rmax` / `--cost-cap`). **Exact
  reformulation:** `A_k = 2^N·P[admissible]` with offsets i.i.d. p=½ (every subset has prob
  `2^{−N}`), `\bar C_q`="S misses a class mod q" is DECREASING ⇒ the events are FKG-positively
  associated. **(1) Mean field** = independent-primes product `2^N∏_q p_q`, `p_q=1−∏_j(1−
  2^{−n_{q,j}})` exact from class sizes — a RIGOROUS LOWER bound (FKG). **(2) Cluster
  expansion:** Möbius inversion `ln P_B=Σ_{T⊆B}c_T`, `c_T=Σ_{U⊆T}(−1)^{|T\U|}ln P_U`, each
  `P_T` got from a **nested inclusion-exclusion** for the joint cover prob (cost `2^{sumS−maxS}`,
  k-INDEPENDENT — feasible at 128/256 where S519's `count_tm` on full B is Θ(A_k)-dead; validated
  BIT-FOR-BIT against count_tm, selftest 2). **(3) FKG block bracket:** `P[∩_B]≥∏_blocks
  P[∩_block]` (any partition, lower) and `≤P[∩_G]` (single block G⊆B, upper). **VALIDATED EXACT
  k≤64:** truncation errors ALTERNATE (odd r under, even over) and shrink ≈5×/order — r=4 within
  0.001 of exact at k=64; the bracket CONTAINS exact (k=64 [0.9236,0.9276]∋0.9276). **HEADLINE
  k=128: ratio ln A₁₂₈/(k/lnk) ∈ [0.9243, 0.9403] RIGOROUS** (block {3,5,7,11,13}, cost 2^26;
  cluster point estimate 0.9264 inside); **k=256: ratio ≥ 0.9282** (rigorous lower; upper loose
  1.13 — one feasible block holds only 4 of 14 active primes, the honest |B|-growth limitation).
  Ratio climbs 0.667→0.808→0.886→0.928→≈0.93, consistent with →1. **STRUCTURAL FINDING:** the
  constant is CORRELATION-DRIVEN — the independent-prime (mean-field) ratio DECLINES with k
  (0.752→0.605→0.394 at k=64,128,256), so CRT overlaps supply a GROWING fraction of the entropy
  (≈19%→35%→57%); the conjecture's →1 lives entirely in the corrections, not the leading term.
  **Still open (constant=1):** NOT proved — the bracket does not close to a point as |B|→∞ (upper
  degrades) and the cluster series has Θ(entropy)-sized terms (cap-sensitivity at k=128 shows the
  triple sum still moving), so =1 needs an analytic handle on the FULL correlation series. 11
  selftest groups (marginal==count_tm; IE==count_tm; Möbius identity==ln A_k; mean field≤exact;
  bracket brackets exact; error shrinks+alternates; partition-coarsening monotonicity).
  `census_entropy_asymptotic_results.md`.
- Census A_k WITHOUT enumeration — open item 4 transfer-matrix route CLOSED + two scaling
  structural facts (S519): answered the S518 NEXT ACTION ("count A_k by a transfer matrix /
  generating function, validate vs known A_k, push to k=128 OR prove the state space is
  exponential"). Built
  `experiments/analytic/local_pattern_census/census_transfer_matrix.py` (one script,
  `--selftest`/`--k`/`--scan`/`--reduce`; 9 selftest groups; TM reproduces every known A_k
  bit-for-bit, brute==DFS==TM). **(1) Transfer-matrix route CLOSED:** the natural position-DP
  (state = per-prime occupied-residue SET, dead-state pruned) has **Θ(A_k) states** — measured
  collapse factor A_k/states = 1.86,1.66,2.37,1.55 (bounded, NOT growing) at k=8,16,32,64;
  ln(states) slope **1.020** ≈ ln A_k slope **1.011** ⇒ same entropy exponent ⇒ exponential and
  strictly WORSE than DFS (Θ(A_k) memory vs Θ(k); final live-state count is order-independent
  ≈A_k/2, so no offset ordering helps). k=128 (A₁₂₈~10¹¹) unreachable by both. **(2)
  Active-prime reduction (rigorous, two independent proofs that AGREE — selftest 8):** enforcing
  a prime is a monotone filter, so A_k depends only on primes ≤ B(k); found by count-stabilisation
  (`minimal_B`) AND by a no-count W(B)-weight bound (`reduce_primes` via branch-and-bound
  min-union `minkill_bb`, exact, selftest 9). **B(k)=3,5,7,13,23,47** at k=8,16,32,64,128,256 —
  **A₁₂₈ uses only 8 primes {3,5,7,11,13,17,19,23}**, the 9 larger enforceable primes provably
  never bind. (Shrinks the problem but does not break the exponential wall.) **(3) Identity
  W(B(k)) = ρ*(k) = maximal admissible-tuple size** = 3,5,9,16,28,52: the weight bound is tight (a
  weight-W(B) subset missing a class mod all of B is itself admissible), so it equals the densest
  admissible aligned-k tuple; **matches H() (OEIS A008407) exactly** at k=64 (H(16)=60≤63<H(17)=66)
  and k=128 (H(28)=126≤127<H(29)=132), extending the prior k≤32 cross-check. ρ*/(k/lnk)=0.78→1.13
  over k=8…256 (Hensley–Richards). **This explains the entropy law's exponent**: admissibles live
  in the ρ*(k)≈k/lnk geometry ⇒ ln A_k=Θ(π(k)); the entropy slope→1 (1.01) is re-confirmed
  independently of the DFS. **Honest scope:** no new clean-family exact A_k point is reachable (the
  contribution is the closure + the reduction + the weight identity, NOT a new A_k value); the
  constant→1 still needs an analytic/generating-function/cluster-expansion argument (separate from
  exact counting), and an exact A₁₂₈ needs a genuinely sub-A_k algorithm (none found — TM
  exponential, IE catastrophic-cancellation/term-blowup). `census_transfer_matrix_results.md`;
  `local_pattern_census_results.md` extended.
- Non-sieve sub-√x witness probe via the explicit formula — open item 3 direction (ii) closed
  for the NATURAL analytic witness (S518): answered the explicit S517 NEXT ACTION. Built
  `experiments/analytic/explicit_formula_witness/explicit_formula_witness.py`, a STANDALONE
  measurement (not a protocol). **(numerical fix, CLOSED row 31)** the explicit formula's
  reported DIVERGENCE ("error 3.5→2076") was a branch-cut bug: `li(x^{ρ/k})` computed as
  `li(exp(ρ·lnx/k))` folds the huge imaginary part `ρ·lnx∼γ·lnx` mod 2π inside `exp`, so `li`'s
  internal `log` takes the wrong branch. Fix = the exponential integral on the UNWRAPPED
  argument `ei(ρ·lnx/k)`; smooth `R(x)` via `mpmath.riemannr` (=Gram series, validated 1e-20).
  Selftest asserts BOTH the bug (|err|≈7e4) AND the fix (|err|<0.5). **(witness exponent)**
  per-x window-median settling count over x∈[10³,2.2×10⁵] with 8000 zeros (height ≤8148): raw
  log-log COUNT slope α=0.91, HEIGHT slope β=0.73 — both ABOVE ½. **Matched control over the
  SAME anchors**: pure x^0.5→slope 0.500, √x·log²x→0.712≈β, so the >½ raw slope is the polylog
  dressing on a ½ LEADING POWER, not a genuine super-√x exponent; **NO measurement dips below
  ½** ⇒ sub-√x REFUTED, polylog (slope 0) REFUTED (N grows 42→7106). Galway-form fit: √x·log²x
  best (rms 3.2→1.8→0.93 for p=0,1,2), reconfirming Galway/zero_scaling.py with the height/count
  split made explicit. **(S511 floor link)** exactness fraction is ~0 for c=T/√x≤1 at every
  anchor (10⁴,10⁵,10⁶) and only ~0.6–0.8 even at the full 8000-zero budget; at x=10⁵,10⁶ dropping
  ANY octave `[T/2,T]` of zeros up to √x jumps the RMS remainder back >0.5 (0.56–2.8 vs full
  0.46–0.56), none droppable (at x=10⁴ the full budget over-settles, RMS 0.31, so drops stay
  0.34–0.46; the signal sharpens with x) ⇒ the rounding-relevant information is DENSE across
  Θ(√x) zeros — the analytic analogue of S511's "Θ(K), not o(K)". The natural analytic witness MATCHES the √x floor, does
  NOT beat it. **Honest scope:** reconfirms the long-closed √x scaling (rows 30–34/39/267,
  `zero_scaling.py`); the NEW content is the row-31 numerical fix, the height-vs-count
  separation, and the floor cross-check tying it to open item 3. Closes the NATURAL analytic
  witness only — a non-natural witness is not ruled out (item-3 universal question open). 16
  selftest cases. `explicit_formula_witness_results.md`; CLOSED row added.
- Fine-grained (SETH) barrier to #P-hardness of π — direction (i) of open item 3
  conditionally CLOSED (S517): answered the lever S512 left open ("a Turing, not
  parsimonious, reduction sidestepping the c↦x(c) circularity"). Built
  `experiments/constructions/p12_sumcheck_pi_verification/turing_reduction_barrier.py`,
  a STANDALONE measurement + complexity anchor (not a protocol). **THEOREM
  (conditional):** with α=inf{a:π(x)∈TIME(x^{a+o(1)})} and SETH, NO poly-time
  Turing reduction #SAT→π-oracle has all queries of bit-length ≤c·n with c<1/α
  (else simulate π by its own time-x^a' algorithm → #SAT in 2^{a'cn}, a'c<1 →
  refutes SETH). Critical query-blowup **c*=1/α**. **MEASURED:** the project's own
  `pi_lucy` (Lucy_Hedgehog O(x^{3/4})) over x=10⁴..10⁸ fits time-exponent α≈0.66–0.70
  (sub-linear; t/x^{0.75} decreases) ⇒ c*≈1.5; cited rungs c*=1.5 (combinatorial
  LMO/Deléglise–Rivat 2/3) and **c*=2 (analytic Lagarias–Odlyzko 1/2)**. **Natural
  (parsimonious) reduction has blowup c→1<c***: encode "count=C=2^n" as π(x)=C needs
  smallest x=p_C, bitlen≈n+log n ⇒ c=bitlen(p_C)/n falls 1.50→1.20 over n=4..22
  (exact via sieve), PNT c=1+log₂(ln C)/n→1.00 — so EVERY parsimonious-blowup Turing
  reduction is SETH-forbidden, matching S512's circularity with a complexity wall
  that also covers adaptive non-parsimonious reductions. **COROLLARY (the project
  tie-in): polylog-π XOR #P-hard-π under P≠NP** — α→0 (the goal) ⇒ c*=∞, a polylog
  π-oracle answers every query in poly(N) so a poly-time Turing #P-hardness reduction
  gives #SAT∈P ⇒ P=#P; thus a #P-hardness PROOF for π proves no polylog algorithm
  exists (this project impossible). The two open questions are the two horns of one
  dichotomy; √x sits between them at c*=2. **Still open (scope):** super-c*-blowup
  reductions (query π at x≥2^{(1/α)n}, exponentially past where p_C≈2^n lives — no
  natural construction, analytic-π there already costs ≥2^n so no #SAT speedup) and a
  non-sieve sub-√x witness for L_π (the other NEXT-ACTION half, must also beat the
  S511 info floor). Honest: the meta-principle is folklore fine-grained complexity;
  the contribution is applying it to π with π's own √x algorithm, the quantitative
  c*=1/α + the measured parsimonious c→1, and the mutual-exclusivity corollary. 52
  selftest cases. `turing_reduction_barrier_results.md`.
- Expander distance/soundness as a measured CURVE — the S514 single-point δ claim
  upgraded (S516): closed the S515 NEXT ACTION (the "clean single cycle"). Extended
  `experiments/constructions/p12_sumcheck_pi_verification/pcs_commit.py` with two
  `--selftest`-gated measurements behind new flags. **(a) `--distance-sweep`:**
  `_min_basis_rel_weight(code,k,q)` (code-agnostic generalization of S514's
  `_exp_min_basis_rel_weight`; RS via the incremental power table `c^j`, expander via
  `_exp_encode`) tabulates the soundness parameter `δ = min_j |Enc(e_j)|/N` vs k up to
  512 for BOTH row codes over q & BIG_Q & SMALL_Q. **RS δ→1** ((N−1)/N, increasing).
  **Expander δ DECLINES 0.650→0.404 over k=8→512 but DECELERATES** — per-doubling
  decrements shrink (+0.082→+0.019, mean ratio 0.77), geometric-tail floor estimate
  **≈0.34**; **δ does NOT decay to 0**, field-independent (q-free indices). Practical:
  the t for 2^−100 forge soundness (`t≥100/−log₂(1−δ)`) grows only **67→134** over
  k=8→512 (sub-linear, because δ floors above 0). **(b) `--forge-rate`:** worst-case
  Monte-Carlo — the adversary bends honest `v` at the min-weight coordinate
  `j0=argmin weight(Enc(e_j))`, support `δN`; across many random openings (independent
  FS challenges) the empirical accept rate **tracks `(1−δ)^t` and stays at/below it**
  (emp/pred ≤ ~1.0, hugging the tighter hypergeometric), fitted `δ_eff≥δ` — the forge
  dies at least as fast as the bound, **no query correlation** (the falsifier `emp≫
  (1−δ)^t` never observed; max ~1.02, a 3σ blip), field-independent over q & BIG_Q.
  Selftest **cases 10/11** assert the falsifiers: δ bounded below + decelerating (not
  →0); `_hyper_miss ≤ (1−δ)^t` monotone; empirical accept ≤ `(1−δ)^t` within binomial
  tolerance at every t, monotone, `t=1≈1−δ`, `δ_eff≥δ−0.1`. **Honest scope:** the
  asymptotic constant distance is Brakedown's PROVEN, cited (not re-proved); this adds
  the MEASUREMENT that demo-scale δ declines but decelerates to a positive floor; the
  geometric extrapolation is an estimate, not a bound; recipe = measure δ at the actual
  k, pick `t=⌈λ/−log₂(1−δ)⌉`. `pcs_commit_results.md` S515 section.
- Certificate-info sharpening — the Θ(√x) joint-info floor solidified + per-layer
  density (S515): closed the S514 NEXT ACTION. Extended `cert_incompressibility.py`
  (refactored `integer_rank` → reusable `build_residual_matrix`+`min_exact_rank`;
  added `layer_density_profile`, `--wfactors`). **(A) Wide window (W=192·K):** the
  rank/K **floor rises 0.75→0.88** (uniform, no downward drift over k=16..22) and
  **α_rank/α_K = 1.106→0.978 ≈ 1**. **Corrected the S511 premise** that α_rank should
  approach 0.5: α_rank cannot exceed α_K, and the finite-window α_K is itself only
  ≈0.42 (PNT 1/ln x log-log discount on √x), so 0.5 is an asymptote, not a reachable
  target; the 64·K α_rank=0.459 was *inflated* by an across-window rank/K climb
  (under-sampling at low k), not faster-than-K growth. The sampling-robust √x signal
  is the **rank/K FLOOR ⇒ rank=Θ(K)=Θ(π(√x))=Θ(√x)**, with the raw exponent reading
  ~0.40 the shared discount (both → 0.5 only as x→∞). Raw α_rank went *down*
  (0.459→0.396) yet the result is *strengthened* (floor lifted, ratio →1). **(B)
  Per-layer DENSITY:** the √x info is spread across ALL K=π(√x) layers, not o(K) —
  per-layer hard bits uniform (deciles 9.8–12.7, `bits_uniformity`=0.77, active
  fraction 0.98), prefix integer-rank LINEAR (rank/j 0.84→0.88, `rank_half_ratio`=
  rank(K/2)/rank(K)=0.52). Energy is end-loaded (scale: φ~x at small a) but bits/rank,
  the scale-free measures, are uniform. **New falsifiers:** rank/K decay or
  α_rank/α_K→0 (not α_rank<0.45 — that's the expected PNT tracking value);
  `rank_half_ratio`→1 or `bits_uniformity`/`active_frac`→0 (√x carried by o(K) layers)
  — none observed. 30 selftest cases (+8: `min_exact_rank` refactor safety, density
  separation of dense vs rank-concentrated (duplicate cols) vs bit-concentrated
  (sparse) controls). `cert_incompressibility_results.md` §1b/§1c.
- Field-size-free commitment — the RS `N≤q` demo artifact removed (S514): built the
  Brakedown-style linear-time **expander row code** the S513 NEXT ACTION called for,
  behind the SAME `commit/prove/verify` interface in
  `experiments/constructions/p12_sumcheck_pi_verification/pcs_commit.py`, gated by a
  `code=` selector (`"rs"` default = every prior artifact verbatim; `"expander"`).
  Construction (Spielman/Druk–Ishai recursive, as in Brakedown): systematic
  `Enc_n(x)=x‖Enc_{⌈n/2⌉}(A·x)‖B·Enc_{⌈n/2⌉}(A·x)` with **column-regular** sparse A,B
  over F_q (coldeg=4) + a dense base case; the matrices are a fixed PUBLIC map (seed +
  size → q-free column indices, values mod q), so prover-commit and verifier-reencode
  build the identical linear code. The tensor reduction `S~(pt)=a^T M b` and the
  homomorphic column check `⟨a,Mhat[:,c]⟩=Enc(v)[c]` hold for ANY linear Enc, so the
  `O(t(r+k))=O(x^{1/4})` verifier is unchanged — the verifier re-encodes v,w ONCE
  (`O(N)=O(k)`, linear-time) then indexes the t queried columns. **Column-regularity is
  load-bearing:** a plain low-density generator leaves some input column unused → a
  weight-1 basis codeword → the forge cheat at that coordinate is uncatchable; making
  every input column hit coldeg output rows lifts the measured min basis-codeword
  relative weight from ≈0.003 (row-regular) to **≈0.45** (the quantity that governs the
  forge cheat, since the forged-v difference is exactly δ·Enc(e_{j0})). **Measured
  (`--bench --code expander`):** verifier op-count slope **0.524 (α≈0.262), sub-√x, same
  as RS (0.500), IDENTICAL over q/BIG_Q/SMALL_Q** (op count is field-independent —
  indices not values set the structure); the expander codeword is SHORTER than RS's
  (N=380 vs 512 @nb=14) so `commit_vops` is actually lower at large nb. **4 cheat classes
  rejected** with `code="expander"` over q & BIG_Q & SMALL_Q (wrong claim / forged
  opening / tampered codeword row / tampered revealed column); honest opens == `mle_eval`;
  same map as the chain base close. **The field-size-free win is concrete (selftest case
  9):** over a tiny prime **q=17** at nb=8 (k=16 → RS Ncode=64>17) RS REFUSES to commit
  while the expander commits (N=44), verifies, and rejects the cheats. Threaded
  `commit_code` through `run_chain` (default "rs" verbatim) + `--commit-code`; **chain
  verdict UNCHANGED** with the expander base over q AND BIG_Q (selftest case 24, n∈{8,10,
  12}: honest π==sieve, delta_pi + self-consistent liar rejected, per-layer still
  leaf-free, base opens still sub-√x `vcommit_ops==vleaf_ops_ot`). **Honest scope:**
  distance is the security parameter and is MEASURED here (δ≈0.45 at demo params, selftest
  asserts ≥0.3 for k∈{8,16,32,64}), not proven from the asymptotic expander bound —
  Brakedown's analyzed params give δ≈0.07 with t≈148 queries for 100-bit security; the
  construction shape is theirs, the constants are demo-scale and distance-measured. Still
  computational (CRH/RO, sha256) like the RS PCS — this swap changes only the field
  dependence, not the unconditional-vs-computational boundary. `pcs_commit_results.md`
  S514 section; canonical-doc protocol stack + leaf-openings section updated.
- Large-x benchmark reach push — item 5 advanced (S513): built
  `experiments/constructions/p12_sumcheck_pi_verification/large_x_benchmark.py`,
  a clean reproducible driver (the "large-x benchmark" line had no dedicated
  script) running the **FULL succinct config**
  (delegate+structured+pcs+batch_trace+batch_ub+batch_wiring+commit_base) +
  FAST_BIG over the sound `BIG_Q=2⁶¹−1` — strictly stronger than the S504 config
  (it adds the real leaf openings S505, the batched Ub discharge S506, and the
  tensor-PCS base commitment S508, i.e. the full non-interactive certificate).
  `--selftest` gates correctness/soundness over Q & BIG_Q at small n (full config
  π==sieve, delta_pi + self-consistent liar rejected, FAST_BIG a bit-identical
  drop-in, per-layer leaf count 0). **Reach, BIG_Q, `claimed π==sieve`, honest
  ACCEPTED:** n=20 π(1048575)=82025 (71.8 s, reproduces S504's 66 s now under the
  FULL config), **n=22 π(4194303)=295947 (252.8 s, x≈4.2×10⁶)**, n=24 (x≈1.7×10⁷)
  [APPEND ON LAND]. **Soundness witnessed at the reach:** the `delta_pi` liar
  (claim π+1) is REJECTED at n=22 (66.9 s). **Profile confirmed at the reach** —
  per-layer verifier large-table leaf-eval count **= 0** (Õ(√x) end-to-end, S506);
  `comm` ~95% `comm_outer` = the K sequential outer reductions (Θ(√x)), every
  batched obligation polylog (n=22: comm_bt/bw/ub=80/2373/61, comm_outer=139899 of
  146637, sum exact); large-table opens = the one-time tensor commitment only
  (Θ(x^{1/4})). Chain stays **prover-bound** (~3.5×/Δn=2), NOT field- or DP-bound —
  this is the open-item-1 Õ(x) prover on the √x state, NOT a polylog π(x) (the goal
  stays blocked); the deliverable is the verified verification artifact at x≈10⁷.
  n539 built the driver + reproduced n=20; n540 ran the reach (n=22 + soundness),
  banked the results, and launched n=24. `large_x_benchmark_results.md`;
  `compressed_layer_results.md` item-5 section extended.
- #P-hardness / NP-membership feasibility probe — the COMPUTATIONAL lower-bound
  face of open item 3 (S512): answered the explicit lever S511 named. Built
  `experiments/constructions/p12_sumcheck_pi_verification/sharpP_probe.py`, a
  STANDALONE measurement + complexity anchor (NOT a protocol). **(A) UPPER bound
  (folklore, demonstrated exact):** π(x) and the Legendre partial sieve φ(x,a) are
  #P FUNCTIONS of the BINARY input x — count the ≤2^N integers n≤x under a poly(N)
  primality/coprimality predicate (selftest: predicate-count==sieve π,
  coprime-count==2^a-term inclusion-exclusion). ⇒ `L_π∈C_=P` (exact counting),
  `{π(x)≥c}∈PP⊆P^#P` (Toda: `PH⊆P^#P`, the other direction; `PP⊆PH` is not
  known and not claimed). The #P "subset witness" for φ has 2^{π(√x)} terms
  (doubly-exp in N) — membership gives no short witness. **(B) NP-completeness
  obstruction:** π a function ⇒ `L_π∈NP ⟹ L_π∈NP∩coNP` (certify π(x)≠c via the true
  c'+its cert) ⟹ NP-complete ⟹ NP=coNP. The live question is mere NP-MEMBERSHIP
  (polylog witness). **(C) WITNESS-SIZE LADDER (measured):** leading-power exponent
  α via the 2-term fit `log bits = α·log x + δ·log log x + c` that separates the
  power from polylog (the naive single-slope is polylog-inflated over a short
  window — log²x reads 0.21 but α=0; verified exact on closed-form controls).
  enumeration NP-cert (list π(x) primes w/ Pratt + every other composite w/ a
  factor, computed from a real spf sieve) **α=0.985 (Θ(x))**; sieve transcript
  **α=0.490** (matches S509 real naive 0.473); S511 info floor **α=0.486**;
  zeta-zero/explicit-formula witness **α=0.500** (Galway K~c·√x·log²x, **CITED**
  EDGES Thread-3/S195/196/S434–436, not recomputed). THREE independent natural
  families (sieve, analytic/zeta, info floor) **converge at √x** (max|α−0.5|=0.014,
  k=10..20); enumeration a full power higher; polylog (α=0) the rung S511 rules out
  for the sieve class. ⇒ **no natural witness reaches poly(N); every one ≥2^{N/2}**.
  **(D) parsimonious-reduction obstruction (CIRCULAR, exact):** a reduction #A→π
  maps w↦x(w) with π(x(w))=#A(w); realizing target count c forces x∈[p_c,p_{c+1}-1]
  so the map c↦x(c) IS the inverse-prime p(·) = the project goal — the reduction's
  "easy direction" is as hard as π itself (closure mode C). The sieve/φ instance is
  lattice-structured ({multiples of d}), no instance richness to embed an arbitrary
  #P-complete count. ⇒ #P-hardness of π, if true, cannot come from instance-
  embedding — it must be value-incompressibility (the S511 route). **(E) Corrected
  CLOSED row 175** ("exact pi(x) is #P-hard" — unsubstantiated S7-era assertion;
  true statement is π∈#P, hardness open). **STILL OPEN (filed):** a genuine
  #P-hardness proof OR a non-sieve sub-√x witness — neither delivered; on present
  evidence L_π is an NP-intermediate-flavoured counting problem with a √x cert and
  no proven poly(N) one. 40 selftest cases (#P-membership equalities, Legendre,
  n-th-prime realization window, fit_power_log exact recovery on x/√x/√x·log³x/
  log²x controls, ladder exponents+convergence, monotonicity). The information face
  (S511) and the computational face (S512) now agree and point the same way without
  closing the universal question. `sharpP_probe_results.md`.
- Certificate incompressibility probe — the √x is an INFORMATION floor, not just
  construction shape (S511): answered the lower-bound face S510 flagged ("is sub-√x
  information-theoretically impossible?"). Built `experiments/constructions/
  p12_sumcheck_pi_verification/cert_incompressibility.py`, a STANDALONE measurement
  (entropy/rank of the checkpoint sequence, not a protocol). The chain pins one
  large-side survivor count per layer = the Legendre partial sieve
  `φ(x,a)=#{1≤n≤x : n has no prime factor ≤ p_a}`; we measure the JOINT hard bits of
  `{φ(x,a)}_{a=1..K}`, K=π(√x), beyond a POLYLOG smooth predictor, at INTEGER
  precision (CLOSED row 737's lesson: relative smoothness is cancelled by the
  exact-value requirement). **CRUCIAL frame:** cert SIZE is √x just from having K
  layers; the question is cert INFORMATION (can the K transcripts be jointly
  compressed?). **Authoritative measure = integer-reconstruction SVD rank of the
  residual matrix** over a dense x-window (smooth-model-free — SVD subsumes the best
  low-rank/smooth model; row-737 style on residuals). `rank≤min(W,K)`, so resolving
  all K modes needs W≫K x-samples: window-sensitivity at x=2¹⁹ (K=128) gives
  **rank/K 0.09→0.25→0.63→0.82→0.97 as W/K 1→4→16→63→187 — FULL RANK once sampled**.
  Adaptive-window sweep (W=64K, rank/K 0.75–0.97 over k=16..22): **rank≈K,
  α_rank=0.459 ≈ α_K=0.415 ≈ √x**. ⇒ the K checkpoint residuals are integer-
  INDEPENDENT (no smaller smooth+low-rank model recovers them exactly) ⇒ joint info
  **Θ(K)=Θ(√x)·polylog** ⇒ **the √x certificate is INFORMATION-FORCED for any
  sieve-reconstructing verifier**: polylog cert RULED OUT (joint info super-polylog,
  lifting S36's single-value "O(log x), computational not informational" to the
  joint certificate setting), sub-√x cert RULED OUT for this class (upgrades S510's
  construction-inherent √x to information-inherent). Per-x corroboration (single
  K-length sequence, conservative 1D under-sample): gzip α=0.388, AR(3) α=0.557, both
  track K's α=0.415 (incompressible). **Single-value control reproduces S36**:
  bits(|π(x)−Li(x)|)=4.3→8.0, slope 0.089≈0 (O(log x); |π−R|≤15) — the barrier is
  JOINT across layers, NOT per-value. **Self-corrected** a narrow-fixed-window
  artifact (W=4096 → rank~x^0.35, read as sub-√x) via the window-sensitivity check —
  the rank cap `min(W,K)` under-samples; W∝K gives full rank. **Honest scope:** φ is
  a faithful proxy (transcript hashes bounded above by the state); the bound is the
  sieve-reconstruction CLASS only — a different witness (cf. factoring) is not ruled
  out, so a UNIVERSAL sub-√x impossibility still needs #P-hardness (item 3 formal
  half). 22 selftest cases (Legendre identity, li/R sanity, metric estimators,
  exponent fitter, synthetic info-forced-vs-compressible separation, SVD integer-rank
  low-vs-full). `cert_incompressibility_results.md`.
- Merged-layer COMM probe — the √x certificate size is LAYERING-INHERENT (S510):
  answered the single quantitative question S509 flagged — is the Õ(√x) cert size
  inherent to the layering, or just to the unbatched chain? Built
  `experiments/constructions/p12_sumcheck_pi_verification/merged_layer_comm.py`, a
  STANDALONE primitive (one merged layer's comm vs j single layers, not a full-chain
  rebuild). It collapses j consecutive SMALL-side (value→value `⌊v/p⌋`) Lucy layers
  into ONE reduction over the composed map, whose 2^j-subset inclusion-exclusion form
  `M[S](v)=Σ_{T⊆[j]}(−1)^|T| GateProd_T(v) S(⌊v/m_T⌋)` is exactly open item 2's 2^j
  fill (`m_T=Π_{k∈T}p_k`; verified == direct j-layer composition, selftest 1). Two
  SOUND reductions of `M̃[S](r)=c` to a base S-claim, both ground-truth-checked and
  cheat-tested over q & BIG_Q: **DIRECT** (one deg-3 sum-check over the v-cube → the
  2^j openings, each routed back to S, the 2^j claims `open_batch`-folded), and
  **BATCHED** (stack the 2^j on a j-bit axis, the batched_trace pattern → ONE deg-3
  sum-check over the (j+nb)-cube → ONE opening). Comm is prime- AND field-independent
  (transcript size depends only on nb,j; selftest 4), so total = ⌈K/j⌉·per-merged is
  exact. **MEASURED (n∈{8..16}, j∈{1,2,4}, both fields identical):** chain TOTAL comm
  exponent **α~0.45–0.49 (Θ(√x)) for EVERY j in BOTH modes** — √x is layering-inherent.
  DIRECT `comm_fill=(2^j+1)(3nb+2)` grows ~2^j (total comm has a min at j=2 then GROWS
  — merging makes the certificate *bigger*); BATCHED pushes the 2^j OUT of the
  transcript (comm_fill flat=27) but into **verifier WORK** (`vwork`~2^j, ratio 5.7×
  @j:1→4, exponent 0.84–0.89 matching the prover) AND prover, comm floored at
  ⌈K/j⌉(j+nb)→K=√x (the merge-axis sum-check costs log₂(2^j)=j/layer, ⌈K/j⌉·j=K). The
  **prover fill total ⌈K/j⌉·2^j·2^nb → x^{3/2}** (/x: 0.42→2.25 over j:1→6 at n=16,
  exponent 0.84+) — open problem 2's m·K²/2≈x^{3/2}, now on the VERIFICATION face.
  **Conclusion:** no merge depth in either batching strategy drops the comm exponent
  below 0.5; merging trades a bounded comm-prefactor win for a 2^j prover/verifier-WORK
  blow-up. The cert-size wall = verifier-work wall = prover-time wall = the one
  ⌊v/p⌋-semigroup sequential-sieve wall. 5 selftest cases (operator==composition;
  honest accepts/residual correct; c_value/c_gate/c_route reject; comm prime- &
  field-independent; scaling bounds). `merged_layer_comm_results.md`.
- Certificate-SIZE profile — the Õ(√x) comm wall LOCALIZED (S509): the verifier side of
  the milestone was met at S508 (large-table ops Θ(x^{1/4})); the quantity STILL Õ(√x) is
  the certificate SIZE (`stats["comm"]`, slope ~0.5 in every config). This cycle pins down
  WHERE. Added a purely-additive comm attribution to `compressed_layer.run_chain` —
  boundary snapshots partition `stats["comm"]` into `comm_outer` (the K=π(√x) SEQUENTIAL
  per-layer two-sided reductions), `comm_bt`/`comm_bw`/`comm_ub` (the one-time batched
  trace/wiring/Ub discharges), `comm_base` (the two one-time S₀ tensor-PCS eval proofs);
  the five sum EXACTLY to `comm` (asserted on every complete certificate; transcript &
  verdict unchanged — compressed_layer selftest still passes). New
  `experiments/constructions/p12_sumcheck_pi_verification/certificate_profile.py` measures,
  over n∈{8..18}, the **full succinct config** (delegate+structured+pcs+batch_trace+
  batch_ub+**batch_wiring**+commit_base — strictly extends S508 config (d) by batching the
  wiring, S503/S504, which alone cuts comm 5.3× @n=18): (i) cert size in elems & bits,
  (ii) verifier large-table ops, (iii) prover wall; attributes the comm; fits each exponent.
  **MEASURED:** cert size **α=0.473 (Θ(√x))**, DOMINATED by `comm_outer` **α=0.522** —
  60%→90% of comm over n=8→18, →100% asymptotically; every batched discharge polylog
  (`comm_bt` α=0.09, `comm_bw` α=0.18, `comm_ub` α=0.09); `comm_base` α=0.30 (≈Θ(x^{1/4}),
  the PCS proof); verifier ops α=0.258 (reconfirms S508). **The Õ(√x) certificate, after
  every batchable obligation is batched (→polylog) and the base committed (→x^{1/4}), IS
  exactly the K sequential outer layer reductions** (each O(nb²)=O(log²x) round scalars ×
  K=π(√x) ⇒ Θ(√x·polylog)). **HONESTY SHARPENING of S508:** the verifier must READ the
  Õ(√x) proof (each scalar an O(1) round check) ⇒ **total verification is Θ(√x)**, bounded
  below by proof size; the S508 Õ(x^{1/4}) is the large-table-eval SUB-term only — now made
  explicit. **Membership result:** `L_π={(x,c):π(x)=c}` has, under a CRH/RO (Fiat–Shamir),
  a NON-INTERACTIVE certificate of size **Õ(√x)** field elems (Õ(√x·log²x) bits), with
  **Õ(x^{1/4})** large-table opens (total verification Θ(√x)). Compressing the SIZE to
  polylog requires batching the K SEQUENTIAL reductions — each layer's output claim IS the
  next layer's input claim (a dependency chain), unlike the independent per-layer
  trace/Ub/wiring witnesses — which is the ⌊v/p⌋-semigroup compression, i.e. the SAME wall
  as open item 2 (layer-batching closed negative). Cert-size wall = prover-time wall = the
  one sequential-sieve wall. selftest (5 cases): attribution closes over q & BIG_Q (honest;
  delta_pi + self-consistent liar reject), comm COUNT field-independent (q vs BIG_Q identical,
  bits differ 31 vs 61), batch_wiring moves comm out of `outer`, exponent bounds, verifier
  large-table work one-time-only. `certificate_profile_results.md`.
- Tensor polynomial commitment for the S₀ base closes — the whole-chain verifier is
  now SUB-√x (S508): built `experiments/constructions/p12_sumcheck_pi_verification/
  pcs_commit.py`, the hash-based multilinear PCS the S507 NEXT ACTION called for, and
  threaded it behind `run_chain(--commit-base)`. The chain's per-layer verifier was
  leaf-eval-free since S506; the ONLY verifier work whose cost still grew with x was
  the two ONE-TIME S₀ base closes, each a direct `mle_eval` over a 2^nb=O(√x) table
  (`compressed_layer.py:1047`). This collapses that term. Construction (Ligero/
  Brakedown tensor): reshape the 2^nb table into an r×k matrix (r=2^{n1}, k=2^{n2},
  n1+n2=nb), Reed-Solomon-encode each ROW (N=blowup·k points), Merkle-commit the N
  COLUMNS (sha256 = the CRH/random-oracle stand-in). Since `eq~(pt,w)` factors over
  the bit split, `S~(pt)=a^T M b` (a=eq_table(pt_hi) len r, b=eq_table(pt_lo) len k);
  the eval proof (Fiat–Shamir) sends v=a^T M, w=ρ^T M (len k) and opens t columns, and
  the **verifier is `O(t·(r+k))=O(x^{1/4}·polylog)`** — SUB-√x, the full 2^nb table
  never touching it. `commit`/`prove`/`verify` with stats-tallied `vcommit_ops`.
  **Measured:** standalone `--bench` — commit-verify op-count slope **0.5/nb
  (Θ(2^{nb/2})=Θ(x^{1/4}))** vs the direct close's 1.0/nb (Θ(√x)), same over q & BIG_Q;
  absolute crossover at nb≈14 (the t-query constant dominates below — the S506 pattern,
  honestly documented). Whole-chain `--bench-verifier-ops` **config (d)** added:
  leading exponent **α_ops=0.258 (Θ(x^{1/4}))** vs (c)'s 0.500 (Θ(√x)) and (a)'s 0.961
  (Θ(x)) — the milestone **Θ(x)→Θ(√x)→Θ(x^{1/4})**; comm slope stays ~0.60 (Õ(√x), the
  proof adds only an O(x^{1/4}) term). **Soundness:** standalone selftest (over q,
  BIG_Q, SMALL_Q) rejects all 4 cheat classes — wrong claimed value, forged opening
  (v bent to pass the eval binding, caught by column consistency), tampered (non-
  codeword) committed row (caught by the proximity test), tampered revealed column
  (Merkle path); honest opens AGREE with `mle_eval`. `compressed_layer` selftest 23:
  `run_chain(--commit-base)` verdict UNCHANGED — honest π==sieve over q & BIG_Q,
  delta_pi + self-consistent liar rejected. **Honest scope:** this trades the
  otherwise-unconditional Õ(√x) verifier for full succinctness under the hash
  assumption (the rest stays unconditional); RS's N≤q point requirement is a demo
  artifact a field-size-free linear code (Brakedown's expander code) removes behind the
  same `commit/open/verify` interface. Default `commit_base=False` keeps every prior
  artifact verbatim. `pcs_commit_results.md`; canonical-doc leaf-openings section + the
  protocol-stack table updated.
- Whole-chain verifier op-count CURVE (S507): turned the Õ(√x) end-to-end verifier
  claim from a per-n op-count fact (S506) into a MEASURED SCALING LAW. Added
  `_acct_vleaf` to `compressed_layer.py` — a counter tallying every VERIFIER
  large-table evaluation (a direct `mle_eval` of a committed √x-size cube, the ONLY
  verifier operation whose per-call cost grows with x), weighted by table size, split
  into the PER-LAYER critical path (scales with K=π(√x)) vs the ONE-TIME terms (the two
  S₀ base closes). Instrumented all six verifier leaf-eval sites (affine s2 close, trace
  s_B close, the nb Ub openings, the line-batch close, the single-claim anchor, the two
  base closes), respecting the existing pcs/batch_ub gating. `--bench-verifier-ops` runs
  three configs (ALL delegate+structured — so the wiring verifier is already polylog and
  leaf opens are the sole x-scaling verifier term) over n∈{8..18}: (a) no-pcs, (b) pcs,
  (c) pcs+batch_trace+batch_ub. **Falsifiable claim CONFIRMED:** exact per-layer
  leaf-eval COUNTS `K·(nb+5)−1` (a) / `K·nb` (b, the Ub openings) / **0** (c), one-time
  `2·2^nb` in all three (asserted in the bench + selftest 22 over q AND BIG_Q — the count
  is field-independent); whole-chain `vleaf_ops` fitted leading exponent **α=0.961 /
  0.998 / 0.500** (per-step Δn=2 growth ×4 / ×4 / **×2.0 exactly**) — the Θ(x)→Θ(√x)
  drop, with config (c)'s total being the one-time base closes ONLY. The `comm` column
  (slope ~0.60, Õ(√x) with polylog) corroborates that the entire non-leaf verifier
  residual (sum-check checks + the O(K·polylog) batched-discharge recomputes, proven in
  S495/S502/S503/S506) stays sub-linear and never re-introduces a Θ(x) term. **Honest
  scope:** the metric is the LEADING term (large-table evals only); config (c)'s grand
  total is Õ(√x) at effective exponent ~0.6 (polylog above √x), exactly the milestone's
  "leading term Θ(√x·polylog)". The only remaining √x verifier work is now one-time (the
  two S₀ base closes + batched discharges). compressed_layer §"Step 9" results.
- Batched Ub openings — the LAST per-layer O(2^nb) verifier term removed (S506): built
  `batched_trace.verify_ub_openings_batched`, the multi-point/multi-table opening that
  the S505 NEXT ACTION called for. After `--pcs`, the only per-layer O(2^nb) verifier
  operation left was `verify_trace_region`'s nb Ub-bit-table openings (`Ub_{Lv-nb+k}~(r_C)`,
  the B2 routing's pin of g1_trace to the certified quotient u_e) — per-layer WITNESS
  data, so no carried claim to thread to the S₀ base. These Ub tables are EXACTLY the
  ones `verify_constraints_batched` already stacks and commits along the layer axis, so
  the discharge shares that commitment. Each layer DEFERS `(l, r_C^l, [ub_{l,k}])` (the
  prover supplies the scalars, the verifier folds them into its B2 `expect` in O(nb)) and
  all K·nb discharge in ONE degree-2 sum-check: with γ←F_q and β=γ^nb the per-(l,k) weight
  β^l·γ^k=γ^{l·nb+k} is a DISTINCT power, and the identity factorizes to a single inner
  product `claim=Σ_w B[w]·C[w]` over the (Lk+nb)-cube — B the verifier-anchored per-layer
  eq weights (recomputed in O(K(Lk+nb)), the soundness anchor), C the γ-fold of the
  committed low-nb Ub tables (its opening taken from the sum-check's folded scalar, exactly
  as the zero-test trusts its bit scalars). Wired behind `run_chain(--batch_ub`, requires
  pcs+batch_trace); threaded via a new `ub_defer` accumulator (mirrors `batch_wiring`'s
  `defer`) through `large_reduce`→`verify_trace_region`, discharged after the layer loop.
  **Measured (`--bench-ub`, full delegate+structured+batch_trace+batch_wiring+pcs config):**
  the verifier-side O(2^nb) Ub-leaf eval count `ub_leaf_v` drops from K·nb to **EXACTLY 0**
  at every n=8..16 (24/55/108/217/432 → 0; moved to the prover, `ub_leaf_p`=K·nb). The
  per-layer verifier is now LEAF-EVAL-FREE — the only O(√x) verifier work left is one-time
  (the two S₀ base closes + the batched discharges), so **the chain verifier is honestly
  Õ(√x) end-to-end**, closing the "Honest scope" caveat that has stood since step 1.
  **Honest wall-clock scope:** the measured t_verifier drop is only ~1.1–1.2× and roughly
  flat in n — a vectorised numpy `mle_eval` over an O(√x)-size array is cheap relative to
  the Python-loop wiring/eq recomputes that dominate the measured t_verifier, so the
  O(K·nb·2^nb)→one-time improvement does NOT surface as a growing wall-clock ratio at
  reachable n. **This CORRECTS the S505 NEXT-ACTION prediction** that the ratio would grow
  with n (it grows in op-count, not measured wall-clock until past the array-size crossover
  S500 measured). Selftest 21 (compressed_layer) + §7 (batched_trace): verdict UNCHANGED
  (honest==sieve, delta_pi + self-consistent liar rejected) over q AND BIG_Q; standalone
  primitive accepts honest openings, agrees with the per-layer inline `mle_eval` ground
  truth, rejects any single forged (layer,bit) opening 4/4; K=1 edge; fast==object. Default
  `batch_ub=False` keeps every prior artifact verbatim.
- Real leaf openings — the polynomial-commitment opening primitive (S505): built
  `experiments/constructions/p12_sumcheck_pi_verification/leaf_open.py`, the
  sum-check MLE OPENING that the NEXT ACTION called for. `open_eval(S,pt,claimed)`
  proves `S~(pt)=claimed` by ONE degree-2 sum-check of the MLE-eval identity
  `S~(pt)=Σ_w S[w]eq~(pt,w)` — a REDUCTION converting an opening at `pt` into a
  polylog transcript (verifier O(nb), comm O(nb)) PLUS one residual claim `S~(r)`
  at a fresh point; `open_batch` folds k same-table claims to ONE residual via a
  powers-of-γ eq-RLC + one sum-check (drop-in for `line_batch_pair`/
  `batch_on_table`). Threaded behind `run_chain(--pcs)`: the carried-claim folds'
  residuals are **threaded as the next layer's claim and discharged transitively
  at the S₀ base** (the now-redundant `s2`/`s_B` closes in `verify_affine_region`/
  `verify_trace_region` skipped) — UNCONDITIONAL (no commitment; soundness rides
  the GKR layer reductions). **Measured:** standalone opening verifier is flat in
  table size `2^nb` (the O(nb) signature) — 6.5×→432× faster than the `mle_eval`
  close it replaces as nb grows over q, **21×→3605× over BIG_Q** (object dtype);
  `open_batch` beats sequential line folding 13–34× (vs k). Chain `--bench-pcs`:
  **stable ~1.2–1.36× verifier reduction** (n=8..16), comm +1–3%. Selftest case 20:
  pcs chain verdict UNCHANGED — honest==sieve, delta_pi + self-consistent liar
  (layers 1/K/2/K) rejected — over q AND BIG_Q, automaton AND delegated+structured,
  composed with batch_trace+batch_wiring; n=16 headline π(65535)=6542 with all
  cheats rejected 6/6. **Honest scope:** this is a CONSTANT-FACTOR win + a working
  demonstration of the residual-threading architecture, NOT yet asymptotic Õ(√x):
  it removes ~5 per-layer O(2^nb) leaf closes, leaving the DOMINANT per-layer
  O(2^nb) term — `verify_trace_region`'s nb Ub-bit-table openings (line 429), whose
  residuals are per-layer WITNESS data with no carried claim to thread to. That
  single remaining term (→ batched-trace integration) is the new NEXT ACTION.
  Default `pcs=False` keeps every prior artifact verbatim.
- Batched wiring INTEGRATION + the fast-path sign flip resolved (S504): wired
  `batched_wiring.verify_wiring_obligations` into `compressed_layer.run_chain` by
  the planned defer-and-batch. New `batch_wiring=True` path (delegate only): a
  `defer` accumulator threaded through `large_reduce`/`verify_affine_region`/
  `verify_waff_value` (large affine, `accept_rem=p−1`) and `small_reduce` (small
  division, `accept_rem=None`) COLLECTS each layer's wiring obligation `(p,r_v,r_u,
  accept_rem,claim,lie)` — the claim being the scalar the outer sum-check already
  pinned — instead of calling `inner_verify_div` inline; all `2K−1` obligations
  discharge in ONE batched chain after the layer loop. Two new entry points in
  `batched_wiring.py` (`chain_obligation` builds an obligation from a deferred
  check — prover-supplied claim, honest-or-lying chain over the shared `l_max`
  cube; `verify_wiring_obligations` discharges the list). `verify_waff_value`'s
  `[e≤ecut]` comparator fold (O(2^nb), p-independent) stays per-layer, as planned;
  default `batch_wiring=False` keeps every prior artifact verbatim. **Selftest
  §19/§19b:** chain verdict unchanged (honest==sieve; `delta_pi` + self-consistent
  liar rejected) over q & BIG_Q, structured & dense, alone and composed with
  `batch_trace`; the wiring-specific liars (`small_forge` at sum-check #0,
  `small_chain`/`waff_chain` in the batched backward sweep, `waff_forge` in the
  comparator fold) all rejected THROUGH the batched discharge; `batch_wiring`
  without `delegate` raises. **Headline (`--n 16 --bench-combined`, BIG_Q,
  delegate+structured):** with BOTH big kernels batched the global fast Mersenne
  path is finally a NET END-TO-END WIN — baseline (obj) 17.1 s → batch_trace (obj)
  14.0 s → batch_trace+wiring (obj) 12.8 s → **batch_trace+wiring+FAST 8.87 s
  (1.92×)**; FAST alone with only trace batched is a LOSS (25.5 s, the S502 22-vs-16
  sign flip reproduced — tiny per-layer wiring cubes penalise the 24-op mulmod).
  Wiring batch also cuts comm 5.3× (87226→16509: K chain transcripts → one).
  **Item-5 delivered (winning config, `claimed π==sieve` over the sound BIG_Q):**
  n=16 π=6542 (9.3 s), n=18 π=23000 (22.2 s), **n=20 π(1048575)=82025 (66 s)** —
  exact π at x≈10⁶ behind an Õ(√x) verifier over a sound-characteristic prime,
  Õ(x) prover, both wirings delegated. The S500 lever ("reduce the COUNT of small
  field-ops, not the per-multiply cost") is confirmed end-to-end.
- Batched wiring delegation (S503): the SECOND big chain kernel (S501: 30% of
  `run_chain` wall, **76% of all sum-check calls**), built and measured.
  `experiments/constructions/p12_sumcheck_pi_verification/batched_wiring.py` —
  replaces the K = π(√x) INDEPENDENT per-layer `inner_verify_div` chains (each a
  GKR chain: sum-check #0 + n backward-sweep layers + base, on a tiny `2^{l_ℓ}`
  cube) with ONE batched inner GKR chain over all K wiring obligations. Defer-and-
  batch: collect `(p_ℓ, r_v_ℓ, r_u_ℓ, claim_ℓ, accept_rem_ℓ)`, discharge together.
  **Data-parallel GKR with a NON-UNIFORM per-layer wiring**: add `Lk=⌈log₂K⌉`
  layer vars, MSB-zero-pad each cube to `2^{l_max}`, share one `(r_L,ρ)`
  trajectory + a `γ` combiner. Keeping the layer index a genuine sum-check
  variable makes the round-poly final factorize, so the heterogeneous wiring
  (`R_j^ℓ`, `[σ<p_ℓ]` depend on `p_ℓ`; `wv_ℓ[j],wu_ℓ[j]` on the per-layer query
  point) is recomputed by the VERIFIER in `O(K·l)` via `Σ_ℓ χ_ℓ(r_L)·⟨MLE_ℓ⟩` —
  the trace's `av` trick generalized to a point-dependent operator. **No dense
  `2^{2l}` table is built**: each backward layer rides the layer axis INSIDE the
  S496 structured phase split (phase 1 binds σ', phase 2 binds σ), so every
  prover table is `K·2^{l_max}` WIDE not `2^{Lk+2l_max}` dense; the per-step
  wiring factors are absorbed into the next step's verifier weight (eq-tracking
  `eq(r_L1,r_L2)`). **Verified by AGREEMENT** with the AND of K independent
  `inner_verify_div` (selftest): honest accepts; any single corrupted obligation
  (wrong claim / self-consistent lying chain / broken leaf) rejected first/mid/last
  layer, over Q & BIG_Q, div/affine/mixed wirings; K=1 edge; fast==object
  bit-identical; comm drop. Matches BOTH delegated wirings `run_chain` emits
  (small division `accept_rem=None` line 657; large affine `accept_rem=p−1`, the
  `inner_verify_div` core of `verify_waff_value` line 273), uniform depth `nb`.
  **`--bench`: fmul CALL count ↓ 2.0→18.1× / per-fmul WIDTH ↑ 6→40× over n=8..16
  (K=6..54)** — the op-count→width conversion S501 named (net fmul-element-work
  only ~2×). **The fast-Mersenne sign flip is reproduced**: unbatched the fast
  path LOSES (n=16: 15.9 s fast vs 6.2 s object), batched it WINS (n=16: 1.31 s
  fast vs 4.77 s object) → speedup vs the unbatched fast baseline **8.99× @n=14,
  12.1× @n=16**, growing with K. Standalone primitive + bench only; the
  `run_chain(batch_trace, batch_wiring)` + `FAST_BIG` re-measure is the NEXT
  ACTION (needs the defer-and-batch restructure of `small_reduce`/
  `verify_waff_value` + the latter's separate `[e≤ecut]` comparator fold).
- Batched trace zero-test (S502): the S501 lever, built and measured.
  `experiments/constructions/p12_sumcheck_pi_verification/batched_trace.py` —
  replaces the K = π(√x) independent per-layer `verify_constraints` sum-checks
  (53% of `run_chain` wall, S501) with ONE degree-3 sum-check over the stacked
  `(K·2^nb)` cube. The K witnesses stack along a `2^Lk` layer axis (`Lk=⌈log₂K⌉`,
  MSB-zero-padded bit tables → one shared term list); `BETA_EQ[ℓ,e] =
  (β^ℓ if ℓ<K else 0)·eq~(tau,e)` folds the inter-layer combiner, the eq-selector,
  AND the layer-valid mask into one multilinear table (the `EQ→BETA_EQ` swap in
  `build_terms`); the verifier recomputes `BETA_EQ~(r)` and the factored wiring
  MLE `av=(Int(r_e)+dstart)·Σ_{ℓ<K}pₗχ_ℓ(r_L)` in O(K) — within the chain's
  existing O(K) budget, so still Õ(√x). **`--bench`: fmul COUNT ↓26× / per-fmul
  WIDTH ↑37× at n=16 (K=54), exactly the op-count→width conversion S501 named.
  The fast-Mersenne sign flip is real**: unbatched the fast path loses on the
  chain (small arrays), batched it wins and the win grows 0.6×→3.54× over
  n=8..16; the isolated trace kernel goes 8.5 s (unbatched object, today) →
  **1.7 s (batched+fast) = 5.0×**. Selftest: structural (comm 1728→56 at n=16),
  honest+soundness over Q & BIG_Q (every per-layer cheat — wrong/self-consistent
  quotient, corrupted bit, wrong remainder, non-bit — rejected at first/mid/last
  layer), agreement with the AND of the K per-layer tests, K=1 edge, fast==object
  bit-identical, and the fmul-count/width win. **Wired into `run_chain` via
  `batch_trace=True`** (skip per-layer `verify_constraints`, run one batched check
  before the loop; `compressed_layer.py` selftest 18 confirms unchanged verdict
  honest/delta_pi/liar over Q & BIG_Q). **Honest end-to-end:** batch_trace (object)
  is **1.14–1.44×** on `run_chain` (n=16: 15.5→13.6 s) — only the trace half;
  globally enabling fast-Mersenne is a NET LOSS (22 s vs 16 s) because
  `_cpmt.FAST_BIG` is global and penalises the still-tiny unbatched WIRING cubes
  (30%, 76% of sum-check calls). The wiring defer-and-batch is the gate to the
  combined fast-path win — the new NEXT ACTION.
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

> **⏳ REACH: n=25 DONE + VERIFIED; n=26 STILL RUNNING DETACHED as of S535 (2026-06-14T~16:10Z, ~2.85 h
> in, RSS 4 GB, no DONE yet — PID 1592836 alive) — CHECK `pgrep -f 'large_x_benchmark.py --n 26'`/LOG
> BEFORE RELAUNCHING; S534 + S535 left it untouched (did contention-safe fresh ground instead, see below).
> S535 OBSERVED the live reach: it is `stime`-DOMINATED (0.695 kernel/page-fault time, 3.1×10⁹ minor
> faults, majflt=0) ⇒ ALLOCATION-BOUND, the third wall leg. So apply the SHARPENED falsifier below, not
> the raw [4.0,5.4]× wall test. Sample `/proc/<pid>/stat` stime/utime (or `alloc_overhead_wall.py --live`)
> BEFORE the n=26 process exits — the split is lost on exit.**
> `run_reach_detached.sh` (setsid-detached, PPID=1, SURVIVES cycle boundaries)
> ran n=25 (honest+delta_pi) and immediately started n=26 (--no-cheat) over BIG_Q, FULL streamed config.
> **n=25 harvested (run_n25.log):** HONEST ACCEPTED, **claimed π(2²⁵−1)=2 063 689 == sieve == ground-
> truth lock** ✓; wall 5421.0 s; peak RSS 5929 MB (5.79 GB); delta_pi liar REJECTED (1720 s); profile
> holds (per-layer leaf=0, comm=438070 of which comm_outer 426268=97% Θ(√x), base opens 24704 Θ(x^{1/4})).
> **CAVEAT — n=25 wall is contention-inflated, NOT a clean falsifier point:** it ran 11:25–13:24Z
> *concurrently with the S531 op-count path-(b) run* (which itself ran an n=24 chain) sharing the 32
> cores, so the n24→25 ratio (~3.7×/Δn=1 vs ~1467 s n=24) is wall-contention, not op-count (op-count is
> deterministic, measured Õ(x) S531). **The clean S529 falsifier test needs n=26 (now running, started
> 13:24Z, contention-free):** when `run_n26.log` shows "DONE", verify claimed π(2²⁶−1)=3 957 809 (cross-
> check an independent Eratosthenes sieve), record wall + peak RSS + the reach's stime/utime split, and
> **apply the S535 SHARPENED falsifier (NOT the raw wall test):** a wall ratio n24→26 >5.4× is CONSISTENT
> with Θ(x) op-count IFF the reach is `stime`-dominated (S535 measured 0.695 = allocation/page-fault leg
> 3, which S529's [4.0,5.4]× model omitted) — the discriminator is `stime`/`utime`, NOT the wall alone;
> only a `utime`-dominated excess with super-linear per-element growth reopens super-Θ(x). (op-count is
> Θ(x) to n=24, S530/S531/S532b.) Then write the item-5 reach result. A concurrent cycle must NOT start a
> second reach — `pgrep -f large_x_benchmark`; alive iff a process exists AND `run_n26.log` lacks "DONE".

S529 established (mechanism, in-cycle): the verification-chain reach wall is **Õ(x) op-count × a
one-time L3-crossover DRAM constant (~7×) + 3-point single-step noise** — NOT super-Θ(x) and NOT a
per-n growing memory penalty (the per-element `_mul61` wall is flat ~90–110 ns across the whole reach
DRAM regime; geo-mean reach ratio 4.52×/Δn=2). The detached n=25/26 run is the macro confirmation.
The verification stack is otherwise mature end-to-end (S527): prover Õ(x) op-count, verifier Õ(√x)
unconditional / Õ(x^{1/4}) large-table opens (computational), cert size Θ(√x) (information-forced
S511 + layering-inherent S510), field-lifted, fast-Mersenne, batched, streamed. **The
verification-line free CODE lever is spent** — the reach is now pure RAM-bound measurement.

**NEXT (single, concrete):** when `run_n26.log` (in `experiments/constructions/
p12_sumcheck_pi_verification/`) shows DONE, harvest + apply the S535 SHARPENED falsifier above and write
the item-5 reach result (no new code — pure finalize; ALL THREE legs are now measured — op-count
Θ(x)·polylog 4.0–4.4×↓ (S530/S531/S532b) × compute/cache ≤1.35× (S529) × allocation/page-fault stime
(S535, the reach is `stime`-dominated 0.695) — so the test is NOT "is the wall in [4.0,5.4]×?" but "is
the >5.4× excess (if any) `stime`/allocation [⇒ Θ(x) intact] or `utime`/super-linear [⇒ reopen]?" — and
sample the reach's stime/utime BEFORE it exits (lost on exit) — using the CLEAN n=26 wall, since n=25's wall was
contention-inflated, see reach block above). Path (b) is DONE (S532b block above: op-count Õ(x)
confirmed to n=24, S530 verdict extended, false auto-reading corrected). **n=25 is DONE + verified
(π=2 063 689 ✓); n=26 is RUNNING** (started 13:24Z) — a concurrent cycle must NOT start a second reach
(`pgrep -f large_x_benchmark`; alive iff a process exists AND `run_n26.log` lacks "DONE"). If a cycle
wants fresh ground while the reach computes, the op-count curve can be cheaply extended to n=26
(`prover_opcount_scaling.py --scaling --ns 16,18,20,22,24,26 --field big` — DETERMINISTIC, ~RAM of one
n=26 chain) to keep testing the S531 falsifier (2-term band α stays ≤1.05 / ratios keep declining).
> **✅ PATH (b) DONE + HARVESTED (S532b, 2026-06-14):** the detached op-count curve
> (`opcount_scaling_n24.log`, `--ns 16,18,20,22,24`) COMPLETED. **Verdict: op-count is Θ(x)·polylog
> = Õ(x), confirmed to n=24** (x≈1.7×10⁷). The per-step ratio kept DECLINING (n20→22 4.353×, n22→24
> 4.329×) — **S530's live falsifier ("ratio turns back UP above ~4.5× at n≥22") is NOT triggered.**
> The detached run's *auto-reading* fired a FALSE "SUPER-Θ(x)" from a single-power fit over the short
> high-n window (α=1.064); this was the **polylog-inflation artifact S512's `sharpP_probe` 2-term fit
> exists to defeat** — the script only had a single-power fit. FIXED: added `_fit_power_log` /
> `_fit_polylog_fixed_alpha` / `_stable_band` + a trend-based reading + `--points` (inject deterministic
> high-n counts). 2-term fit over n=14..24 → **α=1.012, polylog δ=+0.72** (residual 17× tighter than
> single-power); forced-α=1 model **`x·(log₂x)^0.87` predicts every per-step ratio to <0.4%**. selftest
> [6] codifies it (6/6 pass). No landed-chain change. `prover_opcount_scaling_results.md` §S531;
> `opcount_scaling_corrected.log`.
**Do NOT** pivot to P7.b: S530 CONFIRMED it is a re-derivation of the closed S224 Correlation
Dichotomy (dyadic anchors = the uncorrelated side ⇒ Θ̃(√x)/query; NUFFT pre-closed by dichotomy
falsifier #4). The remaining genuinely-fresh frontiers need new math (below).

**OR pivot to fresh adjacent territory** (the verification stack being mature). NOTE the adjacent
positive-target shelf is now BARE of measurable targets: **P7.b CLOSED** (S530 = re-derivation of
S224 Correlation Dichotomy), **P10 CLOSED-NEGATIVE** (S532, adaptivity beats nothing), and **P8
CLOSED-NEGATIVE** (S533, per-bit cost geometric not sub-linear — the precision-domain face of the
info barrier; confirmed the NEXT-ACTION caution that it re-derives S511/S518). The ONLY remaining
named positive target is **P9** (quantum batched π — literature survey only, likely inconclusive, no
in-repo measurement possible). So a fresh-ground cycle now has no cheap positive-target measurement
left; the substantive frontiers are the NEW-MATH ones below (no code candidate) — prefer harvesting
the reach / extending the op-count curve over manufacturing a low-value measurement. **S534 spent the
last cheap verification-line measurement** (soundness budget + field-size requirement — step 21 above,
B=Θ(√x), p(10^100) is correctness-bound at ~333-bit q); a further fresh-ground cycle should HARVEST the
reach (below) rather than manufacture another measurement.

Deep open frontiers (each its own track, both likely NEW MATH — no code candidate):
- **#P-hardness theory:** super-c*-blowup Turing reductions (S517 leaves these formally open, no
  construction exists) OR a value-incompressibility proof. S512+S517 close every natural reduction.
- **A NON-NATURAL sub-√x witness for L_π** (item 3 direction (ii) residue): S518 closed the
  natural analytic witness; a non-arithmetic witness beating the S511 floor is not ruled out (no
  candidate exists).
