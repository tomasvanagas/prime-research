# Novel Finding: Exact π(x) Is Succinctly Verifiable — the Computation/Verification Dichotomy

**Status:** Novel construction + working protocols (S491). First attack
in the project (and, per S491 + S48 literature searches, first in the
literature) on the VERIFICATION complexity of the prime-counting
function, as opposed to its computation complexity. Implemented,
tested, adversarially exercised:
`experiments/constructions/p12_sumcheck_pi_verification/`.

## The protocol stack (all unconditionally sound; all measured)

| protocol | verifies | verifier | prover | status |
|---|---|---|---|---|
| wheel sum-check | Φ(x, polylog wheel) | polylog | Θ̃(x) | built; 100/100 cheaters rejected |
| Lucy base | exact π(x) | Õ(x/log x) | Õ(x^{3/2}) | built; π(2²⁰−1) = 82025 verified; 50/50 rejected |
| Lucy + delegated wiring | exact π(x) | **Õ(√x)** | Õ(x^{3/2}) | built; crossover p ≈ 130 measured |
| GKR-over-AKS | exact π(x) | polylog | poly(x) | literature corollary (GKR'08 + AKS); impractical constants |

Headline measured row: exact π(1048575) = 82025 verified through 172
sieve layers, verifier < 1 s, transcript ~110 KB, soundness
information-theoretic (no crypto, no setup), with self-consistent
cheating provers — including DP-table corruption recomputed coherently
downstream, and false wiring claims — rejected in every trial.

## Summary

Every barrier catalogued by this project — the √x hard zone, the
GUE-random incompressible digits, the 730+ closed computation routes —
constrains *computing* π(x) from scratch. None of it constrains
*checking* a claimed value produced by a prover who did the work. That
asymmetry is real and exploitable:

> **Protocol (S491).** An untrusted prover who runs the Lucy_Hedgehog
> sieve DP can convince a verifier of the EXACT value of π(x) with
> unconditional (information-theoretic) soundness — no hardness
> assumptions, no trusted setup. The verifier checks each of the
> K = π(√x) sieve layers by a degree-≤3 sum-check and never touches
> any Θ(x)-size object. Verifier cost: Õ(n·Σ_{p≤√x} p) field ops with
> automaton wiring (designed improvement to Õ(√x) via carry-chain
> wiring); communication Õ(K·n); soundness error ≤ ~7nK/q.

Measured (q = 2³¹−1): π(65535) = 6542, π(262143) = 23000,
π(1048575) = 82025 all verified exactly; verifier 69 ms → 0.9 s while
x grows 16×; transcripts 28–110 KB; **50/50 cheating provers rejected**,
including self-consistent liars that corrupt one DP entry and recompute
all downstream layers coherently.

## The two structural observations that make it work

1. **The sieve recursion is made of small-automaton predicates.** The
   verifier's bottleneck in any sum-check is evaluating the multilinear
   extension (MLE) of the summand at a random field point. The MLE of
   the *primality indicator* is dense — 0.71 nonzero-coefficient
   fraction, 2^{n/2} coefficients (S48 / E2.1) — which is precisely why
   sum-check-for-computing-π(x) died at S18/S48. But the Lucy recurrence
   `S_i(v) = S_{i−1}(v) − [v ≥ p_i²]·(S_{i−1}(⌊v/p_i⌋) − (i−1))`
   contains no primality predicate: floor-division ⌊v/p⌋ is a p-state
   long-division automaton, the threshold [v ≥ p²] a 3-state comparator,
   and any width-W bit-automaton predicate has an MLE evaluable at an
   arbitrary point in O(n·W) by a transfer-matrix DP. **The sieve is
   verification-friendly because it decomposes primality into automaton
   steps; primality as a monolithic predicate is not.** This is a new
   counterpoint to E2.1: the same density wall that blocks computing
   blocks nothing once the predicate is factored through the recursion.
2. **The Lucy correction scalar is public.** S_{i−1}(p_i − 1) = i − 1
   (the primes below p_i), so each layer needs no extra prover claim
   for its correction term — the layer reductions chain with only two
   point-claims (merged by line restriction) and terminate in the
   closed form S̃_0(z) = Int(z) − 1 + Π(1−z_j).

## Corollaries

- **Succinct certificates for p(n).** "p(n) = y" is certified by (i) an
  ECPP/Pratt primality certificate for y, plus (ii) this protocol run
  on π(y) = n. Kilobyte-scale, unconditionally sound up to the field
  soundness error. Today's trust mechanism for record π(x) values is
  independent recomputation by a different algorithm; this replaces
  months of recomputation with a transcript.
- **Certified comparisons in the hard zone.** The decision problem
  [p(n) ≤ x?] ⟺ [π(x) ≥ n?] is worst-case as hard as computing p(n)
  (binary search, E6.6), and inside the √x-wide hard window around
  R⁻¹(n) each comparison costs a full π(x) computation. But it is
  cheaply *verifiable*: a prover can certify each comparison bit.
  Empirical geography of the decision problem measured in
  `experiments/analytic/guess_comparison_oracle/`: against the biased
  guess li⁻¹(n) the comparison bit is constant (30/30, the
  Rubinstein–Sarnak bias, predictable up to the Skewes scale ~10³¹⁶);
  against the centered guess R⁻¹(n) it is an empirical coin flip
  (10/30 vs 20/30) — the better the guess, the closer the bit is to
  fair, which is the information barrier restated.
- **The verification-complexity landscape (corrected S491-late).** In
  principle, GKR 2008 over the data-parallel uniform primality circuit
  (Σ over w ≤ x of an AKS copy — polylog depth, log-space uniform)
  already gives an unconditionally-sound interactive proof with polylog
  verifier and poly(x) prover — a literature corollary, NOT novel, and
  with astronomically large polylog constants and a nontrivial
  wiring-uniformity verification. (An earlier draft of this file
  posed "polylog or √x floor?" as open; the in-principle answer is
  polylog, by GKR.) The honest open questions are: (a) a *practical*
  polylog-verifier protocol — the S491 Lucy protocol runs in
  milliseconds but its wiring totals Õ(x/log x); **the delegated-wiring
  upgrade (built and measured S491-late,
  `lucy_dp_delegated_wiring.py`) reaches Õ(√x)**: the division
  automaton's DP is itself GKR-delegated, the verifier evaluating the
  transition relation s′ = 2s + a − b·p via a width-4 carry automaton
  in O(log p) — measured 3.2× verifier growth over an 84× p-range (vs
  34× for the automaton baseline, whose linear-in-p cost the
  width-dichotomy theorem proves unavoidable without delegation);
  crossover at p ≈ 130. The remaining √x → polylog gap is now exactly
  the K = π(√x) sequential layer structure of the sieve itself;
  (b) succinct *non-interactive* certificates unconditionally — made
  precise (S491-late): a polylog-size unconditionally-checkable
  certificate for π(x) = c is exactly the statement
  **L_π = {(x, c) : π(x) = c} ∈ NP** in input length n = log x.
  No construction is known (even the "≥ c" direction lacks a succinct
  witness — listing primes costs ~x bits; Pratt/ECPP certify one
  prime, not a count), and no barrier either: for #P-complete
  functions, value-certificates in NP would force NP = coNP, but π(x)
  is not known #P-hard (proving THAT would itself be a major result,
  and is a cleaner formal target than the algorithmic question).
  L_π ∈ NP is the proof-theoretic shadow of the headline problem:
  polylog computation would trivially imply polylog certificates, so
  certificate lower bounds are the easier-to-attack face of the same
  wall; (c) the constructibility dichotomy below.

## The width spectrum (S491 brainstorm follow-up, measured)

`experiments/constructions/automaton_width_dichotomy/` upgrades the
project's OBDD-width edge from a computation lower bound to an exact
verification-cost measure: given an explicit read-once branching
program with cut-widths (W_j), the MLE is evaluable at any field point
in O(Σ_j W_j) (transfer DP — Theorem T1, implemented). Measured at
x = 2²⁴ with matched-random control:

- **Division wiring [u = ⌊v/p⌋] has Nerode width exactly 2p+1** —
  the S491 protocol's O(n·p) wiring check is optimal among explicit
  automatons; prover-supplied carry witnesses are *necessary* (not just
  sufficient) to go below.
- **χ_P's width profile = incompressible plateau + admissibility
  cliff**: W_j = 2^j exactly through j = 16 (every aligned 256-window
  of primes distinct — the global pseudorandom side), then a crash to
  44× below random (107 vs 4707 distinct aligned 16-windows — wheel/HL
  local rigidity made quantitative). The cut-profile invariant
  D_k(x) = #distinct aligned k-windows of χ_P is new to the project
  and bridges its global-randomness and local-structure edge families
  in one curve.
- **Middle-cut width is maximal (2^{n/2} = 4096 > E2.1 rank 2050)** —
  refines the rank edge: the factor-2 deficiency is linear-algebraic,
  not informational.
- **One-shot sum-check verification of π(x) is dead twice over**:
  Σ_j W_j(χ_P) ~ x^c with c drifting 0.73 → 0.77 upward (peak-crossover
  argument gives c → 1), AND the minimal program is non-constructible
  without the primes. Verification must consume primality through
  explicit small-width sieve predicates — layering is structurally
  forced, which is exactly what the Lucy protocol does.
- **The cliff is exactly Hardy-Littlewood (S491 continuation,
  `experiments/analytic/local_pattern_census/`).** W_j = D_k(2ⁿ)
  (k = 2^{n−j}) by definition, and the census saturates at
  **D_k(x) = A_k + 1** — the exactly-enumerable count of HL-admissible
  aligned patterns, plus the single m = 0 window containing the small
  primes. Verified to closure: k = 8 (14 = 13+1), k = 16 (107 = 106+1,
  all admissible realized by x = 2¹⁸, stable through 2²⁸); k = 32
  converging (3385 of A₃₂+1 = 3574 at 2²⁶; the 189 missing are exactly
  the weight ≥ 6 dense constellations — the census deficit is a
  finite-x probe of k-tuple uniformity). Max weights 3/5/9 match the
  classical H(k) narrowest-tuple function. **The prime indicator's
  width spectrum = global-pseudorandom plateau + admissible-pattern
  entropy curve**, with the crossover set by the (open) entropy
  constant lim log₂(A_k)/k (measured 0.46 → 0.37 at k = 8 → 32).

## Honest scope

- No speedup for computing π(x); the prover pays the full sieve cost.
  This is an adjacent-problem partial-positive in the S224 mold —
  computation unchanged, trust exponentially cheaper.
- The protocol composes standard tools (LFKN 1992 sum-check; GKR-style
  layer reductions; CMT 2012 prover techniques). The novelty is the
  application target (no prior verification protocol for π(x) exists),
  the automaton-MLE treatment of the sieve wiring, and the public-
  scalar observation that makes the Lucy chain close without auxiliary
  claims.
- Demo field q = 2³¹−1 (soundness ~10⁻⁶ at n=16). **Field lift (S498,
  corrected):** the multiply-trace identity U·a+R−X=0 is sound only while
  the field's CHARACTERISTIC exceeds the integers it relates (U·a ~ x),
  i.e. q > x, roughly n ≲ 30. The lift is therefore to a larger-
  characteristic PRIME (q ~ 2⁶¹−1 for n ≲ 60, generic ~128-bit prime for
  cryptographic x) — NOT an extension field. An extension field F_{q²} has
  q² *elements* but characteristic q, so an aliasing forgery (off by k·q)
  embeds as the SAME field element as the truth and is undetectable; the
  larger element count only shrinks Schwartz–Zippel error, an orthogonal
  concern. Demonstrated end-to-end in `compressed_prover_mult_trace.py`
  (`--field`, `--alias-demo`, `--refute-q2`). The genuine alternative to a
  bigger prime is a small-field schoolbook carry trace. **The lift is now
  threaded through the WHOLE compressed chain (S499):** `compressed_layer.py`
  `run_chain` and the delegated wiring (`lucy_dp_delegated_wiring.py`) carry a
  prime-modulus `q` (`--field {q,big,small}`), default `q=2³¹−1` bit-identical,
  `q=2⁶¹−1` over exact-int object arrays — the chain over `big` gives
  `claimed π==sieve` with all cheats rejected, and the wrap-around alias the demo
  prime admits is rejected by the lift in the chain's exact trace config. The
  remaining gap to a large-n run is speed, not soundness — but (**S500**) the
  obvious speed fix is not the bottleneck: the numpy `2⁶¹−1` mulmod
  (`_mul61`/`_sum61`, `2⁶¹≡1` fold) is built and threaded chain-wide
  bit-identically, yet it is a *large-array-only* win (the single-`√x`-cube
  certifier gains 1.7–3.8× at n=24–28; the **full chain is slower**, being
  `π(√x)` layers of many small cubes — *operation-count-bound*, not
  width-bound). The pure-Python Lucy DP is also negligible (`<0.1%` of
  `run_chain`). So the real lever for a fast large-x chain run is reducing the
  *count* of small per-layer field-ops (cross-cube batching), kept opt-in
  (`--fast-big`) with BIG_Q on the object reference by default.
  Demo prover is table-based Θ(2ⁿ/layer); production prover works the O(√x)
  compressed Lucy state at Õ(x) total.
- **Leaf openings (S505).** Each layer's claim is grounded by opening a committed
  O(√x)-size MLE at a point. The demo closed these by direct `mle_eval`
  (O(√x) each → O(K·√x) ≈ Õ(x) over the chain, the last √x-/p-linear verifier
  term). `leaf_open.py` builds the standard **sum-check opening** — `S~(pt)=Σ_w
  S[w]eq~(pt,w)` is itself a degree-2 sum-check — reducing an opening at `pt` to a
  polylog transcript (verifier O(nb)) plus ONE residual claim at a fresh point.
  `run_chain(--pcs)` threads each fold's residual as the next layer's carried
  claim, discharged transitively at the S₀ base — **unconditional** (no
  commitment; soundness rides the GKR layer reductions). Measured: the standalone
  opening verifier is flat in 2^nb (6.5×→432× over `mle_eval` as nb grows;
  21×→3605× over BIG_Q); the chain integration is a **stable ~1.2–1.36× verifier
  win**, NOT yet asymptotic — the remaining per-layer O(2^nb) term is the trace
  region's nb Ub-bit-table openings (per-layer witness data, no carried claim to
  thread to), which want the batched-trace integration. A real multilinear
  commitment (Brakedown/FRI-free) would discharge all residuals succinctly under a
  hashing assumption, trading unconditionality for full succinctness.
- **Batched Ub openings (S506) — the per-layer verifier is now leaf-eval-free.**
  That one remaining per-layer O(2^nb) term is closed. The nb Ub-bit-table openings
  `Ub_{Lv−nb+k}~(r_C)` (the B2 routing that pins g1_trace to the certified quotient
  u_e) are EXACTLY the tables the batched trace zero-test (`batched_trace`,
  `verify_constraints_batched`) already stacks and commits along the layer axis, so
  the discharge SHARES that commitment. Each layer defers `(l, r_C^l, [ub_{l,k}])`
  (the prover supplies the scalars; the verifier folds them into B2's `expect` in
  O(nb)) and all K·nb discharge in **ONE degree-2 sum-check**: with γ←F_q, β=γ^nb the
  weight β^l·γ^k=γ^{l·nb+k} is a distinct power, and the identity factorizes to a
  single inner product `Σ_w B[w]·C[w]` over the (Lk+nb)-cube — B the verifier-anchored
  per-layer eq weights (the soundness anchor, recomputed in O(K(Lk+nb))), C the γ-fold
  of the committed low-nb Ub tables (its opening taken from the sum-check's folded
  scalar, exactly as the zero-test trusts its bit scalars). `verify_ub_openings_batched`;
  wired behind `run_chain(--pcs --batch_trace --batch_ub)`. **Measured: the verifier-
  side O(2^nb) Ub-leaf eval count drops K·nb→0 at every n** (moved to the prover); the
  per-layer verifier is leaf-eval-free, so the only O(√x) verifier work left is
  one-time — the two S₀ base closes plus the batched discharges. **The unconditional
  Õ(√x) end-to-end verifier is now structurally real.** Honest scope: the wall-clock
  t_verifier win is a modest ~1.1–1.2× and flat in n (a vectorised numpy `mle_eval`
  over a √x-size array is cheap relative to the Python-loop wiring/eq recomputes that
  dominate the measured t_verifier), so the asymptotic O(K·nb·2^nb)→one-time gain is in
  **op-count, not measured wall-clock** at reachable n — correcting the S505 prediction
  that the verifier ratio would grow with n. Only a real outer commitment for the two
  one-time S₀ base closes remains (computational, not unconditional).
- **Tensor polynomial commitment for the base closes (S508) — the whole-chain
  verifier is now sub-√x.** That last one-time term is closed. `pcs_commit.py` is a
  hash-based multilinear PCS (Ligero/Brakedown tensor style): reshape the 2^nb table
  into an r×k matrix (r=2^{n1}, k=2^{n2}, n1+n2=nb), RS-encode each row, Merkle-commit
  the N columns (sha256 = CRH/RO stand-in). Using `eq~(pt,w)=eq~(pt_hi,w_hi)·eq~(pt_lo,
  w_lo)` so `S~(pt)=a^T M b` (a=eq_table(pt_hi), b=eq_table(pt_lo)), the **evaluation-
  proof verifier is O(t·(r+k)) = O(x^{1/4}·polylog)** — sub-√x: it reads the two
  combined messages v=a^T M, w=ρ^T M (len k), checks the eval binding ⟨v,b⟩=claimed,
  and at t Fiat–Shamir columns checks Merkle membership + the linear-code consistency
  Enc(v)[c]=⟨a,Mhat[:,c]⟩, Enc(w)[c]=⟨ρ,Mhat[:,c]⟩. The full 2^nb table never touches
  the verifier. Wired behind `run_chain(--commit-base)`; with `--pcs --batch_trace
  --batch_ub --commit-base` the **whole-chain verifier is sub-√x**. Measured
  (`--bench-verifier-ops`, config d): leading-term exponent **α_ops = 0.258 ≈ 0.25
  (Θ(x^{1/4}))**, vs config c's 0.500 (Θ(√x)) — the milestone drop Θ(x)→Θ(√x)→Θ(x^{1/4}).
  Verdict unchanged over q AND BIG_Q (honest π==sieve, delta_pi + self-consistent liar
  rejected); cheating-prover tests reject wrong claim / forged opening / tampered
  codeword row / tampered revealed column. **Honest scope:** this trades the otherwise-
  unconditional Õ(√x) verifier for full succinctness under the hash assumption (the rest
  of the verifier stays unconditional); the win is the SLOPE — the absolute op count
  only drops below the direct 2^nb close at nb≈14 (n≈28), the t-query constant dominating
  below that (the S506 crossover pattern); and RS's N≤q point requirement is a demo
  artifact that a field-size-free linear code (Brakedown's expander code) removes behind
  the same interface.
- **Certificate SIZE characterized — the Õ(√x) comm wall localized (S509).** With the
  verifier's large-table work down to Õ(x^{1/4}), the binding Õ(√x) quantity is the
  certificate SIZE itself (the Fiat–Shamir transcript, `stats["comm"]`, slope ~0.5 in
  every config). `certificate_profile.py` measures the full succinct config
  (delegate+structured+pcs+batch_trace+batch_ub+**batch_wiring**+commit_base) over
  n∈{8..18}: cert size **α=0.473 (Θ(√x))**, verifier ops **α=0.258 (Θ(x^{1/4}))**, prover
  wall α=0.41 (overhead-bound at reachable n; asymptotic element-work Õ(x)). `run_chain`
  now partitions `comm` by source (additive snapshots, sum==comm asserted); the Õ(√x) is
  **DOMINATED by `comm_outer` (α=0.522, 60%→90%→100% of comm) = the K=π(√x) SEQUENTIAL
  outer two-sided layer reductions** (each O(nb²)=O(log²x) round scalars × K=π(√x) ⇒
  Θ(√x·polylog)); every batched discharge (trace/wiring/Ub) is polylog (α≤0.18), the base
  commit Θ(x^{1/4}). Batching the wiring (S503/S504) cuts comm 5.3× at n=18 — the reason it
  is in the headline config. **Honest sharpening of S508:** the verifier must READ the
  Õ(√x) proof (each transcript scalar is an O(1) round check), so **total verification is
  Θ(√x)** — bounded below by the proof size; the S508 Õ(x^{1/4}) is the large-table-eval
  SUB-term only. **Membership result:** `L_π = {(x,c) : π(x)=c}` has, under a CRH/RO
  (Fiat–Shamir of the chain), a NON-INTERACTIVE certificate of size **Õ(√x)** field
  elements (Õ(√x·log²x) bits) with **Õ(x^{1/4})** large-table opens. Compressing the size
  (hence total verification) below √x to polylog requires batching the K SEQUENTIAL outer
  reductions — but each layer's OUTPUT claim is the next layer's INPUT claim (a dependency
  chain), unlike the independent per-layer trace/Ub/wiring witnesses that WERE batchable.
  That is exactly the ⌊v/p⌋-semigroup compression — the layer-batching closed negative
  (open problem 2). **The certificate-size wall and the prover-time wall are the same
  sequential-sieve wall, on the two faces of the problem.**
- **Merged-layer comm probe — the √x is layering-inherent (S510).** S509 left one
  quantitative question: is the Õ(√x) certificate size inherent to the layering, or
  just to the unbatched chain? `merged_layer_comm.py` answers it with a sound,
  cheat-tested depth-j merged-layer reduction (the 2^j-subset inclusion-exclusion
  operator, == direct j-layer composition), measured over n∈{8..16}, j∈{1,2,4}, both
  q and BIG_Q (identical comm). **The chain total comm exponent stays ~0.5 (Θ(√x)) for
  every j in both a DIRECT mode (comm carries the 2^j fill; total comm even grows past
  j=2) and a fully-BATCHED mode (the 2^j leaves the transcript → comm polylog, but
  moves into verifier WORK ~2^j and the prover, comm floored at ⌈K/j⌉(j+nb)→K=√x).**
  The prover fill-table total ⌈K/j⌉·2^j·2^nb trends toward x^{3/2} as j→log K — open
  problem 2's m·K²/2 reproduced on this construction. So no merge depth drops the comm
  exponent below 0.5: the Õ(√x) certificate is layering-inherent, and the cert-size,
  verifier-work, and prover-time walls are all the one ⌊v/p⌋-semigroup wall. This
  strengthens the membership result to "the sieve certificate is Θ(√x) at every merge
  depth"; a sub-√x certificate, if one exists, cannot come from merging layers.

## Open problems (the verification program, consolidated S491-final)

1. **Compressed prover** — make the prover pay Õ(x) (polylog overhead
   over plain Lucy) instead of Õ(x^{3/2}). The large-value side forces
   d-addressing whose cross-map u = ⌊x/(dp)⌋ is variable×variable
   multiplication ⟹ a trace-table (AIR-style) GKR. **Step 1 built and
   tested (S492):** `compressed_prover_mult_trace.py` — a batched
   primitive verifying C = Σ_d eq(ρ,d)·S̃(⌊X/(dp)⌋) with every quotient
   certified through prover trace tables (a degree-3 constraint
   zero-test — multiply+remainder identity, bit-range recompositions,
   remainder bound — plus a routing lookup), verifier polylog O(m·Lv)
   and never dividing, prover on the D≈√x cube; 6 cheat classes
   rejected, verifier flat in batch size. Remaining: integrate into the
   layer protocol (small/large-side dispatch; the affine index map for
   dp≤√x; chaining and line-batching the S̃ claims).
2. **Layer batching (closed negative, S491-final; COMM face measured
   S510).** Merging the K = π(√x) sequential layers by a balanced
   product tree fails on fill-in: each layer matrix has fill ≤ 3 per row
   on the compressed state space, but a j-layer merge has fill ~2^j per
   row (one entry per subset-product ⌊v/(p_{i₁}···p_{iₛ})⌋; K < |V| so
   fill never saturates). Sparse-multiplication cost at tree level ℓ is
   K·m·2^{ℓ−1} — geometric, top-dominated: **m·K²/2 ≈ x^{3/2}/ln²x**,
   no better than the unbatched prover. **S510
   (`merged_layer_comm.py`) measured the COMM of a sound depth-j
   merged-layer reduction directly** — collapsing j consecutive
   small-side (value→value ⌊v/p⌋) Lucy layers into one reduction over
   the 2^j-subset inclusion-exclusion operator
   M[S](v)=Σ_{T⊆[j]}(−1)^|T| GateProd_T(v) S(⌊v/m_T⌋) (== direct
   composition, m_T=Π_{k∈T}p_k), in two sound modes: DIRECT (the 2^j
   subset openings) and BATCHED (the 2^j stacked on a j-bit axis, one
   sum-check). **The chain total comm exponent stays ~0.5 (Θ(√x)) for
   every merge depth j∈{1,2,4} in BOTH modes** — DIRECT comm carries the
   2^j fill explicitly (and the total comm *grows* past j=2), BATCHED
   pushes the 2^j out of the transcript (comm polylog) but into verifier
   WORK (measured ~2^j) and the prover, the comm floored at
   ⌈K/j⌉(j+nb)→K=√x; the prover fill total still → x^{3/2}. So polylog
   verification still requires algebraic compression of the semigroup
   generated by v ↦ ⌊v/p⌋ — open — and the √x is **inherent to the
   layering on all three faces (certificate size, verifier work, prover
   time), not an artifact of the unbatched chain**.
3. **L_π ∈ NP?** (see above) — succinct unconditional non-interactive
   certificates; no construction, no barrier. The complementary formal
   target: **is π(x) #P-hard?** Either answer would be major; hardness
   would explain the project's whole closure landscape. **Concrete upper
   bound (S509):** L_π has a *computational* (CRH/RO) non-interactive
   certificate of size **Õ(√x)** — far short of polylog, and the Õ(√x) is
   provably (measured) the K=π(√x) sequential layer reductions, i.e. the
   open-problem-2 semigroup wall. So the unconditional-polylog question and
   the layer-batching question are the same wall; a sub-√x certificate (even
   computational) would already be new. **S510 sharpens this: merging j
   consecutive layers does NOT drop the comm exponent below 0.5 (measured,
   both a direct and a fully-batched merged-layer reduction), so a sub-√x
   certificate cannot come from layer merging — the construction face of this
   upper bound is exhausted. The complementary lower-bound question (is sub-√x
   information-theoretically impossible?) is the natural successor.**
   **S511 answers the information half of that successor.** Measuring the joint
   hard information of the K=π(√x) sieve checkpoints `{φ(x,a)}` (the Legendre
   partial sieve, the chain's per-layer survivor count) at integer precision
   against a polylog smooth predictor: the residual matrix is FULL integer-rank
   once adequately sampled (rank/K → 0.97 at window W≫K; α_rank ≈ α_K ≈ √x), so
   the checkpoints are integer-INDEPENDENT and the joint information is
   **Θ(√x)·polylog**. Hence the √x certificate is **information-forced for any
   sieve-reconstructing verifier** — polylog and sub-√x certs are both ruled out
   for that class, upgrading S510's *construction*-inherent √x to an
   *information*-inherent one. (A single π(x) value remains O(log x) bits, S36:
   the barrier is JOINT across layers, not per-value.) **Scope:** this binds the
   sieve-reconstruction class only; a fundamentally different witness (cf.
   factoring) is not ruled out, so a UNIVERSAL sub-√x lower bound still routes
   through the **#P-hardness** question — now the sole remaining lower-bound
   lever. See `experiments/constructions/p12_sumcheck_pi_verification/
   cert_incompressibility.py`.
4. **A_k entropy law at k = 128** — the census conjecture
   A_k = e^{(1+o(1))π(k)} is enumeration-verified to k = 64 (clean
   2-power family, geometric convergence to constant ≈ 1); k = 128 has
   A ~ e²⁵ and needs an analytic count (the surjection inclusion-
   exclusion suffers catastrophic cancellation; the DFS is
   enumeration-bound). A transfer-matrix or generating-function attack
   is the well-posed successor.
5. **Large-x benchmark** — once (1) lands: verify π(10⁹⁺) and measure
   the verifier-vs-recomputation crossover against primecount-class
   baselines; the asymptotic separation is already measured at the
   wiring level (3.2× vs 34× over an 84× p-range).

## Cross-references

E2.1 (MLE density — circumvented, not contradicted), E6.6 (binary
search), S18 + S48 CLOSED_PATHS rows (scope-corrected: they close
computation, not verification), S224 (partial-positive prototype),
OPEN_POSITIVE_TARGETS.md §P12 / ATTACK_VECTORS.md §H5 (Thread 12),
`experiments/constructions/p12_sumcheck_pi_verification/` (code, logs,
falsifiers), `experiments/analytic/guess_comparison_oracle/` (decision-
problem geography).
