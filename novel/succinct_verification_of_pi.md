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
2. **Layer batching (closed negative, S491-final).** Merging the
   K = π(√x) sequential layers by a balanced product tree fails on
   fill-in: each layer matrix has fill ≤ 3 per row on the compressed
   state space, but a j-layer merge has fill ~2^j per row (one entry
   per subset-product ⌊v/(p_{i₁}···p_{iₛ})⌋; K < |V| so fill never
   saturates). Sparse-multiplication cost at tree level ℓ is
   K·m·2^{ℓ−1} — geometric, top-dominated: **m·K²/2 ≈ x^{3/2}/ln²x**,
   no better than the unbatched prover. Polylog-depth verification
   with Õ(x) prover therefore requires algebraic compression of the
   semigroup generated by the maps v ↦ ⌊v/p⌋ — open, and the precise
   form of the project's √x sequential-sieve wall on the verification
   side.
3. **L_π ∈ NP?** (see above) — succinct unconditional non-interactive
   certificates; no construction, no barrier. The complementary formal
   target: **is π(x) #P-hard?** Either answer would be major; hardness
   would explain the project's whole closure landscape.
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
