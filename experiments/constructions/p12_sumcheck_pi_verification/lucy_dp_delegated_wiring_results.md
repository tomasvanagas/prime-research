# Delegated-Wiring Upgrade: the Õ(√x) Verifier for Exact π(x) (S491)

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

## Cross-references

`lucy_dp_verification.py` (base protocol; this experiment imports its
primitives), `../automaton_width_dichotomy/` (the 2p+1 tightness
theorem that makes delegation *necessary*),
novel/succinct_verification_of_pi.md (updated landscape: practical
verifier exponent now √x, matching the layer count — the remaining gap
to the in-principle GKR polylog verifier is the K = π(√x) sequential
layer structure itself), GKR 2008 / Thaler PAZK (the layer-reduction
form with eq̃-selector).
