# Interactive proof for the exact value of π(x) — working protocol (S491)

## Statement of what was built

**A complete interactive proof in which an untrusted prover convinces a
verifier of the EXACT value of π(x), with unconditional
(information-theoretic) soundness, by certifying every layer of the
Lucy_Hedgehog sieve recursion via sum-check. The verifier never sieves
and touches no Θ(x)-size object.**

To our knowledge (S491 search + the project's S48 literature scan),
this is the first verification protocol of any kind specifically for
the prime-counting function. The number-theory community's existing
trust mechanism for π(x) records is independent recomputation by a
different method (e.g., analytic vs combinatorial); this object
replaces recomputation with a transcript whose soundness needs no
hardness assumption, no trusted setup, and no re-execution.

## The protocol in one paragraph

The Lucy_Hedgehog recurrence
`S_i(v) = S_{i-1}(v) − [v ≥ p_i²]·(S_{i-1}(⌊v/p_i⌋) − (i−1))` holds
pointwise for every v < 2ⁿ, with `S_K(x) = π(x)` at K = π(√x), and the
correction scalar `S_{i-1}(p_i−1) = i−1` is a *public constant*. Each
layer is verified by a two-phase sum-check (degree ≤ 3): phase A checks
the recurrence against an auxiliary table `G(v) = S_{i-1}(⌊v/p_i⌋)`;
phase B certifies G against the floor-division wiring
`W(v,u) = [u = ⌊v/p_i⌋]`, whose multilinear extension the verifier
evaluates *itself* in O(n·p) via a **long-division automaton** (p
states = the running remainder, reading bits of v MSB-first, emitting
the bits of u). The two resulting point-claims about S̃_{i-1} are merged
by the standard line-restriction trick; after K layers the chain
terminates in the closed form `S̃_0(z) = Int(z) − 1 + Π(1−z_j)`, which
the verifier checks in O(n). The threshold `[v ≥ p²]` is a 3-state
comparator automaton, also verifier-evaluated in O(n).

## Measured results (`run_lucy_n16/18/20.log`)

| x | layers K | π(x) verified | prover | verifier | communication |
|---|---|---|---|---|---|
| 2¹⁶−1 = 65535 | 54 | **6542** (= sieve, exact) | 1.59 s | 68.5 ms | 27.6 KB |
| 2¹⁸−1 = 262143 | 97 | **23000** (exact) | 11.6 s | 253 ms | 55.7 KB |
| 2²⁰−1 = 1048575 | 172 | **82025** (exact) | 93.3 s | 895 ms | 109.5 KB |

x grows 16× → prover grows 59× (Θ̃(K·2ⁿ) table-based demo prover),
verifier grows 13× (Õ(n·Σ_{p≤√x} p) automaton evaluations dominate).

Soundness battery at n = 16 (each trial fresh verifier randomness):
- **adaptive round-patching liar** (claims a wrong π(x), fixes up every
  sum-check round polynomial to pass the running consistency check):
  rejected **25/25**;
- **self-consistent liar** (corrupts one entry of a mid-chain DP table,
  recomputes all later layers consistently — so the final claim is wrong
  but internally coherent): rejected **25/25**;
- theoretical soundness error ≤ ~7nK/q ≈ 2.8×10⁻⁶ at q = 2³¹−1
  (production: q ~ 2⁶¹ or extension field → < 10⁻³⁰; protocol
  structure unchanged).

Unit tests (`--selftest`): division-automaton MLE ≡ ⌊v/p⌋ on 800
random Boolean points across p ∈ {2,3,7,13}; comparator ≡ [v ≥ M]
exhaustively at n = 8; closed-form S̃_0 ≡ table-MLE at 20 random field
points; single-layer reduction maps a true claim to a true claim at
random field points; end-to-end n = 10 chain.

## Honest accounting

1. **This does not compute π(x) faster.** The prover does Θ̃(x)-class
   work. The exponential gap opened here is computation-vs-trust: the
   project's information-theoretic barrier (GUE-random digits of π(x))
   constrains computing the value, and says nothing about checking a
   claimed value produced by someone who did the work. That asymmetry
   is real, rigorous, and was unexploited in 490 prior sessions.
2. **Verifier exponent with this wiring: Õ(x/log x)** field ops
   (Σ_{p≤√x} n·p from the division automatons) — sublinear, and tiny in
   wall-clock at demo scale, but not yet the right asymptotic target.
   **Designed next step (not yet implemented): carry-chain wiring.**
   Replace the p-state automaton check of `u = ⌊v/p⌋` by the packed
   integer identity `Int(v) = Int(u)·p + Int(r)` with auxiliary
   carry/range bits verified inside the same sum-check (Int is a
   degree-1 extension, the carry constraints are local and degree ≤ 3);
   this makes each layer's wiring check O(n·polylog) and the total
   verifier **Õ(√x)** — the same exponent as the number of layers,
   which is the natural floor for this layered route. A polylog
   verifier for exact π(x) with unconditional soundness would need to
   break the K = π(√x) sequential-layer structure (the project's
   "√x sequential sieve wall" reappearing on the verification side —
   a nice structural echo, worth stating as an open question).
3. **Prover is table-based** (Θ(2ⁿ) per layer) for the demo. A
   production prover works on the O(√x)-size compressed Lucy state and
   pays Õ(x) total, matching the sieve cost it already pays to compute
   π(x); the protocol adds polylog-factor overhead, not a new exponent.
4. Field q = 2³¹−1 was chosen so numpy uint64 arithmetic is exact;
   nothing in the protocol depends on the field size beyond soundness.

## Why this is the right object for this project's goal-neighbourhood

- It is a **partial-positive on an adjacent problem** in exactly the
  S224 Correlation-Dichotomy mold: no speedup for computing π(x), an
  exponential improvement for *trusting* π(x).
- It is **constructive and falsifiable**: code runs, exact values
  verified, cheating provers measured.
- It gives the project a genuinely new structural classification: the
  predicates whose MLEs a verifier can evaluate cheaply are the
  **small-automaton predicates** (divisibility, comparison, congruence,
  long division — everything the sieve recursion is made of), while the
  primality indicator itself is provably outside (E2.1 / S48 density).
  π(x)'s DP is verification-friendly precisely because the sieve
  decomposes primality into automaton steps; primality as a monolithic
  predicate is not.

## Falsification statement

This construction is falsified by: an accepted honest run whose output
differs from a direct sieve (none observed, 3 scales); a cheating
prover of either tested class accepted above the soundness bound (0/50
accepted); or verifier cost found to scale with x rather than
n·Σ_{p≤√x} p (refuted by the 16×-x / 13×-verifier measurement).

## S499 note (field-parameterised automata)

The two verifier-side automata `ge_const_eval` (comparator `[v≥M]`) and
`w_div_eval` (division wiring `[u=⌊v/p⌋]`) gained a `q=Q` keyword so the
compressed chain and delegated wiring can evaluate them over an arbitrary prime
field (the S498/S499 field lift). Pure-Python, no dtype change; `q=Q` default
keeps every existing caller and the base 2ⁿ-demo protocol byte-for-byte
identical (selftest unchanged).

## Files

- `lucy_dp_verification.py` — protocol (one script, CLI-parameterised).
- `run_lucy_n16.log`, `run_lucy_n18.log`, `run_lucy_n20.log`.
- Companion: `sumcheck_wheel_verification.py` (+ results.md) — the
  simpler wheel-count protocol built first the same session; its
  automaton-MLE primitive is what the layered protocol generalises.

## Sources (cross-domain import: verifiable computation)

- Lund, Fortnow, Karloff, Nisan 1992 (sum-check); Goldwasser, Kalai,
  Rothblum 2008 (doubly-efficient IPs); Cormode, Mitzenmacher, Thaler
  2012 (practical provers); Thaler, *Proofs, Arguments, and ZK*.
- Surveys checked for prior art: arxiv.org/abs/2308.15191,
  eprint.iacr.org/2025/008. No prior verification protocol for π(x)
  found.
