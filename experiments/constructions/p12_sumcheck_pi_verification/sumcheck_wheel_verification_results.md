# P12 PoC — Sum-check verification of partial-sieve counts (S491)

## What was built

A complete, working interactive proof (LFKN 1992 sum-check) in which an
untrusted prover convinces a verifier of the exact value of

    Phi(2^n, P) = #{ 0 <= w < 2^n : p ∤ w  for all p in P }    (odd wheel P)

with **unconditional (information-theoretic) soundness** — no hardness
assumptions, no trusted setup. The verifier runs in
`O(n · Σ_{p∈P} p)` field operations (polylog in x = 2^n for polylog-sized
wheels); the prover does the Θ(2^n) work.

One script, CLI-parameterised: `sumcheck_wheel_verification.py`
(`--selftest`, `--n/--primes/--trials`, `--scaling`). Logs:
`run_n20.log`, `run_scaling.log`.

## The new primitive

The verifier's final check requires evaluating the multilinear extension
(MLE) of the summand at a random field point. For the **primality**
indicator χ_P this is exactly what S48 closed: the MLE of χ_P is dense
(0.71 nonzero-coefficient fraction, 2^{n/2}-size coefficients), so the
final check costs Ω(2^{n/2}) — sum-check dies for *computing* π(x).

For the **fixed-modulus divisibility** indicator δ_p(w) = [w ≡ 0 mod p]
the situation is structurally different: δ_p is recognised by a p-state
automaton reading bits MSB-first, so its MLE is evaluable at an
*arbitrary* field point by a transfer-matrix DP over residue states in
`O(n·p)` time (function `delta_mle_eval`). The summand
`g(z) = Π_{p∈P}(1 − δ̃_p(z))` has per-variable degree |P| and its
restriction to the Boolean cube is the coprime-to-wheel indicator.

A second trick makes the prover near-linear instead of naively
quadratic: for a Boolean suffix of value t (m bits remaining), the
automaton from state s ends at `(s·2^m + t) mod p`, so the suffix
δ̃-value collapses to a single gather `u[s*]`, `s* = −t·2^{−m} mod p`
(`prover_round_poly`).

**General statement of the primitive (new to the project):** any
predicate computed by a width-W bit-reading automaton admits an MLE
evaluable in `O(n·W)` per point, hence sum-check counting with verifier
cost `O(n·(k + Σ W_i))` for products of k such predicates. Divisibility
by fixed p (W = p), congruence classes, digit predicates qualify;
primality itself does NOT (S48 / E2.1 density — that obstruction is
confirmed, not contradicted, by this construction).

## Results

### Correctness + soundness (n = 20, P = {3,5,7,11,13}, q = 2³¹−1)

| check | outcome |
|---|---|
| honest run | ACCEPTED; claimed Φ = 402250 = direct sieve count (exact) |
| adaptive cheating prover (claim Φ+1, per-round additive patch) | rejected **100/100** trials |
| theoretical soundness error | ≤ n·k/q = 4.66×10⁻⁸ |
| communication | 120 field elements = 480 bytes |

### Scaling sweep (`run_scaling.log`)

| n | k (wheel) | prover | verifier | direct sieve | comm |
|---|---|---|---|---|---|
| 12 | 3 | 0.013 s | 0.86 ms | ~0 s | 192 B |
| 16 | 5 | 0.080 s | 1.89 ms | ~0 s | 384 B |
| 20 | 5 | 0.522 s | 2.70 ms | 0.001 s | 480 B |
| 22 | 10 | 6.329 s | 5.64 ms | 0.004 s | 968 B |

Prover tracks Θ(2^n·k²); verifier tracks O(n·Σp) — across a 1024×
range in x the verifier grows ~3× (driven by n and the wheel, not x).

## Honest cost accounting — what this does and does not give

1. **At PoC scale the direct sieve is faster than the verifier.** The
   verification win is architectural, not a small-x wall-clock win: the
   verifier never touches the 2^n integers, so at large x (where
   sieving takes hours/days, or where the count is produced by an
   expensive remote computation) verification stays at milliseconds.
2. **This does NOT verify π(x) itself yet.** π(x) = Φ(x, P_√x) + |P_√x| − 1
   needs the full wheel P = {odd p ≤ √x}, and the verifier's final
   check then costs `Õ(Σ_{p≤√x} p) = Õ(x/log x)` — worse than just
   computing π(x). The polylog-verifier regime of THIS protocol is
   exactly the polylog-wheel partial-sieve counts Φ(x, polylog-wheel).
   Closing the gap from polylog-wheel Φ to full π(x) — via GKR-style
   layered circuits over the product tree, or a native protocol for the
   Meissel-Lehmer / Lucy_Hedgehog DP — is the content of Thread 12
   slots 2–3 (see `.commit_state`).
3. **Soundness is unconditional**, unlike SNARK/STARK pipelines
   (computational soundness + Fiat-Shamir heuristic). This matters for
   the project's standards: the verified count is a theorem modulo a
   4.7×10⁻⁸ random-challenge failure probability, not modulo a
   cryptographic assumption. (Production: q ~ 2¹²⁸ extension field.)

## What would falsify this

- An accepted honest run whose output differs from the direct sieve
  count (completeness bug) — tested, none.
- A cheating-prover acceptance rate materially above n·k/q — tested at
  100 trials, 0 acceptances.
- Verifier cost empirically growing with 2^n rather than n — refuted by
  the scaling table.

## Edges / closures cited

- **E2.1** (rank(π_N) = 2^{N/2−1}+2; MLE density of χ_P): the
  obstruction that kills sum-check for *computing* π(x); this
  construction routes around it by changing the predicate class
  (automaton sieve predicates), not by contradicting it.
- **CLOSED_PATHS S18 row** ("Sumcheck / arithmetization of pi(x)",
  FAIL/E) and **S48 row** ("Multilinear extension of χ_P…", FAIL/I):
  both closures attack the truth-table-MLE layer with the *computation*
  goal. Scope correction filed at S491: they do not close the
  *verification* frame (prover allowed Θ̃(x) work), which is this
  thread.
- **S224 Correlation Dichotomy / Thread 5**: the partial-positive
  prototype shape this thread follows — no speedup for computing π(x),
  but a measurable capability on an adjacent problem (here: trust at
  polylog cost instead of recomputation cost).
- **S48 wheel-factorization row** (CRT rank-1 wheel indicator): the
  same wheel structure that gave only constant-factor *computation*
  speedup gives an exponential *verification* separation here.

## Cross-domain sources (verifiable computing / interactive proofs)

- Lund, Fortnow, Karloff, Nisan 1992 — sum-check protocol (the
  #P-counting IP).
- Goldwasser, Kalai, Rothblum 2008 — doubly-efficient delegation (GKR).
- Cormode, Mitzenmacher, Thaler 2012 — practical sum-check provers.
- Thaler, *Proofs, Arguments, and Zero-Knowledge* (survey book).
- Surveys located during S491 search: arxiv.org/pdf/2308.15191 (verified
  computation SoTA), eprint.iacr.org/2025/008 (low-degree-polynomial
  interactive verifiable computing survey). No prior art found on
  succinct verification of prime-counting specifically (consistent with
  the S48 literature scan).
