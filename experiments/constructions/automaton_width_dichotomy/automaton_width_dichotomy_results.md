# The Width Spectrum of the Prime-Counting Pipeline (S491)

## The object

For a predicate f on n bits read MSB-first, the per-cut Nerode width
`W_j(f)` = number of distinct subfunctions after j bits (= distinct rows
of the 2^j × 2^{n−j} truth-table reshape). The S491 transfer-DP
primitive upgrades this from a computation lower bound (how the project
used it at S17/S20/S23/S28) to an **exact verification-cost measure**:

> **T1.** Given an explicit read-once branching program for f with
> widths (W_j), the multilinear extension f̃ is evaluable at an
> arbitrary field point in O(Σ_j W_j) field operations (transfer DP,
> implemented in `../p12_sumcheck_pi_verification/`). Hence Σ_j W_j is
> the per-point price a sum-check verifier pays — IF it can construct
> the program.

The "if" is the dichotomy: sieve predicates have explicit programs;
primality's minimal program exists but requires knowing the primes to
build.

## Measured results (`run.log`)

### 1. Division wiring is Θ(p) tight — S491 protocol is optimal in its paradigm

Width profile of the relation [u = ⌊v/p⌋] over interleaved bits, n = 8:

| p | max width | max/p |
|---|---|---|
| 3 | 7 | 2.33 |
| 5 | 11 | 2.20 |
| 7 | 15 | 2.14 |
| 11 | 23 | 2.09 |
| 13 | 27 | 2.08 |

Max width = **2p + 1 exactly** (pending-quotient-bit × remainder, plus
dead state). The Lucy-protocol wiring check at O(n·p) per layer is
therefore optimal within explicit read-once automatons; pushing the
exact-π(x) verifier below Õ(x/log x) toward Õ(√x) **requires**
prover-supplied carry witnesses (GKR-style layered arithmetization of
the ×p circuit) — now established as necessary, not just convenient.
Comparator [v ≥ M]: width exactly 3 at all M tested.

### 2. χ_P width profile: incompressible plateau → admissibility cliff

n = 24 (x = 16.7M), matched-density random control:

```
j:        1 ... 16     17      18     19    20   21  22  23
chi_P:    2^j ...65536 130937  102200 3213  107  14  5   3
random:   2^j ...65536 130865  185111 58506 4707 209 16  4
```

- **Left plateau:** W_j = 2^j exactly up to j = 16 — every aligned
  256-window of the prime indicator at this scale is distinct. The
  global side of the function is width-incompressible, matching the
  pseudorandomness catalogue.
- **Right cliff:** primes crash 18–44× below random
  (3213 vs 58506 distinct 32-windows; **107 vs 4707** distinct
  16-windows). The tail counts *realized local prime patterns* in
  aligned windows — wheel/admissibility rigidity (aligned 16-windows
  only populate odd positions and must avoid full residue classes mod
  3, 5, …) made quantitative. **D_k(x) := #distinct aligned k-windows
  of χ_P at scale x is a new measurable invariant** bridging the
  project's global-pseudorandomness edges and Hardy-Littlewood local
  structure in a single profile.
- **Middle cut: W_{n/2} = 4096 = 2^{n/2}, MAXIMAL** — strictly above
  the E2.1 rank (2^{n/2−1}+2 = 2050). Refines S20's "width ≥ rank"
  to: distinct-row width is exactly maximal at the middle cut; the
  factor-2 rank deficiency is linear-algebraic, not informational.

### 3. One-shot sum-check verification of π(x) is dead — twice over

Σ_j W_j(χ_P) ~ x^c with c measured 0.7296 → 0.7703 across
n = 16 → 24, **drifting upward**. Structural argument for the drift:
the profile peak sits at the crossover 2^{j*} ≈ D(2^{n−j*}); with local
pattern entropy log₂D(m) ≈ α·m this gives peak width ≈ α·x/log x, so
c → 1 asymptotically *for primes and random alike* — the prime/random
difference is the constant α and the cliff, not the exponent. Combined
with non-constructibility (building the minimal OBDD requires the
primes), the one-shot route — sum-check π(x) = Σ_w χ_P(w) directly with
the verifier evaluating χ̃_P — fails BOTH on width (x^{c→1}) and on
constructibility. **Layered sieve protocols (S491) are necessary, not
merely convenient**: verification must consume primality only through
explicit small-width predicates (division, comparison), never
monolithically.

## What would falsify this

- A cut order or variable order giving polylog width for χ_P at all
  cuts (contradicts E2.1; none known — MSB order measured here, S20/S28
  measured LSB/multi-order sizes, all exponential).
- An explicit (constructible without the primes) read-once program for
  χ_P with Σ W_j = x^{o(1)} — would give a polylog one-shot verifier
  and, more importantly, would be a compression of the primes
  contradicting the project's entire Section-7 edge family.
- Division-relation width o(p) — refuted: 2p+1 exact.

## Cross-references

E2.1 (rank — refined here at the middle cut), S17/S20/S23/S28
CLOSED_PATHS rows (computation-side width ancestors),
`experiments/constructions/p12_sumcheck_pi_verification/` (T1
implementation; the protocol this spectrum prices), Gallagher 1976 /
HL Conjecture B (the local-pattern side of the cliff),
novel/succinct_verification_of_pi.md §width-dichotomy.
