# The Local Pattern Census of Primes: the Width Cliff Is Exactly Hardy-Littlewood (S491)

## Headline

The width profile W_j of the prime indicator at scale 2ⁿ (measured in
`experiments/constructions/automaton_width_dichotomy/`) is, *by
definition*, the census curve W_j = D_k(2ⁿ) with k = 2^{n−j}, where
D_k(x) = #distinct exact occupancy patterns of aligned k-windows of
primes ≤ x. This experiment characterises that census exactly:

> **D_k(x) = A_k + 1 for all x beyond a saturation threshold x₀(k)**,
> where A_k is the exactly-enumerable count of Hardy-Littlewood
> admissible aligned patterns and the single "+1" is the m = 0 window
> containing the small primes themselves (the only inadmissible pattern
> that ever occurs).

So the two-phase structure of the prime width spectrum is now fully
explained: **plateau = census unsaturated (all windows distinct,
global pseudorandomness); cliff = census saturated at the HL admissible
count (local rigidity)**. The prime indicator's read-once width profile
*is* the admissible-pattern entropy curve.

## Admissibility criterion (exact, aligned windows, proof = divisibility)

S ⊆ {0..k−1} is admissible iff (i) for every prime q | k, no o ∈ S has
q | o (else k·m + o ≡ 0 mod q in *every* window); (ii) for every prime
q ≤ k, q ∤ k, the residues {o mod q : o ∈ S} do not cover Z_q (else
every window has a q-divisible occupied offset). Patterns violating
(i)/(ii) can only occur in windows where the obstructing element *is*
the prime q — i.e., finitely many small-m windows. Under the prime
k-tuple conjecture every admissible S occurs as an *exact* pattern
(prime on S, composite off S) infinitely often.

## Measured (runs `run_k8.log`, `run_k16.log`, `run_k32.log`)

| k | A_k (exact) | weights | D_k saturates at | value | exceptional |
|---|---|---|---|---|---|
| 8 | 13 | ≤ 3 | < 2¹⁶ | **14 = 13+1** | [0,8): {2,3,5,7} |
| 16 | 106 | ≤ 5 | 2¹⁸ | **107 = 106+1** | [0,16): {2,...,13} |
| 32 | 3573 | ≤ 9 | > 2²⁶ (converging) | 3385 at 2²⁶ | [0,32): {2,...,31} |

- **k = 16: all 106 admissible patterns realized by x = 262144**,
  including all six weight-5 constellations; D₁₆ then constant through
  2²⁸. The width-profile value 107 at (n=24, j=20) is exactly A₁₆+1.
- **k = 32 convergence curve:** D = 1093 → 1825 → 2452 → 2940 → 3213 →
  3385 across x = 2¹⁶..2²⁶ (the 2²⁴ value matches the width profile's
  3213 exactly). **189 admissible patterns remain unrealized at
  x = 2²⁶ — every one of weight ≥ 6** (23 of weight 6, 115 of 7, 47 of
  8, 4 of 9). The census deficit is precisely the
  first-occurrence-of-dense-constellations tail, HL-predicted to decay
  like the rarest pattern's density ~ S_pattern·x/log^w x. D_k(x)−A_k
  is thus a direct finite-x probe of k-tuple uniformity.
- **External cross-validation:** maximal realized/admissible weights
  3, 5, 9 at k = 8, 16, 32 agree with the classical narrowest-
  admissible-tuple function (H(9) = 30 ≤ 31 < H(10) = 32 — a 10-tuple
  cannot fit in a 32-window). The four admissible weight-9 patterns are
  the classical densest 9-tuples translated to odd offsets.
- **Admissible entropy law (S491-late, `--entropy-scan`, `run_entropy.log`):**
  the brute-force enumeration was replaced by a pruned DFS exploiting
  downward-closure of admissibility (a subset of an admissible set is
  admissible ⟹ a partial pattern covering all residues mod any prime
  kills its supertree); node count scales with A_k itself, reaching
  k = 80 in seconds. Exact counts: A₈=13, A₁₂=16, A₁₆=106, A₂₀=121,
  A₂₄=227, A₂₈=640, A₃₂=3573, A₃₆=2704, A₄₀=5704, A₄₄=19825,
  A₄₈=29002, A₅₂=87438, A₅₆=93751, A₆₀=53602, A₆₄=1581920,
  A₈₀=5777381. Two structural facts:
  1. **A_k is NOT monotone in k** — the dips (A₃₆ < A₃₂, A₆₀ ≪ A₅₆)
     are the divisor structure of k: each q | k removes the 1/q
     fraction of offsets divisible by q (condition (i)), so highly
     composite k have far fewer allowed offsets (k = 60: only
     φ-fraction 4/15). Cross-check: A₈₀ (allowed = 32 offsets, same as
     k = 64, but weaker mod-q binding) = 5.78M > A₆₄ = 1.58M. ✓
  2. **On the clean family k = 2^m** (only the parity divisor):
     ln A_k / (k/ln k) = 0.667, 0.808, 0.886, 0.928 at k = 8, 16, 32,
     64 — increments halving at each doubling (geometric convergence),
     extrapolating to **1.0 ± 0.03**. Conjecture (falsifiable at
     k = 128 with an analytic/transfer-matrix count — DFS is
     enumeration-bound there): **ln A_k ~ k/ln k**, i.e.
     **A_k = e^{(1+o(1))·π(k)} — the admissible-pattern entropy of a
     window equals its prime count.** Lower bound is rigorous:
     A_k ≥ 2^{ρ*(k)} = e^{Θ(k/ln k)} (subsets of a maximal admissible
     tuple, Hensley-Richards); the binomial upper bound gives
     e^{O(k·ln ln k/ln k)} — the data says the lower bound's exponent
     order is the truth with constant 1. This constant sets the
     plateau→cliff corner of the χ_P width profile asymptotically.

## Consequences for the verification thread

The cliff height being A_k + 1 (not polylog) re-confirms quantitatively
that no cut of χ_P is narrow enough for one-shot sum-check
verification: even the maximally-compressed right tail is exponential
in the window size, and the left plateau is maximal. Verification must
(and does — the S491 Lucy protocol) route around χ_P entirely.

## What would falsify

- A second inadmissible pattern at m > 0 (impossible by divisibility —
  bug check); none found through 2²⁸.
- An admissible pattern of weight ≤ 5 (k=32) missing at 2²⁶ while
  higher-weight cousins appear — would contradict HL density ordering;
  not observed (missing set is exactly the weight ≥ 6 stratum).
- D₁₆ rising above 107 at larger x — impossible unless the census code
  or the admissibility proof is wrong.

## Files

`local_pattern_census.py` (one script, `--k`, `--nmax`),
`run_k8.log`, `run_k16.log`, `run_k32.log`.

## Cross-references

`experiments/constructions/automaton_width_dichotomy/` (the width
profile this closes), E2.1 (middle-cut maximality — the plateau side),
novel/pseudorandomness_of_pi.md (global-randomness edges — the plateau
is their width-spectrum shadow), Hardy-Littlewood Conjecture B /
Gallagher 1976 (the cliff side), H(k) narrowest-tuple tables
(external validation), novel/succinct_verification_of_pi.md
§width-dichotomy.
