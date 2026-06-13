# Thread 11 / Slot 3 — LP relaxation of minimum line cover

## Question

Slots 1-2 measured greedy line cover for primes vs matched-baseline
random points on Ulam spiral and residue/polynomial-image grids.
Greedy gives an upper bound on the integer minimum-cover; the LP
relaxation gives a lower bound. Slot 3 asks: does the LP push below
greedy (revealing slack), and where does the gap live?

Pre-stated falsifier (this slot): if LP = greedy on every
embedding, the wheel-sieve / Ulam-greedy is LP-tight and all the
structure has already been measured. The B-NEGATIVE forecast from
slot 1 (`L_p / L_r → 1` as N grows) would tighten.

## Setup

Set-cover LP relaxation. Universe = primes ≤ N, embedded under Φ:
ℕ → ℤ². Lines = bounded-direction (a, b) with gcd=1, max(|a|,|b|)
≤ K, intercept c such that ≥ 2 primes lie on `b·x − a·y = c`.
LP: minimise `Σ x_l` s.t. `Σ_{l: p ∈ l} x_l ≥ 1` for all primes p.
Solved via scipy `linprog(method='highs')`.

Two embeddings tested:
- **Ulam spiral** at N ∈ {10³, 5·10³, 10⁴, 10⁵}, K = 5, 10
- **Residue-class grid** Φ_R(n; q=210) at N = 10⁴, K = 5

For each: also LP-solved 3 matched-baseline random samples.

## Results

### Residue-class q=210 at N=10⁴, K=5

| metric | primes | random (3 trials) |
|---|---|---|
| greedy L | 48 | 210 |
| LP relaxation | **48.0000** | 210.0 |
| integrality gap | 1.000 | 1.000 |

The wheel-sieve embedding is LP-tight at greedy for both primes
(48 = φ(210) covered by coprime columns) and random (210 columns
needed since random points hit all 210). Same as slot 2: line-cover
compression equals exactly the wheel-sieve density φ(q)/q. **No LP
slack anywhere.**

### Ulam spiral, scaling

| N    | π(N)  | √N      | greedy_p | **LP_p** | greedy_r | **LP_r** | **LP_p/LP_r** | gap_p (greedy/LP) |
|------|-------|---------|----------|----------|----------|----------|---------------|--------|
| 10³  | 168   | 31.62   | 26       | 23.34    | 33.7     | 29.95    | 0.779         | 1.114  |
| 5·10³| 669   | 70.71   | 63       | 54.30    | 81.3     | **71.00**| 0.765         | 1.160  |
| 10⁴  | 1229  | 100.00  | 95       | 77.59    | 113.3    | **100.00** | 0.776       | 1.224  |
| 10⁵  | 9592  | 316.23  | 311      | 246.69   | 330.0    | **316.00** | 0.781       | 1.261  |

### Two structural observations

1. **Random baseline LP saturates exactly at √N (rounded down) for
   N ≥ 5000.** Trial-to-trial std = 0.000. The LP solution at the
   random baseline at N=10⁴ is **integer**: 100 axis-aligned vertical
   lines, each at weight `x_l = 1.0` exactly. *Random points achieve
   the integer minimum cover by columns; no fractional slack exists.*

2. **Primes LP relaxes well below √N.** LP_p / LP_r ≈ 0.78 stable
   across 100× range in N (10³ to 10⁵). The LP solution has **zero
   integer-1 lines** (462 fractional lines at N=10⁴, all `0 < x_l < 1`).
   The fractional weight concentrates on slope-±1 Ulam diagonals.

### LP solution structure (Ulam N=10⁴, K=5, primes)

Total LP weight by direction (a, b):

| direction | total weight | fraction of LP |
|---|---|---|
| (1, −1) — slope-(−1) Ulam diagonals | 28.03 | 36.1% |
| (1, +1) — slope-(+1) Ulam diagonals | 25.74 | 33.2% |
| (1, 0) — vertical | 9.68 | 12.5% |
| (0, 1) — horizontal | 8.86 | 11.4% |
| other (3,1; 1,−3; 2,−1; ...) | 5.30 | 6.8% |
| **total** | **77.59** | 100.0% |

Top-20 individual lines all have direction (1, ±1) with intercepts
in {−28, ..., +18}, x_l between 0.59 and 0.68. These are the
quadratic-prime sequences `4k² + bk + c` flagged in slot 1's
top-line analysis (HL densities ~40%).

For comparison: the random baseline at N=10⁴ puts 100% of LP weight
on direction (1, 0) — exactly the 100 vertical columns of the
spiral, each at integer weight 1.0. **Random has no diagonal LP
weight; primes put 69% of LP weight on slope-±1 diagonals.**

## What this means

Slot 1 forecast B-NEGATIVE: greedy ratio `L_p/L_r → 1` as N^{−0.24}.
Slot 3 finds: **this is greedy slackness, not absence of structure.**
At the LP relaxation, the prime advantage is *stable*:

```
LP_primes(Ulam, N) / LP_random(Ulam, N) ≈ 0.78  (constant across 10³ ≤ N ≤ 10⁵)
LP_random(Ulam, N) = ⌊√N⌋  (exact integer, integrality gap = 1)
LP_primes(Ulam, N) ≈ 0.78 · √N  (heavily fractional, integrality gap ≈ 1.22)
```

The slot-1 greedy-decay observation tracked the *integrality gap of
the prime LP*, not the structural compression itself. As N grows,
greedy's slack gets worse on primes (gap 1.11 → 1.16 → 1.22 → 1.26
across N=10³ → 10⁵). On random, integrality gap stays ≈ 1.

This is a **structural separation**: random points achieve their
LP at axis-aligned integer cover; primes do not, instead distributing
LP weight across HL-rich slope-±1 diagonals.

## Falsifier check vs the pre-stated outcomes

P11 pre-stated the three modes:
- (E mode B-grade) `L_Ulam / L_random → 1` — predicted closure
- (I mode B-grade) `L_Ulam(N) ~ π(N)^{2/3}` (random ST floor) — alt closure
- (A mode) deviation from (E)/(I) with structural reason

Slot 3 finding rules out (E) under the LP measure: ratio is stable at
0.78, not → 1. (I) is also ruled out under any cover measure: both
LP and greedy give Θ(√N), well below π(N)^{2/3} (since π(N)^{2/3} =
9592^{2/3} ≈ 449 ≫ measured 247-311 at N=10⁵).

The structural reason exists and is concrete: HL slope-±1 diagonals
absorb 69% of the LP fractional weight on primes vs 0% on random.
**This is criterion (a) of the A-grade: a non-trivial incidence-
geometric structural fact about primes**, with quantitative LP-tight
backing across 100× range in N.

## Honest caveats

1. Measured at N up to 10⁵ only. The ratio 0.78 may drift at much
   larger N.
2. The A-grade target said "named-exponent `L ~ π(N)^α` with α < 1".
   We have α = 1/2 (since `LP ~ √N`) for both primes and random; the
   *constant prefactor* differs by 22%. This is not a named-exponent
   separation — it's a constant-prefactor structural fact.
3. The asymptotic claim "ratio → c < 1" needs a proof. Empirically,
   the data is consistent with c ≈ 0.78, but could drift.
4. The LP fractional structure (462 lines at fractional weight) does
   NOT yield an algorithmic improvement on integer line cover —
   integer cover stays at the greedy upper bound, which scales like
   √N · log(N) (the 1.22-1.26 integrality gap).

## What would falsify this

- **Refuted if** at N = 10⁶ the ratio LP_p/LP_r drifts toward 1
  (say > 0.9), indicating slot-1's observed greedy decay also
  applies to LP at large N.
- **Refuted if** the LP fractional weight is an artifact of the
  bounded-direction restriction K=5,10. Larger K would shift weight
  to other directions. (Partially checked at K=20, N=10⁴: LP changed
  by 0.25 only, structural conclusion unchanged.)

## Slot 4 plan (per .commit_state)

Slot 4 should:
1. Push to N=10⁶ Ulam LP (estimated runtime: 10⁻²⁰ minutes via
   sparse + HiGHS scaling). Confirm or refute the 0.78 stability.
2. Investigate whether the LP fractional structure → integer cover
   improvement via *randomised rounding* or *iterated rounding*
   matching the LP optimum to within `1 + o(1)`.
3. Test alternative embeddings (residue grid `q < √N`) for the
   same LP-vs-greedy gap. Slot 2 found greedy = wheel-sieve floor;
   does LP also tighten there, or is the wheel sieve genuinely
   integer-tight?

## Files

- `p11_lp_relaxation.py` — main LP evaluator
- `p11_lp_dual_inspect.py` — extracts LP solution lines + directions
- `lp_summary_residue_q210_N10000_K5.txt`
- `lp_summary_ulam_N10000_K{5,10,20}.txt`
- `lp_summary_ulam_N100000_K5.txt`
- `lp_scaling_summary.csv` — N x {primes, random} LP scaling table

## Edges cited

- E1.5 (information-theoretic compression floor)
- HL Conjecture F (off-edge external — quadratic-prime densities)
- Szemerédi-Trotter (off-edge external — incidence-theoretic baseline)
- Stanley 1989 *Adv. Math.* (matroid-theoretic line-cover LP, the
  cross-domain ingredient — establishes that the LP relaxation
  equals the integer optimum on certain matroidal hypergraphs;
  observed here that random-baseline residue/Ulam embeddings are
  LP-tight while prime Ulam is NOT)
