# Thread 11 / Slot 2 — alternative 2D embeddings

## Goal (slot 2)

Slot 1 (S484) measured `L_Ulam(primes, N)` on the Ulam spiral and
found both primes and matched-baseline random points scale as ~√N
to leading order, with the prime advantage `L_p/L_r` shrinking as
`N^{-0.24}` (B-NEGATIVE forecast for Ulam alone).

Slot 2 asks whether an *alternative* 2D embedding gives qualitatively
different scaling. Two embedding families are tested:

- **Residue-class grid** `Φ_R(n; q) = (n mod q, ⌊n/q⌋)`
  — grid is `q` columns × `N/q` rows.
- **Polynomial-image grid** `Φ_Q(n; q) = (n, n² mod q)`
  — strip is `N` wide × `q` tall.

For each, the bounded-direction greedy line cover is computed at
`N ∈ {10⁴, 10⁵, 10⁶}`, varying the modulus `q`.

## What's new in this session

- `p11_alt_embeddings.py` — generalises the slot-1 evaluator to
  accept different embeddings. Imports the bounded-greedy machinery
  from `p11_ulam_bounded.py`; supports `--p-list` (any modulus list)
  and toggles between residue / polynomial / both.
- Two embedding families × two modulus regimes (prime moduli
  `p ∈ {2,3,5,7,11}` and primorial moduli `q ∈ {6, 30, 210, 2310}`)
  × three N values × matched-baseline trials.

## Results

### Prime modulus p ∈ {2, 3, 5, 7, 11}

| q | N=10⁴ | N=10⁵ | N=10⁶ |
|---|---|---|---|
| residue_p=2 | L_p=2, L_r=2 (1.00) | L_p=2, L_r=2 (1.00) | L_p=2, L_r=2 (1.00) |
| residue_p=3 | L_p=3, L_r=3 (1.00) | L_p=3, L_r=3 (1.00) | L_p=3, L_r=3 (1.00) |
| residue_p=5 | L_p=5, L_r=5 (1.00) | L_p=5, L_r=5 (1.00) | L_p=5, L_r=5 (1.00) |
| residue_p=7 | L_p=7, L_r=7 (1.00) | L_p=7, L_r=7 (1.00) | L_p=7, L_r=7 (1.00) |
| residue_p=11 | L_p=11, L_r=11 (1.00) | L_p=11, L_r=11 (1.00) | L_p=11, L_r=11 (1.00) |
| poly_p=2 | L_p=2, L_r=2 (1.00) | L_p=2, L_r=2 (1.00) | L_p=2, L_r=2 (1.00) |
| poly_p=3 | L_p=2, L_r=2 (1.00) | L_p=2, L_r=2 (1.00) | L_p=2, L_r=2 (1.00) |
| poly_p=5 | L_p=3, L_r=3 (1.00) | L_p=3, L_r=3 (1.00) | L_p=3, L_r=3 (1.00) |
| poly_p=7 | L_p=4, L_r=4 (1.00) | L_p=4, L_r=4 (1.00) | L_p=4, L_r=4 (1.00) |
| poly_p=11 | L_p=6, L_r=6 (1.00) | L_p=6, L_r=6 (1.00) | L_p=6, L_r=6 (1.00) |

**Observation 1.** For prime modulus p, both embeddings give a CONSTANT
line cover (independent of N) for both primes and matched-random.
- Residue-class: `L_R(p) = p`. Reason: vertical cover by p columns is
  optimal — each column = an arithmetic progression mod p — and
  primes occupy all p classes (one is `{p}` itself, the others are
  φ(p) = p−1 coprime classes). Random samples also occupy all p.
- Polynomial-image: `L_Q(p) = (p+1)/2` for p odd, 2 for p=2. Reason:
  `n² mod p` takes (p+1)/2 distinct values (the QR set ∪ {0}); both
  primes and random fill all of them.

**Compression ratio is exactly 1.00 in all 30 prime-modulus cells.**
For prime modulus, primes provide no advantage over matched-random:
the embedding's quotient structure is the only constraint, and primes
saturate it just as random samples do.

### Primorial modulus q ∈ {6, 30, 210, 2310}

| q | embedding | N=10⁵ L_p | L_r | ratio | N=10⁶ L_p | L_r | ratio |
|---|---|---|---|---|---|---|---|
| 6=2·3 | residue | 3 | 6 | 0.500 | 3 | 6 | 0.500 |
| 6 | poly | 2 | 4 | 0.500 | 2 | 4 | 0.500 |
| 30=2·3·5 | residue | 9 | 30 | 0.300 | 9 | 30 | 0.300 |
| 30 | poly | 4 | 12 | 0.333 | 4 | 12 | 0.333 |
| 210=2·3·5·7 | residue | 49 | 210 | 0.233 | 49 | 210 | 0.233 |
| 210 | poly | 9 | 48 | 0.188 | 9 | 48 | 0.188 |
| 2310=2·3·5·7·11 | residue | 44 | 44 | 1.000 | 433 | 433 | 1.000 |
| 2310 | poly | 34 | 288 | 0.118 | 34 | 288 | 0.118 |

**Observation 2 (column-dominated regime, q² < N).** When `q² < N`
(rows = `⌊N/q⌋ > q` = columns), vertical-line cover wins for the
greedy. Prime advantage is rigid:

- residue_q : `L_p ≈ φ(q) + ω(q)`, `L_r = q` ⇒ ratio ≈ `φ(q)/q`
  (Mertens-like: 0.500, 0.300, 0.233 for q = 6, 30, 210).
- poly_q : `L_p ≈ #{n² mod q : n coprime to q} + ω(q)` (small primes
  p|q each give unique non-coprime y), `L_r = #{n² mod q : n ∈ ℤ}`.
  Asymptotic ratio = `(1/2) · ∏_{p odd, p|q} (p−1)/(p+1)`. Predicted:
  q=6: 0.5·0.5 = 0.250 (data 0.500, dominated by ω(q) correction);
  q=30: 0.5·0.5·0.667 = 0.167 (data 0.333);
  q=210: 0.5·0.5·0.667·0.75 = 0.125 (data 0.188);
  q=2310: 0.5·0.5·0.667·0.75·0.833 = 0.104 (data 0.118).
  The asymptotic prediction is the limit at q → ∞; the data is +
  the ω(q)/L_r additive boost from non-coprime small-prime images.

**Observation 3 (row-dominated regime, q² > N).** When `q > √N`
(columns > rows), the greedy switches to horizontal lines: one line
per row. Both primes and random need ⌊N/q⌋ + 1 horizontal lines —
the embedding's quotient structure is no longer used. Compression
collapses to ratio 1.

This is exactly what residue_q=2310 shows: at N=10⁵ rows = 44 ≈ q
boundary, ratio = 1.0. At N=10⁶ rows = 433, still < q = 2310, and
again `L_p = L_r = 433` (vertical lines `φ(q)+ω(q) = 485` would cost
more than 433 horizontals, so greedy picks horizontals).

The polynomial-image embedding has q-rows × N-columns, so its rows
never run out — `q < N` always — and the column-direction (slope-0,
horizontal lines) compression persists across q. Hence poly_q=2310
at N=10⁵ retains ratio 0.118 even though residue_q=2310 has lost
compression.

### Verifying wheel-sieve decomposition (top-line CSV, N=10⁵)

```
residue_q=30 top lines: all (a=0, b=1) [vertical]
  intercepts: {1, 7, 11, 13, 17, 19, 23, 29}  ← exactly φ(30)=8 coprime classes
  + (a=1, b=0, c=0) horizontal y=0 picking up {2, 3, 5}        ← ω(30)=3 small primes
  L_p = 8 + 1 = 9 ✓

residue_q=210 top lines: all (a=0, b=1) [vertical]
  top-10 intercepts: {73, 197, 107, 179, 187, 17, 53, 139, 11, 19}
  ← all coprime to 210 (each is in φ(210)=48 coprime residues)
  L_p = 48 + 1 = 49 ✓

poly_q=30 top lines: all (a=1, b=0) [horizontal]
  intercepts c = -y where y ∈ {1, 19} for coprime primes
  ← only 2 distinct n² mod 30 values for n coprime to 30  
  L_p = 2 + 2 = 4 ✓
```

### Closed-form predictions for L_Φ at primorial q

Let `q = ∏_{p ∈ S} p` for a set S of distinct primes. Define
- `φ(q) = ∏_{p∈S}(p−1)` (Euler totient)
- `ω(q) = |S|`
- `Q(q) = ∏_{p∈S, p>2} ((p−1)/2 + 1) · 2` if 2 ∈ S else `∏ ((p−1)/2 + 1)`
  — counting all squares mod q via CRT (numerator of poly L_r)

Then for N >> q²:
```
L_R(primes, N; q) = φ(q) + ω(q) − k    [k = singleton-merge savings, small]
L_R(random, N; q) = q                  [random hits all q residues]
L_R-ratio          = (φ(q)+ω(q))/q  ≈ φ(q)/q ≈ ∏ (1 − 1/p)
```

For N ≪ q² (row-dominated):
```
L_R(primes, N; q) = ⌊N/q⌋ + 1   (horizontal cover)
L_R(random, N; q) = ⌊N/q⌋ + 1
L_R-ratio          = 1
```

### What this means for the thread

The slot 2 measurement gives a **complete map of where line-cover
compression exists** under the natural arithmetic 2D embeddings:

1. **Prime modulus p**: no compression (rigid `L_p = L_r`).
2. **Composite primorial q with q² < N**: constant-factor compression
   ratio ≈ φ(q)/q (the Mertens factor). For optimal `q ≈ √N`,
   `ratio ~ 2e^{-γ}/log(N)` (slow log-decay).
3. **Composite primorial q with q² > N**: compression collapses
   (row-dominated, embedding-rigid).
4. **Polynomial-image at composite q**: best ratio observed,
   `~ ∏_{p|q, p>2} (p−1)/(p+1)` (also slow log-decay).

**The compression is real but is exactly the classical wheel sieve
encoded geometrically.** No incidence-geometric compression beyond
wheel-sieve / quadratic-residue structure was observed across this
slot's embeddings.

This is consistent with the slot-1 Ulam result: under any embedding,
prime line-cover compression is bounded by classical density factors
(HL for Ulam diagonals; wheel-sieve for residue-grids; QR-density
for polynomial-image). No "new" structure appears.

## What would falsify this

- Find an embedding where `L_p/L_r → 0` faster than `1/log(q)` (i.e.,
  beats the Mertens floor). None known; it would require structural
  prime-clustering that wheel-sieve does not predict.
- An LP relaxation or branch-and-bound improvement showing the greedy
  is loose by more than the current `ω(q)`-merge savings (~ 1-3
  fewer lines).

## Edges cited / composed

- E1.5 (information-theoretic bound — primes' bit-content matches
  random subject to wheel-sieve density correction). Confirmed:
  geometric line-cover compression matches the wheel-sieve ratio,
  no further compression observable.
- HL Conjecture F (off-edge external) — irrelevant to non-Ulam
  embeddings tested here; the structure observed here is residue-
  class only, no quadratic-form coupling.
- Szemerédi-Trotter (1983) — random-points line-cover lower bound,
  used as theoretical floor for `L_random`.
- Wheel-sieve density — classical (Mertens 1874): `∏_{p≤y} (1−1/p)
  ~ e^{-γ}/log(y)`. Geometrically realised here as line-cover ratio.

## Files

- `p11_alt_embeddings.py` — multi-embedding evaluator (this slot)
- `summary_alt_N{10000,100000,1000000}_K5.csv` — scaling data
- `top_lines_alt_N{10000,100000,1000000}_K5.csv` — dominant lines

## Slot 2 verdict

**B-grade.** Substantive empirical map across 5 prime moduli + 4
primorial moduli × 2 embeddings × 3 N values. Found:
- Prime modulus: no compression (`L_p/L_r = 1.000` in all 30 cells).
- Primorial modulus: classical wheel-sieve compression, ratio
  `~ φ(q)/q` (residue) or `~ ∏ (p−1)/(p+1)` (polynomial-image),
  bounded away from 1 only when `q² ≪ N`.
- Optimal q for compression: `q ≈ √N`. At this regime, ratio
  decays as `1/log(N)` — slow but real.

**No A-grade.** The compression is the wheel sieve; no
embedding-specific incidence-geometric structure beyond what
classical density theory predicts. The result extends slot 1's
B-NEGATIVE forecast: line-cover compression in 2D embeddings is
universally bounded by classical density factors.

**Next action (slot 3):** theoretical shape. Either prove a
Szemerédi-Trotter-style lower bound `L_Φ(N) ≥ c · π(N)^{2/3}` in
the column-dominated regime (matching the random points lower bound,
extending to all "narrow" embeddings), OR investigate whether an
LP-relaxation gives a strict improvement over greedy.

## Self-evaluation

1. **What did I produce that was not in the project before this
   session?**
   - First multi-embedding line-cover map. Established that the
     compression exists *only* under composite (primorial) modulus
     in the column-dominated regime, with ratio bounded by Mertens
     factor.
   - Confirmed top-line decomposition: residue-class lines are
     coprime residues; polynomial-image lines are coprime QR
     residues. Pure wheel-sieve structure.

2. **What edges did my work compose or cite?** E1.5; Szemerédi-Trotter
   (off-edge external); Mertens 1874 wheel-sieve.

3. **If duplicate closures, why?** N/A — slot 2 is empirical setup,
   not closure.

4. **Next action:** see Slot 2 verdict above.
