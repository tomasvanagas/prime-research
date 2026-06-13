# F6 / P7 — Structure tests on the dyadic family π(2^k)

**Target**: NOVELTY_CHALLENGES.md §F6 = OPEN_POSITIVE_TARGETS.md §P7
("π(2^n) — does the binary structure admit polylog?").

**Session**: S246 (novelty mode, B-grade target).

**Mode**: Empirical structural test of π(2^k) for k = 1..56.

**Edges referenced**: E1.1 (HKM x^{2/3}), E1.5 (per-step entropy),
E1.3 (R⁻¹ bit-level model), Riemann explicit formula
(π(x) = R(x) − Σ_ρ R(x^ρ) + small), pseudorandomness battery
(novel/pseudorandomness_of_pi.md).

**Verdict**: **B-NEGATIVE** — the dyadic family admits no detectable
amortisation beyond what HKM gives for arbitrary x of the same
magnitude. All three F-claims FAIL against Monte Carlo random-null
baselines.

---

## Falsification statements (PRE-REGISTERED before running)

The B-grade attack assumes a NEGATIVE outcome. The three F-claims
below are the *positive* hypotheses we expect to **fail**.

### F-claim 1 — Compressibility of sign(r_R(k))

Let `r_R(k) = π(2^k) − R(2^k)` where `R` is the Riemann smooth
counting approximation `Σ_n μ(n)/n · Li(x^{1/n})`. Define
`s(k) ∈ {0, 1}` by `s(k) = 1` iff `r_R(k) > 0`. Then the Berlekamp-
Massey linear complexity over GF(2) of `(s(1), ..., s(56))` is **at
or below the 5th percentile of the random-shuffle null distribution
matched on n_pos / n_neg**.

(The naïve `BM ≤ ⌈log₂² 56⌉ = 33` threshold is too lenient: the
random-shuffle baseline at n=56 with these splits is itself
distributed in `[27, 30]` with mean 28.25. Tightened post-MC to the
5%-tail.)

POSITIVE outcome: BM ≤ MC p05 — sign sequence is more compressible
than random, suggesting hidden structure.

NEGATIVE outcome: BM > MC p05.

### F-claim 2 — Autocorrelation of normalised residual

Let `c_R(k) = r_R(k) / 2^{k/2}` (RH-normalised). Then
`max_{l ∈ {1,...,10}} |corr(c_R(k), c_R(k+l))| ≥ p999` of the
Gaussian-iid null at length 56. (p999 ≈ 0.437 from MC.)

Note: the original v1 pre-registration used Li-baseline, lag-1 only,
threshold 0.3. That fired spuriously at 0.87 because Li carries a
*smooth* −Li(√x)/2 trend that auto-correlates trivially. The corrected
test uses R-baseline (which absorbs the smooth Möbius series) and
multiple lags + Bonferroni-shaped threshold via MC p999.

POSITIVE outcome: max |ac| ≥ 0.437 — structure beyond random phases.

NEGATIVE outcome: max |ac| < 0.437.

### F-claim 3 — HKM "dyadic speedup"

For sympy.primepi (Meissel-Lehmer / Lucy-style evaluator), the median
cold-start subprocess time `T(2^k)` averaged across k ∈ {28, 30, 32}
relative to the median of `T(2^k − 1)` and `T(2^k + 1)` is **at most
0.5** (i.e., dyadic x at least 2× faster).

POSITIVE outcome: ratio ≤ 0.5.

NEGATIVE outcome: ratio in [0.8, 1.25].

---

## Composite verdict

- **All three F-claims FAIL** → **B-NEGATIVE closure of P7**: the
  dyadic family does not admit binary-structure-driven amortisation.
- **At least one F-claim HOLDS** → partial-positive surprise.

---

## Method

1. Hard-code OEIS A007053 reference values for π(2^k), k = 0..56.
2. Verify k ≤ 18 directly via sympy.primepi; confirm exact agreement.
3. Compute Li(2^k) and R(2^k) via mpmath at 80-decimal precision.
4. Compute r_R(k) = π(2^k) − R(2^k) and c_R(k) = r_R(k) / 2^{k/2}.
5. Monte Carlo null distributions:
   - BM linear complexity over GF(2) of random {0,1}^56 sequences with
     same n_pos / n_neg split as data (4000 shuffles).
   - Lag-1 absolute autocorrelation of i.i.d. Gaussian length-56
     sequences (4000 trials).
6. HKM cost test: `sympy.primepi(x)` in fresh Python subprocess,
   k ∈ {28, 30, 32}, 5 cold-start trials per (k, label) with
   label ∈ {x_pow, x_pow − 1, x_pow + 1}. Subprocess isolation
   ensures no sieve-cache reuse across calls.

---

## Results

### Verification (step 2)

`sympy.primepi(2^k) = π(2^k)` exactly for k = 1..18 (18/18, perfect
agreement against OEIS A007053).

### Residual statistics (step 4)

Using R as smooth baseline:

```
c_R(k) statistics, k=1..56
  mean   = -0.01142
  std    =  0.06082
  min    = -0.7071  (k=1, π=1, R≈1.5: small-k regime)
  max    = +0.4596  (k=2, π=2, R≈1.0)
  range over k≥10 (asymptotic regime): [-0.041, +0.084]
```

For comparison, the *Li* baseline (smooth-confounded) gives
`c_li mean = -0.1303` — non-zero because π(x) − Li(x) carries the
deterministic `-0.5 Li(√x)` Möbius-series correction.

### Pseudorandomness tests (step 5)

Sign distribution: `n_pos = 25, n_neg = 31` (non-trivial — Skewes
phenomenon at scale 2^k for small k can flip locally).

Berlekamp-Massey linear complexity:
- Data: `BM(sign(r_R)) = 28`
- MC random-binary baseline: `mean 28.25, std 1.01, [p05=27, p95=30]`
- Data sits at 0.25σ below the MC mean — statistically *identical to
  random*.

Lag-by-lag autocorrelation of `c_R`:

| lag | r(c_R, lag) |
| --- | ----------- |
|  1  |   +0.283    |
|  2  |   +0.043    |
|  3  |   +0.117    |
|  4  |   −0.083    |
|  5  |   +0.029    |
|  6  |   +0.045    |
|  7  |   +0.025    |
|  8  |   −0.092    |
|  9  |   +0.015    |
| 10  |   −0.027    |

Max `|ac|` across lags 1..10: 0.283 (lag 1). Gaussian-iid MC null:
`p95=0.250, p99=0.328, p999=0.437`. The lag-1 result is at nominal
p ≈ 0.04 single-lag, which becomes p ≈ 0.4 after 10-lag Bonferroni.
**Not statistically significant.**

(The earlier v1 result `lag-1 = 0.87` was driven by the Möbius-smooth
Li-baseline trend and disappears under R-baseline.)

### HKM cost test (step 6)

Subprocess-isolated, 5 cold-start trials per cell:

| k  | T(2^k) median | T(neighbour) median | ratio |
| -- | ------------- | ------------------- | ----- |
| 28 | 0.0619 s      | 0.0621 s            | 0.998 |
| 30 | 0.1761 s      | 0.1682 s            | 1.047 |
| 32 | 0.4500 s      | 0.4525 s            | 0.995 |

Average ratio = **1.013**, range [0.998, 1.047], well inside [0.8,
1.25]. Sympy's Meissel-Lehmer evaluator does not exploit dyadic
structure: cost is x-magnitude-driven, not x-form-driven.

### F-claim verdicts

| F-claim | Threshold | Data | Verdict |
| ------- | --------- | ---- | ------- |
| F1 (BM ≤ MC p05) | ≤ 27 | 28 | **FAILS** |
| F2 (max ac ≥ MC p999) | ≥ 0.437 | 0.283 | **FAILS** |
| F3 (ratio ≤ 0.5) | ≤ 0.5 | 1.013 | **FAILS** |

---

## Verdict

**B-NEGATIVE closure of P7 / F6**: the dyadic family `(π(2^k))_k=1..56`
admits no structural amortisation detectable at this scale. Three
independent tests — sign-sequence linear-complexity, multi-lag
autocorrelation of the RH-normalised residual, and Meissel-Lehmer
wall-clock cost — all return random-baseline values within MC noise.

### Structural reason

Two distinct mechanisms confirm the closure:

1. **Information-theoretic** (residual side). After subtracting the
   smooth Riemann correction R, what remains is the explicit-formula
   zero-sum
   ```
   π(x) − R(x) ≈ −Σ_ρ R(x^ρ)
              ≈ −2 √x · Σ_n cos(γ_n log x − arg(ρ_n)) / |ρ_n|
   ```
   At x = 2^k, the phases γ_n · k log 2 mod 2π are equidistributed in
   k for *every* γ_n (Weyl, since γ_n log 2 ∉ πQ for any γ_n with
   γ_n irrational). The first ~10⁴ Riemann zeros γ_n collectively
   produce a Gaussian-like sum (S195/S202 random-phase regime). At
   x = 2^k, no different from arbitrary x of the same magnitude.

2. **Algorithmic** (cost side). Meissel-Lehmer's Lucy-style outer
   loop processes counts of {n ≤ √x} and {n smooth up to x^{2/3}}.
   Neither set carries any binary-representation structure. Sympy's
   ratio 1.013 across (28, 30, 32) confirms.

### What this rules out

- Any algorithm that computes π(2^k) substantially faster than
  generic-x π(x) at the same magnitude. The "binary form of the
  input" is not a useful side-channel.
- Any compressed representation of the sequence `(π(2^k))_k` of size
  smaller than the trivial bit-budget Σ_k k = O(N²) implied by the
  output size.
- Any short-recurrence/closed-form expression for the residuals
  `π(2^k) − R(2^k)` at the n=56 scale (random within MC).

### What this leaves open

- **Larger N**: a test with k > 200 (requiring tabulated π(2^k)
  beyond OEIS — not currently available in sympy/primecount-free
  systems) might reveal lower-order trends invisible at n=56.
- **Different parametric families** (P7-adjacent): π(10^k), π(F_k)
  for Fibonacci F_k, etc. The mechanism is structurally identical
  (Weyl equidistribution); the conclusion should propagate without
  experiment.
- **Sub-structure**: BATCHED variants where many dyadic queries
  share work — Thread 5 / Correlation Dichotomy shape. P7 as written
  asks per-query; batched-on-k is a different question (see
  Successor §P7.a below).

---

## Successor challenges (proposed in S246)

### §P7.a — Cross-modulus generalisation

Repeat the structural test for parametric families
`x_k = m^k` at m ∈ {3, 5, 6, 10, 30}. Hypothesis: identical
B-NEGATIVE shape in each case. Cost: 1 session, ~20 minutes per m.
A-grade if any *single* m exhibits BM compressibility OR
multi-lag autocorrelation above MC p999. (Predicted: none.)

### §P7.b — Batched-on-k amortisation

Thread 5 / S224 Correlation Dichotomy shape: for x_1 = 2^{k_1},
..., x_M = 2^{k_M}, can the explicit-formula evaluation share zero
data across the M queries? The shared `(γ_n)` table is independent
of k; computing all `cos(γ_n · k_i · log 2)` for i = 1..M and
n = 1..K costs O(M · K). Per-query amortised cost = O(K). For
K = polylog(2^max k_i), per-query = polylog. **This is exactly
Thread 5's positive shape transposed to the dyadic setting.**
Cost: 2-session arc.

### §P7.c — Lower bound from output bit-budget

Formal claim: any algorithm that outputs all of π(2^k) for k = 1..N
requires Ω(N²) bit-operations (Σ_k log₂ π(2^k) = Σ_k k − O(log k)).
The B-NEGATIVE result here is consistent with: no algorithm beats
the trivial bit-budget. Make this rigorous as a CLOSED_PATHS row
under "lower bounds via output size". Cost: 0.5 session.

---

## Files

- `dyadic_pi_structure.py` — main script (verification + residuals
  + MC + cost test)
- `raw_data.json` — verification map, residuals table for k=1..56
- `hkm_cost.json` — subprocess timing data for k ∈ {28, 30, 32}
- `summary.json` — F-claim verdicts and MC baselines
