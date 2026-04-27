# C7 — Calibrated 1-bit-bias random control: definition

## Object

A randomised Boolean function `f_cal : {0,1}^N → {0,1}` whose
*joint* distribution with the parity bit `x_0` matches that of
PRIMES at N=6 exactly, but is otherwise unstructured.

## Signature

For N = 6, π_table(x) := 1 iff x ∈ Primes ∩ {0..63}, with weight 18:
- `n_p_odd  := |{x : x odd, π_table(x)=1}|  = 17`
- `n_p_even := |{x : x even, π_table(x)=1}| = 1` (the prime 2)

**Stratified mode (primary):**
```
f_cal(x) := indicator{x ∈ S_odd ∪ S_even}
where  S_odd  ⊆ {odd   x ∈ [0,63]}, |S_odd|  = 17, S_odd uniformly random
       S_even ⊆ {even  x ∈ [0,63]}, |S_even| = 1,  S_even uniformly random
```
Always: `weight(f_cal) = 18`, `bit_0_accuracy(f_cal) = 0.75`.

**Bernoulli mode (secondary):**
```
f_cal(x) ~ Bernoulli(p_odd)   if x odd, with p_odd = 17/32
f_cal(x) ~ Bernoulli(p_even)  if x even, with p_even = 1/32
```
Mean weight 18, mean bit_0_accuracy 0.75, but each individual sample
fluctuates.

## Intended relationship to π(x)

`f_cal` is a class-conditioned matched-density random control for
PRIMES, designed to isolate the bit_0 ("oddness") component of any
PRIMES-vs-random gap.

In particular: the depth-2 sign-threshold W=1 circuit complexity of
`f_cal` should equal that of PRIMES *if and only if* the S84 gap is
entirely driven by the oddness predictor.

## Edges composed

- **E1.10** — pseudorandomness of π(x) mod 2 across 35 measures. The
  S84 depth-2 sign-threshold result is the **36th** such measure; it
  is the first that DEVIATED from matched-density random behaviour.
  This composition partitions the deviation into "explained-by-bit_0"
  and "residual" components.
- **S84 result** — `experiments/circuit_complexity/sat_tc0_primes_n8/`:
  PRIMES min M = 6, unbiased random min M ∈ {7,8} over 10 seeds at N=6.
  This composition uses S84's enumeration + ILP harness verbatim with a
  modified target-table generator.

## Output reported

For each of 20 stratified + 20 bernoulli samples:
- weight (count of 1s in the truth table)
- bit_0 accuracy (fraction of x where (1 if x odd else 0) == f_cal(x))
- min_M = smallest M such that depth2_search(f_cal, N=6, K=1458, M)
  returns a verified SAT solution

Output histogram of min_M for each mode and compare to:
- PRIMES baseline: min_M = 6
- Unbiased matched-density random baseline: min_M ∈ {7, 8} (S84)

## Why this is novel

- No previous project session has measured the depth-2 sign-threshold
  size of any class-conditional matched random function (calibrated or
  otherwise). The S84 result was a single-mechanism conjecture; this
  experiment tests it directly.
- The construction mechanically separates "structure-from-oddness" from
  "structure-beyond-oddness" in the PRIMES sample, using a control
  whose bit_0 statistic is matched exactly. No published work performs
  this separation for the depth-2 sign-threshold metric.

## Falsification

See `calibrated_d2_results.md` for the four pre-stated falsifiers
(F1-F4).
