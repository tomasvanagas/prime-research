# Soundness budget & field-size requirement of the π(x) verification chain (S534)

**File:** `soundness_field_size.py` (`--selftest` 9 groups, `--budget`,
`--sz-atomic`, `--sz-chain`, `--extrapolate`, `--all`).
**Run log:** `soundness_field_size_run.log`.

## Question

The compressed π(x) verification chain (`compressed_layer.run_chain`) is a
polynomial IOP over a prime field **F_q**, Fiat–Shamir'd to non-interactive.
Two *distinct* lower bounds on the prime q decide whether it works at a given x;
only the first was previously pinned down:

1. **Correctness / integer-range.** The trace identities `U·a + R − X = 0` and
   the bit recompositions carry integers up to ~X; if `q ≤ X` a wrong quotient
   can *alias* the truth mod q. **Already CLOSED:** `status/CLOSED_PATHS.md`
   row 498 (the F_{q²} refutation) — the fix is a larger-**characteristic**
   prime, `q > X`, i.e. `BIG_Q=2⁶¹−1` is sound to n≲60. Row 498 explicitly
   notes the element count "only shrinks Schwartz–Zippel error, **orthogonal**
   to the integer-range hole."

2. **Soundness / Schwartz–Zippel.** That orthogonal axis. Every random-challenge
   round binds a polynomial identity of degree ≤ d; a cheating prover
   false-accepts that round with prob ≤ d/q. Over the whole chain there are
   R(x) field challenges carrying a **degree budget B(x) = Σ per-round degrees**.
   Union bound: total field-soundness error ≤ **B(x)/q**. For λ-bit security:
   **log₂ q > λ + log₂ B(x)**.

This cycle measures axis (2) — never previously characterized — and locates which
axis binds at the project's `p(10^100)` target.

## What was measured

All RAM-light / cache-resident (n ≤ 16, working set ≤ a few MB) — deliberately
chosen to run beside the detached large-x reach without DRAM-bandwidth contention.

### (A) Soundness budget B(x), R(x), d_max — deterministic instrumentation

A counting rng (`CountRng`, bit-identical drop-in — selftest [2]) tallies every
field challenge R; a `sumcheck` monkeypatch (per-module, mirroring
`prover_opcount_scaling.py`) tallies the degree budget B = Σ_rounds (nvars·deg)
and d_max. Run on the **real** honest `run_chain`, default config, over BIG_Q.
B, R, d_max are **field-independent** (q=257 == BIG_Q — selftest [4]):

| n | x | K = π(√x) | R (challenges) | B (degree budget) | d_max | B/K = c_BK |
|---|---|---|---|---|---|---|
| 8  | 255   | 6  | 206  | 530   | 6  | 88.3 |
| 10 | 1023  | 11 | 472  | 1282  | 7  | 116.5 |
| 12 | 4095  | 18 | 922  | 2632  | 8  | 146.2 |
| 14 | 16383 | 31 | 1844 | 5512  | 9  | 177.8 |
| 16 | 65535 | 54 | 3654 | 11406 | 10 | 211.2 |

**Structure.** `K = layers = π(√x)` exactly (6,11,18,31,54 = π(15),…,π(255)) —
the K sequential outer layer reductions. `d_max = n−2` (the B2 routing-lookup
degree, polylog = O(log x)). `c_BK = B/K ≈ 15.35·n − 36.2` grows **linearly in
n = log x** (polylog).

**HEADLINE — the polylogs cancel:** `B = K · c_BK ~ (√x/log x)·(log x) =
Θ(√x)`. The `1/log x` in the prime count K cancels the `log x` per-layer round
budget. Measured single-power exponent **α = 0.548 ≈ 0.5** over n=8..16. So the
**soundness budget is Θ(√x)** — the same √x as the cert size (S509) and the
information floor (S511), now on the *field-soundness* face.

> *Honest note on the fit:* a 2-term `α·n + δ·log₂n + c` fit reads α=0.27 at
> small n — the **S511/S515 PNT 1/log x deflation** (overfitting a polylog term
> on few points). We do **not** extrapolate B by a small-n power; we anchor on
> the **exact** π(√x) backbone × the measured c_BK polylog (next section). The
> honest single-power α=0.55 is the √x reading.

### (B) Per-round Schwartz–Zippel law — atomic (deg-3 trace zero-test)

Build the witness once, corrupt one quotient, re-verify with many fresh
challenge sets at primes q>X (the only soundness-bearing knob). The deg-3 trace
zero-test false-accepts at **rate ≈ c/q** with `c = O(deg·rounds)`.

Clean (high count) — n=8, q just over X, ~44–88 accepts/3000 (selftest [6]):
```
cheat=u_value   rate(q=257)=0.0147  rate(q=521)=0.0073   rate·q ≈ 3.8 (const)   slope ≈ −1
```
n=10 sweep, 6000 seeds/q (run log) — false-accept clearly **drops ~1/q** (an
order of magnitude over q=1031→8191), bounded by `c/q` with `c ≈ 4–9`:
```
cheat=u_value      q=1031 rate=0.0043 (26)   q=2053 rate=0.0030 (18)   q=4093 (2)   q=8191 (4)
cheat=r_value      q=1031 rate=0.0060 (36)   q=2053 rate=0.0045 (27)   q=4093 (2)   q=8191 (1)
```
**Honest caveat:** the *fitted* slopes over the full sweep read −1.1…−1.9, but
the high-q points carry only **1–4 accepts/6000** (expected count < 5 ⇒ Poisson
noise dominates and over-steepens the fit). The robust statement is from the
high-count points (≥15 accepts): `rate·q` is bounded (3.8–9.2) and `rate` falls
~1/q; no point decays *slower* than 1/q. A range-violating corruption (`nonbit`,
a 'bit' = 2) is caught at rate ≈ 0 even at the smallest q (it fails the bit-range
part w.h.p. — selftest [7]).

### (C) Chain composition — two soundness regimes, no amplification

Run the **real** chain with a lie over many seeds at q>X. Two regimes:

| lie | mechanism | false-accept (n=8, 800 seeds/q) |
|---|---|---|
| `delta_pi` (claim π(x)+1) | output-value lie, caught **deterministically** at the final base equality | **0/800** at q=257, 1021, 4091 (no field dependence) |
| `corrupt_layer` (self-consistent internal DP lie) | caught only by **Schwartz–Zippel** | 0.0088 (7) → 0.0038 (3) → 0.0025 (2); `rate·q` 2.3, 3.8, 10.2; always **≤ B/q** |

(The q=4091 `rate·q=10.2` is 2/800 Poisson noise; the trend is a clear drop with q.)

So lying about the *answer* is free to catch; a prover must lie *internally and
self-consistently* to get any Schwartz–Zippel chance, and even then the measured
false-accept stays **under the union bound B/q** — no soundness amplification
past the union bound (falsifier F3 un-triggered). This is the first empirical
chain-soundness audit at **non-negligible** error; prior cheat tests are
single-shot at BIG_Q (error ≈ 0, where a 10× inflation would be invisible).

### (D) Which axis binds — field size vs x

Anchoring on the exact π(√x) backbone × measured c_BK (λ=100):

- **Soundness curve:** `log₂ q > λ + log₂ K(x) + log₂ c_BK ~ λ + n/2 − log₂ n + O(log n)`.
- **Correctness curve:** `log₂ q > n + slack`, measured **slack ≈ 0–1 bit**
  (q just above X accepts honestly: q=257 for X=255).

| | crossover | at p(10^100), n = 332 |
|---|---|---|
| | **n ≈ 209 (x ≈ 10^63)** | log₂K(√x)=159, log₂B=172 |
| soundness needs | — | **q ≈ 272 bits** |
| correctness needs | — | **q ≈ 333 bits** |
| **binds** | — | **CORRECTNESS** |

**Crossover at n ≈ 2λ.** Below it (x ≲ 10^63) the λ-bit security floor binds;
above it — **including the `p(10^100)` target** — the integer-range/correctness
requirement binds at `log₂ q ≈ n ≈ 333 bits`, and that field is **automatically
more than large enough** for 100-bit soundness (272 < 333). The √x soundness
budget never overtakes the linear-in-n correctness requirement at large x.

## Takeaway

- The π-verification chain's **field-soundness budget is Θ(√x)** (the K=π(√x)
  layers × a polylog per-layer round budget; the prime-count 1/log x cancels the
  per-layer log x). It is the *same* √x as the cert size (S509) and info floor
  (S511), now confirmed on the field-soundness axis.
- For `p(10^100)`, the binding field-size constraint is **correctness** (q > X,
  ~333 bits), not soundness (~272 bits for 100-bit security). They cross at
  x ≈ 10^63 (n ≈ 2λ).
- This sharpens the "Õ(√x) verifier" honesty: the field elements are
  `Θ(log x)`-bit (correctness-forced), so each field op is `O((log x)^{1+ε})`
  bit-ops — the polylog in `Õ(√x)` concretely carries an `~n^{1+ε} = (log x)^{1+ε}`
  factor, set by correctness, **independent of the security parameter** once
  x ≳ 10^63.
- Empirically, the protocol's soundness is **real and tight**: per-round error
  = deg/q (slope −1), chain false-accept ≤ B/q with no amplification, and
  output-value lies are caught deterministically.

## What would falsify this (pre-registered)

- **F1** B(x) single-power exponent reads clearly > 0.6 (with the π(√x)-anchored
  c_BK growing super-linearly in n) → soundness budget is not Θ(√x); the
  crossover analysis is wrong. *(Observed α=0.55, c_BK linear → un-triggered.)*
- **F2** atomic accept-rate does not halve when q doubles (slope far from −1)
  → the per-round Schwartz–Zippel model fails in this implementation.
  *(Observed slope ≈ −1, rate·q const → un-triggered.)*
- **F3** a chain cheat's measured false-accept exceeds B/q at some q>X
  → union bound violated (amplification / skipped check). *(All ≤ B/q → un-triggered.)*
- **F4** honest `run_chain` rejects at some prime q ∈ (X, BIG_Q] → correctness
  needs q strictly larger than X (slack > the measured ~1 bit), moving the
  crossover. *(q=257>X=255 accepts → un-triggered.)*

## Honest scope

This is an **adjacent / verification-line** result, **not** progress on the
polylog π(x) goal — the goal stays blocked. It is a quantitative *honesty
sharpening* of the verification stack's "Õ(√x)" claim (the hidden field-size
factor) plus the first empirical soundness audit of the chain at non-negligible
error, and it answers the concrete "what field does `p(10^100)` verification
need" question (≈ 333-bit prime, correctness-bound). The soundness budget being
Θ(√x) is a *new measurement* of a known √x wall from the field-soundness angle —
a measurement, not novelty (per the contract). The union bound is a standard
(loose) composition; the empirical audit confirms the implementation does not
violate it, it does not prove a *matching* soundness lower bound.
