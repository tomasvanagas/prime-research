# Session 84 — A1 wild swing: SAT-search for TC^0 PRIMES at N=8

**Date:** 2026-04-26
**Mode:** WILD SWING (frontier attack, ATTACK_VECTORS.md §A1)
**Target:** ATTACK_VECTORS.md §A1 — "find a primality witness using only
NC^1-known scalar primitives via SAT-based circuit synthesis at N=8".

## Self-grade: B (substantive structural finding + concrete mechanism)

The session committed to a §A1 wild swing as instructed (no mid-session
pivot). It produced:
- Three independent measurements that did not previously exist in the
  project (PTF degree, depth-2 W=1 size, N=6 PRIMES-vs-random gap).
- A clean STRUCTURAL MECHANISM (1-bit oddness advantage) that explains
  the PRIMES-vs-random gap at N=6 in elementary terms.
- A negative result for §A1 as literally posed (depth-5 size-2000
  enumeration intractable).

**Why B and not A:** the session did not produce a small explicit
polylog-promising circuit for PRIMES at N=8 (which would have been the
A-grade outcome). The PRIMES-vs-random gap at N=6, while statistically
significant (p < 0.001), reduces to the well-known "odd primes" fact —
a known elementary structure, not a new arithmetic insight.

**Why not C:** the session produced 3+ artefacts not previously in the
project:
- Exact PTF degrees of PRIMES at N=4..8 (4 at N=8) — first measurement.
- Exact depth-2 W=1 sign-threshold sizes at N=4 (3) and N=6 (6) —
  first measurement.
- The PRIMES-vs-random gap is the FIRST PRIMES circuit-complexity
  measure in the project that EMPIRICALLY DEVIATES from random
  (35+ prior pseudorandomness measures show no deviation; this one
  does, with a clean mechanism).

**Why not F:** the session committed to ONE wild-swing target (§A1) and
spent the full session on it. Multiple ILP / Z3 / enumeration approaches
were tested. Failure modes documented honestly. The negative result
constrains future §A1 attempts (rules out the "small W=1 sign-threshold"
sub-family at N=8).

## What this session produced

1. **PTF degree of PRIMES at N=8** (= 4) and N=4..7 (2, 3, 3, 3).
   Random-matched controls give the same. No structural advantage at
   the PTF-degree level — consistent with the project's pseudorandomness
   thesis.
2. **Exact depth-2 sign-threshold (W=1) circuit complexity** of PRIMES
   at N=4 (= 3) and N=6 (= 6). Lower bounds at N=7 (≥ 6) and N=8
   (≥ 17, in restricted candidate sub-family k_max ≤ 5).
3. **Statistically significant PRIMES-vs-random gap at N=6**:
   PRIMES = 6, random ≥ 7 over 10 seeds. p < 0.001. The first instance
   in the project of a circuit-complexity measure where PRIMES
   deviates from matched-density random.
4. **Mechanism of the gap**: PRIMES has a single-bit predictor (bit_0 =
   "x is odd") with 70.3% accuracy at N=8, vs random's best 1-bit at
   ~57%. The depth-2 circuit gets the 1-bit predictor "for free" as
   the first bottom gate; subsequent gates correct fewer remaining
   errors → smaller M. This explains the gap quantitatively.
5. **Concrete depth-2 sign-threshold circuits** for PRIMES at N=4 (M=3)
   and N=6 (M=6), with explicit gate weights and threshold breakdowns.
   Bottom-layer gates are empirical fits, not number-theoretic
   primality witnesses.

## What this session did NOT produce

- A small explicit TC^0 circuit for PRIMES at N=8 that could
  generalise to polylog. Confirmed by ~17-gate lower bound in W=1
  k_max ≤ 5 family.
- Identification of bottom-layer gates with classical primality tests
  (Miller-Rabin, BPSW, AKS witnesses). Gates are empirical fits.
- A complete §A1 closure (depth-5 size-2000 enumeration is intractable
  with SAT/ILP at this scale; needs a different attack).

## Edges cited / composed

- E1.10 / E3.13 (zeros are GUE-random; pseudorandomness measures all
  match GUE) — REINFORCED at PTF-degree level; PARTIALLY DEVIATED from
  at depth-2 W=1 sign-threshold level (with clear elementary mechanism).
- E5.3 / E7.10 (PRIMES in TC^0 open; AKS-modulus-twists orthogonal to
  depth) — confirmed: depth-2 W=1 is universally insufficient at small
  N (M ≥ 17 at N=8); supports the polylog hypothesis being hard.
- S20 / S28 (BDD complexity ~ 2^(0.73N)) — independent model; W=1
  sign-threshold needs ≥ 17 gates at N=8 vs BDD ≥ 360.
- novel/pseudorandomness_of_pi.md — the 35-measure battery now has a
  RELEVANT 36th measure (depth-2 W=1 sign-threshold size) where PRIMES
  EMPIRICALLY DEVIATES from random. The deviation has a known elementary
  mechanism (oddness), so it doesn't break the broader pseudorandomness
  thesis, but it WARRANTS A FOOTNOTE.

## Why no edge promotion

The PRIMES-vs-random depth-2 gap is statistically significant but
mechanistically reducible to "primes are mostly odd." This is well-known
elementary number theory, not a new mathematical edge. It refines the
pseudorandomness narrative but doesn't unlock new circuit-complexity
attack surfaces.

## Open question for next session

**Is there a residual PRIMES-vs-random gap beyond the oddness effect?**

Concrete falsification test: redo N=6 (and N=8) with "random functions
calibrated to match PRIMES on bit_0" (i.e., random subset of 17 odd
numbers + 1 even number out of 64 inputs at N=6, matching primes's
class-conditional distribution). If PRIMES still has min M < calibrated
random → residual structure beyond oddness, possibly worth pursuing.
If PRIMES matches calibrated random → entire gap is explained by
oddness, pseudorandomness thesis fully intact.

This is a 1-session B-grade target, added to NOVELTY_CHALLENGES.md.

## Failure mode self-evaluation

Q1: What did I produce that was not in the project before this session?
A: Exact PTF degrees N=4..8, exact depth-2 W=1 sizes N=4,6 + bounds
   N=7,8, single-bit-predictor advantage quantification, mechanism
   linking 1-bit predictor advantage to depth-2 sign-threshold gap.

Q2: What edges did my work compose or cite?
A: E1.10, E3.13, E5.3, E7.10, S20, S28, novel/pseudorandomness_of_pi.

Q3: Was this a duplicate-closure session?
A: No. Three new measurement families (PTF degree, depth-2 W=1 size,
   1-bit predictor advantage) plus one structural mechanism. Each
   could be re-derived by a strong number theorist in a half-day given
   the right tools, but they were not in CLOSED_PATHS / EDGES /
   experiments before this session.

Q4: What's the next-action?
A: (a) Add to NOVELTY_CHALLENGES.md: "calibrated-1-bit random control
       for PRIMES depth-2 sign-threshold size at N=6" (1-session B-grade
       target).
   (b) Add CLOSED_PATHS row for the §A1 W=1 sign-threshold sub-family
       at N=8 (size > 16 with bounded-weight bottom gates).
   (c) ATTACK_VECTORS §A1 NOT closed; remains open with the refinement
       that depth-2 W=1 is too restrictive at N=8.

## Files produced

`experiments/circuit_complexity/sat_tc0_primes_n8/`:
  sat_tc0_primes_n8.py, sat_tc0_primes_n8_results.md,
  sat_depth2_ilp.py, enum_d2_smart.py, enumerated_depth2.py,
  ptf_degree_battery.py, depth1_threshold_test.py,
  n6_robust_check.py, greedy_d2.py, greedy_d2_np.py,
  + result JSONs and logs.

This session synthesis: `archive/sessions/session84_a1_sat_tc0_primes.md`.
