# Space-Time Tradeoff Lower Bounds for pi(x): Results

**Date:** 2026-04-05 (Session 20)
**Code:** `experiments/circuit_complexity/space_time_tradeoff.py`

---

## Executive Summary

We investigated whether formal space-time tradeoff lower bounds can be established
for computing pi(x). **The answer is NO with current techniques**, but we clarified
exactly why and obtained several useful intermediate results.

**Key finding:** The strongest provable bound is T * S >= Omega(log^2(x)), which
does NOT rule out polylog(x) algorithms. Proving T >= x^{Omega(1)} would require
circuit lower bounds beyond the current frontier of complexity theory.

---

## Experiment 1: OBDD Sizes

OBDD (Ordered Binary Decision Diagram) sizes for pi(x) output bits, under the
best variable ordering found:

| N  | max pi(x) | output bits | Total OBDD | Random baseline | Ratio  |
|----|-----------|-------------|------------|-----------------|--------|
| 4  | 6         | 3           | 18         | 12              | 1.500  |
| 6  | 18        | 5           | 58         | 53              | 1.088  |
| 8  | 54        | 6           | 194        | 192             | 1.010  |
| 10 | 172       | 8           | 533        | 819             | 0.651  |
| 12 | 564       | 10          | 1442       | 3413            | 0.422  |

**Growth fit:** log2(OBDD) ~ 0.792 * N => OBDD size ~ 2^{0.79*N} = x^{0.79}

**Interpretation:**
- OBDD size grows exponentially in N = log(x), i.e., polynomially in x.
- For N >= 10, pi(x) OBDDs are SMALLER than random functions of the same density.
  This suggests some exploitable structure in the variable ordering.
- However, OBDD is a restricted model. General branching programs can be
  exponentially smaller than OBDDs.

## Experiment 2: Nechiporuk Lower Bound

Using rank(M_pi, k-bit partition) = 2^{k-1} + 2 (verified Session 17):

The Nechiporuk bound partitions N input bits into blocks of size s and sums
subfunctions per block. The optimal block size is s* = 2/ln(2) ~ 3, giving:

**L(pi) >= Omega(N) = Omega(log x)**

This is a trivially weak bound. The Nechiporuk method is inherently limited to
O(N^2) formula size bounds for ANY function, and for pi(x) it gives only O(N).
The exponential rank 2^{s/2} cannot be leveraged because it grows with block size,
not with the number of blocks.

**Note:** When block size = N (single block), Nechiporuk gives L >= 2^{N/2}/N,
but this is just the communication complexity in disguise, not a formula size bound.

## Experiment 3: Abrahamson-Beame Time-Space Tradeoff

From D(pi) = N/2 - O(1) bits of communication complexity:

| Method | Bound | Sufficient for polylog? |
|--------|-------|------------------------|
| Crossing sequence (Cobham) | T >= N + (N/2)/S | Trivially polylog |
| R-round comm (Beame) | T * S >= Omega(N^2) | YES -- polylog works |
| Multi-partition (Borodin-Cook) | T * S >= Omega(N^2) | YES -- polylog works |
| OBDD size from rank | OBDD >= 2^{N/2-1} = sqrt(x) | Only for OBDDs |

**Critical insight:** Communication complexity of f: {0,1}^N -> Z is bounded by
N bits (the input length). This means D(f) = O(N) always, so any time-space
tradeoff derived from communication complexity gives T * S >= poly(N) at best.
Since poly(N) = polylog(x), **communication complexity can never rule out
polylog(x) algorithms**.

To get T >= 2^{Omega(N)} = x^{Omega(1)}, we need circuit lower bounds.

## Experiment 4: Pebbling the Meissel-Lehmer DAG

DAG parameters (Lucy DP formulation):

| x       | floor vals | primes | nodes  | depth | width |
|---------|-----------|--------|--------|-------|-------|
| 100     | 19        | 4      | 95     | 4     | 19    |
| 1,000   | 62        | 11     | 744    | 11    | 62    |
| 10,000  | 199       | 25     | 5,174  | 25    | 199   |
| 100,000 | 631       | 65     | 41,646 | 65    | 631   |

Depth = pi(sqrt(x)), confirming Session 12.

**Exact pebbling results (small DAGs):**
- x=10: 4 pebbles (15 nodes, depth 2, width 5)
- x=15: 5 pebbles (18 nodes, depth 2, width 6)
- x=20: 5 pebbles (24 nodes, depth 2, width 8)

**Pebbling space-time tradeoff (M-L DAG only):**
T * S >= Omega(x^{5/6} / ln x) for the specific Meissel-Lehmer computation DAG.

This is ALGORITHM-SPECIFIC. It says: if you compute pi(x) using the
Meissel-Lehmer method, you cannot beat T * S ~ x^{5/6}. But a completely
different algorithm could potentially bypass this DAG entirely.

## Experiment 5: BDD Growth Comparison

| N  | isPrime | pi mod 2 | pi total | random | 2^{N/2} |
|----|---------|----------|----------|--------|---------|
| 4  | 6       | 8        | 19       | 6      | 4       |
| 6  | 16      | 19       | 60       | 17     | 8       |
| 8  | 43      | 55       | 197      | 56     | 16      |
| 10 | 122     | 147      | 533      | 163    | 32      |
| 12 | 344     | 421      | 1442     | 484    | 64      |

Growth exponents (log2(OBDD) ~ c * N):
- isPrime: c = 0.731
- pi mod 2: c = 0.719
- pi(x) total: c = 0.782
- random: c = 0.796

**All grow as 2^{~0.7-0.8 * N}**, confirming that pi(x) behaves like a
pseudo-random function for OBDD complexity. No significant structural advantage
over random functions.

---

## Formal Results Summary

### What IS proven:
1. **OBDD lower bound:** Any OBDD for pi(x) has size >= 2^{N/2-1} = Omega(sqrt(x))
2. **Time-space (branching programs):** T * S >= Omega(N^2) = Omega(log^2 x)
3. **Pebbling (M-L DAG only):** T * S >= Omega(x^{5/6} / ln x) for Meissel-Lehmer
4. **Formula size:** L(pi) >= Omega(N) = Omega(log x) [trivial]

### What CANNOT be proven (with current techniques):
1. T >= x^{Omega(1)} for general algorithms
2. T * S >= x^{Omega(1)} for general algorithms
3. pi(x) is not in NC (polylog depth, poly size circuits)

### Why:
The fundamental obstacle is the **circuit complexity barrier**. Communication
complexity arguments are bounded by the input length N = log(x), so they can
only give poly(N) = polylog(x) bounds. To prove super-polynomial (in N) bounds
would require proving circuit lower bounds, which faces:
- The **Natural Proofs barrier** (Razborov-Rudich 1997)
- The **relativization barrier** (Baker-Gill-Solovay 1975)
- The **algebrization barrier** (Aaronson-Wigderson 2009)

### Status update for OPEN_PROBLEMS:
The space-time tradeoff problem for pi(x) should be reclassified:
- **General algorithms:** OPEN, and likely as hard as P vs NC.
- **OBDD model:** CLOSED -- Omega(sqrt(x)) proven from communication rank.
- **Meissel-Lehmer specifically:** CLOSED -- T*S >= Omega(x^{5/6}/ln x) from pebbling.
- **Communication complexity route:** CLOSED as a route to super-polylog bounds.
  Maximum achievable: T*S >= Omega(log^2 x), which is too weak.
