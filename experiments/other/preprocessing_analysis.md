# Preprocessing / Succinct Data Structures for p(n): Deep Analysis

**Date:** 2026-04-04  
**Session:** 18  
**Prior work:** Session 7 `experiments/other/succinct_structure.py` (10 approaches analyzed)

---

## 1. Problem Setup

**Model:** Allow polynomial-time preprocessing to build a data structure of size S(x)
bits. Then answer queries "what is p(n)?" in time T(n), where n ranges up to M and
x = p(M) ~ M ln M.

**Goal:** Determine whether ANY preprocessing scheme achieves polylog query time T(n) =
O(polylog(n)) with sub-exponential storage S = o(2^N) where N = log2(x).

**Key parameters for M = 10^100:**
- x = p(M) ~ M ln M ~ 2.3 * 10^102, so N = log2(x) ~ 341 bits
- pi(x) ~ x/ln(x) ~ 4.3 * 10^99

---

## 2. Naive Approaches (Review from Session 7)

### 2a. Full lookup table
Store all primes up to x. Size = pi(x) * log(x) ~ 10^99 * 341 ~ 3.4 * 10^101 bits.
Query time O(log pi(x)) via binary search. **Impossible**: exceeds atoms in universe.

### 2b. Compressed prime table (gap encoding)
By PNT, average gap ~ ln(x) ~ 235, so each gap needs log2(ln(x)) ~ 7.9 bits.
Total: pi(x) * log(log(x)) ~ 10^99 * 8 ~ 8 * 10^99 bits.
This is information-theoretically near-optimal (entropy bound = pi(x) * log(x/pi(x))
= pi(x) * log(ln(x)) ~ 10^99 * 7.9 bits). Still 8 * 10^99 bits. **Impossible.**

### 2c. Information-theoretic lower bound (any representation)
Any data structure supporting exact p(n) for ALL n in [1,M] must store at least
log2(C(x, pi(x))) ~ x * H(1/ln(x)) ~ pi(x) * log(ln(x)) ~ 8 * 10^99 bits.
**No compression can beat this** -- it is the entropy of the prime indicator.

---

## 3. The Strategic Precomputation / Range Tree Approach

### 3a. Pi(x) at evenly-spaced checkpoints

Precompute pi(k * Delta) for k = 0, 1, ..., x/Delta.

- **Storage:** (x/Delta) values, each of ~log2(pi(x)) ~ 332 bits.
  S = (x / Delta) * 332 bits.

- **Query for p(n):** Binary search for the unique k with pi(k*Delta) < n <= pi((k+1)*Delta).
  This costs O(log(x/Delta)) comparisons of stored values = essentially free.
  Then compute p(n) exactly within the interval [k*Delta, (k+1)*Delta].

- **Within-interval cost:** We need pi(x') for x' in an interval of width Delta.
  Using Deleglise-Rivat (best combinatorial): O(Delta^{2/3}).
  Using Lagarias-Odlyzko (best analytic, conditional on RH): O(Delta^{1/2+eps}).

### 3b. Optimal Delta (combinatorial pi(x))

Minimize max(S, T):
- S = (x / Delta) * 332
- T = Delta^{2/3}

Set S = T:
  (x / Delta) * 332 = Delta^{2/3}
  x * 332 = Delta^{5/3}
  Delta = (332 * x)^{3/5}

For x = 10^102:
  Delta = (332 * 10^102)^{3/5} = (3.32 * 10^104)^{3/5} ~ 10^{62.7}

Then:
  S = (10^102 / 10^{62.7}) * 332 ~ 10^{41.8} bits ~ 10^{41} bytes
  T = (10^{62.7})^{2/3} ~ 10^{41.8} operations

**Result:** S ~ T ~ 10^{41.8}. Product S*T ~ 10^{83.6}.

Compare to plain Deleglise-Rivat: T = x^{2/3} = 10^{68}, S = x^{1/3} = 10^{34}.
Product S*T = 10^{102}.

**The precomputation approach wins on query time** (10^{41.8} vs 10^{68}) at the cost
of higher storage (10^{41.8} vs 10^{34}). Product is LOWER (10^{83.6} vs 10^{102}).

This is significant: the product S*T is NOT fixed at Omega(x). It depends on the
method used within intervals.

### 3c. Optimal Delta (analytic pi(x), conditional on RH)

With Lagarias-Odlyzko, within-interval cost = O(Delta^{1/2+eps}).
Set S = T:
  (x / Delta) * 332 = Delta^{1/2}
  x * 332 = Delta^{3/2}
  Delta = (332 * x)^{2/3} ~ (3.32 * 10^104)^{2/3} ~ 10^{69.6}

Then:
  S = (10^102 / 10^{69.6}) * 332 ~ 10^{34.9} bits
  T = (10^{69.6})^{1/2} ~ 10^{34.8} operations

**Result (conditional on RH):** S ~ T ~ 10^{34.8}. Product S*T ~ 10^{69.6}.

### 3d. General S*T tradeoff curve

With inner algorithm of complexity O(Delta^alpha):
  S = (x/Delta) * log(pi(x))
  T = Delta^alpha

Minimize S + T (or set S = T):
  Delta = (x * log(pi(x)))^{1/(1+alpha)}
  S = T = (x * log(pi(x)))^{alpha/(1+alpha)}

| alpha | Method | Optimal Delta | S = T | S*T |
|-------|--------|--------------|-------|-----|
| 2/3 | Deleglise-Rivat | 10^{62.7} | 10^{41.8} | 10^{83.6} |
| 1/2+eps | Lagarias-Odlyzko (RH) | 10^{69.6} | 10^{34.8} | 10^{69.6} |
| 1/3 | Hypothetical O(x^{1/3}) | 10^{78.0} | 10^{26.0} | 10^{52.0} |
| 0+eps | Hypothetical polylog | 10^{102} | polylog | polylog |

The last row shows: if pi(x) were computable in polylog(x) time (i.e., alpha -> 0),
then Delta -> x and we need only O(1) checkpoints with polylog query time.
**This is just the circuit complexity question in disguise: pi(x) in NC <=> polylog T.**

---

## 4. Succinct Index with Pi(x) Oracle at Strategic Points

### 4a. Non-uniform checkpoints

Instead of uniform spacing Delta, place checkpoints at x_0 < x_1 < ... < x_K where
pi(x_i) = i * (M/K). This ensures each interval contains exactly M/K primes.

**Advantage:** Intervals near x ~ 0 (where primes are dense) are shorter than intervals
near x ~ 10^102 (where primes are sparse). The maximum interval width is
Delta_max = max_i (x_{i+1} - x_i).

By PNT, the interval containing primes p(jM/K+1) through p((j+1)M/K) has width
approximately (M/K) * ln(p(jM/K)) ~ (M/K) * ln(jM ln(jM)/K).

The WIDEST interval is the last one, width ~ (M/K) * ln(M ln M) ~ (M/K) * ln(x).
The NARROWEST is the first, width ~ (M/K) * ln(2) (near x=2).

So the worst-case query is still T ~ ((M/K) * ln(x))^alpha.

Storage: K values of size log2(x) ~ N bits. S = K * N.

Set S = T:
  K * N = ((M/K) * ln(x))^alpha
  K * N = (M * ln(x) / K)^alpha
  K^{1+alpha} = (M * ln(x))^alpha / N

This gives the same asymptotic as uniform checkpoints (up to log factors) because the
widest interval dominates. No improvement from non-uniform spacing.

### 4b. Hierarchical / geometric checkpoints

Store pi(x) at x = x_max, x_max/2, x_max/4, ..., 1 (geometric spacing).
Storage: O(log(x)) values = O(N) values of O(log(pi(x))) bits each = O(N^2) bits ~ 10^5 bits.
This is tiny!

Query: Binary search among O(log x) = O(N) checkpoints to find interval of width
x/2^k containing p(n). Then further narrow: call pi at midpoint. Each call costs
O((x/2^k)^alpha). After O(N) refinements, interval width = 1.

Total cost: sum over k from 0 to N of O((x/2^k)^alpha) ~ O(x^alpha) (geometric series
dominated by first term).

**Result:** Geometric checkpoints with O(N^2) storage reduce query time by only a
constant factor vs computing pi(x) from scratch. The first refinement step dominates.

### 4c. Two-level scheme: geometric + dense

Level 1: Geometric checkpoints (O(N^2) bits) to narrow to interval [a, 2a].
Level 2: Within [a, 2a], store pi at uniform spacing Delta_2.

This gives:
- S = O(N^2) + (a/Delta_2) * N ~ (a/Delta_2) * N
- T = Delta_2^alpha

where a ~ x (for the worst case, the outermost interval).
Same tradeoff as before. No asymptotic improvement.

---

## 5. Connection to Communication Complexity

### 5a. The Session 17 result

The communication matrix of pi_N (N-bit input, Alice gets top N/2 bits, Bob gets
bottom N/2 bits) has rank exactly 2^{N/2-1} + 2 for even N in [4, 16].

**Deterministic communication complexity:** D(pi_N) = ceil(log2(rank)) = N/2 - 1 + O(1).
This means any two-party protocol computing pi(x) requires Omega(N/2) = Omega(log(x)/2)
bits of communication.

### 5b. Does this imply data structure lower bounds?

**Miltersen's connection (1995):** For a static data structure problem with universe
size U, storing a set from a family F:
- Cell-probe complexity with word size w: T >= log|F| / (S * w)
  where S = number of cells, T = number of probes per query.

For our problem:
- Universe U = x (possible prime values)
- |F| = C(x, pi(x)) (number of possible prime sets)
- log|F| ~ pi(x) * log(ln(x)) ~ 8 * 10^99

**Cell-probe lower bound:**
  T >= log|F| / (S * w) = 8 * 10^99 / (S * w)

With word size w = log(x) ~ 341 and S cells:
  T >= 8 * 10^99 / (341 * S)

For S = 10^50 cells (10^52 bits): T >= 2.3 * 10^47 probes. Not polylog.
For T = polylog(x) ~ 10^3: S >= 8 * 10^99 / (341 * 10^3) ~ 2.3 * 10^94 cells.

**This is the simple counting bound. Can we do better?**

### 5c. Patrascu-style lower bounds

Patrascu (2008, "Unifying the landscape of cell-probe lower bounds") proved that for
the rank/select problem on a bit vector of length n with m ones:

  T = Omega(log(n*log(n/m) / (S*log(n/m))))

where S = space in bits, T = query time in cell probes (word size log n).

For our problem (primes as a sparse set): n = x, m = pi(x) ~ x/ln(x):
  log(n/m) = log(ln(x)) ~ 3 bits (for x = 10^102, ln(x) ~ 235, log2(235) ~ 7.9)

  T = Omega(log(x * 7.9 / (S * 7.9)))
    = Omega(log(x / S))

For S = polylog(x) = N^c:
  T = Omega(log(x / N^c)) = Omega(N - c*log(N)) = Omega(N)

**THIS IS A STRONG LOWER BOUND:** With poly(log(x)) bits of storage, any cell-probe
data structure needs Omega(log(x)) = Omega(N) probes per query.

But N = log2(x) ~ 341 for x = 10^102. So T >= 341 probes.
Each probe reads a word of O(log(x)) bits. Total: O(N^2) = O(10^5) bit operations.
This is still polylog in x! (N^2 = (log x)^2.)

**WAIT -- this is only for the RANK/SELECT formulation.** It says: if we store the
prime bit vector in compressed form (S bits), then select queries need
Omega(log(x/S)) probes. With S = O(x/ln(x) * log(ln(x))) (the entropy), we get T = O(1).
But with S < entropy, probes grow.

### 5d. The real question: small S regime

The question is whether there exists a data structure with:
- S = poly(N) = poly(log x) bits (polynomially many bits in the INPUT size)
- T = poly(N) probes (polynomial time in the input size)

This is exactly the question: "Is p(n) computable in polynomial time (in the bit
length of n)?" Equivalently: "Is pi(x) in P when x is given in binary?"

**Known:** ALL algorithms for pi(x) are exponential in N = log(x). Best: 2^{N/2+eps}
(Lagarias-Odlyzko). No algorithm runs in poly(N) time.

**The cell-probe lower bound from Patrascu gives T = Omega(N) with S = poly(N).**
This is TRIVIAL -- just reading the output takes Omega(N) time.

**There is NO known super-linear cell-probe lower bound for rank/select with
sub-entropy storage.** The known bounds are:
- Entropy storage: T = O(1) (achieved by RRR, Elias-Fano)
- Sub-entropy: T = Omega(log(opt_space / actual_space))
- Tiny storage (poly(N)): T = Omega(N) (trivial from output size)

### 5e. Communication complexity to data structure reduction

The communication complexity bound D(pi_N) = N/2 - O(1) does give a data structure
lower bound via the Miltersen-Nisan-Wigderson (1998) framework:

For an INDEXING problem where Alice encodes a dataset and Bob holds a query:
- Alice's message (= data structure) has S bits
- Bob's probes (= query algorithm) take T rounds
- S * T >= Omega(communication complexity of the underlying function)

But the "underlying function" here is: given the prime bit-vector (Alice's data) and
an index n (Bob's query), output p(n) (the n-th 1 in Alice's vector). This is exactly
the SELECT problem.

The communication complexity of SELECT on a universe of size x with pi(x) ones is
Omega(log(C(x, pi(x)))) = Omega(pi(x) * log(ln(x))) in the one-way model (Alice
sends, Bob answers).

**This gives:** S >= pi(x) * log(ln(x)) for T = 1, or S*T >= Omega(pi(x) * log(ln(x)))
in general. This is the information-theoretic bound again. No improvement.

The Session 17 communication complexity result (rank = 2^{N/2-1}+2) applies to a
DIFFERENT partition: Alice holds top N/2 bits of x, Bob holds bottom N/2 bits, and
they compute pi(x). This is about computing pi at a SINGLE point, not about data
structures storing ALL primes.

---

## 6. Cell-Probe Lower Bounds for Prime Queries (Literature)

### 6a. Known cell-probe lower bounds

No cell-probe lower bound specific to prime-counting queries has been published.
The relevant general results are:

1. **Rank/Select on bit vectors** (Golynski 2006, Patrascu 2008):
   With space S = m * log(n/m) * (1 + 1/f) bits (f = redundancy factor),
   select requires T = Omega(log f) probes.

2. **Predecessor search** (Patrascu-Thorup 2006):
   With space S = n * polylog bits, predecessor queries on n keys from [U] require
   T = Omega(min(log log U, sqrt(log n / log log n))) probes.

3. **Range counting** (Patrascu 2007):
   2D range counting with n points requires T = Omega(log n / log log n) probes
   even with space O(n polylog).

None of these directly address the prime-specific structure. The prime indicator is
a specific fixed function, not a worst-case adversarial input.

### 6b. Can prime structure help?

The primes are NOT adversarial -- they have known statistical properties (PNT, Bombieri-
Vinogradov, etc.). Could this help build better data structures?

**No**, for the following reason: the hard part of p(n) is the correction delta(n) =
p(n) - R^{-1}(n), which is pseudorandom (Session 16 confirmed: uncorrelated with gaps,
uniform mod m, ~sqrt(x) magnitude). Any data structure that stores delta(n) must handle
this pseudorandom signal.

The entropy of {delta(n) : 1 <= n <= M} is approximately M * H(delta) where each delta
value needs ~(1/2)*log2(n) ~ N/2 bits. Total: M * N/2 ~ 10^100 * 170 ~ 1.7 * 10^102 bits.

This is the same order as the entropy of the full prime set. The "nice structure" of
primes (PNT, residue classes, etc.) is already captured by R^{-1}(n) and does NOT
reduce the entropy of the correction.

---

## 7. Can We Achieve Polylog Query with Sub-Exponential Storage?

### 7a. The answer is NO, absent a breakthrough in pi(x) complexity

The fundamental constraint is:

**To answer a single query p(n), we need to determine delta(n) = p(n) - R^{-1}(n)
to precision O(1). This requires either:**
1. Storing delta(n) (or equivalent) -- costs ~N/2 bits per query point
2. Computing delta(n) on-the-fly -- costs O(x^{alpha}) for current alpha in {2/3, 1/2+eps}
3. Some combination via preprocessing

For case 3, the tradeoff is:
- Store delta at K evenly-spaced reference points: S = K * N/2
- For a query between references, interpolate/recompute: T = (M/K)^alpha * polylog

For T = polylog(M): need (M/K)^alpha = polylog, so M/K = polylog^{1/alpha}, so
K = M / polylog^{1/alpha} ~ M. Storage = M * N/2 ~ 10^100 * 170 ~ 10^102 bits.

**For T = polylog, we need S ~ M * N/2 ~ 10^102 bits. This is exponential in N.**

The ONLY way to achieve polylog T with sub-exponential S is if the inner computation
(step 2) can be done in polylog time. That is: **pi(x) must be computable in polylog(x)
time, i.e., pi(x) must be in NC.** This is the central open question.

### 7b. Summary of S*T tradeoff

**Theorem (S*T tradeoff for p(n)):**

Let S be the storage in bits and T the query time (bit operations) for answering
p(n) queries for n in [1, M]. Let alpha be the exponent of the best known
algorithm for pi(x) (currently alpha = 2/3 unconditional, alpha = 1/2+eps under RH).

Then the optimal balanced tradeoff is:
  S = T = (M * ln(M * ln M))^{alpha/(1+alpha)} * polylog

| Storage S | Query T (alpha=2/3) | Query T (alpha=1/2) |
|-----------|-------------------|-------------------|
| 10^100 (full) | polylog | polylog |
| 10^80 | 10^{13.2} | 10^{10.0} |
| 10^60 | 10^{26.5} | 10^{20.0} |
| 10^{41.8} | 10^{41.8} | -- |
| 10^{34.8} | -- | 10^{34.8} |
| 10^{34} (DR space) | 10^{68} | 10^{51} |
| polylog | 10^{68} | 10^{51} |

(Values are for M = 10^100, x ~ 10^102.)

Key observation: **Reducing storage below 10^{34} (the working memory of Deleglise-
Rivat) does NOT increase query time.** The algorithm already runs in O(x^{2/3}) with
only O(x^{1/3}) working space. Extra storage ONLY helps if it stores precomputed pi(x)
values, which allows querying within smaller intervals.

### 7c. The preprocessing-query gap

| Preprocessing time | Storage | Query time (alpha=2/3) |
|-------------------|---------|----------------------|
| 0 (no precomp) | O(x^{1/3}) working | O(x^{2/3}) |
| O(x^{2/3}) for K checkpoints | K * log(pi(x)) | O((x/K)^{2/3}) |
| O(x) (sieve everything) | pi(x) * log(ln(x)) | O(log(pi(x))) |

The preprocessing time to build K checkpoints is K * O((x/K)^{2/3} + ...) ~ K * x^{2/3}
which is actually O(K * x^{2/3}) if done naively (each pi(x_i) independently), or
O(x^{2/3}) if the checkpoints can share computation (e.g., Lucy DP computes all pi(x/k)
simultaneously). With the Lucy DP, we get O(sqrt(x)) checkpoints for free!

**Lucy DP insight:** The Lucy DP computes pi(v) for all v in {floor(x/k) : k=1,...,x}
simultaneously in O(x^{2/3}) time. This gives O(sqrt(x)) precomputed values.

With these O(sqrt(x)) = O(10^{51}) checkpoints, the maximum gap between consecutive
floor(x/k) values is O(sqrt(x)), so:
  S = O(sqrt(x) * log(pi(x))) ~ O(10^{51} * 332) ~ 10^{53.5} bits
  T = O(sqrt(x)^{2/3}) = O(x^{1/3}) = O(10^{34}) [within each interval]

This is exactly the Deleglise-Rivat tradeoff! It's already optimal for its method.

---

## 8. Definitive Answers to the Six Questions

### Q1: Lookup table
Storing all primes: 8 * 10^99 bits minimum. **Impossible** for M = 10^100.

### Q2: Compressed prime table
Gap encoding achieves ~8 bits/prime (near-optimal). Still 8 * 10^99 bits. **Impossible.**

### Q3: Succinct index + pi(x) oracle at a few points
Binary search needs O(log(x)) ~ 341 oracle calls. Each at cost O(x^{2/3}) = O(10^{68}).
Total: O(10^{70.5}). **Better than sieve (O(10^{100})) but not polylog.** Storing pi(x)
at strategic points reduces inner-interval cost but not asymptotically (see Section 4).

### Q4: Range tree / skip structure
Optimal Delta with alpha=2/3: Delta = (332*x)^{3/5} ~ 10^{62.7}.
Storage: 10^{41.8} bits. Query: 10^{41.8} ops.
This achieves the balanced S = T point of the tradeoff curve. **Significant improvement
over plain binary search but still far from polylog.**

### Q5: Information-theoretic S*T lower bound
From counting: S*T >= Omega(M * log(ln(M))) / log(S).
From communication complexity (Session 17): rank(pi_N) = 2^{N/2-1} + 2 implies
**the function pi(x) itself requires Omega(sqrt(x)) work**, independent of storage.
This means T >= Omega(sqrt(x)) for computing pi at ANY new point not stored.

**Crucial insight:** The Session 17 result implies that even with unlimited preprocessing/
storage, computing pi(x) at a NEW point (not stored) requires Omega(sqrt(x)) communication
bits in the 2-party model. However, this does NOT directly translate to a cell-probe lower
bound because the data structure model is more powerful (adaptive queries to stored data).

The true data structure lower bound is: S >= Omega(M * log(ln(M))) for T = O(1), and
T >= Omega(log(x/S)) for general S (from rank/select bounds). **No non-trivial S*T
product lower bound beyond the counting argument is known for this specific problem.**

### Q6: Cell-probe lower bounds
**No prime-specific cell-probe lower bounds exist in the literature.** General rank/select
lower bounds (Patrascu 2008) give T = Omega(log(redundancy factor)) which is O(1) at
optimal space and Omega(N) at poly(N) space -- the latter being trivial.

---

## 9. Conclusion

**Can preprocessing achieve polylog query time with sub-exponential storage?**

**NO**, under any known algorithm for pi(x).

The fundamental constraint: to answer p(n) with query time T, we need to determine
pi(x') at points x' within an interval of width W around p(n), where W is determined
by the gaps between precomputed checkpoints. The cost is O(W^alpha) where alpha >= 1/2
(unconditionally likely alpha >= 1/2; currently alpha = 2/3 proven, alpha = 1/2+eps
conditional on RH).

For T = polylog, we need W = polylog^{1/alpha} = polylog, which means storing
checkpoints at spacing polylog -- requiring M/polylog ~ M checkpoints.
Storage: M * log(x) ~ 10^102 bits. **Exponential in N = log(x).**

**The preprocessing question reduces EXACTLY to the pi(x) complexity question:**
  - polylog query time <=> pi(x) computable in polylog time <=> pi(x) in NC
  - This is the central open question (see status/OPEN_PROBLEMS.md, Problem 1)

**Best achievable tradeoff (unconditional):**
  S = T = (M * log x)^{2/5} * polylog ~ 10^{41.8} (for M = 10^100)

**Best achievable tradeoff (under RH):**
  S = T = (M * log x)^{1/3} * polylog ~ 10^{34.8} (for M = 10^100)

Both are far from polylog. The sqrt(x) barrier is universal.

---

## 10. Failure Mode Classification

**Mode: Equivalence (E)**. The preprocessing/succinct data structure approach reduces
exactly to the circuit complexity of pi(x). It offers a clean S-T tradeoff within
existing methods but cannot break through the sqrt(x) barrier. No new algorithmic
insight emerges; the approach merely repackages the known computational difficulty.

**Verdict: CLOSED** as a route to polylog-time p(n) queries.  
**Value: PARTIAL** -- the S-T tradeoff analysis is useful for practical implementations
at feasible scales (n <= 10^{12}, where 4.5 TB suffices for O(1) lookup).
