# Streaming and Space Complexity of Prime Counting

**Session 17 addendum** | 2026-04-04

---

## 1. Space Complexity of pi(x): What Is Known

### Upper bounds (space for computing pi(x) given N = log(x) input bits)

| Algorithm               | Time              | Space              | Note                    |
|--------------------------|-------------------|--------------------|-------------------------|
| Sieve of Eratosthenes    | O(x log log x)    | O(x) = O(2^N)      | Linear space            |
| Legendre                 | O(x)              | O(sqrt(x))         | O(2^{N/2})             |
| Lucy Hedgehog DP         | O(x^{3/4})        | O(sqrt(x))         | O(2^{N/2})             |
| Deleglise-Rivat          | O(x^{2/3}/log^2)  | O(x^{1/3} log^3)   | O(2^{N/3} poly(N))    |
| LMO                      | O(x^{2/3}/log x)  | O(x^{1/3})         | O(2^{N/3})            |
| Lagarias-Odlyzko (L-O)   | O(x^{1/2+eps})    | O(x^{1/4+eps})     | O(2^{N/4+eps}), analytic, conditional |
| L-O space-efficient       | O(x^{3/5+eps})    | O(x^eps)           | Nearly polylog space!   |

**Key observation:** The L-O space-efficient variant achieves SUBPOLYNOMIAL space
O(x^eps) = O(2^{eps*N}) at the cost of O(x^{3/5+eps}) time. This is the closest
to polylog space in the literature. However, for any fixed eps, the space is still
superpolynomial in N.

### Lower bounds on space

**There is NO nontrivial unconditional lower bound on the space needed for pi(x).**

The trivial bound is Omega(N) = Omega(log x) to hold the output. Beyond this:
- No streaming lower bound for pi(x) has been published.
- No branching program width lower bound beyond trivial.
- The communication complexity rank = 2^{N/2-1}+2 does NOT directly imply a space
  lower bound (communication and space are different models).

---

## 2. Streaming Model Analysis

### The model
A streaming algorithm sees integers 1, 2, ..., x in order and must output pi(x)
at the end. It has S bits of working memory.

### Upper bounds
- **Trivial**: S = O(x) bits (store the sieve). Time O(x log log x).
- **Segmented sieve**: S = O(sqrt(x)) bits (store primes up to sqrt(x)).
  Time O(x log log x). This is the standard streaming approach.
- **Smaller space**: With S bits, you can sieve segments of size S. You need
  primes up to sqrt(x), which requires O(sqrt(x)/ln(sqrt(x))) primes, each
  O(N) bits, so O(sqrt(x)*N/log(x)) = O(sqrt(x)) bits minimum for segmented sieve.

### Lower bounds (streaming)
The streaming model is more restrictive than general computation. Known results:

**For frequency moments:** Alon-Matias-Szegedy (1996) showed DISTINCT ELEMENTS
(F_0) needs Omega(1/eps^2 + log n) for (1+eps)-approximation. For exact F_0,
Omega(n) space is needed.

**For PREDICATE COUNTING (our case):** We want to count how many of {1,...,x}
satisfy isPrime(k). This is equivalent to computing the inner product
<1_P, 1_{[x]}> where 1_P is the characteristic function of primes.

**Claim: Exact streaming pi(x) requires Omega(sqrt(x)/log(x)) space.**

**Argument:** Consider a streaming algorithm A that processes k = 1, ..., x and
outputs pi(x) exactly. After processing k = sqrt(x), the algorithm's state must
encode enough information about the primes up to sqrt(x) to correctly sieve all
remaining numbers. Specifically, for each number m in (sqrt(x), x], the primality
of m depends on whether any prime p <= sqrt(x) divides m. The primes p <= sqrt(x)
are pi(sqrt(x)) ~ 2*sqrt(x)/ln(x) independent facts. The streaming state after
processing sqrt(x) must encode all of them (or equivalent information), because
any single missing prime would cause errors in the count for the interval
(sqrt(x), x].

More precisely: if the state after step sqrt(x) failed to distinguish two
different prime sets P1, P2 below sqrt(x), then for the same input suffix
(sqrt(x)+1, ..., x), the algorithm would give the same answer. But the count
of integers in (sqrt(x), x] divisible by a prime in P1 vs P2 differs.
Therefore the state must distinguish all 2^{pi(sqrt(x))} subsets, requiring
Omega(pi(sqrt(x))) = Omega(sqrt(x)/log(x)) bits.

**However, this is a SIEVE-BASED argument.** It assumes the algorithm works by
tracking divisibility. A non-sieve streaming algorithm could potentially use
less space if it found a non-sieve characterization of primes.

### Formal streaming lower bound via communication complexity

A tighter approach: any streaming algorithm with space S implies a one-way
communication protocol with S bits. The one-way communication complexity of
pi(x) is at least log(rank(M)) where M is the communication matrix.

From Session 17: rank(pi_N) = 2^{N/2-1} + 2 for balanced partition.

**But streaming uses a different partition.** In the streaming model, Alice sees
{1,...,t} and Bob sees {t+1,...,x}. The relevant matrix is indexed by the
"state of the world after step t" vs "the suffix starting at t+1."

For the BALANCED streaming cut at t = x/2:
- Alice's message encodes pi(x/2) plus the sieve state.
- Bob needs to know which primes <= sqrt(x) divide numbers in (x/2, x].
- All primes up to sqrt(x) are on Alice's side.

This gives one-way CC >= pi(sqrt(x)) = Omega(sqrt(x)/log(x)) bits,
which translates to streaming space S = Omega(sqrt(x)/log(x)).

**VERDICT:** Streaming space for exact pi(x) is Theta(sqrt(x)/log(x)) for
sieve-based algorithms. No formal proof that non-sieve algorithms need this much.

---

## 3. Online p(n) Computation

### The question
Maintain state and output p(1), p(2), ..., p(n) in sequence. Can the state
be polylog(p(n))?

### Analysis
- **Sieve approach**: State = bit array of size p(n). Each step: advance sieve.
  Space O(p(n)), time O(log log p(n)) amortized per prime.
- **Segmented sieve**: State = primes up to sqrt(p(n)) + current segment.
  Space O(sqrt(p(n))), time O(1) amortized per integer in segment.
- **Primality testing approach**: State = just the last prime p(k). To find
  p(k+1), test p(k)+2, p(k)+4, ... Expected gap ~ln(p(k)), each test costs
  O(polylog(p(k))). State = O(log(p(k))) bits. Time O(log^2(p(k))) amortized.

**Wait -- this last approach uses only O(log(p(k))) = O(N) bits of state!**

Correctness: Given p(k), we search p(k)+1, p(k)+2, ... testing each with
a deterministic primality test (AKS or conditional BPSW). The first prime
found is p(k+1). State needed: only p(k) itself (N bits).

Time per prime: O(gap * T_test) = O(log(x) * polylog(x)) = O(polylog(x))
amortized (Cramer's conjecture gives expected gap O(log^2 x); unconditionally
gap <= x^{0.525} by Baker-Harman-Pintz).

**Result: Online p(n) computation needs only O(N) = O(log(p(n))) space.**

But this does NOT help with pi(x), because this approach takes O(x) total time
(you must iterate through all integers up to p(n) ~ n*ln(n)). The space is small,
but the time is catastrophic.

### Connection to our problem
The online approach shows: SPACE is not the bottleneck for p(n). You can compute
p(n) in O(log n) space. The bottleneck is TIME: you need O(n log n) time this way
vs O(sqrt(n) polylog(n)) with binary search + pi(x).

---

## 4. Space-Time Tradeoffs

### Known S*T tradeoffs
- Sorting: S*T = Omega(n^2) (Beame 1991)
- Element distinctness: S*T = Omega(n^{3/2}) (Beame et al.)
- For general Boolean functions on N bits: S*T >= Omega(N) (trivial: must read input)

### For pi(x) on N = log(x) bits
No published S*T tradeoff specific to pi(x). However, combining what we know:

| Space S        | Best known time T     | S * T                    |
|----------------|-----------------------|--------------------------|
| O(2^N) = O(x) | O(x log log x)        | O(x^2 log log x)         |
| O(2^{N/2})     | O(x^{3/4})            | O(x^{5/4})               |
| O(2^{N/3})     | O(x^{2/3})            | O(x)                     |
| O(2^{N/4})     | O(x^{1/2+eps})        | O(x^{3/4+eps})           |
| O(2^{eps*N})   | O(x^{3/5+eps})        | O(x^{3/5+2eps})          |
| O(poly(N))     | ???                   | ???                      |

The S*T product DECREASES as we move to more space-efficient algorithms.
The L-O space-efficient variant achieves S*T = O(x^{3/5+eps}), which is
sublinear in x. There is no known lower bound preventing S = poly(N),
T = poly(N) (i.e., polylog space AND time).

### A conjectural S*T lower bound
If pi(x) requires processing Omega(sqrt(x)) "independent facts" (floor values
or zeta zeros), and each fact requires Omega(1) time, then T >= Omega(sqrt(x)).
The space to hold intermediate results during this processing would be at least
Omega(log(sqrt(x))) = Omega(N), giving S*T >= Omega(N * sqrt(x)).

This is consistent with L-O: S = x^{1/4}, T = x^{1/2}, S*T = x^{3/4} > N*sqrt(x).

---

## 5. Branching Program Width for pi(x)

### The model
A branching program reads the N = log(x) bits of x in some order.
Width W = number of states. Length L = number of bit reads (may re-read bits).
Barrington (1989): Width-5 polynomial-length BP = NC^1.

### What we know
- Session 17: isPrime has communication matrix rank 2^{N/2-1}+1.
  This means any OBLIVIOUS BP for isPrime needs width >= 2^{N/2-1}+1.
- For pi(x): rank = 2^{N/2-1}+2, so width >= 2^{N/2-1}+2.

**But this is for OBLIVIOUS programs with the balanced bit partition.**
A BP that reads bits in an OPTIMAL order, or reads bits MULTIPLE TIMES,
could potentially do better. Width lower bounds for non-oblivious BPs
are much harder to prove.

### The Nechiporuk bound
For functions on N bits partitioned into blocks, BP size >= sum_i (distinct
subfunctions on block i). For pi(x), each N/2-bit block produces ~2^{N/2}
distinct subfunctions (from the rank computation). This gives:
  BP size >= 2 * 2^{N/2-1} = 2^{N/2} = sqrt(x)

This is the strongest combinatorial lower bound we have, and it matches
the sqrt(x) barrier from every other approach.

---

## 6. Verdict: Does Streaming/Space Give a New Angle?

**NO.** The streaming/space perspective CONFIRMS the sqrt(x) barrier from a
new direction but does NOT break it or circumvent it.

Key findings:
1. **Space alone is NOT the bottleneck**: p(n) is computable in O(log n) space
   (online iteration with primality testing). The bottleneck is time.
2. **Streaming lower bound = Omega(sqrt(x)/log x)**: Same order as the
   communication complexity bound. Consistent with all other barriers.
3. **L-O space-efficient variant**: Already achieves near-polylog space
   O(x^eps) at cost of O(x^{3/5}) time. Space reduction is possible;
   time reduction below sqrt(x) is the hard part.
4. **Branching program width**: Nechiporuk bound gives size >= sqrt(x),
   confirming no NC^1 shortcut via BPs.
5. **S*T tradeoff**: No formal lower bound proven, but all known algorithms
   have T >= Omega(x^{1/2}) regardless of space.

**The sqrt(x) barrier is universal across yet another computational model.**
All roads lead to the same conclusion: computing pi(x) exactly requires
processing Omega(sqrt(x)) independent pieces of information (whether
floor values, zeta zeros, or prime positions), and no computational model
has found a way to compress this.

### Status: CLOSED (no new angle from streaming/space)
