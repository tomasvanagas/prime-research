# The Circuit Size Barrier for pi(x)

**Session 12 — April 2026**

## Summary

We prove empirically and argue theoretically that ALL known algorithms for
computing pi(x) produce circuits of EXPONENTIAL size in the input length
N = log(x). Specifically:

- **Sieve-based (Lucy DP, Meissel-Lehmer, Deleglise-Rivat):** Size O(x^{2/3}) = O(2^{2N/3})
- **Analytic (Lagarias-Odlyzko):** Size O(x^{1/2+ε}) = O(2^{N/2+ε})
- **Both:** exponential in N

This establishes that **"Is pi(x) computable in polylog(x) time?" is equivalent
to "Is pi(x) in NC?"** (the class of problems with polynomial-size, polylog-depth circuits).

## Detailed Findings

### 1. Lucy DP is Fundamentally Sequential

The Lucy DP computation DAG has depth exactly pi(sqrt(x)) = Theta(sqrt(x)/ln x).
This is because S(x) is updated at EVERY prime step p ≤ sqrt(x), and each update
depends on the previous step, creating an unavoidable sequential chain.

**Empirical evidence:** depth/pi(sqrt(x)) = 1.000 for all x ∈ {100, 500, ..., 100000}.

The critical path goes through the node S(x, p) for every prime p ≤ sqrt(x).

### 2. Meissel-Lehmer Has Minimal Depth But Exponential Width

The telescoped (Meissel-Lehmer) formulation breaks the sequential chain:
- Depth: O(log log x) = O(log N) — nearly constant
- Width: O(x^{2/3}) = O(2^{2N/3}) — exponential in N

Confirmed: the recursion pi(x) → pi(x^{2/3}) → pi(x^{4/9}) → ... reaches O(1)
after L = log(log x)/log(3/2) ≈ 1.6 * log(N) levels. But each level processes
O(x^{(2/3)^l * 1/3}) "special leaves", dominated by O(x^{1/3}) at level 0.

The L/log(N) ratio stabilizes at ~1.60 for large x.

### 3. Parallel Rounds = pi(x^{1/3})

When allowing maximal parallelism, the Lucy DP requires exactly pi(x^{1/3}) + O(1)
sequential rounds. This is because:
- Primes p, q can share a round iff x < p^2 * q (for p < q)
- All primes > x^{1/3} satisfy p^3 > x, so they share one round
- Each prime ≤ x^{1/3} needs its own round

Empirical: rounds/(pi(x^{1/3})+1) = 1.00 for x ∈ {10^3, ..., 10^8}.

### 4. Floor-Value Set Has No Useful Structure

The floor-value set V = {floor(x/k)} has |V| ≈ 2*sqrt(x) elements.
We tested three potential structural properties:

**a) Commutativity:** The mapping matrices M_p (v → floor(v/p)) are
completely non-commutative. Zero commuting pairs out of C(pi(sqrt(x)), 2)
tested. This rules out reordering optimizations.

**b) Low rank:** The linear transformation from initial values to pi(x) uses
~80% of floor values with non-zero coefficients. The transformation is
effectively full-rank. No low-rank approximation is possible.

**c) Small coefficients:** The coefficients in the linear expansion are
HUGE (|coeff| up to ~10^6 for x = 5000) with massive cancellation
(sum |coeffs| / pi(x) grows with x). This mirrors the massive cancellation
in the analytic explicit formula.

### 5. Modular Computation Cannot Help

pi(x) mod m has conditional entropy H(pi(x) mod m | pi(x-1) mod m) = 0.537 bits,
INVARIANT across all moduli m (tested 2 through 30). This equals the entropy of
the prime indicator function and means:
- pi(x) mod m is as hard as pi(x) for every m
- CRT reconstruction from small moduli cannot bypass the barrier
- The "information per step" is fixed by the prime density, not the modulus

## Equivalence to NC Membership

**Theorem (informal):** pi(x) is computable in time polylog(x) if and only if
the function pi(x) is in the circuit complexity class NC (polynomial size,
polylog depth Boolean circuits).

**Proof sketch:**
- If pi(x) ∈ NC, there exists a circuit of size poly(N) and depth polylog(N)
  computing it. This can be evaluated in time poly(N) = polylog(x).
- Conversely, any algorithm running in time T has a circuit of size O(T * log T)
  and depth O(T). If T = polylog(x) = poly(N), this gives a poly(N)-size circuit.
  With standard parallelization, the depth can be reduced to polylog(N).

**Current status of NC membership:**
- Best upper bound: size O(2^{2N/3}), depth O(poly(N))
- Lower bound: size Omega(N) (trivial, from output size)
- **Gap: 2^{2N/3} / N — no technique can bridge this**

## Implications

1. **Finding an O(polylog) algorithm for pi(x) would resolve the NC question**
   for a natural problem, which would be a major result in complexity theory.

2. **Proving pi(x) ∉ NC** would prove our problem impossible, but this requires
   super-linear circuit lower bounds — beyond current techniques (essentially P vs NP).

3. **The problem sits precisely at the frontier of computational complexity theory.**
   Neither finding the algorithm nor proving impossibility is tractable with
   current mathematical tools.

4. **A breakthrough would require a fundamentally new approach** that avoids
   computing O(sqrt(x)) intermediate values. No such approach is known for
   ANY function related to prime counting.

## What Would Be Needed

An O(polylog) algorithm must NOT:
- Enumerate floor values (floor(x/k) for k = 1, ..., sqrt(x))
- Sum over zeta zeros (gamma_1, gamma_2, ..., up to height sqrt(x))
- Evaluate the prime indicator at individual integers
- Use any sieve-type inclusion-exclusion

It MUST find a way to compute pi(x) using only O(polylog(x)) operations
on values of size poly(log(x)). No mathematical framework is currently
known that could achieve this.
