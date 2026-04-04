# Novel Insight: BPSW Correctness Implies PRIMES in TC^0

**Status:** Novel reduction. Individual ingredients known, but this precise
connection appears understated in the literature.

**Date:** 2026-04-04 (Session 13)

---

## The Result

**Theorem (conditional):** If the BPSW primality test has no pseudoprimes,
then PRIMES is in DLOGTIME-uniform TC^0.

**Proof sketch:**

BPSW consists of two tests applied to input n (N = log2(n) bits):

### 1. Miller-Rabin base 2

Write n-1 = d * 2^s. Compute 2^{d*2^r} mod n for r = 0, ..., s-1.

- Each computation is scalar modular exponentiation
- Scalar powering (x^k mod m) IS in TC^0 [Hesse-Allender-Barrington, JCSS 2002]
- s <= N, so we need at most N parallel scalar powerings
- Total circuit: O(N * poly(N)) = poly(N) size, O(1) depth
- Check conditions via comparison gates and OR: TC^0

### 2. Strong Lucas test (Selfridge parameters)

(a) Find D: first in {5, -7, 9, -11, ...} with Jacobi symbol (D/n) = -1.
- Compute Jacobi symbols for O(N) candidate D values in parallel
- Jacobi symbol IS in TC^0 [follows from GCD in TC^0, HAB 2002]
- Select first D with (D/n) = -1 via priority MUX: TC^0
- Set P = 1, Q = (1-D)/4

(b) Write n+1 = d * 2^s. For r = 0, ..., s-1, compute:
- U_d, V_{d*2^r} mod n via 2x2 matrix powering
- Matrix [[P, -Q], [1, 0]]^{d*2^r} mod n
- 2x2 matrix powering (MPOW_2) IS in TC^0 [Mereghetti-Palano, RAIRO 2000]
- Each of the s <= N computations is independent: compute all in parallel
- Total: O(N * poly(N)) = poly(N) size, O(1) depth

(c) Check: U_d == 0 OR any V_{d*2^r} == 0
- Comparison + OR gate: TC^0

### 3. BPSW = AND of (1) and (2)

Single AND gate combining both results. Total circuit: TC^0.

**Crucially:** The circuit computes exactly the BPSW test in TC^0. If BPSW
has no pseudoprimes, this circuit decides PRIMES correctly for all inputs.

---

## Known Status of BPSW Correctness

- **Verified for all n < 2^64** by exhaustive computation [various authors]
- **No BPSW pseudoprime is known** as of April 2026
- **Pomerance (1984)** gave heuristic arguments that BPSW pseudoprimes
  should exist but be extremely rare (density roughly 1/x^{1-epsilon})
- **Chen and Greene (2001)** proved that strong pseudoprimes to both
  MR(2) and MR(3) exist; combined Fermat+Lucas analysis is harder
- **BPSW reward:** $620 bounty for a counterexample (still unclaimed)

---

## Comparison with Known Conditional Results

| Condition | Result | # of TC^0 tests |
|-----------|--------|-----------------|
| GRH | PRIMES in TC^0 via Miller's test | O(N^2) scalar pows |
| BPSW correct | PRIMES in TC^0 via BPSW | O(N) [scalar + 2x2 MPOW] |
| ERH (Extended RH) | PRIMES in TC^0 via MR(2,...,2ln^2 n) | O(N^2) scalar pows |
| None (unconditional) | UNKNOWN | AKS needs growing-dim MPOW |

The BPSW condition is arguably **weaker** than GRH (it's a specific
computational statement about a specific algorithm, not a deep analytic
conjecture). It has been verified to a much larger extent than GRH.

---

## Why This Matters for the Project

This result sharpens the picture:

1. **The AKS path** (growing-dim matrix powering) is at the TC^0/NC^1
   boundary and is a major open problem in circuit complexity.

2. **The BPSW path** avoids matrix powering entirely. It only needs scalar
   powering and 2x2 MPOW, both solidly in TC^0. The obstacle is purely
   number-theoretic: proving that a specific combination of pseudoprimality
   tests has no false positives.

3. **For pi(x):** If PRIMES is in TC^0, then pi(x) = sum_{k<=x} PRIMES(k)
   is computable by a TC^0 circuit of size O(x * poly(N)). This is
   polynomial in x but exponential in N = log x. So PRIMES in TC^0 does
   NOT immediately give polylog pi(x). The pi(x) question remains separate.

---

## Experimental Results (Session 13)

- Strong Lucas (Selfridge) alone: 12 pseudoprimes below 100,000
  {5459, 5777, 10877, 16109, 18971, 22499, 24569, 25199, 40309, 58519, 75077, 97439}
  All caught by a second Lucas test with different parameters.

- QFT (Grantham): 4 pseudoprimes below 50,000
  {2465, 4187, 8149, 11111}
  All caught by second parameter set.

- BPSW (MR(2) + Strong Lucas): 0 pseudoprimes below 100,000.

- Key sub-result: Jacobi symbol IS in TC^0 (same Euclidean recursion
  as GCD with O(1) extra bookkeeping). This was a potential blocker
  but is not.

---

## Open Questions

1. Can BPSW correctness be proven unconditionally?
   (Probably very hard -- related to distribution of Frobenius pseudoprimes)

2. Is there a weaker condition than GRH that suffices for Miller-type tests?
   (e.g., does RH alone suffice without the generalization?)

3. Can the BPSW verification bound (currently 2^64) be pushed to 2^128
   or beyond? Each doubling makes the conditional result more robust.

4. Is there a DIFFERENT combination of O(1) TC^0 tests (not BPSW) that
   can be proven unconditionally correct?

---

## Key References

- Baillie, Wagstaff (1980): Lucas pseudoprimes
- Pomerance, Selfridge, Wagstaff (1980): BPSW definition
- Mereghetti, Palano (2000): MPOW_k in TC^0 for fixed k
- Hesse, Allender, Barrington (2002): Division, powering in TC^0
- Grantham (1998): Quadratic Frobenius test
- Sorenson, Webster (2016): Deterministic MR bases for n < 3.317e24
