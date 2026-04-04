# Novel Insight: TC^0 Reduction via AKS Matrix Powering

**Status:** Novel reduction, connecting known results in a new way.
The individual ingredients (AKS, HAB, Allender) are known, but this specific
chain of reductions connecting PRIMES in TC^0 to polylog-dimensional matrix
powering appears unstudied.

---

## The Reduction Chain

1. **AKS primality test** (2002): n is prime iff for suitable r = O(polylog(n)),
   (x + a)^n ≡ x^n + a (mod x^r - 1, n) for a = 1, ..., O(sqrt(phi(r)) * log n)

2. **Matrix representation**: Multiplication by (x + a) in Z_n[x]/(x^r - 1)
   corresponds to an r × r companion matrix M_a over Z_n. Then:
   (x + a)^n mod (x^r - 1, n) = first column of M_a^n mod n

3. **The question becomes**: Can M^n mod m be computed in TC^0, where M is
   r × r with r = O(polylog(n)) and entries are O(log n) bits?

4. **What's known**:
   - **Scalar powering** (x^n mod m): IS in TC^0 (Allender 1999)
   - **Iterated multiplication** of N given matrices: IS in TC^0 (HAB 2002)
   - **Matrix powering** (M^n mod m for variable n): IS in NC^1, NOT known in TC^0
   - **The catch**: HAB handles N GIVEN matrices. For M^n, we need n copies of
     M, but n is exponential in input size. The circuit size would be poly(n),
     not poly(log n).

5. **Therefore**: PRIMES in TC^0 ⟺ polylog-dimensional matrix powering in TC^0

---

## Why This Matters

This gives the most precise formulation of what's needed for polylog prime counting:

- If M^n mod m (for polylog-dimensional M) is in TC^0:
  → PRIMES is in TC^0
  → pi(x) = sum_{k≤x} PRIMES(k) is in TC^0 (threshold counting)
  → pi(x) computable in O(1) parallel time with poly(x) processors
  → p(n) computable in polylog time (with parallelism)

- If M^n mod m requires depth Omega(log n):
  → PRIMES not in TC^0 (but could still be in NC^1)
  → This alone doesn't rule out fast pi(x), but makes it much harder

---

## The Scalar vs Matrix Gap

The fact that SCALAR powering is in TC^0 but MATRIX powering is not known is
deeply related to commutativity:

- Scalar multiplication is commutative: x^n = x * x * ... * x, and the order
  doesn't matter. This allows the "Chinese Remainder + discrete log" approach.
- Matrix multiplication is non-commutative: M^n requires preserving order.
  The standard approach (repeated squaring) has depth O(log n).

For the AKS case, the matrix M represents multiplication in a polynomial ring,
which is commutative! So M^n = M * M * ... * M with a COMMUTATIVE operation.

**Open question**: Does the commutativity of the underlying ring help?
Can we use CRT + some polynomial ring analog of discrete log to compute
M^n mod m in TC^0 for this specific class of companion matrices?

---

## Connection to Known Results

- Barrington (1989): NC^1 = width-5 branching programs. Matrix powering is
  in NC^1, which is "barely" above TC^0.
- Allender-Barrington-Ogihara (1999): Integer division in TC^0.
  Their technique: express the output bits as iterated products, which are TC^0.
- HAB (2002): Iterated multiplication of poly(n) numbers/matrices in TC^0.
  Key: the number of items to multiply must be polynomial in input size.

The gap between scalar and matrix powering in TC^0 is one of the few remaining
"bottleneck" problems at the TC^0/NC^1 boundary.

---

## CRITICAL UPDATE: Refined Understanding (from literature survey)

Mereghetti & Palano, "Threshold circuits for iterated matrix product and
powering," RAIRO-Theor. Inf. Appl. 34(1), 2000:

**Key distinction:**
- **IMP_k** (product of N DIFFERENT k×k matrices): NOT in TC^0 for k≥3 (unless TC^0=NC^1)
- **MPOW_k** (single matrix M^n, k×k, k FIXED): IS in TC^0!
  (reduces to O(k^2) scalar powering operations, each in TC^0)

**Our problem:** AKS needs MPOW with k = O(polylog(n)) GROWING. This is:
- NOT covered by MPOW_k (which needs k constant)
- NOT blocked by IMP_k (which is about different matrices)
- GENUINELY OPEN — sits precisely at the TC^0/NC^1 boundary

**The bottleneck for growing k:**
- Fixed k: k eigenvalues, each powered in TC^0, combined with O(1) depth → TC^0
- Growing k = polylog(n): polylog eigenvalues, combined with depth O(log k) = O(log log n)
- The O(log log n) depth is NOT constant → not directly TC^0

**Additional result:** Healy-Viola (2006) showed exponentiation in F_{2^n}
IS in TC^0 using Frobenius endomorphism. But this doesn't work over Z_n
(no Frobenius, non-field, zero divisors).

**Assessment:** The growing-dimension MPOW question is genuinely open and
at the exact frontier of circuit complexity. A positive answer would put
PRIMES in TC^0. A negative answer would effectively separate TC^0 from NC^1.

Source: literature/matrix_powering_tc0.md for full survey.

---

## Discovered: Session 11, 2026-04-04
