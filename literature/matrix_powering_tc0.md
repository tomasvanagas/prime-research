# Matrix Powering in TC^0: Literature Survey

Last updated: 2026-04-04 (Session 11)

## Context

We showed: PRIMES in TC^0 iff polylog-dimensional matrix powering (M^n mod m,
where M is O(polylog(n)) x O(polylog(n)), the companion matrix of (x+a)^n in
Z_n[x]/(x^r-1)) is in TC^0.

This survey collects what is known about this question.

---

## 1. Scalar Powering IS in TC^0

**Hesse, Allender, Barrington.** "Uniform constant-depth threshold circuits for
division and iterated multiplication." JCSS 65, 2002.
[ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0022000002000259)

**Key result:** Integer division, scalar powering (x^k mod m), and iterated
multiplication of n given n-bit integers are all in DLOGTIME-uniform TC^0.

**Technique:** Chinese Remainder Representation (CRR). An n-bit number is stored
as its residues modulo polynomially many O(log n)-bit primes. The key insight:
raising x to the nth power modulo a small prime p can be expressed as a
first-order formula (FOM), which places it in uniform AC^0 -- even simpler than
TC^0. This is the POW predicate.

**Proof structure:**
1. Division is in FOM+POW (FOM augmented with a powering-modulo-small-primes predicate)
2. POW is in FOM (actually in FO = DLOGTIME-uniform AC^0)
3. Therefore division (and powering) are in DLOGTIME-uniform TC^0

**Why it works for scalars:** Scalar powering mod a small prime p reduces to
discrete logarithm in Z_p^*, which is a cyclic group of order p-1. The key is
that p has only O(log n) bits, so the entire group is small enough to handle
in constant depth.

---

## 2. Iterated Multiplication of N GIVEN Matrices IS in TC^0

**Immerman, Landau.** "The complexity of iterated multiplication." Information
and Computation 116(1), 1995.
[PDF](https://people.cs.umass.edu/~immerman/pub/mult.pdf)

**Key result:** Iterated multiplication of n given integers (a_1 * a_2 * ... * a_n)
is in TC^0. Uses CRR with discrete logarithms.

**Hesse, Allender, Barrington (2002):** Same paper as above, tightened to
DLOGTIME-uniform TC^0.

**IMPORTANT DISTINCTION:** This is for scalar integers only. For matrices, the
situation is different -- see section 3.

---

## 3. Iterated Matrix Product (IMP_k): Mostly NOT in TC^0

**Mereghetti, Palano.** "Threshold circuits for iterated matrix product and
powering." RAIRO - Theoretical Informatics and Applications 34(1), 39-46, 2000.
[NUMDAM](https://www.numdam.org/article/ITA_2000__34_1_39_0.pdf)
[RAIRO](https://www.rairo-ita.org/articles/ita/abs/2000/01/ita9910/ita9910.html)

**Key results:**
- IMP_k (iterated product of N given k x k matrices) is in NC^1 for all k
- For k >= 3: IMP_k is NOT in TC^0, unless TC^0 = NC^1
- For stochastic 2x2 matrices: IMP_2 IS in TC^0
- For k >= 4 stochastic: IMP_k is NOT in TC^0, unless TC^0 = NC^1

**Technique:** Uses Barrington-Therien algebraic characterization of TC^0 via
finite monoids/groups. The word problem over a finite monoid M is in TC^0 iff
every simple group dividing M's syntactic monoid is abelian (i.e., the monoid is
"solvable"). For k >= 3, the matrix group contains non-solvable subgroups
(like PSL(2,5) = A_5), making IMP_k NC^1-complete under this characterization.

**Key reference for the algebraic characterization:**

**Barrington, Therien.** "Finite monoids and the fine structure of NC^1."
JACM 35(4), 1988.
[ACM](https://dl.acm.org/doi/10.1145/48014.63138)

---

## 4. Fixed-Dimension Matrix POWERING (MPOW_k) IS in TC^0

**Mereghetti, Palano (2000).** Same paper as above.

**Key result:** For any FIXED k, MPOW_k (computing M^n for a single k x k matrix
M raised to the nth power) IS in TC^0.

**Why powering is easier than iterated product:** When computing M^n, the matrix M
is fixed. The computation reduces to: for each entry (i,j) of M^n, compute a
degree-k^2 polynomial in the eigenvalues of M raised to the nth power. Since k is
a CONSTANT, these are a constant number of scalar powering operations, each of
which is in TC^0 by Hesse-Allender-Barrington.

**CRITICAL CAVEAT:** "For any k" means for any FIXED constant k. The circuit family
depends on k. When k = O(polylog(n)) GROWS with input size, this result does NOT
apply. There is a different TC^0 circuit for k=2, k=3, k=4, etc., but no single
TC^0 circuit family that handles k growing with n.

---

## 5. The Critical Gap: Growing-Dimension Matrix Powering

**Our problem:** M is an r x r companion matrix where r = O(polylog(n)). This is
NOT covered by Mereghetti-Palano's MPOW_k result because k is not constant.

**What IS known for growing dimension:**
- M^n for d x d matrices (d growing) is in NC^1 by Barrington's theorem
  (branching programs of width 5 simulate log-depth circuits)
- M^n for d x d matrices is in NC^2 (direct repeated squaring: O(log n) matrix
  multiplications, each doable in NC^1)
- M^n is in L (logspace) for any polynomial d (sequential repeated squaring)
- NOT known to be in TC^0 for any d = omega(1)

**The bottleneck:** The Mereghetti-Palano technique reduces MPOW_k to O(k^2)
scalar powering operations. When k = polylog(n), this is polylog(n)^2 = polylog(n)
operations. The problem is that these operations interact: each entry of M^n is
a polynomial combination of eigenvalue powers, and the number of eigenvalues
grows with k. The circuit depth needed to combine them is O(log k) = O(log log n),
which is NOT constant depth.

---

## 6. Barrington's Theorem and NC^1

**Barrington, D.A.M.** "Bounded-width polynomial-size branching programs recognize
exactly those languages in NC^1." JCSS 38(1), 1989.
[ScienceDirect](https://www.sciencedirect.com/science/article/pii/0022000089900378)

**Key result:** NC^1 = width-5 polynomial-size branching programs. The proof uses
the non-solvability of S_5 (the symmetric group on 5 elements).

**Relevance:** This places matrix powering (for growing dimension) firmly in NC^1,
since a branching program of width d that reads bits of n can simulate M^n.
The question is whether NC^1 operations can be pushed down to TC^0.

---

## 7. Finite Field Exponentiation in TC^0

**Healy, Viola.** "Constant-Depth Circuits for Arithmetic in Finite Fields of
Characteristic Two." STACS 2006.
[PDF](https://www.ccs.neu.edu/home/viola/papers/FieldOps.pdf)

**Key result:** For the explicit realization F_{2^n} = F_2[x]/(x^{2*3^l} + x^{3^l} + 1),
exponentiation alpha^k is computable by DLOGTIME-uniform poly(n,t)-size TC^0
circuits (where t = |k| is the bit-length of the exponent).

**Why this matters:** Exponentiation in F_{2^n} IS powering in a polynomial ring
of degree n. This shows that for specific polynomial rings over F_2, powering IS
in TC^0 even for growing dimension (the polynomial has degree n, corresponding to
an n x n companion matrix).

**Why it doesn't directly help:** This works over F_2, exploiting the Frobenius
endomorphism (x -> x^2 is the identity on coefficients in F_2). Our problem is
over Z_n[x]/(x^r - 1), which has different structure. The Chinese Remainder
approach in F_{2^n} uses that F_{2^n}^* is cyclic and its discrete log is
computable in TC^0 for specific irreducible polynomials. Over Z_n, the
multiplicative structure is more complex.

---

## 8. Algebraic Characterization of TC^0

**Barrington, Therien (1988):** A regular language is in (nonuniform) TC^0 iff
its syntactic monoid is solvable (every simple group divisor is abelian).

**Krebs, Lange, Reifferscheid.** "Characterizing TC^0 in Terms of Infinite Groups."
STACS 2005.
[Springer](https://link.springer.com/chapter/10.1007/978-3-540-31856-9_41)

**Key result:** Extended the algebraic characterization of TC^0 to infinite groups
using "typed monoids."

**Relevance:** The algebraic characterization suggests that the TC^0/NC^1 boundary
is controlled by the group theory of the monoid generated by the matrices. For
our companion matrix in Z_n[x]/(x^r-1), the relevant monoid structure depends on
the factorization of n, which is what we're trying to determine (circularity risk).

---

## 9. Recent Results (2024-2026)

**Andrews, Wigderson.** "Constant-Depth Arithmetic Circuits for Linear Algebra
Problems." FOCS 2024. [arXiv:2404.10839](https://arxiv.org/abs/2404.10839)

**Key result:** Polynomial GCD, discriminant, resultant, Bezout coefficients,
squarefree decomposition, and inversion of Sylvester/Bezout matrices are all
computable by constant-depth polynomial-size arithmetic circuits (AC^0_F) over
fields. Does NOT include general matrix powering.

**Andrews, Wigderson + Bhattacharjee, Kumar, Rai, Ramanathan, Saptharishi, Saraf.**
"Constant-depth circuits for polynomial GCD over any characteristic." 2025.
[arXiv:2506.23220](https://arxiv.org/abs/2506.23220)

Extended to any characteristic. Again, does not address matrix powering.

**Transformer circuit complexity (2024-2025):** Multiple papers showed that
constant-precision transformers are in TC^0, and cannot solve NC^1-complete
problems like balanced parentheses unless TC^0 = NC^1. This is relevant because
it shows the TC^0/NC^1 question has practical implications for AI architectures.

---

## 10. Summary: What This Means for Our Project

### The Precise Question
Can we compute M^n mod composite-m, where M is an r x r matrix with r = O(polylog(n)),
in DLOGTIME-uniform TC^0?

### What's Known
| Problem | In TC^0? | Reference |
|---------|----------|-----------|
| Scalar x^n mod m | YES | Hesse-Allender-Barrington 2002 |
| Product of N given integers | YES | Immerman-Landau 1995, HAB 2002 |
| Product of N given k x k matrices (k fixed, k >= 3) | NO (unless TC^0=NC^1) | Mereghetti-Palano 2000 |
| M^n for fixed k x k matrix (k constant) | YES | Mereghetti-Palano 2000 |
| M^n for d x d matrix (d growing with input) | UNKNOWN (in NC^1) | -- |
| alpha^k in F_{2^n} (specific irreducible) | YES | Healy-Viola 2006 |
| Polynomial GCD, discriminant, etc. | YES (algebraic AC^0) | Andrews-Wigderson 2024 |
| Determinant of n x n matrix | In NC^2, not known in TC^0 | Csanky 1976, Berkowitz 1984 |

### The Gap
Our problem falls in the "M^n for d x d matrix (d growing with input)" row.
This is the unresolved case. The key obstacle:

1. **Fixed-k MPOW uses scalar powering of eigenvalues.** When k is constant,
   there are k eigenvalues, each powered in TC^0. The combination is a
   constant-size computation.

2. **Growing-k MPOW requires combining polylog(n) eigenvalue powers.** This
   combination step requires depth O(log k) = O(log log n), which is NOT constant.

3. **The Healy-Viola trick for F_{2^n} exploits Frobenius.** Over F_2, the map
   x -> x^2 permutes coefficients. This eliminates the eigenvalue combination
   step entirely. No analog exists over Z_n.

4. **The algebraic obstruction:** If the matrix generates a non-solvable monoid,
   IMP is NC^1-complete and NOT in TC^0 (assuming TC^0 != NC^1). But MPOW
   (powering a single matrix) avoids this because powers of a single matrix
   generate a cyclic (hence abelian, hence solvable) subgroup. The question is
   whether this cyclic structure can be exploited when the dimension grows.

### The Commutative Case (Our Specific Problem)
Our companion matrix in Z_n[x]/(x^r-1) generates a COMMUTATIVE ring. This is
potentially easier than general matrix powering because:
- Commutative matrix algebras are simultaneously diagonalizable (over algebraic closure)
- The ring Z_n[x]/(x^r-1) decomposes via CRT into a product of Z_n[x]/(f_i(x))
  for cyclotomic factors f_i
- Each factor is a quotient of Z_n by a cyclotomic polynomial

However, this decomposition requires knowing the factorization structure of x^r-1
mod n, which depends on the multiplicative order of primes dividing n modulo r.
This creates a potential circularity with the primality-testing goal.

### Assessment
**The question "is polylog-dimensional matrix powering in TC^0?" is GENUINELY OPEN.**
It is not covered by any known result. It sits precisely at the TC^0/NC^1 boundary.
A positive answer would put PRIMES in TC^0. A negative answer (showing it requires
depth omega(1)) would separate TC^0 from NC^1, which is a major open problem.

This confirms that our reduction (PRIMES in TC^0 iff polylog matrix powering in TC^0)
has identified a problem at the exact frontier of circuit complexity.

---

## Key References (Chronological)

1. Barrington (1989): NC^1 = width-5 branching programs
2. Barrington, Therien (1988): Algebraic characterization of TC^0 via solvable monoids
3. Immerman, Landau (1995): Iterated integer multiplication in TC^0
4. Mereghetti, Palano (2000): MPOW_k in TC^0 for fixed k; IMP_k not in TC^0 for k>=3
5. Hesse, Allender, Barrington (2002): Division and powering in uniform TC^0
6. Healy, Viola (2006): Exponentiation in F_{2^n} in TC^0
7. Andrews, Wigderson (2024): Linear algebra problems in constant-depth arithmetic circuits
8. Bhattacharjee et al. (2025): Polynomial GCD in constant depth over any characteristic
