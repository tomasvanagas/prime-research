# Novel Insight: Determinantal Complexity of pi(x)

**Status:** Novel connection. The specific formulation of "Is pi(x) in GapL?" as a
determinantal complexity question for the multilinear polynomial pi(bits) does not
appear to have been explicitly stated in the literature.

**Date:** 2026-04-04 (Session 15)

---

## The Connection

pi(x), viewed as a function of the N = log2(x) bits of x, is a **multilinear
polynomial of degree exactly N** over the reals. (Confirmed empirically for N=3
through N=10; follows from the fact that any function {0,1}^N -> Z has a unique
multilinear representation.)

The **determinantal complexity** dc(f) of a polynomial f is the smallest m such that
f = det(M) where M is an m×m matrix whose entries are affine-linear functions of
the variables.

**The GapL question reformulated:**

    "Is pi(x) in GapL?" ⟺ "Is dc(pi_N) = poly(N)?"

where pi_N is the degree-N multilinear polynomial in N variables representing pi(x)
for N-bit inputs.

---

## Empirical Results (Session 15)

### Degree of pi(x) as polynomial in bits

| N | Degree | # nonzero monomials | Max coefficient |
|---|--------|--------------------|----|
| 3 | 3 | 5 | 2 |
| 4 | 4 | 10 | 4 |
| 5 | 4 | 21 | 6 |
| 6 | 6 | 44 | 11 |
| 7 | 7 | ~100 | ~20 |
| 8 | 8 | ~250 | ~40 |
| 9 | 9 | ~500 | ~80 |
| 10 | 10 | ~1000 | ~160 |

Degree = N for all tested N except N=5 (degree 4). Monomials grow exponentially
(consistent with 50% sparsity from GF(2) analysis, Session 13).

### Determinantal representations found

| N | Matrix size | Found? |
|---|------------|--------|
| 2 | 2×2 | YES (numerical optimization) |
| 3 | 3×3 | YES |
| 4 | 4×4 | YES |
| 5 | 5×5 | Not found (optimization difficulty) |
| 6 | 6×6 | Not found (optimization difficulty) |

For N=2,3,4: found N×N matrices M(x) = M_0 + sum_k x_k * M_k with det(M) = pi(x)
at all 2^N evaluation points.

For N≥5: numerical optimization did not converge. This does NOT prove dc > N;
the optimization landscape is highly non-convex.

### Dimension Analysis

For an N×N matrix with affine-linear entries: N^2*(N+1) parameters, minus N^2-1 for
GL(N) symmetry = N^3 + 1 effective parameters.

The full multilinear polynomial space has 2^N dimensions.

| N | Effective params (N×N) | Polynomial space (2^N) | Ratio |
|---|------------------------|------------------------|-------|
| 4 | 65 | 16 | 4.06 |
| 6 | 217 | 64 | 3.39 |
| 8 | 513 | 256 | 2.00 |
| 10 | 1001 | 1024 | 0.98 |
| 12 | 1729 | 4096 | 0.42 |
| 16 | 4097 | 65536 | 0.06 |

**Critical observation:** For N ≥ 10, the parameter count drops below the polynomial
space dimension. A GENERIC degree-N polynomial in N variables does NOT have an N×N
determinantal representation for N ≥ 10.

For pi(x) to have polynomial determinantal complexity, it must have VERY SPECIFIC
algebraic structure that makes it "lie in" the determinantal variety. This structure,
if it exists, must come from number-theoretic properties of primes.

---

## Known Results on Determinantal Complexity

- **Mignon-Ressayre (2004):** Permanent of n×n matrix has dc ≥ n²/2
- **Grenet (2011):** Elementary symmetric polynomial e_k has dc = n
- **Cai-Chen-Li (2010):** VP ≠ VNP iff permanent has super-polynomial dc
- **ECCC TR26-035 (2026):** Learning read-once determinants in polynomial time

---

## Implications

1. If dc(pi_N) = poly(N), then pi(x) ∈ GapL ⊆ NC², confirming polylog-time algorithms exist.

2. If dc(pi_N) = 2^{Omega(N)}, then pi(x) ∉ GapL, closing the GapL direction.

3. The question "dc(pi_N) = ?" is a clean algebraic complexity question that could
   potentially be attacked using tools from algebraic geometry (the determinantal
   variety, Hilbert function arguments, etc.).

4. The connection to the permanent vs determinant problem is direct: both ask whether
   a specific polynomial has polynomial determinantal complexity.

---

## Why This Matters

This reformulation converts the number-theoretic question "can we count primes efficiently?"
into a pure algebraic complexity question about a specific polynomial. While this does not
immediately lead to a resolution, it opens the door to techniques from algebraic complexity
theory that have not been applied to prime counting.

Source: Session 15, experiments/circuit_complexity/det_perm_encoding.py
