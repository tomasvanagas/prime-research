# Three Failure Modes: A Taxonomy

Every one of 380+ failed approaches falls into exactly one of three modes.
This classification was developed empirically across Sessions 1-10.

---

## Mode C: Circularity (~60 approaches)

**Definition:** The method requires knowing primes to compute primes.

**Examples:**
- Mills' constant: the constant itself encodes all primes
- Generating functions: need primes to evaluate the EGF/OGF
- Chebotarev density: need p to compute Frobenius at p
- Euler product inversion: product is over primes
- Modular forms: Euler factors are indexed by primes
- Graph spectrum approaches: graph construction needs primes
- Expander graphs: Ramanujan graphs built from primes

**Pattern:** The formula or structure LOOKS like it avoids primes, but
some intermediate step requires prime input. The circularity is often
hidden behind several layers of abstraction.

**Detection:** Check whether any step requires evaluating a function AT
or OVER primes. If yes, circularity is present.

---

## Mode E: Equivalence (~160 approaches)

**Definition:** The method reduces to the explicit formula / zeta zero sum,
possibly in disguise.

**Examples:**
- All spectral methods (Selberg trace, heat kernel, Jacobi matrix)
- Perron contour integral = explicit formula
- Resurgent trans-series: non-perturbative sectors ARE zeta zeros
- Carlson's theorem: analytic continuation = explicit formula
- Geometric Langlands: same complexity as explicit formula
- Prismatic cohomology: O(x^{1/2+eps}) = analytic method
- Mertens function: O(x^{2/3}) = combinatorial method
- All Mobius inversion variants

**Pattern:** Different mathematical frameworks provide different
REPRESENTATIONS of the same underlying information. Converting between
"spectral side" (zeros) and "geometric side" (primes) always requires
O(sqrt(x)) or O(x^{2/3}) work.

**Detection:** Check whether the method's core computation involves
summing oscillatory terms indexed by zeta zeros (possibly implicit).

---

## Mode I: Information Loss (~150 approaches)

**Definition:** The method produces a smooth approximation that lacks
the ~170 bits needed to pinpoint the exact prime.

**Examples:**
- All ML/neural approaches (best 1.1%)
- All interpolation/regression on delta(n)
- Prime number theorem and its refinements
- Cramer model corrections
- Gap prediction (random walk, 5.04 bits/prime irreducible)
- Quantum holographic principle (volume-law, not area-law)
- p-adic analysis (not continuous in p-adic topology)

**Pattern:** The method captures the SMOOTH part of prime distribution
(which is O(polylog) to compute) but cannot capture the OSCILLATORY part
(which requires O(x^{2/3}) to compute). The gap between approximate and
exact is where all the difficulty lies.

**Detection:** Check whether the output is a continuous/smooth function
of the input. If yes, it loses the discontinuous prime structure.

---

## Why Only Three Modes?

Computing p(n) requires two things:
1. Accessing the right mathematical structure (primes, zeta zeros, etc.)
2. Extracting sufficient information from that structure

- If you can't access the structure without already knowing the answer: **C**
- If you access equivalent structure but at the same computational cost: **E**  
- If you access a simplified version that loses critical information: **I**

These are exhaustive because any approach must either:
- Use primes as input (C), or
- Use some transform of prime data, which is either equivalent (E) or lossy (I)

---

## Developed: Sessions 1-10, refined in Session 9
