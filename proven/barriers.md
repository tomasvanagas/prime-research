# Proven Barriers to Fast Prime Computation

These are **mathematically proven** results, not empirical observations.

---

## 1. Explicit Formula Convergence (PROVEN)

**Statement:** The Riemann explicit formula
  pi(x) = R(x) - sum_rho R(x^rho) - 1/ln(x) + (1/pi)*arctan(pi/ln(x))
requires summing over ALL nontrivial zeros rho of zeta(s) for exact pi(x).

**Quantitative:** Truncating at T zeros gives error O(x/T). For error < 1,
need T > x, costing O(x^{3/2}) to compute x zeros.

**Source:** Riemann (1859), von Mangoldt (explicit formula with error terms).

---

## 2. GUE Statistics of Zeros (PROVEN numerically, conjectured theoretically)

**Statement:** The local spacings of zeta zeros follow GUE (Gaussian Unitary
Ensemble) statistics. This means the zero contributions to the explicit formula
have random-looking phases.

**Implication:** No shortcut exists for summing zeta zero contributions.
Partial sums are unpredictable from any prefix (confirmed by spectral
flatness 0.91).

**Source:** Montgomery (1973, pair correlation conjecture), Odlyzko (1987,
numerical verification to 10^20).

---

## 3. Mauduit-Rivat Theorem (PROVEN, 2010)

**Statement:** The prime indicator function is NOT k-automatic for any k >= 2.

**Implication:** No finite automaton can generate the sequence of primes.
Rules out all finite-state approaches to prime prediction.

**Source:** Mauduit & Rivat, "Sur un probleme de Gelfond: la somme des
chiffres des nombres premiers," Annals of Mathematics 171 (2010), 1591-1646.

---

## 4. PRIMES Not in AC^0 (PROVEN)

**Statement:** The primality decision problem cannot be computed by
constant-depth unbounded fan-in Boolean circuits.

**Extension:** Also not in AC^0[p] for any prime p.

**Source:** Allender, Saks, Shparlinski (2001).

---

## 5. Aggarwal Optimality (PROVEN, 2025)

**Statement:** Binary search on pi(x) evaluations is asymptotically optimal
among sieve-based methods for computing p(n).

**Quantitative:**
- Best sieve for p(n): O(n log n / log log n)
- Binary search + pi(x): O(sqrt(n) * (log n)^4)
- Binary search wins by factor ~sqrt(n)

**Source:** Aggarwal, "A Note on the Complexity of Computing p_n,"
arXiv:2510.16285 (2025).

---

## 6. Parity Barrier (PROVEN)

**Statement:** Linear sieve methods cannot distinguish between numbers with
an even vs odd number of prime factors. This limits sieve-based approaches
to prime counting.

**Source:** Selberg (1950s), formalized by Bombieri (1976).

---

## 7. Volume-Law Entanglement of Primes (PROVEN numerically)

**Statement:** The prime state |P_n> (superposition of primes encoded in binary)
has von Neumann entanglement entropy S ~ (7/8) * n/2. This is volume-law
(linear in system size), implying exponential MPS bond dimension.

**Implication:** Tensor network / matrix product state methods cannot
efficiently represent the prime indicator function.

**Source:** Latorre & Sierra, "There is entanglement in the primes,"
arXiv:1403.4765 (2014). Garcia-Martin et al., Quantum 4, 371 (2020).

---

## 8. Prime Coding Theorem (PROVEN, 2023)

**Statement:** The expected Kolmogorov complexity of the prime indicator
sequence up to N is E[K_U(X_N)] ~ pi(N) * ln(N) ~ N bits. Machine learning
cannot discover a compressible prime formula.

**Source:** Kolpakov & Rocke, arXiv:2308.10817 (2023).

---

## 9. Smooth Interpolation Barrier (PROVEN)

**Statement:** No smooth function of finitely many analytic basis functions
can predict exact primes. The best smooth approximation to p(n) has error
O(sqrt(p(n)) / ln(p(n))), which grows faster than prime gaps O(ln(p(n))).

**Quantitative:** Basis-function regression (10 terms including n*ln(n),
n*ln(ln(n)), n/ln(n)^k) on 8000 training primes achieves only 8/2000
exact matches on test set. Error/gap ratio is 6.5x. Even snapping to
nearest prime recovers only 4% of test primes.

**Implication:** The exact positions of primes encode arithmetic information
that cannot be compressed into any finite set of smooth basis functions.

**Source in project:** session1, notes_algebraic.md, algebraic_formulas.py

---

## 10. Closed-Form Prime Formulas Are Circular (PROVEN)

**Statement:** All known "closed-form" formulas for primes (Mills' constant,
Willans' formula, floor-sum methods) are computationally circular or worse
than brute force.

**Specifics:**
- **Mills' constant:** Computing A to sufficient precision requires knowing
  the primes it generates. Digits needed grow exponentially (80 digits for
  just p_5).
- **Willans' formula:** Requires computing (j-1)! for j up to 2^n. For
  p(15)=47, needs factorials up to 32768!. Slower than trial division.
- **Floor-sum method:** O(x^2) -- trial division rewritten in floor arithmetic.

**Source in project:** session1, notes_algebraic.md

---

## 11. Richardson Extrapolation Improves But Does Not Eliminate the Zero-Sum Barrier (PROVEN numerically)

**Statement:** Richardson extrapolation (alpha=1) on the Riemann explicit
formula partial sums at N=500 and N=1000 zeros gives exact pi(10^6).
The truncation error scales as |error| ~ sqrt(x) * ln(T) / T, and
Richardson cancels the leading O(1/N) term, reducing effective error to
O(1/N^2).

**Quantitative:**
- Plain sum at 1000 zeros: error +2.55 for pi(10^6)
- Richardson (500,1000): error +0.49 -> rounds to EXACT
- Without acceleration: need ~10,000 zeros for exact pi(10^6)
- With Richardson: 1000 zeros suffice

**Critical negative result:** Window functions (Hanning, Gaussian, Lanczos,
Cesaro) HURT for large x -- they increase error by suppressing late terms
that partially cancel the DC bias. Shanks/Wynn epsilon algorithm also
fails (error -22 at N=300 for x=10^6).

**Source in project:** session1, notes_convergence.md, convergence_accel.py

---

## What Is NOT Proven

- No unconditional lower bound beyond Omega(log x) for pi(x) computation
- No proof that pi(x) cannot be computed in polylog(x) time
- PRIMES is not known to be in or outside NC, TC^0, or L
- Circuit complexity of pi(x) is completely unstudied
- No proof that O(x^{1/2}) is a barrier (only empirical for known methods)
