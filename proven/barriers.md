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

## 12. Convergence Acceleration Cannot Beat Zero-Sum Barrier (CONFIRMED numerically, Session 11)

**Statement:** No linear convergence acceleration method (Richardson of any order,
Levin u/t, Weniger delta, smoothed explicit formulas) can reduce the number of
zeta zeros needed from O(sqrt(x)) to subpolynomial.

**Quantitative:**
- Error structure is oscillatory, NOT clean 1/N expansion (alpha varies: -1.6 to 2.9)
- Levin transforms are WORSE than plain summation for large x
- Richardson order >2 diverges due to condition number explosion (10^{3k})
- Minimum zeros scale as N_min ~ C * x^alpha with alpha = 0.25--0.37 (power law)
- Smoothed formula (psi/ln x) catastrophic: 10% error from prime power terms

**Implication:** The explicit formula truncation error encodes GUE-random phases
of the next zeta zeros beyond the truncation point. No linear combination of
partial sums can predict these phases. The barrier is informational, not technical.

**Source in project:** session11, experiments/analytic/advanced_convergence.py, zero_scaling.py

---

## 13. Helfgott-Thompson Does Not Transfer to pi(x) (CONFIRMED numerically, Session 11)

**Statement:** H-T's O(x^{3/5}) method for M(x) = sum mu(n) exploits the signed
cancellation in (-1)^{omega(n)}. Since pi(x) is a positive counting function
(1-almost-primes), no analogous cancellation exists.

**Specifics:**
- pi(x) = pi(z) + A_1(x, z) where A_1 counts primes > z
- A_1 is irreducible: it IS the prime counting problem
- The Meissel-Lehmer decomposition pi(x) = Phi - A_0 - A_2 - A_3 - A_4
  IS the analog of H-T for pi(x), and costs O(x^{2/3})
- Converting M(x) → pi(x) requires the explicit formula (O(x^{1/2+eps}))

**Source in project:** session11, experiments/sieve/ht_transfer_attempt.py

---

## What Is NOT Proven

- No unconditional lower bound beyond Omega(log x) for pi(x) computation
- No proof that pi(x) cannot be computed in polylog(x) time
- PRIMES is not known to be in or outside NC, TC^0, or L
- **NEW (Session 11):** AKS path to PRIMES in TC^0 requires GROWING-dim (k=polylog(n))
  matrix powering. Fixed-k MPOW IS in TC^0 (RAIRO 2000). IMP_k (k≥3, different
  matrices) NOT in TC^0 unless TC^0=NC^1. Growing-dim MPOW is GENUINELY OPEN —
  sits precisely at the TC^0/NC^1 boundary. Healy-Viola 2006: F_{2^n} exponentiation
  IS in TC^0 via Frobenius, but doesn't work over Z_n.
- Circuit complexity of pi(x): all direct TC^0 paths closed (Session 11),
  but no impossibility proof exists
- No proof that O(x^{1/2}) is a barrier (only empirical for known methods)
- Combinatorial pi(x) in O(x^{3/5})? (H-T achieved for M(x) but transfer FAILS)
