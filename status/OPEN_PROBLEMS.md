# Open Problems: Viable Research Directions

Last updated: 2026-04-04

These are the ONLY directions not yet proven closed. Everything else has been
tested (380+ approaches) and confirmed to hit one of three failure modes:
Circularity, Equivalence, or Information Loss.

---

## 1. Circuit Complexity of pi(x) [MOST PROMISING]

**Question:** What is the circuit complexity of the prime-counting function?

**What's known:**
- PRIMES (decision) is NOT in AC^0, NOT in AC^0[p] (Allender-Saks-Shparlinski 2001)
- PRIMES is in P (AKS 2002), but NOT known to be in NC, TC^0, or L
- pi(x) (counting) has NO circuit lower bounds beyond trivial Omega(log x)
- Division and iterated multiplication ARE in uniform TC^0 (Hesse-Allender-Barrington 2002)

**Why it matters:** If pi(x) is computable by polylog-depth, poly-size circuits,
then polylog time is possible (with parallelism). If not, it proves impossibility.
This is genuinely unstudied territory.

**The gap:** Omega(log x) proven lower bound vs O(x^{1/2+epsilon}) best upper bound.
This is "one of the least-explored gaps in complexity theory" for a natural problem.

**Approach:** Try to show pi(x) mod 2 requires super-constant depth, or conversely
find a TC^0 reduction from pi(x) to known TC^0-computable functions.

---

## 2. Time-Bounded Kolmogorov Complexity of delta(n)

**Question:** What is Kt(delta(n) | n) -- the time-bounded conditional complexity
of the prime correction term?

**What's known:**
- Unbounded: K(p(n)|n) = O(1) trivially (program that sieves)
- delta(n) = p(n) - R^{-1}(n) has O(log n) bits of information
- Computing those bits takes O(x^{2/3}) time with best known methods
- Oliveira (2019) connected Kt complexity to circuit lower bounds

**Why it matters:** If Kt(delta(n)|n) = O(polylog(n)), that would imply
small circuits exist, connecting to Problem 1. If Kt(delta(n)|n) = omega(polylog),
that would prove no fast algorithm exists.

---

## 3. Zeta Zero Compressibility

**Question:** Do the Riemann zeta zeros have exploitable global structure
that allows fast summation of sum_rho R(x^rho)?

**What's known:**
- Zeros follow GUE statistics (Montgomery-Odlyzko)
- Individual zeros are computable in O(t^{1/3+epsilon}) time (Turing-type methods)
- Spectral flatness of zero sequence: 0.91 (high but not 1.0 = white noise)
- No FMM-type acceleration works due to incommensurability of zero heights

**Why it matters:** The entire barrier rests on needing O(sqrt(x)) zeta zeros.
If the sum could be compressed/approximated using structure in the zeros,
the barrier collapses.

**Status:** R(n) is 24% more compressible than random, but residual incompressible.

---

## 4. Berry-Keating / Hilbert-Polya Hamiltonian

**Question:** Does there exist a concrete self-adjoint operator H whose eigenvalues
are the zeta zeros, AND can its spectrum be computed efficiently?

**What's known:**
- The Hilbert-Polya conjecture posits such an operator exists
- Berry-Keating proposed H = xp + px (not rigorous)
- Connes' 2026 paper advances the Weil quadratic form approach
- Even if H exists, QPE requires 10^51 zeros for p(10^100)

**Why it matters:** A quantum-simulable Hamiltonian with zeta zero eigenvalues
would reduce the problem to quantum eigenvalue estimation.

---

## 5. Novel Number-Theoretic Identity

**Question:** Is there an identity that relates sum_rho R(x^rho) to a
computable function of x and n, without enumerating zeros?

**What's known:**
- All known identities (Weil explicit, Selberg trace, etc.) are transformations
  of the same underlying information
- No "shortcut identity" has been found in 165+ years of analytic number theory

**Why it matters:** This would be the most direct path -- a formula that bypasses
the zero sum entirely. The least likely to succeed but highest impact if found.

---

## Priority Assessment

| Direction | Feasibility | Impact | Recommended Effort |
|-----------|-------------|--------|-------------------|
| Circuit complexity | Medium | High | PRIMARY focus |
| Kt complexity | Medium | High | Theoretical exploration |
| Zero compressibility | Low | Very High | Numerical experiments |
| Berry-Keating | Very Low | Very High | Literature monitoring |
| Novel identity | Very Low | Maximal | Serendipity only |
