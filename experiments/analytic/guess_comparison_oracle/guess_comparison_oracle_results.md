# Guess-comparison oracle for p(n): measured geography (S491)

**Question (user, S491):** if exact p(n) is out of reach, can we at
least take a guess x and decide whether x is below or above the real
p(n)?

**Formal content:** the comparison [p(n) ≤ x?] ⟺ [π(x) ≥ n?] is the
*decision version* of prime counting. Binary search recovers p(n)
exactly from ~log x comparisons (E6.6 / Aggarwal), so a cheap
worst-case comparison oracle would be the full breakthrough — the
decision version is not an escape hatch. What IS true is that the
difficulty has a sharp geography, measured here at n up to 10⁶
(30 log-spaced samples, `run.log`):

1. **Far zone — polylog-decidable today.** If the guess is farther from
   the truth than the maximal oscillation (≈ √x·log x under RH), the
   sign of R(x) − n (computable in polylog) decides the comparison
   rigorously. This is exactly why ~the leading half of p(n)'s digits
   are already polylog-computable via R⁻¹(n).
2. **Biased guess — predictable bit.** Against x = li⁻¹(n): p(n) >
   li⁻¹(n) in **30/30** cases. The li-vs-π bias (Rubinstein–Sarnak
   1994: under RH+LI the log-density of π(x) > li(x) is ≈ 2.6×10⁻⁷;
   first sign change near 1.4×10³¹⁶) makes this single comparison bit
   essentially free at any reachable scale — including n = 10¹⁰⁰,
   whose p(n) ≈ 2.5×10¹⁰² sits far below the Skewes region.
3. **Centered guess — coin flip.** Against x = R⁻¹(n): above in 10/30,
   below in 20/30 — no usable bias (consistent with the mean-zero
   GUE-random oscillation of π − R; the li bias is exactly the
   ½·li(√x) term that R subtracts). **The better the guess, the closer
   the comparison bit is to a fair coin** — this is the project's
   information barrier restated as a decision problem.
4. **Hard-zone width = √p(n).** |p(n) − R⁻¹(n)|/√p(n): median 0.13,
   max 0.53 across the sweep. Comparisons inside this window each cost
   a full π(x) computation (best known x^{1/2+ε}); outside it they are
   polylog.

**Verification twist (S491 protocols):** although a comparison inside
the hard zone cannot be *computed* cheaply, it can be *certified*
cheaply — the Lucy-DP interactive proof
(`experiments/constructions/p12_sumcheck_pi_verification/`) verifies
exact π(x), hence the bit [π(x) ≥ n], with a ~100 KB unconditionally-
sound transcript; adding an ECPP primality certificate for y yields a
succinct certificate for the exact statement p(n) = y.

**Falsifiers:** a stable predictable bias in the R-guess comparison bit
(would imply structure beyond R — none seen); a li-guess sign flip at
reachable scale (contradicts Skewes localisation — none seen).

**Edges/closures cited:** E6.6 (Aggarwal binary search), the
barrier paragraph (√x hard zone), S224 (adjacent-problem framing),
novel/succinct_verification_of_pi.md (certified comparisons).
